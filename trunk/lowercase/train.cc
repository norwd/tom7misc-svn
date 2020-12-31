// This code was forked from ../chess/blind/unblind.cc, which
// came from ../../mtoz, which came from ../../redi,
// so check that for some development history / thoughts.
// (In this version, I cleaned up the code a bunch, added
// native support for dense layers, and made several performance
// improvements. This version should definitely supersede the
// previous!)

// This first experiment tries to predict the letter's shape as
// Bezier (only) curves. I don't really expect this to work but
// it may be funny/interesting.



// TODO: Output error history during training. Could just concatenate
// it to a file. Would be nice to get it over a large number of examples
// (more than the batch size) to reduce variance.

// TODO: Fix error display.

// TODO: 2x/3x for RGB/FLAT.

// TODO: Show timer breakdown in GUI.

// TODO: Now we pass stuff to the video thread at multiple different moments,
// so they can be desynchronized. Should do something to synchronize this?

// TODO: In preparation for having multiple models that we care about,
// would be good to have the network configuration completely stored within
// the serialized format, so that we can mix and match models (at least
// multiple ones in the same process, if not spread out over code versions).
//   - TODO: Network now stores width/height/channels; restructure so
//     that rendering uses these.
//   - TODO: Had to get rid of const field in Stimulation and Errors,
//     because we copy assign them to video. Maybe better to heap
//     allocate; it's fine to use the copy constructor.

#include "SDL.h"
#include "SDL_main.h"
#include "sdl/sdlutil.h"
#include "sdl/font.h"

#include <CL/cl.h>

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <cmath>
#include <chrono>
#include <algorithm>
#include <tuple>
#include <utility>
#include <set>
#include <vector>
#include <map>
#include <unordered_set>
#include <deque>
#include <shared_mutex>

#include "base/stringprintf.h"
#include "base/logging.h"
#include "arcfour.h"
#include "util.h"
#include "stb_image.h"
#include "stb_image_write.h"
#include "vector-util.h"
#include "threadutil.h"
#include "randutil.h"
#include "base/macros.h"
#include "color-util.h"
#include "image.h"
#include "lines.h"

#include "loadfonts.h"
#include "network.h"

#include "clutil.h"
#include "timer.h"
#include "top.h"
#include "font-problem.h"
#include "autoparallel.h"

#include "../bit7/embed9x9.h"

#define FONTCHARS " ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789`-=[]\\;',./~!@#$%^&*()_+{}|:\"<>?" /* removed icons */
#define FONTSTYLES 7

static constexpr int VERBOSE = 2;
// Perform somewhat expensive sanity checking for NaNs.
// (Beware: Also enables some other really expensive diagnostics.)
// XXX PERF turn off once it's working!
static constexpr bool CHECK_NANS = false;

using namespace std;

using uint8 = uint8_t;
// Better compatibility with CL.
using uchar = uint8_t;

using uint32 = uint32_t;
using uint64 = uint64_t;

using Contour = TTF::Contour;

// TODO: These should be owned by network.h and exposed as
// constants; we want the CPU implementation to agree with
// the GPU #defines!

// We use a different compiled kernel depending on the layer's
// transfer function (these are in the inner loops, so we want
// to avoid overhead). The implementation of the transfer
// function is straightforward. The derivative is given in
// terms of the *output* of the transfer function, because
// this is the most natural/efficient for the sigmoid, and
// can be done (a bit less naturally) for ReLU.
static const char *SIGMOID_FN =
  "#define FORWARD(potential) (1.0f / (1.0f + exp(-potential)))\n"
  // This wants to be given the actual output value f(potential).
  "#define DERIVATIVE(fx) (fx * (1.0f - fx))\n";

static const char *RELU_FN =
  "#define FORWARD(potential) ((potential < 0.0f) ? 0.0f : potential)\n"
  // This is normally given as x < 0 ? 0 : 1, but note that f(x)
  // tells us which side of 0 the input is on (retaining the
  // not-very-important ambiguity at exactly 0), anyway. So we define
  // it in terms of f(x) to maintain the same interface we use for
  // sigmoid.
  "#define DERIVATIVE(fx) ((fx < 0.0f) ? 0.0f : 1.0f)\n";

// Like RELU but slight slope in the "zero" region.
static const char *LEAKY_RELU_FN =
  "#define FORWARD(potential) ((potential < 0.0f) ? potential * 0.01f : potential)\n"
  // See note above.
  "#define DERIVATIVE(fx) ((fx < 0.0f) ? 0.01f : 1.0f)\n";

// Defined at the bottom.
static std::optional<string> GetExclusiveApp();

// For checks in performance-critical code that should be skipped when
// we're confident the code is working and want speed.
#if 0
# define ECHECK(a)
# define ECHECK_EQ(a, b)
# define ECHECK_LT(a, b)
# define ECHECK_GT(a, b)
# define ECHECK_LE(a, b)
# define ECHECK_GE(a, b)
#else
# define ECHECK(a) CHECK(a)
# define ECHECK_EQ(a, b) CHECK_EQ(a, b)
# define ECHECK_LT(a, b) CHECK_LT(a, b)
# define ECHECK_GT(a, b) CHECK_GT(a, b)
# define ECHECK_LE(a, b) CHECK_LE(a, b)
# define ECHECK_GE(a, b) CHECK_GE(a, b)
#endif

// Graphics.
#define FONTWIDTH 9
#define FONTHEIGHT 16
static Font *font = nullptr;
#define SCREENW 1920
#define SCREENH 1280
static SDL_Surface *screen = nullptr;

static constexpr int MAX_FONTS = 100'000;
static constexpr int ENOUGH_FONTS = 100;

// Thread-safe, so shared between train and ui threads.
static CL *global_cl = nullptr;

std::shared_mutex print_mutex;
#define Printf(fmt, ...) do {				\
    WriteMutexLock Printf_ml(&print_mutex);		\
    printf(fmt, ##__VA_ARGS__);				\
    fflush(stdout);					\
  } while (0);

template<class C>
static void DeleteElements(C *cont) {
  for (auto &elt : *cont) {
    delete elt;
  }
  cont->clear();
}

// Communication between threads.
static bool train_should_die = false;
std::shared_mutex train_should_die_m;
static bool train_done = false;
std::shared_mutex train_done_m;

static uint8 FloatByte(float f) {
  if (f <= 0.0f) return 0;
  if (f >= 1.0f) return 255;
  else return f * 255.0;
}

static constexpr float ByteFloat(uint8 b) {
  return b * (1.0 / 255.0f);
}

// The input is the coordinates of each contour.
//  - We normalize these so that they are nominally in [0,1].
//    (See ttf.h).
//  - To keep a regular structure, a contour starts with
//    the start coordinates, then consists only of cubic Bezier
//    splines. Each one is given as its control point, then as
//    the end point.
//  - Since a letter shape can have multiple contours, we need
//    some way to delimit them. Simplest way is to fix a maximum
//    number of contours, and points per contour, so that it's
//    representable in a matrix.

// Almost all the fonts (9100 out of 9135 in the first data set)
// have at most 3 contours per letter (and 3 is by far the most
// common maximum value, presumably for letters like B). So it
// seems like a good choice for this experiment.
//
// The three contours don't have to be the same length. If sorted
// by length, the 8795/9135 have 100 or fewer points in their
// longest path, 8770 have 25 or fewer in the second longest,
// and 8930 have 16 or fewer in the third.
//
// Each contour/row has:
//   2 floats for start pos
//   2*2 floats for each bezier curve
//   so (100 * 4 + 2) + (25 * 4 + 2) + (16 * 4 + 2) =
//      402 + 102 + 66 = 570
//
// Note that 0,0 is a valid coordinate and not that weird. What do
// we do when the path is shorter than the row allows? Note we
// also need to be able to recognize such a thing to render the
// output.
// Couple options:
//   - Just pad with zeroes or -1 or whatever.
//   - Duplicate points (making zero-length beziers) so that the
//     row is filled.
//       - could just duplicate the last point, but maybe better
//         to spread it out.
//   - Make the path longer by splitting beziers into multiple
//     segments.

// The output is the same shape.
// Additionally, consider adding a one-hot prediction of the letter
// being lowercased, to coax the network to distinguish between
// letters.

static constexpr int ROW0_MAX_PTS = 38;
static constexpr int ROW1_MAX_PTS = 14;
static constexpr int ROW2_MAX_PTS = 10;

// XXX constexpr c++20
static const std::vector<int> ROW_MAX_POINTS = {
  ROW0_MAX_PTS, ROW1_MAX_PTS, ROW2_MAX_PTS,
};

static constexpr int INPUT_LAYER_SIZE =
  // start point, then beziers as control point, end point.
  2 + ROW0_MAX_PTS * (2 + 2) +
  2 + ROW1_MAX_PTS * (2 + 2) +
  2 + ROW2_MAX_PTS * (2 + 2);

static constexpr int OUTPUT_LAYER_SIZE =
  INPUT_LAYER_SIZE +
  // one-hot hint for what letter is this?
  26;

static constexpr int NEIGHBORHOOD = 1;

enum UserRenderStyle : uint32_t {
  RENDERSTYLE_INPUTXY = RENDERSTYLE_USER + 1000,
  RENDERSTYLE_OUTPUTXY = RENDERSTYLE_USER + 1001,
};

// Nominal pixel height of a character when rendering.
// Character can exceed these bounds; it's even normal
// for this to happen width-wise for wide characters.
static constexpr int NOMINAL_CHAR_SIZE = 400;


// Network that lives entirely on the GPU, but can be copied back to
// the Network object.
struct NetworkGPU {
  // XXX DESTRUCTOR.
  NetworkGPU(CL *cl, Network *net) : cl(cl), net(net) {

    layers.resize(net->layers.size());
    for (int layer = 0; layer < net->layers.size(); layer++) {
      layers[layer].indices =
	MoveMemoryToGPU(cl->context, cl->queue, true,
			&net->layers[layer].indices);
      layers[layer].weights =
	MoveMemoryToGPU(cl->context, cl->queue, false,
			&net->layers[layer].weights);
      layers[layer].biases =
	MoveMemoryToGPU(cl->context, cl->queue, false,
			&net->layers[layer].biases);
    }

    inverted_indices.resize(net->inverted_indices.size());
    for (int layer = 0; layer < net->layers.size(); layer++) {
      inverted_indices[layer].start =
	MoveMemoryToGPUConst(cl->context, cl->queue,
			     net->inverted_indices[layer].start);
      inverted_indices[layer].length =
	MoveMemoryToGPUConst(cl->context, cl->queue,
			     net->inverted_indices[layer].length);
      inverted_indices[layer].output_indices =
	MoveMemoryToGPUConst(cl->context, cl->queue,
			     net->inverted_indices[layer].output_indices);
    }

    clFinish(cl->queue);
  }

  // Read the weights and biases (which is the only thing that can change) from
  // GPU back to the Network object. Not thread safe!
  void ReadFromGPU() {
    for (int layer = 0; layer < net->layers.size(); layer++) {
      ReadTo(layers[layer].weights, &net->layers[layer].weights);
      ReadTo(layers[layer].biases, &net->layers[layer].biases);
    }
    clFinish(cl->queue);
  }

  // Like CopyBufferFromGPUTo, but don't wait for the command to finish.
  template<class T>
  void ReadTo(cl_mem buf, vector<T> *vec) {
    CHECK_SUCCESS(
	clEnqueueReadBuffer(cl->queue, buf, CL_TRUE, 0,
			    sizeof (T) * vec->size(),
			    vec->data(),
			    // No wait-list or event.
			    0, nullptr,
			    nullptr));
  }

  struct Layer {
    // Const
    cl_mem indices;
    cl_mem weights;
    cl_mem biases;
  };

  struct InvertedIndices {
    // Const
    cl_mem start;
    // Const
    cl_mem length;
    // Const
    cl_mem output_indices;
  };

  vector<Layer> layers;
  vector<InvertedIndices> inverted_indices;

  CL *cl;
  Network *net;
 private:
  DISALLOW_COPY_AND_ASSIGN(NetworkGPU);
};

// Data on the GPU for a single example in a single training round. Can
// be reused across rounds.
struct TrainingRoundGPU {
  TrainingRoundGPU(CL *cl, const Network &net) : cl(cl), net(&net) {
    for (int i = 0; i < net.num_layers + 1; i++) {
      stimulations.push_back(
	  CreateUninitializedGPUMemory<float>(cl->context, net.num_nodes[i]));
    }

    for (int i = 0; i < net.num_layers; i++) {
      errors.push_back(
	  CreateUninitializedGPUMemory<float>(cl->context,
					      net.num_nodes[i + 1]));
    }

    expected =
      CreateUninitializedGPUMemory<float>(cl->context,
					  net.num_nodes[net.num_layers]);
  }

  void LoadInput(const vector<float> &inputs) {
    CHECK_EQ(inputs.size(), net->num_nodes[0]);
    if (CHECK_NANS) { for (float f : inputs) CHECK(!std::isnan(f)); }
    CopyBufferToGPU(cl->queue, inputs, stimulations[0]);
  }

  void LoadExpected(const vector<float> &values) {
    CHECK_EQ(values.size(), net->num_nodes[net->num_layers]);
    if (CHECK_NANS) { for (float f : values) CHECK(!std::isnan(f)); }
    CopyBufferToGPU(cl->queue, values, expected);
  }

  void ExportStimulation(Stimulation *stim) {
    CHECK_EQ(stim->values.size(), stimulations.size());
    for (int i = 0; i < stim->values.size(); i++) {
      CopyBufferFromGPUTo(cl->queue, stimulations[i], &stim->values[i]);
    }
  }

  void ExportErrors(Errors *err) {
    CHECK_EQ(err->error.size(), errors.size());
    for (int i = 0; i < err->error.size(); i++) {
      CopyBufferFromGPUTo(cl->queue, errors[i], &err->error[i]);
    }
  }

  // Copy (only) the final layer of the stimulation back to main memory.
  void ExportOutput(vector<float> *out) {
    CopyBufferFromGPUTo(cl->queue, stimulations.back(), out);
  }

  // num_nodes + 1 layers. 0th is input, final is the output.
  vector<cl_mem> stimulations;
  // num_nodes layers.
  vector<cl_mem> errors;
  // Size of final stimulation.
  cl_mem expected;

  ~TrainingRoundGPU() {
    for (cl_mem m : stimulations) {
      CHECK_SUCCESS(clReleaseMemObject(m));
    }
    for (cl_mem m : errors) {
      CHECK_SUCCESS(clReleaseMemObject(m));
    }
    CHECK_SUCCESS(clReleaseMemObject(expected));
  }

  CL *cl;
  const Network *net;
 private:
  DISALLOW_COPY_AND_ASSIGN(TrainingRoundGPU);
};

// XXX this is probably obsolete -- now we do more specializations
// of constants in each wrapper's constructor.
static std::pair<std::vector<cl_program>, std::vector<cl_kernel>>
MakeTransferKernels(CL *cl, const char *base_file, const char *function_name) {
  std::vector<cl_program> programs;
  std::vector<cl_kernel> kernels;
  string base_src = Util::ReadFile(base_file);
  for (int tf = 0; tf < NUM_TRANSFER_FUNCTIONS; tf++) {
    cl_program program;
    cl_kernel kernel;
    string kernel_src;
    switch (tf) {
    case SIGMOID: kernel_src += SIGMOID_FN; break;
    case RELU: kernel_src += RELU_FN; break;
    case LEAKY_RELU: kernel_src += LEAKY_RELU_FN; break;
    default:
      CHECK(false) << "Invalid transfer function " << tf;
    }
    kernel_src += base_src;
    std::tie(program, kernel) = cl->BuildOneKernel(kernel_src, function_name);
    programs.push_back(program);
    kernels.push_back(kernel);
  }
  return make_pair(programs, kernels);
}

struct ForwardLayerCL {
  explicit ForwardLayerCL(CL *cl, const Network &net) : cl(cl) {
    string base_src = Util::ReadFile("forwardlayer.cl");
    for (int layer = 0; layer < net.layers.size(); layer++) {
      const TransferFunction transfer_function =
	net.layers[layer].transfer_function;
      const int indices_per_node = net.layers[layer].indices_per_node;

      const bool dense = net.layers[layer].type == LAYER_DENSE;

      string kernel_src;
      switch (transfer_function) {
      case SIGMOID: kernel_src += SIGMOID_FN; break;
      case RELU: kernel_src += RELU_FN; break;
      case LEAKY_RELU: kernel_src += LEAKY_RELU_FN; break;
      default:
	CHECK(false) << "Invalid transfer function " << transfer_function;
      }

      StringAppendF(&kernel_src, "\n#define INDICES_PER_NODE %d\n",
		    indices_per_node);

      kernel_src += base_src;
      auto [program, kernel] =
	cl->BuildOneKernel(kernel_src,
			   dense ?
			   "ForwardLayerDense" :
			   "ForwardLayerSparse");
      // XXX don't save debugging stuff
      layer_kernels.emplace_back(program, kernel,
				 dense, indices_per_node, transfer_function);
    }

  }

  struct ForwardContext {
    ForwardContext(ForwardLayerCL *parent, NetworkGPU *net_gpu, int layer) :
      parent(parent), net_gpu(net_gpu), layer(layer) {
      indices = net_gpu->layers[layer].indices;
      weights = net_gpu->layers[layer].weights;
      biases = net_gpu->layers[layer].biases;
    }

    // TODO: Do we really want to share the same command queue across
    // threads? Presumably clFinish can't tell "this thread's
    // commands" apart from others, so we may be prematurely
    // waiting/running other thread's work.
    void Forward(TrainingRoundGPU *train) {
      ECHECK_LT(layer + 1, train->stimulations.size());

      CL *cl = parent->cl;

      cl_mem src_values = train->stimulations[layer];
      cl_mem dst_values = train->stimulations[layer + 1];

      Network *net = net_gpu->net;

      auto [program, kernel, dense, kernel_ipn, kernel_tf] =
	parent->layer_kernels[layer];

      // Sanity check we have the right kernel
      CHECK(dense == (net->layers[layer].type == LAYER_DENSE));
      CHECK(kernel_ipn == net->layers[layer].indices_per_node);
      CHECK(kernel_tf == net->layers[layer].transfer_function);

      // Can't have multiple threads setting a kernel's argument at one time.
      {
	WriteMutexLock ml(&parent->m);

	if (dense) {
	  // No indices parameter for dense layers.
	  CHECK_SUCCESS(clSetKernelArg(kernel, 0, sizeof (cl_mem),
				       (void *)&src_values));
	  CHECK_SUCCESS(clSetKernelArg(kernel, 1, sizeof (cl_mem),
				       (void *)&weights));
	  CHECK_SUCCESS(clSetKernelArg(kernel, 2, sizeof (cl_mem),
				       (void *)&biases));
	  CHECK_SUCCESS(clSetKernelArg(kernel, 3, sizeof (cl_mem),
				       (void *)&dst_values));
	} else {
	  CHECK_SUCCESS(clSetKernelArg(kernel, 0, sizeof (cl_mem),
				       (void *)&src_values));
	  CHECK_SUCCESS(clSetKernelArg(kernel, 1, sizeof (cl_mem),
				       (void *)&indices));
	  CHECK_SUCCESS(clSetKernelArg(kernel, 2, sizeof (cl_mem),
				       (void *)&weights));
	  CHECK_SUCCESS(clSetKernelArg(kernel, 3, sizeof (cl_mem),
				       (void *)&biases));
	  CHECK_SUCCESS(clSetKernelArg(kernel, 4, sizeof (cl_mem),
				       (void *)&dst_values));
	}

	size_t global_work_offset[] = { 0 };
        size_t global_work_size[] = { (size_t)(net->num_nodes[layer + 1]) };
	// Printf("Run FL Kernel.\n");
	Timer kernel_timer;
	CHECK_SUCCESS(clEnqueueNDRangeKernel(cl->queue, kernel,
					     // work dimensions
					     1,
					     // global work offset
					     global_work_offset,
					     // global work size
					     global_work_size,
					     // local work size
					     nullptr,
					     // no wait list
					     0, nullptr,
					     // no event
					     nullptr));
	clFinish(cl->queue);
	kernel_ms += kernel_timer.MS();
      }
    }

    ~ForwardContext() {
    }

    cl_mem indices;
    cl_mem weights;
    cl_mem biases;
    ForwardLayerCL *parent = nullptr;
    NetworkGPU *net_gpu = nullptr;
    const int layer;
    double kernel_ms = 0.0;
  };

  ~ForwardLayerCL() {
    for (auto &[p, k, dense_unused, ipn_unused, tf_unused] : layer_kernels) {
      CHECK_SUCCESS(clReleaseKernel(k));
      CHECK_SUCCESS(clReleaseProgram(p));
    }
  }

  CL *cl = nullptr;
  // Owned. Indexed by layer id.
  std::vector<std::tuple<cl_program, cl_kernel, bool, int, TransferFunction>>
  layer_kernels;

  std::shared_mutex m;
};

// Set the error values; this is almost just a memcpy.
struct SetOutputErrorCL {
  explicit SetOutputErrorCL(CL *cl) : cl(cl) {
    std::tie(programs, kernels) =
      MakeTransferKernels(cl, "setoutputerror.cl", "SetOutputError");
  }

  struct Context {
    Context(SetOutputErrorCL *parent, NetworkGPU *net_gpu) :
      parent(parent), net_gpu(net_gpu) {}

    void SetOutputError(TrainingRoundGPU *train) {
      CL *cl = parent->cl;

      const Network *net = net_gpu->net;
      // Errors are always computed for the output layer, which is the last one.
      const TransferFunction transfer_function =
	net->layers.back().transfer_function;

      cl_kernel kernel = parent->kernels[transfer_function];

      // All three memories here have num_nodes floats.
      int num_nodes = net->num_nodes[net->num_layers];
      cl_mem actual_outputs = train->stimulations.back();
      cl_mem expected = train->expected;
      cl_mem output_error = train->errors.back();

      // Can't have multiple threads setting a kernel's argument at one time.
      {
	WriteMutexLock ml(&parent->m);

	CHECK_SUCCESS(clSetKernelArg(kernel, 0, sizeof (cl_mem),
				     (void *)&actual_outputs));
	CHECK_SUCCESS(clSetKernelArg(kernel, 1, sizeof (cl_mem),
				     (void *)&expected));
	CHECK_SUCCESS(clSetKernelArg(kernel, 2, sizeof (cl_mem),
				     (void *)&output_error));

	size_t global_work_offset[] = { 0 };
        size_t global_work_size[] = { (size_t)(num_nodes) };

	CHECK_SUCCESS(clEnqueueNDRangeKernel(cl->queue, kernel,
					     // work dimensions
					     1,
					     // global work offset
					     global_work_offset,
					     // global work size
					     global_work_size,
					     // local work size
					     nullptr,
					     // no wait list
					     0, nullptr,
					     // no event
					     nullptr));
	clFinish(cl->queue);
      }
    }

   private:
    SetOutputErrorCL *parent = nullptr;
    NetworkGPU *net_gpu = nullptr;
  };

  ~SetOutputErrorCL() {
    for (auto &k : kernels)
      CHECK_SUCCESS(clReleaseKernel(k));
    for (auto &p : programs)
      CHECK_SUCCESS(clReleaseProgram(p));
  }

 private:
  CL *cl = nullptr;
  // Owned:
  std::vector<cl_program> programs;
  std::vector<cl_kernel> kernels;

  std::shared_mutex m;

  DISALLOW_COPY_AND_ASSIGN(SetOutputErrorCL);
};

// Propagate errors backwards. Note that errors flow from "dst" to "src".
struct BackwardLayerCL {
  explicit BackwardLayerCL(CL *cl, const Network &net) : cl(cl) {
    string base_src = Util::ReadFile("backwardlayer.cl");
    for (int src_layer = 0; src_layer < net.layers.size() - 1; src_layer++) {
      int dst_layer = src_layer + 1;
      const TransferFunction transfer_function =
	net.layers[src_layer].transfer_function;
      const int dst_indices_per_node = net.layers[dst_layer].indices_per_node;
      // The dest layer, but num_nodes offset by 1.
      const int dst_num_nodes = net.num_nodes[dst_layer + 1];
      const bool dst_dense = net.layers[dst_layer].type == LAYER_DENSE;

      string kernel_src;
      switch (transfer_function) {
      case SIGMOID: kernel_src += SIGMOID_FN; break;
      case RELU: kernel_src += RELU_FN; break;
      case LEAKY_RELU: kernel_src += LEAKY_RELU_FN; break;
      default:
	CHECK(false) << "Invalid transfer function " << transfer_function;
      }

      StringAppendF(&kernel_src,
		    "\n"
		    "#define DST_INDICES_PER_NODE %d\n"
		    "#define DST_NUM_NODES %d\n",
		    dst_indices_per_node,
		    dst_num_nodes);

      kernel_src += base_src;
      auto [program, kernel] =
	cl->BuildOneKernel(kernel_src,
			   dst_dense ?
			   "BackwardLayerDense" :
			   "BackwardLayerSparse");
      // XXX don't save debugging stuff
      layer_kernels.emplace_back(program, kernel,
				 dst_dense,
				 dst_indices_per_node, transfer_function);
    }
  }

  struct Context {
    Context(BackwardLayerCL *parent, NetworkGPU *net_gpu, int dst_layer) :
      parent(parent), net_gpu(net_gpu), dst_layer(dst_layer) {

      const int gap = dst_layer;
      // const int src_layer = dst_layer - 1;

      starts = net_gpu->inverted_indices[gap].start;
      lengths = net_gpu->inverted_indices[gap].length;
      inverted_index = net_gpu->inverted_indices[gap].output_indices;
      dst_weights = net_gpu->layers[dst_layer].weights;
    }

    void Backward(TrainingRoundGPU *train) {
      CL *cl = parent->cl;
      const Network *net = net_gpu->net;
      const int gap = dst_layer;
      const int src_layer = dst_layer - 1;

      auto [program_, kernel, dst_dense, dst_ipn, tf] =
	parent->layer_kernels[src_layer];
      // Sanity check that this was compiled with the right ipn / tf.
      CHECK(dst_ipn == net->layers[dst_layer].indices_per_node);
      CHECK(tf == net->layers[src_layer].transfer_function);
      CHECK(dst_dense == (net->layers[dst_layer].type == LAYER_DENSE));

      cl_mem src_output = train->stimulations[src_layer + 1];
      cl_mem dst_error = train->errors[dst_layer];

      // This is the source layer, but num_nodes is offset by one
      // since it includes the size of the input layer as element 0.
      int src_num_nodes = net->num_nodes[src_layer + 1];
      cl_mem src_error = train->errors[src_layer];

      CHECK_EQ(src_num_nodes, net->inverted_indices[gap].start.size());

      // Can't have multiple threads setting a kernel's argument at one time.
      {
	WriteMutexLock ml(&parent->m);

	// XXX don't set unused arguments in dense mode
	CHECK_SUCCESS(clSetKernelArg(kernel, 0, sizeof (cl_mem),
				     (void *)&starts));
	CHECK_SUCCESS(clSetKernelArg(kernel, 1, sizeof (cl_mem),
				     (void *)&lengths));
	CHECK_SUCCESS(clSetKernelArg(kernel, 2, sizeof (cl_mem),
				     (void *)&inverted_index));
	CHECK_SUCCESS(clSetKernelArg(kernel, 3, sizeof (cl_mem),
				     (void *)&dst_weights));
	CHECK_SUCCESS(clSetKernelArg(kernel, 4, sizeof (cl_mem),
				     (void *)&src_output));
	CHECK_SUCCESS(clSetKernelArg(kernel, 5, sizeof (cl_mem),
				     (void *)&dst_error));
	CHECK_SUCCESS(clSetKernelArg(kernel, 6, sizeof (cl_mem),
				     (void *)&src_error));

	size_t global_work_offset[] = { 0 };
	size_t global_work_size[] = { (size_t)src_num_nodes };
	Timer kernel_timer;
	CHECK_SUCCESS(clEnqueueNDRangeKernel(cl->queue, kernel,
					     // work dimensions
					     1,
					     // global work offset
					     global_work_offset,
					     // global work size
					     global_work_size,
					     // local work size
					     nullptr,
					     // no wait list
					     0, nullptr,
					     // no event
					     nullptr));
	clFinish(cl->queue);
	kernel_ms += kernel_timer.MS();
      }
    }

    ~Context() {}

    cl_mem starts, lengths, inverted_index, dst_weights;
    BackwardLayerCL *parent = nullptr;
    NetworkGPU *net_gpu = nullptr;
    const int dst_layer;
    double kernel_ms = 0.0;
  };

  ~BackwardLayerCL() {
    for (auto &[p, k, deleteme0, deleteme1, deleteme2] : layer_kernels) {
      CHECK_SUCCESS(clReleaseKernel(k));
      CHECK_SUCCESS(clReleaseProgram(p));
    }
  }

  CL *cl = nullptr;
  // Indexed by source layer index.
  // Note: Unlike others, these have the layer's parameters
  // baked in.
  std::vector<std::tuple<cl_program, cl_kernel, bool, int, TransferFunction>>
    layer_kernels;

  std::shared_mutex m;
};

struct UpdateWeightsCL {
  explicit UpdateWeightsCL(CL *cl, const Network &net) : cl(cl) {
    // Note that this one doesn't depend on the transfer function/derivative.

    string base_src = Util::ReadFile("updateweights.cl");
    for (int layer = 0; layer < net.layers.size(); layer++) {
      const int indices_per_node = net.layers[layer].indices_per_node;

      const bool dense = net.layers[layer].type == LAYER_DENSE;

      string kernel_src;
      StringAppendF(&kernel_src, "\n#define INDICES_PER_NODE %d\n",
		    indices_per_node);

      kernel_src += base_src;
      auto [program, kernel] =
	cl->BuildOneKernel(kernel_src,
			   dense ?
			   "UpdateWeightsDense" :
			   "UpdateWeightsSparse");
      // XXX don't save debugging stuff
      layer_kernels.emplace_back(program, kernel,
				 dense, indices_per_node);
    }
  }

  struct Context {
    Context(UpdateWeightsCL *parent, NetworkGPU *net_gpu, int layer) :
      parent(parent), net_gpu(net_gpu), layer(layer) {

      layer_indices = net_gpu->layers[layer].indices;
      layer_weights = net_gpu->layers[layer].weights;
      layer_biases = net_gpu->layers[layer].biases;
    }

    void Update(float learning_rate, TrainingRoundGPU *train, int layer) {
      CL *cl = parent->cl;

      // Really can't run these in parallel because of concurrent writes to net.
      WriteMutexLock ml(&parent->m);

      cl_mem layer_error = train->errors[layer];
      cl_mem layer_values = train->stimulations[layer];

      auto [program_, kernel, dense, ipn] =
	parent->layer_kernels[layer];
      // Sanity check that this was compiled with the right constants.
      CHECK(ipn == net_gpu->net->layers[layer].indices_per_node);
      CHECK(dense == (net_gpu->net->layers[layer].type == LAYER_DENSE));

      const int num_nodes = net_gpu->net->num_nodes[layer + 1];

      if (dense) {
	// Dense version does not need indices arg.
	CHECK_SUCCESS(clSetKernelArg(kernel, 0, sizeof (cl_float),
				     (void *)&learning_rate));
	CHECK_SUCCESS(clSetKernelArg(kernel, 1, sizeof (cl_mem),
				     (void *)&layer_error));
	CHECK_SUCCESS(clSetKernelArg(kernel, 2, sizeof (cl_mem),
				     (void *)&layer_values));
	CHECK_SUCCESS(clSetKernelArg(kernel, 3, sizeof (cl_mem),
				     (void *)&layer_weights));
	CHECK_SUCCESS(clSetKernelArg(kernel, 4, sizeof (cl_mem),
				     (void *)&layer_biases));
      } else {
	CHECK_SUCCESS(clSetKernelArg(kernel, 0, sizeof (cl_float),
				     (void *)&learning_rate));
	CHECK_SUCCESS(clSetKernelArg(kernel, 1, sizeof (cl_mem),
				     (void *)&layer_error));
	CHECK_SUCCESS(clSetKernelArg(kernel, 2, sizeof (cl_mem),
				     (void *)&layer_indices));
	CHECK_SUCCESS(clSetKernelArg(kernel, 3, sizeof (cl_mem),
				     (void *)&layer_values));
	CHECK_SUCCESS(clSetKernelArg(kernel, 4, sizeof (cl_mem),
				     (void *)&layer_weights));
	CHECK_SUCCESS(clSetKernelArg(kernel, 5, sizeof (cl_mem),
				     (void *)&layer_biases));
      }

      size_t global_work_offset[] = { 0 };
      size_t global_work_size[] = { (size_t)num_nodes };
      CHECK_SUCCESS(clEnqueueNDRangeKernel(cl->queue, kernel,
					   // work dimensions
					   1,
					   // global work offset
					   global_work_offset,
					   // global work size
					   global_work_size,
					   // local work size
					   nullptr,
					   // no wait list
					   0, nullptr,
					   // no event
					   nullptr));
      clFinish(cl->queue);
    }

    void Finish() {
      CL *cl = parent->cl;
      clFinish(cl->queue);
    }

    cl_mem layer_indices, layer_weights, layer_biases;
    UpdateWeightsCL *parent = nullptr;
    NetworkGPU *net_gpu;
    const int layer;
    double kernel_ms = 0.0;
  };


  ~UpdateWeightsCL() {
    for (auto &[p, k, dense_unused, ipn_unused] : layer_kernels) {
      CHECK_SUCCESS(clReleaseKernel(k));
      CHECK_SUCCESS(clReleaseProgram(p));
    }
  }

  CL *cl = nullptr;
  // Owned. Indexed by layer id. (is_dense, indices_per_node)
  std::vector<std::tuple<cl_program, cl_kernel, bool, int>>
  layer_kernels;

  std::shared_mutex m;
};

// Make indices. This assumes that nodes are 2D "pixel" data, where on
// each layer we have width[l] * height[l] pixels, with channels[l]
// nodes per pixel. Row-major order.
//
// We mostly sample from a Gaussian near each pixel, but:
//  - we reject duplicates (inefficient),
//  - we reject pixels off the image (doesn't make sense; wrapping around
//      would work but an image is not a torus)
//  - we require that a small neighborhood around the pixel is mapped
//      directly (especially the pixel itself; this preserves spatial
//      locality and makes sure we don't have any statically dead nodes).
//
// TODO: In addition to the spatial channels, it might be useful to
// have some region that is densely connected. For example, say we
// specify DENSE_PREFIX_ROWS, which are then subtracted from the
// implied "image", and each pixel is assigned indices from the
// neighborhood, then from dense_prefix_rows, then from a gaussian.
//
// TODO: Allow specifying a strategy for index assignment. For chess,
// taking full rows, columns, and diagonals is probably better than
// random gaussians!
//
// TODO: If the number of indices requested is close to the number
// available (or equal to it), then this approach is super inefficient.
// Could instead randomly delete indices. But at least we only do it
// once.
static void MakeIndices(ArcFour *rc, Network *net) {
  static_assert(NEIGHBORHOOD >= 0, "must include the pixel itself.");
  auto OneNode = [](ArcFour *rc, RandomGaussian *gauss,
		    int64 *rejected, int64 *duplicate,
		    int indices_per_node,
		    int src_width, int src_height, int src_channels,
		    int dst_width, int dst_height, int dst_channels,
		    int idx) -> vector<uint32> {

    CHECK(indices_per_node <= src_width * src_height * src_channels) <<
    "Can't get " << indices_per_node
    << " distinct indices from a layer with " <<
    src_width << " x " << src_height << " x " << src_channels <<
    " = " << (src_width * src_height * src_channels) << " sources";

    // Whenever we read the neighborhood, we include all source channels.
    CHECK((NEIGHBORHOOD * 2 + 1) * (NEIGHBORHOOD * 2 + 1) *
	  src_channels <= indices_per_node) <<
    "neighborhood doesn't fit in indices!";
    // Which pixel is this?
    const int dst_nodes_per_row = dst_width * dst_channels;
    [[maybe_unused]]
    const int c = idx % dst_channels;
    const int x = (idx % dst_nodes_per_row) / dst_channels;
    const int y = idx / dst_nodes_per_row;

    const double xf = x / (double)dst_width;
    const double yf = y / (double)dst_height;

    // Use hash set for deduplication; we re-sort for locality of access later.
    unordered_set<int> indices;
    // clips xx,yy if they are out of the image.
    // cc must be a valid channel index.
    auto AddNodeByCoordinates = [src_width, src_height, src_channels, &indices,
				 rejected, duplicate](int xx, int yy, int cc) {
      ECHECK_GE(cc, 0);
      ECHECK_LT(cc, src_channels);
      if (xx < 0 || yy < 0 || xx >= src_width || yy >= src_height) {
	++*rejected;
	return;
      }
      int idx = (yy * src_width * src_channels) + xx * src_channels + cc;
      ECHECK_GE(idx, 0);
      ECHECK_LT(idx, src_width * src_height * src_channels);
      auto p = indices.insert(idx);
      if (!p.second) ++*duplicate;
    };

    // Find the closest corresponding pixel in the src layer; add all
    // its channels.
    const int cx = round(xf * src_width);
    const int cy = round(yf * src_height);
    for (int ny = -NEIGHBORHOOD; ny <= NEIGHBORHOOD; ny++) {
      for (int nx = -NEIGHBORHOOD; nx <= NEIGHBORHOOD; nx++) {
        // Note that the pixel may be clipped.
	for (int nc = 0; nc < src_channels; nc++) {
	  AddNodeByCoordinates(cx + nx, cy + ny, nc);
	}
      }
    }

    CHECK_LE(indices.size(), indices_per_node);

    // XXX Select this dynamically based on how many unused nodes
    // are even left?
    #if 0
    static constexpr double stddev = 1 / 16.0;

    // Sample gaussian pixels.
    while (indices.size() < indices_per_node) {
      double dx = gauss->Next() * stddev;
      double dy = gauss->Next() * stddev;

      AddNodeByCoordinates((int)round((xf + dx) * src_width),
			   (int)round((yf + dy) * src_height));
    }
    #else

    // XXXXX
    int hood = NEIGHBORHOOD;
    while (indices.size() < indices_per_node) {
      hood++;
      const int cx = round(xf * src_width);
      const int cy = round(yf * src_height);
      for (int ny = -hood; ny <= hood; ny++) {
	for (int nx = -hood; nx <= hood; nx++) {
	  // In the interests of getting more spatial
	  // dispersion, only add one channel at random. As
	  // we expand we can try these multiple times, so
	  // pixels closer to the center are more likely to
	  // have all channels used.
	  int nc = RandTo(rc, src_channels);
	  AddNodeByCoordinates(cx + nx, cy + ny, nc);
	  if (indices.size() == indices_per_node) goto done;
	}
      }
    }
  done:;
    #endif

    CHECK_EQ(indices_per_node, indices.size());
    vector<uint32> ret;
    ret.reserve(indices_per_node);
    for (int idx : indices) {
      CHECK_GE(idx, 0);
      CHECK_LT(idx, src_width * src_height * src_channels);
      ret.push_back(idx);
    }
    return ret;
  };

  // This must access rc serially.
  vector<ArcFour *> rcs;
  for (int i = 0; i < net->num_layers; i++) rcs.push_back(Substream(rc, i));

  auto OneLayer =
    [&rcs, &OneNode, &net](int layer) {
      if (net->layers[layer].type == LAYER_DENSE) {
	// Dense layers have a predetermined structure.
	const int prev_num_nodes = net->num_nodes[layer];
	const int this_num_nodes = net->num_nodes[layer + 1];
	CHECK(net->layers[layer].indices_per_node == prev_num_nodes);
	vector<uint32> *layer_indices = &net->layers[layer].indices;
	CHECK(layer_indices->size() == this_num_nodes * prev_num_nodes);
	for (int n = 0; n < this_num_nodes; n++) {
	  for (int p = 0; p < prev_num_nodes; p++) {
	    (*layer_indices)[n * prev_num_nodes + p] = p;
	  }
	}

      } else {
	// Assign sparse layers randomly.
	CHECK(net->layers[layer].type == LAYER_SPARSE);

	const int indices_per_node = net->layers[layer].indices_per_node;
	Printf("Intializing %d indices for layer %d...\n",
	       indices_per_node, layer);
	vector<uint32> *layer_indices = &net->layers[layer].indices;
	CHECK_LT(layer + 1, net->width.size());
	CHECK_LT(layer + 1, net->channels.size());
	CHECK_LT(layer + 1, net->num_nodes.size());
	const int src_width = net->width[layer];
	const int src_channels = net->channels[layer];
	CHECK_EQ(0, net->num_nodes[layer] % (src_width * src_channels));
	const int src_height =
	  net->num_nodes[layer] / (src_width * src_channels);
	const int dst_width = net->width[layer + 1];
	const int dst_channels = net->channels[layer + 1];
	CHECK_EQ(0, net->num_nodes[layer + 1] % (dst_width * dst_channels));
	const int dst_height = net->num_nodes[layer + 1] /
	  (dst_width * dst_channels);
	RandomGaussian gauss{rcs[layer]};
	int64 rejected = 0LL, duplicate = 0LL;
	for (int node_idx = 0;
	     node_idx < dst_height * dst_width * dst_channels;
	     node_idx++) {
	  vector<uint32> indices =
	    OneNode(rcs[layer], &gauss, &rejected, &duplicate,
		    indices_per_node,
		    src_width, src_height, src_channels,
		    dst_width, dst_height, dst_channels,
		    node_idx);
	  // Sort them, for better locality of access later.
	  std::sort(indices.begin(), indices.end());
	  CHECK_EQ(indices_per_node, indices.size());
	  const int start_idx = node_idx * indices_per_node;
	  for (int i = 0; i < indices_per_node; i++) {
	    ECHECK_LT(i, indices.size());
	    ECHECK_LT(start_idx + i, layer_indices->size())
	      << "start " << start_idx
	      << " i " << i
	      << " indices size " << layer_indices->size()
	      << " indices per node "
	      << indices_per_node;
	    (*layer_indices)[start_idx + i] = indices[i];
	  }
	  if (node_idx % 1000 == 0) {
	    Printf("  %d. [%d/%d] %.1f%% (%lld rejected %lld dupe)\n",
		   layer,
		   node_idx, dst_height * dst_width * dst_channels,
		   (100.0 * node_idx) / (dst_height * dst_width * dst_channels),
		   rejected, duplicate);
	  }
	}
      }
      Printf("... done with layer %d.\n", layer);
    };


  ParallelComp(net->num_layers, OneLayer, 12);

  Printf("DeleteElements:\n");
  DeleteElements(&rcs);
  Printf("Exiting MakeIndices.\n");
}

// Randomize the weights in a network. Doesn't do anything to indices.
static void RandomizeNetwork(ArcFour *rc, Network *net) {
  auto RandomizeFloats = [](float mag, ArcFour *rc, vector<float> *vec) {
    RandomGaussian gauss{rc};
    for (int i = 0; i < vec->size(); i++) {
      (*vec)[i] = mag * gauss.Next();
    }
  };

  // This must access rc serially.
  vector<ArcFour *> rcs;
  for (int i = 0; i < net->num_layers; i++) rcs.push_back(Substream(rc, i));

  // But now we can do all layers in parallel.
  CHECK_EQ(net->num_layers, net->layers.size());
  ParallelComp(net->num_layers, [rcs, &RandomizeFloats, &net](int layer) {
    // XXX such hacks. How to best initialize?

    /*
    RandomizeFloats(powf(0.025f, (layer / 3.0f) + 1.0), rcs[layer], &net->layers[layer].biases);
    RandomizeFloats(1.0f / (net->layers[layer].indices_per_node * ((layer / 3.0f) + 5)),
		    rcs[layer], &net->layers[layer].weights);
    */

    for (float &f : net->layers[layer].biases) f = 0.0f;
    // RandomizeFloats(0.000025f, rcs[layer], &net->layers[layer].biases);
    // RandomizeFloats(0.025f, rcs[layer], &net->layers[layer].weights);

    // The more indices we have, the smaller initial weights we should
    // use.
    float base_weight = 1.0f;
    float gaussian_width = base_weight /
      (float)net->layers[layer].indices_per_node;
    RandomizeFloats(gaussian_width, rcs[layer], &net->layers[layer].weights);
  }, 12);

  DeleteElements(&rcs);
}


// These must be initialized before starting the UI thread!
static constexpr int NUM_VIDEO_STIMULATIONS = 6;
static constexpr int EXPORT_EVERY = 4;

static constexpr bool DRAW_ERRORS = false;

// XXX this needs to support multiple models somehow.
// right now it just overwrites with whatever the last call was !!
struct UI {
  UI() {
    current_stimulations.resize(NUM_VIDEO_STIMULATIONS);
    current_errors.resize(NUM_VIDEO_STIMULATIONS);
    current_expected.resize(NUM_VIDEO_STIMULATIONS);
  }

  std::shared_mutex video_export_m;
  int current_round = 0;
  double examples_per_second = 0.0;
  vector<Stimulation> current_stimulations;
  vector<Errors> current_errors;
  vector<vector<float>> current_expected;
  Network *current_network = nullptr;
  bool allow_updates = true;
  bool dirty = true;
  bool take_screenshot = false;
  double current_learning_rate = 0.0;
  double current_total_error = 0.0;

  void SetTakeScreenshot(int model_idx) {
    WriteMutexLock ml(&video_export_m);
    take_screenshot = true;
    dirty = true;
  }

  void ExportRound(int model_idx, int r) {
    WriteMutexLock ml(&video_export_m);
    if (allow_updates) {
      current_round = r;
      // dirty = true;
    }
  }

  void ExportExamplesPerSec(int model_idx, double eps) {
    WriteMutexLock ml(&video_export_m);
    if (allow_updates) {
      examples_per_second = eps;
      // dirty = true;
    }
  }

  void ExportLearningRate(int model_idx, double rl) {
    WriteMutexLock ml(&video_export_m);
    if (allow_updates) {
      current_learning_rate = rl;
      // dirty = true;
    }
  }

  void ExportNetworkToVideo(int model_idx, const Network &net) {
    WriteMutexLock ml(&video_export_m);
    if (allow_updates) {
      if (current_network == nullptr) {
	current_network = Network::Clone(net);
      } else {
	current_network->CopyFrom(net);
      }
      dirty = true;
    }
  }

  void ExportStimulusToVideo(int model_idx, int example_id, const Stimulation &stim) {
    WriteMutexLock ml(&video_export_m);
    if (allow_updates) {
      CHECK_GE(example_id, 0);
      CHECK_LT(example_id, current_stimulations.size());
      current_stimulations[example_id] = stim;
      dirty = true;
    }
  }

  void ExportExpectedToVideo(int model_idx, int example_id,
			     const vector<float> &expected) {
    WriteMutexLock ml(&video_export_m);
    if (allow_updates) {
      CHECK_GE(example_id, 0);
      CHECK_LT(example_id, current_expected.size());
      current_expected[example_id] = expected;
      dirty = true;
    }
  }

  void ExportErrorsToVideo(int model_idx, int example_id, const Errors &err) {
    WriteMutexLock ml(&video_export_m);
    if (allow_updates) {
      CHECK_GE(example_id, 0);
      CHECK_LT(example_id, current_errors.size());
      current_errors[example_id] = err;
      dirty = true;
    }
  }

  void ExportTotalErrorToVideo(int model_idx, double t) {
    WriteMutexLock ml(&video_export_m);
    if (allow_updates) {
      current_total_error = t;
    }
  }

  // TODO: For this drawing stuff, maybe do it in terms of headless
  // ImageRGBA, and then UI thread just blits that stuff to SDL screen.
  // This lets us export training images for movie, and maybe make the
  // rendering code more portable (e.g. could have web-based UI rather
  // than SDL)?

  // Returns unscaled width, height, and scale (= number of screen
  // pixels per layer pixel)
  static std::tuple<int, int, int>
  PixelSize(const Network &net, int layer) {
    auto MakeScale =
      [](int w, int h) {
	int d = std::max(w, h);
	return (d <= 16) ? 3 : (d <= 128) ? 2 : 1;
      };
    switch (net.renderstyle[layer]) {
    default:
    case RENDERSTYLE_RGB: {
      // Could do a hybrid where we use 3 channels per pixel, but still
      // show them all
      int w = net.width[layer], h = net.height[layer];
      return {w, h, MakeScale(w, h)};
    }
    case RENDERSTYLE_FLAT: {
      int w = net.width[layer] * net.channels[layer];
      int h = net.height[layer];
      return {w, h, MakeScale(w, h)};
    }

    case RENDERSTYLE_INPUTXY:
      return {NOMINAL_CHAR_SIZE, NOMINAL_CHAR_SIZE, 1};
    case RENDERSTYLE_OUTPUTXY:
      return {std::max(NOMINAL_CHAR_SIZE, EmbeddedFont::CHAR_WIDTH * 26),
	      NOMINAL_CHAR_SIZE + 1 + EmbeddedFont::CHAR_HEIGHT,
	      1};
    }
  }

  static int ErrorWidth(const Network &net, int layer) {
    switch (net.renderstyle[layer]) {
    default:
      return net.width[layer];
    }
  }

  // Draw error pixels to the screen.
  template<int SCALE>
  static void Error2D(const vector<float> &errs,
		      // Screen coordinates.
		      int xpos, int ypos,
		      // Unscaled width/height.
		      int width, int height) {
    static_assert(SCALE >= 1);
    for (int y = 0; y < height; y++) {
      const int ystart = ypos + y * SCALE;
      for (int x = 0; x < width; x++) {
	const int idx = y * width + x;
	const int xstart = xpos + x * SCALE;
	// Allow this to not be rectangular, e.g. chessboard.
	uint8 r = ((x + y) & 1) ? 0xFF : 0x00;
	uint8 g = 0x00;
	uint8 b = ((x + y) & 1) ? 0xFF : 0x00;
	if (idx < errs.size()) {
	  const float f = errs[idx];
	  // Use different colors for negative/positive.
	  // XXX: Show errors greater than magnitude 1!
	  if (f < 0.0f) {
	    r = FloatByte(-f);
	    g = 0;
	    b = 0;
	  } else {
	    r = 0;
	    g = FloatByte(f);
	    b = 0;
	  }
	}

	// Hopefully this gets unrolled.
	for (int yy = 0; yy < SCALE; yy++) {
	  for (int xx = 0; xx < SCALE; xx++) {
	    sdlutil::drawpixel(screen, xstart + xx, ystart + yy, r, g, b);
	  }
	}
      }
    }
  }

  // Supposedly SDL prefers this to be called from the main thread.
  void Loop() {
    [[maybe_unused]]
    int mousex = 0, mousey = 0;
    int vlayer = 0;

    for (;;) {
      // int round = ReadWithLock(&video_export_m, &current_round);
      {
	ReadMutexLock ml(&video_export_m);
	if (dirty) {
	  sdlutil::clearsurface(screen, 0x0);
	  const char *paused_msg = allow_updates ? "" : " [^2VIDEO PAUSED]";
	  string menu = StringPrintf(
	      "  round ^3%d ^1|  ^3%0.4f^0 eps    "
	      "^1%.6f^< learning rate   ^1%.6f^< total err%s",
	      current_round,
	      examples_per_second,
	      current_learning_rate,
	      current_total_error,
	      paused_msg);
	  font->draw(2, 2, menu);

	  if (current_network == nullptr) {
	    font->draw(50, 50, "[ ^2No network yet^< ]");
	    SDL_Flip(screen);
	    continue;
	  }

	  int max_width = 0;
	  for (int layer = 0;
	       layer < current_network->num_layers + 1;
	       layer++) {
	    int w = std::get<0>(PixelSize(*current_network, layer));
	    if (DRAW_ERRORS)
	      w += ErrorWidth(*current_network, layer) * 2 + 2;
	    max_width = std::max(w, max_width);
	  }
	  max_width += 4;

	  CHECK(current_stimulations.size() == NUM_VIDEO_STIMULATIONS);
	  CHECK(current_errors.size() == NUM_VIDEO_STIMULATIONS);
	  CHECK(current_expected.size() == NUM_VIDEO_STIMULATIONS);
	  for (int s = 0; s < NUM_VIDEO_STIMULATIONS; s++) {
	    const Stimulation &stim = current_stimulations[s];
	    const Errors &err = current_errors[s];
	    const vector<float> &expected = current_expected[s];
	    // These are not filled in atomically, so it's possible
	    // that we have stimulations but not expected values yet.
	    if (expected.empty()) continue;
	    
	    // Skip if not yet set. These get initialized to empty
	    // sentinels and then updated by the training thread
	    // periodically.
	    if (stim.num_layers == 0 ||
		err.num_layers == 0)
	      continue;

	    CHECK(stim.values.size() == current_network->num_layers + 1);
	    CHECK(stim.values.size() == current_network->renderstyle.size());
	    CHECK(err.error.size() == current_network->num_layers);
	    CHECK(expected.size() == current_network->num_nodes.back());
	    
	    const int xstart = 4 + s * max_width;
	    if (xstart >= SCREENW)
	      break;


	    int ystart = 24;
	    for (int l = 0; l < stim.values.size(); l++) {

	      const auto render_style = current_network->renderstyle[l];
	      switch (render_style) {

	      case RENDERSTYLE_OUTPUTXY: {
		// XXX delete, but save this piece of character prediction code
		const vector<float> &outs = stim.values[l];
		for (int i = 0; i < 26; i++) {
		  const uint8 v = FloatByte(outs[INPUT_LAYER_SIZE + i]);
		  EmbeddedFont::Blit('A' + i,
				     xstart + EmbeddedFont::CHAR_HEIGHT * i,
				     ystart + NOMINAL_CHAR_SIZE + 1,
				     [this, v](int x, int y) {
				       sdlutil::drawclippixel(
					   screen, x, y,
					   v, v, 0xFF);
				     },
				     [this](int x, int y) {
				       sdlutil::drawclippixel(
					   screen, x, y,
					   0, 0, 0);
				     });
		}
	      }

	      case RENDERSTYLE_RGB:
		for (int y = 0; y < current_network->height[l]; y++) {
		  int yy = ystart + y;
		  if (yy >= SCREENH) break;

		  for (int x = 0; x < current_network->width[l]; x++) {

		    int xx = xstart + x;
		    if (x >= SCREENW) break;

		    int cidx = y * current_network->width[l] *
		      current_network->channels[l] +
		      x * current_network->channels[l];

		    // XXX! Need to support more than 3 channels for
		    // chess. Color is dubious anyway?
		    switch (current_network->channels[l]) {
		    case 0: break;
		    case 1: {
		      const uint8 v = FloatByte(stim.values[l][cidx]);
		      sdlutil::drawpixel(screen, xx, yy, v, v, v);
		      break;
		    }
		    case 2: {
		      const uint8 r = FloatByte(stim.values[l][cidx + 0]);
		      const uint8 g = FloatByte(stim.values[l][cidx + 1]);
		      sdlutil::drawpixel(screen, xx, yy, r, g, 0);
		      break;
		    }
		    default:
		      // If more than 3, later ones are just ignored.
		    case 3: {
		      const uint8 r = FloatByte(stim.values[l][cidx + 0]);
		      const uint8 g = FloatByte(stim.values[l][cidx + 1]);
		      const uint8 b = FloatByte(stim.values[l][cidx + 2]);
		      sdlutil::drawpixel(screen, xx, yy, r, g, b);
		      break;
		    }
		    }
		  }
		}
		break;

	      case RENDERSTYLE_FLAT:
		for (int y = 0; y < current_network->height[l]; y++) {
		  const int yy = ystart + y;
		  for (int x = 0;
		       x < current_network->width[l] *
			 current_network->channels[l];
		       x++) {
		    const int xx = xstart + x;
		    if (x >= SCREENW) break;
		    const int cidx =
		      y * current_network->width[l] *
		      current_network->channels[l] + x;
		    const uint8 v = FloatByte(stim.values[l][cidx]);
		    sdlutil::drawpixel(screen, xx, yy, v, v, v);
		  }
		}
		break;
	      }

	      // Now errors.
	      // (XXX Something wrong with this? It draws on top of
	      // the output)
	      if (DRAW_ERRORS && l > 0) {
		const int exstart = xstart +
		  std::get<0>(PixelSize(*current_network, l));
		const vector<float> &errs = err.error[l - 1];
		switch (current_network->renderstyle[l]) {
		default:
		case RENDERSTYLE_FLAT:
		case RENDERSTYLE_RGB:
		  Error2D<2>(errs, exstart, ystart,
			     current_network->width[l] *
			     current_network->channels[l],
			     current_network->height[l]);
		}
	      }

	      ystart += std::get<1>(PixelSize(*current_network, l)) + 4;
	    }


	    if (vlayer >= 0) {
	      double tot = 0.0;
	      int yz = ystart + 4;
	      const vector<float> &vals = stim.values[vlayer];
	      for (int i = 0; i < vals.size(); i++) {
		tot += vals[i];
		if (i < 32) {
		  // Font color.
		  const int c = vals[i] > 0.5f ? 0 : 1;
		  if (vlayer > 0) {
		    const float e = vlayer > 0 ? err.error[vlayer - 1][i] : 0.0;
		    // Font color.
		    const int sign = (e == 0.0) ? 4 : (e < 0.0f) ? 2 : 5;
		    font->draw(xstart, yz,
			       StringPrintf("^%d%.9f ^%d%.12f",
					    c, vals[i], sign, e));
		  } else {
		    font->draw(xstart, yz,
			       StringPrintf("^%d%.9f", c, vals[i]));
		  }
		  yz += FONTHEIGHT;
		}
	      }
	      font->draw(xstart, yz,
			 StringPrintf("[%d] tot: %.9f", vlayer, tot));
	    }
	  }

	  SDL_Flip(screen);

	  if (take_screenshot) {
	    string scfile = StringPrintf("video/train%d.png", current_round);
	    sdlutil::SavePNG(scfile, screen);
	    Printf("Wrote screenshot to %s\n", scfile.c_str());
	    take_screenshot = false;
	  }

	  dirty = false;
	}
      }

      if (ReadWithLock(&train_done_m, &train_done)) {
	Printf("UI thread saw that training finished.\n");
	return;
      }

      SDL_Event event;
      if (SDL_PollEvent(&event)) {
	if (event.type == SDL_QUIT) {
	  Printf("QUIT.\n");
	  return;

	} else if (event.type == SDL_MOUSEMOTION) {
	  SDL_MouseMotionEvent *e = (SDL_MouseMotionEvent*)&event;

	  mousex = e->x;
	  mousey = e->y;
	  // (and immediately redraw)

	} else if (event.type == SDL_KEYDOWN) {
	  switch (event.key.keysym.sym) {
	  case SDLK_ESCAPE:
	    Printf("ESCAPE.\n");
	    return;
	  case SDLK_SPACE:
	    {
	      WriteMutexLock ml(&video_export_m);
	      allow_updates = !allow_updates;
	    }
	    break;

	  case SDLK_v:
	    {
	      ReadMutexLock ml(&video_export_m);
	      // Allow vlayer to go to -1 (off), but wrap around
	      // below that.

	      if (current_network != nullptr) {
		vlayer++;
		if (vlayer < -1) {
		  // There are actually num_layers + 1 stimulations.
		  vlayer = current_network->num_layers;
		} else if (vlayer > current_network->num_layers) {
		  vlayer = -1;
		}
		dirty = true;
	      }
	    }
	    break;
	  default:;
	  }
	}
      } else {
	SDL_Delay(1000);
      }
    }
  }

};


static std::unique_ptr<Network> CreateInitialNetwork(ArcFour *rc) {

  const int num_layers = 4;
  const vector<int> width_config =
    { INPUT_LAYER_SIZE,
      INPUT_LAYER_SIZE * 2,
      INPUT_LAYER_SIZE * 2,
      OUTPUT_LAYER_SIZE * 2,
      OUTPUT_LAYER_SIZE, };

  // If zero, automatically factor to make square-ish.
  const vector<int> height_config =
    { 0, 0, 0, 0, 0, };

  // input layer could maybe be thought of as two-channel (x,y),
  // but output has extra stuff (predicted character). Instead
  // we just think of these as flat and do custom rendering and
  // index assignment.
  const vector<int> channels = { 1, 1, 1, 1, 1, };

  // When zero, create a dense layer.
  const vector<int> indices_per_node_config =
    {
      // Start dense
      0,
      0,
      // Then half of the points
      INPUT_LAYER_SIZE,
      OUTPUT_LAYER_SIZE,
    };

  const vector<uint32_t> renderstyle = {
    RENDERSTYLE_INPUTXY,
    RENDERSTYLE_FLAT,
    RENDERSTYLE_FLAT,
    RENDERSTYLE_FLAT,
    RENDERSTYLE_OUTPUTXY,
  };

  const vector<TransferFunction> transfer_functions = {
    LEAKY_RELU,
    LEAKY_RELU,
    LEAKY_RELU,
    LEAKY_RELU,
  };

  vector<int> height, width;
  CHECK(height_config.size() == width_config.size());
  for (int i = 0; i < height_config.size(); i++) {
    int w = width_config[i];
    int h = height_config[i];
    CHECK(w > 0);
    CHECK(h >= 0);
    if (h == 0) {
      if (w == 1) {
	width.push_back(1);
	height.push_back(1);
      } else {
	// Try to make a rectangle that's squareish.
	vector<int> factors = Util::Factorize(w);

	CHECK(!factors.empty()) << w << " has no factors??";

	// XXX Does this greedy approach produce good results?
	int ww = factors.back(), hh = 1;
	factors.pop_back();

	for (int f : factors) {
	  if (ww < hh)
	    ww *= f;
	  else
	    hh *= f;
	}

	CHECK(ww * hh == w);
	width.push_back(ww);
	height.push_back(hh);
      }
    } else {
      width.push_back(w);
      height.push_back(h);
    }
  }

  // num_nodes = width * height * channels
  // //  indices_per_node = indices_per_channel * channels
  // //  vector<int> indices_per_node;
  vector<int> num_nodes;
  CHECK_EQ(width.size(), height.size());
  CHECK_EQ(width.size(), channels.size());
  CHECK_EQ(width.size(), renderstyle.size());
  CHECK_EQ(num_layers + 1, height.size());
  CHECK_EQ(num_layers, indices_per_node_config.size());
  CHECK_EQ(num_layers, transfer_functions.size());
  for (int i = 0; i < num_layers + 1; i++) {
    CHECK(width[i] >= 1);
    CHECK(height[i] >= 1);
    CHECK(channels[i] >= 1);
    num_nodes.push_back(width[i] * height[i] * channels[i]);
  }

  vector<int> indices_per_node;
  vector<LayerType> layer_type;
  for (int i = 0; i < indices_per_node_config.size(); i++) {
    int ipc = indices_per_node_config[i];
    if (ipc == 0) {
      layer_type.push_back(LAYER_DENSE);
      indices_per_node.push_back(width[i] * height[i] * channels[i]);
    } else {
      CHECK(ipc > 0);
      CHECK(ipc <= width[i] * height[i] * channels[i]);
      layer_type.push_back(LAYER_SPARSE);
      indices_per_node.push_back(ipc);
    }
  }

  std::unique_ptr<Network> net =
    std::make_unique<Network>(
	num_nodes, indices_per_node, transfer_functions);
  net->width = width;
  net->height = height;
  net->channels = channels;
  net->renderstyle = renderstyle;

  // XXX should probably be ctor arg?
  CHECK(net->layers.size() == layer_type.size());
  for (int i = 0; i < net->layers.size(); i++) {
    net->layers[i].type = layer_type[i];
  }

  Printf("Randomize weights:\n");
  RandomizeNetwork(rc, net.get());
  Printf("Gen indices:\n");
  MakeIndices(rc, net.get());
  Printf("Invert indices:\n");
  Network::ComputeInvertedIndices(net.get());

  // Should be well-formed now.
  net->StructuralCheck();
  net->NaNCheck("initial network");

  return net;
}

static UI *ui = nullptr;
static VectorLoadFonts *load_fonts = nullptr;

struct TrainingExample {
  vector<float> input;
  vector<float> output;
};

// Encapsulates training of a single network. Allows for multiple
// simultaneous instances, allowing networks to be cotrained.
struct Training {

  // Number of examples per round of training. Includes eval
  // examples.
  static constexpr int EXAMPLES_PER_ROUND = 1024;
  // Number of examples that are eval inputs (not trained); the
  // remainder are training examples.
  static constexpr int EVAL_INPUTS_PER_ROUND = EXAMPLES_PER_ROUND / 4;
  static constexpr int TRAINING_PER_ROUND = EXAMPLES_PER_ROUND - EVAL_INPUTS_PER_ROUND;
  static_assert(TRAINING_PER_ROUND > 0);
  static_assert(TRAINING_PER_ROUND <= EXAMPLES_PER_ROUND);
  static_assert(EVAL_INPUTS_PER_ROUND < EXAMPLES_PER_ROUND);
  static_assert(EVAL_INPUTS_PER_ROUND >= 0);
  static_assert(EXAMPLES_PER_ROUND > 0);
  
  // Write a screenshot of the UI (to show training progress for
  // videos, etc.) every time the network does this many rounds.
  static constexpr int SCREENSHOT_ROUND_EVERY = 50;

  // On a verbose round, we write a network checkpoint and maybe some
  // other stuff to disk. XXX: Do this based on time, since round
  // speed can vary a lot based on other parameters!
  static constexpr int VERBOSE_ROUND_EVERY = 250;
  
  
  explicit Training(int model_index) :
    model_index(model_index),
    model_filename(StringPrintf("net%d.val", model_index)),
    rc(GetRandomSeed(model_index)) {
    Timer setup_timer;

    CHECK(ui != nullptr) << "Must be created first.";
    CHECK(load_fonts != nullptr) << "Must be created first.";

    rc.Discard(2000);

    // Random state that can be accessed in parallel for each example.
    // Note, 256 * 4096 is kinda big (one megabyte). We could get away
    // with capping this at the number of threads, if threadutil (really
    // autoparallel) had a way of passing thread local data.
    printf("Initialize example-local random streams...\n");
    vector<ArcFour> example_rc;
    for (int i = 0; i < EXAMPLES_PER_ROUND; i++) {
      std::vector<uint8> init;
      init.reserve(64);
      for (int j = 0; j < 64; j++) init.push_back(rc.Byte());
      example_rc.emplace_back(init);
    }

    // Load the existing network from disk or create the initial one.
    Timer initialize_network_timer;
    // Try loading from disk; null on failure.
    net.reset(Network::ReadNetworkBinary(model_filename));

    if (net.get() == nullptr) {
      Printf("Initializing new network...\n");
      net = CreateInitialNetwork(&rc);
      CHECK(net.get() != nullptr);

      Printf("Writing network so we don't have to do that again...\n");
      Network::SaveNetworkBinary(*net, model_filename);
    }

    Printf("Initialized network in %.1fms.\n", initialize_network_timer.MS());

    // Create kernels right away so that we get any compilation errors early.
    forwardlayer = make_unique<ForwardLayerCL>(global_cl, *net);
    setoutputerror = make_unique<SetOutputErrorCL>(global_cl);
    backwardlayer = make_unique<BackwardLayerCL>(global_cl, *net);
    updateweights = make_unique<UpdateWeightsCL>(global_cl, *net);
    
    net_gpu = make_unique<NetworkGPU>(global_cl, net.get());

    Printf("Network uses %.2fMB of storage (without overhead).\n",
	   net->Bytes() / (1024.0 * 1024.0));
    {
      Stimulation tmp(*net);
      int64 stim_bytes = tmp.Bytes();
      Printf("A stimulation is %.2fMB, so for %d examples we need %.2fMB\n",
	     stim_bytes / (1024.0 * 1024.0), EXAMPLES_PER_ROUND,
	     (stim_bytes * EXAMPLES_PER_ROUND) / (1024.0 * 1024.0));
    }

    for (int i = 0; i < EXAMPLES_PER_ROUND; i++) {
      training.push_back(new TrainingRoundGPU{global_cl, *net});
    }

    // Automatically tune parallelism for some loops, caching the results
    // on disk. The experiment string should change (or cache files deleted)
    // when significant parameters change (not just these!).
    // Ideally these should be shared between the two models since they should
    // have the same performance, but it's not that wasteful to duplicate
    // the sampling and maybe would be worse to have lock contention.
    const string experiment =
      StringPrintf("sdf-%d.%d", EXAMPLES_PER_ROUND, model_index);

    stim_init_comp = std::make_unique<AutoParallelComp>(32, 50, false,
							StringPrintf("autoparallel.%s.stim.txt",
								     experiment.c_str()));

    forward_comp = std::make_unique<AutoParallelComp>(32, 50, false,
						      StringPrintf("autoparallel.%s.fwd.txt",
								   experiment.c_str()));
    
    error_comp = std::make_unique<AutoParallelComp>(32, 50, false,
						    StringPrintf("autoparallel.%s.err.txt",
								 experiment.c_str()));
    
    backward_comp = std::make_unique<AutoParallelComp>(32, 50, false,
						       StringPrintf("autoparallel.%s.bwd.txt",
								    experiment.c_str()));
    
  }

  std::unique_ptr<AutoParallelComp> stim_init_comp, forward_comp,
    error_comp, backward_comp;
  
  // Run one training round. Might stall if starved for examples.
  // Will exit early if global train_should_die becomes true.
  void RunRound() {

    // XXX members?
    double setup_ms = 0.0, stimulation_init_ms = 0.0, forward_ms = 0.0,
      fc_init_ms = 0.0, bc_init_ms = 0.0, kernel_ms = 0.0, backward_ms = 0.0,
      output_error_ms = 0.0, update_ms = 0.0,
      eval_ms = 0.0;

    Timer round_timer;

    if (VERBOSE > 2) Printf("\n\n");
    Printf("[%d] ** NET ROUND %d (%d in this process) **\n",
	   model_index, net->rounds, rounds_executed);

    // The learning rate should maybe depend on the number of examples
    // per round, since we integrate over lots of them. We could end
    // up having a total error for a single node of like
    // +EXAMPLES_PER_ROUND or -EXAMPLES_PER_ROUND, which could yield
    // an unrecoverable-sized update. We now divide the round learning rate
    // to an example learning rate, below. UpdateWeights also caps the
    // maximum increment to +/- 1.0f, which is not particularly principled
    // but does seem to robustly prevent runaway.

    constexpr int TARGET_ROUNDS = 50000;
    auto Linear =
      [](double start, double end, double round_target, double input) {
	if (input < 0.0) return start;
	if (input > round_target) return end;
	double height = end - start;
	double f = input / round_target;
	return start + f * height;
      };
    const float round_learning_rate =
      Linear(0.10, 0.002, TARGET_ROUNDS, net->rounds);

    CHECK(!std::isnan(round_learning_rate));
    if (VERBOSE > 2) Printf("Learning rate: %.4f\n", round_learning_rate);

    const float example_learning_rate =
      round_learning_rate / (double)TRAINING_PER_ROUND;

    const bool is_verbose_round =
      0 == ((rounds_executed /* + 1 */) % VERBOSE_ROUND_EVERY);
    if (is_verbose_round) {
      Printf("Writing network:\n");
      net_gpu->ReadFromGPU();
      Network::SaveNetworkBinary(*net,
				 StringPrintf("network-%d-checkpoint.bin",
					      model_index));
    }

    if (ShouldDie()) return;
    
    if (VERBOSE > 2) Printf("Export network:\n");
    ui->ExportRound(model_index, net->rounds);
    ui->ExportLearningRate(model_index, round_learning_rate);

    // XXX do this in video?
    const bool take_screenshot =
      net->rounds % SCREENSHOT_ROUND_EVERY == 0;
    if (rounds_executed % EXPORT_EVERY == 0 ||
	take_screenshot) {
      net_gpu->ReadFromGPU();
      ui->ExportNetworkToVideo(model_index, *net);
    }

    if (take_screenshot) {
      // Training screenshot.
      ui->SetTakeScreenshot(model_index);

      // Stable eval screenshot on Helvetica (iterative).
      FontProblem::RenderVector("helvetica.ttf",
				*net,
				ROW_MAX_POINTS,
				StringPrintf("eval%d/eval%lld.png",
					     model_index,
					     net->rounds));
    }

    if (CHECK_NANS) {
      net_gpu->ReadFromGPU();
      net->NaNCheck(StringPrintf("round start model %d", model_index));
    }

    Timer setup_timer;
    if (VERBOSE > 2) Printf("Setting up batch:\n");

    // examples include training examples and "eval" examples
    vector<TrainingExample> examples;
    examples.reserve(EXAMPLES_PER_ROUND);

    auto GetExamples = [&examples](std::shared_mutex *mut,
				   deque<TrainingExample> *queue,
				   int target_num,
				   const char *type) {
	for (;;) {
	  {
	    WriteMutexLock ml{mut};
	    while (examples.size() < target_num &&
		   !queue->empty()) {
	      examples.push_back(std::move(queue->front()));
	      queue->pop_front();
	    }
	  }

	  if (examples.size() >= target_num)
	    return;
	  
	  if (VERBOSE > 0)
	    Printf("Blocked grabbing %s examples (still need %d)...\n",
		   type,
		   target_num - examples.size());
	  std::this_thread::sleep_for(100ms);
	}
      };

    GetExamples(&example_queue_m, &example_queue,
		TRAINING_PER_ROUND, "training");
    GetExamples(&eval_queue_m, &eval_queue,
		EXAMPLES_PER_ROUND, "eval");

    CHECK(examples.size() == EXAMPLES_PER_ROUND);
    
    setup_ms += setup_timer.MS();

    if (false && CHECK_NANS /* && net->rounds == 3 */) {
      // Actually run the forward pass on CPU, trying to find the
      // computation that results in nan...
      net_gpu->ReadFromGPU(); // should already be here, but...
      CHECK(!examples.empty());
      const TrainingExample &example = examples[0];
      Stimulation stim{*net};
      CHECK_EQ(stim.values[0].size(), example.input.size());
      for (int i = 0; i < example.input.size(); i++)
	stim.values[0][i] = example.input[i];
      net->RunForwardVerbose(&stim);
    }

    // TODO: may make sense to pipeline this loop somehow, so that we
    // can parallelize CPU/GPU duties?

    // Run a batch of images all the way through. (Each layer requires
    // significant setup.)
    Timer stimulation_init_timer;

    if (VERBOSE > 2) Printf("Setting input layer of Stimulations...\n");
    // These are just memory copies; easy to do in parallel.
    CHECK_EQ(examples.size(), training.size());
    stim_init_comp->ParallelComp(
	examples.size(),
	[this, &examples](int i) {
	  training[i]->LoadInput(examples[i].input);
	});
    stimulation_init_ms += stimulation_init_timer.MS();

    if (ShouldDie()) return;

    // We run the forward pass for both training and eval examples.
    // The loop over layers must be in serial.
    for (int src = 0; src < net->num_layers; src++) {
      if (VERBOSE > 2) Printf("FWD Layer %d: ", src);
      Timer fc_init_timer;
      ForwardLayerCL::ForwardContext fc(forwardlayer.get(), net_gpu.get(), src);
      fc_init_ms += fc_init_timer.MS();

      // Can be parallel, but watch out about loading the GPU with
      // too many simultaneous value src/dst buffers.
      Timer forward_timer;
      if (VERBOSE > 3) Printf("Parallelcomp...\n");
      forward_comp->ParallelComp(
	  examples.size(),
	  [this, num_examples = examples.size(), &fc](int example_idx) {
	    fc.Forward(training[example_idx]);

	    if (rounds_executed % EXPORT_EVERY == 0 &&
		example_idx < NUM_VIDEO_STIMULATIONS) {
	      // XXX this uses unintialized/stale memory btw
	      Stimulation stim{*net};
	      training[example_idx]->ExportStimulation(&stim);
	      // Copy to screen.
	      ui->ExportStimulusToVideo(model_index, example_idx, stim);
	    }
	  });
      forward_ms += forward_timer.MS();
      kernel_ms += fc.kernel_ms;
      if (VERBOSE > 2) Printf("\n");
    }

    if (CHECK_NANS) {
      for (int example_idx = 0; example_idx < training.size(); example_idx++) {
	Stimulation stim{*net};
	training[example_idx]->ExportStimulation(&stim);
	stim.NaNCheck(StringPrintf("forward pass model %d", model_index));
      }
    }

    // TODO HERE: Publish the eval results (or return them from the function call?)

    // Compute total error (average per training example).
    if (rounds_executed % EXPORT_EVERY == 0) {
      CHECK_EQ(examples.size(), training.size());
      // We don't use the stimulus, because we want the total over all
      // examples (but we only export enough for the video above), and
      // only need the final output values, not internal layers.
      double total_error = 0.0;
      vector<float> values;
      values.resize(net->num_nodes[net->num_layers]);
      CHECK(TRAINING_PER_ROUND <= examples.size());
      for (int i = 0; i < TRAINING_PER_ROUND; i++) {
	training[i]->ExportOutput(&values);
	CHECK(examples[i].output.size() == values.size());
	for (int j = 0; j < values.size(); j++) {
	  float d = values[j] - examples[i].output[j];
	  total_error += fabs(d);
	}
      }
      ui->ExportTotalErrorToVideo(model_index, total_error / (double)TRAINING_PER_ROUND);
    }

    if (ShouldDie()) return;

    // For error computation (etc.) we only use the training examples.
    // Compute expected.
    if (VERBOSE > 2) Printf("Error calc.\n");
    CHECK(TRAINING_PER_ROUND <= examples.size());
    Timer output_error_timer;
    error_comp->ParallelComp(
	TRAINING_PER_ROUND,
	[this, &examples](int example_idx) {
	  // (Do this after finalizing the expected vector above.)
	  if (rounds_executed % EXPORT_EVERY == 0 &&
	      example_idx < NUM_VIDEO_STIMULATIONS) {
	    // Copy to screen.
	    ui->ExportExpectedToVideo(model_index,
				      example_idx,
				      examples[example_idx].output);
	  }
	  
	  // PERF we could have loaded this a lot earlier
	  training[example_idx]->LoadExpected(examples[example_idx].output);
	  SetOutputErrorCL::Context sc{setoutputerror.get(), net_gpu.get()};
	  sc.SetOutputError(training[example_idx]);
	});
    output_error_ms += output_error_timer.MS();
    if (VERBOSE > 2) Printf("\n");

    if (ShouldDie()) return;
    if (VERBOSE > 2) Printf("Backwards:\n");
    // Also serial, but in reverse.
    Timer backward_timer;
    // We do NOT propagate errors to the input layer, so dst is
    // strictly greater than 0.
    for (int dst = net->num_layers - 1; dst > 0; dst--) {
      if (VERBOSE > 2) Printf("BWD Layer %d: ", dst);

      Timer bc_init_timer;
      BackwardLayerCL::Context bc{backwardlayer.get(), net_gpu.get(), dst};
      bc_init_ms += bc_init_timer.MS();

      backward_comp->ParallelComp(
	  TRAINING_PER_ROUND,
	  [this, &bc](int example) {
	    bc.Backward(training[example]);
	  });
      if (VERBOSE > 2) Printf("\n");
    }
    backward_ms += backward_timer.MS();

    if (rounds_executed % EXPORT_EVERY == 0) {
      for (int example_idx = 0;
	   example_idx < NUM_VIDEO_STIMULATIONS;
	   example_idx++) {
	Errors err{*net};
	training[example_idx]->ExportErrors(&err);
	ui->ExportErrorsToVideo(model_index, example_idx, err);
      }
    }

    if (ShouldDie()) return;
    if (VERBOSE > 2) Printf("Update weights:\n");
    Timer update_timer;

    // Don't parallelize! These are all writing to the same network
    // weights. Each call is parallelized, though.
    for (int layer = 0; layer < net->num_layers; layer++) {
      UpdateWeightsCL::Context uc{updateweights.get(), net_gpu.get(), layer};

      // PERF Faster to try to run these in parallel (maybe
      // parallelizing memory traffic with kernel execution -- but we
      // can't run the kernels at the same time).
      for (int example = 0; example < TRAINING_PER_ROUND; example++) {
	uc.Update(example_learning_rate, training[example], layer);
      }

      // Now we leave the network on the GPU, and the version in the
      // Network object will be out of date. But flush the command
      // queue. (why? I guess make sure that we're totally done
      // writing since other parts of the code assume concurrent reads
      // are ok?)
      uc.Finish();
      /*
	Printf("[%d/%d] = (%.2f%%) ", layer, net->num_layers,
	layer * 100.0 / net->num_layers);
      */
    }
    update_ms += update_timer.MS();
    if (VERBOSE > 2) Printf("\n");

    if (CHECK_NANS) {
      net_gpu->ReadFromGPU();
      net->NaNCheck(StringPrintf("updated weights model %d", model_index));
    }

    if (ShouldDie()) return;

    net->rounds++;
    net->examples += TRAINING_PER_ROUND;

    double round_ms = round_timer.MS();
    auto Pct = [round_ms](double d) { return (100.0 * d) / round_ms; };
    double denom = rounds_executed + 1;

    // TODO: Would be nice to average this over the last few rounds.
    const double round_eps = EXAMPLES_PER_ROUND / (round_ms / 1000.0);
    // (EXAMPLES_PER_ROUND * denom) / (total_ms / 1000.0)
    ui->ExportExamplesPerSec(model_index, round_eps);
    if (VERBOSE > 1)
      Printf("Round time: %.1fs. (%.1f eps)\n"
	     "We spent %.1fms in setup (%.1f%%),\n"
	     "%.1fms in stimulation init (%.1f%%),\n"
	     "%.1fms in eval (main thread; amortized) (%.1f%%),\n"
	     "%.1fms in forward layer (%.1f%%),\n"
	     "%.1fms in fc init (%.1f%%),\n"
	     "%.1fms in forward layer kernel (at most; %.1f%%).\n"
	     "%.1fms in error for output layer (%.1f%%),\n"
	     "%.1fms in bc init (%.1f%%),\n"
	     "%.1fms in backwards pass (%.1f%%),\n"
	     "%.1fms in updating weights (%.1f%%),\n",
	     round_ms / 1000.0, round_eps,
	     setup_ms / denom, Pct(setup_ms),
	     stimulation_init_ms / denom, Pct(stimulation_init_ms),
	     eval_ms / denom, Pct(eval_ms),
	     forward_ms / denom, Pct(forward_ms),
	     fc_init_ms / denom, Pct(fc_init_ms),
	     kernel_ms / denom, Pct(kernel_ms),
	     output_error_ms / denom, Pct(output_error_ms),
	     bc_init_ms / denom, Pct(bc_init_ms),
	     backward_ms / denom, Pct(backward_ms),
	     update_ms / denom, Pct(update_ms));
    
  }
  
  // Generate training examples in another thread and feed them to AddExampleToQueue.
  // To avoid overloading, throttle while WantsExamples is false.
  bool WantsExamples() {
    ReadMutexLock ml(&example_queue_m);
    return example_queue.size() < EXAMPLE_QUEUE_TARGET;
  }
  void AddExampleToQueue(TrainingExample example) {
    WriteMutexLock ml(&example_queue_m);
    example_queue.emplace_back(std::move(example));
    // PERF start moving it to GPU?
  }

  // Like the above, but for "eval" inputs, which do not have expected
  // outputs and are not used for training. This is particular to the
  // SDF font problem where we generate inverse examples for cotraining.
  bool WantsEvalInputs() {
    ReadMutexLock ml(&eval_queue_m);
    return eval_queue.size() < EVAL_QUEUE_TARGET;
  }
  // XXX use trainingexample or just vector<float>?
  void AddEvalInputToQueue(TrainingExample example) {
    WriteMutexLock ml(&eval_queue_m);
    eval_queue.emplace_back(std::move(example));
  }

  ~Training() {
    for (TrainingRoundGPU *trg : training)
      delete trg;
    training.clear();
  }
  
private:
  // Try to keep twice that in the queue all the time.
  static constexpr int EXAMPLE_QUEUE_TARGET =
    std::max(TRAINING_PER_ROUND * 3, 256);
  static constexpr int EVAL_QUEUE_TARGET =
    EVAL_INPUTS_PER_ROUND * 3;
  
  static string GetRandomSeed(int idx) {
    const string start_seed = StringPrintf("%d  %lld  %d",
					   getpid(),
					   (int64)time(nullptr),
					   idx);
    Printf("Start seed: [%s]\n", start_seed.c_str());
    return start_seed;
  }

  bool ShouldDie() {
    CHECK(net.get()) << "Only call this after the network has been loaded";
    const bool should_die = ReadWithLock(&train_should_die_m, &train_should_die);
    if (should_die) {
      Printf("Train thread signaled death.\n");
      Printf("Saving to %s...\n", model_filename.c_str());
      Network::SaveNetworkBinary(*net, model_filename);
    }
    return should_die;
  }

  
  const int model_index = 0;
  const string model_filename;
  ArcFour rc;
  vector<ArcFour> example_rc;
  std::unique_ptr<Network> net;
  int64 rounds_executed = 0;
  
  // Separate kernels per training instance, since they can be specialized
  // to the network's parameters.
  std::unique_ptr<ForwardLayerCL> forwardlayer;
  std::unique_ptr<SetOutputErrorCL> setoutputerror;
  std::unique_ptr<BackwardLayerCL> backwardlayer;
  std::unique_ptr<UpdateWeightsCL> updateweights;

  // The network's presence on the GPU. It can be out of date with the
  // CPU copy in the net variable.
  std::unique_ptr<NetworkGPU> net_gpu;

  // We use the same structures to hold all the stimulations and errors
  // on the GPU. Size EXAMPLES_PER_ROUND.
  vector<TrainingRoundGPU *> training;

  // Protects example_queue.
  std::shared_mutex example_queue_m;
  // (Note: This could just be a vector and use stack ordering, but
  // we want to make sure that we get to every training example for cases
  // that they are not randomly sampled (e.g. co-training).)
  // (Actually since these are separate, it probably could just be a vector now.)
  deque<TrainingExample> example_queue;

  // Same idea but eval examples, for which we don't have an "expected output".
  // Used for co-training.
  std::shared_mutex eval_queue_m;
  deque<TrainingExample> eval_queue;
};

// Inserts examples into the Training objects, and reads their eval outputs
// to generate examples for the co-trained pair.
// FIXME this function is fuct
static void MakeTrainingExamplesThread(
    Training *make_lowercase,
    Training *make_uppercase) {

  Printf("Training example thread startup.\n");
  string seed = StringPrintf("make ex %lld", (int64)time(nullptr));
  ArcFour rc(seed);
  rc.Discard(2000);

  #if 0
  // Returns true if successful.
  // XXXX olde
  auto PopulateExampleFromFont =
    [](ArcFour *rc, const TTF *ttf, TrainingExample *example) -> bool {

      const int char_offset = RandTo(rc, 26);
      const int upper_char = 'A' + char_offset;
      const int lower_char = 'a' + char_offset;

      example->input.resize(INPUT_LAYER_SIZE);
      if (!FontProblem::FillVector(ttf, upper_char, ROW_MAX_POINTS,
				   example->input.data())) {
	return false;
      }

      example->output.resize(OUTPUT_LAYER_SIZE, 0.0f);
      if (!FontProblem::GetRows(ttf, lower_char, ROW_MAX_POINTS,
				&example->expected_contours)) {
	return false;
      }

      #if 0
      // filled later.
      if (!FontProblem::FillVector(ttf, lower_char, ROW_MAX_POINTS,
				   example->output.data())) {
	return false;
      }
      #endif

      // Output is same as input size, plus 26 letters. So use the
      // length of the input to know where to start putting letters.
      int next_idx = example->input.size();
      for (int i = 0; i < 26; i++) {
	CHECK(next_idx < example->output.size());
	example->output[next_idx] = i == char_offset ? 1.0f : 0.0f;
	next_idx++;
      }
      CHECK(next_idx == example->output.size());

      if (CHECK_NANS) {
	for (float f : example->input) { CHECK(!std::isnan(f)); }
	for (float f : example->output) { CHECK(!std::isnan(f)); }
      }

      return true;
    };

  auto MakeTrainingExamplesThread = [&training_examples_m,
				     &training_examples,
				     &PopulateExampleFromFont]() {

    for (;;) {
      if (ReadWithLock(&train_should_die_m, &train_should_die)) {
	return;
      }

      training_examples_m.lock_shared();
      // Make sure we have plenty of examples so that learning doesn't stall.
      if (training_examples.size() < EXAMPLE_QUEUE_TARGET) {
	training_examples_m.unlock_shared();

	TrainingExample example;
	{
	  load_fonts->fonts_m.lock_shared();

	  if (load_fonts->fonts.size() < ENOUGH_FONTS) {
	    load_fonts->fonts_m.unlock_shared();
	    Printf("Not enough training data loaded yet!\n");
	    std::this_thread::sleep_for(1s);
	    continue;
	  } else {
	    const int idx = RandTo(&rc, load_fonts->fonts.size());

	    bool ok = PopulateExampleFromFont(&rc,
					      load_fonts->fonts[idx],
					      &example);
	    load_fonts->fonts_m.unlock_shared();

	    if (ok) {
	      WriteMutexLock ml(&training_examples_m);
	      training_examples.push_back(std::move(example));
	    }

	    // PERF: If not ok, we could mark this as a bad training
	    // example so we don't keep trying (or just remove it from
	    // the vector?) ... but these should have been pre-filtered.
	  }
	}

      } else {
	training_examples_m.unlock_shared();
	std::this_thread::sleep_for(10ms);
      }
    }
    Printf("Training example generator exiting.\n");
  };
  #endif
  
}

// XXX it ded
#if 0
static void TrainThread() {

  if (ReadWithLock(&train_should_die_m, &train_should_die))
    return;


  if (ShouldDie()) return;

  // (Note that the below actually works while the examples are still
  // being loaded, but maybe overkill since this dataset should load
  // pretty fast? So syncing here.)
  Printf("Waiting for fonts...\n");
  CHECK(load_fonts != nullptr);
  load_fonts->Sync();
  Printf("Fonts all loaded.\n");




  if (ShouldDie()) return;



  Printf(" ** Done. **");

  WriteWithLock(&train_done_m, &train_done, true);
}
#endif

// Periodically we check to see if any process name matches something
// in this function. If so, we pause training.
static std::optional<string> GetExclusiveApp() {
  vector<string> procs = Top::Enumerate();

  for (const string &proc : procs) {
    string match = Util::lcase(proc);
    // Now you can see what games I'm playing in December 2020!
    if (match == "spel2.exe") return {proc};
    if (match == "disc room.exe") return {proc};
    if (match == "superliminalsteam.exe") return {proc};
    // Can add more here, including regexes etc...
  }

  return nullopt;
}

int SDL_main(int argc, char **argv) {
  // XXX This is specific to my machine. You probably want to remove it.
  // Assumes that processors 0-16 are available.
  // CHECK(SetProcessAffinityMask(GetCurrentProcess(), 0xF));

  if (!SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS)) {
    LOG(FATAL) << "Unable to go to BELOW_NORMAL priority.\n";
  }

  /* Initialize SDL and network, if we're using it. */
  CHECK(SDL_Init(SDL_INIT_VIDEO |
		 SDL_INIT_TIMER |
		 SDL_INIT_AUDIO) >= 0);
  fprintf(stderr, "SDL initialized OK.\n");

  SDL_EnableKeyRepeat(SDL_DEFAULT_REPEAT_DELAY,
                      SDL_DEFAULT_REPEAT_INTERVAL);

  SDL_EnableUNICODE(1);

  SDL_Surface *icon = SDL_LoadBMP("lowercase-icon.bmp");
  if (icon != nullptr) {
    SDL_WM_SetIcon(icon, nullptr);
  }

  screen = sdlutil::makescreen(SCREENW, SCREENH);
  CHECK(screen);

  font = Font::create(screen,
		      "font.png",
		      FONTCHARS,
		      FONTWIDTH, FONTHEIGHT, FONTSTYLES, 1, 3);
  CHECK(font != nullptr) << "Couldn't load font.";

  global_cl = new CL;

  // Start loading fonts in background.
  load_fonts = new VectorLoadFonts(
      []() { return ReadWithLock(&train_should_die_m, &train_should_die); },
      {ROW0_MAX_PTS, ROW1_MAX_PTS, ROW2_MAX_PTS},
      12,
      MAX_FONTS);

  ui = new UI;

  vector<Training *> training = {nullptr, nullptr};
  ParallelComp(2, [&training](int idx) { training[idx] = new Training(idx); }, 2);
  // Aliases for convenience.
  Training *make_lowercase = training[0];
  Training *make_uppercase = training[1];
  
  // Start generating examples in another thread and feeding them to
  // the two training instances.
  std::thread examples_thread{
    [make_lowercase, make_uppercase]() {
      MakeTrainingExamplesThread(make_lowercase, make_uppercase);
    }};

  // FIXME call RunRound, but throttle if necessary
  std::thread train_thread{
    [make_lowercase, make_uppercase]() {

      while (!ReadWithLock(&train_should_die_m, &train_should_die)) {

	// Pause if exclusive app is running.
	while (std::optional<string> excl = GetExclusiveApp()) {
	  Printf("(Sleeping because of exclusive app %s)\n",
		 excl.value().c_str());
	  std::this_thread::sleep_for(5000ms);
	}

	// XXX in parallel!
	make_lowercase->RunRound();
	make_uppercase->RunRound();	
	
      }
    }};
  
  ui->Loop();

  Printf("Killing train thread (might need to wait for round to finish)...\n");
  WriteWithLock(&train_should_die_m, &train_should_die, true);
  for (Training *t : training) delete t;
  examples_thread.join();
  
  Printf("Train is dead; now UI exiting.\n");

  SDL_Quit();
  return 0;
}
