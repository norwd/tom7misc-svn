
#include "network-gpu.h"

#include <string>
#include <optional>
#include <vector>
#include <mutex>
#include <cmath>

#include "base/logging.h"
#include "base/stringprintf.h"

#include "util.h"
#include "network.h"
#include "timer.h"
#include "clutil.h"
#include "threadutil.h"

using namespace std;


NetworkGPU::NetworkGPU(CL *cl, Network *net) : cl(cl), net(net) {
  layers.resize(net->layers.size());
  for (int layer = 0; layer < net->layers.size(); layer++) {
    const Layer &cpu_layer = net->layers[layer];
    GPULayer *gpu_layer = &layers[layer];

    auto CopyMemoryAllowEmpty = [&](auto vec, bool readonly) -> cl_mem {
        if (vec.empty()) return 0;
        else return CopyMemoryToGPU(cl->context, cl->queue, vec, readonly);
      };

    gpu_layer->chunks.resize(cpu_layer.chunks.size());
    for (int chunk = 0; chunk < cpu_layer.chunks.size(); chunk++) {
      CHECK(chunk < cpu_layer.chunks.size());
      CHECK(chunk < gpu_layer->chunks.size());
      const Chunk &cpu_chunk = cpu_layer.chunks[chunk];
      GPUChunk *gpu_chunk = &gpu_layer->chunks[chunk];

      // Normal for this to be empty for dense chunks.
      gpu_chunk->indices = CopyMemoryAllowEmpty(cpu_chunk.indices, true);
      // These two should only be empty for the (token) input layer,
      // which we will never use, but we still build the gpu copies
      // for uniformity.
      gpu_chunk->weights = CopyMemoryAllowEmpty(cpu_chunk.weights, false);
      gpu_chunk->biases = CopyMemoryAllowEmpty(cpu_chunk.biases, false);

      // Empty when weight_update is SGD.
      gpu_chunk->weights_aux =
        CopyMemoryAllowEmpty(cpu_chunk.weights_aux, false);
      gpu_chunk->biases_aux =
        CopyMemoryAllowEmpty(cpu_chunk.biases_aux, false);

      InvertedIndices inverted = net->ComputeInvertedIndices(layer, chunk);
      gpu_chunk->ii_start = CopyMemoryAllowEmpty(inverted.start, true);
      gpu_chunk->ii_length = CopyMemoryAllowEmpty(inverted.length, true);
      gpu_chunk->ii_indices =
        CopyMemoryAllowEmpty(inverted.output_indices, true);
    }
  }

  clFinish(cl->queue);
}

NetworkGPU::~NetworkGPU() {
  for (GPULayer &gpu_layer : layers) {
    for (GPUChunk &gpu_chunk : gpu_layer.chunks) {
      if (gpu_chunk.indices != 0)
        CHECK_SUCCESS(clReleaseMemObject(gpu_chunk.indices));
      if (gpu_chunk.weights != 0)
        CHECK_SUCCESS(clReleaseMemObject(gpu_chunk.weights));
      if (gpu_chunk.biases != 0)
        CHECK_SUCCESS(clReleaseMemObject(gpu_chunk.biases));
      if (gpu_chunk.weights_aux != 0)
        CHECK_SUCCESS(clReleaseMemObject(gpu_chunk.weights_aux));
      if (gpu_chunk.biases_aux != 0)
        CHECK_SUCCESS(clReleaseMemObject(gpu_chunk.biases_aux));
      if (gpu_chunk.ii_start != 0)
        CHECK_SUCCESS(clReleaseMemObject(gpu_chunk.ii_start));
      if (gpu_chunk.ii_length != 0)
        CHECK_SUCCESS(clReleaseMemObject(gpu_chunk.ii_length));
      if (gpu_chunk.ii_indices != 0)
        CHECK_SUCCESS(clReleaseMemObject(gpu_chunk.ii_indices));
    }
  }
}


void NetworkGPU::ReadFromGPU() {
  CHECK(net->layers.size() == layers.size());
  for (int layer = 0; layer < net->layers.size(); layer++) {
    Layer *cpu_layer = &net->layers[layer];
    GPULayer *gpu_layer = &layers[layer];
    CHECK(cpu_layer->chunks.size() == gpu_layer->chunks.size());
    for (int chunk = 0; chunk < cpu_layer->chunks.size(); chunk++) {
      Chunk *cpu_chunk = &cpu_layer->chunks[chunk];
      GPUChunk *gpu_chunk = &gpu_layer->chunks[chunk];
      ReadToZeroOk(gpu_chunk->weights, &cpu_chunk->weights);
      ReadToZeroOk(gpu_chunk->biases, &cpu_chunk->biases);
      ReadToZeroOk(gpu_chunk->weights_aux, &cpu_chunk->weights_aux);
      ReadToZeroOk(gpu_chunk->biases_aux, &cpu_chunk->biases_aux);
    }
  }
  clFinish(cl->queue);
}


static string ForwardKernelName(ChunkType ct) {
  switch (ct) {
  case CHUNK_DENSE: return "ForwardChunkDense";
  case CHUNK_SPARSE: return "ForwardChunkSparse";
  case CHUNK_CONVOLUTION_ARRAY: return "ForwardChunkConvolutional";
  case CHUNK_INPUT:
    LOG(FATAL) << "Can't run INPUT chunks forward.";
  default:
    CHECK(false) << "Unsupported chunk type "
                 << ChunkTypeName(ct);
  }
  return "";
}


ForwardLayerCL::ForwardLayerCL(CL *cl, const Network &net) : cl(cl) {
  // Compile the appropriate kernel with baked in constants for
  // each chunk in the network.
  string base_src = Util::ReadFile("forwardchunk.cl");
  for (int layer = 0; layer < net.layers.size(); layer++) {
    if (layer == 0) {
      // No forward kernels for input layer.
      CHECK(net.layers[layer].chunks.size() == 1);
      layer_kernels.push_back({ChunkKernel()});
      continue;
    }

    const int src_layer_size = net.layers[layer - 1].num_nodes;
    const int dst_layer_size = net.layers[layer].num_nodes;

    std::vector<ChunkKernel> chunk_kernels;
    int chunk_start = 0;
    for (const Chunk &chunk : net.layers[layer].chunks) {
      string kernel_src =
        Network::TransferFunctionDefines(chunk.transfer_function);

      StringAppendF(&kernel_src,
                    "\n"
                    "#define INDICES_PER_NODE %d\n"
                    "#define NUM_FEATURES %d\n"
                    "#define SPAN_START %d\n"
                    "#define SPAN_SIZE %d\n"
                    "#define CHUNK_START %d\n"
                    "#define SRC_LAYER_SIZE %d\n"
                    "#define DST_LAYER_SIZE %d\n",
                    chunk.indices_per_node,
                    chunk.num_features,
                    chunk.span_start,
                    chunk.span_size,
                    chunk_start,
                    src_layer_size,
                    dst_layer_size);

      const string kernel_name = ForwardKernelName(chunk.type);
      kernel_src += base_src;
      auto [program, kernel] = cl->BuildOneKernel(kernel_src, kernel_name);
      CHECK(program != 0 && kernel != 0);

      ChunkKernel ck;
      ck.program = program;
      ck.kernel = kernel;
      chunk_kernels.push_back(ck);

      chunk_start += chunk.num_nodes;
    }
    layer_kernels.push_back(std::move(chunk_kernels));
  }
}

void ForwardLayerCL::RunForward(
    NetworkGPU *net_gpu, TrainingRoundGPU *train, int src_layer) {
  const int dst_layer = src_layer + 1;

  // TODO: Do we really want to share the same command queue across
  // threads? Presumably clFinish can't tell "this thread's
  // commands" apart from others, so we may be prematurely
  // waiting/running other thread's work.

  CHECK_LT(dst_layer, train->stimulations.size());

  // TODO: could keep net_gpu and net as members?
  const Network &net = *net_gpu->net;

  // We do the chunks in serial.
  const Layer &layer = net.layers[dst_layer];
  for (int chunk_idx = 0; chunk_idx < layer.chunks.size(); chunk_idx++) {

    const Chunk &chunk = layer.chunks[chunk_idx];
    CHECK(dst_layer < layer_kernels.size());
    CHECK(chunk_idx < layer_kernels[dst_layer].size());
    const ChunkKernel &ck = layer_kernels[dst_layer][chunk_idx];

    CHECK(chunk_idx < net_gpu->layers[dst_layer].chunks.size());
    NetworkGPU::GPUChunk &gpu_chunk =
      net_gpu->layers[dst_layer].chunks[chunk_idx];

    cl_mem indices = gpu_chunk.indices;
    CHECK(chunk.type == CHUNK_DENSE || indices != 0);
    cl_mem weights = gpu_chunk.weights;
    CHECK(weights != 0);
    cl_mem biases = gpu_chunk.biases;
    CHECK(biases != 0);

    // const int src_num_nodes = net.layers[src_layer].num_nodes;
    // size src_num_nodes * num_examples
    cl_mem all_src_values = train->stimulations[src_layer];


    // const int dst_num_nodes = net.layers[dst_layer].num_nodes;
    // size dst_num_nodes * num_examples
    cl_mem all_dst_values = train->stimulations[dst_layer];

    // Can't have multiple threads setting a kernel's argument at one time.
    {
      MutexLock ml(&m);

      // All the kernels take the same args.
      CHECK_SUCCESS(clSetKernelArg(ck.kernel, 0, sizeof (cl_mem),
                                   (void *)&all_src_values));
      // (Note that indices can be 0, an invalid cl_mem, which
      // we use to represent an empty memory. The memory will not
      // be accessed in this case (dense kernel). Not totally sure
      // that it is defined behavior to pass 0.)
      CHECK_SUCCESS(clSetKernelArg(ck.kernel, 1, sizeof (cl_mem),
                                   (void *)&indices));
      CHECK_SUCCESS(clSetKernelArg(ck.kernel, 2, sizeof (cl_mem),
                                   (void *)&weights));
      CHECK_SUCCESS(clSetKernelArg(ck.kernel, 3, sizeof (cl_mem),
                                   (void *)&biases));
      CHECK_SUCCESS(clSetKernelArg(ck.kernel, 4, sizeof (cl_mem),
                                   (void *)&all_dst_values));

      // TODO: Break into smaller chunks if nodes * examples is
      // too large?

      // Arg 0 is node, Arg 1 is example
      size_t global_work_offset[] = { 0, 0, };
      size_t global_work_size[] = {
        (size_t)(chunk.num_nodes),
        (size_t)(train->num_examples),
      };

      CHECK(ck.kernel != 0);
      CHECK_SUCCESS(clEnqueueNDRangeKernel(cl->queue, ck.kernel,
                                           // work dimensions
                                           2,
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
}

ForwardLayerCL::~ForwardLayerCL() {
  for (std::vector<ChunkKernel> &ckv : layer_kernels) {
    for (ChunkKernel &ck : ckv) {
      // These will be 0 for the input layer.
      if (ck.kernel != 0) CHECK_SUCCESS(clReleaseKernel(ck.kernel));
      if (ck.program != 0) CHECK_SUCCESS(clReleaseProgram(ck.program));
    }
  }
  layer_kernels.clear();
}



SetOutputErrorCL::SetOutputErrorCL(
    CL *cl, const Network &net,
    const std::optional<std::string> remap_define) : cl(cl) {
  // This only runs on one layer, the output. But we do need to have the
  // transfer function's derivative.
  string base_src = Util::ReadFile("setoutputerror.cl");

  const Layer &layer = net.layers.back();
  int chunk_start = 0;
  for (int chunk_idx = 0; chunk_idx < layer.chunks.size(); chunk_idx++) {
    const Chunk &chunk = layer.chunks[chunk_idx];
    const TransferFunction transfer_function =
      chunk.transfer_function;

    string kernel_src =
      Network::TransferFunctionDefines(transfer_function);

    // Add remapping function or fill in identity if disabled.
    kernel_src += "\n";
    if (remap_define.has_value()) {
      kernel_src += remap_define.value();
    } else {
      kernel_src += "#define REMAP(c, i, x) x";
    }
    kernel_src += "\n";

    StringAppendF(&kernel_src,
                  "\n"
                  "#define CHUNK_START %d\n"
                  "#define CHUNK_IDX %d\n"
                  "#define LAYER_SIZE %d\n",
                  chunk_start,
                  chunk_idx,
                  layer.num_nodes);

    kernel_src += base_src;

    auto pk = cl->BuildOneKernel(kernel_src, "SetOutputError");
    kernels.push_back(ChunkKernel{.program = pk.first,
                                  .kernel = pk.second});
    chunk_start += chunk.num_nodes;
  }
  CHECK(chunk_start == layer.num_nodes);
}

void SetOutputErrorCL::SetOutputError(
    NetworkGPU *net_gpu, TrainingRoundGPU *train) {
  // TODO: Could keep alias to this?
  const Network *net = net_gpu->net;

  // Full buffers for all examples.
  cl_mem all_actual_outputs = train->stimulations.back();
  cl_mem all_expected = train->expected;
  cl_mem all_output_error = train->errors.back();

  const Layer &layer = net->layers.back();
  for (int chunk_idx = 0; chunk_idx < layer.chunks.size(); chunk_idx++) {
    CHECK(chunk_idx < kernels.size());
    const Chunk &chunk = layer.chunks[chunk_idx];
    const ChunkKernel &ck = kernels[chunk_idx];

    {
      MutexLock ml(&m);

      CHECK_SUCCESS(clSetKernelArg(ck.kernel, 0, sizeof (cl_mem),
                                   (void *)&all_actual_outputs));
      CHECK_SUCCESS(clSetKernelArg(ck.kernel, 1, sizeof (cl_mem),
                                   (void *)&all_expected));
      CHECK_SUCCESS(clSetKernelArg(ck.kernel, 2, sizeof (cl_mem),
                                   (void *)&all_output_error));

      size_t global_work_offset[] = { 0, 0, };
      size_t global_work_size[] = {
        (size_t)chunk.num_nodes,
        (size_t)train->num_examples,
      };

      CHECK_SUCCESS(clEnqueueNDRangeKernel(cl->queue, ck.kernel,
                                           // work dimensions
                                           2,
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
}


static std::string Backward1KernelName(ChunkType ct) {
  switch (ct) {
  case CHUNK_DENSE: return "BackwardChunkDense";
  case CHUNK_SPARSE: return "BackwardChunkSparse";
  case CHUNK_CONVOLUTION_ARRAY: return "BackwardChunkConvolutional";
  default:
    CHECK(false) << "Unsupported chunk type for BackwardLayer1 "
                 << ChunkTypeName(ct);
  }
  return "ERROR";
}

BackwardLayerCL::BackwardLayerCL(CL *cl, const Network &net) : cl(cl) {
  string base1_src = Util::ReadFile("backwardchunk.cl");

  // Dummy kernels for input layer, which can't be a destination layer.
  CHECK(net.layers[0].chunks.size() == 1);
  layer_kernels.push_back(std::vector<ChunkKernel>{ChunkKernel()});

  for (int dst_layer_idx = 1;
       dst_layer_idx < net.layers.size();
       dst_layer_idx++) {
    const Layer &dst_layer = net.layers[dst_layer_idx];

    // First pass is chunk-by-chunk in the destination layer.
    std::vector<ChunkKernel> chunk_kernels;
    int out_idx = 0;
    for (int chunk_idx = 0; chunk_idx < dst_layer.chunks.size(); chunk_idx++) {
      const Chunk &chunk = dst_layer.chunks[chunk_idx];
      string kernel_src =
        StringPrintf("#define SRC_LAYER_SIZE %d\n"
                     "#define DST_LAYER_SIZE %d\n"
                     "#define CHUNK_START %d\n"
                     "#define SPAN_START %d\n"
                     "#define SPAN_SIZE %d\n"
                     "#define DST_INDICES_PER_NODE %d\n"
                     "#define DST_NUM_NODES %d\n"
                     "#define DST_NUM_FEATURES %d\n"
                     "#define SRC_SPAN_IS_ZERO %d\n",
                     net.layers[dst_layer_idx - 1].num_nodes,
                     net.layers[dst_layer_idx].num_nodes,
                     out_idx,
                     chunk.span_start,
                     chunk.span_size,
                     chunk.indices_per_node,
                     chunk.num_nodes,
                     chunk.num_features,
                     // PERF: set this to true for first chunk,
                     // or (better) for all chunks writing into
                     // a span that hasn't yet been written by
                     // a previous chunk
                     false);

      kernel_src += base1_src;

      auto [program, kernel] =
        cl->BuildOneKernel(kernel_src, Backward1KernelName(chunk.type));

      ChunkKernel ck{
        .program1 = program,
        .kernel1 = kernel
        // program2, kernel2 filled in below.
      };
      chunk_kernels.push_back(ck);

      out_idx += chunk.num_nodes;
    }
    CHECK(out_idx == dst_layer.num_nodes);
    layer_kernels.push_back(std::move(chunk_kernels));
  }

  // And the second pass. This is associated with the source layer.
  string base2_src = Util::ReadFile("backwardsecondpass.cl");

  // We don't actually need this for the first layer (as we don't
  // propagate error to the input normally) nor the last (it
  // cannot be a source layer) but we generate it for each chunk
  // nonetheless, for uniformity.
  for (int layer_idx = 0;
       layer_idx < net.layers.size();
       layer_idx++) {
    const Layer &layer = net.layers[layer_idx];
    // We modify these in place so we expect them to have already been
    // created above.
    CHECK(layer_idx < layer_kernels.size());

    // Second pass is on the source layer.
    const int src_layer_size = net.layers[layer_idx].num_nodes;

    int out_idx = 0;
    for (int chunk_idx = 0; chunk_idx < layer.chunks.size(); chunk_idx++) {
      const Chunk &chunk = layer.chunks[chunk_idx];
      CHECK(chunk_idx < layer_kernels[layer_idx].size());

      string kernel_src =
        Network::TransferFunctionDefines(chunk.transfer_function);

      StringAppendF(&kernel_src,
                    "#define SRC_LAYER_SIZE %d\n"
                    "#define CHUNK_START %d\n"
                    "#define CLIP_ERROR %s\n"
                    "#define LARGE_ERROR %0.8ff\n",
                    src_layer_size,
                    out_idx,
                    CLIP_ERROR ? "true" : "false",
                    LARGE_ERROR);

      kernel_src += base2_src;

      // One kernel for all chunk types (but it does depend on
      // the chunk's transfer function).
      auto [program, kernel] =
        cl->BuildOneKernel(kernel_src, "BackwardSecondPass");
      layer_kernels[layer_idx][chunk_idx].program2 = program;
      layer_kernels[layer_idx][chunk_idx].kernel2 = kernel;

      out_idx += chunk.num_nodes;
    }
    CHECK(out_idx == layer.num_nodes);
  }
}

BackwardLayerCL::~BackwardLayerCL() {
  for (auto &v : layer_kernels) {
    for (auto &ck : v) {
      if (ck.kernel1 != 0) CHECK_SUCCESS(clReleaseKernel(ck.kernel1));
      if (ck.program1 != 0) CHECK_SUCCESS(clReleaseProgram(ck.program1));
      if (ck.kernel2 != 0) CHECK_SUCCESS(clReleaseKernel(ck.kernel2));
      if (ck.program2 != 0) CHECK_SUCCESS(clReleaseProgram(ck.program2));
    }
  }
}

void BackwardLayerCL::BackwardLayer(NetworkGPU *net_gpu,
                                    TrainingRoundGPU *train,
                                    int dst_layer) {
  const Network &net = *net_gpu->net;
  CHECK(dst_layer > 0);
  CHECK(dst_layer < net.layers.size());

  const int src_layer = dst_layer - 1;

  {
    // Full source error for all examples.
    cl_mem all_src_error = train->errors[src_layer];

    // In the general case we are accumulating the weighted error sum
    // (+=) rather than writing it once (=), so unlike the other
    // kernels we need to start by clearing the src_error.
    //
    // PERF we could consider doing this only for chunks where we have
    // SRC_SPAN_IS_ZERO false, although it gets pretty complicated with
    // multiple training examples. Better perhaps to just skip when ALL
    // the chunks are SRC_SPAN_IS_ZERO, which would happen in the common
    // case that there is just one chunk.
    cl_float zero = 0.0f;
    CHECK_SUCCESS(
        clEnqueueFillBuffer(cl->queue,
                            all_src_error,
                            // pattern and its size in bytes
                            &zero, sizeof (cl_float),
                            // offset and size to fill (in BYTES)
                            0, (size_t)(net.layers[src_layer].num_nodes *
                                        train->num_examples *
                                        sizeof (cl_float)),
                            // no wait list or event
                            0, nullptr, nullptr));
    // This needs to be done before the kernel runs below. PERF that
    // if we are running many examples in parallel, we'll be synchronizing
    // on each other's writes, and possibly starving the kernel below. (This
    // is probably no longer an issue since the training round is managed
    // internally?)
    // This would be a good place to use wait lists!
    clFinish(cl->queue);
  }

  // First pass over chunks in the DEST layer.
  const Layer &layer = net.layers[dst_layer];
  for (int chunk_idx = 0; chunk_idx < layer.chunks.size(); chunk_idx++) {
    const Chunk &chunk = layer.chunks[chunk_idx];
    NetworkGPU::GPUChunk &gpu_chunk =
      net_gpu->layers[dst_layer].chunks[chunk_idx];

    cl_mem ii_start = gpu_chunk.ii_start;
    cl_mem ii_length = gpu_chunk.ii_length;
    cl_mem ii_indices = gpu_chunk.ii_indices;
    cl_mem dst_weights = gpu_chunk.weights;

    CHECK(dst_layer < layer_kernels.size() &&
          chunk_idx < layer_kernels[dst_layer].size());
    const ChunkKernel &ck = layer_kernels[dst_layer][chunk_idx];

    cl_mem all_dst_error = train->errors[dst_layer];
    cl_mem all_src_error = train->errors[src_layer];

    {
      MutexLock ml(&m);

      CHECK_SUCCESS(clSetKernelArg(ck.kernel1, 0, sizeof (cl_mem),
                                   (void *)&ii_start));
      CHECK_SUCCESS(clSetKernelArg(ck.kernel1, 1, sizeof (cl_mem),
                                   (void *)&ii_length));
      CHECK_SUCCESS(clSetKernelArg(ck.kernel1, 2, sizeof (cl_mem),
                                   (void *)&ii_indices));
      CHECK_SUCCESS(clSetKernelArg(ck.kernel1, 3, sizeof (cl_mem),
                                   (void *)&dst_weights));
      CHECK_SUCCESS(clSetKernelArg(ck.kernel1, 4, sizeof (cl_mem),
                                   (void *)&all_dst_error));
      CHECK_SUCCESS(clSetKernelArg(ck.kernel1, 5, sizeof (cl_mem),
                                   (void *)&all_src_error));

      size_t global_work_offset[] = { 0, 0, };
      size_t global_work_size[] = {
        (size_t)chunk.span_size,
        (size_t)train->num_examples,
      };
      Timer kernel_timer;
      CHECK_SUCCESS(clEnqueueNDRangeKernel(cl->queue, ck.kernel1,
                                           // work dimensions
                                           2,
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

  // Second pass over chunks in the SOURCE layer.
  const Layer &source_layer = net.layers[src_layer];
  for (int chunk_idx = 0; chunk_idx < source_layer.chunks.size(); chunk_idx++) {
    const Chunk &chunk = source_layer.chunks[chunk_idx];

    CHECK(src_layer < layer_kernels.size() &&
          chunk_idx < layer_kernels[src_layer].size());
    const ChunkKernel &ck = layer_kernels[src_layer][chunk_idx];

    // Full source error for this example.
    cl_mem all_src_output = train->stimulations[src_layer];
    cl_mem all_src_error = train->errors[src_layer];

    {
      MutexLock ml(&m);

      CHECK_SUCCESS(clSetKernelArg(ck.kernel2, 0, sizeof (cl_mem),
                                   (void *)&all_src_output));
      CHECK_SUCCESS(clSetKernelArg(ck.kernel2, 1, sizeof (cl_mem),
                                   (void *)&all_src_error));

      size_t global_work_offset[] = { 0, 0 };
      size_t global_work_size[] = {
        (size_t)chunk.num_nodes,
        (size_t)train->num_examples,
      };
      Timer kernel_timer;
      CHECK_SUCCESS(clEnqueueNDRangeKernel(cl->queue, ck.kernel2,
                                           // work dimensions
                                           2,
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
}

DecayWeightsCL::DecayWeightsCL(CL *cl, const Network &net,
                               float decay_factor) : cl(cl) {
  string base_src = Util::ReadFile("decayweights.cl");

  string kernel_src;
  StringAppendF(&kernel_src, "\n#define DECAY_FACTOR %.9ff\n",
                decay_factor);
  kernel_src += base_src;
  auto p = cl->BuildOneKernel(kernel_src, "DecayWeights");
  program = p.first;
  kernel = p.second;
}

DecayWeightsCL::~DecayWeightsCL() {
  CHECK_SUCCESS(clReleaseKernel(kernel));
  CHECK_SUCCESS(clReleaseProgram(program));
}

void DecayWeightsCL::Decay(NetworkGPU *net_gpu, int layer_idx) {
  const Network &net = *net_gpu->net;
  // PERF: Should actually be able to run in parallel across the entire
  // network if we weren't sharing a single kernel. Every weight
  // just scaled independently.
  CHECK(layer_idx >= 0 && layer_idx < net.layers.size());
  for (int chunk_idx = 0;
       chunk_idx < net.layers[layer_idx].chunks.size();
       chunk_idx++) {
    const Chunk &cpu_chunk = net.layers[layer_idx].chunks[chunk_idx];
    if (cpu_chunk.fixed) continue;

    NetworkGPU::GPUChunk &gpu_chunk =
      net_gpu->layers[layer_idx].chunks[chunk_idx];
    const int num_weights =
      net.layers[layer_idx].chunks[chunk_idx].weights.size();

    {
      MutexLock ml(&m);
      CHECK_SUCCESS(clSetKernelArg(kernel, 0, sizeof (cl_mem),
                                   (void *)&gpu_chunk.weights));

      size_t global_work_offset[] = { 0 };
      // Total number of weights in this chunk.
      size_t global_work_size[] = { (size_t)num_weights };
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
}

static string UpdateWeights1KernelName(ChunkType ct) {
  switch (ct) {
    case CHUNK_DENSE: return "UpdateWeightsDense";
    case CHUNK_SPARSE: return "UpdateWeightsSparse";
    case CHUNK_CONVOLUTION_ARRAY: return "UpdateWeightsConvolutional";
    case CHUNK_INPUT:
      CHECK(false) << "Can't update weights for input layer.";
      return "";
    default:
      CHECK(false) << "Unsupported chunk type "
                   << ChunkTypeName(ct);
      return "";
  }
}


UpdateWeightsCL::UpdateWeightsCL(int examples_per_round,
                                 CL *cl, const Network &net) :
  examples_per_round(examples_per_round), cl(cl) {
  // Note that this one doesn't depend on the transfer function/derivative.

  // Also unlike others, we actually invoke the kernel differently in
  // the convolution case.

  const string base_src1 = Util::ReadFile("updateweights.cl");
  const string base_src2 = Util::ReadFile("updatesecondpass.cl");

  // Not possible on the input layer, but we have placeholder
  // ChunkKernels for uniformity with other layer arrays.
  CHECK(net.layers[0].chunks.size() == 1);
  layer_kernels.push_back(std::vector<ChunkKernel>{ChunkKernel()});

  for (int layer_idx = 1; layer_idx < net.layers.size(); layer_idx++) {

    std::vector<ChunkKernel> chunk_kernels;
    int out_idx = 0;
    for (int chunk_idx = 0;
         chunk_idx < net.layers[layer_idx].chunks.size();
         chunk_idx++) {
      const Chunk &chunk = net.layers[layer_idx].chunks[chunk_idx];
      ChunkKernel ck;

      const int num_occurrences =
        chunk.num_occurrences_across *
        chunk.num_occurrences_down;

      // PERF: Don't even compile a kernel if the chunk is fixed.

      string kernel1_src;

      StringAppendF(&kernel1_src,
                    "#define SRC_LAYER_SIZE %d\n"
                    "#define DST_LAYER_SIZE %d\n"
                    "#define CHUNK_START %d\n"
                    "#define SPAN_START %d\n"
                    "#define SPAN_SIZE %d\n"
                    "#define INDICES_PER_NODE %d\n"
                    "#define NUM_OCCURRENCES %d\n"
                    "#define NUM_FEATURES %d\n",
                    net.layers[layer_idx - 1].num_nodes,
                    net.layers[layer_idx].num_nodes,
                    out_idx,
                    chunk.span_start,
                    chunk.span_size,
                    chunk.indices_per_node,
                    num_occurrences,
                    chunk.num_features);

      const string kernel1_name = UpdateWeights1KernelName(chunk.type);

      kernel1_src += base_src1;
      auto pk1 = cl->BuildOneKernel(kernel1_src, kernel1_name);
      ck.program1 = pk1.first;
      ck.kernel1 = pk1.second;

      // Second pass.
      string kernel2_src;
      switch (chunk.weight_update) {
      case SGD:
        kernel2_src += "#define WEIGHT_UPDATE_SGD 1\n";
        break;
      case ADAM:
        kernel2_src += "#define WEIGHT_UPDATE_ADAM 1\n";
        break;
      default:
        LOG(FATAL) << "Unsupported weight update type " <<
          WeightUpdateName(chunk.weight_update);
        break;
      }

      StringAppendF(&kernel2_src,
                    "#define EXAMPLES_PER_ROUND %d\n"
                    "#define CLIPPING %s\n"
                    "#define CONSTRAIN %s\n",
                    examples_per_round,
                    CLIPPING ? "true" : "false",
                    CONSTRAIN ? "true" : "false");


      kernel2_src += base_src2;
      auto pk2 = cl->BuildOneKernel(kernel2_src,
                                   "UpdateWeightsSecondPass");
      ck.program2 = pk2.first;
      ck.kernel2 = pk2.second;


      // Allocate temporary buffer.
      // TODO: This is the minimum size. Allow for multiples up
      // to some configured budget, then run in parallel.
      const int num_weights =
        net.layers[layer_idx].chunks[chunk_idx].weights.size();
      const int num_biases =
        net.layers[layer_idx].chunks[chunk_idx].biases.size();
      ck.weight_grad_tmp =
        CreateUninitializedGPUMemory<float>(cl->context, num_weights);
      ck.bias_grad_tmp =
        CreateUninitializedGPUMemory<float>(cl->context, num_biases);

      chunk_kernels.push_back(ck);
      out_idx += chunk.num_nodes;
    }
    CHECK(out_idx == net.layers[layer_idx].num_nodes);

    layer_kernels.push_back(std::move(chunk_kernels));
  }
}

UpdateWeightsCL::~UpdateWeightsCL() {
  for (auto &vec : layer_kernels) {
    for (auto &ck : vec) {
      if (ck.kernel1 != 0) CHECK_SUCCESS(clReleaseKernel(ck.kernel1));
      if (ck.program1 != 0) CHECK_SUCCESS(clReleaseProgram(ck.program1));
      if (ck.kernel2 != 0) CHECK_SUCCESS(clReleaseKernel(ck.kernel2));
      if (ck.program2 != 0) CHECK_SUCCESS(clReleaseProgram(ck.program2));

      if (ck.weight_grad_tmp != 0)
        CHECK_SUCCESS(clReleaseMemObject(ck.weight_grad_tmp));
      if (ck.bias_grad_tmp != 0)
        CHECK_SUCCESS(clReleaseMemObject(ck.bias_grad_tmp));
    }
  }
}

void UpdateWeightsCL::Update(NetworkGPU *net_gpu, TrainingRoundGPU *train,
                             float learning_rate, int layer_idx) {
  // new chunks
  CHECK(layer_idx > 0) << "Can't update the weights for the input layer, "
    "which would not be useful anyway since there aren't any.";
  CHECK(train->num_examples == examples_per_round) << "Must match the "
    "configured constant.";

  const Network &net = *net_gpu->net;

  const Layer &layer = net.layers[layer_idx];
  NetworkGPU::GPULayer &gpu_layer = net_gpu->layers[layer_idx];

  // XXX overflow is possible here.
  // PERF: After the round is suitably high, we should be ignoring this.
  const cl_int round_number = net.rounds;

  for (int chunk_idx = 0; chunk_idx < layer.chunks.size(); chunk_idx++) {
    const Chunk &chunk = layer.chunks[chunk_idx];
    // For fixed chunks, just skip the update step.
    if (chunk.fixed) continue;

    NetworkGPU::GPUChunk &gpu_chunk = gpu_layer.chunks[chunk_idx];
    ChunkKernel &ck = layer_kernels[layer_idx][chunk_idx];

    cl_mem all_prev_layer_output = train->stimulations[layer_idx - 1];
    cl_mem all_layer_error = train->errors[layer_idx];

    // TODO PERF: These all write to the same weights, so it is much
    // trickier to parallelize over examples than other passes. We
    // could run each example to a temporary weight vector
    // and then sum them at the end. Probably we do
    // not want to use num_examples * num_weights in temporary memory
    // for this. But we could pick some width W and do W examples in
    // parallel, each writing to their own shards, and then sum those,
    // perhaps into an accumulator of the same size. Here we would
    // only need W * num_weights memory, and this is fine because as
    // num_weights gets high, we no longer benefit from internal
    // parallelism.
    // If we have W buffers, it could work like this:
    //   clear buffer[0].
    //   while there are examples left:
    //     clear buffer[1]...buffer[W-1]  (but not the 0th)
    //     run the next W examples, writing with += to buffer[0]...buffer[W-1]
    //     sum buffer[0]...buffer[W-1] into buffer[0]
    //   // now buffer[0] holds the sum of weight updates
    //   apply clipping, adam etc., write to weights_aux and weights.
    // We do need W>0, so this requires at least that extra memory
    // overhead. It could be shared between chunks/layers if we're willing
    // to give up on parallelism over those.

    // (XXX BTW, it now occurs to me that I'm currently just doing
    // Adam wrong, because I run it 1000 times per round (and update
    // the moving averages) -- but I should only be taking one
    // timestep. Right? So this would recommend a two-pass approach,
    // computing the sums and adding at the end. We would not need to
    // keep copies of weight_aux, since we only do this in the second
    // pass.)

    // TODO: We want to use more than one temp buffer for parallelism!
    const int num_weights =
      net.layers[layer_idx].chunks[chunk_idx].weights.size();
    const int num_biases =
      net.layers[layer_idx].chunks[chunk_idx].biases.size();

    // Accumulate all weight grads into these temp buffers, then
    // apply them below.
    {
      // Clear temp memory.
      cl_float zero = 0.0f;
      CHECK_SUCCESS(
          clEnqueueFillBuffer(cl->queue,
                              ck.weight_grad_tmp,
                              // pattern and its size in bytes
                              &zero, sizeof (cl_float),
                              // offset and size to fill (in BYTES)
                              0, (size_t)(num_weights * sizeof (cl_float)),
                              // no wait list or event
                              0, nullptr, nullptr));
      CHECK_SUCCESS(
          clEnqueueFillBuffer(cl->queue,
                              ck.bias_grad_tmp,
                              // pattern and its size in bytes
                              &zero, sizeof (cl_float),
                              // offset and size to fill (in BYTES)
                              0, (size_t)(num_biases * sizeof (cl_float)),
                              // no wait list or event
                              0, nullptr, nullptr));
      clFinish(cl->queue);
    }

    for (int example_num = 0;
         example_num < train->num_examples;
         example_num++) {

      cl_int example_num_cl = example_num;

      {
        // PERF: Could run the chunks and/or layers in parallel.
        MutexLock ml(&m);

        CHECK_SUCCESS(clSetKernelArg(ck.kernel1, 0, sizeof (cl_mem),
                                     (void *)&all_prev_layer_output));
        CHECK_SUCCESS(clSetKernelArg(ck.kernel1, 1, sizeof (cl_mem),
                                     (void *)&all_layer_error));
        CHECK_SUCCESS(clSetKernelArg(ck.kernel1, 2, sizeof (cl_mem),
                                     (void *)&gpu_chunk.indices));
        CHECK_SUCCESS(clSetKernelArg(ck.kernel1, 3, sizeof (cl_mem),
                                     (void *)&ck.weight_grad_tmp));
        CHECK_SUCCESS(clSetKernelArg(ck.kernel1, 4, sizeof (cl_mem),
                                     (void *)&ck.bias_grad_tmp));
        CHECK_SUCCESS(clSetKernelArg(ck.kernel1, 5, sizeof (cl_int),
                                     (void *)&example_num_cl));

        // Arguments are the same for the different chunk types,
        // but for convolution array it's actually a 2D kernel
        // over features and weights.

        if (chunk.type == CHUNK_CONVOLUTION_ARRAY) {
          size_t global_work_offset[] = { 0, 0 };
          size_t global_work_size[] = { (size_t)chunk.num_features,
                                        (size_t)chunk.indices_per_node };
          CHECK_SUCCESS(clEnqueueNDRangeKernel(cl->queue, ck.kernel1,
                                               // work dimensions
                                               2,
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
        } else {
          size_t global_work_offset[] = { 0 };
          size_t global_work_size[] = { (size_t)chunk.num_nodes };
          CHECK_SUCCESS(clEnqueueNDRangeKernel(cl->queue, ck.kernel1,
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
        }
        clFinish(cl->queue);
      }
    }

    // Second pass for the chunk here.
    {
      MutexLock ml(&m);

      auto SecondPass = [&](int num_weights,
                            cl_mem grad_sums,
                            cl_mem chunk_weights,
                            cl_mem chunk_weights_aux,
                            cl_float constrain_max) {
          CHECK_SUCCESS(clSetKernelArg(ck.kernel2, 0, sizeof (cl_int),
                                       (void *)&round_number));
          CHECK_SUCCESS(clSetKernelArg(ck.kernel2, 1, sizeof (cl_float),
                                       (void *)&learning_rate));
          CHECK_SUCCESS(clSetKernelArg(ck.kernel2, 2, sizeof (cl_float),
                                       (void *)&constrain_max));
          CHECK_SUCCESS(clSetKernelArg(ck.kernel2, 3, sizeof (cl_mem),
                                       (void *)&grad_sums));
          CHECK_SUCCESS(clSetKernelArg(ck.kernel2, 4, sizeof (cl_mem),
                                       (void *)&chunk_weights));
          CHECK_SUCCESS(clSetKernelArg(ck.kernel2, 5, sizeof (cl_mem),
                                       (void *)&chunk_weights_aux));

          size_t global_work_offset[] = { 0 };
          size_t global_work_size[] = { (size_t)num_weights };
          CHECK_SUCCESS(clEnqueueNDRangeKernel(cl->queue, ck.kernel2,
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
        };

      SecondPass(num_weights,
                 ck.weight_grad_tmp,
                 gpu_chunk.weights,
                 gpu_chunk.weights_aux,
                 WEIGHT_CONSTRAIN_MAX);

      SecondPass(num_biases,
                 ck.bias_grad_tmp,
                 gpu_chunk.biases,
                 gpu_chunk.biases_aux,
                 BIAS_CONSTRAIN_MAX);
    }
  }
}
