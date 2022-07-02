
#include "network-gpu.h"

#include <cmath>
#include <memory>
#include <vector>
#include <functional>
#include <string>
#include <ctype.h>
#include <chrono>
#include <thread>
#include <deque>
#include <numbers>

#include "network.h"
#include "network-test-util.h"
#include "clutil.h"
#include "base/logging.h"
#include "base/stringprintf.h"
#include "arcfour.h"
#include "randutil.h"
#include "threadutil.h"
#include "image.h"
#include "util.h"
#include "train-util.h"
#include "periodically.h"
#include "timer.h"
#include "error-history.h"

#include "mnist.h"

using namespace std;

using TestNet = NetworkTestUtil::TestNet;
using TrainNet = NetworkTestUtil::TrainNet;
using TestExample = NetworkTestUtil::TestExample;

static CL *cl = nullptr;

using int64 = int64_t;

static constexpr WeightUpdate WEIGHT_UPDATE = ADAM;

static constexpr float TWO_PI = 2.0f * std::numbers::pi_v<float>;

#define MODEL_BASE "grad"
#define MODEL_NAME MODEL_BASE ".val"

// Known dimensions for mnist training data.
static constexpr int IMG_WIDTH = 28;
static constexpr int IMG_HEIGHT = 28;

constexpr int INPUT_SIZE = IMG_WIDTH * IMG_HEIGHT;
constexpr int OUTPUT_SIZE = 10;

// Very small examples -- this can probably be much higher.
static constexpr int EXAMPLES_PER_ROUND = 1000;

struct ExampleThread {
  // A full round's worth of examples.
  // First element is INPUT_SIZE * EXAMPLES_PER_ROUND,
  // then OUTPUT_SIZE * EXAMPLES_PER_ROUND.
  using example_type = std::pair<std::vector<float>, std::vector<float>>;

  example_type GetExamples(int64 round_num) {
    for (;;) {
      {
        MutexLock ml(&m);
        auto it = q.find(round_num);
        if (it != q.end()) {
          example_type ret = std::move(it->second);
          q.erase(it);
          return ret;
        }
      }

      std::this_thread::sleep_for(std::chrono::microseconds(100));
    }
  }

  ExampleThread(int64 next_round) : next_example(next_round) {
    mnist.reset(new MNIST("mnist/train"));
    CHECK(mnist.get() != nullptr);
    CHECK(mnist->width == IMG_WIDTH);
    CHECK(mnist->height == IMG_HEIGHT);

    work_thread.reset(new std::thread(&Generate, this, 1));
  }

  ~ExampleThread() {
    LOG(FATAL) << "unimplemented";
  }

private:
  static constexpr int TARGET_SIZE = 10;

  void Generate(int id) {
    printf("Started example thread %d.\n", id);

    // Not deterministic.
    ArcFour rc(StringPrintf("%lld.ex", time(nullptr)));

    for (;;) {
      // If we need to generate an example, claim its number.
      // Otherwise, return -1.
      const int64 next = [&]() -> int64 {
          MutexLock ml(&m);
          if (q.size() < TARGET_SIZE) {
            int64 ex = next_example;
            next_example++;
            return ex;
          } else {
            return -1;
          }
        }();

      if (next >= 0) {
        // examples_per_round * INPUT_SIZE
        std::vector<float> examples;
        // examples_per_round * OUTPUT_SIZE
        std::vector<float> outputs;
        outputs.reserve(EXAMPLES_PER_ROUND * OUTPUT_SIZE);

        // Fill 'em
        for (int i = 0; i < EXAMPLES_PER_ROUND; i++) {

          const int idx = RandTo(&rc, mnist->Num());

          // one-hot output label
          int lab = mnist->labels[idx];
          for (int o = 0; o < OUTPUT_SIZE; o++) {
            outputs.push_back(lab == o ? 1.0f : 0.0f);
          }

          // image data.
          // TODO: Randomly add noise, squish, pan, etc.
          const ImageA &img = mnist->images[idx];
          for (int y = 0; y < IMG_HEIGHT; y++) {
            for (int x = 0; x < IMG_WIDTH; x++) {
              float f = (float)img.GetPixel(x, y) / 255.0f;
              examples.push_back(f);
            }
          }
        }

        CHECK(examples.size() == EXAMPLES_PER_ROUND * INPUT_SIZE);
        CHECK(outputs.size() == EXAMPLES_PER_ROUND * OUTPUT_SIZE);

        {
          MutexLock ml(&m);
          CHECK(q.find(next) == q.end()) << next;
          q[next] = make_pair(std::move(examples), std::move(outputs));
        }
      } else {
        std::this_thread::sleep_for(std::chrono::milliseconds(25));
      }
    }
  }

  std::unique_ptr<MNIST> mnist;

  std::mutex m;
  // Round number for the next example to generate.
  int64 next_example = 0;
  std::map<int64, example_type> q;

  std::unique_ptr<std::thread> work_thread;
};

static void Train(Network *net) {
  ExampleThread example_thread(net->rounds);

  ErrorHistory error_history("error-history.tsv");

  static constexpr int max_parallelism = 4;
  // 0, 1, 2
  static constexpr int VERBOSE = 1;
  static constexpr bool SAVE_INTERMEDIATE = true;

  UpdateWeightsCL::UpdateConfig update_config;
  update_config.base_learning_rate = 0.1f;
  // This is conservative, but with larger exponents I would
  // get divergence after hundreds of thousands of rounds.
  // This happened again with the plugin parameter predictor
  // with a value of 1e-6!
  update_config.adam_epsilon = 1.0e-4;

  // XXX this should probably depend on the learning rate; if the
  // learning rate is too small, it won't even be able to overcome
  // the decay
  static constexpr float DECAY_RATE = 0.999999f;
  static constexpr bool DO_DECAY = false;

  // On a verbose round we compute training error and print out
  // examples.
  constexpr int VERBOSE_EVERY = 1000;
  int64 total_verbose = 0;

  // How often to also save error history to disk (if we run a
  // verbose round).
  static constexpr double HISTORY_EVERY_SEC = 120.0;
  Periodically history_per(HISTORY_EVERY_SEC);

  static constexpr double SAVE_ERROR_IMAGES_EVERY = 130.0;
  Periodically error_images_per(SAVE_ERROR_IMAGES_EVERY);

  static constexpr double TIMING_EVERY_SEC = 20.0;
  Periodically timing_per(TIMING_EVERY_SEC);


  static constexpr int64 CHECKPOINT_EVERY_ROUNDS = 100000;

  constexpr int IMAGE_EVERY = 100;
  TrainingImages images(*net, "train", MODEL_NAME, IMAGE_EVERY);

  printf("Training!\n");

  constexpr int NUM_COLUMNS = 10;
  std::vector<std::unique_ptr<ErrorImage>> error_images;
  for (int i = 0; i < NUM_COLUMNS; i++) {
    string filename = StringPrintf("error-%d.png", i);
    error_images.emplace_back(
        std::make_unique<ErrorImage>(2000, EXAMPLES_PER_ROUND, filename));
  }

  auto net_gpu = make_unique<NetworkGPU>(cl, net);

  std::unique_ptr<ForwardLayerCL> forward_cl =
    std::make_unique<ForwardLayerCL>(cl, net_gpu.get());
  std::unique_ptr<SetOutputErrorCL> error_cl =
    std::make_unique<SetOutputErrorCL>(cl, net_gpu.get());
  std::unique_ptr<BackwardLayerCL> backward_cl =
    std::make_unique<BackwardLayerCL>(cl, net_gpu.get());
  [[maybe_unused]]
  std::unique_ptr<DecayWeightsCL> decay_cl =
    std::make_unique<DecayWeightsCL>(cl, net_gpu.get(), DECAY_RATE);
  std::unique_ptr<UpdateWeightsCL> update_cl =
    std::make_unique<UpdateWeightsCL>(cl, net_gpu.get(),
                                      EXAMPLES_PER_ROUND,
                                      update_config);

  // Uninitialized training examples on GPU.
  std::unique_ptr<TrainingRoundGPU> training(
      new TrainingRoundGPU(EXAMPLES_PER_ROUND, cl, *net));

  CHECK(!net->layers.empty());
  CHECK(net->layers[0].num_nodes == INPUT_SIZE);
  CHECK(net->layers.back().num_nodes == OUTPUT_SIZE);

  double round_ms = 0.0;
  double getexample_ms = 0.0;
  double example_ms = 0.0;
  double forward_ms = 0.0;
  double error_ms = 0.0;
  double decay_ms = 0.0;
  double backward_ms = 0.0;
  double update_ms = 0.0;
  double loss_ms = 0.0;
  double image_ms = 0.0;

  Timer train_timer;
  int64 total_examples = 0LL;
  // seconds since timer started
  double last_save = 0.0;
  for (int iter = 0; true; iter++) {
    Timer round_timer;

    const bool verbose_round = (iter % VERBOSE_EVERY) == 0;

    // Initialize training examples.
    // (PERF: parallelize?)

    Timer getexample_timer;
    // All examples for the round, flat.
    auto [inputs, expecteds] = example_thread.GetExamples(net->rounds);

    getexample_ms += getexample_timer.MS();

    {
      Timer example_timer;
      CHECK(inputs.size() == INPUT_SIZE * EXAMPLES_PER_ROUND);
      CHECK(expecteds.size() == OUTPUT_SIZE * EXAMPLES_PER_ROUND);
      training->LoadInputs(inputs);
      training->LoadExpecteds(expecteds);
      example_ms += example_timer.MS();
    }

    if (VERBOSE > 1)
      printf("Prepped examples.\n");

    {
      Timer forward_timer;
      for (int src_layer = 0;
           src_layer < net->layers.size() - 1;
           src_layer++) {
        forward_cl->RunForward(training.get(), src_layer);
      }
      forward_ms += forward_timer.MS();
    }

    if (VERBOSE > 1)
      printf("Forward done.\n");

    {
      Timer error_timer;
      error_cl->SetOutputError(training.get());
      error_ms += error_timer.MS();
    }

    if (VERBOSE > 1)
      printf("Set error.\n");

    {
      Timer backward_timer;
      for (int dst_layer = net->layers.size() - 1;
           // Don't propagate to input.
           dst_layer > 1;
           dst_layer--) {
        backward_cl->BackwardLayer(training.get(), dst_layer);
      }
      backward_ms += backward_timer.MS();
    }

    if (VERBOSE > 1)
      printf("Backward pass.\n");

    if (DO_DECAY) {
      Timer decay_timer;
      for (int layer_idx = 0; layer_idx < net->layers.size(); layer_idx++) {
        decay_cl->Decay(layer_idx);
      }
      decay_ms += decay_timer.MS();
    }

    {
      Timer update_timer;
      // Can't run training examples in parallel because these all write
      // to the same network. But each layer is independent.
      UnParallelComp(net->layers.size() - 1,
                     [&](int layer_minus_1) {
                       const int layer_idx = layer_minus_1 + 1;
                       update_cl->Update(training.get(),
                                         layer_idx);
                     },
                     max_parallelism);
      update_ms += update_timer.MS();
    }

    if (VERBOSE > 1)
      printf("Updated weights.\n");

    total_examples += EXAMPLES_PER_ROUND;

    net->examples += EXAMPLES_PER_ROUND;
    net->rounds++;

    // (currently no way to actually finish, but we could set a
    // training error threshold below.)
    const bool finished = false;

    if (verbose_round) {
      Timer loss_timer;
      if (VERBOSE > 1)
        printf("Verbose round...\n");

      // Get loss as abs distance, plus number of incorrect (as booleans).

      // PERF could do this on the flat vector, but we only need to
      // run it for verbose rounds
      std::vector<std::vector<float>> expected;
      expected.reserve(EXAMPLES_PER_ROUND);
      for (int i = 0; i < EXAMPLES_PER_ROUND; i++) {
        std::vector<float> one;
        one.resize(OUTPUT_SIZE, 0.0f);
        for (int j = 0; j < OUTPUT_SIZE; j++) {
          one[j] = expecteds[i * OUTPUT_SIZE + j];
        }
        expected.emplace_back(std::move(one));
      }

      if (VERBOSE > 1)
        printf("Got expected\n");

      // Get loss as abs distance, plus number of incorrect (as booleans).
      std::vector<float> actuals(EXAMPLES_PER_ROUND * OUTPUT_SIZE);
      training->ExportOutputs(&actuals);

      const bool save_history = history_per.ShouldRun();
      const bool save_error_images = error_images_per.ShouldRun();
      // compute loss for each digit separately

      printf("%lld rounds [%.2f eps]:\n",
             net->rounds,
             (iter * EXAMPLES_PER_ROUND) / train_timer.Seconds());

      for (int col = 0; col < NUM_COLUMNS; col++) {

        std::vector<std::pair<float, float>> ex(EXAMPLES_PER_ROUND);
        float average_loss = 0, min_loss = 1.0f / 0.0f, max_loss = 0.0f;
        for (int idx = 0; idx < EXAMPLES_PER_ROUND; idx++) {
          float expected = expecteds[idx * NUM_COLUMNS + col];
          float actual = actuals[idx * NUM_COLUMNS + col];
          ex[idx] = std::make_pair(expected, actual);
          float loss = fabsf(expected - actual);
          min_loss = std::min(min_loss, loss);
          max_loss = std::max(max_loss, loss);
          average_loss += loss;
        }
        average_loss /= EXAMPLES_PER_ROUND;

        if (save_history) {
          error_history.Add(net->rounds, average_loss, col);
        }

        error_images[col]->Add(std::move(ex));

        if (save_error_images) {
          error_images[col]->Save();
          printf("Wrote %s\n", error_images[col]->Filename().c_str());
        }

        printf("   %c: %.3f<%.3f<%.3f\n",
               '0' + col,
               min_loss, average_loss, max_loss);
      }

      total_verbose++;
      loss_ms += loss_timer.MS();
    }

    // TODO: Should probably synchronize saving images with saving
    // the model. Otherwise when we continue, we lose pixels that
    // were written to the image but not saved to disk.
    if ((iter % IMAGE_EVERY) == 0) {
      Timer image_timer;
      net_gpu->ReadFromGPU();
      images.Sample(*net);
      image_ms += image_timer.MS();
    }

    static constexpr double SAVE_EVERY_SEC = 180.0;
    bool save_timeout = false;
    if ((train_timer.MS() / 1000.0) > last_save + SAVE_EVERY_SEC) {
      save_timeout = true;
      last_save = train_timer.MS() / 1000.0;
    }

    // Checkpoint (no overwrite) every X rounds.
    bool checkpoint_timeout = (net->rounds % CHECKPOINT_EVERY_ROUNDS) == 0;

    if (SAVE_INTERMEDIATE && (save_timeout || checkpoint_timeout || finished ||
                              iter % 5000 == 0)) {
      net_gpu->ReadFromGPU();
      // Note that if we write a checkpoint, we don't also overwrite
      // the main model, which is less redundant but might be surprising?
      const string filename = checkpoint_timeout ?
        StringPrintf(MODEL_BASE ".%lld.val", net->rounds) :
        (string)MODEL_NAME;
      net->SaveToFile(filename);
      if (VERBOSE)
        printf("Wrote %s\n", filename.c_str());
      error_history.Save();
    }

    // Parameter for average_loss termination?
    if (finished) {
      printf("Successfully trained!\n");
      return;
    }

    round_ms += round_timer.MS();

    if (timing_per.ShouldRun()) {
      double accounted_ms = getexample_ms + example_ms + forward_ms +
        error_ms + decay_ms + backward_ms + update_ms + loss_ms +
        image_ms;
      double other_ms = round_ms - accounted_ms;
      double pct = 100.0 / round_ms;
      printf("%.1f%% ge  "
             "%.1f%% ex  "
             "%.1f%% fwd  "
             "%.1f%% err  "
             "%.1f%% dec  "
             "%.1f%% bwd  "
             "%.1f%% up  "
             "%.1f%% loss "
             "%.1f%% img "
             "%.1f%% other\n",
             getexample_ms * pct,
             example_ms * pct,
             forward_ms * pct,
             error_ms * pct,
             decay_ms * pct,
             backward_ms * pct,
             update_ms * pct,
             loss_ms * pct,
             image_ms * pct,
             other_ms * pct);
      double msr = 1.0 / (iter + 1);
      printf("%.1fms ge  "
             "%.1fms ex  "
             "%.1fms fwd  "
             "%.1fms err  "
             "%.1fms dec  "
             "%.1fms bwd  "
             "%.1fms up  "
             "%.2fms loss  "
             "%.2fms img  "
             "%.1fms other\n",
             getexample_ms * msr,
             example_ms * msr,
             forward_ms * msr,
             error_ms * msr,
             decay_ms * msr,
             backward_ms * msr,
             update_ms * msr,
             loss_ms * msr,
             image_ms * msr,
             other_ms * msr);
    }
  }
}

static unique_ptr<Network> NewDigitsNetwork() {
  // Deterministic!
  ArcFour rc("learn-digits-network");

  std::vector<Layer> layers;

  Chunk input_chunk;
  input_chunk.type = CHUNK_INPUT;
  input_chunk.num_nodes = INPUT_SIZE;
  input_chunk.width = IMG_WIDTH;
  input_chunk.height = IMG_HEIGHT;
  input_chunk.channels = 1;

  layers.push_back(Network::LayerFromChunks(input_chunk));

  const int CONV1_SIZE = 3;
  const int CONV1_FEATURES = 32;

  Chunk first_conv_chunk =
    Network::Make2DConvolutionChunk(
        // Entire image
        0, IMG_WIDTH, IMG_HEIGHT,
        CONV1_FEATURES, CONV1_SIZE, CONV1_SIZE,
        // full overlap
        1, 1,
        LEAKY_RELU, WEIGHT_UPDATE);

  const int CONV2_SIZE = 8;
  const int CONV2_FEATURES = 32;

  Chunk second_conv_chunk =
    Network::Make2DConvolutionChunk(
        // Entire image
        0, IMG_WIDTH, IMG_HEIGHT,
        CONV2_FEATURES, CONV2_SIZE, CONV2_SIZE,
        // full overlap
        1, 1,
        LEAKY_RELU, WEIGHT_UPDATE);

  int SPARSE1_SIZE = 256;
  Chunk sparse_chunk =
    Network::MakeRandomSparseChunk(
        &rc, SPARSE1_SIZE,
        {{.span_start = 0, .span_size = INPUT_SIZE, .ipn = 64}},
        LEAKY_RELU, WEIGHT_UPDATE);

  layers.push_back(Network::LayerFromChunks(first_conv_chunk,
                                            second_conv_chunk,
                                            sparse_chunk));

  int num_conv1 =
    CONV1_FEATURES * first_conv_chunk.num_occurrences_across *
    first_conv_chunk.num_occurrences_down;
  Chunk next_conv1 =
    Network::Make1DConvolutionChunk(
        0, num_conv1,
        32, CONV1_FEATURES,
        // process the features as a block,
        CONV1_FEATURES,
        LEAKY_RELU, WEIGHT_UPDATE);

  int num_conv2 =
    CONV2_FEATURES * second_conv_chunk.num_occurrences_across *
    second_conv_chunk.num_occurrences_down;
  Chunk next_conv2 =
    Network::Make1DConvolutionChunk(
        num_conv1, num_conv2,
        32, CONV2_FEATURES,
        CONV2_FEATURES,
        LEAKY_RELU, WEIGHT_UPDATE);

  Chunk sparse2 =
    Network::MakeRandomSparseChunk(
        &rc, 256, {{.span_start = num_conv1 + num_conv2,
                    .span_size = SPARSE1_SIZE, .ipn = 64}},
        LEAKY_RELU, WEIGHT_UPDATE);

  layers.push_back(Network::LayerFromChunks(next_conv1,
                                            next_conv2,
                                            sparse2));

  Chunk sparse3 =
    Network::MakeRandomSparseChunk(
        &rc, 512, {{.span_start = 0, .span_size = layers.back().num_nodes,
                    .ipn = 64}},
        LEAKY_RELU, WEIGHT_UPDATE);

  layers.push_back(Network::LayerFromChunks(sparse3));

  // more layers?

  // Dense?
  Chunk sparse_out =
    Network::MakeRandomSparseChunk(
        &rc, OUTPUT_SIZE,
        {{.span_start = 0,
              .span_size = layers.back().num_nodes,
          .ipn = 64}},
        LEAKY_RELU, WEIGHT_UPDATE);

  layers.push_back(Network::LayerFromChunks(sparse_out));

  CHECK(layers.back().num_nodes > 0);

  auto net = std::make_unique<Network>(layers);

  printf("Randomize..\n");
  RandomizeNetwork(&rc, net.get(), 2);
  printf("New network with %lld parameters\n", net->TotalParameters());
  return net;
}

int main(int argc, char **argv) {
  cl = new CL;

  std::unique_ptr<Network> net(
      Network::ReadFromFile(MODEL_NAME));

  if (net.get() == nullptr) {
    net = NewDigitsNetwork();
    CHECK(net.get() != nullptr);
    net->SaveToFile(MODEL_NAME);
    printf("Wrote to %s\n", MODEL_NAME);
  }

  net->StructuralCheck();
  net->NaNCheck(MODEL_NAME);

  Train(net.get());

  printf("OK\n");
  return 0;
}

