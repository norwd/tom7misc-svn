#include "network-gpu.h"

#include <cmath>
#include <memory>
#include <vector>
#include <functional>
#include <string>
#include <ctype.h>

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
#include "wikipedia.h"
#include "error-history.h"

using namespace std;

using TestNet = NetworkTestUtil::TestNet;
using TrainNet = NetworkTestUtil::TrainNet;
using TestExample = NetworkTestUtil::TestExample;

static CL *cl = nullptr;

using int64 = int64_t;

constexpr int MAX_WORD_LEN = 18;
constexpr int RADIX = 27;

static bool AllLetters(const string &s) {
  for (int i = 0; i < s.size(); i++) {
    if (s[i] < 'a' || s[i] > 'z') return false;
  }
  return true;
}

struct Wikibits {
  static constexpr int NUM_SHARDS = 128;

  Wikibits() : rc("wikibits" + StringPrintf("%lld", time(nullptr))) {
    std::vector<string> filenames;
    for (int i = 0; i < NUM_SHARDS; i++)
      filenames.push_back(StringPrintf("wikibits/wiki-%d.txt", i));
    std::vector<std::unordered_map<string, int64>> countvec =
      ParallelMap(filenames,
                  [](const std::string &filename) {
                    printf("Reading %s...\n", filename.c_str());
                    std::unordered_map<string, int64> counts;
                    int64 not_letters = 0;
                    auto co = Util::ReadFileOpt(filename);
                    CHECK(co.has_value());
                    string contents = std::move(co.value());
                    int64 pos = 0;
                    auto ReadByte = [&filename, &contents, &pos]() ->
                      uint8_t {
                        CHECK (pos < contents.size())
                          << filename << " @ " << pos << " of "
                          << contents.size();
                        return (uint8_t)contents[pos++];
                      };

                    auto Read32 = [&ReadByte]() {
                        uint32_t a = ReadByte();
                        uint32_t b = ReadByte();
                        uint32_t c = ReadByte();
                        uint32_t d = ReadByte();
                        return (a << 24) | (b << 16) | (c << 8) | d;
                      };

                    auto ReadString = [&contents, &pos, &ReadByte](int len) {
                        // printf("Read string of length %d:\n", len);
                        string s;
                        s.reserve(len);
                        for (int i = 0; i < len; i++) {
                          CHECK(pos < contents.size())
                            << i << "/" << len;
                          char c = ReadByte();
                          // printf("[%c]", c);
                          s.push_back(c);
                        }
                        return s;
                      };

                    int64 num_articles = 0;
                    auto Write = [&](uint32_t x) {
                      return StringPrintf("(%lld) %d [%c%c%c%c]",
                                          num_articles,
                                          x,
                                          (x >> 24) & 0xFF,
                                          (x >> 16) & 0xFF,
                                          (x >>  8) & 0xFF,
                                          x         & 0xFF);
                      };

                    while (pos < contents.size()) {
                      Wikipedia::Article art;
                      uint32_t t_len = Read32();
                      CHECK(t_len < 1000) << Write(t_len);
                      art.title = ReadString(t_len);
                      uint32_t b_len = Read32();
                      CHECK(b_len < 10000000) << Write(b_len);
                      art.body = ReadString(b_len);

                      for (int i = 0; i < art.title.size(); i++) {
                        CHECK(art.title[i] != 0);
                      }

                      for (int i = 0; i < art.body.size(); i++) {
                        CHECK(art.body[i] != 0);
                      }

                      // but convert article to words.
                      std::vector<string> tokens =
                        Util::Tokens(art.body,
                                     [](char c) { return isspace(c); });
                      for (const string &token : tokens) {
                        const string ltoken = Util::lcase(token);
                        if (AllLetters(token)) {
                          counts[ltoken]++;
                        } else {
                          not_letters++;
                        }
                      }
                      num_articles++;
                    }

                    printf("Distinct words: %lld, Not letters: %lld\n",
                           counts.size(), not_letters);
                    return counts;
                  }, 8);

    printf("Now build all map:\n");
    for (const auto &m : countvec) {
      for (const auto &[s, c] : m) {
        counts[s] += c;
      }
    }

    printf("Done. All distinct words: %lld\n", counts.size());
    distinct.reserve(counts.size());
    for (const auto &[s, c] : counts) {
      distinct.emplace_back(s, c);
    }
  }

  // XXX not thread-safe
  string RandomWord() {
    // XXX This should be weighted by frequency!
    int64 pos = RandTo(&rc, distinct.size());
    return distinct[pos].first;
  }

  std::unordered_map<string, int64> counts;
  std::vector<std::pair<string, int64>> distinct;
  ArcFour rc;
};


static void Train(Network *net) {

  Wikibits wikibits;

  ErrorHistory error_history("learn-words-error-history.tsv");

  // 0, 1, 2
  static constexpr int max_parallelism = 4;
  static constexpr int VERBOSE = 1;
  static constexpr bool SAVE_INTERMEDIATE = true;
  // Very small examples; could easily do 100x this...
  static constexpr int EXAMPLES_PER_ROUND = 1000;
  // XXX need to reduce this over time
  static constexpr float EXAMPLE_LEARNING_RATE =
    0.001f / (float)EXAMPLES_PER_ROUND;

  // XXX!
  std::vector<std::unique_ptr<ImageRGBA>> images;
  constexpr int IMAGE_WIDTH = 3000;
  constexpr int IMAGE_HEIGHT = 1000;
  constexpr int IMAGE_EVERY = 5;
  int image_x = 0;
  for (int i = 0; i < net->layers.size(); i++) {
    // XXX skip input layer?
    images.emplace_back(new ImageRGBA(IMAGE_WIDTH, IMAGE_HEIGHT));
    images.back()->Clear32(0x000000FF);
    images.back()->BlendText2x32(
        2, 2, 0x9999AAFF,
        StringPrintf("Layer %d | Start Round %lld | 1 px = %d rounds ",
                     i, net->rounds, IMAGE_EVERY));
  }

  printf("Training!\n");
  ArcFour rc("training");
  RandomGaussian gauss(&rc);

  auto net_gpu = make_unique<NetworkGPU>(cl, net);

  std::unique_ptr<ForwardLayerCL> forward_cl =
    std::make_unique<ForwardLayerCL>(cl, *net);
  std::unique_ptr<SetOutputErrorCL> error_cl =
    std::make_unique<SetOutputErrorCL>(cl, *net);
  std::unique_ptr<BackwardLayerCL> backward_cl =
    std::make_unique<BackwardLayerCL>(cl, *net);
  [[maybe_unused]]
  std::unique_ptr<DecayWeightsCL> decay_cl =
    std::make_unique<DecayWeightsCL>(cl, *net, 0.99999f);
  std::unique_ptr<UpdateWeightsCL> update_cl =
    std::make_unique<UpdateWeightsCL>(cl, *net);

  // Uninitialized training examples on GPU.
  std::vector<std::unique_ptr<TrainingRoundGPU>> training;
  for (int i = 0; i < EXAMPLES_PER_ROUND; i++)
    training.emplace_back(new TrainingRoundGPU(cl, *net));

  // Used to compute loss.
  std::vector<std::vector<float>> expected;
  expected.resize(training.size());

  Timer train_timer;
  int64 total_examples = 0LL;
  // seconds since timer started
  double last_save = 0.0;
  for (int iter = 0; true; iter++) {

    // Initialize training examples.
    // (PERF: parallelize?)
    for (int i = 0; i < training.size(); i++) {
      // One-hot
      std::vector<float> inputs(MAX_WORD_LEN * RADIX, 0.0f);

      string word = wikibits.RandomWord();
      for (int j = 0; j < MAX_WORD_LEN; j++) {
        const int c = j < word.size() ? word[j] - 'a' + 1 : 0;
        inputs[j * RADIX + c] = 1.0f;
      }
      training[i]->LoadInput(inputs);
      training[i]->LoadExpected(inputs);
      expected[i] = std::move(inputs);
    }

    if (VERBOSE > 1)
      printf("Prepped examples.\n");

    for (int src_layer = 0;
         src_layer < net->layers.size() - 1;
         src_layer++) {
      ParallelComp(
          training.size(),
          [&](int idx) {
            forward_cl->RunForward(
                net_gpu.get(), training[idx].get(), src_layer);
          },
          max_parallelism);
    }

    if (VERBOSE > 1)
      printf("Forward done.\n");

    ParallelComp(
        training.size(),
        [&](int idx) {
          error_cl->SetOutputError(net_gpu.get(), training[idx].get());
        },
        max_parallelism);

    if (VERBOSE > 1)
      printf("Set error.\n");

    for (int dst_layer = net->layers.size() - 1;
         // Don't propagate to input.
         dst_layer > 1;
         dst_layer--) {
      ParallelComp(
          training.size(),
          [&](int idx) {
            backward_cl->BackwardLayer(net_gpu.get(),
                                       training[idx].get(),
                                       dst_layer);
          },
          max_parallelism);
    }

    if (VERBOSE > 1)
      printf("Backward pass.\n");

    for (int layer_idx = 0; layer_idx < net->layers.size(); layer_idx++) {
      decay_cl->Decay(net_gpu.get(), layer_idx);
    }

    // Can't run training examples in parallel because these all write
    // to the same network. But each layer is independent.
    ParallelComp(net->layers.size() - 1,
                 [&](int layer_minus_1) {
                   const int layer_idx = layer_minus_1 + 1;
                   for (int i = 0; i < training.size(); i++) {
                     update_cl->Update(net_gpu.get(), training[i].get(),
                                       EXAMPLE_LEARNING_RATE, layer_idx);
                   }
                 },
                 max_parallelism);

    if (VERBOSE > 1)
      printf("Updated errors.\n");

    // Get loss as abs distance, plus number of incorrect (as booleans).
    // Size of examples = Number of training instances.
    string example_correct, example_predicted;
    std::vector<std::pair<float, int>> losses =
      ParallelMapi(expected,
                   [&](int idx, const std::vector<float> &exp) {
                     std::vector<float> got;
                     got.resize(exp.size());
                     training[idx]->ExportOutput(&got);

                     float loss = 0.0f;
                     for (int i = 0; i < exp.size(); i++) {
                       loss += fabsf(exp[i] - got[i]);
                     }

                     auto MaxChar = [](const std::vector<float> &v) {
                         CHECK(v.size() == MAX_WORD_LEN * RADIX);
                         string s;
                         s.reserve(MAX_WORD_LEN);
                         for (int c = 0; c < MAX_WORD_LEN; c++) {
                           int maxi = 0;
                           float maxv = -999999.0f;
                           for (int a = 0; a < RADIX; a++) {
                             int idx = c * RADIX + a;
                             if (v[idx] > maxv) {
                               maxi = a;
                               maxv = v[idx];
                             }
                           }
                           if (maxi == 0) s.push_back('_');
                           else s.push_back('a' + (maxi - 1));
                         }
                         return s;
                       };
                     string sexp = MaxChar(exp);
                     string sgot = MaxChar(got);

                     if (idx == 0) {
                       // careful about thread safety
                       example_correct = sexp;
                       example_predicted = sgot;
                     }

                     CHECK(sexp.size() == sgot.size());
                     int incorrect = 0;
                     for (int i = 0; i < sexp.size(); i++)
                       if (sexp[i] != sgot[i]) incorrect++;

                     return std::make_pair(loss, incorrect);
                   }, max_parallelism);

    if (VERBOSE > 1)
      printf("Got losses.\n");

    float min_loss = 1.0f / 0.0f, average_loss = 0.0f, max_loss = 0.0f;
    int min_inc = net->layers.back().num_nodes + 1, max_inc = 0;
    float average_inc = 0.0f;
    for (auto [loss_dist, loss_incorrect] : losses) {
      min_loss = std::min(loss_dist, min_loss);
      max_loss = std::max(loss_dist, max_loss);
      average_loss += loss_dist;

      min_inc = std::min(loss_incorrect, min_inc);
      max_inc = std::max(loss_incorrect, max_inc);
      average_inc += loss_incorrect;
    }
    average_loss /= losses.size();
    average_inc /= losses.size();

    total_examples += EXAMPLES_PER_ROUND;
    const double total_sec = train_timer.MS() / 1000.0;
    const double eps = total_examples / total_sec;

    constexpr int HISTORY_EVERY = 5;
    if ((iter % HISTORY_EVERY) == 0) {
      error_history.Add(net->rounds, average_loss, false);
    }

    net->examples += EXAMPLES_PER_ROUND;
    net->rounds++;

    if ((iter % IMAGE_EVERY) == 0) {

      // XXX would be better if this was more accurate,
      // but we only want to read from GPU if we're going to
      // actually do anything below
      if (images.size() > 2 &&
          images[1].get() != nullptr &&
          image_x < images[1]->Width()) {

        net_gpu->ReadFromGPU();

        for (int target_layer = 1; target_layer < net->layers.size();
             target_layer++) {
          ImageRGBA *image = images[target_layer].get();
          if (image == nullptr) continue;
          if (image_x >= image->Width()) continue;

          CHECK(net->layers.size() > 0);
          CHECK(target_layer < net->layers.size());
          const Layer &layer = net->layers[target_layer];
          CHECK(layer.chunks.size() > 0);
          const Chunk &chunk = layer.chunks[0];
          // x axis
          auto ToScreenY = [](float w) {
              int yrev = w * float(IMAGE_HEIGHT / 4) + (IMAGE_HEIGHT / 2);
              int y = IMAGE_HEIGHT - yrev;
              // Always draw on-screen.
              return std::clamp(y, 0, IMAGE_HEIGHT - 1);
            };
          if (image_x & 1) {
            image->BlendPixel32(image_x, ToScreenY(1), 0xCCFFCC40);
            image->BlendPixel32(image_x, ToScreenY(0), 0xCCCCFFFF);
            image->BlendPixel32(image_x, ToScreenY(-1), 0xFFCCCC40);
          }

          uint8 weight_alpha =
            std::clamp((255.0f / sqrtf(chunk.weights.size())), 10.0f, 240.0f);

          for (float w : chunk.weights) {
            // maybe better to AA this?
            image->BlendPixel32(image_x, ToScreenY(w),
                                0xFFFFFF00 | weight_alpha);
          }

          uint8 bias_alpha =
            std::clamp((255.0f / sqrtf(chunk.biases.size())), 10.0f, 240.0f);

          for (float b : chunk.biases) {
            image->BlendPixel32(image_x, ToScreenY(b),
                                0xFF777700 | bias_alpha);
          }

          if (chunk.weight_update == ADAM) {
            CHECK(chunk.weights_aux.size() == 2 * chunk.weights.size());
            CHECK(chunk.biases_aux.size() == 2 * chunk.biases.size());
            for (int idx = 0; idx < chunk.weights.size(); idx++) {
              const float m = chunk.weights_aux[idx * 2 + 0];
              const float v = sqrtf(chunk.weights_aux[idx * 2 + 1]);

              image->BlendPixel32(image_x, ToScreenY(m),
                                  0xFFFF0000 | weight_alpha);
              image->BlendPixel32(image_x, ToScreenY(v),
                                  0xFF00FF00 | weight_alpha);
            }
            // Also bias aux?
          }


          if ((image_x % 100 == 0) || image_x == image->Width()) {
            string filename = StringPrintf("train-image-%d.png",
                                           target_layer);
            image->Save(filename);
            printf("Wrote %s\n", filename.c_str());
          }
        }
        image_x++;
      }
    }

    const bool finished = max_inc == 0;

    static constexpr double SAVE_EVERY_SEC = 120.0;
    bool save_timeout = false;
    if ((train_timer.MS() / 1000.0) > last_save + SAVE_EVERY_SEC) {
      save_timeout = true;
      last_save = train_timer.MS() / 1000.0;
    }

    if (SAVE_INTERMEDIATE && (save_timeout || finished ||
                              iter == 1000 || iter % 5000 == 0)) {
      net_gpu->ReadFromGPU();
      const string file = StringPrintf("words.val", iter);
      net->SaveToFile(file);
      if (VERBOSE)
        printf("Wrote %s\n", file.c_str());
      error_history.Save();
    }

    // Parameter for average_loss termination?
    if (finished) {
      printf("Successfully trained!\n");
      return;
    } else {
      if (VERBOSE || (iter % 100 == 0)) {
        printf("%d: %.3f<%.3f<%.3f", iter,
               min_loss, average_loss, max_loss);
        printf(" | %d<%.3f<%d",
               min_inc, average_inc, max_inc);
        printf(" (%.2f eps)\n", eps);
        printf("   [%s] got [%s]\n",
               example_correct.c_str(), example_predicted.c_str());
      }
    }
  }


}


// Create an Nx1 convolutional chunk. It reads the entire previous
// layer of size prev_size.
static Chunk ConvolutionalChunk1D(int prev_size,
                                  int num_features,
                                  int x_stride,
                                  int pattern_width) {
  Chunk chunk;
  chunk.type = CHUNK_CONVOLUTION_ARRAY;
  chunk.num_features = num_features;
  chunk.occurrence_x_stride = x_stride;
  chunk.occurrence_y_stride = 1;
  chunk.pattern_width = pattern_width;
  chunk.pattern_height = 1;
  chunk.src_width = prev_size;
  chunk.src_height = 1;
  chunk.transfer_function = LEAKY_RELU;
  chunk.span_start = 0;
  chunk.span_size = prev_size;
  chunk.indices_per_node = pattern_width;

  {
    const auto [indices, this_num_nodes,
                num_occurrences_across, num_occurrences_down] =
      Network::MakeConvolutionArrayIndices(0, prev_size,
                                           chunk.num_features,
                                           chunk.pattern_width,
                                           chunk.pattern_height,
                                           chunk.src_width,
                                           chunk.src_height,
                                           chunk.occurrence_x_stride,
                                           chunk.occurrence_y_stride);
    chunk.num_nodes = this_num_nodes;
    chunk.width = chunk.num_nodes;
    chunk.height = 1;
    chunk.channels = 1;

    chunk.num_occurrences_across = num_occurrences_across;
    CHECK(num_occurrences_down == 1);
    chunk.num_occurrences_down = num_occurrences_down;
    chunk.indices = indices;

    chunk.weights = std::vector<float>(
        chunk.indices_per_node * chunk.num_features,
        0.0f);
    chunk.biases = std::vector<float>(chunk.num_features, 0.0f);
  }

  chunk.weight_update = ADAM;
  chunk.weights_aux.resize(chunk.weights.size() * 2, 0.0f);
  chunk.biases_aux.resize(chunk.biases.size() * 2, 0.0f);

  return chunk;
}

static Chunk DenseChunk(int prev_size,
                        int this_size) {
  Chunk chunk;
  chunk.type = CHUNK_DENSE;
  chunk.num_nodes = this_size;
  chunk.transfer_function = LEAKY_RELU;
  chunk.width = this_size;
  chunk.height = 1;
  chunk.channels = 1;
  chunk.span_start = 0;
  chunk.span_size = prev_size;
  chunk.indices_per_node = prev_size;
  chunk.indices = {};
  chunk.weights = std::vector<float>(
      chunk.num_nodes * chunk.indices_per_node, 0.0f);
  chunk.biases = std::vector<float>(chunk.num_nodes, 0.0f);

  chunk.weight_update = ADAM;
  chunk.weights_aux.resize(chunk.weights.size() * 2, 0.0f);
  chunk.biases_aux.resize(chunk.biases.size() * 2, 0.0f);

  return chunk;
}

static Network *NewNetwork() {
  auto L = [&](const Chunk &chunk) {
      return Layer{.num_nodes = chunk.num_nodes, .chunks = {chunk}};
    };

  constexpr int INPUT_SIZE = MAX_WORD_LEN * RADIX;
  Chunk input_chunk;
  input_chunk.type = CHUNK_INPUT;
  input_chunk.num_nodes = INPUT_SIZE;
  input_chunk.width = INPUT_SIZE;
  input_chunk.height = 1;
  input_chunk.channels = 1;

  // Convolve to 10 bits/char. We should only need 5...
  constexpr int BITS_PER_CHAR = 10;
  Chunk conv_chunk1 = ConvolutionalChunk1D(INPUT_SIZE,
                                           BITS_PER_CHAR,
                                           // non-overlapping
                                           RADIX, RADIX);

  // Then bigrams.
  // It would be useful to support overlap here. But I don't know
  // how to expand the thing afterwards, since we'd end up with
  // MAX_WORD_LEN - 1 bigrams in that case?

  // probably also more than we should need
  constexpr int BITS_PER_BIGRAM = 15;
  Chunk conv_chunk2 = ConvolutionalChunk1D(BITS_PER_CHAR * MAX_WORD_LEN,
                                           BITS_PER_BIGRAM,
                                           BITS_PER_CHAR * 2,
                                           BITS_PER_CHAR * 2);

  // Number of nodes in output of previous
  static_assert((MAX_WORD_LEN & 1) == 0);
  constexpr int PRE_ENCODED_WIDTH = (MAX_WORD_LEN / 2) * BITS_PER_BIGRAM;
  CHECK(conv_chunk2.num_nodes == PRE_ENCODED_WIDTH);

  constexpr int ENCODED = 64;
  Chunk dense_encode = DenseChunk(PRE_ENCODED_WIDTH, ENCODED);

  // And now the reverse.
  Chunk dense_decode = DenseChunk(ENCODED, PRE_ENCODED_WIDTH);

  Chunk unconv_chunk2 = ConvolutionalChunk1D(PRE_ENCODED_WIDTH,
                                             BITS_PER_CHAR * 2,
                                             // non-overlapping
                                             BITS_PER_BIGRAM,
                                             BITS_PER_BIGRAM);
  Chunk unconv_chunk1 = ConvolutionalChunk1D(BITS_PER_CHAR * MAX_WORD_LEN,
                                             RADIX,
                                             BITS_PER_CHAR,
                                             BITS_PER_CHAR);


  return new Network(vector<Layer>{
                         L(input_chunk),
                         L(conv_chunk1),
                         L(conv_chunk2),
                         L(dense_encode),
                         L(dense_decode),
                         L(unconv_chunk2),
                         L(unconv_chunk1)});

#if 0
  Chunk unconv_chunk1 = ConvolutionalChunk1D(BITS_PER_CHAR * MAX_WORD_LEN,
                                             RADIX,
                                             BITS_PER_CHAR,
                                             BITS_PER_CHAR);

  return new Network(vector<Layer>{
      L(input_chunk),
      L(conv_chunk1),
      // ... TODO: grow here ...
      L(unconv_chunk1)
        });
#endif
}


int main(int argc, char **argv) {
  cl = new CL;

  std::unique_ptr<Network> net(
      Network::ReadFromFile("words.val"));

  if (net.get() == nullptr) {
    net.reset(NewNetwork());
    net->StructuralCheck();
    ArcFour rc("new");
    RandomizeNetwork(&rc, net.get(), 2);
    printf("New network with %lld parameters\n", net->TotalParameters());
  }

  Train(net.get());

  printf("OK\n");
  return 0;
}
