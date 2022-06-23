
#ifndef _NETWORK_TEST_UTIL_H
#define _NETWORK_TEST_UTIL_H

#include <vector>
#include <functional>
#include <string>

#include "network.h"

// Note: Evaluates a and b a second time if check fails!
#define CHECK_FEQ(a, b) CHECK(fabs((a) - (b)) < 0.00001) \
  << #a " = " << (a) << " vs " #b " = " << (b)

// Same, but for vectors of floats.
#define CHECK_FEQV(aav, bbv) do {                           \
    auto av = (aav), bv = (bbv);                            \
    CHECK_EQ(av.size(), bv.size()) <<                       \
      av.size() << " vs " << bv.size();                     \
    for (size_t i = 0; i < av.size(); i++)                  \
      CHECK_FEQ(av[i], bv[i]) <<                            \
        "\n" #aav "[" << i << "] vs " #bbv "[" << i << "]"; \
  } while (0)

struct NetworkTestUtil {

  struct TestExample {
    std::string name;
    std::vector<float> input;
    std::vector<float> output;
  };

  struct TestNet {
    std::string name;
    Network net;
    std::vector<TestExample> examples;
  };

  // A network structure that should be able to learn the given
  // function. (And its input/output shape should match the function.)
  struct TrainNet {
    std::string name;
    // Weights and biases initialized to zero.
    Network net;
    std::function<std::vector<float>(const std::vector<float>&)> f;
    // If true, inputs should just be 0.0 or 1.0.
    bool boolean_input = false;
    // If true, outputs should be treated as false/true with a cutoff of 0.5.
    bool boolean_output = false;

    int NumInputs() const;
    int NumOutputs() const;
  };

  // Convert all chunks to ADAM weight_update.
  static TrainNet ForceAdam(TrainNet net);
  // Convert all chunks to YOGI weight_update.
  static TrainNet ForceYogi(TrainNet net);

  // Trivial network with just one node, sparse chunk.
  static TestNet SingleSparse();
  // Trivial network with just one node, dense chunk.
  static TestNet SingleDense();
  // Trivial network with one 1x1 "convolutional" chunk.
  static TestNet SingleConvolution();

  // Input size 2; output is a single node 2a + 3b + 1.
  static TestNet TwoInputSparse();

  // One input, two outputs in separate dense chunks.
  // computes 5a - 7a.
  static TestNet TwoDenseChunks();

  // 3x2 convolution on 4x3 input.
  // One feature:
  // 3  2 1
  // 4 -5 6
  //  -100
  static TestNet SimpleConv();

  // Small network that computes a nontrivial function.
  // This one also has rounds and examples set > 2^31.
  static TestNet Net1();

  // Copies input to output, using MakeCopyChunk.
  static TestNet Copy();

  // With a fixed chunk.
  static TestNet FixedSingle();

  // Nontrivial network with two dense layers (one chunk each).
  // input  a0, b0
  // hidden a1 = leaky(.5 + 2a0 - 3b0)   b1 = leaky(-1 - 1a0 + .25b0)
  // out    a2 = leaky(0 + 1.5a1 + 1b0)  b2 = leaky(1 - 0.25a1 + 3b1)
  static TestNet TwoDenseLayers();

  // From an input of size 8, counts internal 0-1 and 1-0 transitions,
  // with a convolution and a dense node. (Output is not exact.)
  static TestNet CountInternalEdges();

  // TODO: Need some tests with convolutions!

  // F(x) = x. One sparse node.
  static TrainNet LearnTrivialIdentitySparse();
  static TrainNet LearnTrivialIdentityDense();
  static TrainNet LearnTrivialIdentityConvolution();

  // All Boolean functions of three variables (2^(2^3) many), with two
  // dense layers.
  static TrainNet LearnBoolean();

  // TODO: More simple, fast-converging tests that have multiple chunks
  // or layers.

  // Count ones with a single dense layer.
  static TrainNet LearnCountOnesDense();

  // First a 4x1 convolution to sum, then dense over the 20 remaining.
  static TrainNet LearnCountOnesConvDense();

  // 80 through 4x1 -> 20 through 5x1 -> 4 through dense -> 1
  // (With the arg set to true, the dense layer is fixed (to the correct
  // answer), which is kinda "cheating".)
  static TrainNet LearnCountOnesConvConvDense(bool fix_dense_layer);

  // 2x1 convolution that should be able to learn to count 0-1 and 1-0
  // edges in a bit string.
  static TrainNet LearnCountEdges();

  // Uses Adam weight update.
  // --dim--
  // c a a a |
  // b c a a dim
  // b b c a |
  // b b b c |
  //
  // 1.0 if the sum of a is greater than b,
  // or sum(c) > sum(b) && sum(c) < sum(a). 0 otherwise.
  //
  // Last layer is 1 dense node (sigmoid). depth gives the number of
  // dim*dim dense hidden layers preceding that.
  static TrainNet TriangleSumsAdam(int depth, int dim = 4);

  // The input defines a random point (x, y) and a circle (cx, cy, radius),
  // and the output is 1.0 if inside the circle, otherwise 0.0. This
  // requires approximating the square (or square root).
  //
  // All hidden layers (depth) are dense, and of size width.
  static TrainNet InCircleAdam(int width, int depth);
  // Same but with sparse layers using the given indices per node.
  static TrainNet SparseInCircleAdam(int width, int ipn, int depth);

  // The input defines two random line segments (p0x, p0y, p1x, p1y)
  // and (p2x, p2y, p3x, p3y).
  // If the segments do not intersect, the output is (0.0, 0.0, 0.0).
  // Otherwise, the output is (1.0, x, y) where (x,y) is the point of
  // intersection.
  static TrainNet SparseLineIntersectionAdam(int width, int ipn, int depth,
                                             int64_t seed = 1);

  // The input defines a random line segment and a random point. The
  // output is the point reflected about that line.
  static TrainNet ReflectAdam(int width, int ipn, int depth);

  // The atan2 function of two variables.
  static TrainNet Atan2Adam(int width, int ipn, int depth);

  // Takes two circular objects. An object is (x, y, dx, dy, radius).
  // The arguments are scaled internally so that gaussian inputs produce
  // reasonably interesting results (with about 12% colliding, 5pp of
  // which because they are already overlapping in the initial state).
  // Model starts with convolutions, since the objects are symmetric.
  // Predict whether the two objects will collide within 100 timesteps,
  // and if so return the first colliding timestep/100.
  static TrainNet DodgeballAdam(int width, int ipn, int depth,
                                int64_t seed = 1);

  // Tests a very simple case for the tanh transfer function. The
  // sign of the sum of inputs (-1, 0, 1).
  static TrainNet TanhSignFunctionAdam();
};


#endif
