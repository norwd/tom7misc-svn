
#include "network-test-util.h"

#include "network.h"
#include "base/logging.h"

NetworkTestUtil::TestNet NetworkTestUtil::SingleSparse() {
  Chunk input_chunk;
  input_chunk.type = CHUNK_INPUT;
  input_chunk.num_nodes = 1;
  input_chunk.width = 1;
  input_chunk.height = 1;
  input_chunk.channels = 1;

  Chunk sparse_chunk;
  sparse_chunk.type = CHUNK_SPARSE;
  sparse_chunk.num_nodes = 1;
  sparse_chunk.transfer_function = IDENTITY;
  sparse_chunk.width = 1;
  sparse_chunk.height = 1;
  sparse_chunk.channels = 1;
  sparse_chunk.span_start = 0;
  sparse_chunk.span_size = 1;
  sparse_chunk.indices_per_node = 1;
  sparse_chunk.indices = {0};
  sparse_chunk.weights = {1.0};
  sparse_chunk.biases = {0.0};

  Layer input_layer;
  input_layer.num_nodes = 1;
  input_layer.chunks = {input_chunk};

  Layer real_layer;
  real_layer.num_nodes = 1;
  real_layer.chunks = {sparse_chunk};

  Network net({input_layer, real_layer});
  net.NaNCheck(__func__);

  CHECK(net.layers.size() == 2);
  CHECK(net.layers[0].chunks.size() == 1);
  CHECK(net.layers[1].chunks.size() == 1);

  TestExample example1{
    .name = "five",
    .input = {5.0},
    .output = {5.0},
  };

  return TestNet{
    .name = "one real layer with one sparse node, computes identity",
    .net = net,
    .examples = {example1},
  };
}

NetworkTestUtil::TestNet NetworkTestUtil::SingleDense() {
  Chunk input_chunk;
  input_chunk.type = CHUNK_INPUT;
  input_chunk.num_nodes = 1;
  input_chunk.width = 1;
  input_chunk.height = 1;
  input_chunk.channels = 1;

  Chunk dense_chunk;
  dense_chunk.type = CHUNK_DENSE;
  dense_chunk.num_nodes = 1;
  dense_chunk.transfer_function = IDENTITY;
  dense_chunk.width = 1;
  dense_chunk.height = 1;
  dense_chunk.channels = 1;
  dense_chunk.span_start = 0;
  dense_chunk.span_size = 1;
  dense_chunk.indices_per_node = 1;
  // indices not stored for dense chunks
  dense_chunk.indices = {};
  dense_chunk.weights = {1.0};
  dense_chunk.biases = {0.0};

  Layer input_layer;
  input_layer.num_nodes = 1;
  input_layer.chunks = {input_chunk};

  Layer real_layer;
  real_layer.num_nodes = 1;
  real_layer.chunks = {dense_chunk};

  Network net({input_layer, real_layer});
  net.NaNCheck(__func__);

  CHECK(net.layers.size() == 2);
  CHECK(net.layers[0].chunks.size() == 1);
  CHECK(net.layers[1].chunks.size() == 1);

  TestExample example1{
    .name = "seven",
    .input = {7.0},
    .output = {7.0},
  };

  return TestNet{
    .name = "one real layer with one dense node, computes identity",
    .net = net,
    .examples = {example1},
  };
}

NetworkTestUtil::TestNet NetworkTestUtil::TwoInputSparse() {
  Chunk input_chunk;
  input_chunk.type = CHUNK_INPUT;
  input_chunk.num_nodes = 2;
  input_chunk.width = 2;
  input_chunk.height = 1;
  input_chunk.channels = 1;

  Chunk sparse_chunk;
  sparse_chunk.type = CHUNK_SPARSE;
  sparse_chunk.num_nodes = 1;
  sparse_chunk.transfer_function = IDENTITY;
  sparse_chunk.width = 1;
  sparse_chunk.height = 1;
  sparse_chunk.channels = 1;
  sparse_chunk.span_start = 0;
  sparse_chunk.span_size = 2;
  sparse_chunk.indices_per_node = 2;
  sparse_chunk.indices = {0, 1};
  sparse_chunk.weights = {2.0, 3.0};
  sparse_chunk.biases = {1.0};

  Layer input_layer;
  input_layer.num_nodes = 2;
  input_layer.chunks = {input_chunk};

  Layer real_layer;
  real_layer.num_nodes = 1;
  real_layer.chunks = {sparse_chunk};

  Network net({input_layer, real_layer});
  net.NaNCheck(__func__);

  CHECK(net.layers.size() == 2);
  CHECK(net.layers[0].chunks.size() == 1);
  CHECK(net.layers[1].chunks.size() == 1);

  TestExample example1{
    .name = "eight-nine",
    .input = {8.0f, 9.0f},
    .output = {2.0f * 8.0f + 3.0f * 9.0f + 1.0f},
  };

  return TestNet{
    .name = "one node computing 2a + 3b + 1",
    .net = net,
    .examples = {example1},
  };
}


NetworkTestUtil::TestNet NetworkTestUtil::TwoDenseChunks() {
  Chunk input_chunk;
  input_chunk.type = CHUNK_INPUT;
  input_chunk.num_nodes = 1;
  input_chunk.width = 1;
  input_chunk.height = 1;
  input_chunk.channels = 1;

  Chunk dense_chunk1;
  dense_chunk1.type = CHUNK_DENSE;
  dense_chunk1.num_nodes = 1;
  dense_chunk1.transfer_function = IDENTITY;
  dense_chunk1.width = 1;
  dense_chunk1.height = 1;
  dense_chunk1.channels = 1;
  dense_chunk1.span_start = 0;
  dense_chunk1.span_size = 1;
  dense_chunk1.indices_per_node = 1;
  // indices not stored for dense chunks
  dense_chunk1.indices = {};
  dense_chunk1.weights = {5.0f};
  dense_chunk1.biases = {1.0f};

  Chunk dense_chunk2 = dense_chunk1;
  dense_chunk2.weights = {-7.0f};
  dense_chunk2.biases = {2.0f};

  Layer input_layer;
  input_layer.num_nodes = 1;
  input_layer.chunks = {input_chunk};

  Layer real_layer;
  real_layer.num_nodes = 2;
  real_layer.chunks = {dense_chunk1, dense_chunk2};

  Network net({input_layer, real_layer});
  net.NaNCheck(__func__);

  CHECK(net.layers.size() == 2);
  CHECK(net.layers[0].chunks.size() == 1);
  CHECK(net.layers[1].chunks.size() == 2);

  TestExample example1{
    .name = "three",
    .input = {3.0f},
    .output = {5.0f * 3.0f + 1.0f, -7.0f * 3.0f + 2.0f},
  };

  return TestNet{
    .name = "two dense chunks, computing (5a + 1, -7a + 2)",
    .net = net,
    .examples = {example1},
  };
}


NetworkTestUtil::TestNet NetworkTestUtil::Net1() {
  Chunk input_chunk;
  input_chunk.type = CHUNK_INPUT;
  input_chunk.num_nodes = 3;
  input_chunk.width = 3;
  input_chunk.height = 1;
  input_chunk.channels = 1;

  Chunk dense_chunk =
    Network::MakeDenseChunk(2,
                            // Span
                            0, 2,
                            //
                            IDENTITY);
  dense_chunk.weights = {2.0, 3.0, 4.0, 5.0};
  dense_chunk.biases = {-100.0, -200.0};

  Chunk sparse_chunk;
  sparse_chunk.type = CHUNK_SPARSE;
  sparse_chunk.num_nodes = 2;
  sparse_chunk.transfer_function = LEAKY_RELU;
  sparse_chunk.width = 2;
  sparse_chunk.height = 1;
  sparse_chunk.channels = 1;
  sparse_chunk.span_start = 1;
  sparse_chunk.span_size = 2;
  sparse_chunk.indices_per_node = 1;
  sparse_chunk.indices = {2, 1};
  sparse_chunk.weights = {10.0, 70.0};
  sparse_chunk.biases = {-1000.0, -2000.0};

  Layer input_layer;
  input_layer.num_nodes = 3;
  input_layer.chunks = {input_chunk};

  Layer real_layer;
  real_layer.num_nodes = 4;
  real_layer.chunks = {dense_chunk, sparse_chunk};

  Network net({input_layer, real_layer});
  net.NaNCheck(__func__);

  CHECK(net.layers.size() == 2);
  CHECK(net.layers[0].chunks.size() == 1);
  CHECK(net.layers[1].chunks.size() == 2);

  TestExample example1{
    .name = "one",
    .input = {3.0, 5.0, 7.0},
    .output = {
      // dense chunk (identity transfer function)
      -100.0 + 2.0 * 3.0 + 3.0 * 5.0,
      -200.0 + 4.0 * 3.0 + 5.0 * 5.0,
      // sparse chunk (leaky relu)
      0.01f * (-1000.0 + 10.0 * 7.0),
      0.01f * (-2000.0 + 70.0 * 5.0),
    },
  };

  return TestNet{
    .name = "one real layer, dense id and sparse leaky relu chunks",
    .net = net,
    .examples = {example1},
  };
}
