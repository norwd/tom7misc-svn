
#include "modelinfo.h"

#include <memory>
#include <string>
#include <cstdint>
#include <cmath>
#include <time.h>

#include "network.h"
#include "base/logging.h"
#include "base/stringprintf.h"

using namespace std;

int main(int argc, char **argv) {
  const string modelfile = argc > 1 ? (string)argv[1] : "net0.val";

  // Try loading from disk; null on failure.
  printf("Load network from %s...\n", modelfile.c_str());
  std::unique_ptr<Network> net(Network::ReadFromFile(modelfile));

  CHECK(net.get() != nullptr) << modelfile;

  // Get histogram of weights per layer.
  // Only layers with inputs (i.e. not the input layer) have weights.
  const int num_histos = net->layers.size() - 1;

  static constexpr int HISTOW = 800;
  static constexpr int HISTOH = 220;
  static constexpr int MARGIN = 4;

  const int WIDTH = HISTOW * 2 + MARGIN;
  const int HEIGHT = HISTOH * num_histos;

  const ImageRGBA histos = ModelInfo::Histogram(*net, WIDTH, HEIGHT,
                                                // {-0.0000001f},
                                                // {+0.0000001f},
                                                nullopt,
                                                nullopt,
                                                nullopt,
                                                nullopt);

  char dates[128] = {};
  time_t tt = time(nullptr);
  // XXX even though TZ is set and 'date +"%H:%M"' works as expected
  // in mingw, this reports UTC time?
  tt -= 3600 * 5;

  strftime(dates, 127, "%d %b %Y  %H:%M", localtime(&tt));
  vector<string> lines = {
    StringPrintf("%s  round %lld   examples %lld   bytes %lld   layers %d",
                 dates,
                 net->rounds, net->examples, net->Bytes(), net->layers.size()),
    StringPrintf("                                      %lld total params",
                 net->TotalParameters()),
  };

  for (int layer_idx = 0; layer_idx < net->layers.size(); layer_idx++) {
    const Layer &layer = net->layers[layer_idx];
    string line = StringPrintf("%d: %d nodes; ", layer_idx, layer.num_nodes);
    for (int chunk_idx = 0; chunk_idx < layer.chunks.size(); chunk_idx++) {
      const Chunk &chunk = layer.chunks[chunk_idx];
      StringAppendF(&line, "%sx%d ",
                    ChunkTypeName(chunk.type), chunk.num_nodes);

      /*
      string types =
        layer.type == LAYER_DENSE ? (string)"DENSE" :
        layer.type == LAYER_SPARSE ? (string)"SPARSE" :
        layer.type == LAYER_CONVOLUTION_ARRAY ?
        StringPrintf("CONVx%d from %dx%d pat %dx%d +%dx%d",
                     layer.num_features,
                     layer.src_width, layer.src_height,
                     layer.pattern_width, layer.pattern_height,
                     layer.occurrence_x_stride, layer.occurrence_y_stride) :
        "???";
      */
    }
    lines.push_back(std::move(line));
  }

  const int TOP = 20 * lines.size();
  ImageRGBA img(histos.Width(), histos.Height() + TOP);
  img.Clear32(0x000000FF);
  img.BlendImage(0, TOP, histos);
  for (int i = 0; i < lines.size(); i++) {
    img.BlendText2x32(0, i * 20, 0xCCCCCCFF, lines[i]);
  }

  const string outfile = argc > 2 ? (string)argv[2] : "modelinfo.png";
  img.Save(outfile);

  printf("Wrote %s.\n", outfile.c_str());
  return 0;
}
