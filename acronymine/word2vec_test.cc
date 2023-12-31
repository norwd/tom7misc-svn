
#include "word2vec.h"

#include <cmath>
#include <cstdint>
#include <memory>
#include <vector>
#include <string>
#include <unordered_set>

#include "base/logging.h"
#include "image.h"

#define DATA_FILE "c:\\code\\word2vec\\GoogleNews-vectors-negative300.bin"
#define FILL_FILE "word2vecfill.txt"

using namespace std;

static void CheckSimple(bool use_filter) {
  std::unordered_set<string> filter;
  if (use_filter) {
    filter.insert("miles");
    filter.insert("kilometers");
    filter.insert("sadness");
    filter.insert("melancholy");
    filter.insert("test");
    filter.insert("color");
    filter.insert("colour");
  }

  std::unique_ptr<Word2Vec> wv(
      Word2Vec::Load(DATA_FILE, FILL_FILE, filter));

  int miles = wv->GetWord("miles");
  int kilometers = wv->GetWord("kilometers");
  int sadness = wv->GetWord("sadness");
  int melancholy = wv->GetWord("melancholy");
  int color = wv->GetWord("color");
  int colour = wv->GetWord("colour");
  CHECK(miles != -1);
  CHECK(kilometers != -1);
  CHECK(sadness != -1);
  CHECK(melancholy != -1);
  CHECK(color != -1);
  CHECK(colour != -1) << "This word comes from fill file.";

  CHECK(miles != kilometers);
  CHECK(miles != sadness);
  CHECK(miles != melancholy);

  // Should not be exactly the same.
  CHECK(color != colour);

  {
    std::vector<int> testwords =
      {miles, kilometers, sadness, melancholy, color, colour};
    ImageRGBA img(wv->Size(), testwords.size());
    img.Clear32(0x000000FF);

    int y = 0;
    for (int word : testwords) {
      std::vector<float> vec = wv->NormVec(word);
      for (int x = 0; x < vec.size(); x++) {
        float f = vec[x];
        uint8_t r = f < 0.0 ? sqrtf(sqrtf(-f)) * 255.0f : 0.0;
        uint8_t g = f > 0.0 ? sqrtf(sqrtf(f)) * 255.0f : 0.0;
        img.SetPixel(x, y, r, g, 0, 0xFF);
      }
      y++;
    }

    img.ScaleBy(4).Save("word2vectest.png");
  }

  for (int word : {miles, kilometers, sadness, melancholy}) {
    std::vector<float> vec = wv->NormVec(word);
    printf("%d:", word);
    for (float f : vec) printf(" %.3f", f);
    printf("\n");
  }

  const float distdist = wv->Similarity(miles, kilometers);
  const float saddist = wv->Similarity(sadness, melancholy);

  const float ms = wv->Similarity(miles, sadness);
  const float mm = wv->Similarity(miles, melancholy);

  const float ks = wv->Similarity(kilometers, sadness);
  const float km = wv->Similarity(kilometers, melancholy);

  const float usuk = wv->Similarity(color, colour);

  printf("distdist %.11g\n"
         "saddist  %.11g\n"
         "ms mm    %.11g %.11g\n"
         "ks km    %.11g %.11g\n"
         "usuk     %.11g\n",
         distdist, saddist,
         ms, mm,
         ks, km,
         usuk);

  // These should be similar (i.e. positive)
  CHECK(distdist > 0.0f) << distdist;
  CHECK(saddist > 0.0f) << saddist;

  // Check that distances between similar words is smaller
  // than distances between dissimilar ones.
  CHECK(distdist > ms) << distdist << " " << ms;
  CHECK(distdist > mm) << distdist << " " << mm;

  CHECK(distdist > ks) << distdist << " " << ks;
  CHECK(distdist > km) << distdist << " " << km;

  CHECK(saddist > ms) << saddist << " " << ms;
  CHECK(saddist > mm) << saddist << " " << mm;

  CHECK(saddist > ks) << saddist << " " << ks;
  CHECK(saddist > km) << saddist << " " << km;
}

int main(int argc, char **argv) {
  CheckSimple(true );
  return 0;
}
