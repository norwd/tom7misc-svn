#ifndef _ACRONYMY_FREQ_H
#define _ACRONYMY_FREQ_H

#include <cstdint>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

struct Freq {
  static Freq *Load(const std::string &filename,
                    const std::unordered_set<std::string> &domain);

  // Returns 0 for a word not in the list.
  // Otherwise a word has a count of at least 1, even if it was
  // never actually seen. (Laplace smoothing)
  int64_t RawFreq(const std::string &word) const {
    auto it = smoothed_counts.find(word);
    if (it == smoothed_counts.end()) return 0;
    else return it->second;
  }

  // A probability distribution over all words in the domain.
  // Returns zero for words outside the domain.
  double Probability(const std::string &word) const {
    return RawFreq(word) * probability_scale;
  }

  // Where the most common word is 1.0.
  double NormalizedFreq(const std::string &word) const {
    return RawFreq(word) * norm_scale;
  }

  // With the most frequent words first.
  std::vector<std::string> SortedWords() const { return sorted; }

  int64_t TotalCount() const { return total_count; }
  // For sample_idx < total_count, returns the weighted sample
  // at that position. (e.g. ByWeightedSampleIndex(RandTo(TotalCount()))
  // gives you a random word according to the frequency distribution.)
  std::string ByWeightedSampleIndex(int64_t sample_idx) const;

 private:
  // Use Load.
  Freq();
  // These two already account for Laplace smoothing.
  std::unordered_map<std::string, int64_t> smoothed_counts;
  // Already accounting for Laplace smoothing.
  int64_t total_count = 0;
  // 1 / total_count
  double probability_scale = 1.0;
  // 1 / max_count
  double norm_scale = 1.0;

  // Parallel, sorted by descending frequency.
  std::vector<std::string> sorted;
  std::vector<int64_t> counts;
};

#endif
