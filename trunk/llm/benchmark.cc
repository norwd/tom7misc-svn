
#include "llama.h"

#include <algorithm>
#include <cassert>
#include <cinttypes>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "base/logging.h"
#include "base/stringprintf.h"
#include "ansi.h"
#include "timer.h"
#include "util.h"
#include "vector-util.h"
#include "arcfour.h"
#include "randutil.h"

#include "llm.h"
#include "llm-util.h"

// Prevent runaway (no correct answer can be longer).
static constexpr int MAX_ANSWER = 80;
static constexpr int MAX_THOUGHT = 160;

static constexpr bool USE_THOUGHT = true;

using namespace std;

static inline std::string CleanWhitespace(std::string s) {
  return Util::LoseWhiteL(Util::LoseWhiteR(std::move(s)));
}

// TODO: Add progress bar.

// Question text; expected answer(s).
struct Question {
  string statement;
  // "Show your work" for the example problem.
  string thought;
  std::vector<string> answers;
};

struct Problem {
  // The problem name; a-z0-9_.
  string name;
  // Introductory prompt describing the problem.
  string prompt;
  // Unscored example input/output pair.
  Question example;
  std::vector<Question> questions;

  string question_before = "Question: ";
  string question_after = "\n";

  string answer_before = "Answer: [";
  string answer_after = "]\n";
};

// Consumes lines from the vector. Returns the first one that
// is not empty nor a comment (#), after trimming outer whitespace.
string GetNonCommentLine(const string &error_context,
                         std::vector<string> *lines) {
  for (;;) {
    CHECK(!lines->empty()) << error_context;
    string line = CleanWhitespace((*lines)[0]);
    lines->erase(lines->begin());
    if (line.empty()) continue;
    if (line[0] == '#') continue;
    return line;
  }
}

// Read the header of a problem file. This is currently the name and
// prompt. Consumes lines from the vector.
void ReadHeader(const string &error_context,
                std::vector<string> *lines, Problem *problem) {
  problem->name = GetNonCommentLine(error_context, lines);
  problem->prompt = GetNonCommentLine(error_context, lines);
}

Problem LoadFreeformProblem(const string &filename) {
  std::vector<string> lines = Util::ReadFileToLines(filename);
  for (string &line : lines) {
    line = CleanWhitespace(std::move(line));
  }

  FilterVector(&lines, [](const string &line) {
      if (line.empty()) return false;
      if (line[0] == '#') return false;
      return true;
    });

  Problem problem;
  ReadHeader(filename, &lines, &problem);

  for (int idx = 0; idx < (int)lines.size(); idx++) {
    const string &line = lines[idx];
    vector<string> v = Util::Split(line, '|');
    CHECK(v.size() >= 2) << filename << ": " << line;
    // Strip all surrounding whitespace.
    for (string &s : v) s = CleanWhitespace(s);
    Question question;
    question.statement = std::move(v[0]);
    v.erase(v.begin());

    // If the first "answer" is surrounded by [square brackets], it is
    // a thought (for examples).

    if (v[0][0] == '[') {
      string thought = std::move(v[0]);
      v.erase(v.begin());

      // Check and strip [brackets].
      CHECK(thought[thought.size() - 1] == ']') << filename << " : "
                                                << thought;
      thought = thought.substr(1, thought.size() - 2);
      CHECK(!thought.empty()) << filename;
      question.thought = thought;
    }

    // Answers are all the remaining fields
    question.answers = std::move(v);
    problem.questions.push_back(std::move(question));
  }
  CHECK(!problem.questions.empty());

  // The first question is the example.
  problem.example = std::move(problem.questions[0]);
  problem.questions.erase(problem.questions.begin());

  return problem;
}

// Format is like this
// LANG: C++
// THOUGHT: This is a loop that only terminates when the variable x is zero.
// OUTPUT: aaaz
// <code>
// for (int x = 3; x != 0; x--) {
//   printf("a");
// }
// printf("z");
// </code>
Problem LoadCodeProblem(const string &filename) {
  std::vector<string> lines = Util::ReadFileToLines(filename);

  Problem problem;
  ReadHeader(filename, &lines, &problem);

  // This is sort of a hack, but since the question itself has
  // "Language: ", we don't need a prefix.
  problem.question_before = "";

  // Make it clearer we're looking for the program's output, not
  // some other "answer".
  problem.answer_before = "Output: [";

  string lang;
  string thought;
  string output;
  string code;
  std::map<string, string> parts;
  auto EmitQuestion = [&filename,
                       &problem, &lang, &thought, &output, &code]() {
      // But thought is optional.
      CHECK(!lang.empty());
      CHECK(!output.empty());
      CHECK(!code.empty());

      Question question;
      question.statement =
        StringPrintf("Language: %s\n"
                     "<code>\n"
                     "%s"
                     "</code>\n", lang.c_str(), code.c_str());
      question.thought = std::move(thought);
      question.answers.push_back(std::move(output));
      problem.questions.emplace_back(std::move(question));

      lang.clear();
      output.clear();
      code.clear();
      thought.clear();
    };

  bool in_code = false;
  for (int i = 0; i < (int)lines.size(); i++) {
    string line = lines[i];

    if (in_code) {
      if (line == "</code>") {
        in_code = false;
        // Done with question.
        EmitQuestion();
      } else {
        code += line;
        code += "\n";
        // StringAppendF(&code, "%s\n", line.c_str());
      }
    } else {
      // Normalize whitespace only outside of code block.
      line = CleanWhitespace(std::move(line));
      if (line.empty()) continue;
      if (line[0] == '#') continue;

      if (Util::TryStripPrefix("<code>", &line)) {
        in_code = true;
      } else if (Util::TryStripPrefix("LANG: ", &line)) {
        lang = line;
      } else if (Util::TryStripPrefix("THOUGHT: ", &line)) {
        thought = line;
      } else if (Util::TryStripPrefix("OUTPUT: ", &line)) {
        output = line;
      } else {
        CHECK(false) << filename << ": Invalid line: " << line;
      }
    }
  }

  // The first question is the example.
  CHECK(!problem.questions.empty()) << filename;
  problem.example = std::move(problem.questions[0]);
  problem.questions.erase(problem.questions.begin());

  return problem;
}

Problem LoadMultipleChoiceProblem(const string &filename) {
  ArcFour rc(StringPrintf("%s.%lld", filename.c_str(), time(nullptr)));
  std::vector<string> lines = Util::ReadFileToLines(filename);
  for (string &line : lines) {
    line = CleanWhitespace(std::move(line));
  }

  FilterVector(&lines, [](const string &line) {
      if (line.empty()) return false;
      if (line[0] == '#') return false;
      return true;
    });

  Problem problem;
  ReadHeader(filename, &lines, &problem);

  // State of current question; we read line-by-line.
  string statement;
  string thought;
  // Exactly one choice should be marked correct.
  vector<std::pair<bool, string>> choices;

  auto EmitQuestion = [&filename, &rc,
                       &problem, &statement, &thought, &choices]() {
      CHECK(!statement.empty());
      CHECK(!choices.empty());
      Shuffle(&rc, &choices);

      int correct_idx = -1;
      string rendered = StringPrintf("%s\n", statement.c_str());
      for (int i = 0; i < (int)choices.size(); i++) {
        if (choices[i].first) correct_idx = i;
        StringAppendF(&rendered, " (%c) %s", 'a' + i,
                      choices[i].second.c_str());
        if (i < (int)choices.size() - 1) StringAppendF(&rendered, "\n");
      }
      CHECK(correct_idx >= 0) << filename
                              << "\nNo correct answer:\n"
                              << rendered;

      Question question;
      question.statement = rendered;
      question.thought = std::move(thought);
      question.answers.push_back(StringPrintf("%c", 'a' + correct_idx));
      problem.questions.emplace_back(std::move(question));

      statement.clear();
      thought.clear();
      choices.clear();
  };
  for (int idx = 0; idx < (int)lines.size(); idx++) {
    string &line = lines[idx];
    bool is_thought = false, is_correct = false, is_choice = false;
    if (Util::TryStripPrefix("(*)", &line)) is_correct = is_choice = true;
    else if (Util::TryStripPrefix("*", &line)) is_choice = true;
    else if (Util::TryStripPrefix("[", &line)) is_thought = true;

    if (is_choice) {
      // Asterisk has already been stripped.
      choices.emplace_back(is_correct,
                           CleanWhitespace(line));
    } else if (is_thought) {
      CHECK(thought.empty()) << filename
                             << " multiple thoughts not implemented: "
                             << line;
      CHECK(Util::TryStripSuffix("]", &line)) << filename;
      thought = CleanWhitespace(line);

    } else {
      if (!statement.empty()) EmitQuestion();
      statement = CleanWhitespace(line);
    }
  }

  // Emit final question.
  if (!statement.empty()) EmitQuestion();

  // The first question is the example.
  CHECK(!problem.questions.empty());
  problem.example = std::move(problem.questions[0]);
  problem.questions.erase(problem.questions.begin());

  return problem;
}

#define THOUGHT_BEFORE "Thought: ["
#define THOUGHT_AFTER "]\n"

struct Result {
  Result() {}
  Result(string name, int correct, int total, double sec) :
    name(name), correct(correct), total(total), total_sec(sec) {}
  string name;
  int correct = 0;
  int total = 0;
  double total_sec = 0.0;
};

static string AnsiResultString(const Result &result) {
  const double pct = (result.correct * 100.0) / result.total;
  return StringPrintf(
      APURPLE("%s") ": "
      "Scored " AYELLOW("%d") "/" AWHITE("%d")
      " (" AGREEN("%.2f%%") ") in %s",
      result.name.c_str(),
      result.correct, result.total,
      pct, AnsiTime(result.total_sec).c_str());
}


static Result RunProblem(const Problem &problem,
                         LLM *llm) {
  string prompt = problem.prompt;
  if (USE_THOUGHT) {
    StringAppendF(&prompt,
                  " Before answering, show your thought process, "
                  "also surrounded by square brackets.");
  }
  StringAppendF(&prompt, "\n");
  printf(AWHITE(" == ") APURPLE("%s") AWHITE(" == ") "\n",
         problem.name.c_str());
  printf(AGREY("%s"), prompt.c_str());
  Timer startup_timer;
  llm->Reset();
  llm->DoPrompt(prompt);

  string example_string =
    StringPrintf("%s" "%s" "%s",
                 problem.question_before.c_str(),
                 problem.example.statement.c_str(),
                 problem.question_after.c_str());
  if (USE_THOUGHT) {
    StringAppendF(&example_string,
                  THOUGHT_BEFORE "%s" THOUGHT_AFTER,
                  problem.example.thought.c_str());
  }

  StringAppendF(&example_string,
                "%s" "%s" "%s",
                problem.answer_before.c_str(),
                problem.example.answers[0].c_str(),
                problem.answer_after.c_str());

  printf(AGREY("%s"), example_string.c_str());
  llm->InsertString(example_string);
  const LLM::State start_state = llm->SaveState();
  EmitTimer("started problem", startup_timer);

  int num_correct = 0;
  Timer answer_timer;
  bool first = true;
  for (const Question &question : problem.questions) {
    if (!first) {
      llm->LoadState(start_state);
    }
    first = false;

    string qp =
      StringPrintf("%s" "%s" "%s",
                   problem.question_before.c_str(),
                   question.statement.c_str(),
                   problem.question_after.c_str());

    // Put the cursor after "Answer: ", maybe consuming and printing
    // a thought if that is enabled.
    if (USE_THOUGHT) {
      StringAppendF(&qp, THOUGHT_BEFORE);
      printf(AGREY("%s"), qp.c_str());
      llm->InsertString(qp);
      string thought =
        llm->GenerateUntilEx("]", MAX_THOUGHT, true).first;
      printf(ABLUE("%s") "]", thought.c_str());
      // XXX THOUGHT_AFTER?
      string qa = "\n" + problem.answer_before;
      printf(AGREY("%s"), qa.c_str());
      llm->InsertString(qa);
    } else {
      qp += problem.answer_before;
      printf(AGREY("%s"), qp.c_str());
      llm->InsertString(qp);
    }

    string answer =
      llm->GenerateUntilEx("]", MAX_ANSWER, true).first;
    const bool is_correct =
      [&]() {
        const string lanswer = Util::lcase(answer);
        for (const auto &a : question.answers)
          if (lanswer == Util::lcase(a))
            return true;
        return false;
      }();
    if (answer.empty()) {
      printf(ARED("INVALID") "\n");
    } else if (is_correct) {
      printf(AGREEN("%s") "]\n", answer.c_str());
      num_correct++;
    } else {
      // (Only shows one of the correct answers)
      const string first_answer = question.answers[0];
      printf(ARED("%s") "] (want %s)\n",
             answer.c_str(), first_answer.c_str());
    }
  }

  Result result(problem.name,
                num_correct, (int)problem.questions.size(),
                answer_timer.Seconds());
  printf("== Done ==\n"
         "%s\n", AnsiResultString(result).c_str());
  return result;
}

int main(int argc, char ** argv) {
  AnsiInit();
  Timer model_timer;
  Timer total_timer;

  std::vector<Problem> problems;

  // Load problems first so that parse errors don't cause us to
  // waste a huge amount of time!

  problems.emplace_back(LoadCodeProblem("code.txt"));
  problems.emplace_back(LoadFreeformProblem("blank.txt"));
  problems.emplace_back(LoadFreeformProblem("classical-logic.txt"));
  problems.emplace_back(LoadMultipleChoiceProblem("multiple.txt"));
  problems.emplace_back(LoadFreeformProblem("trivia.txt"));

  problems.emplace_back(LoadFreeformProblem("digits.txt"));
  // Would prefer for this to parse and validate PGN, as well as
  // prompt a different way (no need for Question and Answer).
  problems.emplace_back(LoadFreeformProblem("chess.txt"));

  problems.emplace_back(LoadFreeformProblem("trivia-personal.txt"));
  printf("Loaded " APURPLE("%d") " problems.\n", (int)problems.size());

  ContextParams cparams;
  // cparams.model = "../llama/models/7B/ggml-model-q4_0.bin";
  // cparams.model = "../llama/models/7B/ggml-model-f16.bin";
  // cparams.model = "../llama/models/7B/ggml-model-q8_0.bin";
  // cparams.model = "../llama/models/65B/ggml-model-q4_0.bin";
  // cparams.model = "../llama/models/65B/ggml-model-q8_0.bin";
  // cparams.model = "llama2/7b/ggml-model-f16.gguf";

  cparams.model = "llama2/70b/ggml-model-f16.gguf";

  SamplerParams sparams;
  sparams.type = SampleType::MIROSTAT_2;

  LLM llm(cparams, sparams);
  EmitTimer("Loaded model", model_timer);

  std::map<string, Result> results;

  for (const Problem &problem : problems) {
    Result result = RunProblem(problem, &llm);
    results[result.name] = result;

    // Always print all, since this takes a long time.
    printf("\n\nResults so far:\n");
    for (const auto &[name_, res] : results) {
      printf("%s\n", AnsiResultString(res).c_str());
    }
  }

  llama_print_timings(llm.context.lctx);

  // TODO: Dump full context/sampler params, like to a spreadsheet.
  printf("Model: " AWHITE("%s") "\n", cparams.model.c_str());
  printf("Sample type: " AYELLOW("%s") "\n",
         Sampler::SampleTypeString(sparams.type));
  printf("Using thoughts: " AYELLOW("%s") "\n", USE_THOUGHT ? "YES" : "NO");
  printf("\nTotal benchmark time: %s\n",
         AnsiTime(total_timer.Seconds()).c_str());

  // Save to file.
  std::string out = "\n== Result ==\n";
  StringAppendF(&out, "Finished: %lld\n", (int64_t)time(nullptr));
  StringAppendF(&out, "Took: %lld (%s)\n",
                (int64_t)total_timer.Seconds(),
                ANSI::StripCodes(AnsiTime(total_timer.Seconds())).c_str());
  StringAppendF(&out, "Model: %s\n", cparams.model.c_str());
  StringAppendF(&out, "Sampler: %s\n",
                Sampler::SampleTypeString(sparams.type));
  StringAppendF(&out, "Thoughts: %s\n",
                USE_THOUGHT ? "YES" : "NO");
  for (const auto &[name_, result] : results) {
    const double pct = (result.correct * 100.0) / result.total;
    StringAppendF(&out, "%s: %d/%d in %lld (%.2f%%; %s)\n",
                  result.name.c_str(),
                  result.correct,
                  result.total,
                  (int64_t)result.total_sec,
                  pct,
                  ANSI::StripCodes(AnsiTime(result.total_sec)).c_str());
  }

  static constexpr const char *RESULT_FILE = "benchmark-results.txt";
  {
    FILE *f = fopen(RESULT_FILE, "ab");
    fprintf(f, "%s", out.c_str());
    fclose(f);
  }
  printf("Appended to " AGREEN("%s") ".\n", RESULT_FILE);

  return 0;
}
