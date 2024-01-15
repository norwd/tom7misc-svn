
#include "lex.h"

#include <cstdint>
#include <string>
#include <vector>

#include "ansi.h"
#include "re2/re2.h"
#include "base/logging.h"
#include "base/stringprintf.h"

const char *TokenTypeString(TokenType tok) {
  switch (tok) {
  case LPAREN: return "LPAREN";
  case RPAREN: return "RPAREN";
  case LBRACKET: return "LBRACKET";
  case RBRACKET: return "RBRACKET";
  case COMMA: return "COMMA";
  case PERIOD: return "PERIOD";

  case DIGITS: return "DIGITS";

  // Keywords.
  case FN: return "FN";

  // Identifier.
  case ID: return "ID";

  case LAYOUT_LIT: return "LAYOUT_LIT";
  case STR_LIT: return "STR_LIT";

  case INVALID: return "INVALID";

  default: return "???UNHANDLED???";
  }
}


// Background colors to highlight tokens.
static std::array<uint32_t, 6> COLORS = {
  0x000066FF,
  0x006600FF,
  0x660000FF,
  0x660066FF,
  0x006666FF,
  0x666600FF,
};

std::pair<std::string, std::string> ColorTokens(
    const std::string &input_string,
    const std::vector<Token> &tokens) {
  std::string source, ctokens;
  int color_idx = 0;
  size_t input_pos = 0;
  for (int t = 0; t < (int)tokens.size(); t++) {
    const Token &tok = tokens[t];
    /*
    printf("pos %zu. tok %d from %zu for %zu\n", input_pos,
           t, tok.start, tok.length);
    */
    CHECK(input_pos <= tok.start);
    CHECK(tok.start < input_string.size());
    CHECK(tok.start + tok.length <= input_string.size());
    // If there is whitespace to skip, emit it with black bg.
    if (input_pos < tok.start) {
      source += ANSI::BackgroundRGB(0x22, 0x22, 0x22);
      // StringAppendF(&source, "%d", tok.start - input_pos);
      while (input_pos < tok.start) {
        source.push_back(input_string[input_pos]);
        input_pos++;
      }
    }

    CHECK(input_pos == tok.start);
    const uint32_t c = COLORS[color_idx];
    color_idx++;
    color_idx %= COLORS.size();

    std::string bc = ANSI::BackgroundRGB32(c);
    std::string fc = ANSI::ForegroundRGB32(c | 0x808080FF);
    std::string in = input_string.substr(tok.start, tok.length);
    StringAppendF(&source, "%s%s", bc.c_str(), in.c_str());
    StringAppendF(&ctokens, "%s%s" ANSI_RESET " ",
                  fc.c_str(), TokenTypeString(tok.type));
    input_pos += tok.length;
  }
  source += ANSI_RESET;
  ctokens += ANSI_RESET;
  return std::make_pair(source, ctokens);
}

// Lexing.
std::vector<Token> Lex(const std::string &input_string) {
  static const RE2 whitespace("[ \r\n\t]+");
  static const RE2 digits("[0-9]+");
  // TODO: Allow UTF-8
  static const RE2 ident("([A-Za-z][-'_A-Za-z0-9]*)");
  static const RE2 strlit(
    // Starting with double quote
    "\""
    // Then some number of characters.
    "(?:"
    // Any character other than " or newline or backslash.
    R"([^"\\\r\n])"
    "|"
    // Or an escaped character. This is permissive; we parse
    // the string literals and interpret escape characters
    // separately.
    R"(\\.)"
    ")*"
    "\"");

  // Regex for the interior of a [layout literal].
  // When we use antiquote for a nested expression, like
  // [aaa[exp]bbb], aaa and bbb are lexed as separate tokens.
  // Unlike string literals, which include their quotation
  // marks in the token, when we lex a layout literal we
  // also output tokens for LBRACKET and RBRACKET as appropriate.
  //
  // We don't actually have an escape character inside layout.
  // The way to write a ] or [ is to use expression mode, like ["["].
#define ANY_BRACKET R"([\[\]])"
  static const RE2 layoutlit(
      ANY_BRACKET
      "(?:"
      R"([^\[\]])" "*"
      ")*"
      ANY_BRACKET);

  static const std::unordered_map<std::string, TokenType> keywords = {
    {"fn", FN},
  };

  re2::StringPiece input(input_string);

  // Get the current offset of the stringpiece.
  auto Pos = [&input_string, &input]() {
      CHECK(input.data() >= input_string.data());
      CHECK(input.data() <= input_string.data() + input_string.size());
      return input.data() - input_string.data();
    };

  std::vector<Token> ret;

  while (!input.empty()) {
    const size_t start = Pos();
    std::string match;
    if (RE2::Consume(&input, whitespace)) {
      /* no tokens. */
      // printf("Saw whitespace at %zu for %zu\n", start, Pos() - start);
    } else if (RE2::Consume(&input, digits)) {
      ret.emplace_back(DIGITS, start, Pos() - start);
    } else if (RE2::Consume(&input, ident, &match)) {
      const auto kit = keywords.find(match);
      if (kit != keywords.end()) {
        ret.emplace_back(kit->second, start, Pos() - start);
      } else {
        ret.emplace_back(ID, start, Pos() - start);
      }
    } else if (RE2::Consume(&input, strlit)) {
      ret.emplace_back(STR_LIT, start, Pos() - start);
    } else if (RE2::Consume(&input, layoutlit)) {
      const size_t len = Pos() - start;
      // Must match open and closing bracket, at least.
      CHECK(len >= 2);
      TokenType brack1 = INVALID, brack2 = INVALID;
      char c1 = input_string[start];
      char c2 = input_string[start + len - 1];
      switch (c1) {
      case '[': brack1 = LBRACKET; break;
      case ']': brack1 = RBRACKET; break;
      default: CHECK(false) << c1;
      }
      switch (c2) {
      case '[': brack2 = LBRACKET; break;
      case ']': brack2 = RBRACKET; break;
      default: CHECK(false) << c2;
      }

      ret.emplace_back(brack1, start, 1);
      ret.emplace_back(LAYOUT_LIT, start + 1, len - 2);
      ret.emplace_back(brack2, start + len - 1, 1);
    } else {
      char c = input[0];
      switch (c) {
        // Could use a table for this...
      case '(':
        ret.emplace_back(LPAREN, start, 1);
        input.remove_prefix(1);
        continue;
      case ')':
        ret.emplace_back(RPAREN, start, 1);
        input.remove_prefix(1);
        continue;
      case ',':
        ret.emplace_back(COMMA, start, 1);
        input.remove_prefix(1);
        continue;
      case '.':
        ret.emplace_back(PERIOD, start, 1);
        input.remove_prefix(1);
        continue;

      default: {
        const char *msg = "Unexpected character";
        switch (c) {
        case '\"': msg = "Invalid string literal starting with"; break;
        case '[':
        case ']': msg = "Invalid layout literal starting with"; break;
        default: break;
        }
        // XXX get line info
        CHECK(false) <<
          StringPrintf("%s '%c' (0x%02x) at offset %zu\n",
                       msg, c, c, start) << "\nIn input: "
                     << input_string;
      }
      }
    }
  }
  return ret;
}
