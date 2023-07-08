#include "karate.h"

#include <string>

#include "base/logging.h"
#include "base/stringprintf.h"

#include "re2/re2.h"

using namespace std;

using Aim = Karate::Aim;
using Dir = Karate::Dir;
using Attack = Karate::Attack;
using enum Karate::Aim;
using enum Karate::Dir;
using enum Karate::Attack;

Karate::Karate(const std::string &fighter1, const std::string &fighter2) :
  fighter1(fighter1), fighter2(fighter2) {
}

void Karate::Reset() {
  x1 = 1;
  x2 = 3;
  jump1 = jump2 = nullopt;
  parry1 = parry2 = false;
}

string PosString(int x, std::optional<Dir> jump) {
  if (!jump.has_value()) {
    return StringPrintf("standing on the %dm marker", x);
  } else {
    return StringPrintf("in the air above the %dm marker", x);
  }
}

string Karate::GetStatus() {
  string ret = StringPrintf(
      "Score: %s (%d)  %s (%d)\n"
      "%s is %s, facing right.\n"
      "%s is %s, facing left.\n",
      fighter1.c_str(), score1, fighter2.c_str(), score2,
      fighter1.c_str(), PosString(x1, jump1).c_str(),
      fighter2.c_str(), PosString(x2, jump2).c_str());
  // At most one of these should be true.
  if (parry1)
    StringAppendF(&ret, "%s just had a successful parry.\n", fighter1.c_str());
  if (parry2)
    StringAppendF(&ret, "%s just had a successful parry.\n", fighter2.c_str());

  return ret;
}

#define AIM_RE2 "HIGH|MEDIUM|LOW"
#define DIR_RE2 "LEFT|NEUTRAL|RIGHT"

static std::optional<Aim> GetAim(const std::string &line) {
  if (line == "HIGH") return {HIGH};
  else if (line == "MEDIUM") return {MEDIUM};
  else if (line == "LOW") return {LOW};
  return nullopt;
}

static std::optional<Dir> GetDir(const std::string &line) {
  if (line == "LEFT") return {LEFT};
  else if (line == "NEUTRAL") return {NEUTRAL};
  else if (line == "RIGHT") return {RIGHT};
  return nullopt;
}

std::optional<Aim> Karate::GetBlock(const std::string &line) {
  string aim;
  if (RE2::FullMatch(line, "BLOCK (" AIM_RE2 ")", &aim))
    return GetAim(aim);
  return nullopt;
}

std::optional<Aim> Karate::GetPunch(const std::string &line) {
  string aim;
  if (RE2::FullMatch(line, "PUNCH (" AIM_RE2 ")", &aim))
    return GetAim(aim);
  return nullopt;
}

std::optional<Aim> Karate::GetKick(const std::string &line) {
  string aim;
  if (RE2::FullMatch(line, "KICK (" AIM_RE2 ")", &aim))
    return GetAim(aim);
  return nullopt;
}

std::optional<std::pair<Attack, Aim>> Karate::GetAttack(const std::string &line) {
  if (auto po = GetPunch(line)) return {make_pair(PUNCH, po.value())};
  else if (auto ko = GetKick(line)) return {make_pair(KICK, ko.value())};
  else return nullopt;
}

std::optional<Dir> Karate::GetMove(const std::string &line) {
  string dir;
  if (RE2::FullMatch(line, "MOVE (" DIR_RE2 ")", &dir))
    return GetDir(dir);
  return nullopt;
}

std::optional<Dir> Karate::GetJump(const std::string &line) {
  string dir;
  if (RE2::FullMatch(line, "JUMP (" DIR_RE2 ")", &dir))
    return GetDir(dir);
  return nullopt;
}

static int DirValue(std::optional<Dir> od) {
  if (!od.has_value()) return 0;
  switch (od.value()) {
  default:
  case LEFT: return -1;
  case NEUTRAL: return 0;
  case RIGHT: return +1;
  }
}

string Karate::Process(const std::string &say1, const std::string &say2,
                       const std::string &do1, const std::string &do2) {
  bool prev_parry1 = parry1, prev_parry2 = parry2;
  parry1 = parry2 = false;

  optional<Dir> old_jump1 = jump1, old_jump2 = jump2;
  jump1 = nullopt;
  jump2 = nullopt;

  printf("Process %s|%s|%s|%s\n", say1.c_str(), say2.c_str(), do1.c_str(), do2.c_str());

  // First, move players
  auto mo1 = GetMove(do1);
  auto mo2 = GetMove(do2);

  // If there is a jump in progress, it overrides any explicit move.
  if (old_jump1.has_value()) mo1 = old_jump1;
  if (old_jump2.has_value()) mo2 = old_jump2;

  // A new jump is essentially like a movement, but it also puts
  // the player in the air, and queues up an attack for the next
  // turn.
  auto jo1 = GetJump(do1);
  auto jo2 = GetJump(do2);
  // Jump is ignored if already in the air.
  if (old_jump1.has_value()) jo1 = nullopt;
  if (old_jump2.has_value()) jo2 = nullopt;
  // Otherwise, treat it like a movement.
  if (jo1.has_value()) mo1 = jo1;
  if (jo2.has_value()) mo2 = jo2;

  // Process movement first, possibly ending due to ringout.
  if (mo1.has_value() || mo2.has_value()) {
    // Simultaneous movement.
    // If the players walk into each other,
    // the action fails.
    int nx1 = x1 + DirValue(mo1);
    int nx2 = x2 + DirValue(mo2);
    if (nx1 < nx2) {
      // Success
      x1 = nx1;
      x2 = nx2;
      jump1 = jo1;
      jump2 = jo2;
    }

    // XXX should still be able to jump if the space is occupied.
    // in fact, if just one player is jumping, they should win the race.

    if (x1 < 0 && x2 >= 5) {
      Reset();
      return "Clash: Fighters stepped out of the ring simultaneously.";
    }

    if (x1 < 0) {
      score2++;
      Reset();
      return "Fighter 1 stepped out of the ring.";
    }
    if (x2 >= 5) {
      score1++;
      Reset();
      return "Fighter 2 stepped out of the ring.";
    }

    // Otherwise, continue in new position.
  }

  auto bo1 = GetBlock(do1);
  auto bo2 = GetBlock(do2);

  // Mid-air blocking not allowed.
  if (old_jump1.has_value()) bo1 = nullopt;
  if (old_jump2.has_value()) bo2 = nullopt;

  auto ao1 = GetAttack(do1);
  auto ao2 = GetAttack(do2);

  // No mid-air punching.
  // XXX what to do about kick aims in air?
  if (old_jump1.has_value() &&
      ao1.has_value() && std::get<0>(ao1.value()) != KICK)
    ao1 = nullopt;
  if (old_jump2.has_value() &&
      ao2.has_value() && std::get<0>(ao2.value()) != KICK)
    ao2 = nullopt;

  // If the players are not within striking distance, attacks do nothing.
  if (std::abs(x1 - x2) <= 1) {
    if (ao1.has_value() && ao2.has_value()) {
      // Simultaneous attack.
      auto [a1, aim1] = ao1.value();
      auto [a2, aim2] = ao2.value();

      // If one player is in the air, their kick beats any
      // attack.
      if (old_jump1.has_value() &&
          !old_jump2.has_value()) {
        CHECK(a1 == KICK);
        score1 += 2;
        return "Fighter 1 lands a jump kick.";
      }
      if (old_jump2.has_value() &&
          !old_jump1.has_value()) {
        CHECK(a2 == KICK);
        score2 += 2;
        return "Fighter 2 lands a jump kick.";
      }

      // Otherwise, both players are in the air or both
      // are on the ground, so regular rules apply. Some
      // situations will be impossible in the air, but
      // it doesn't matter.

      // Punch wins.
      if (a1 == PUNCH && a2 == KICK) {
        score1++;
        return "Fighter 1 lands a quick punch.";
      }
      if (a2 == PUNCH && a1 == KICK) {
        score2++;
        return "Fighter 2 lands a quick punch.";
      }

      CHECK(a1 == a2);
      CHECK(!(prev_parry2 && prev_parry2));
      if (prev_parry1) {
        if (a1 == PUNCH) {
          score1++;
          return "Fighter 1 lands a parry-punch combo.";
        } else {
          score1 += 2;
          return "Fighter 1 lands a parry-kick combo.";
        }
      }
      if (prev_parry2) {
        if (a2 == PUNCH) {
          score2++;
          return "Fighter 2 lands a parry-punch combo.";
        } else {
          score2 += 2;
          return "Fighter 2 lands a parry-kick combo.";
        }
      }

      return "Clash: Both fighters attacked simultaneously.";
    }

    // XXX below, we should interpret a low kick for a jumping
    // player as high (and discard other attacks, which inherently
    // miss). We should also discard medium and low attacks for
    // the player on the ground.

    // Now we have some asymmetric attack (if any).
    if (ao1.has_value()) {
      const auto &[a1, aim1] = ao1.value();
      CHECK(!ao2.has_value());

      if (bo2.has_value() && aim1 == bo2.value()) {
        parry2 = true;
        return "Fighter 2 successfully parried.";
      }

      if (a1 == PUNCH) {
        score1++;
        return "Fighter 1 lands a punch.";
      } else {
        score1 += 2;
        return "Fighter 1 lands a kick.";
      }
    }

    if (ao2.has_value()) {
      const auto &[a2, aim2] = ao2.value();
      CHECK(!ao1.has_value());

      if (bo1.has_value() && aim2 == bo1.value()) {
        parry1 = true;
        return "Fighter 1 successfully parried.";
      }

      if (a2 == PUNCH) {
        score2++;
        return "Fighter 2 lands a punch.";
      } else {
        score2 += 2;
        return "Fighter 2 lands a kick.";
      }
    }

  }

  return "";
}
