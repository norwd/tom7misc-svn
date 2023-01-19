#include <string.h>
#include "suborkb.h"

#define FKB_ESCAPE 0x01
#define FKB_F1 0x02
#define FKB_F2 0x03
#define FKB_F3 0x04
#define FKB_F4 0x05
#define FKB_F5 0x06
#define FKB_F6 0x07
#define FKB_F7 0x08
#define FKB_F8 0x09
#define FKB_F9 0x0A
#define FKB_F10 0x0B
#define FKB_F11 0x0C
#define FKB_F12 0x0D
#define FKB_PAUSE 0x0E
#define FKB_GRAVE 0x0F
#define FKB_1 0x10
#define FKB_2 0x11
#define FKB_3 0x12
#define FKB_4 0x13
#define FKB_5 0x14
#define FKB_6 0x15
#define FKB_7 0x16
#define FKB_8 0x17
#define FKB_9 0x18
#define FKB_0 0x19
#define FKB_MINUS 0x1A
#define FKB_EQUALS 0x1B
#define FKB_BACK 0x1C
#define FKB_INSERT 0x1D
#define FKB_HOME 0x1E
#define FKB_PRIOR 0x1F
#define FKB_NUMLOCK 0x20
#define FKB_DIVIDE 0x21
#define FKB_MULTIPLY 0x22
#define FKB_SUBTRACT 0x23
#define FKB_TAB 0x24
#define FKB_Q 0x25
#define FKB_W 0x26
#define FKB_E 0x27
#define FKB_R 0x28
#define FKB_T 0x29
#define FKB_Y 0x2A
#define FKB_U 0x2B
#define FKB_I 0x2C
#define FKB_O 0x2D
#define FKB_P 0x2E
#define FKB_LBRACKET 0x2F
#define FKB_RBRACKET 0x30
#define FKB_RETURN 0x31
#define FKB_DELETE 0x32
#define FKB_END 0x33
#define FKB_NEXT 0x34
#define FKB_NUMPAD7 0x35
#define FKB_NUMPAD8 0x36
#define FKB_NUMPAD9 0x37
#define FKB_ADD 0x38
#define FKB_CAPITAL 0x39
#define FKB_A 0x3A
#define FKB_S 0x3B
#define FKB_D 0x3C
#define FKB_F 0x3D
#define FKB_G 0x3E
#define FKB_H 0x3F
#define FKB_J 0x40
#define FKB_K 0x41
#define FKB_L 0x42
#define FKB_SEMICOLON 0x43
#define FKB_APOSTROPHE 0x44
#define FKB_NUMPAD4 0x45
#define FKB_NUMPAD5 0x46
#define FKB_NUMPAD6 0x47
#define FKB_LSHIFT 0x48
#define FKB_Z 0x49
#define FKB_X 0x4A
#define FKB_C 0x4B
#define FKB_V 0x4C
#define FKB_B 0x4D
#define FKB_N 0x4E
#define FKB_M 0x4F
#define FKB_COMMA 0x50
#define FKB_PERIOD 0x51
#define FKB_SLASH 0x52
#define FKB_BACKSLASH 0x53
#define FKB_UP 0x54
#define FKB_NUMPAD1 0x55
#define FKB_NUMPAD2 0x56
#define FKB_NUMPAD3 0x57
#define FKB_LCONTROL 0x58
#define FKB_LMENU 0x59
#define FKB_SPACE 0x5A
#define FKB_LEFT 0x5B
#define FKB_DOWN 0x5C
#define FKB_RIGHT 0x5D
#define FKB_NUMPAD0 0x5E
#define FKB_DECIMAL 0x5F

#define AK2(x, y) ((FKB_##x) | (FKB_##y << 8))
#define AK(x) FKB_##x

static constexpr uint16 const matrix[13][2][4] = {
    {{AK(4), AK(G), AK(F), AK(C)}, {AK(F2), AK(E), AK(5), AK(V)}},
    {{AK(2), AK(D), AK(S), AK(END)}, {AK(F1), AK(W), AK(3), AK(X)}},
    {{AK(INSERT), AK(BACK), AK(NEXT), AK(RIGHT)},
     {AK(F8), AK(PRIOR), AK(DELETE), AK(HOME)}},
    {{AK(9), AK(I), AK(L), AK(COMMA)}, {AK(F5), AK(O), AK(0), AK(PERIOD)}},
    {{AK(RBRACKET), AK(RETURN), AK(UP), AK(LEFT)},
     {AK(F7), AK(LBRACKET), AK(BACKSLASH), AK(DOWN)}},
    {{AK(Q), AK(CAPITAL), AK(Z), AK(TAB)},
     {AK(ESCAPE), AK(A), AK(1), AK(LCONTROL)}},
    {{AK(7), AK(Y), AK(K), AK(M)}, {AK(F4), AK(U), AK(8), AK(J)}},
    {{AK(MINUS), AK(SEMICOLON), AK(APOSTROPHE), AK(SLASH)},
     {AK(F6), AK(P), AK(EQUALS), AK(LSHIFT)}},
    {{AK(T), AK(H), AK(N), AK(SPACE)}, {AK(F3), AK(R), AK(6), AK(B)}},
    {{0, 0, 0, 0}, {0, 0, 0, 0}},
    {{AK(LMENU), AK(NUMPAD4), AK(NUMPAD7), AK(F11)},
     {AK(F12), AK(NUMPAD1), AK(NUMPAD2), AK(NUMPAD8)}},
    {{AK(SUBTRACT), AK(ADD), AK(MULTIPLY), AK(NUMPAD9)},
     {AK(F10), AK(NUMPAD5), AK(DIVIDE), AK(NUMLOCK)}},
    {{AK(GRAVE), AK(NUMPAD6), AK(PAUSE), AK(SPACE)},
     {AK(F9), AK(NUMPAD3), AK(DECIMAL), AK(NUMPAD0)}},
};

namespace {
struct SuborKB final : public InputCFC {
  using InputCFC::InputCFC;

  void Write(uint8 v) override {
    v >>= 1;
    if (v & 2) {
      if ((ksmode & 1) && !(v & 1)) ksindex = (ksindex + 1) % 13;
    }
    ksmode = v;
  }

  uint8 Read(int w, uint8 ret) override {
    if (w) {
      ret &= ~0x1E;
      //  if(ksindex==9)
      //  {
      //     if(ksmode&1)
      //        ret|=2;
      //  }
      //  else
      //  {
      for (int x = 0; x < 4; x++)
        if (bufit[matrix[ksindex][ksmode & 1][x] & 0xFF] ||
            bufit[matrix[ksindex][ksmode & 1][x] >> 8])
          ret |= 1 << (x + 1);
      //  }
      ret ^= 0x1E;
    }
    return ret;
  }

  void Strobe() override {
    ksmode = 0;
    ksindex = 0;
  }

  void Update(void *data, int arg) override {
    memcpy(bufit + 1, data, 0x60);
  }

  uint8 bufit[0x61] = {};
  uint8 ksmode = 0;
  uint8 ksindex = 0;
};
}  // namespace


InputCFC *CreateSuborKB(FC *fc) {
  return new SuborKB(fc);
}
