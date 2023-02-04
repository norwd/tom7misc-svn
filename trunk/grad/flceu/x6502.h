/* FCE Ultra - NES/Famicom Emulator
 *
 * Copyright notice for this file:
 *  Copyright (C) 2002 Xodnizel
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

#ifndef _X6502_H
#define _X6502_H

#include <cstdint>
#include <bit>

#include "tracing.h"
#include "fceu.h"
#include "fc.h"

#include "fluint8.h"
#include "fluint16.h"

struct X6502 {
  // Initialize with fc pointer, since memory reads/writes
  // trigger callbacks.
  explicit X6502(FC *fc);

  // PERF: Cheap but unnecessary histogram of instructions
  // used.
  int64 inst_histo[256] = {};
  void ClearInstHisto() {
    for (int i = 0; i < 256; i++) inst_histo[i] = 0;
  }

  /* Temporary cycle counter */
  int32 tcount;

  /* I'll change this to uint32 later...
     I'll need to AND PC after increments to 0xFFFF
     when I do, though.  Perhaps an IPC() macro? */

  // XXX make sure if you change something you also change its
  // size in state.cc.

  // Program counter.
  Fluint16 reg_PC;
  // A is accumulator, X and Y other general purpose registes.
  // S is the stack pointer. Stack grows "down" (push decrements).
  // P is processor flags (from msb to lsb: NVUBDIZC).
  // PI is another copy of the processor flags, maybe having
  // to do with interrupt handling, which I don't understand.
  Fluint8 reg_A, reg_X, reg_Y, reg_S, reg_P, reg_PI;
  // I think this is set for some instructions that halt the
  // processor until an interrupt.
  uint8 jammed;

  uint16_t GetPC() const { return reg_PC.ToInt(); }
  uint8 GetA() const { return reg_A.ToInt(); }
  uint8 GetX() const { return reg_X.ToInt(); }
  uint8 GetY() const { return reg_Y.ToInt(); }
  uint8 GetS() const { return reg_S.ToInt(); }
  uint8 GetP() const { return reg_P.ToInt(); }
  uint8 GetPI() const { return reg_PI.ToInt(); }

  int32 count;
  /* Simulated IRQ pin held low (or is it high?).
     And other junk hooked on for speed reasons. */
  uint32 IRQlow;
  /* Data bus "cache" for reads from certain areas */
  uint8 DB;

  void Run(int32 cycles);
  void RunLoop();

  void Init();
  void Reset();
  void Power();

  void TriggerNMI();
  void TriggerNMI2();

  void IRQBegin(int w);
  void IRQEnd(int w);

  // DMA read and write, I think. Like normal memory read/write,
  // but consumes a CPU cycle.
  uint8 DMR(uint32 A);
  void DMW(uint32 A, uint8 V);

  uint32 timestamp = 0;

  void (*MapIRQHook)(FC *, int) = nullptr;

  // TODO: compute some hash of memory reads/writes so that we can
  // ensure that we're preserving the order for them
  uint64_t mem_trace = 0;
  inline void ClearMemTrace() { mem_trace = 0; }
  inline void TraceRead(uint16_t AA) {
    // XXX improve this hash!
    mem_trace ^= AA;
    mem_trace = std::rotr(mem_trace, 31);
    mem_trace = mem_trace * 5 + 0xe6546b64;
  }

  inline void TraceWrite(uint16_t AA, uint8_t VV) {
    // XXX improve this hash!
    mem_trace ^= AA;
    mem_trace = std::rotr(mem_trace, 41);
    mem_trace += VV;
    mem_trace = mem_trace * 5 + 0xe6546b64;
  }

private:
  inline void PUSH(Fluint8 v) {
    WrRAM(0x100 + reg_S.ToInt(), v.ToInt());
    reg_S--;
  }

  inline void PUSH16(Fluint16 vv) {
    WrRAM(0x100 + reg_S.ToInt(), vv.Hi().ToInt());
    reg_S--;
    WrRAM(0x100 + reg_S.ToInt(), vv.Lo().ToInt());
    reg_S--;
  }

  // XXX Fluint8
  inline uint8 POP() {
    return RdRAM(0x100 + (++reg_S).ToInt());
  }

  inline Fluint16 POP16() {
    Fluint8 lo(POP());
    Fluint8 hi(POP());
    return Fluint16(hi, lo);
  }

  /* Indexed Indirect */
  Fluint16 GetIX() {
    Fluint8 tmp = RdMem(reg_PC);
    reg_PC++;
    tmp += reg_X;
    Fluint8 lo = RdRAM(Fluint16(tmp));
    tmp++;
    Fluint8 hi = RdRAM(Fluint16(tmp));
    return Fluint16(hi, lo);
  }

  // Zero Page
  Fluint8 GetZP() {
    Fluint8 ret = RdMem(reg_PC);
    reg_PC++;
    return ret;
  }

  /* Zero Page Indexed */
  Fluint8 GetZPI(Fluint8 i) {
    Fluint8 ret = i + RdMem(reg_PC);
    reg_PC++;
    return ret;
  }

  /* Absolute */
  Fluint16 GetAB() {
    Fluint8 lo = RdMem(reg_PC);
    reg_PC++;
    Fluint8 hi = RdMem(reg_PC);
    reg_PC++;
    return Fluint16(hi, lo);
  }

  /* Absolute Indexed (for writes and rmws) */
  Fluint16 GetABIWR(Fluint8 i) {
    Fluint16 rt = GetAB();
    Fluint16 target = rt;
    target += Fluint16(i);
    (void)RdMem((target & Fluint16(0x00FF)) | (rt & Fluint16(0xFF00)));
    return target;
  }

  /* Absolute Indexed (for reads) */
  Fluint16 GetABIRD(Fluint8 i) {
    Fluint16 tmp = GetAB();
    Fluint16 ret = tmp + Fluint16(i);
    Fluint8::Cheat();
    if (Fluint16::RightShift<8>((ret ^ tmp) & Fluint16(0x100)).ToInt()) {
      (void)RdMem(ret ^ Fluint16(0x100));
      ADDCYC(1);
    }
    return ret;
  }

  /* Indirect Indexed (for reads) */
  Fluint16 GetIYRD() {
    Fluint8 tmp(RdMem(reg_PC));
    reg_PC++;
    Fluint8 lo = RdRAM(Fluint16(tmp));
    Fluint8 hi = RdRAM(Fluint16(tmp + Fluint8(1)));
    Fluint16 rt(hi, lo);
    Fluint16 ret = rt + Fluint16(reg_Y);
    Fluint8::Cheat();
    if (Fluint16::RightShift<8>((ret ^ rt) & Fluint16(0x100)).ToInt()) {
      (void)RdMem(ret ^ Fluint16(0x100));
      ADDCYC(1);
    }
    return ret;
  }


  /* Indirect Indexed(for writes and rmws) */
  Fluint16 GetIYWR() {
    Fluint8 tmp(RdMem(reg_PC));
    reg_PC++;
    Fluint8 lo(RdRAM(Fluint16(tmp)));
    Fluint8 hi(RdRAM(Fluint16(tmp + Fluint8(0x01))));
    Fluint16 rt(hi, lo);
    // PERF directly add Fluint8
    Fluint16 ret = rt + Fluint16(reg_Y);
    (void)RdMem((ret & Fluint16(0x00FF)) | (rt & Fluint16(0xFF00)));
    return ret;
  }

  // Implements the ZNTable (zero flag and negative flag),
  // returning 0, Z_FLAG, or N_FLAG.
  Fluint8 ZnFlags(Fluint8 zort) {
    static_assert(N_FLAG8 == 0x80, "This requires the negative flag "
                  "to be the same as the sign bit.");
    static_assert(Z_FLAG8 == 0x02, "This expects the zero flag to "
                  "have a specific value, although this would be "
                  "easily remedied.");
    Fluint8 zf = Fluint8::LeftShift1Under128(Fluint8::IsZero(zort));
    Fluint8 nf = Fluint8::AndWith<N_FLAG8>(zort);
    // Can't overflow because these are two different bits.
    Fluint8 res = Fluint8::PlusNoOverflow(nf, zf);

    return res;
  }

  void X_ZN(Fluint8 zort) {
    reg_P = Fluint8::AndWith<(uint8_t)~(Z_FLAG8 | N_FLAG8)>(reg_P);
    // We just masked out the bits, so this can't overflow.
    reg_P = Fluint8::PlusNoOverflow(reg_P, ZnFlags(zort));
  }

  void X_ZNT(Fluint8 zort) {
    reg_P |= ZnFlags(zort);
  }

  void LDA(Fluint8 x) {
    reg_A = x;
    X_ZN(reg_A);
  }
  void LDX(Fluint8 x) {
    reg_X = x;
    X_ZN(reg_X);
  }
  void LDY(Fluint8 x) {
    reg_Y = x;
    X_ZN(reg_Y);
  }

  void AND(Fluint8 x) {
    reg_A &= x;
    X_ZN(reg_A);
  }

  void BIT(Fluint8 x) {
    reg_P = Fluint8::AndWith<(uint8_t)~(Z_FLAG8 | V_FLAG8 | N_FLAG8)>(reg_P);
    // PERF: AddNoOverflow
    /* PERF can simplify this ... just use iszero? */
    reg_P |= Fluint8::AndWith<Z_FLAG8>(ZnFlags(x & reg_A));
    reg_P |= Fluint8::AndWith<(uint8_t)(V_FLAG8 | N_FLAG8)>(x);
  }

  void EOR(Fluint8 x) {
    reg_A ^= x;
    X_ZN(reg_A);
  }
  void ORA(Fluint8 x) {
    reg_A |= x;
    X_ZN(reg_A);
  }

  void CMPL(Fluint8 a1, Fluint8 a2) {
    auto [carry, diff] = Fluint8::SubtractWithCarry(a1, a2);
    X_ZN(diff);
    reg_P =
      Fluint8::PlusNoOverflow(
          Fluint8::AndWith<(uint8_t)~C_FLAG8>(reg_P),
          Fluint8::XorWith<C_FLAG8>(
              Fluint8::AndWith<C_FLAG8>(carry)));
  }

  // Input should be 1 or 0.
  void JR(Fluint8 cond) {
    {
      uint8_t cc = cond.ToInt();
      CHECK(cc == 0 || cc == 1) << cc;
    }
    Fluint8::Cheat();
    if (cond.ToInt()) {
      Fluint8::Cheat();
      int32 disp = (int8)RdMem(reg_PC).ToInt();
      reg_PC++;
      ADDCYC(1);
      Fluint16 tmp = reg_PC;
      // Need signed addition
      Fluint8::Cheat();
      if (disp < 0) {
        reg_PC -= Fluint16(-disp);
      } else {
        reg_PC += Fluint16(disp);
      }

      Fluint8::Cheat();
      if (Fluint16::IsntZero((tmp ^ reg_PC) & Fluint16(0x100)).ToInt()) {
        ADDCYC(1);
      }
    } else {
      reg_PC++;
    }
  }

  Fluint8 ASL(Fluint8 x) {
    reg_P = Fluint8::AndWith<(uint8_t)~C_FLAG8>(reg_P);
    reg_P |= Fluint8::RightShift<7>(x);
    x = Fluint8::LeftShift<1>(x);
    X_ZN(x);
    return x;
  }

  Fluint8 LSR(Fluint8 x) {
    reg_P = Fluint8::AndWith<(uint8_t)~(C_FLAG8 | N_FLAG8 | Z_FLAG8)>(reg_P);
    reg_P |= Fluint8::AndWith<1>(x);
    x = Fluint8::RightShift<1>(x);
    X_ZNT(x);
    return x;
  }

  Fluint8 DEC(Fluint8 x) {
    x--;
    X_ZN(x);
    return x;
  }

  Fluint8 INC(Fluint8 x) {
    x++;
    X_ZN(x);
    return x;
  }

  Fluint8 ROL(Fluint8 x) {
    Fluint8 l = Fluint8::RightShift<7>(x);
    x = Fluint8::LeftShift<1>(x);
    x |= Fluint8::AndWith<C_FLAG8>(reg_P);
    reg_P = Fluint8::AndWith<(uint8_t)~(Z_FLAG8 | N_FLAG8 | C_FLAG8)>(reg_P);
    reg_P |= l;
    X_ZNT(x);
    return x;
  }

  Fluint8 ROR(Fluint8 x) {
    Fluint8 l = Fluint8::AndWith<1>(x);
    x = Fluint8::RightShift<1>(x);
    x |= Fluint8::LeftShift<7>(Fluint8::AndWith<C_FLAG8>(reg_P));
    reg_P = Fluint8::AndWith<(uint8_t)~(Z_FLAG8 | N_FLAG8 | C_FLAG8)>(reg_P);
    reg_P |= l;
    X_ZNT(x);
    return x;
  }

  void ADDCYC(int x) {
    this->tcount += x;
    this->count -= x * 48;
    timestamp += x;
  }

  template<class F>
  void ST_ZP(F rf) {
    Fluint16 AA(GetZP());
    WrRAM(AA, rf(AA));
  }

  template<class F>
  void ST_ZPX(F rf) {
    Fluint16 AA(GetZPI(reg_X));
    WrRAM(AA, rf(AA));
  }

  template<class F>
  void ST_ZPY(F rf) {
    Fluint16 AA(GetZPI(reg_Y));
    WrRAM(AA, rf(AA));
  }

  template<class F>
  void ST_AB(F rf) {
    Fluint16 AA = GetAB();
    WrMem(AA, rf(AA));
  }

  template<class F>
  void ST_ABI(Fluint8 reg, F rf) {
    Fluint16 AA = GetABIWR(reg);
    WrMem(AA, rf(AA));
  }

  template<class F>
  void ST_ABX(F rf) {
    return ST_ABI(reg_X, rf);
  }

  template<class F>
  void ST_ABY(F rf) {
    return ST_ABI(reg_Y, rf);
  }

  template<class F>
  void ST_IX(F rf) {
    Fluint16 AA = GetIX();
    WrMem(AA, rf(AA));
  }

  template<class F>
  void ST_IY(F rf) {
    Fluint16 AA = GetIYWR();
    WrMem(AA, rf(AA));
  }

  void ADC(Fluint8 x) {
    static_assert(C_FLAG8 == 0x01, "we assume this is the one's place");
    const Fluint8 p_carry_bit = Fluint8::AndWith<C_FLAG8>(reg_P);
    auto [carry1, sum1] = Fluint8::AddWithCarry(reg_A, x);
    auto [carry2, sum] = Fluint8::AddWithCarry(sum1, p_carry_bit);

    // Since p_carry_bit is at most 1, these can't both overflow.
    Fluint8 carry = Fluint8::PlusNoOverflow(carry1, carry2);

    // uint32 l = reg_A.ToInt() + (x).ToInt() + p_carry_bit.ToInt();
    reg_P = Fluint8::AndWith<
      (uint8_t)~(Z_FLAG8 | C_FLAG8 | N_FLAG8 | V_FLAG8)>(reg_P);
    // The overflow is for signed arithmetic. It tells us if we've
    // added two positive numbers but got a negative one, or added two
    // negative numbers but got a positive one. (If the signs are
    // different, overflow is not possible.) This is computed from the
    // sign bits.
    Fluint8 aaa = Fluint8::XorWith<0x80>(
        Fluint8::AndWith<0x80>(reg_A ^ x));
    Fluint8 bbb = Fluint8::AndWith<0x80>(reg_A ^ sum);
    static_assert(V_FLAG8 == 0x40);

    CHECK((reg_P.ToInt() & (V_FLAG8 | C_FLAG8)) == 0);
    // Sets overflow bit, which was cleared above.
    reg_P = Fluint8::PlusNoOverflow(reg_P, Fluint8::RightShift<1>(aaa & bbb));
    // Sets carry bit, which was cleared above.
    reg_P = Fluint8::PlusNoOverflow(reg_P, carry);
    reg_A = sum;
    // PERF since we already cleared Z and N flags, can use
    // PlusNoOverflow
    X_ZNT(reg_A);
  }

  void SBC(Fluint8 x) {
    static_assert(C_FLAG8 == 0x01, "we assume this is the one's place");
    // On 6502, the borrow flag is !Carry.
    Fluint8 p_ncarry_bit = Fluint8::XorWith<C_FLAG8>(
        Fluint8::AndWith<C_FLAG8>(reg_P));

    auto [carry1, diff1] = Fluint8::SubtractWithCarry(reg_A, x);
    auto [carry2, diff] = Fluint8::SubtractWithCarry(diff1, p_ncarry_bit);

    // As in ADC.
    Fluint8 carry = Fluint8::PlusNoOverflow(carry1, carry2);

    // uint32 l = reg_A.ToInt() - x.ToInt() - p_ncarry_bit.ToInt();
    reg_P = Fluint8::AndWith<
      (uint8_t)~(Z_FLAG8 | C_FLAG8 | N_FLAG8 | V_FLAG8)>(reg_P);
    // As above, detect overflow by looking at sign bits.
    Fluint8 aaa = reg_A ^ diff;
    Fluint8 bbb = reg_A ^ x;
    Fluint8 overflow = Fluint8::AndWith<0x80>(aaa & bbb);
    static_assert(V_FLAG8 == 0x40);

    CHECK((reg_P.ToInt() & (V_FLAG8 | C_FLAG8)) == 0);
    // V_FLAG8 bit is cleared above.
    reg_P = Fluint8::PlusNoOverflow(reg_P, Fluint8::RightShift<1>(overflow));
    // C_FLAG8 bit is cleared above.
    reg_P = Fluint8::PlusNoOverflow(
        reg_P,
        Fluint8::XorWith<C_FLAG8>(Fluint8::AndWith<C_FLAG8>(carry)));
    reg_A = diff;
    // PERF since we already cleared Z and N flags, can use
    // PlusNoOverflow here too
    X_ZNT(reg_A);
  }

  void LSRA() {
    /* For undocumented instructions, maybe for other things later... */
    reg_P = Fluint8::AndWith<(uint8_t)~(C_FLAG8 | N_FLAG8 | Z_FLAG8)>(reg_P);
    reg_P |= Fluint8::AndWith<1>(reg_A);
    reg_A = Fluint8::RightShift<1>(reg_A);
    X_ZNT(reg_A);
  }

  /* Special undocumented operation.  Very similar to CMP. */
  void AXS(Fluint8 x) {
    // TODO: Should be easy with SubtractWithCarry, but we have no
    // test coverage for this instruction :/
    Fluint8::Cheat();
    uint32 t = (reg_A & reg_X).ToInt() - x.ToInt();
    X_ZN(Fluint8(t));
    reg_P =
      Fluint8::PlusNoOverflow(
          Fluint8::AndWith<(uint8_t)~C_FLAG8>(reg_P),
          Fluint8::XorWith<C_FLAG8>(
            Fluint8::AndWith<C_FLAG8>(Fluint8((uint8)(t >> 8)))));
    reg_X = Fluint8((uint8)t);
  }

  // normal memory read, which calls hooks. We trace the sequence
  // of these to make sure we're not dropping or reordering them
  // (which often does not matter to games, but could).
  inline Fluint8 RdMem(Fluint16 A) {
    uint16_t AA = A.ToInt();
    TraceRead(AA);
    return Fluint8(DB = fc->fceu->ARead[AA](fc, AA));
  }

  // normal memory write
  inline void WrMem(Fluint16 A, Fluint8 V) {
    uint16_t AA = A.ToInt();
    uint8_t VV = V.ToInt();
    TraceWrite(AA, VV);
    fc->fceu->BWrite[AA](fc, AA, VV);
  }

  inline uint8 RdRAM(unsigned int A) {
    return (DB = fc->fceu->RAM[A]);
  }

  inline Fluint8 RdRAM(Fluint16 A) {
    return Fluint8(DB = fc->fceu->RAM[A.ToInt()]);
  }

  inline void WrRAM(unsigned int A, uint8 V) {
    fc->fceu->RAM[A] = V;
  }

  inline void WrRAM(Fluint16 A, Fluint8 V) {
    fc->fceu->RAM[A.ToInt()] = V.ToInt();
  }

  // Commonly we do bitwise ops with compile-time constants,
  // which can be faster than fully dynamic fluint operations.
  static constexpr uint8_t N_FLAG8{0x80};
  static constexpr uint8_t V_FLAG8{0x40};
  static constexpr uint8_t U_FLAG8{0x20};
  static constexpr uint8_t B_FLAG8{0x10};
  static constexpr uint8_t D_FLAG8{0x08};
  static constexpr uint8_t I_FLAG8{0x04};
  static constexpr uint8_t Z_FLAG8{0x02};
  static constexpr uint8_t C_FLAG8{0x01};

  // But the constants are also available as fluints.
  static constexpr Fluint8 N_FLAG{N_FLAG8};
  static constexpr Fluint8 V_FLAG{V_FLAG8};
  static constexpr Fluint8 U_FLAG{U_FLAG8};
  static constexpr Fluint8 B_FLAG{B_FLAG8};
  static constexpr Fluint8 D_FLAG{D_FLAG8};
  static constexpr Fluint8 I_FLAG{I_FLAG8};
  static constexpr Fluint8 Z_FLAG{Z_FLAG8};
  static constexpr Fluint8 C_FLAG{C_FLAG8};

  FC *fc = nullptr;
  DISALLOW_COPY_AND_ASSIGN(X6502);
};

#define NTSC_CPU 1789772.7272727272727272
#define PAL_CPU  1662607.125

#define FCEU_IQEXT      0x001
#define FCEU_IQEXT2     0x002
/* ... */
#define FCEU_IQRESET    0x020
#define FCEU_IQNMI2     0x040  // Delayed NMI, gets converted to *_IQNMI
#define FCEU_IQNMI      0x080
#define FCEU_IQDPCM     0x100
#define FCEU_IQFCOUNT   0x200
#define FCEU_IQTEMP     0x800

#endif
