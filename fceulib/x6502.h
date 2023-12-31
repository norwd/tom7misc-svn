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

#ifndef _FCEULIB_X6502_H
#define _FCEULIB_X6502_H

#include "tracing.h"
#include "fceu.h"
#include "fc.h"

// XXX
#include "base/logging.h"

struct X6502 {
  // Initialize with fc pointer, since memory reads/writes
  // trigger callbacks.
  explicit X6502(FC *fc);

  #ifdef AOT_INSTRUMENTATION
  // Number of times the PC had the corresponding value.
  int64 pc_histo[0x10000] = {};
  // Number of times Run was called with the given number of cycles
  // (prior to multiplication for NTSC/PAL). Last element means that
  // number or greater.
  int64 cycles_histo[1024] = {};
  #endif
  // int64 entered_aot[0x10000] = {};

  /* Temporary cycle counter */
  int32 tcount;

  /* I'll change this to uint32 later...
     I'll need to AND PC after increments to 0xFFFF
     when I do, though.  Perhaps an IPC() macro? */

  // Program counter.
  uint16 reg_PC;
  // A is accumulator, X and Y other general purpose registes.
  // S is the stack pointer. Stack grows "down" (push decrements).
  // P is processor flags (from msb to lsb: NVUBDIZC).
  // PI is another copy of the processor flags, maybe having
  // to do with interrupt handling, which I don't understand.
  uint8 reg_A, reg_X, reg_Y, reg_S, reg_P, reg_PI;
  // I think this is set for some instructions that halt the
  // processor until an interrupt.
  uint8 jammed;

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

private:
  // normal memory read
  inline uint8 RdMem(unsigned int A) {
    return DB = fc->fceu->ARead[A](fc, A);
  }

  // normal memory write
  inline void WrMem(unsigned int A, uint8 V) {
    fc->fceu->BWrite[A](fc, A, V);
  }

  inline uint8 RdRAM(unsigned int A) {
    // PERF: We should read directly from ram in this case (and
    // see what other ones are possible); cheats at this level
    // are not important. -tom7
    //bbit edited: this was changed so cheat substitution would work
    return (DB = fc->fceu->ARead[A](fc, A));
    // return (DB=RAM[A]);
  }

  inline void WrRAM(unsigned int A, uint8 V) {
    fc->fceu->RAM[A] = V;
  }

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
