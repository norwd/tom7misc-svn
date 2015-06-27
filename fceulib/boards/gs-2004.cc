/* FCE Ultra - NES/Famicom Emulator
 *
 * Copyright notice for this file:
 *  Copyright (C) 2007 CaH4e3
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

#include "mapinc.h"

static uint8 reg, mirr;
static SFORMAT StateRegs[]=
{
  {&reg, 1, "REGS"},
  {&mirr, 1, "MIRR"},
  {0}
};

static void Sync(void)
{
  fceulib__cart.setprg8r(1,0x6000,0);
  fceulib__cart.setprg32(0x8000,reg);
  fceulib__cart.setchr8(0);
}

static DECLFW(BMCGS2004Write)
{
  reg=V;
  Sync();
}

static void BMCGS2004Power(void)
{
  reg=~0;
  Sync();
  fceulib__fceu.SetReadHandler(0x6000,0x7FFF,Cart::CartBR);
  fceulib__fceu.SetReadHandler(0x8000,0xFFFF,Cart::CartBR);
  fceulib__fceu.SetWriteHandler(0x8000,0xFFFF,BMCGS2004Write);
}

static void BMCGS2004Reset(void)
{
  reg=~0;
}

static void StateRestore(int version)
{
  Sync();
}

void BMCGS2004_Init(CartInfo *info)
{
  info->Reset=BMCGS2004Reset;
  info->Power=BMCGS2004Power;
  fceulib__fceu.GameStateRestore=StateRestore;
  AddExState(&StateRegs, ~0, 0, 0);
}
