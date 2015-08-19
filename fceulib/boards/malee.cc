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
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA
 */

#include "mapinc.h"

namespace {
struct MALEE : public CartInterface {
  uint8 WRAM[2048] = {};

  void Power() override {
    fceulib__.cart->setprg2r(0x10, 0x7000, 0);
    fceulib__.fceu->SetReadHandler(0x8000, 0xFFFF, Cart::CartBR);
    fceulib__.fceu->SetReadHandler(0x6000, 0x67FF, Cart::CartBR);
    fceulib__.fceu->SetReadHandler(0x7000, 0x77FF, Cart::CartBR);
    fceulib__.fceu->SetWriteHandler(0x7000, 0x77FF, Cart::CartBW);
    fceulib__.cart->setprg2r(1, 0x6000, 0);
    fceulib__.cart->setprg32(0x8000, 0);
    fceulib__.cart->setchr8(0);
  }

  MALEE(FC *fc, CartInfo *info) : CartInterface(fc) {
    fceulib__.cart->SetupCartPRGMapping(0x10, WRAM, 2048, 1);
    fceulib__.state->AddExState(WRAM, 2048, 0, "WRAM");
  }

};
}

CartInterface *MALEE_Init(FC *fc, CartInfo *info) {
  return new MALEE(fc, info);
}
