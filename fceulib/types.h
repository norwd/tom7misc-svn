/* FCE Ultra - NES/Famicom Emulator
 *
 * Copyright notice for this file:
 *  Copyright (C) 2001 Aaron Oneal
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

#ifndef _FCEULIB_TYPES_H
#define _FCEULIB_TYPES_H

#include <cstdint>
#include <sys/stat.h>
#include <unistd.h>

typedef int8_t int8;
typedef int16_t int16;
typedef int32_t int32;
typedef int64_t int64;

typedef uint8_t uint8;
typedef uint16_t uint16;
typedef uint32_t uint32;
typedef uint64_t uint64;

#ifndef DISABLE_SOUND
#define DISABLE_SOUND 0
#endif

#ifndef DISABLE_VIDEO
#define DISABLE_VIDEO 0
#endif

// XXX why are styles 1 and 3 the same?
#if PSS_STYLE==2

#define PSS "\\"
#define PS '\\'

#elif PSS_STYLE==1

#define PSS "/"
#define PS '/'

#elif PSS_STYLE==3

#define PSS "\\"
#define PS '\\'

#elif PSS_STYLE==4

#define PSS ":"
#define PS ':'

#else

#error PSS_STYLE undefined or invalid; see "types.h" for possible values, and add as compile-time option.

#endif

#include "utils/endian.h"

#endif
