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

#ifndef _FCEULIB_VERSION_H
#define _FCEULIB_VERSION_H

#include <string>

#define FCEU_NAME "FCEUX"

#define FCEU_VERSION_NUMERIC 21060
#define FCEU_VERSION_STRING "2.1.6-fceulib"
#define FCEU_NAME_AND_VERSION FCEU_NAME " " FCEU_VERSION_STRING

std::string FCEUI_GetAboutString();

#endif
