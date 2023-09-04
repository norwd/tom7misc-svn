//
// This file is part of Alpertron Calculators.
//
// Copyright 2019-2021 Dario Alejandro Alpern
//
// Alpertron Calculators is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Alpertron Calculators is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Alpertron Calculators.  If not, see <http://www.gnu.org/licenses/>.

#ifndef _GLOBALS_H
#define _GLOBALS_H

#include "bignbr.h"

// Stuff used from afar with "extern" moved here.
// Probably we should get rid of as much of this as possible.

// extern char output[3000000];
extern char *output;
extern bool lang;
extern int groupLen;

extern bool hexadecimal;

void beginLine(char** pptrOutput);
void finishLine(char** pptrOutput);

#endif