/*
 *
 * Copyright (C) 2022 Juan Domingo (Juan.Domingo@uv.es)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "../headers/debugpar_ppam.h"

// This is the only place where this variable is declared. It is global for the full package and should only be changed by ParallelpamSetDebug
unsigned char DEB=NODEBUG;

//' ParellelpamSetDebug
//'
//' Sets debugging in parallelpam package to ON (with TRUE) or OFF (with FALSE) for both parts of it.\cr
//' If this function is not called the default status of debug will be OFF.\cr
//' Setting debugging of any part to ON shows a message. Setting to OFF does not show anything (since debugging is OFF...)
//'
//' @param deb     boolean, TRUE to generate debug messages for the PAM algorithm and FALSE to turn them off.
//' @param debjmat boolean, TRUE to generate debug messages for the jmatrix part inside this package and FALSE to turn them off.
void ParallelpamSetDebug(bool deb,bool debjmat)
{ 
 if (deb)
 {
  DEB |= DEBPP;
  std::cout << "Debugging for PAM algorithm set to ON.\n";
 }
 else
  DEB &= (~DEBPP);
  
 if (debjmat)
 {
  DEB |= DEBJM;
  std::cout << "Debugging for jmatrix inside parallelpam library set to ON.\n";
 }
 else
  DEB &= (~DEBJM);
}

void ParallelpamStop(std::string errortext)
{
 std::cerr << "Error message from the parallelpam library:\n";
 std::cerr << "   " << errortext;
 exit(1);
}

void ParallelpamWarning(std::string warntext)
{
 std::cout << "Warning message from the parallelpam library:\n";
 std::cout << "   " << warntext;
}
