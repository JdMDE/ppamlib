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

#include "../headers/diftimehelper.h"

DifftimeHelper::DifftimeHelper()
{
 tp.clear();
 messages.clear();
}

void DifftimeHelper::StartClock(std::string message)
{
 std::chrono::steady_clock::time_point newtp=std::chrono::steady_clock::now();
 tp.push_back(newtp);
 messages.push_back(message);
}

double DifftimeHelper::EndClock(bool deb)
{
 std::chrono::steady_clock::time_point just_now=std::chrono::steady_clock::now();
 if (tp.size()==0)
 {
     if (deb)
      std::cout << "Error: unmatched call to EndClock()\n";
     return(0);
     
 }
 std::chrono::steady_clock::time_point before = tp.back();
 tp.pop_back();
 std::chrono::duration<double> elapsed = just_now - before;
 std::string message=messages.back();
 messages.pop_back();
 
 double ret=elapsed.count();
 if (deb)
 {
  std::cout << message << " " << "Elapsed time: " << elapsed.count() << " s\n";
  std::cout.flush();
 }
 return(ret);
}
