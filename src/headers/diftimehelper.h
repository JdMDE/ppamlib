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

#ifndef _DIFTIME_HELPER_H
#define _DIFTIME_HELPER_H

#include <iostream>
#include <chrono>
#include <vector>
#include <string>

/// @file diftimehelper.h

/**
 * @Difftimehelper class to help in measuring and printing time spent by parts of the programs
 */
class DifftimeHelper
{
 public:
    /**
     * Default constructor
     */
    DifftimeHelper();

    /**
     * Function to start counting of time. It can be called nested inside other call made before; when it finishes each pair of StartClock/EndClock calls will show its message
     *
     * @param[in] message Message that will be printed in the console when the corresponding EndClock that matches this call be called.
     */
    void StartClock(std::string message);

    /**
     * Function to end counting of time elapsed from the last call to StartClock. It prints the message with which StartClock was called if requested.
     *
     * @param[in] deb Boolean value to print or not the message stored by the last call to StartClock
     */
    double EndClock(bool deb);
 private:
    std::vector<std::chrono::steady_clock::time_point> tp;
    std::vector<std::string> messages;
};

#endif
