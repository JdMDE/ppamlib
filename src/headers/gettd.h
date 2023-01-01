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

#ifndef _GETTD_H
#define _GETTD_H

#include <jmatrixlib/symmetricmatrix.h>

/// @file gettd.h

/**
 * Function to get the value of the metrics usually employed in PAM minimization: sum of distances of each point to its closest medoid divided by number of points.
 *
 * @param[in] Lmed    A vector with the indices of the points which are medoids. These indices refer to the order of points in the distance/dissimilarity matrix
 * @param[in] Lclasif A vector with the index (as position in Lmed) of the medoid closest to each point
 * @param[in] D       A reference to the dissimilariry matrix, as a SymmetricMatrix
 *
 * @return The value of the total sum of distances divided by the number of points
 */
template<typename disttype> double GetTD(std::vector<indextype> Lmed,std::vector<indextype> Lclasif,SymmetricMatrix<disttype> &D);

#endif
