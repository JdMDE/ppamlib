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

#ifndef _SILHOUETTE_H
#define _SILHOUETTE_H

#include <thread>

#include <jmatrixlib/fullmatrix.h>
#include <jmatrixlib/sparsematrix.h>
#include <jmatrixlib/symmetricmatrix.h>
#include <jmatrixlib/memhelper.h>

/// @file silhouette.h

/**
 * The values of the silhouette will be stored as double. Since its number is always linear
 * with the number of points, this should not increase too much memory usage and is simpler.
 */
typedef double siltype;

#ifndef DOXYGEN_SHOULD_SKIP_THIS
// A vector of structures like the following one will be created, one structure per point
// Its objective is to be able to return to R a matrix compatible with package cluster which allow representation in standard form.
typedef struct
    {
     indextype pnum;                    // The original point number (to keep it after sorting)
     indextype ownclus;                 // The cluster the point belongs to (C++-numbering, 0..nclus-1)
     indextype neiclus;                 // The cluster which is at minimal average distance (except the own cluster), i.e.: the closest neighbour (C++-numbering)
     siltype  silvalue;                 // The silhouette value
    } silinfo;
    
template <typename disttype>
struct SilhoutteThread_IO
    {
     indextype num_obs;
     indextype nmed;
     const std::vector<indextype> *nearest;
     std::vector<siltype> *current_sil;
     std::vector<unsigned long> *hist;
     std::vector<silinfo> *silres;
     SymmetricMatrix<disttype> *D;
    };
#endif

/**
 * Function to calculate in parallel the silhouette of each point after a clustering has been done\n
 * disttype is the value type used to represent distances in the dissimilarity matrix, either float or double\n
 * siltype is the value type used to store the silhouette, here defined as double
 *
 * @param[in] cl A vector with the class each point belong to, as a number in [0..(num_classes-1)]. Its length must be the number of points, which is the number of rows (and of columns) of the dissimilarity matrix
 * @param[in]  D A reference to the dissimilariry matrix, as a SymmetricMatrix
 * @param[in] nt Number of threads to be opened. Normally, use the result of function ChooseNumThreads(AS_MANY_AS_POSSIBLE) to get this parameter
 *
 * @return A vector with as many components as points containing the silhouette value of each one. Order of points is as in the dissimilarity matrix.
 */
template <typename disttype> std::vector<siltype> CalculateSilhouette(std::vector<indextype> cl,SymmetricMatrix<disttype> &D,unsigned int nt);

/**
 * Function to calculate in parallel the mean values of the silhouette of all points after a clustering has been done\n
 * disttype is the value type used to represent distances in the dissimilarity matrix, either float or double\n
 * siltype is the value type used to store the silhouette, here defined as double
 *
 * @param[in] cl A vector with the class each point belong to, as a number in [0..(num_classes-1)]. Its length must be the number of points, which is the number of rows (and of columns) of the dissimilarity matrix
 * @param[in]  D A reference to the dissimilariry matrix, as a SymmetricMatrix
 * @param[in] nt Number of threads to be opened. Normally, use the result of function ChooseNumThreads(AS_MANY_AS_POSSIBLE) to get this parameter
 *
 * @return The mean value of the silhouette of all points.
 */
template <typename disttype> siltype CalculateMeanSilhouette(std::vector<indextype> cl,indextype nmed,SymmetricMatrix<disttype> *D,unsigned int nt);

#endif
