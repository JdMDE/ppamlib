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

#ifndef _DISSIMMAT_H
#define _DISSIMMAT_H

#include <iostream>
#include <sstream>
#include <string>
#include <unistd.h>
#include <cmath>

#include <jmatrixlib/fullmatrix.h>
#include <jmatrixlib/sparsematrix.h>
#include <jmatrixlib/symmetricmatrix.h>
#include <jmatrixlib/memhelper.h>

/// @file dissimmat.h

const unsigned char DL1=0x0;   // L1 distance
const unsigned char DL2=0x1;   // L2 distance
const unsigned char DPe=0x2;   // Pearson dissimilarity coefficient

#ifndef DOXYGEN_SHOULD_SKIP_THIS

const unsigned char EMPTY=0x0;
const unsigned char IN_FIRST=0x01;
const unsigned char IN_SECOND=0x02;   // This is IN_FIRST << 1

// From now on, and in the .cpp files, counttype is the value type of the input files and disttype the value type of the dissimilarity matrix (our output)
template <typename counttype,typename disttype>
struct args_to_sp_thread
{
    indextype initial_row1;
    indextype final_row1;
    indextype initial_row2;
    indextype final_row2;
    SparseMatrix<counttype> *M;
    SymmetricMatrix<disttype> *D;
    std::vector<counttype> *mu;
    unsigned char dtype;
};

template <typename counttype,typename disttype>
struct args_to_full_thread
{
    unsigned long initial_row1;
    unsigned long final_row1;
    indextype initial_row2;
    indextype final_row2;
    FullMatrix<counttype> *M;
    SymmetricMatrix<disttype> *D;
    std::vector<counttype> *mu;
    unsigned char dtype;
};
#endif

/**
 * Function to calculate the distance matrix from the data matrix if such matrix is a FullMatrix (in the terminology of the JMatrix library)\n
 * counttype is the data type of the data matrix\n
 * disttype is the data type of the dissimilarity matrix to be returned (use float or double)
 *
 * @param[in] M     The FullMatrix with the data where rows represent individuals (points) and columns are characteristics (dimensions)
 * @param[in] dtype Distance type. Use one of the constants DL1 for Manhattan/City block distance, DL2 for Euclidean distance and Dpe for Pearson dissimilarity coefficient
 * @param[in] nthr  Number of threads to be opened. Normally, use the result of function ChooseNumThreads(AS_MANY_AS_POSSIBLE) to get this parameter.
 */
template <typename counttype,typename disttype>
SymmetricMatrix<disttype> &CalcDistFromFull(FullMatrix<counttype> &M,unsigned char dtype, unsigned int nthr);

/**
 * Function to calculate the distance matrix from the data matrix if such matrix is a SparseMatrix (in the terminology of the JMatrix library)\n
 * counttype is the data type of the data matrix\n
 * disttype is the data type of the dissimilarity matrix to be returned (use float or double)
 *
 * @param[in] M     The SparseMatrix with the data where rows represent individuals (points) and columns are characteristics (dimensions)
 * @param[in] dtype Distance type. Use one of the constants DL1 for Manhattan/City block distance, DL2 for Euclidean distance and Dpe for Pearson dissimilarity coefficient
 * @param[in] nthr  Number of threads to be opened. Normally, use the result of function ChooseNumThreads(AS_MANY_AS_POSSIBLE) to get this parameter.
 */
template <typename counttype,typename disttype>
SymmetricMatrix<disttype> &CalcDistFromSparse(SparseMatrix<counttype> &M,unsigned char dtype,unsigned int nthr);

#endif
