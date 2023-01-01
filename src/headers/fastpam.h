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

#ifndef _FASTPAM_H
#define _FASTPAM_H

#include <utility>    // For pair(,)
#include <map>        // For map(,), in our case, map<pair<unsigned int,unsigned int>,struct stnode>
#include <vector>
#include <typeinfo>

#include <jmatrixlib/fullmatrix.h>
#include <jmatrixlib/symmetricmatrix.h>
#include <jmatrixlib/memhelper.h>

/// @file fastpam.h

///@{
/**
 * The values of this constants is arbitrary, just a mark to distinguish it as a different initalization method\n
 * If you add other initialization method, do it at the end and increase the NUM_INIT_METHODS constant.
 */
const unsigned char INIT_METHOD_PREVIOUS=0;
const unsigned char INIT_METHOD_BUILD=1;
const unsigned char INIT_METHOD_LAB=2;
const unsigned char NUM_INIT_METHODS=3;
///@}

/**
 * Names of the initialization methods. Their positions in the array must coincide with its constant.
 */
const std::string init_method_names[NUM_INIT_METHODS]={"PREV","BUILD","LAB"};

///@{
/**
 * Arbitrary constant, just a mark to distinguish the different algorithms for the optimizacion phase
 */
const unsigned char OPT_METHOD_FASTPAM1=0;
const unsigned char OPT_METHOD_FASTPAMBSIL=1;
const unsigned char NUM_OPT_METHODS=2;
///@}

/**
 * Names of the optimization methods. Their positions in the array must coincide with its constant.
 */
const std::string opt_method_names[NUM_OPT_METHODS]={"FASTPAM1","TWOBRANCH"};

/**
 * The maximum number of iterations we will allow
 */
const unsigned int  MAX_ITER=1001;

/**
 * The maximum number of medoids we allow.
 */
const indextype MAX_MEDOIDS=std::numeric_limits<indextype>::max()-1;

/**
 * A convenience constant to indicate a point is not currently assigned to any medoid
 */
const indextype NO_CLUSTER=MAX_MEDOIDS;

/**
 * @class FastPAM
 * A class to implement the Partitioning Around Medoids (PAM) clustering metho described in\n
 * \n
 * Schubert, E. and Rousseeuw, P.J.: "Fast and eager k-medoids clustering: O(k) runtime improvement of the PAM, CLARA, and CLARANS algorithms."\n
 * Information Systems, vol. 101, p. 101804, 2021.\n
 * doi: https://doi.org/10.1016/j.is.2021.101804\n
 * \n
 * Notice that the actual values of the vectors (instances) are not needed. To recover them, look at the data matrix
 * used to generate the distance matrix.\n
 * The number of instances, N, is never passed since dissimilarity matrix is NxN and therefore its size indicates the N value.\n
 * With respect to the calculated value, it consists on two vectors. The first one has as many components as requested medoids and
 * the second has as many components as instances.\n
 * Medoids are expressed in the first one by its number in the array of points (row in the dissimilarity matrix) starting at 0 (C++ convention).\n
 * The second vector contains the number of the medoid (i.e.: the cluster) to which each instance has been assigned, according to their order in the first vector (also from 0).\n
 * These vectors are returned by the functions GetMedoids and GetAssign (see their respective documentation)
 */
template <typename disttype>
class FastPAM
{
 public:
  /**
   * Default (and only available) constructor
   *
   * @param[in] Dm          A pointer to a SymmetricMatrix which is the distance/dissimilarity matrix
   * @param[in] num_medois  The number of medoids to be found
   * @param[in] initmet     Initialization method (one of the constants INIT_METHOD_PREVIOUS, INIT_METHOD_BUILD or INIT_METHOD_LAB)
   * @param[in] limiter     Maximum number of iterations allowed in the optimization phase. Use 0 to perform only initialization.
   * @param[in] nthreads    Number of threads to be opened. Normally, use the result of function ChooseNumThreads(AS_MANY_AS_POSSIBLE) to get this parameter.
   */
  FastPAM(SymmetricMatrix<disttype> *Dm,indextype num_medoids,unsigned char inimet,int limiter,int nthreads);

  /**
   * This function performs the initialization according to the method set at the class constructor
   *
   * @param[in] initmedoids A vector with the indices of the points that are considered as medoids after the initialization phase.\n
   *            This parameter makes sense (and it is used) ONLY for the initialization method PREV and is probably the result of a previous application of the algorithm, possibly with limiter=0. For other methods it is ignored; just pass an empty vector
   * @param[in] nt          Number of threads to be opened. Normally, use the result of function ChooseNumThreads(AS_MANY_AS_POSSIBLE) to get this parameter.
   */
  void Init(std::vector<indextype> initmedoids, unsigned int nt);

  /**
   * This function runs the optimization phase according to the chosen optimization method
   *
   * @param[in] opt_method Optimization method (one of the constants OPT_METHOD_FASTPAM1 or OPT_METHOD_FASTPAMBSIL)
   * @param[in] nt          Number of threads to be opened. Normally, use the result of function ChooseNumThreads(AS_MANY_AS_POSSIBLE) to get this parameter.
   */
  void Run(unsigned char opt_method,unsigned int nt);
 
  /**
   * This function gets the medoids as a FullMatrix of dimension (num_medoids x 1), i.e. a column vector
   *
   * @return The column vector (as a FullMatrix) with the indices of the medoids in the order they appear in the dissimilarity matrix
   */
  FullMatrix<indextype> &GetMedoids();

   /**
   * This function gets the medoids as a FullMatrix of dimension (num_medoids x 1), i.e. a column vector with point names
   *
   * @param[in] rownames The names of all points as they are stored in the dissimilarity matrix, if it has names.\n
   *                     The function selects specifically those which are medoids and uses them as names for the returned vector.\n
   *                     These parameter can be obtained from the dissimilarity matrix with D.GetRowNames()
   *
   * @return The column vector (as a FullMatrix) with the indices of the medoids in the order they appear in the dissimilarity matrix
   */
  FullMatrix<indextype> &GetMedoids(std::vector<std::string> rownames);

  /**
   * This function gets the medoid to which each point is closest to as a FullMatrix of dimension (num_points x 1), i.e. a column vector
   *
   * @return The column vector (as a FullMatrix) with the indices of the medoids in the vector of medoids (as returned by GetMedoids())\n
   *         Obviously, and since this index is in [0..(num_medoids-1)], is also a class label.
   */
  FullMatrix<indextype> &GetAssign();

  /**
   * This function gets the medoid to which each point is closest to as a FullMatrix of dimension (num_medoids x 1), i.e. a column vector with point names
   *
   * @param[in] rownames The names of all points as they are stored in the dissimilarity matrix, if it has names.\n
   *                     These names are simply attached in the same order to the returned vector.\n
   *                     These parameter can be obtained from the dissimilarity matrix with D.GetRowNames()
   *
   * @return The column vector (as a FullMatrix) with the indices of the medoids in the vector of medoids (as returned by GetMedoids())\n
   *         Obviously, and since this index is in [0..(num_medoids-1)], is also a class label.
   */
  FullMatrix<indextype> &GetAssign(std::vector<std::string> rownames);
  
  /**
   * This function returns the values of the optimization metrics TD (i.e.: the sum of distances of each point to its closest medoid, divided by the number of points) along the succesive optimization iterations.
   *
   * @return The vector with the values of TD for each optimization step
   */
  std::vector<disttype>  GetTDHistory()       { return TDkeep; };

  /**
   * This function returns the number of points that have been swapped between tow clusters along the succesive optimization iterations.
   *
   * @return The vector with the number of swapped points for each optimization step
   */
  std::vector<indextype> GetReassignHistory() { return NpointsChangekeep; };
  
  /**
   * This function returns the total time (in seconds) used for the initialization phase.
   *
   * @return Time (in seconds) used for initalization
   */
  double GetInTime()  { return(time_in_initialization); };

  /**
   * This function returns the total time (in seconds) used for the optimization phase.
   *
   * @return Time (in seconds) used for optimization
   */
  double GetOptTime() { return(time_in_optimization); };

  /**
   * This function returns the number of iterations done in the optimization phase until convergence (or limiter of no convergence is reached)
   *
   * @return Number of iterations used in optimization
   */
  unsigned int GetNumIter() { return(num_iterations_in_opt); };

 private:
  // The multiplicative factor to calculate the threshold for stopping.
  // If the change of TD between consecutive iterations is less than this factor multiplied by the initial TD value we will stop
  // This is to prevent for numerical unstability (very few points move from cluster to cluster getting the algorithm stuck in
  // an infinite loop. This should never happen, but sometimes when working with floats it is possible.
  // It this happens with your data, increase the value of this constant.
  const disttype tlimit=1e-6;

  ///@{
  /**
   * The possible minimum and maximum values the dissimilarity can take. To initialize before loops.
   */
  #define MIND std::numeric_limits<disttype>::min()
  #define MAXD std::numeric_limits<disttype>::max()
  ///@}

  SymmetricMatrix<disttype> *D;  // The dissimilarity matrix
  indextype nmed;             // The number of medoids we want to find
  indextype num_obs;             // The number of observations; it is equal to the number of rows of D, but just for convenience/clarity.
  unsigned char method;          // The initialization method (see constants to codify methods a few lines up)
  unsigned int maxiter;          // Maximum number of iterations we allow
  unsigned int nt;               // Number of threads the user asks for (it may be changed if there are few points)
  bool is_initialized;           // To mark if the chosen initalization algorithm has already be executed.
  
  double time_in_initialization;  // Time in seconds used in the initalization phase (BUILD, ParBUILD or LAB).
  double time_in_optimization;    // Time in seconds used in the optimization phase (FastPAM1 or ParallelFastPAM1).
  unsigned int num_iterations_in_opt;  // Number of itereations used in the optimization phase, never more than maxiter.
  
  // The next fields are filled by initialization (whatever method) and updated by Run
  std::vector<indextype> medoids;     // The current medoids (point index of each one). This is the vector to be returned at the end.
  std::vector<bool>      ismedoid;    // A vector of marks with true for the medoids and false for the others. It is just to accelerate,
                                      // since the medoids vector has already such information.
  std::vector<indextype> nearest;     // The index of the medoid _in_the_array_of_medoids closest to each point
  std::vector<disttype>  dnearest;    // The dissimilarity of every point to its current closest medoid. It plays as a cache
  std::vector<disttype>  dsecond;     // The dissimilarity of every point to its current second closest medoid. It plays as a cache
  
  // These vectors and values are for statistics/information and measures
  disttype               currentTD;          // The value of the optimization function at the current iteration
  std::vector<disttype>  TDkeep;             // Value of TD at each iteration
  indextype              current_npch;       // The value of number of points that have changed cluster at the current iteration
  std::vector<indextype> NpointsChangekeep;  // Number of points that change class at each iteration
 
  // 1) Initialization of the variables; valid for serial and parallel versions
  void InitializeInternals();
  // end 1)
  
  // 2) A function to fill the medoids vector from a previous list; valid for serial and parallel versions
  void InitFromPreviousSet(std::vector<indextype> medoidslist);
  // end 2)
   
  // 3) Internal functions to help in random samples generation
  std::vector<indextype> randomSample(indextype samplesize, indextype n);
  std::vector<indextype> randomSampleExc(indextype samplesize, indextype n,std::vector<bool> &toexclude);

  // 4) Initalization algorithms
  // 4.1) Brute-force initialization algorithm (BUILD), serial version
  void BUILD();

  // 4.2) Brute-force initialization algorithm, parallel version
  void ParBUILD(unsigned int nt);
       
  // 4.2.1) A structure needed for parallel implementation of BUILD
  struct BUILDThread_IO
  {
      FastPAM *FPp;
      indextype *foundmed;
      disttype *DeltaT;
  };
  // end 4.2.1)
  // 4.2.2) Threads needed for parallel implementation of BUILD
  static void *FindFirstMedoidBUILDThread(void *arg);
  static void *FindSuccessiveMedoidBUILDThread(void *arg);
  // end 4.2.2)
            
  // 4.3) Linear approximative build (LAB), serial version
  void LAB();
  // end 4.3)
  // end 4)
  
  // 5) Optimization phase
  // 5.1) Serial version, improved implementation as described in Schubert & Rousseeuw 2021, Algorithm 3
  void RunImprovedFastPAM1();
  // end 5.1)
  
  // 5.1.1) Serial version, improved implementation with my variant, B branches explored simultaneously.
  const unsigned int NBRANCHES=4;

  // This structure contains all necessary data to define a exchange between a medoid and other point
  struct exchange_struct
  {
    disttype DeltaTDst;   // TD improvement that this exchange would provoke
    indextype mst;        // Number of the medoid to be swapped
    indextype xst;        // Number of the point that will be swapped with the medoid
    indextype imst;       // Index in the array of medoids where new point will be put (the place of the current medoid)
  };
  typedef struct exchange_struct exchange;

  void ExploreBranches(disttype *DeltaTDminusm,disttype *DeltaTD,std::vector<exchange> &xcg);
  void ChooseExchange(std::vector<exchange> &xcg,exchange &best_xcg,unsigned int nt);

  void RunImprovedFastPAMMultiBranch(unsigned int branching_index,unsigned int nt);

  // 5.2) Parallel version, improved implementation as described in Schubert & Rousseeuw 2021, Algorithm 3
  void RunParallelImprovedFastPAM1(unsigned int nt);
  // 5.2.1) A structure needed for parallel implementation of optimization
  struct FastPAM1Thread_IO
  {
  	FastPAM *FPp;
  	indextype *mst;
  	indextype *xst;
  	indextype *imst;
  	disttype *DeltaTDst;
  	disttype *DeltaTDminusm;         // This will work as an array
  };
  // end 5.2.1)
  // 5.2.2) Thread needed for parallel implementation of optimization
  static void *FastPAM1InternalThread(void *arg);
  // end 5.2.2)
  // end 5.2)
  // 5.3) Parallel version, implementation of my variation.
  // 5.3.1) A structure needed for parallel implementation of my version of optimization
  struct ExploreBranchesThread_IO
  {
  	FastPAM *FPp;
  	indextype *mst;
  	indextype *xst;
  	indextype *imst;
  	disttype *DeltaTDst;
  	disttype *DeltaTDminusm;         // This will work as an array
  	exchange *xcg;                   // This will work as an array, too.
  };
  // end 5.3.1)
  // 5.3.2) Thread needed for parallel implementation of my variation of optimization
  static void *ExploreBranchesInternalThread(void *arg);
  // end 5.3.2)
  void ExploreBranchesParallel(disttype *DeltaTDminusm,disttype *DeltaTD,std::vector<exchange> &xcg,unsigned int nt);
  // end 5.3)
  // end 5)
  
  // 6) Auxiliary functions used inside all versions of optimization
  void FillSecond();
  void SwapRolesAndUpdate(indextype mst,indextype xst,indextype i);
  // end 6)
};

#endif
