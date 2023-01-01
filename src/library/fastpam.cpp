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

#include <unordered_map>
#include <random>
#include "../headers/fastpam.h"
#include "../headers/threadhelper.h"
#include "../headers/diftimehelper.h"
#include "../headers/debugpar_ppam.h"
#include "../headers/silhouette.h"

extern unsigned char DEB;

using namespace std;

// Returns a vector with a random sample of samplesize numbers uniformly chosen from range 0..n-1
template <typename disttype>
vector<indextype> FastPAM<disttype>::randomSample(indextype samplesize, indextype n)
{
        vector<indextype> samples(samplesize);

        unordered_map<indextype, bool> sample_dict;
        //unordered_map<int, bool>::iterator it;

	    // Only for R
        // GetRNGstate();

        std::random_device r;
        std::seed_seq seed{r(), r(), r(), r(), r(), r(), r(), r()};
        std::mt19937 eng(seed);
        std::uniform_int_distribution<> dist{1,int(n)-1};

        indextype i=0;
        while (i < samplesize)
        {
            // Only for R
            // indextype rnd = indextype(double(n)*unif_rand());
            indextype rnd = dist(eng);
            if (sample_dict.find(rnd) == sample_dict.end())
            {
                samples[i] = rnd;
                sample_dict[rnd] = true;
                i++;
            }
        }
        // Only for R
        // PutRNGstate();

        return samples;
}

// Returns a vector random sample of samplesize numbers uniformly chosen from range 0..n-1 but which have no 'true' mark in the n-sized array of booleans
template<typename disttype>
vector<indextype> FastPAM<disttype>::randomSampleExc(indextype samplesize, indextype n,vector<bool> &toexclude)
{
        vector<indextype> samples(samplesize);

        unordered_map<indextype, bool> sample_dict;
        //std::unordered_map<int, bool>::iterator it;

        for (indextype e=0;e<n;e++)
         if (toexclude[e])
          sample_dict[e]=true;

        // Only for R
	    // GetRNGstate();

        indextype i=0;

        std::random_device r;
        std::seed_seq seed{r(), r(), r(), r(), r(), r(), r(), r()};
        std::mt19937 eng(seed);
        std::uniform_int_distribution<> dist{1,int(n)-1};

        while (i < samplesize)
        {
            // Only for R
            // indextype rnd = indextype(double(n)*unif_rand());
            indextype rnd = dist(eng);
            if (sample_dict.find(rnd) == sample_dict.end())
            {
                samples[i] = rnd;
                sample_dict[rnd] = true;
                i++;
            }
        }
        // Only for R
        // PutRNGstate();

        return samples;
}

template <typename disttype>
FastPAM<disttype>::FastPAM(SymmetricMatrix<disttype> *Dm, indextype num_medoids, unsigned char inimet, int limiter, int nthreads)
{
 D = Dm;
 nmed = num_medoids;
 is_initialized = false;

 if (!D->TestDistDisMat())
 {
  ostringstream errst;
  errst << "  Sorry, the matrix is not a distance/dissimilarity matrix.\n";
  errst << "  It has either non-zero elements in the main diagonal or null or negative elements outside it.\n";
  errst << "  The PAM algorithm does not work with this type of matrices.\n";
  ParallelpamStop(errst.str());
 }
 if (DEB & DEBPP)
  std::cout << "  Matrix is a correct distance/dissimilarity matrix.\n";

 indextype npoints=Dm->GetNRows();
 if (npoints<1000)
 {
  nt=1;
  if (DEB & DEBPP)
   std::cout << "Calculating with a single thread, since you have only " << npoints << " vectors and the overhead of using threads would be excessive.\n";
 }
 else
  nt=ChooseNumThreads(nthreads);

 time_in_initialization = 0.0;
 time_in_optimization = 0.0;
 num_iterations_in_opt = 0;
  
 // Even not strictly needed, this variable makes the code clearer
 num_obs=D->GetNRows();
 
 if (inimet >= NUM_INIT_METHODS)
 {
     ParallelpamStop("Error: unknown method passed to FastPAM constructor.\n");
     return;
 }
 method=inimet;
 
 if ((unsigned int)limiter>MAX_ITER)
 {
     ostringstream errst;
     errst << "Error: maximum number of iteration limited to " << MAX_ITER << ".\n";
     errst << "If you need more, change the constant MAX_ITER at fastpam.h and recompile.\n";
     ParallelpamStop(errst.str());
     return;
 }
 if (limiter==0)
  maxiter=0;
 else
  maxiter=(unsigned int)limiter-1;
  
 // No medoids found yet...
 medoids.clear();
 
 // All this vectors have as length the number of points
 ismedoid.resize(num_obs);
 nearest.resize(num_obs);
 dnearest.resize(num_obs);
 
 for (indextype q=0; q<num_obs; q++)
 {
  // The vector of marks is initialized to all-false.
  ismedoid[q]=false;
  // Initial assignment: all points are not yet assigned to any cluster
  nearest[q]=NO_CLUSTER;
  // Therefore, its distance to closest medoid is infinite
  dnearest[q]=MAXD;
 }
 
 // The vectors of TD data as long as the current values are cleared
 TDkeep.clear();
 currentTD=MAXD;
 NpointsChangekeep.clear();
 current_npch=0;
 time_in_initialization=time_in_optimization=0.0;
}

template FastPAM<float>::FastPAM(SymmetricMatrix<float> *Dm,indextype num_medoids,unsigned char imet,int miter,int nthreads);
template FastPAM<double>::FastPAM(SymmetricMatrix<double> *Dm,indextype num_medoids,unsigned char imet,int miter,int nthreads);

/********* InitializeInternals **************/
template <typename disttype>
void FastPAM<disttype>::InitializeInternals()
{
 // This function is called when the medoids vector has been populated with the chosen number of medoids, i.e. after initialization with BUILD, LAB or PREV.
 
 // Apart from the medoids' vector, there are other things we must fill: ismedoid, nearest, dnearest and TD
 
 // No point is a medoid, except those explictly stated in the array medoids
 for (indextype q=0; q<num_obs; q++)
  ismedoid[q]=false;
 for (indextype m=0; m<nmed; m++)
  ismedoid[medoids[m]]=true;
 
 disttype d,mindist;
 indextype index_of_mindist;
 currentTD = disttype(0);
 for (indextype q=0; q<num_obs; q++)
 {
     // A simple sequential search of the minimum distance (and to which medoid) along the array of medoids
     index_of_mindist=nmed+1;
     mindist=MAXD;
     for (indextype m=0; m<nmed; m++)
      if ( (d=D->Get(q,medoids[m])) < mindist )
      {
       mindist=d;
       index_of_mindist=m;
      }
    
     if (index_of_mindist>nmed)
     {
      ostringstream errst;
      errst << "Point " << q << " does not seem to have a closest medoid. Unexpected error.\n";
      ParallelpamStop(errst.str());
      return;
     } 
     
     // The nearest array contain the index in the array of medoids of the medoid closest to point q
     nearest[q]=index_of_mindist;
     // dnearest contains the distance to such medoid...
     dnearest[q]=mindist;
     // and the TD function is updated adding up such distance.
     currentTD += mindist;
 }  
}

template void FastPAM<float>::InitializeInternals();
template void FastPAM<double>::InitializeInternals();

/**************** Init **********************/
template <typename disttype>
void FastPAM<disttype>::Init(std::vector< indextype > initmedoids, unsigned int nt)
{
 switch (method)
 {
     case INIT_METHOD_PREVIOUS: InitFromPreviousSet(initmedoids); break; 
     case INIT_METHOD_BUILD:
     {
         DifftimeHelper Dt;
         if ( nt==1 || D->GetNRows()<1000)
         {
             Dt.StartClock("BUILD initialization method (serial version) finished.");
             BUILD();
             time_in_initialization=Dt.EndClock(DEB & DEBPP);
         }
         else
         {
             Dt.StartClock("BUILD initialization method (parallel version) finished.");
             ParBUILD(nt);
             time_in_initialization=Dt.EndClock(DEB & DEBPP);
         }
         break;
     }
     case INIT_METHOD_LAB: 
     {
         // LAB is not yet implemented in parallel, so we don't check the number of threads.
     	 DifftimeHelper Dt;
     	 Dt.StartClock("LAB initialization method (serial version) finished.");
     	 LAB();
     	 time_in_initialization=Dt.EndClock(DEB & DEBPP);
         break;
     }
     default: ParallelpamStop("Unknown initialization method.\n"); break;
 }
 
 // Mark to indicate that initalization phase, whatever variant, is finished. Checked by swap algorithm before start
 is_initialized=true;
 
 // Called to initialize some internal variables. See actual function code up.
 InitializeInternals();
}

template void FastPAM<float>::Init(std::vector<indextype> initmedoids,unsigned int nt);
template void FastPAM<double>::Init(std::vector<indextype> initmedoids,unsigned int nt);

/************ Run ******************/
template <typename disttype>
void FastPAM<disttype>::Run(unsigned char opt_method,unsigned int nt)
{
    if (!is_initialized)
    {
        ParallelpamStop("Function FastPAM::Run(int nthreads) called before calling FastPAM::Init()\n");
        return;
    }
    if (maxiter==0)
     return;

    DifftimeHelper Dt;
    if (nt==1)
    {
     switch (opt_method)
     {
         case OPT_METHOD_FASTPAM1:
             Dt.StartClock("Optimization method FASTPAM1 (serial version) finished.");
             RunImprovedFastPAM1();
             break;
         case OPT_METHOD_FASTPAMBSIL:
             Dt.StartClock("Optimization method TWOBRANCH (serial version) finished.");
             RunImprovedFastPAMMultiBranch(NBRANCHES,nt);
             break;
         default: ParallelpamStop("Unexpected error in Run: unknonw optimization method.\n"); break;
     }
     time_in_optimization=Dt.EndClock(DEB & DEBPP);
    }
    else
    {
     switch (opt_method)
     {
         case OPT_METHOD_FASTPAM1:
             Dt.StartClock("Optimization method FASTPAM1 (parallel version) finished.");
             RunParallelImprovedFastPAM1(nt);
             break;
         case OPT_METHOD_FASTPAMBSIL:
             // Yes, the same fuction as in the serial version is called, but the value of nt here will be bigger than 1.
             RunImprovedFastPAMMultiBranch(NBRANCHES,nt);
             break;
         default: ParallelpamStop("Unexpected error in Run: unknonw optimization method.\n"); break;
     }
     time_in_optimization=Dt.EndClock(DEB & DEBPP);
    }

    if (DEB & DEBPP)
    {
     std::cout << "Time summary ";
     if (nt==1)
      std::cout << " (serial implementation).\n";
     else
      std::cout << " (parallel implementation with " << nt << " threads).\n";
     std::cout << "   Initalization: " << GetInTime() << " s (method " << init_method_names[method] << ").\n";
     std::cout << "   Optimization:  " << GetOptTime() << " s in " << GetNumIter() << " iterations";
     if (GetNumIter()!=0)
      std::cout << " (" << GetOptTime()/double(GetNumIter()) << " seconds/iteration).\n";
     else
      std::cout << ".\n";

     double tt=GetInTime()+GetOptTime();
     std::cout << "   Total time:    " << tt << " s (" << int(tt/60.0) << " minutes, " << tt-60.0*int(tt/60.0) << " seconds).\n";
   }
}

template void FastPAM<float>::Run(unsigned char opt_method,unsigned int nt);
template void FastPAM<double>::Run(unsigned char opt_method,unsigned int nt);

// FROM NOW ON, INITIALIZATION ALGORTIHMS: form given set of medoids, BUILD (serial and parallel versions) and LAB.

/****************** InitFromPreviousSet ******************************/
template <typename disttype>
void FastPAM<disttype>::InitFromPreviousSet(std::vector<indextype> inlist)
{
 if (inlist.size() != nmed)
 {
     ostringstream errst;
     errst << "Error reading initial medoids file: passed list with " << inlist.size() << " medoids. We expected " << nmed << "\n";
     ParallelpamStop(errst.str());
 }

 // Differently to the R version, the read indices of medoids are 0-based, so we simply store it as we read them. 
 for (size_t t=0; t<size_t(inlist.size()); t++)
  medoids.push_back(inlist[t]);
}

template void FastPAM<float>::InitFromPreviousSet(std::vector<indextype> initmedlist);
template void FastPAM<double>::InitFromPreviousSet(std::vector<indextype> initmedlist);


/********************** BUILD (serial) ****************/
template <typename disttype>
void FastPAM<disttype>::BUILD()
{
    if (DEB & DEBPP)
    {
        std::cout << "Starting BUILD initialization method, serial version\n";
        //std::cout << "WARNING: all successive messages use R-numbering (from 1) for points and medoids. Substract 1 to get the internal C-numbers.\n";
        std::cout << "Looking for medoid 0. ";
        std::cout.flush();
    }

    // Find the first medoid: the point with minimal sum of distances to the rest.
    indextype initial_best=num_obs+1;
    disttype dbest=MAXD;
    disttype sumofrow;
    for (indextype r=0;r<num_obs;r++)
    {
        sumofrow=(disttype)0;
        for (indextype c=0;c<num_obs;c++)
            sumofrow += D->Get(r,c);
        if (sumofrow<dbest)
        {
            dbest = sumofrow;
            initial_best = r;
        }
    }
    if (initial_best>num_obs)
    {
        ParallelpamStop("No best medoid found. Unexpected error.\n");
        return;
    }
    // First, the total distance is the sum of distances of the best to all others
    currentTD=dbest;
    
    if (DEB & DEBPP)
    {
        std::cout << "Medoid 0 found. Point " << initial_best << ". TD=" << std::fixed << currentTD/float(num_obs) << "\n";
        std::cout.flush();
    }
    medoids.push_back(initial_best);
    for (indextype m=1; m<nmed; m++)
     medoids.push_back(NO_CLUSTER);
     
    
    // Initialize the arrays of assignments and closest dissimilarities...
    // To start, all points are in the cluster of the first medoid, so distance to closest medoid (dlcmed) is distance to it.
    for (indextype q=0; q<num_obs; q++)
    {
        nearest[q] = 0;   // The only medoid now is initial_best, at place 0 of medoids' vector
        dnearest[q] = D->Get(q,initial_best);
    }
    
    // The only medoid which is such now is signalled in the array of marks
    ismedoid[initial_best]=true;
    // and the distance to itself is 0. Strictly, this should have been already done in the former loop, but just for clarity...
    dnearest[initial_best]=disttype(0);
    
    // Now, the rest of medoids
   
    disttype d;
    double tdchange,most_negative_tdchange;
    indextype best_up_to_now;
    vector<indextype> reassign_list;
    
    for (indextype nextmed=1; nextmed<nmed; nextmed++)
    {
     if (DEB & DEBPP)
     {
         std::cout << "Looking for medoid " << nextmed << ". ";
         std::cout.flush();
     }
     
     // The maximum decrease in TD is initialized to the lowest possible number
     most_negative_tdchange = numeric_limits<double>::max();
     best_up_to_now = num_obs+1;
     
     // For each point, it is a candidate to be a new medoid...   
     for (indextype cand=0; cand<num_obs; cand++)
     {
         // ...unless it is one of the already found medoids.
         if ( !ismedoid[cand] )
         {
             // The total change in TD is initialized to 0
             tdchange=0.0;
             
             // Now, let's look at each of the other points...
             for (indextype other=0; other<num_obs; other++)
                // .. as said before, 'other' points different from cand and its dissimilarity with its
                // current closest medoid is bigger than the one to me...
                if ( (other!=cand) && ((d=D->Get(cand,other))<dnearest[other]) )
                  // Then, this point should be assigned to the cluster leaded by cand, if it effectively ends up being a medoid...
                  // We will decide on that based on the sum of distances to all these "adherent" points
                  tdchange += double(d-dnearest[other]);
                
             // This is because the distance of the prospective medoid to its closest medoid must be diminished
             // from TD, too, since the cand would be a medoid and the distance to its closest medoid (itself) will become 0.
             // This was not counted before due to the condition (other!=cand)
             tdchange -= dnearest[cand];
             
             // This is to retain the best candidate, i.e.: that which makes tdchange as much negative as possible...
             if ((tdchange<0) && (tdchange<most_negative_tdchange))
             {
                 // We take note of the change in TD to go on comparing. 
                 most_negative_tdchange = tdchange;
                 // Also of who is this candidate..
                 best_up_to_now = cand;
             }
         }
     }
     
     if (best_up_to_now>num_obs)
     {
         ostringstream errst;
         errst << "Error: medoid number " << nextmed << " has not been found. Unexpected error.\n";
         ParallelpamStop(errst.str());
         return;
     }
     // The best candidate is admitted as initial medoid
     medoids[nextmed]=best_up_to_now;
     ismedoid[best_up_to_now]=true;
     dnearest[best_up_to_now]=disttype(0);
     
     // TD is updated:
     if (most_negative_tdchange < -currentTD)
     {
         ParallelpamStop("Error: TD canot become negative.\n");
         return;
     }
     currentTD += most_negative_tdchange;
     
     // Update assignments and closests dissimilarities
     indextype num_updated=0;
     for (indextype q=0; q<num_obs; q++)
     {
         d=D->Get(q,best_up_to_now);
        
         if (d<dnearest[q])
         {
          dnearest[q]=d;
          nearest[q]=nextmed;
          num_updated++;
         }
     }
     
     // The medoid itself is of course in its own cluster, and its dissimilarity with the "closest" (ifself) is obviously 0
     // This has been probably updated in the former loop, but ... 
     nearest[best_up_to_now]=best_up_to_now;
     dnearest[best_up_to_now]=disttype(0);
     
     if (DEB & DEBPP)
     {
         std::cout << "Medoid " << nextmed << " found. Point " << best_up_to_now << ". " << num_updated << " reassigned points. TD=" << std::fixed << currentTD/float(num_obs) << "\n";
         std::cout.flush();
     }
     
    }
    
    if (DEB & DEBPP)
      std::cout << "Current TD: " << std::fixed << currentTD/float(num_obs) << "\n";
}

template void FastPAM<float>::BUILD();
template void FastPAM<double>::BUILD();

/***************** FindFirstMedoidBUILDThread (first thread for parallel BUILD) ****************/
template <typename disttype>
void *FastPAM<disttype>::FindFirstMedoidBUILDThread(void *arg)
{
 unsigned int numthreads = GetNumThreads(arg);
 unsigned int current_thread_num = GetThisThreadNumber(arg);
  
 FastPAM *FPp = GetField(arg,BUILDThread_IO,FPp);
 indextype *bestmed = GetField(arg,BUILDThread_IO,foundmed);
 disttype *DeltaT = GetField(arg,BUILDThread_IO,DeltaT);
 
 indextype start;
 
 // The number of points to be analyzed by each thread, rounded down
 unsigned int num_points_per_th = FPp->num_obs/numthreads;
 // If the division is not exact, some medoids remain 
 unsigned int remaining_points = FPp->num_obs % numthreads;
 // If this is the case, each one of the remaining medoids will be assigned to each of the first threads.
 if ( remaining_points != 0 && current_thread_num < remaining_points)
  num_points_per_th++;
 
 start = current_thread_num*num_points_per_th;
 // The latests threads start higher, since remaining_medoids have been assigned to the initial threads.
 if (current_thread_num >= remaining_points)
  start += (FPp->num_obs % numthreads);
 
 // The end of the points assigned to this thread. The point with this number will NOT be managed by this thread;
 // it is the first point of the next one. Notice that the loops will run from point=start to point<end, not to point<=end.
 indextype end = start + num_points_per_th;
 if (end>FPp->num_obs)
  end=FPp->num_obs;
  
 // Find the first best medoid among the points assigned to this thread: the point with minimal sum of distances to _all_ others
 indextype initial_best=FPp->num_obs+1;
 disttype dbest=MAXD;
 disttype sumofrow;
 for (indextype r=start;r<end;r++)
 {
    sumofrow=(disttype)0;
    for (indextype c=0;c<FPp->num_obs;c++)
      sumofrow += FPp->D->Get(r,c);
    if (sumofrow<dbest)
    {
        dbest = sumofrow;
        initial_best = r;
    }
 }
 
 *bestmed = initial_best;
 *DeltaT = dbest;
 
 pthread_exit(nullptr); 
}

template void *FastPAM<float>::FindFirstMedoidBUILDThread(void *arg);
template void *FastPAM<double>::FindFirstMedoidBUILDThread(void *arg);

/***************** FindSuccessiveMedoidBUILDThread (second thread for parallel BUILD) ****************/
template <typename disttype>
void *FastPAM<disttype>::FindSuccessiveMedoidBUILDThread(void *arg)
{
 unsigned int numthreads = GetNumThreads(arg);
 unsigned int current_thread_num = GetThisThreadNumber(arg);
  
 FastPAM *FPp = GetField(arg,BUILDThread_IO,FPp);
 indextype *returnedmed = GetField(arg,BUILDThread_IO,foundmed);
 disttype *returnedchange = GetField(arg,BUILDThread_IO,DeltaT);
 
 indextype start;
 
 // The number of points to be analyzed by each thread, rounded down
 unsigned int num_points_per_th = FPp->num_obs/numthreads;
 // If the division is not exact, some medoids remain 
 unsigned int remaining_points = FPp->num_obs % numthreads;
 // If this is the case, each one of the remaining medoids will be assigned to each of the first threads.
 if ( remaining_points != 0 && current_thread_num < remaining_points)
  num_points_per_th++;
 
 start = current_thread_num*num_points_per_th;
 // The latests threads start higher, since remaining_medoids have been assigned to the initial threads.
 if (current_thread_num >= remaining_points)
  start += (FPp->num_obs % numthreads);
 
 // The end of the points assigned to this thread. The point with this number will NOT be managed by this thread;
 // it is the first point of the next one. Notice that the loops will run from point=start to point<end, not to point<=end.
 indextype end = start + num_points_per_th;
 if (end>FPp->num_obs)
  end=FPp->num_obs;
  
 disttype d;
 // The maximum decrease in TD is initialized to the lowest possible number
 disttype most_negative_tdchange = numeric_limits<double>::max();
 disttype tdchange;
 indextype best_up_to_now = FPp->num_obs+1;
 // For each point, it is a candidate to be a new medoid...   
 for (indextype cand=start; cand<end; cand++)
 {
  // ...unless it is one of the already found medoids.
  if ( !(FPp->ismedoid)[cand] )
  {
   // The total change in TD is initialized to 0
   tdchange=0.0;
             
   // Now, let's look at each of the other points...
   for (indextype other=0; other<FPp->num_obs; other++)
    // .. as said before, 'other' points different from cand and its dissimilarity with its
    // current closest medoid is bigger than the one to me...
    if ( (other!=cand) && ((d=FPp->D->Get(cand,other))<(FPp->dnearest)[other]) )
     // Then, this point should be assigned to the cluster leaded by cand, if it effectively ends up being a medoid...
     // We will decide on that based on the sum of distances to all these "adherent" points
     tdchange += double(d-FPp->dnearest[other]);
                
   // This is because the distance of the prospective medoid to its closest medoid must be diminished
   // from TD, too, since the cand would be a medoid and the distance to its closest medoid (itself) will become 0.
   // This was not counted before due to the condition (other!=cand)
   tdchange -= (FPp->dnearest)[cand];
             
   // This is to retain the best candidate, i.e.: that which makes tdchange as much negative as possible...
   if ((tdchange<0) && (tdchange<most_negative_tdchange))
   {
    // We take note of the change in TD to go on comparing. 
    most_negative_tdchange = tdchange;
    // Also of who is this candidate..
    best_up_to_now = cand;
   }
  }
 }
 
 *returnedmed=best_up_to_now;
 *returnedchange=most_negative_tdchange;
 
 pthread_exit(nullptr);    
}

template void *FastPAM<float>::FindSuccessiveMedoidBUILDThread(void *arg);
template void *FastPAM<double>::FindSuccessiveMedoidBUILDThread(void *arg);

/***************** ParBUILD (BUILD in parallel version) ***********************/
template <typename disttype>
void FastPAM<disttype>::ParBUILD(unsigned int nt)
{
    if (DEB & DEBPP)
    {
        std::cout << "Starting BUILD initialization method, parallel version with " << nt << " threads.\n";
        //std::cout << "WARNING: all successive messages use R-numbering (from 1) for points and medoids. Substract 1 to get the internal C-numbers.\n";
        std::cout << "Looking for medoid 0. ";
        std::cout.flush();
    }
    
    disttype dbest=MAXD;
    BUILDThread_IO *BUILDargs = new BUILDThread_IO [nt];
    
    indextype initial_best=num_obs+1;
    {
     indextype *initial_bestTh = new indextype [nt];
     disttype *dbestTh = new disttype [nt];
    
     for (unsigned int t=0; t<nt; t++)
     {
        BUILDargs[t].FPp = this;
        BUILDargs[t].foundmed = &initial_bestTh[t];
        BUILDargs[t].DeltaT = &dbestTh[t];
     }
       
     CreateAndRunThreadsWithDifferentArgs(nt,FindFirstMedoidBUILDThread,BUILDargs,sizeof(BUILDThread_IO));

     for (unsigned int t=0; t<nt; t++)
        if (dbestTh[t] < dbest)
        {
            dbest = dbestTh[t];
            initial_best = initial_bestTh[t];
        }
        
     if (initial_best>num_obs)
       ParallelpamStop("Error: no best medoid found. Unexpected error.\n");
    
     delete[] initial_bestTh;
     delete[] dbestTh;
    }
    
    // First, the total distance is the sum of distances of the best to all others
    currentTD=dbest;
    
    medoids.resize(nmed,num_obs+1);
    medoids[0]=initial_best;
    
    if (DEB & DEBPP)
    {
        std::cout << "Medoid 0 found. Point " << initial_best << ". TD=" << std::fixed << currentTD/float(num_obs) << "\n";
        std::cout.flush();
    }

    // Initialize the arrays of assignments and closest dissimilarities...
    // To start, all points are in the cluster of the first medoid, so distance to closest medoid (dlcmed) is distance to it.
    for (indextype r=0;r<num_obs;r++)
    {
        nearest[r] = 0;   // The only medoid now is initial_best, at place 0 of medoids' vector
        dnearest[r] = D->Get(initial_best,r);
    }
    
    ismedoid[initial_best]=true;
    dnearest[initial_best]=disttype(0);
    
    // Now, the rest of medoids
   
    disttype d;
    double most_negative_tdchange;
    indextype best_up_to_now;
    vector<indextype> reassign_list;
    
    for (indextype nextmed=1; nextmed<nmed; nextmed++)
    {
     if (DEB & DEBPP)
     {
         std::cout << "Looking for medoid " << nextmed << ". ";
         std::cout.flush();
     }
     
     {
      indextype *best_up_to_nowTh = new indextype [nt];
      disttype *most_negative_tdchangeTh = new disttype [nt];
      for (unsigned int t=0; t<nt; t++)
      {
         BUILDargs[t].FPp=this;
         BUILDargs[t].foundmed=&best_up_to_nowTh[t];
         BUILDargs[t].DeltaT=&most_negative_tdchangeTh[t];
      }
    
      CreateAndRunThreadsWithDifferentArgs(nt,FindSuccessiveMedoidBUILDThread,BUILDargs,sizeof(BUILDThread_IO));
     
      // The maximum decrease in TD is initialized to the lowest possible number
      most_negative_tdchange = numeric_limits<double>::max();
      best_up_to_now = num_obs+1;
     
      for (unsigned int t=0; t<nt; t++)
        if (most_negative_tdchangeTh[t]<most_negative_tdchange)
        {
            most_negative_tdchange = most_negative_tdchangeTh[t];
            best_up_to_now = best_up_to_nowTh[t];
        }
        
      delete[] best_up_to_nowTh;
      delete[] most_negative_tdchangeTh;
     }
     
     if (best_up_to_now>num_obs)
     {
     	ostringstream errst;
        errst << "Error: medoid number " << nextmed << " has not been found. Unexpected error.\n";
        ParallelpamStop(errst.str());
        return;
     }
     // The best candidate is admitted as initial medoid
     medoids[nextmed]=best_up_to_now;
     ismedoid[best_up_to_now]=true;
     dnearest[best_up_to_now]=disttype(0);
     
     // TD is updated:
     if (most_negative_tdchange < -currentTD)
     {
         ParallelpamStop("Error: TD canot become negative.\n");
         return;
     }
     currentTD += most_negative_tdchange;
     
     // Update assignations and closests dissimilarities
     indextype num_updated=0;
     for (indextype q=0;q<num_obs;q++)
     {
         d=D->Get(q,best_up_to_now);
        
         if (d<dnearest[q])
         {
          dnearest[q]=d;
          nearest[q]=nextmed;
          num_updated++;
         }
     }
     
     // The medoid itself is of course in its own cluster, and its dissimilarity with the "closest" (ifself) is obviously 0
     // This has been probably updated in the former loop, but ... 
     nearest[best_up_to_now]=best_up_to_now;
     dnearest[best_up_to_now]=disttype(0);
     
     if (DEB & DEBPP)
     {
         std::cout << "Medoid " << nextmed << " found. Point " << best_up_to_now << ". " << num_updated << " reassigned points. TD=" << std::fixed << currentTD/float(num_obs) << "\n";
         std::cout.flush();
     }
    }
    
    if (DEB & DEBPP)
     std::cout << "Current TD: " << std::fixed << currentTD/float(num_obs) << "\n";
       
    delete[] BUILDargs;
}

template void FastPAM<float>::ParBUILD(unsigned int nt);
template void FastPAM<double>::ParBUILD(unsigned int nt);

/*********************** LAB (serial version) **********************************/
template <typename disttype>
void FastPAM<disttype>::LAB()
{
    if (DEB & DEBPP)
    {
        std::cout << "Starting LAB initialization method, serial version.\n";
        //Rcpp::Rcout << "WARNING: all successive messages use R-numbering (from 1) for points and medoids. Substract 1 to get the internal C-numbers.\n";
        std::cout << "Looking for medoid 0. ";
        std::cout.flush();
    }
    
    // First, we get a subsample
    size_t samplesize = 20 + 2*ceil(sqrt(double(num_obs)));
    // This check is to prevent the special case of a really low number of observations...
    if (samplesize>num_obs)
        samplesize=num_obs;
    
    // A random sample is chosen:
    vector<indextype> S=randomSample(samplesize,num_obs);   
    
    // Find the first medoid in this sample: the point with minimal sum of distances to the rest.
    
    indextype initial_best;
    disttype dbest=MAXD;
    disttype sumofdist;
    for (indextype r=0;r<S.size();r++)
    {
        sumofdist=disttype(0);
        for (unsigned c=0;c<S.size();c++)
          if (r!=c)
              sumofdist += D->Get(S[r],S[c]);
          
        if (sumofdist<dbest)
        {
            dbest = sumofdist;
            initial_best = S[r];
        }
    }
    
    medoids.clear();
    medoids.push_back(initial_best);
    
    // Initialize the arrays of assignments and closest dissimilarities...
    // To start, all points are in the cluster of the first medoid, so distance to closest medoid (dlcmed) is to it.
    // Also, the total distance is the sum of distances of the best of this subsample to all others
    currentTD=disttype(0);
    for (indextype r=0;r<num_obs;r++)
    {
        nearest[r] = 0;        // The only medoid now is initial_best, at place 0 of medoids' vector
        dnearest[r] = D->Get(initial_best,r);
        currentTD += dnearest[r];
    }
    
    if (DEB & DEBPP)
    {
        std::cout << "Medoid 0 found. Point " << initial_best << ". TD=" << std::fixed << currentTD/float(num_obs) << "\n";
        std::cout.flush();
    } 
    
    // The array of marks to check easily if a point has been found as a medoid is updated with the first medoid
    ismedoid[initial_best]=true;
    dnearest[initial_best]=disttype(0);
    
    // Now, the rest of medoids
    disttype DeltaTD,DeltaTDstar,delta;
    indextype xstar;
    vector<indextype> reassign_list;
    
    for (indextype nextmed=1; nextmed<nmed; nextmed++)
    {
     if (DEB & DEBPP)
     {
         std::cout << "Looking for medoid " << nextmed << ". ";
         std::cout.flush();
     }
     
     DeltaTDstar = MAXD;
     xstar = num_obs+1;
     
     S = randomSampleExc(samplesize,num_obs,ismedoid);
     
     for (indextype j=0; j<samplesize; j++)
     {
         
         DeltaTD=0;
         
         for (indextype xz=0; xz < samplesize; xz++)
         {
          if (S[xz] != S[j])
          {
           delta = D->Get(S[xz],S[j])-dnearest[S[xz]];
           if (delta<0)
            DeltaTD += delta;   
          }
         }
         
         if (DeltaTD < DeltaTDstar)
         {
             DeltaTDstar = DeltaTD;
             xstar = S[j];
         }
      }    
      
      medoids.push_back(xstar);
      ismedoid[xstar]=true;
      
      // Update assignations and closests dissimilarities
      indextype num_updated=0;
      disttype d,dummy;
      for (indextype q=0;q<num_obs;q++)
      {
         d=D->Get(q,xstar);
         if (d<dnearest[q])
         {
             dummy=dnearest[q];
             dnearest[q]=d;
             currentTD -= dummy;
             currentTD += d;
             nearest[q]=medoids.size()-1;   // The new medoid has just appended to the medoid's vector, so it is at the last position.
             num_updated++;
         }
      }
      if (currentTD<0)
      {
          ParallelpamStop("Error: TD cannot be negative.\n");
          return;
      }
      // The medoid itself is of course in its own cluster, and its dissimilarity with the "closest" (ifself) is obviously 0
      // This has been probably updated in the former loop, but ... 
      nearest[xstar]=medoids.size()-1;
      dnearest[xstar]=disttype(0);
     
      if (DEB & DEBPP)
      {
         std::cout << "Medoid " << nextmed << " found. Point " << xstar << ". " << num_updated << " reassigned points. TD=" << std::fixed << currentTD/float(num_obs) << "\n";
         std::cout.flush();
      }
    }
    if (DEB & DEBPP)
     std::cout << "Current TD: " << std::fixed << currentTD/float(num_obs) << "\n";
}

template void FastPAM<float>::LAB();
template void FastPAM<double>::LAB();

// FROM HERE, ONE OF THE ALGORITHMS FOR THE OPTIMIZATION PHASE, FastPAM1, in serial and parallel version

/**************************** RunImprovedFastPAM1 (optimization phase, serial version) *****************/
// This function closely follows the notation in the original work (Schubert and Rousseauw 2021)
// Comments with Ln refer to line n of Algortihm 3 in such paper.
template <typename disttype>
void FastPAM<disttype>::RunImprovedFastPAM1()
{
 if (DEB & DEBPP)
 {
  std::cout << "Starting improved FastPAM1 method in serial implementation...\n";
  //Rcpp::Rcout << "WARNING: all successive messages use R-numbering (from 1) for points and medoids. Substract 1 to get the internal C-numbers.\n";
  std::cout.flush();
 }
 
 // dsecond is to be filled in advance, mostly as cache.
 FillSecond();
 
 // The threshold that will stop the algorithm if TD changes less than this value at any iteration.
 disttype tol_limit=currentTD*tlimit;
 
 // Now, local variables used in the paper's algorithm. Ths star (*) is translated as st so m* will be named mst
 disttype *DeltaTDminusm = new disttype [nmed];
 disttype DeltaTDplusxc,DeltaTDst,d0j;
 disttype *DeltaTD = new disttype [nmed];
 
 // I take these as number of the point and number of the medoid
 indextype xst,mst;
 // Nevertheless, these are indexes in the vector of medoids. imst is not explictly named in the original work.
 indextype i,imst;
  
 unsigned int iteration=0;
 bool out=false;                         // Used to leave in special case of no TD improvement, i.e., no better solucion exist.
 do                                                          // L2
 {
  if (DEB & DEBPP)
  {
   std::cout << "Iteration " << iteration << ". ";
   std::cout.flush();
  }
  
  for (indextype m=0; m<nmed; m++)                            // L3
  {
    DeltaTDminusm[m]=disttype(0);
    for (indextype q=0; q<num_obs; q++)
     if (nearest[q]==m)
      DeltaTDminusm[m] += (dsecond[q]-dnearest[q]);   // Since dsecond[q] is always > dnearest[q], DeltaTDminus[m] will always be possitive for all m
  }
 
  DeltaTDst = disttype(0);                                     // L4
  mst = num_obs+1;                       // This is our 'null'
  xst = num_obs+1;                       // Same here...
  i = nmed+1;                            // Just as a check to be tested outside the next loop, to see that i has been changed.
  imst = nmed+1;
   
  for (indextype xc=0; xc<num_obs; xc++)                        // L5
  {
    if (!ismedoid[xc])                                          // L5
    {
        
       for (indextype m=0; m<nmed; m++)                         // L6   DeltaTD is initialized to a possitive value for all m, since DeltaTDminusm[m] is possitive
           DeltaTD[m] = DeltaTDminusm[m];
       DeltaTDplusxc = disttype(0);                             // L7
        
       for (indextype x0=0; x0<num_obs; x0++)                   // L8
       {
         d0j = D->Get(x0,xc);                                   // L9
         if (d0j < dnearest[x0])                                // L10
         {
            DeltaTDplusxc += (d0j - dnearest[x0]);              // L11  Since d0j here is smaller then dnearest[x0], DeltaTDplusxc is reduced and becomes more and more negative
            DeltaTD[nearest[x0]] += (dnearest[x0]-dsecond[x0]); // L12     and DeltaTD is reduced, since dnearest[x0] < dsecond[x0]
         }
         else
            if (d0j<dsecond[x0])                                // L13
              DeltaTD[nearest[x0]] += (d0j-dsecond[x0]);        // L14   Here DeltaTD is reduced, too, since in this part of the conditional d0j < dsecond[x0]
       }  // for (indextype x0=...
        
       disttype ddummy=MAXD;                                   // L15
       i=nmed+1;
       for (indextype m=0; m<nmed; m++)
        if (DeltaTD[m]<ddummy)
        {
             ddummy = DeltaTD[m];
             i = m;
        }
         
       DeltaTD[i] += DeltaTDplusxc;                       // L16
        
       if (DeltaTD[i]<DeltaTDst)                          // L17     DeltaTDst is the best improvement up to now. If this loop turn has made it improve...
       {
           DeltaTDst = DeltaTD[i];                        // ... update the best improvement and take note of the swap that has provoked it:
           mst = medoids[i];                              // The number of the point which is the medoid that should be left out
           xst = xc;                                      // The number of the point which is to become the new medoid
           imst = i;                                      // The index in the array of medoids which has the medoid to be swapped
       }
    }           // if (!ismedoid)...
  }   // for (indextype xc=...
    
  if (DeltaTDst>=disttype(0))                        // L18
  {
      if (DEB & DEBPP)
        std::cout << "   Exiting, since DeltaTDst is " << std::fixed << DeltaTDst/float(num_obs) << ". Final value of TD is " << std::fixed << currentTD/float(num_obs) << "\n";
      out=true;
      break;
  }
   
  if ((DEB & DEBPP) && (!out) && (imst<nmed))
   std::cout << "Medoid at place " << imst << " (point " << medoids[imst] << ") swapped with point " << xst << "; ";
  
  // Could imst be left unchanged by the loop? Not except by error. Anyway, let's check it
  if ((imst<nmed) && (!out))
  {   
   SwapRolesAndUpdate(mst,xst,imst);                                // L19-20
  
   currentTD += DeltaTDst;                                          // L21
   
   if (DEB & DEBPP)
    std::cout << "TD-change=" << std::fixed << DeltaTDst/float(num_obs) << "; TD=" << std::fixed << currentTD/float(num_obs) << ". " << current_npch << " reassigned points.\n";
  }
  else
  {
     if (DEB & DEBPP)
     {
      std::cout << "   No exchange of medoid/point found which can improve result. Exact result found?\n";
      std::cout << "   Last TD change has been " << std::fixed << DeltaTDst/float(num_obs) << "\n";
      if (imst>nmed)
          std::cout << "Best medoid has not been updated.\n";
      else
          std::cout << "Nevertheless, best medoid has been updated to " << medoids[i] << ". ????\n";
     }
     out=true;
  }
  
  iteration++;
    
  // This is the only point (apart from the messages in the screen) in which TD is converted from raw sum to sum per point.
  TDkeep.push_back(currentTD/float(num_obs));
  NpointsChangekeep.push_back(current_npch); 
 }
 while ((fabs(DeltaTDst)>tol_limit) && (iteration<maxiter) && (!out));   // fabs because DeltaTDst is negative...
 num_iterations_in_opt=(iteration>0) ? iteration-1 : 0;
 
 delete[] DeltaTDminusm;
 delete[] DeltaTD;
} 

template void FastPAM<float>::RunImprovedFastPAM1();
template void FastPAM<double>::RunImprovedFastPAM1();

/**************** FastPAM1InternalThread (thread for PAM optimization phase) ********************/
// This thread is called by RunParallelImprovedFastPAM1. See comments there on original source and notation.
template <typename disttype>
void *FastPAM<disttype>::FastPAM1InternalThread(void *arg)
{
 unsigned int numthreads = GetNumThreads(arg);
 unsigned int current_thread_num = GetThisThreadNumber(arg);
  
 FastPAM *FPp = GetField(arg,FastPAM1Thread_IO,FPp);
 indextype *mst = GetField(arg,FastPAM1Thread_IO,mst);
 indextype *xst = GetField(arg,FastPAM1Thread_IO,xst);
 indextype *imst = GetField(arg,FastPAM1Thread_IO,imst);
 disttype *DeltaTDst = GetField(arg,FastPAM1Thread_IO,DeltaTDst);
 disttype *DeltaTDminusm = GetField(arg,FastPAM1Thread_IO,DeltaTDminusm);
 
 indextype start;
 
 // The number of points to be analyzed by each thread, rounded down
 unsigned int num_points_per_th = FPp->num_obs/numthreads;
 // If the division is not exact, some medoids remain 
 unsigned int remaining_points = FPp->num_obs % numthreads;
 // If this is the case, each one of the remaining medoids will be assigned to each of the first threads.
 if ( remaining_points != 0 && current_thread_num < remaining_points)
  num_points_per_th++;
 
 start = current_thread_num*num_points_per_th;
 // The latests threads start higher, since remaining_medoids have been assigned to the initial threads.
 if (current_thread_num >= remaining_points)
  start += (FPp->num_obs % numthreads);
 
 // The end of the points assigned to this thread. The point with this number will NOT be managed by this thread;
 // it is the first point of the next one. Notice that the loops will run from point=start to point<end, not to point<=end.
 indextype end = start + num_points_per_th;
 if (end>FPp->num_obs)
  end=FPp->num_obs;
   
 for (indextype xc=start; xc<end; xc++)                                            // L5
 {
  if (!(FPp->ismedoid)[xc])                                                        // L5
  {
    disttype *DeltaTD = new disttype [FPp->nmed];
    for (indextype m=0; m<FPp->nmed; m++)                                          // L6   DeltaTD is initialized to a possitive value for all m, since DeltaTDminusm[m] is possitive
      DeltaTD[m] = DeltaTDminusm[m];
    disttype DeltaTDplusxc = disttype(0);                                          // L7
            
    for (indextype x0=0; x0<FPp->num_obs; x0++)                                    // L8
    {
       disttype d0j = FPp->D->Get(x0,xc);                                          // L9
       if (d0j < (FPp->dnearest)[x0])                                              // L10
       {
          DeltaTDplusxc += (d0j - (FPp->dnearest)[x0]);                            // L11  Since d0j here is smaller then dnearest[x0], DeltaTDplusxc is reduced and become more and more negative
          DeltaTD[(FPp->nearest)[x0]] += ((FPp->dnearest)[x0]-(FPp->dsecond)[x0]); // L12     and DeltaTD is reduced, since dnearest[x0] < dsecond[x0]
       }
       else
          if (d0j<(FPp->dsecond)[x0])                                              // L13
            DeltaTD[(FPp->nearest)[x0]] += (d0j-(FPp->dsecond)[x0]);               // L14   Here DeltaTD is reduced, too, since in this part of the conditional d0j < dsecond[x0]
    }
        
    disttype ddummy=MAXD;                                                          // L15
    indextype i=FPp->nmed+1;
    for (indextype m=0; m<FPp->nmed; m++)
      if (DeltaTD[m]<ddummy)
      {
         ddummy = DeltaTD[m];
         i = m;
      }
      
    if (i>FPp->nmed)
    {
       // This is just a check that should never be true, and moreover, would provoke a crash inside a thread. To be elliminated in last version.
       std::ostringstream errst;
       errst << "In loop with xc=" << xc << ": no closest medoid found. Unexpected error.\n";
       ParallelpamStop(errst.str());
    }
         
    DeltaTD[i] += DeltaTDplusxc;                                                    // L16
        
    if (DeltaTD[i]<(*DeltaTDst))                                                    // L17
    {
       *DeltaTDst = DeltaTD[i];
       *mst = (FPp->medoids)[i];
       *xst = xc;
       *imst = i;
    } 
    
    delete [] DeltaTD;  
  }  // if (!ismedoid)...
 }  // for (indextype x0...
 
 pthread_exit(nullptr);
}

template void *FastPAM<float>::FastPAM1InternalThread(void *arg);
template void *FastPAM<double>::FastPAM1InternalThread(void *arg);

/**************** RunParallelImprovedFastPAM1 (optimization, parallel version) ********************/
// This function closely follows the notation in the original work (Schubert and Rousseauw 2021)
// Comments with Ln refer to line n of Algortihm 3 in such paper.
template <typename disttype>
void FastPAM<disttype>::RunParallelImprovedFastPAM1(unsigned int nt)
{
 if (DEB & DEBPP)
 {
  std::cout << "Starting improved FastPAM1 method in parallel implementation with " << nt << " threads.\n";
  //Rcpp::Rcout << "WARNING: all successive messages use R-numbering (from 1) for points and medoids. Substract 1 to get the internal C-numbers.\n";
  std::cout.flush();
 }
 
 // dsecond is to be filled in advance, mostly as cache.
 FillSecond();
 
 // The threshold that will stop the algorithm if TD changes less than this value at any iteration.
 disttype tol_limit=currentTD*tlimit;
 
 // Now, local variables used in the paper's algorithm. Ths star (*) is translated as st so m* will be named mst
 disttype *DeltaTDstPerTh = new disttype [nt];
 disttype DeltaTDst;
 
 // I take these as number of the point and number of the medoid respectively
 indextype *xstPerTh = new indextype [nt];
 indextype *mstPerTh = new indextype [nt];
 indextype xst,mst;
 // Nevertheless, these are indexes in the vector of medoids. imst is not explictly named in the original work.
 indextype *imstPerTh = new indextype [nt];
 indextype imst;
 
 struct FastPAM1Thread_IO *FastPAM1args = new struct FastPAM1Thread_IO [nt];

 unsigned int iteration=0;
 bool out=false;            // Used to leave in special case of no TD improvement, i.e., no better solucion exist.
 do                                                        // L2
 {
  if (DEB & DEBPP)
  {
   std::cout << "Iteration " << iteration << ". ";
   std::cout.flush();
  }
  
  disttype *DeltaTDminusm = new disttype [nmed];
  for (indextype m=0; m<nmed; m++)                         // L3
  {
    DeltaTDminusm[m]=disttype(0);
    // Since dsecond[q] is always > dnearest[q], DeltaTDminus[m] will always be possitive for all m
    for (indextype q=0; q<num_obs; q++)
     if (nearest[q]==m)
           DeltaTDminusm[m] += (dsecond[q]-dnearest[q]);
  }
 
  // Filling of arguments to each thread
  for (unsigned int t=0; t<nt; t++)
  {
   DeltaTDstPerTh[t] = disttype(0);                        // L4
   mstPerTh[t] = num_obs+1;                                // This is our 'null'
   xstPerTh[t] = num_obs+1;                                // Same here...
   imstPerTh[t] = nmed+1;
   FastPAM1args[t].FPp = this;
   FastPAM1args[t].mst = &mstPerTh[t];
   FastPAM1args[t].xst = &xstPerTh[t];
   FastPAM1args[t].imst = &imstPerTh[t];
   FastPAM1args[t].DeltaTDst = &DeltaTDstPerTh[t];
   FastPAM1args[t].DeltaTDminusm = DeltaTDminusm;         // Notice that here we initialize a pointer as such, not as the indirection of a variable.
                                                          // It is the same value for all threads, since inside the thread this is used as a read-only variable to initialize other variables.                                                  
  }
  
                                                          // Lines 5 to 17 are inside each thread
  CreateAndRunThreadsWithDifferentArgs(nt,FastPAM1InternalThread,FastPAM1args,sizeof(FastPAM1Thread_IO));
  
  // Now, let's fuse the results of each thread keeping only the minimum value of DeltaTDst, which indicates
  // the right values of the associate variables mst, imst and xst.
  DeltaTDst = MAXD;                     
  mst = num_obs+1;                       
  xst = num_obs+1;                                                 
  imst = nmed+1;
  for (unsigned int t=0; t<nt; t++)
  {
   if (DeltaTDstPerTh[t]<DeltaTDst)
   {
    DeltaTDst = DeltaTDstPerTh[t];
    mst = mstPerTh[t];
    xst = xstPerTh[t];
    imst = imstPerTh[t];
   }
  }
   
  if (DeltaTDst>=disttype(0))                               // L18
  {
     if (DEB & DEBPP)
       std::cout << "   Exiting, since DeltaTDst is " << std::fixed << DeltaTDst/float(num_obs) << ". Final value of TD is " << std::fixed << currentTD/float(num_obs) << "\n";
     out=true;
     break;
  }
   
  if ((DEB & DEBPP) && (!out) && (imst<nmed))
    std::cout << "Medoid at place " << imst << " (point " << medoids[imst] << ") swapped with point " << xst << "; ";
  
  // Could imst be left unchanged by the loop? Not except by error. Anyway, let's check it
  if ((imst<nmed) && (!out))
  {   
    SwapRolesAndUpdate(mst,xst,imst);                        // L19-20
  
    currentTD += DeltaTDst;                                  // L21
   
    if (DEB & DEBPP)
     std::cout << "TD-change=" << std::fixed << DeltaTDst/float(num_obs) << "; TD=" << std::fixed << currentTD/float(num_obs) << ". " << current_npch << " reassigned points.\n";
  }
  else
  {
    if (DEB & DEBPP)
    {
      std::cout << "   No exchange of medoid/point found which can improve result. Exact result found?\n";
      std::cout << "   Last TD change has been " << std::fixed << DeltaTDst/float(num_obs) << "\n";
      if (imst>nmed)
        std::cout << "Best medoid has not been updated.\n";
      else
        std::cout << "Nevertheless, best medoid has been updated to " << medoids[imst] << ". ????\n";
    }
    out=true;
  }
  
  iteration++;
    
  // This is the only point (apart from the messages in the screen) in which TD is converted from raw sum to sum per point.
  TDkeep.push_back(currentTD/float(num_obs));
  NpointsChangekeep.push_back(current_npch); 
  
  delete[] DeltaTDminusm;
 }
 while ((fabs(DeltaTDst)>tol_limit) && (iteration<maxiter) && (!out));   // fabs because DeltaTDst is negative...
 num_iterations_in_opt=(iteration>0) ? iteration-1 : 0;
 
 delete[] DeltaTDstPerTh;
 delete[] imstPerTh;
 delete[] xstPerTh;
 delete[] mstPerTh;
 delete[] FastPAM1args;
} 

template void FastPAM<float>::RunParallelImprovedFastPAM1(unsigned int nt);
template void FastPAM<double>::RunParallelImprovedFastPAM1(unsigned int nt);

/*********************************************************************
 * FROM HERE, NEW VARIANT OF FASTPAM1, The MultiBranch version
 *********************************************************************/
/* *******************************************************************
 * Multibranch, serial implementation
 *********************************************************************/
template <typename disttype>
void FastPAM<disttype>::ExploreBranches(disttype *DeltaTDminusm,disttype *DeltaTD,vector<exchange> &xcg)
{
  // Now, local variables used in the paper's algorithm. Ths star (*) is translated as st so m* will be named mst
  disttype DeltaTDplusxc,d0j,DeltaTDst;
  indextype mst,xst;
  // i is an index in the vector of medoids.
  indextype i,imst;

  for (indextype m=0; m<nmed; m++)                            // L3
  {
    DeltaTDminusm[m]=disttype(0);
    for (indextype q=0; q<num_obs; q++)
     if (nearest[q]==m)
      DeltaTDminusm[m] += (dsecond[q]-dnearest[q]);   // Since dsecond[q] is always > dnearest[q], DeltaTDminus[m] will always be possitive for all m
  }

  DeltaTDst = disttype(0);                                     // L4
  mst = num_obs+1;                       // This is our 'null'
  xst = num_obs+1;                       // Same here...
  i = nmed+1;                                // Just as a check to be tested outside the next loop, to see that i has been changed.
  imst = nmed+1;

  size_t B=xcg.size(),k;
  vector<disttype> vmin(B,MAXD);
  bool done_exchange=false;
  for (indextype xc=0; xc<num_obs; xc++)                        // L5
  {
    done_exchange=false;
    if (!ismedoid[xc])                                          // L5
    {

       for (indextype m=0; m<nmed; m++)                         // L6   DeltaTD is initialized to a possitive value for all m, since DeltaTDminusm[m] is possitive
           DeltaTD[m] = DeltaTDminusm[m];
       DeltaTDplusxc = disttype(0);                             // L7

       for (indextype x0=0; x0<num_obs; x0++)                   // L8
       {
         d0j = D->Get(x0,xc);                                   // L9
         if (d0j < dnearest[x0])                                // L10
         {
            DeltaTDplusxc += (d0j - dnearest[x0]);              // L11  Since d0j here is smaller then dnearest[x0], DeltaTDplusxc is reduced (becomes more negative)
            DeltaTD[nearest[x0]] += (dnearest[x0]-dsecond[x0]); // L12  and DeltaTD is reduced, since dnearest[x0] < dsecond[x0]
         }
         else
          if (d0j<dsecond[x0])                                // L13
            DeltaTD[nearest[x0]] += (d0j-dsecond[x0]);        // L14   Here DeltaTD is reduced, too, since in this part of the conditional d0j < dsecond[x0]
       }  // for (indextype x0=...

       disttype ddummy=MAXD;                                   // L15
       i=nmed+1;
       for (indextype m=0; m<nmed; m++)
        if (DeltaTD[m]<ddummy)
        {
            ddummy = DeltaTD[m];
            i = m;
        }

       DeltaTD[i] += DeltaTDplusxc;                       // L16

       if (DeltaTD[i]<DeltaTDst)                      // L17   DeltaTDst is the best improvement up to now. If this loop turn has made it improve...
       {
           DeltaTDst = DeltaTD[i];                    // ... update the best improvement and take note of the swap that has provoked it:
           mst = medoids[i];                          // The number of the point which is the medoid that should be left out
           xst = xc;                                  // The number of the point which is to become the new medoid
           imst = i;                                  // The index in the array of medoids which has the medoid to be swapped
           done_exchange=true;
       }
    }           // if (!ismedoid)...

    if (done_exchange)
    {
     k=0;
     while ( (k<B) && (vmin[k]<DeltaTDst) )
      k++;
     if (k<B)
     {
      for (size_t l=B-1;l>k;l--)
      {
       vmin[l]=vmin[l-1];
       xcg[l]=xcg[l-1];
      }
      vmin[k]=DeltaTDst;
      xcg[k].DeltaTDst=DeltaTDst;
      xcg[k].mst=mst;
      xcg[k].xst=xst;
      xcg[k].imst=imst;
     }
    }
  }   // for (indextype xc=...
}

template void FastPAM<float>::ExploreBranches(float *DeltaTDminusm,float *DeltaTD,vector<exchange> &xcg);
template void FastPAM<double>::ExploreBranches(double *DeltaTDminusm,double *DeltaTD,vector<exchange> &xcg);

/*********************************************************************
 * Multibranch, parallel implementation
 *********************************************************************/
/**************** ExploreBranchesInternalThread (thread for PAM optimization phase in my variant) ********************/
// This thread is called by ExploreBranchesParallel.
/*
template <typename disttype>
void *FastPAM<disttype>::ExploreBranchesInternalThread(void *arg)
{
 unsigned int numthreads = GetNumThreads(arg);
 unsigned int current_thread_num = GetThisThreadNumber(arg);

 FastPAM *FPp = GetField(arg,ExploreBranchesThread_IO,FPp);
 indextype *mst = GetField(arg,ExploreBranchesThread_IO,mst);
 indextype *xst = GetField(arg,ExploreBranchesThread_IO,xst);
 indextype *imst = GetField(arg,ExploreBranchesThread_IO,imst);
 disttype *DeltaTDst = GetField(arg,ExploreBranchesThread_IO,DeltaTDst);
 disttype *DeltaTDminusm = GetField(arg,ExploreBranchesThread_IO,DeltaTDminusm);
 exchange *xcg = GetField(arg,ExploreBranchesThread_IO,xcg);

 indextype start;

 // The number of points to be analyzed by each thread, rounded down
 unsigned int num_points_per_th = FPp->num_obs/numthreads;
 // If the division is not exact, some medoids remain
 unsigned int remaining_points = FPp->num_obs % numthreads;
 // If this is the case, each one of the remaining medoids will be assigned to each of the first threads.
 if ( remaining_points != 0 && current_thread_num < remaining_points)
  num_points_per_th++;

 start = current_thread_num*num_points_per_th;
 // The latests threads start higher, since remaining_medoids have been assigned to the initial threads.
 if (current_thread_num >= remaining_points)
  start += (FPp->num_obs % numthreads);

 // The end of the points assigned to this thread. The point with this number will NOT be managed by this thread;
 // it is the first point of the next one. Notice that the loops will run from point=start to point<end, not to point<=end.
 indextype end = start + num_points_per_th;
 if (end>FPp->num_obs)
  end=FPp->num_obs;

 pthread_exit(nullptr);
}

template void *FastPAM<float>::ExploreBranchesInternalThread(void *arg);
template void *FastPAM<double>::ExploreBranchesInternalThread(void *arg);
*/
template <typename disttype>
void FastPAM<disttype>::ExploreBranchesParallel(disttype *DeltaTDminusm,disttype *DeltaTD,vector<exchange> &xcg,unsigned int nt)
{
  // Now, local variables used in the paper's algorithm. Ths star (*) is translated as st so m* will be named mst
  disttype DeltaTDplusxc,d0j,DeltaTDst;
  indextype mst,xst;
  // i is an index in the vector of medoids.
  indextype i,imst;

  for (indextype m=0; m<nmed; m++)                            // L3
  {
    DeltaTDminusm[m]=disttype(0);
    for (indextype q=0; q<num_obs; q++)
     if (nearest[q]==m)
      DeltaTDminusm[m] += (dsecond[q]-dnearest[q]);   // Since dsecond[q] is always > dnearest[q], DeltaTDminus[m] will always be possitive for all m
  }

  // DIVIDE by threads from here...
  DeltaTDst = disttype(0);                                     // L4
  mst = num_obs+1;                       // This is our 'null'
  xst = num_obs+1;                       // Same here...
  i = nmed+1;                                // Just as a check to be tested outside the next loop, to see that i has been changed.
  imst = nmed+1;

  size_t B=xcg.size(),k;
  vector<disttype> vmin(B,MAXD);
  bool done_exchange=false;
  for (indextype xc=0; xc<num_obs; xc++)                        // L5
  {
    done_exchange=false;
    if (!ismedoid[xc])                                          // L5
    {

       for (indextype m=0; m<nmed; m++)                         // L6   DeltaTD is initialized to a possitive value for all m, since DeltaTDminusm[m] is possitive
           DeltaTD[m] = DeltaTDminusm[m];
       DeltaTDplusxc = disttype(0);                             // L7

       for (indextype x0=0; x0<num_obs; x0++)                   // L8
       {
         d0j = D->Get(x0,xc);                                   // L9
         if (d0j < dnearest[x0])                                // L10
         {
            DeltaTDplusxc += (d0j - dnearest[x0]);              // L11  Since d0j here is smaller then dnearest[x0], DeltaTDplusxc is reduced (becomes more negative)
            DeltaTD[nearest[x0]] += (dnearest[x0]-dsecond[x0]); // L12  and DeltaTD is reduced, since dnearest[x0] < dsecond[x0]
         }
         else
          if (d0j<dsecond[x0])                                // L13
            DeltaTD[nearest[x0]] += (d0j-dsecond[x0]);        // L14   Here DeltaTD is reduced, too, since in this part of the conditional d0j < dsecond[x0]
       }  // for (indextype x0=...

       disttype ddummy=MAXD;                                   // L15
       i=nmed+1;
       for (indextype m=0; m<nmed; m++)
        if (DeltaTD[m]<ddummy)
        {
            ddummy = DeltaTD[m];
            i = m;
        }

       DeltaTD[i] += DeltaTDplusxc;                       // L16

       if (DeltaTD[i]<DeltaTDst)                      // L17   DeltaTDst is the best improvement up to now. If this loop turn has made it improve...
       {
           DeltaTDst = DeltaTD[i];                    // ... update the best improvement and take note of the swap that has provoked it:
           mst = medoids[i];                          // The number of the point which is the medoid that should be left out
           xst = xc;                                  // The number of the point which is to become the new medoid
           imst = i;                                  // The index in the array of medoids which has the medoid to be swapped
           done_exchange=true;
       }
    }           // if (!ismedoid)...

    if (done_exchange)
    {
     k=0;
     while ( (k<B) && (vmin[k]<DeltaTDst) )
      k++;
     if (k<B)
     {
      for (size_t l=B-1;l>k;l--)
      {
       vmin[l]=vmin[l-1];
       xcg[l]=xcg[l-1];
      }
      vmin[k]=DeltaTDst;
      xcg[k].DeltaTDst=DeltaTDst;
      xcg[k].mst=mst;
      xcg[k].xst=xst;
      xcg[k].imst=imst;
     }
    }
  }   // for (indextype xc=...
 // ... TO HERE.
}

template void FastPAM<float>::ExploreBranchesParallel(float *DeltaTDminusm,float *DeltaTD,vector<exchange> &xcg,unsigned int nt);
template void FastPAM<double>::ExploreBranchesParallel(double *DeltaTDminusm,double *DeltaTD,vector<exchange> &xcg,unsigned int nt);

/*
 * Version 1: force increasing of intermedoid distace
template <typename disttype>
void FastPAM<disttype>::ChooseExchange(std::vector<exchange> &xcg,exchange &best_xcg)
{
 vector<disttype> val(xcg.size());
 vector<indextype> newmedoids;

 disttype vinit=disttype(0);
 for (indextype m=0;m<nmed;m++)
  for (indextype m1=0;m1<nmed;m1++)
   vinit += D->Get(medoids[m],medoids[m1]);
 for (size_t ex=0; ex<xcg.size(); ex++)
 {
  newmedoids=medoids;
  newmedoids[xcg[ex].imst]=xcg[ex].xst;
  val[ex]=disttype(0);
  for (indextype m=0;m<nmed;m++)
   for (indextype m1=0;m1<nmed;m1++)
    val[ex] += D->Get(newmedoids[m],newmedoids[m1]);
 }

 size_t best_exc=0;
 while (best_exc<xcg.size() && val[best_exc]<vinit)
  best_exc++;
 best_xcg = xcg[(best_exc==xcg.size()) ? 0 : best_exc];

 if (DEB & DEBPP)
 {
  cout << "\n  Possible exchanges: (vinit=" << vinit << ")\n";
  for (size_t x=0;x<xcg.size();x++)
   cout << "    Point " << xcg[x].mst << " <--> " << xcg[x].xst << " (place " << xcg[x].imst << ", DeltaTDst= " << xcg[x].DeltaTDst << ", val=" << val[x] << ")\n";
  cout << "  Choice: " << best_xcg.mst << " <--> " << best_xcg.xst << "\n";
 }

}
*/

// Version 2: force decreasing of silhouette
template <typename disttype>
void FastPAM<disttype>::ChooseExchange(std::vector<exchange> &xcg,exchange &best_xcg,unsigned int nt)
{
 siltype vinit=CalculateMeanSilhouette(nearest,nmed,D,nt);

 vector<disttype> val(xcg.size());
 vector<indextype> newmedoids;
 vector<indextype> newnearest(num_obs);
 indextype minmed;
 disttype mindist;
 for (size_t ex=0; ex<xcg.size(); ex++)
 {
  if (xcg[ex].DeltaTDst<0)
  {
   newmedoids=medoids;
   newmedoids[xcg[ex].imst]=xcg[ex].xst;
   for (indextype q=0;q<num_obs;q++)
   {
    mindist=MAXD;
    minmed=nmed+1;
    for (indextype m=0;m<nmed;m++)
     if (D->Get(q,newmedoids[m])<mindist)
     {
      mindist=D->Get(q,newmedoids[m]);
      minmed=m;
     }
    if (minmed>nmed)
     ParallelpamStop("Error: incorrect new medoid.\n");
    newnearest[q]=minmed;
   }
   val[ex]=CalculateMeanSilhouette(newnearest,nmed,D,nt);
  }
  else
   val[ex]=-1; // The minimum silhouette value, will never be chosen
 }

 size_t best_exc=0;
 siltype maxsil=-1;
 for (size_t i=0;i<val.size();i++)
  if (val[i]>maxsil)
  {
   maxsil=val[i];
   best_exc=i;
  }
 best_xcg = xcg[best_exc];

 //while (best_exc<xcg.size() && val[best_exc]<vinit)
 // best_exc++;
 // If no new value of the silhouette is bigger than the current silhouett, best_exc has arriave to the and of the array.
 // In such a case we choose exchange number 0, since it is the one with biggest decrease (most negative value) of TD.
 // Otherwise, we choose the biggest of those decreases which simultaneosuly increases silhouette
 //best_xcg = xcg[(best_exc==xcg.size()) ? 0 : best_exc];

 if (DEB & DEBPP)
 {
  cout << "\n  Possible exchanges: (vinit=" << vinit << ")\n";
  for (size_t x=0;x<xcg.size();x++)
   cout << "    Point " << xcg[x].mst << " <--> " << xcg[x].xst << " (place " << xcg[x].imst << ", DeltaTDst= " << xcg[x].DeltaTDst << ", val=" << val[x] << ")\n";
  cout << "  Choice: " << best_xcg.mst << " <--> " << best_xcg.xst << "\n";
 }
}

template void FastPAM<float>::ChooseExchange(std::vector<exchange> &xcg,exchange &best_xcg,unsigned int nt);
template void FastPAM<double>::ChooseExchange(std::vector<exchange> &xcg,exchange &best_xcg,unsigned int nt);

/**************************** RunImprovedFastPAMMultiBranch (optimization phase, serial version) *****************/
template <typename disttype>
void FastPAM<disttype>::RunImprovedFastPAMMultiBranch(unsigned int B,unsigned int nt)
{
 if (DEB & DEBPP)
 {
  std::cout << "Starting improved FastPAM Multibranch method in ";
  if (nt>1)
   std::cout << "parallel implementation with " << nt << " threads...\n";
  else
   std::cout << "serial implementation...\n";
  //Rcpp::Rcout << "WARNING: all successive messages use R-numbering (from 1) for points and medoids. Substract 1 to get the internal C-numbers.\n";
  std::cout.flush();
 }

 // dsecond is to be filled in advance, mostly as cache.
 FillSecond();

 // The threshold that will stop the algorithm if TD changes less than this value at any iteration.
 disttype tol_limit=currentTD*tlimit;

 // Now, local variables used in the paper's algorithm. Ths star (*) is translated as st so m* will be named mst
 disttype *DeltaTDminusm = new disttype [nmed];
 //disttype DeltaTDplusxc,DeltaTDst,d0j;
 disttype *DeltaTD = new disttype [nmed];

 /*
 // I take these as number of the point and number of the medoid
 indextype xst,mst;
 // Nevertheless, these are indexes in the vector of medoids. imst is not explictly named in the original work.
 indextype i,imst;
 */

 unsigned int iteration=0;
 bool out=false;                         // Used to leave in special case of no TD improvement, i.e., no better solucion exist.
 vector<exchange> xcg(B);

 exchange chosen_exchange;
 do                                                          // L2
 {
  if (DEB & DEBPP)
  {
   std::cout << "Iteration " << iteration << ". ";
   std::cout.flush();
  }

  for (size_t i=0;i<B;i++)
  {
   xcg[i].DeltaTDst=disttype(0);
   xcg[i].mst=0;
   xcg[i].xst=0;
   xcg[i].imst=0;
  }
  if (nt>1)
   ExploreBranchesParallel(DeltaTDminusm,DeltaTD,xcg,nt);
  else
   ExploreBranches(DeltaTDminusm,DeltaTD,xcg);

  ChooseExchange(xcg,chosen_exchange,nt);

  if (xcg[0].DeltaTDst>=disttype(0))                        // L18
  {
      if (DEB & DEBPP)
        std::cout << "   Exiting, since DeltaTDst is " << std::fixed << chosen_exchange.DeltaTDst/float(num_obs) << ". Final value of TD is " << std::fixed << currentTD/float(num_obs) << "\n";
      out=true;
      break;
  }

  if ((DEB & DEBPP) && (!out) && (chosen_exchange.imst<nmed))
   std::cout << "Medoid at place " << chosen_exchange.imst << " (point " << medoids[chosen_exchange.imst] << ") swapped with point " << chosen_exchange.xst << "; ";

  // Could imst be left unchanged by the loop? Not except by error. Anyway, let's check it
  if ((chosen_exchange.imst<nmed) && (!out))
  {
   SwapRolesAndUpdate(chosen_exchange.mst,chosen_exchange.xst,chosen_exchange.imst);                                // L19-20

   currentTD += chosen_exchange.DeltaTDst;                                          // L21

   if (DEB & DEBPP)
    std::cout << "TD-change=" << std::fixed << chosen_exchange.DeltaTDst/float(num_obs) << "; TD=" << std::fixed << currentTD/float(num_obs) << ". " << current_npch << " reassigned points.\n";
  }
  else
  {
     if (DEB & DEBPP)
     {
      std::cout << "   No exchange of medoid/point found which can improve result. Exact result found?\n";
      std::cout << "   Last TD change has been " << std::fixed << chosen_exchange.DeltaTDst/float(num_obs) << "\n";
      if (chosen_exchange.imst>nmed)
          std::cout << "Best medoid has not been updated.\n";
      else
          std::cout << "Nevertheless, best medoid has been updated to " << medoids[chosen_exchange.imst] << ". ????\n";
     }
     out=true;
  }

  iteration++;

  // This is the only point (apart from the messages in the screen) in which TD is converted from raw sum to sum per point.
  TDkeep.push_back(currentTD/float(num_obs));
  NpointsChangekeep.push_back(current_npch);
 }
 while ((fabs(chosen_exchange.DeltaTDst)>tol_limit) && (iteration<maxiter) && (!out));   // fabs because DeltaTDst is negative...
 num_iterations_in_opt=(iteration>0) ? iteration-1 : 0;

 delete[] DeltaTDminusm;
 delete[] DeltaTD;
}

template void FastPAM<float>::RunImprovedFastPAMMultiBranch(unsigned int B,unsigned int nt);
template void FastPAM<double>::RunImprovedFastPAMMultiBranch(unsigned int B,unsigned int nt);

// FINALLY, TWO AUXILIARY FUNCTIONS USED BY ALL VERSIONS (serial and parallel) OF FASTPAM1 AND FASTPAM2B

/***************** FillSecond (first auxiliary function) **************************/
template <typename disttype>
void FastPAM<disttype>::FillSecond()
{ 
 dsecond.clear();
 for (indextype q=0; q<num_obs; q++)
  dsecond.push_back(MAXD);
 
 // The dnearest (distance to closest medoid) is already in the class data, since it is used in BUILD/LAB 
 // and also later in the algorithm. But the distance to second-closest medoid (dsecond) is to be used just here.
 // We only need distances. The concrete medoid which is second is not used later, so we don't search for neither store it 
 
 disttype minseconddist,dd;
 for (indextype q=0; q<num_obs; q++)
 {
     minseconddist=MAXD;
     for (indextype m=0; m<nmed; m++)
      if (m!=nearest[q])   // By doing this we are excluding the closest medoid of the search.
                               // Therefore, the minimum will be the second-closest
      {
          dd=D->Get(q,medoids[m]);
          if (dd < minseconddist)
              minseconddist = dd;
      }
    dsecond[q]=minseconddist;
 }
}

template void FastPAM<float>::FillSecond();
template void FastPAM<double>::FillSecond();

/******************** SwapRolesAndUpdate (second auxiliary function) **************************/
template <typename disttype>
void FastPAM<disttype>::SwapRolesAndUpdate(indextype mst,indextype xst,indextype imst)
{
   if (mst!=medoids[imst])
   {
    std::ostringstream errst;
    errst << "Error in SwapRolesAndUpdate: medoid " << mst << " is not at place " << imst << " of medoids array.\n";
    errst << "The medoid at such place is point " << medoids[imst] << ".\n";
    errst << "Unexpected error.\n";
    ParallelpamStop(errst.str());
   }
   
   // L16 (swap roles...) comprises several different tasks: 
   
   // Updates the array of marks:
   ismedoid[mst]=false;             
   ismedoid[xst]=true;
   
   medoids[imst]=xst;

   // Now, update nearest and dnearest 
   current_npch = 0;
   disttype dd,mind;
   indextype closestmed=nmed+1;

   for (indextype q=0; q<num_obs; q++)
   {
    mind=MAXD;
    for (indextype m=0; m<nmed; m++)
    {
     if ((dd=D->Get(q,medoids[m]))<mind)
     {
      mind=dd;
      closestmed=m;
     } 
    }
  
    if (nearest[q]!=closestmed)
     current_npch++;
    nearest[q]=closestmed;
    dnearest[q]=mind;
   }
   
   // and fill the array of distance to the second-closest point
   FillSecond();
}

template void FastPAM<float>::SwapRolesAndUpdate(indextype mst,indextype xst,indextype imst);
template void FastPAM<double>::SwapRolesAndUpdate(indextype mst,indextype xst,indextype imst);

/******************** Functions to return JMatrix from the internal representation *****/
template<typename disttype>
FullMatrix<indextype> & FastPAM<disttype>::GetMedoids()
{
 FullMatrix<indextype> *M = new FullMatrix<indextype>(medoids.size(),1);
 for (indextype m=0;m<medoids.size();m++)
  M->Set(m,0,medoids[m]);
 return(*M);
}

template FullMatrix<indextype> &FastPAM<float>::GetMedoids();
template FullMatrix<indextype> &FastPAM<double>::GetMedoids();

/**********************************/
template<typename disttype>
FullMatrix<indextype> & FastPAM<disttype>::GetMedoids(vector<string> rownames)
{
 FullMatrix<indextype> &M = GetMedoids();
 if (rownames.size()>0)
 {
  vector<string> mednames;
  for (indextype m=0;m<medoids.size();m++)
   if (medoids[m]<rownames.size())
    mednames.push_back(rownames[medoids[m]]);
   else
    ParallelpamStop("In function GetMedoids: number of medoid would be outside the vector or point names. Have you passed a correct vector of names?");
  M.SetRowNames(mednames);
 }
 return(M);
}

template FullMatrix<indextype> &FastPAM<float>::GetMedoids(vector<string> rownames);
template FullMatrix<indextype> &FastPAM<double>::GetMedoids(vector<string> rownames);

/***********************************/
template<typename disttype>
FullMatrix<indextype> & FastPAM<disttype>::GetAssign()
{
 FullMatrix<indextype> *M = new FullMatrix<indextype>(nearest.size(),1);
 for (indextype m=0;m<nearest.size();m++)
  M->Set(m,0,nearest[m]);
 return(*M);
}

template FullMatrix<indextype> &FastPAM<float>::GetAssign();
template FullMatrix<indextype> &FastPAM<double>::GetAssign();

/**********************************/
template<typename disttype>
FullMatrix<indextype> & FastPAM<disttype>::GetAssign(vector<string> rownames)
{
 FullMatrix<indextype> &M = GetAssign();

 if (rownames.size()>0)
 {
  if (rownames.size()!=nearest.size())
   ParallelpamStop("In function GetAssign: length of vector of names is not equal to the number of points. Have you passed a correct vector of names?");
  M.SetRowNames(rownames);
 }
 return(M);
}

template FullMatrix<indextype> &FastPAM<float>::GetAssign(vector<string> rownames);
template FullMatrix<indextype> &FastPAM<double>::GetAssign(vector<string> rownames);
