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

#include "../headers/threadhelper.h"
#include "../headers/debugpar_ppam.h"

// Auxiliary function to decide the number of threads
unsigned int ChooseNumThreads(int nthreads)
{
 // This is the mark to indicate strictly serial computation (i.e.: no threads)
 if (nthreads==NO_THREADS)
  return 1;
  
 unsigned int nthc=std::thread::hardware_concurrency();
 unsigned int nt;
    
 if (nthreads==AS_MANY_AS_POSSIBLE)
  nt=nthc;
 else
 {
  if ((unsigned int)nthreads>nthc)
  {
   std::ostringstream warst;
   warst << "You have requested " << nthreads << " threads but this machine has only " << nthc << " cores. This is allowed, but discouraged.\n";
   ParallelpamWarning(warst.str());
  }
  nt=nthreads;
 }
 
 return nt;
}

void CreateAndRunThreadsWithDifferentArgs(unsigned int numthreads, void *(*ThreadFunction)(void *),void *io,size_t iosize)
{
 pthread_t *thr = new pthread_t [numthreads];
 arg_to_thread *a = new arg_to_thread [numthreads];

 for (unsigned t=0; t<numthreads; t++)
 {
  a[t].numthreads=numthreads;
  a[t].num_this_thread=t;
  a[t].ioargs=(void *)((char *)io+t*iosize);
  pthread_create(&thr[t],NULL,ThreadFunction,(void *)(&(a[t])));
 }
 for (unsigned t=0; t<numthreads; t++)
  pthread_join(thr[t],nullptr);
  
 delete[] thr;
 delete[] a;
}

void CreateAndRunThreadsWithSameArgs(unsigned int numthreads, void *(*ThreadFunction)(void *),void *io)
{
 pthread_t *thr = new pthread_t [numthreads];
 arg_to_thread *a = new arg_to_thread [numthreads];

 for (unsigned t=0; t<numthreads; t++)
 {
  a[t].numthreads=numthreads;
  a[t].num_this_thread=t;
  a[t].ioargs=(void *)((char *)io);
  pthread_create(&thr[t],NULL,ThreadFunction,(void *)(&(a[t])));
 }
 for (unsigned t=0; t<numthreads; t++)
  pthread_join(thr[t],nullptr);
  
 delete[] thr;
 delete[] a;
}

unsigned int GetNumThreads(void *arg)
{
    return ((arg_to_thread *)arg)->numthreads;
}

unsigned int GetThisThreadNumber(void *arg)
{
    return ((arg_to_thread *)arg)->num_this_thread;
}
