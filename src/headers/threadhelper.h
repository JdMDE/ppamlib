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

#ifndef _THREADHELPER_H
#define _THREADHELPER_H

#include <sstream>
#include <cstdlib>
#include <thread>

/// @file threadhelper.h

const int NO_THREADS=-1;
const int AS_MANY_AS_POSSIBLE=0;

/**
 * Auxiliary function to decide the number of threads
 *
 * @param nthreads Number of thread the user wish to launch. It can be:\n
 * -1 to indicate not to use threads (strictly sequential computation)\n
 * 0 to let the program choose according to the available number of cores/threads in the machine\n
 * Other possitive number to ask for such number of threads\n
 *
 * @return The number of threads that will be really used. This is:\n
 * 1 for input value -1\n
 * the number of cores/threads in the machine for input value 0 (if hypertheading exists, it is taken into account)\n
 * the chosen number for any other input value (but a warning message is raised if the chosen number is bigger than the number of cores/threads)
 */
unsigned int ChooseNumThreads(int nthreads);

#ifndef DOXYGEN_SHOULD_SKIP_THIS
// This is a structure that will be passed to every thread, it does not matter which function uses it,
// since it is meant as a kind of "universal" container for thread-related information
// This is really internal, the user do no initialize it. It is filled by one of the helper functions explained immediately after.
// Nevertheless, the user can access to it inside each thread, and he/she should do it exclusively through functions GetNumThreads
// and GetThisThreadNumber and the macro GetField
typedef struct
{
 unsigned int numthreads;         // How many thread we are creating for this task
 unsigned int num_this_thread;    // Which one of them is this particular thread
 void *ioargs;                    // A pointer to an array of arbitrary structures used to pass inputs and return outputs from each thread.
} arg_to_thread;
  
// The helper function to create threads. It will take the number of threads to be created and the function which each thread must run,
// together with an array (pointer) of input/output structures. Obviously, we don't know the size of each structure (it depends on the
// actual application) so we must pass the size of such structure, too.
// There is one obvious limitation: all threads will run the same function, and this function (which will receive a pointer to arg_to_thread
// structure as its only argument) will know which thread it is, and consequently, which part of the total job is assigned to it.
void CreateAndRunThreadsWithDifferentArgs(unsigned int numthreads,void *(*ThreadFunction)(void *),void *io,size_t iosize);

// Another helper function to create threads. It will take the number of threads to be created and the function which each thread must run,
// together with an pointer to the input/output structure. 
// There is one obvious limitation: all threads will run the same function, and this function (which will receive a pointer to arg_to_thread
// structure as its only argument) will know which thread it is, and consequently, which part of the total job is assigned to it.
void CreateAndRunThreadsWithSameArgs(unsigned int numthreads,void *(*ThreadFunction)(void *),void *io);

unsigned int GetNumThreads(void *arg);
unsigned int GetThisThreadNumber(void *arg);

// Convenience macros to get the member of a structure inside a thread to which only a pointer to the structure has been passed (because this is
// the only thing allowed by the pthread library)
// Obtain field F from a struct of type S pointed by the pointer a. S can be a templated structure since mystruct<something> is passed as an argument, in the same
// way it is with a non-templated structure
#define GetField(a,S,F) ((S *)(((arg_to_thread *)a)->ioargs))->F

// The problem is with templated structured that depend on two (or more) template values like mystruct<something,othering> since the macro preprocessor
// considers the comma as an argument separator. Protecting the comma with (,) makes the preprocessor happy but confuses the compiler syntatic analizer,
// so struct name, template instance 1 and instance 2 must be passed separately.
#define GetFieldDT(a,S,T1,T2,F) ((S<T1,T2> *)(((arg_to_thread *)a)->ioargs))->F
#endif

#endif
