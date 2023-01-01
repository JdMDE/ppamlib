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

#include "../headers/dissimmat.h"
#include "../headers/debugpar_ppam.h"
#include "../headers/threadhelper.h"
#include "../headers/diftimehelper.h"

extern unsigned char DEB;

// This function will fill part of the distance matrix D, concretely, lines between initial_row and (but not including) final_row
template <typename counttype,typename disttype>
void FillMetricMatrixFromFull(indextype initial_row,indextype final_row,FullMatrix<counttype> *M,SymmetricMatrix<disttype> *D,bool L1dist)
{
 disttype d,dif;
 indextype ncols=M->GetNCols();
 indextype nrows=D->GetNRows();
 
 // This should not be done from inside a thread. But, if we have failed anyway....
 if ( (initial_row >= nrows) || (final_row > nrows) )
 {
     std::ostringstream errst;
     errst << "Error in FillMetricMatrixFromFull: start of area at " << initial_row << " or end of area at " << final_row << " outside matrix limits.\n";
     ParallelpamStop(errst.str());
     return;
 }
 
 counttype *va = new counttype [ncols];
 counttype *vb = new counttype [ncols];
 unsigned char *mark = new unsigned char [ncols];
 unsigned char *rowmark = new unsigned char [ncols];

 for (indextype rowA=initial_row; rowA<final_row; rowA++)
 {
  // The values of the current row (let's call it rowA) and the places they are, are stored in va and rowmark respectively
  memset((void *)va,0x0,ncols*sizeof(counttype));
  memset((void *)rowmark,EMPTY,ncols);
  
  M->GetFullRow(rowA,rowmark,IN_FIRST,va);
  
  // The next loop calculates the distance between the current row (rowA) and all others with numbers below its own number, let's call rowB to each.
  // (remember that to fill distance matrix we only need to fill the lower-diagonal part...)
  
  for (unsigned rowB=0; rowB<rowA; rowB++)
  {
   // The mark array is initialized to the marks of rowA
   memcpy(mark,rowmark,ncols);
   
   // and vb will contains the values of rowB. First, we clean it...
   memset((void *)vb,0x0,ncols*sizeof(counttype));
   
   // and now, we fill its values and update mark adding 2 to the rigth places. In that way, mark will return with 0, 1, 2 or 3
   // for the cases of "No value at that place in any row", "Value only in first row", "Value only in second row" or "Value in both", respectively.
   M->GetFullRow(rowB,mark,IN_SECOND,vb);
   
   d=0.0;
   for (unsigned col=0; col<ncols; col++)
   {
    // That component of the vector is updated only if there is something to consider...
    if (mark[col]!=EMPTY)
    {
        dif = ((mark[col]==IN_FIRST) ? disttype(va[col]) : ((mark[col]==IN_SECOND) ? -disttype(vb[col]) : disttype(va[col])-disttype(vb[col]) ));
        d += (L1dist ? fabs(dif) : dif*dif);
    }
   } 
   
   D->Set(rowA,rowB,(L1dist ? d : sqrtf(d)));
  }
  
  // This is just to set the main diagonal.
  D->Set(rowA,rowA,disttype(0));
 }
 
 delete[] va;
 delete[] vb;
 delete[] mark;
 delete[] rowmark;
}

template void FillMetricMatrixFromFull(indextype initial_row,indextype final_row,FullMatrix<float> *M,SymmetricMatrix<float> *D,bool L1dist);
template void FillMetricMatrixFromFull(indextype initial_row,indextype final_row,FullMatrix<double> *M,SymmetricMatrix<float> *D,bool L1dist);
template void FillMetricMatrixFromFull(indextype initial_row,indextype final_row,FullMatrix<float> *M,SymmetricMatrix<double> *D,bool L1dist);
template void FillMetricMatrixFromFull(indextype initial_row,indextype final_row,FullMatrix<double> *M,SymmetricMatrix<double> *D,bool L1dist);

template <typename counttype>
void CalculateMeansFromFull(FullMatrix<counttype> &M,std::vector<counttype> &mu)
{
 indextype ncells=M.GetNRows();
 indextype ngenes=M.GetNCols();
 
 counttype s;   
 for (indextype gene=0; gene<ngenes; gene++)
 {
  s=counttype(0); 
  for (indextype cell=0; cell<ncells; cell++)
   s += M.Get(cell,gene);
  mu.push_back(s/counttype(ncells)); 
 }
}

template void CalculateMeansFromFull(FullMatrix<float> &M,std::vector<float> &mu);
template void CalculateMeansFromFull(FullMatrix<double> &M,std::vector<double> &mu);

// This function will fill part of the distance matrix D, concretely, lines between initial_row and (but not including) final_row
template <typename counttype,typename disttype>
void FillPearsonMatrixFromFull(indextype initial_row,indextype final_row,FullMatrix<counttype> *M,std::vector<counttype> *mu,SymmetricMatrix<disttype> *D)
{
 disttype da,db,sxx,syy,sxy,den,pearson;
 disttype dtol=1e-06;
 indextype ncols=M->GetNCols();
 indextype nrows=D->GetNRows();
 
 // This should not be done from inside a thread. But, if we have failed anyway....
 if ( (initial_row >= nrows) || (final_row > nrows) )
 {
     std::ostringstream errst;
     errst << "Error in FillPearsonMatrixFromFull: start of area at " << initial_row << " or end of area at " << final_row << " outside matrix limits.\n";
     ParallelpamStop(errst.str());
     return;
 }
 
 // Memory is booked for all the row. Only part of it (from initial_col to initial_col+extension) will be
 // used, but we have memory enough and it is easier this way. The same applies to the arrays of marks.
 counttype *va = new counttype [ncols];
 counttype *vb = new counttype [ncols];
 
 for (unsigned rowA=initial_row; rowA<final_row; rowA++)
 {
  // The values of the current row (let's call it rowA) are stored in va 
  memset((void *)va,0x0,ncols*sizeof(counttype));
  
  M->GetRow(rowA,va);
  
  // The next loop calculates the distance between the current row (rowA) and all others with numbers below its own number, let's call rowB to each.
  // (remember that to fill distance matrix we only need to fill the lower-diagonal part...)
  
  for (unsigned rowB=0; rowB<rowA; rowB++)
  {
   // bb will contains the values of rowB. First, we clean it...
   memset((void *)vb,0x0,ncols*sizeof(counttype));
   
   // and now, we fill its values
   M->GetRow(rowB,vb);
    
   sxx=sxy=syy=0.0;
   for (unsigned col=0; col<ncols; col++)
   {
    da=va[col]-disttype((*mu)[col]);
    db=vb[col]-disttype((*mu)[col]);
    sxx += (da*da);
    syy += (db*db);
    sxy += (da*db);
   }
   
   // This is to prevent from numerical problems
   den=disttype(sqrt(double(sxx)))*disttype(sqrt(double(syy)));
   if (den==0.0)   // This is the pathological case in which both vectors are the null vector. Then, they are of course completely "similar" (indeed, identical...)
    D->Set(rowA,rowB,disttype(0.0));
   else
   {
    pearson=0.5-(sxy/den/2.0);
    D->Set(rowA,rowB,(fabs(pearson)<dtol) ? disttype(0.0) : pearson);
   }
  }
  // This is just to set the main diagonal.
  D->Set(rowA,rowA,disttype(0));
 }
 delete[] va;
 delete[] vb;
}

template void FillPearsonMatrixFromFull(indextype initial_row,indextype final_row,FullMatrix<float> *M,std::vector<float> *mu,SymmetricMatrix<float> *D);
template void FillPearsonMatrixFromFull(indextype initial_row,indextype final_row,FullMatrix<double> *M,std::vector<double> *mu,SymmetricMatrix<float> *D);
template void FillPearsonMatrixFromFull(indextype initial_row,indextype final_row,FullMatrix<float> *M,std::vector<float> *mu,SymmetricMatrix<double> *D);
template void FillPearsonMatrixFromFull(indextype initial_row,indextype final_row,FullMatrix<double> *M,std::vector<double> *mu,SymmetricMatrix<double> *D);

template <typename counttype,typename disttype>
void *BasicThreadFull(void *arg)
{
 // Commas into () are called protected commas; this is the way to pass to a macro commas which must be part of a macro argument
 indextype initial_row1 = GetFieldDT(arg,args_to_full_thread,counttype,disttype,initial_row1);
 indextype final_row1 = GetFieldDT(arg,args_to_full_thread,counttype,disttype,final_row1);
 indextype initial_row2 = GetFieldDT(arg,args_to_full_thread,counttype,disttype,initial_row2);
 indextype final_row2 = GetFieldDT(arg,args_to_full_thread,counttype,disttype,final_row2);
 FullMatrix<counttype> *M = GetFieldDT(arg,args_to_full_thread,counttype,disttype,M);
 SymmetricMatrix<disttype> *D = GetFieldDT(arg,args_to_full_thread,counttype,disttype,D);
 std::vector<counttype> *mu = GetFieldDT(arg,args_to_full_thread,counttype,disttype,mu);
 unsigned char dtype = GetFieldDT(arg,args_to_full_thread,counttype,disttype,dtype);

 switch (dtype)
 {
  case DL1: {
  		FillMetricMatrixFromFull(initial_row1,final_row1,M,D,true);
            	FillMetricMatrixFromFull(initial_row2,final_row2,M,D,true);
            }
            break;
  case DL2: {
  		FillMetricMatrixFromFull(initial_row1,final_row1,M,D,false);
            	FillMetricMatrixFromFull(initial_row2,final_row2,M,D,false);
            }
            break;
  case DPe: {
  		FillPearsonMatrixFromFull(initial_row1,final_row1,M,mu,D);
            	FillPearsonMatrixFromFull(initial_row2,final_row2,M,mu,D);
            }
            break;
  default: break;
 }
 pthread_exit(NULL);
}

template void *BasicThreadFull<float,float>(void *arg);
template void *BasicThreadFull<float,double>(void *arg);
template void *BasicThreadFull<double,float>(void *arg);
template void *BasicThreadFull<double,double>(void *arg);

template <typename counttype,typename disttype>
SymmetricMatrix<disttype> &CalcDistFromFull(FullMatrix<counttype> &M,unsigned char dtype, unsigned int nthr)
{
 indextype nrows=M.GetNRows();
 if (DEB & DEBPP)
  std::cout << "Creating dissimilarity matrix of size (" << nrows << "x" << nrows << ")\n";
 SymmetricMatrix<disttype> *D = new SymmetricMatrix<disttype>(nrows,true);

 std::vector<counttype> mu;
 if (dtype==DPe)
 {
  if (DEB & DEBPP)
   std::cout << "Calculating vector of means used by the Pearson dissimilarity...\n";
  CalculateMeansFromFull(M,mu);
  if (mu.size()!=M.GetNCols())
   ParallelpamStop("Error from CalcAndWriteAuxSparse: length of vector of means is not the number of columns of the data matrix.\n");
 }
 
 if ((nrows<1000) && (nthr!=1))
 {
  nthr=1;
  if (DEB & DEBPP)
   std::cout << "We will calculate with a single thread, since you have only " << nrows << " vectors and the overhead of using threads would be excessive.\n";
 }
 
 DifftimeHelper Dt;
 
 if (nthr==1)
 {
  Dt.StartClock("End of dissimilarity matrix calculation (serial version)."); 
  switch (dtype)
  {
   case DL1: FillMetricMatrixFromFull(0,D->GetNRows(),&M,&(*D),true); break;
   case DL2: FillMetricMatrixFromFull(0,D->GetNRows(),&M,&(*D),false); break;
   case DPe: FillPearsonMatrixFromFull(0,D->GetNRows(),&M,&mu,&(*D)); break;
   default: break;
  }
  Dt.EndClock(DEB & DEBPP);
 }
 else
 {
  Dt.StartClock("End of dissimilarity matrix calculation (parallel version)."); 

  indextype nrows=D->GetNRows();
 
  args_to_full_thread<counttype,disttype> *fullargs = new args_to_full_thread<counttype,disttype> [nthr];
  // This strange distribution of rows (each thread has two intervals which are symmetric with respect to the middle of the matrix rows)
  // is to balance well the number of distances calculated by each thread. For each 'sort' row (those before the middle) the same thread does also a 'long' row (those after the middle)
  for (unsigned int t=0; t<nthr; t++)
  {
   fullargs[t].initial_row1 = indextype(float(t)*float(nrows)/(2*float(nthr)));
   if (t==nthr-1)
    fullargs[t].final_row1 = nrows/2;
   else
    fullargs[t].final_row1 = indextype(float(t+1)*float(nrows)/(2*float(nthr)));
      
   fullargs[t].initial_row2 = nrows-fullargs[t].final_row1 - 1;
   if ((t==nthr-1) && (!(nrows%2)))
    fullargs[t].initial_row2++;
    
   fullargs[t].final_row2 = nrows-fullargs[t].initial_row1 - 1;
   if (t==0)
    fullargs[t].final_row2++;
    
   fullargs[t].M = &M;
   fullargs[t].D = D;
   fullargs[t].mu = &mu;
   fullargs[t].dtype = dtype;
  }
  if (DEB & DEBPP)
     std::cout << "Launching up to " << nthr << " simultaneous threads.\n";

  CreateAndRunThreadsWithDifferentArgs(nthr,BasicThreadFull<counttype,disttype>,(void *)fullargs,sizeof(args_to_full_thread<counttype,disttype>));
  
  delete[] fullargs;
 
  Dt.EndClock(DEB & DEBPP);
 }
     
 D->SetRowNames(M.GetRowNames());

 return(*D);
}

template SymmetricMatrix<float>  &CalcDistFromFull<float,float>(  FullMatrix<float>  &M, unsigned char dtype, unsigned int nthr);
template SymmetricMatrix<double> &CalcDistFromFull<float,double>( FullMatrix<float>  &M, unsigned char dtype, unsigned int nthr);
template SymmetricMatrix<float>  &CalcDistFromFull<double,float>( FullMatrix<double> &M, unsigned char dtype, unsigned int nthr);
template SymmetricMatrix<double> &CalcDistFromFull<double,double>(FullMatrix<double> &M, unsigned char dtype, unsigned int nthr);

