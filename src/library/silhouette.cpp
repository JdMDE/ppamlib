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

#include "../headers/silhouette.h"
#include "../headers/diftimehelper.h"
#include "../headers/threadhelper.h"
#include "../headers/debugpar_ppam.h"

#include <jmatrixlib/matmetadata.h>

extern unsigned char DEB;

// Thread used to calculate the silhouette in its parallel implementation
template <typename disttype>
void *SilhoutteThread(void *arg)
{
  unsigned int numthreads = GetNumThreads(arg);
  unsigned int current_thread_num = GetThisThreadNumber(arg);
  
  indextype num_obs = GetField(arg,SilhoutteThread_IO<disttype>,num_obs);
  indextype nmed = GetField(arg,SilhoutteThread_IO<disttype>,nmed);
  const std::vector<indextype> *nearest = GetField(arg,SilhoutteThread_IO<disttype>,nearest);
  std::vector<siltype> *current_sil = GetField(arg,SilhoutteThread_IO<disttype>,current_sil);
  std::vector<unsigned long> *hist = GetField(arg,SilhoutteThread_IO<disttype>,hist);
  std::vector<silinfo> *silres = GetField(arg,SilhoutteThread_IO<disttype>,silres);
  SymmetricMatrix<disttype> *D = GetField(arg,SilhoutteThread_IO<disttype>,D);
  
  indextype start = current_thread_num*(num_obs/numthreads);
  indextype end   = (current_thread_num == numthreads-1) ? num_obs : (current_thread_num+1)*(num_obs/numthreads);
  
  // To interpret all variables and operation here, please look at the comments in the serial version below
  siltype *a = new siltype [end-start];
  siltype *b = new siltype [end-start];
  siltype *bav = new siltype [nmed];
  siltype dmin;
  indextype which_neimin=nmed+1;
  for (indextype q=start; q<end; q++)
  {
     if ((*hist)[(*nearest)[q]]==1)
      (*current_sil)[q]=0.0;
     else
     {
      for (indextype m=0; m<nmed; m++)
       bav[m]=0.0;
      
      for (indextype q1=0; q1<num_obs; q1++)
       bav[(*nearest)[q1]] += D->Get(q,q1);
          
      for (indextype m=0; m<nmed; m++)
       if (m==(*nearest)[q])
        bav[m] /= double((*hist)[m]-1);
       else
        bav[m] /= double((*hist)[m]);
      
      a[q-start] = bav[(*nearest)[q]];
        
      dmin=std::numeric_limits<siltype>::max();
      for (indextype m=0; m<nmed; m++)
       if ( (m!=(*nearest)[q]) && (bav[m]<dmin) )
       {
           which_neimin=m;
           dmin = bav[m];
       } 
       
      b[q-start]=dmin;   
      
      (*current_sil)[q]=(b[q-start]-a[q-start])/std::max(a[q-start],b[q-start]);
     }
     
     (*silres)[q].neiclus=which_neimin;
     (*silres)[q].silvalue=(*current_sil)[q];
  }
 
  delete[] a;
  delete[] b;
  delete[] bav;
  
  pthread_exit(nullptr);  
}

template void *SilhoutteThread<float>(void *arg);
template void *SilhoutteThread<double>(void *arg);

// Auxiliary function with the real implementation of silhouette (serial version)
template <typename disttype>
void SilhouetteSerial(indextype num_obs,indextype nmed,
                      std::vector<indextype> &nearest,std::vector<siltype> &current_sil,std::vector<unsigned long> &hist,
                      std::vector<silinfo> &silres,SymmetricMatrix<disttype> &D)
{
 siltype *a = new siltype [num_obs];
 siltype *b = new siltype [num_obs];
 siltype *bav = new siltype [nmed];
 siltype dmin;
    
 indextype which_neimin=nmed+1;
 for (indextype q=0; q<num_obs; q++)
 {
   // Special case: cluster with one isolated point. Silhouette in this case is defined as 0
   if (hist[nearest[q]]==1)
     current_sil[q]=0.0;
   else
   {  
    // bav will contain the average distance between this points and the others in each cluster, including its own cluster
    for (indextype m=0; m<nmed; m++)
     bav[m]=0.0;
         
    for (indextype q1=0; q1<num_obs; q1++)
     bav[nearest[q1]] += D.Get(q,q1);
          
    // bav contains now the sum of distances to point q, by cluster. Let's divide to calculate the average.
    // The 'minus 1' is because the point itself is not counted. This is the definition of silhouette.
    // In this case, the denominator cannot be 0, since the special case of hist[nearest[q]]==1 was managed before
    // Also, it is not possible to have hist==0 for other clusters, since every cluster is requested to have at least one point
    for (indextype m=0; m<nmed; m++)
     if (m==nearest[q])
      bav[m] /= double(hist[m]-1);
     else
      bav[m] /= double(hist[m]);
           
    // Now, for each point q, a will contain the average distance to the other points in its own cluster:
    a[q] = bav[nearest[q]];
         
    // and b will contain the minimal average distance to points in _other_ clusters.
    // This is why we leave out the own cluster of the point
    // We keep also the number of the cluster at minimal distance, in which_neimin
    dmin=std::numeric_limits<siltype>::max();
    for (indextype m=0; m<nmed; m++)
     if ( (m!=nearest[q]) && (bav[m]<dmin) )
     {
      which_neimin=m;
      dmin = bav[m];
     }  
     
    b[q]=dmin;   
          
    current_sil[q]=(b[q]-a[q])/std::max(a[q],b[q]);
   }
        
   silres[q].neiclus=which_neimin;
   silres[q].silvalue=current_sil[q];       
 }
 
 delete[] a;
 delete[] b;
 delete[] bav;
  
}

template <typename disttype>
std::vector<siltype> CalculateSilhouette(std::vector<indextype> cl,SymmetricMatrix<disttype> &D,unsigned int nt)
{
 DifftimeHelper Dt;
 if (nt==1)
 {
   if (DEB & DEBPP)
   {
    std::cout << "   Calculating silhouette (serial implementation)...\n";
    std::cout.flush();
   }
   Dt.StartClock("Finished serial implementation of silhouette (including dissimilarity matrix load).");
 }
 else
 {
   if (DEB & DEBPP)
   {
    std::cout << "   Calculating silhouette (parallel version) with " << nt << " threads.\n";
    std::cout.flush();
   }
   Dt.StartClock("Finished parallel implementation of silhouette (including dissimilarity matrix load).");
 }

 indextype num_obs=D.GetNRows();
 
 // Nearest is the vector of nearest cluster number (in C++ numbering, from 0)
 std::vector<indextype> nearest;

 indextype mincl=num_obs;
 indextype maxcl=0;
 
 // This is just to check cl vector, which is in C numbering (from 0) and set nearest from it
 for (size_t t=0; t<size_t(cl.size()); t++)
 {
  if (cl[t]<0 || cl[t]>=num_obs)
   ParallelpamStop("The clasification array contains at least one invalid value (outside the range 1.. number_of_points).\n");
  if (cl[t] < mincl)
   mincl=cl[t];
  if (cl[t] > maxcl)
   maxcl=cl[t];
  nearest.push_back(cl[t]);
 }
 if (mincl!=0)
 {
  std::cerr << "Minimum cluster number: " << mincl << "\n";
  std::cerr << "Maximumn cluster number: " << maxcl << "\n";
  ParallelpamStop("The classification array has not 0 as minimum value.\n");
 }

 // The number of clusters (number of medoids)
 indextype nmed=maxcl+1;
 
 if (num_obs!=nearest.size())
  ParallelpamStop("Different number of points in the array of classes and in the dissimilarity matrix.\n");
 
 if (DEB & DEBPP)
  std::cout << num_obs << " points classified in " << nmed << " classes.\n";
 
 silinfo sdummy;
 std::vector<silinfo> silres;
 
 // hist is the histogram to know how many individuals are assigned to each cluster
 std::vector<unsigned long> hist;
 hist.resize(nmed,0);
 for (indextype q=0; q<num_obs; q++)
 {
  hist[nearest[q]]++;
  // The structure with the silhouette info for each point is initialized: original point number and its cluster number, which will not be changed.
  sdummy.pnum=q;
  sdummy.ownclus=nearest[q];
  // Also, the nearest (other than its own) cluster and silhouette value are initalized to control values, to be filled.
  sdummy.neiclus=nmed;            // Absurd values to serve as control..
  sdummy.silvalue=std::numeric_limits<siltype>::max();
  silres.push_back(sdummy);
 }
  
 std::vector<siltype> current_sil;
 current_sil.resize(num_obs,siltype(0));
 
 if (nt==1)
    SilhouetteSerial(num_obs,nmed,nearest,current_sil,hist,silres,D);   
 else
 {
    SilhoutteThread_IO<disttype> *silargs = new SilhoutteThread_IO<disttype> [nt];
    for (unsigned int t=0; t<nt; t++)
    {
        silargs[t].num_obs=num_obs;
        silargs[t].nmed=nmed;
        silargs[t].current_sil=&current_sil;
        silargs[t].nearest=&nearest;
        silargs[t].hist=&hist;
        silargs[t].silres=&silres;
        silargs[t].D = &D;
    }
    CreateAndRunThreadsWithDifferentArgs(nt,SilhoutteThread<disttype>,(void *)silargs,sizeof(SilhoutteThread_IO<disttype>));
    
    delete[] silargs;
 }
 Dt.EndClock(DEB & DEBPP); 

 std::vector<siltype> ret(num_obs);

 for (size_t t=0; t<current_sil.size(); t++)
  ret[t]=current_sil[t];
   
 return(ret);
}

template std::vector<siltype> CalculateSilhouette(std::vector<indextype> cl,SymmetricMatrix<float> &D,unsigned int nt);
template std::vector<siltype> CalculateSilhouette(std::vector<indextype> cl,SymmetricMatrix<double> &D,unsigned int nt);

template <typename disttype>
siltype CalculateMeanSilhouette(std::vector<indextype> cl,indextype nmed,SymmetricMatrix<disttype> *D,unsigned int nt)
{
 indextype num_obs=D->GetNRows();

 silinfo sdummy;
 std::vector<silinfo> silres;

 // hist is the histogram to know how many individuals are assigned to each cluster
 std::vector<unsigned long> hist;
 hist.resize(nmed,0);
 for (indextype q=0; q<num_obs; q++)
 {
  hist[cl[q]]++;
  // The structure with the silhouette info for each point is initialized: original point number and its cluster number, which will not be changed.
  sdummy.pnum=q;
  sdummy.ownclus=cl[q];
  // Also, the nearest (other than its own) cluster and silhouette value are initalized to control values, to be filled.
  sdummy.neiclus=nmed;            // Absurd values to serve as control..
  sdummy.silvalue=std::numeric_limits<siltype>::max();
  silres.push_back(sdummy);
 }

 std::vector<siltype> current_sil;
 current_sil.resize(num_obs,siltype(0));

 if (nt==1)
    SilhouetteSerial(num_obs,nmed,cl,current_sil,hist,silres,(*D));
 else
 {
    SilhoutteThread_IO<disttype> *silargs = new SilhoutteThread_IO<disttype> [nt];
    for (unsigned int t=0; t<nt; t++)
    {
        silargs[t].num_obs=num_obs;
        silargs[t].nmed=nmed;
        silargs[t].current_sil=&current_sil;
        silargs[t].nearest=&cl;
        silargs[t].hist=&hist;
        silargs[t].silres=&silres;
        silargs[t].D = D;
    }
    CreateAndRunThreadsWithDifferentArgs(nt,SilhoutteThread<disttype>,(void *)silargs,sizeof(SilhoutteThread_IO<disttype>));

    delete[] silargs;
 }

 siltype ret=0.0;

 for (size_t t=0; t<current_sil.size(); t++)
  ret+=current_sil[t];
 ret/=siltype(num_obs);

 return(ret);
}

template siltype CalculateMeanSilhouette(std::vector<indextype> cl,indextype nmed,SymmetricMatrix<float> *D,unsigned int nt);
template siltype CalculateMeanSilhouette(std::vector<indextype> cl,indextype nmed,SymmetricMatrix<double> *D,unsigned int nt);


