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

#include <iostream>
#include <cstdlib>

/**
 * @file parpam.cpp
 * @brief <h2>parpam</h2>
 *        See program use in the documention to main() below\n
 * \n
 *        NOTE: The includes in this source file are for compilation of this program as an example together with the library,\n
 *        before the library itself is installed. Once you have installed the library (assuming headers are in\n
 *        /usr/local/include, lib is in /usr/local/lib or in other place included in your compiler search path)\n
 *        you should substitute this by\n
 *\n
 *        #include <parallelpam/debugpar_ppam.h>   etc...\n
 *\n
 *        and compile with something like
 *
 *        g++ -Wall parpam.cpp -o parpam -ljmatrix -lppam
 *
*/
#include "../headers/debugpar_ppam.h"
#include "../headers/threadhelper.h"
#include "../headers/fastpam.h"

#ifndef DOXYGEN_SHOULD_SKIP_THIS
extern unsigned char DEB;

using namespace std;

void Usage(char *pname,string error)
{
 cerr << "Usage:\n\n" << "  " << pname << " ds_file k [-imet method (medoids_file)] [-omet method] [-mit max_iter] [-nt numthreads] -o root_file_name\n\n";
 cerr << "  where\n\n";
 cerr << "   ds_file:     File with the dissimilarity matrix in jmatrix format.\n";
 cerr << "                It must be a symmetric matrix of float or double with dimension (n x n).\n";
 cerr << "                This argument is compulsory and must be the first one after the program name.\n";
 cerr << "   k:           Requested number of medoids (possitive integer number, k<n).\n";
 cerr << "                This argument is compulsory and must be the second one after the program name.\n";
 cerr << "   imet:        Initialization method, which must be one of the strings 'BUILD', 'LAB' or 'PREV'\n";
 cerr << "                If you use PREV the file with the initial medoids must be given, too, which must be\n";
 cerr << "                a jmatrix FullMatrix of unsiged int with dimension (n x 1) (as returned by another call to this program)\n";
 cerr << "                If you use BUILD or LAB no initial medoids file should be provided. Default value: BUILD.\n";
 cerr << "   omet:        Optimization method, which must be one of the strings 'FASTPAM1' or 'TWOBRANCH'. Default value: FASTPAM1\n";
 cerr << "   max_iter:    Maximum number of iterations. Set it to 0 to do only the initialization phase (with BUILD or LAB method).\n";
 cerr << "                Default value: " << MAX_ITER << ".\n";
 cerr << "   numthreads:  Requested number of threads.\n";
 cerr << "                Setting it to 0 will make the program to choose according to the number of processors/cores\n";
 cerr << "                of your machine (default value).\n";
 cerr << "                Setting to -1 forces serial implementation (no threads)\n";
 cerr << "   root_fname:  A string used to build root_fname_med.bin and root_fname_clas.bin.\n";
 cerr << "                This argument is compulsory and must be the last one.\n\n";
 cerr << "   Calling this program as parpamd turns on debugging; calling it as parpamdd turns on the jmatrix library debugging, too.\n";
 cerr << "   The output files will contain jmatrix vectors of final medoids and classification, respectively.\n";
 cerr << "   Both are FullMatrix of indextype (unsigned int) with dimensions (k x 1) for med and (n x 1) for clas.\n";
 cerr << "   The first one contains the indices of the found medoids as row indices of the dissimilarity matrix, from 0.\n";
 cerr << "   (i.e.: integers in range [0..n-1])\n";
 cerr << "   The second contains the index in the first one (from 0) of the medoid to which class each point belongs to.\n";
 cerr << "   (i.e.: integers in range [0..k-1])\n";
 cerr << "   If the dissimilarity matrix contained row names (i.e.: point names) the output vectors will keep them, too.\n";
 cerr << "   The program will refuse to load the dissimilarity matrix if not enough RAM is available;\n";
 cerr << "   also, it will show a warning if the required amount of memory to load it is above 75% of the available RAM.\n";
 cerr << "   Remember that using the program 'jmat csvdump ...' you can convert the output files to .csv format.\n\n";

 if (error.length()>0)
  cerr << "Error was: " << error << "\n\n";

 exit(1);
}

void Verifyk(string kv,int &k)
{
 for (unsigned i=0;i<kv.size();i++)
  if ((kv[i]<'0') || (kv[i]>'9'))
   ParallelpamStop("Argument 'k' must be a possitive integer number.");
 k=atol(kv.c_str());
 if ((indextype)k>=MAX_MEDOIDS)
 {
     ostringstream errst;
     errst << "Asking for too many medoids. Maximum is " << MAX_MEDOIDS-1 << ".\n";
     ParallelpamStop(errst.str());
 }
}

void VerifyInitMethod(vector<string> args,unsigned char &init_method,vector<indextype> &inimeds)
{
 inimeds.clear();
 vector<string>::iterator it=find(args.begin(),args.end(),"-imet");

 if (it==args.end())
 {
  init_method=INIT_METHOD_BUILD;
  return;
 }

 string imethod=*(it+1);

 if ((imethod!="BUILD") && (imethod!="LAB") && (imethod!="PREV"))
  ParallelpamStop("Initializetion method must be BUILD, LAB or PREV.");

 if (imethod=="PREV")
 {
  string inimed_file=*(it+2);
  if (inimed_file[0]=='-')
   ParallelpamStop("Initialization method PREV must be followed by a file name (which cannot start with '-').");
  init_method=INIT_METHOD_PREVIOUS;
  unsigned char mtype,ctype,e,md;
  indextype nr,nc;
  MatrixType(inimed_file,mtype,ctype,e,md,nr,nc);
  cout << "**** nr=" << nr << "; nc=" << nc << "\n";
  // WARNING: this might fail if definition of indextype is changed...
  if ((mtype!=MTYPEFULL) || (ctype!=UITYPE) || (nc!=1))
   ParallelpamStop("The file of initial medoids is wrong. It must contain a FullMatrix of unsigned ints with just one column.\n");
  if (nr==0)
   ParallelpamStop("The file of initial medoids is empty. Check how it was created.\n");
  FullMatrix<indextype> V(inimed_file);
  for (size_t i=0;i<V.GetNRows();i++)
   inimeds.push_back(V.Get(i,0));
  if (DEB & DEBPP)
   cout << "Initial medoids loaded from file " << inimed_file << ".\n";
  return;
 }


 init_method = (imethod=="LAB") ? INIT_METHOD_LAB : INIT_METHOD_BUILD;
 return;
}

void VerifyOptMethod(vector<string> args,unsigned char &opt_method)
{
 vector<string>::iterator it=find(args.begin(),args.end(),"-omet");

 if (it==args.end())
 {
  opt_method=OPT_METHOD_FASTPAM1;
  return;
 }

 string omethod=*(it+1);

 if ((omethod!="FASTPAM1") && (omethod!="TWOBRANCH"))
  ParallelpamStop("Method must be FASTPAM1 or TWOBRANCH.");

 opt_method = (omethod=="FASTPAM1") ? OPT_METHOD_FASTPAM1 : OPT_METHOD_FASTPAMBSIL;
}

void VerifyMaxIter(vector<string> args,int &max_iter)
{
 vector<string>::iterator it=find(args.begin(),args.end(),"-mit");
 if (it==args.end())
 {
  max_iter=MAX_ITER-1;
  return;
 }

 string its=*(it+1);
 for (size_t i=0;i<its.length();i++)
  if ((its[i]<'0') || (its[i]>'9'))
   ParallelpamStop("Argument -mit must be followed by a possitive integer number.");

 max_iter=atoi(its.c_str());
 if ((unsigned int)max_iter>MAX_ITER)
 {
     ostringstream errst;
     errst << "Asking for too many limit iterations. Maximum is " << MAX_ITER-1 << ".\n";
     errst << "If you need more, change the constant MAX_ITER at fastpam.h and reinstall the package.\n";
     ParallelpamStop(errst.str());
 }
}

void VerifyNThreads(vector<string> args,unsigned int &nt)
{
 int nthreads;
 vector<string>::iterator it=find(args.begin(),args.end(),"-nt");
 if (it!=args.end())
 {
  string nts=*(it+1);
  for (size_t i=0;i<nts.length();i++)
   if (nts[i]!='-')
    if ((nts[i]<'0') || (nts[i]>'9'))
     ParallelpamStop("Argument -nt must be followed by a number (may be negative for no threads).");
  nthreads=atoi(nts.c_str());
 }
 else
  nthreads=0;

 nt=ChooseNumThreads(nthreads);
 if (DEB & DEBPP)
  std::cout << nt << " threads will be used.\n";
}

void ParseArguments(int argc,char *argv[],
                    string &dissim_file,
                    int &k,
                    unsigned char &init_method,
                    vector<indextype> &inimeds,
                    unsigned char &opt_method,
                    int &max_iter,
                    unsigned int &nt,
                    string &mfile,
                    string &cfile)
{
 if (argc==1)
  Usage(argv[0],"");
 if ((argc<5) || (argc>14))
  Usage(argv[0],"Incorrect number of arguments.");

 dissim_file=string(argv[1]);

 if (string(argv[argc-2])!="-o")
  Usage(argv[0],"Last but one argument must be -o.");

 string res_rname=string(argv[argc-1]);
 size_t wheredot=res_rname.find(".");
 if (wheredot==string::npos)
 {
  mfile=res_rname+"_med.bin";
  cfile=res_rname+"_clas.bin";
 }
 else
 {
  mfile=res_rname.substr(0,wheredot)+"_med"+res_rname.substr(wheredot);
  cfile=res_rname.substr(0,wheredot)+"_clas"+res_rname.substr(wheredot);
 }

 Verifyk(string(argv[2]),k);

 vector<string> args;
 for (int i=3;i<argc-2;i++)
  args.push_back(string(argv[i]));

 VerifyInitMethod(args,init_method,inimeds);

 VerifyOptMethod(args,opt_method);

 VerifyMaxIter(args,max_iter);

 VerifyNThreads(args,nt);
}

void NameChanged(vector<string> ends)
{
 cerr << "You have changed the name of this program. Don't do that. Its name must be (or at least, must end in) ";
 for (size_t j=0;j<ends.size();j++)
  cerr << "'" << ends[j] << "' ";
 cerr << "\n";
 exit(1);
}

int CheckProgName(string pname,vector<string> possible_endings)
{
 sort(possible_endings.begin(),possible_endings.end(),[](string a, string b) { return a.size()<b.size(); });
 size_t i=0;
 while (i<possible_endings.size())
 {
  if (pname.size()<possible_endings[i].size())
   NameChanged(possible_endings);
  if (pname.substr(pname.size()-possible_endings[i].size())==possible_endings[i])
   return i;
  i++;
 }
 NameChanged(possible_endings);
 return -1;  // Just to avoid a warning
}

#endif

/**
 * <h2>parpam</h2>
 * A program to apply the Partitionin Around Medoids (PAM) clustering method to a set of individuals whose
 * dissimilarity matrix is given, in parallel. It implements the FASTPAM1 algorithm described in\n
 * \n
 * Schubert, E. and Rousseeuw, P.J.: "Fast and eager k-medoids clustering: O(k) runtime improvement of the PAM, CLARA, and CLARANS algorithms."\n
 * Information Systems, vol. 101, p. 101804, 2021.\n
 * doi: https://doi.org/10.1016/j.is.2021.101804\n
 * \n
 * See documentation of class FastPAM for more information\n
 *
 * The program must be called as
 *
 * parpam ds_file k [-imet method (medoids_file)] [-omet method] [-mit max_iter] [-nt numthreads] -o root_file_name
 *
 * where\n
 * \n
 * <b>ds_file</b>:     File with the dissimilarity matrix in jmatrix format.\n
 *              It must be a symmetric matrix of float or double with dimension (n x n). You can use program pardis.cpp to generate it.\n
 *              This argument is compulsory and must be the first one after the program name.\n
 * \n
 * <b>k</b>:           Requested number of medoids (possitive integer number, k<n).\n
 *              This argument is compulsory and must be the second one after the program name.\n
 * \n
 * <b>imet</b>:        Initialization method, which must be one of the strings 'BUILD', 'LAB' or 'PREV'\n
 *              If you use PREV the file with the initial medoids must be given, too, which must be\n
 *              a jmatrix FullMatrix of unsiged int with dimension (n x 1) (as returned by another call to this program)\n
 *              If you use BUILD or LAB no initial medoids file should be provided. Default value: BUILD.\n
 * \n
 * <b>omet</b>:        Optimization method, which must be one of the strings 'FASTPAM1' or 'TWOBRANCH'. Default value: FASTPAM1\n
 * \n
 * <b>max_iter</b>:    Maximum number of iterations. Set it to 0 to do only the initialization phase (with BUILD or LAB method).\n
 *              Default value: the value of constant MAX_ITER defined in fastpam.h\n
 * \n
 * <b>numthreads</b>:  Requested number of threads.\n
 *              Setting it to 0 will make the program to choose according to the number of processors/cores of your machine (default value).\n
 *              Setting to -1 forces serial implementation (no threads)\n
 * \n
 * <b>root_fname</b>:  A string used to build root_fname_med.bin and root_fname_clas.bin.\n
 *              This argument is compulsory and must be the last one.\n
 * \n
 * Calling this program as <b>parpamd</b> turns on debugging; calling it as <b>parpamdd</b> turns on the jmatrix library debugging, too.\n
 * The output files will contain jmatrix vectors of final medoids and classification, respectively.\n
 * Both are FullMatrix of indextype (unsigned int) with dimensions (k x 1) for med and (n x 1) for clas.\n
 * The first one contains the indices of the found medoids as row indices of the dissimilarity matrix, from 0 (i.e.: integers in range [0..n-1]).\n
 * The second contains the index in the first one (from 0) of the medoid to which class each point belongs to (i.e.: integers in range [0..k-1]).\n
 * If the dissimilarity matrix contained row names (i.e.: point names) the output vectors will keep them, too.\n
 * The program will refuse to load the dissimilarity matrix if not enough RAM is available;\n
 * also, it will show a warning if the required amount of memory to load it is above 75% of the available RAM.\n
 * Remember that using the program 'jmat csvdump ...' you can convert the output files to .csv format.
 *
 */

int main(int argc,char *argv[])
{
 int call=CheckProgName(string(argv[0]),{"parpam","parpamd","parpamdd"});
 // if call is 0 (pardis) debug is off by default.
 if (call==1)
  ParallelpamSetDebug(true,false);
 if (call==2)
  ParallelpamSetDebug(true,true);

 string dissim_file;
 int k;
 unsigned char init_method;
 vector<indextype> inimeds;
 unsigned char opt_method;
 int max_iter;
 unsigned int nt;
 string mfile,cfile;

 ParseArguments(argc,argv,dissim_file,k,init_method,inimeds,opt_method,max_iter,nt,mfile,cfile);

 if (DEB & DEBPP)
 {
  cout << "Applying PAM with arguments:\n";
  cout << "  Dissimilarity file: " << dissim_file << "\n";
  cout << "  Number of medoids: " << k << "\n";
  cout << "  Intialization method: ";
  switch (init_method)
  {
   case INIT_METHOD_BUILD: cout << "BUILD\n"; break;
   case INIT_METHOD_LAB: cout << "LAB\n"; break;
   case INIT_METHOD_PREVIOUS: cout << "PREV\n";
   default: break;
  }
  cout << "  Optimization method: ";
  switch (opt_method)
  {
   case OPT_METHOD_FASTPAM1: cout << "FASTPAM1\n"; break;
   case OPT_METHOD_FASTPAMBSIL: cout << "FASTPAMBSIL\n"; break;
   default: break;
  }
  cout << "  Maximum number of iterations: " << max_iter << ((max_iter==0) ? " (only initial phase)\n" : "\n");
  cout << "  Number of threads: " << nt;
  cout << "  Medoid indices will be stored in file " << mfile << ".\n";
  cout << "  Clasification will be stored in file " << cfile << ".\n";
 }

 unsigned char mtype,ctype,e,md;
 indextype nr,nc;
 MatrixType(dissim_file,mtype,ctype,e,md,nr,nc);
 if (mtype!=MTYPESYMMETRIC)
  ParallelpamStop("This function can operate only with binary symmetric matrices.\n");
 if ((ctype!=FTYPE) && (ctype!=DTYPE))
  ParallelpamStop("This function can operate only with binary symmetric matrices with float or double elements.\n");

 if (DEB & DEBPP)
 {
  std::cout << "Reading symmetric distance/dissimilarity matrix " << dissim_file << "\n";
  std::cout.flush();
 }

 if (ctype==FTYPE)
 {
  SymmetricMatrix<float> D(dissim_file,true);

  FastPAM<float> FP(&D,k,init_method,max_iter,nt);
  FP.Init(inimeds,nt);
  FP.Run(opt_method,nt);

  FullMatrix<indextype> &Lmed=FP.GetMedoids(D.GetRowNames());
  Lmed.WriteBin(mfile);

  FullMatrix<indextype> &Lclasif=FP.GetAssign(D.GetRowNames());
  Lclasif.WriteBin(cfile);
 }
 else
 {
  SymmetricMatrix<double> D(dissim_file);

  FastPAM<double> FP(&D,k,init_method,max_iter,nt);
  FP.Init(inimeds,nt);
  FP.Run(opt_method,nt);

  FullMatrix<indextype> &Lmed=FP.GetMedoids(D.GetRowNames());
  Lmed.WriteBin(mfile);

  FullMatrix<indextype> &Lclasif=FP.GetAssign(D.GetRowNames());
  Lclasif.WriteBin(cfile);
 }
}

