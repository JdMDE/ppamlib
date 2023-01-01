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

/**
 * @file parsil.cpp
 * @brief <h2>parsil</h2>
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
 *        g++ -Wall parsil.cpp -o parsil -ljmatrix -lppam
 *
*/
#include "../headers/debugpar_ppam.h"
#include "../headers/threadhelper.h"
#include "../headers/fastpam.h"
#include "../headers/silhouette.h"

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include <iostream>
#include <cstdlib>

extern unsigned char DEB;

using namespace std;

void Usage(char *pname,string error)
{
 cerr << "Usage:\n\n" << "  " << pname << " dissim_file clasif_file [-nt numthreads] -o out_file_name\n\n";
 cerr << "  where\n\n";
 cerr << "   dissim_file:    File with the dissimilarity matrix in jmatrix format.\n";
 cerr << "                   It must be a SymmetricMatrix of float or double with dimension (n x n).\n";
 cerr << "   clasif_file:    File with the clasification result, as obtained from program parpam.\n";
 cerr << "                   It must be a (n x 1) matrix (a column vector) of unsigned int values with values in 0..(k-1) being k the number of clusters.\n";
 cerr << "   numthreads:     Requested number of threads.\n";
 cerr << "                   Setting it to 0 will make the program to choose according to the number of processors/cores of your machine (default value).\n";
 cerr << "                   Setting to -1 forces serial implementation (no threads)\n";
 cerr << "   out_file_name:  Name of the file contaning the silhouette. Compulsory.\n\n";
 cerr << "   The output file will be a FullMatrix of double type and dimension (n x 1) (a column vector) with the value of the silhouette for each point.\n";
 cerr << "   Points are assumed to be in the same order in the dissimilarity matrix and the clasification vector, and this is the order in which their\n";
 cerr << "   silhouettes will be written in the output vector. If the matrix has row names, they will be set for the output file. It the clasif vector\n";
 cerr << "   has row names, they will be checked against the row names of the matrix, if both are present. If only clasification vector has names,\n";
 cerr << "   they will be set for the output vector.\n";
 cerr << "   The program will refuse to load the dissimilarity matrix if not enough RAM is available; also, it will show a warning if the required amount\n";
 cerr << "   of memory to load it is above 75% of the available RAM.\n";
 cerr << "   Remember that using the program 'jmat csvdump ...' you can convert the output file to .csv format.\n\n";
 if (error.length()>0)
  cerr << "Error was: " << error << "\n\n";

 exit(1);
}

vector<string> CheckNameConsistency(vector<string> dn,vector<string> cn,size_t n)
{
 vector<string> ret;

 if ( dn.size()==0 && cn.size()==0 )
  return ret;

 if ( dn.size()>0  && cn.size()==0 )
  ret=dn;
 else
 {
  if ( dn.size()==0 && cn.size()>0 )
   ret=cn;
  else
  {
   if ( dn.size()!=cn.size() )
    ParallelpamStop("The lengths of names in dissimilarity matrix and classification file are not the same.\n");
   if ( dn != cn )
    ParallelpamStop("The point names in dissimilarity matrix and in classification file are not equal.\n");
   ret=dn;
  }
 }
 if (ret.size() != n)
  ParallelpamStop("The lengths of names in dissimilarity matrix and classification file are not the length of the returned silhouette vector.\n");

 return ret;
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
 * <h2>parsil</h2>
 * A program to calculate the silhouette of a clustering (usually obtained with parpam) in parallel
 *
 * The program must be called as
 *
 * parsil dissim_file clasif_file [-nt numthreads] -o out_file_name
 *
 *
 * where\n
 * \n
 *  <b>dissim_file</b>:    File with the dissimilarity matrix in jmatrix format.\n
 *                  It must be a SymmetricMatrix of float or double with dimension (n x n).\n
 * \n
 *  <b>clasif_file</b>:    File with the clasification result, as obtained from program parpam.\n
 *                  It must be a (n x 1) matrix (a column vector) of unsigned int values with values in 0..(k-1) being k the number of clusters.\n
 * \n
 *  <b>numthreads</b>:     Requested number of threads.\n
 *                  Setting it to 0 will make the program to choose according to the number of processors/cores of your machine (default value).\n
 *                  Setting to -1 forces serial implementation (no threads)\n
 * \n
 *  <b>out_file_name</b>:  Name of the file contaning the silhouette. Compulsory.\n
 *  The output file will be a FullMatrix of double type and dimension (n x 1) (a column vector) with the value of the silhouette for each point.\n
 * \n
 *  Points are assumed to be in the same order in the dissimilarity matrix and the clasification vector, and this is the order in which their\n
 *  silhouettes will be written in the output vector. If the matrix has row names, they will be set for the output file. It the clasif vector\n
 *  has row names, they will be checked against the row names of the matrix, if both are present. If only clasification vector has names,\n
 *  they will be set for the output vector.\n
 *  The program will refuse to load the dissimilarity matrix if not enough RAM is available; also, it will show a warning if the required amount\n
 *  of memory to load it is above 75% of the available RAM.\n
 *  Remember that using the program 'jmat csvdump ...' you can convert the output file to .csv format.\n
 *
 */
int main(int argc,char *argv[])
{
 int call=CheckProgName(string(argv[0]),{"parsil","parsild","parsildd"});
 // if call is 0 (parsil) debug is off by default.
 if (call==1)
  ParallelpamSetDebug(true,false);
 if (call==2)
  ParallelpamSetDebug(true,true);

 if (argc==1)
  Usage(argv[0],"");
 if ((argc!=5) && (argc!=7))
  Usage(argv[0],"Incorrect number of arguments.");

 string dfile=string(argv[1]);
 string cfile=string(argv[2]);

 if (string(argv[argc-2])!="-o")
  Usage(argv[0],"Last but one argument must be -o.");

 string outname=string(argv[argc-1]);

 int nthreads;
 if (argc==5)
  nthreads=0;
 else
 {
  if (string(argv[3])!="-nt")
   Usage(argv[0],"Using the program with six or seven arguments, but the third or fourth one is not -nt.");
  string nts=string(argv[4]);
  for (size_t i=0;i<nts.length();i++)
   if (nts[i]!='-')
    if ((nts[i]<'0') || (nts[i]>'9'))
     Usage(argv[0],"Argument -nt must be followed by a number (may be negative for no threads).");
  nthreads=atoi(nts.c_str());
 }
 unsigned int nt=ChooseNumThreads(nthreads);

 if (DEB & DEBPP)
 {
  cout << "Calculating silhouette with arguments:\n";
  cout << "  Dissimilarity file: " << dfile << "\n";
  cout << "  Classification file: " << cfile << "\n";
  cout << "  Number of threads: " << nt;
  cout << "  Output file: " << outname << "\n";
 }

 unsigned char mtype,ctype,e,md;
 indextype nr,nc;
 MatrixType(dfile,mtype,ctype,e,md,nr,nc);
 if (mtype!=MTYPESYMMETRIC)
  ParallelpamStop("This program can operate only with binary symmetric matrices as dissimilarity matrices.\n");
 if ((ctype!=FTYPE) && (ctype!=DTYPE))
  ParallelpamStop("This program can operate only with binary symmetric matrices with float or double elements as dissimilarity matrices.\n");

 FullMatrix<indextype> Lclas(cfile);
 if ((Lclas.GetNRows()!=nr) || (Lclas.GetNCols()!=1))
  ParallelpamStop("Inconsistent dimensions in the vector or matrix stored in the input files. Check them with jmatrix info <the_file>\n");
 vector<string> Cnames=Lclas.GetRowNames();

 vector<indextype> Lc;
 for (size_t i=0;i<Lclas.GetNRows();i++)
  Lc.push_back(Lclas.Get(i,0));

 vector<siltype> sil;
 vector<string> Dnames;
 if (ctype==FTYPE)
 {
  SymmetricMatrix<float> D(dfile,true);
  sil=CalculateSilhouette<float>(Lc,D,nt);
  Dnames=D.GetRowNames();
 }
 else
 {
  SymmetricMatrix<double> D(dfile,true);
  sil=CalculateSilhouette<double>(Lc,D,nt);
  Dnames=D.GetRowNames();
 }
 FullMatrix<double> Vsil(sil.size(),1);
 for (size_t i=0;i<sil.size();i++)
  Vsil.Set(i,0,sil[i]);

 vector<string> names=CheckNameConsistency(Dnames,Cnames,sil.size());
 if (names.size()>0)
  Vsil.SetRowNames(names);

 Vsil.WriteBin(outname);

 return 0;
}
