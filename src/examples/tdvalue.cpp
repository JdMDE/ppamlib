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
 * @file tdvalue.cpp
 * @brief <h2>tdvalue</h2>
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
 *        g++ -Wall tdvalue.cpp -o tdvalue -ljmatrix -lppam
 *
*/
#include <jmatrixlib/fullmatrix.h>
#include "../headers/debugpar_ppam.h"
#include "../headers/gettd.h"

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include <iostream>
#include <cstdlib>

extern unsigned char DEB;

using namespace std;

void Usage(char *pname,string error)
{
 cerr << "Usage:\n\n" << "  " << pname << " med_file class_file ds_file\n\n";
 cerr << "  where, if n is the number of points and k the number of medoids,\n\n";
 cerr << "   med_file:    File with the indexes of the medoids in jmatrix format. Compulsory\n";
 cerr << "                It must be a (k x 1) full matrix (column vector) of indextype (unsigned int)\n";
 cerr << "   class_file   File with the number (from 0 to k-1) of the medoid each point is closest to. Compulsory\n";
 cerr << "                It must be a (n x 1) full matrix (column vector) of indextype (unsigned int)\n";
 cerr << "   ds_file:     File with the dissimilarity matrix in jmatrix format. Compulsory\n";
 cerr << "                It must be a symmetric matrix of float or double with dimension (n x n).\n\n";
 cerr << "   The only output will be a double number written in the screen (unless you call the program as tdvalued or tdvaluedd for debugging).\n";
 cerr << "   Points are assumed to be in the same order in the dissimilarity matrix and the classification vector.\n";
 cerr << "   The program will refuse to load the dissimilarity matrix if not enough RAM is available.\n";
 cerr << "   also, it will show a warning if the required amount of memory to load it is above 75% of the available RAM.\n";

 if (error.length()>0)
  cerr << "Error was: " << error << "\n\n";

 exit(1);
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
 * <h2>tdvalue</h2>
 * A program to obtain the value of the TD optimization value of a clustering result.\n
 * TD is defined as the sum of distances of each point to its closest medoid divided by the total number of points.\n
 * This program takes as all its inputs binary files in jmatrix format.
 *
 * The program must be called as
 *
 *     tdvalue med_file clas_file ds_file
 *
 * where, if n is the number of points and k the number of medoids,\n
 * \n
 *  <b>med_file</b>:    File with the indexes of the medoids in jmatrix format. Compulsory\n
 *               It must be a (k x 1) full matrix (column vector) of indextype (unsigned int)\n
 * \n
 *  <b>class_file</b>:  File with the number (from 0 to k-1) of the medoid each point is closest to. Compulsory\n
 *               It must be a (n x 1) full matrix (column vector) of indextype (unsigned int)\n
 * \n
 *  <b>ds_file</b>:     File with the dissimilarity matrix in jmatrix format. Compulsory\n
 *               It must be a symmetric matrix of float or double with dimension (n x n).\n
 * \n
 * The only output will be a double number written in the screen (unless you call the program as <b>tdvalued</b> or <b>tdvaluedd</b> for debugging).\n
 * Points are assumed to be in the same order in the dissimilarity matrix and the classification vector.\n
 * The program will refuse to load the dissimilarity matrix if not enough RAM is available.\n
 * Also, it will show a warning if the required amount of memory to load it is above 75% of the available RAM.\n
 *
 */
int main(int argc,char *argv[])
{
 int call=CheckProgName(string(argv[0]),{"tdvalue","tdvalued","tdvaluedd"});
 // if call is 0 (tdvalue) debug is off by default.
 if (call==1)
  ParallelpamSetDebug(true,false);
 if (call==2)
  ParallelpamSetDebug(true,true);

 if (argc==1)
  Usage(argv[0],"");
 if (argc!=4)
  Usage(argv[0],"Incorrect number of arguments.");

 string mfile=string(argv[1]);
 string cfile=string(argv[2]);
 string dfile=string(argv[3]);

 FullMatrix<indextype> Lmed(mfile);
 FullMatrix<indextype> Lclas(cfile);

 unsigned char mtype,ctype,e,md;
 indextype nr,nc;
 MatrixType(dfile,mtype,ctype,e,md,nr,nc);
 if (mtype!=MTYPESYMMETRIC)
  ParallelpamStop("This program can operate only with binary symmetric matrices as dissimilarity matrices.\n");
 if ((ctype!=FTYPE) && (ctype!=DTYPE))
  ParallelpamStop("This program can operate only with binary symmetric matrices with float or double elements as dissimilarity matrices.\n");

 if ( (Lmed.GetNCols()!=1) || (Lclas.GetNCols()!=1) || (Lclas.GetNRows()!=nr))
  ParallelpamStop("Inconsistent dimensions in the vectors or matrix stored in the input files. Check them with jmatrix info <the_file>\n");

 vector<indextype> Lv;
 for (size_t i=0;i<Lmed.GetNRows();i++)
  Lv.push_back(Lmed.Get(i,0));

 vector<indextype> Lc;
 for (size_t i=0;i<Lclas.GetNRows();i++)
  Lc.push_back(Lclas.Get(i,0));

 double td;
 if (ctype==FTYPE)
 {
  SymmetricMatrix<float> D(dfile,true);
  td=GetTD<float>(Lv,Lc,D);
 }
 else
 {
  SymmetricMatrix<double> D(dfile,true);
  td=GetTD<double>(Lv,Lc,D);
 }

 cout << td << "\n";

 return 0;
}

