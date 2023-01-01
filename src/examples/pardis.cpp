/* Copyright (C) 2022 Juan Domingo (Juan.Domingo@uv.es)
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
 * @file pardis.cpp
 * @brief <h2>pardis</h2>
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
 *        g++ -Wall pardis.cpp -o pardis -ljmatrix -lppam
 *
*/
#include "../headers/fastpam.h"
#include "../headers/debugpar_ppam.h"
#include "../headers/threadhelper.h"
#include "../headers/dissimmat.h"

#ifndef DOXYGEN_SHOULD_SKIP_THIS
extern unsigned char DEB;

using namespace std;

void Usage(char *pname,string error)
{
 cerr << "Usage:\n\n" << "  " << pname << " input_file [-dis distype] [-vtype valuetype] [-nt numthreads] [-com comment] -o out_file_name\n\n";
 cerr << "  where\n\n";
 cerr << "   input_file:     File with the input matrix in jmatrix format.\n";
 cerr << "                   It must be a matrix of float or double with dimension (n x p) where the individuals (points/vectors,\n";
 cerr << "                   which are n) must be the rows and components/dimensions (which are p) must be the columns.\n";
 cerr << "                   This argument is compulsory and must be immediately after the program name.\n";
 cerr << "   dis:            Type of metrics/dissimilarity, which must be one of the strings 'L1' (Manhattan), 'L2' (Euclidean)\n";
 cerr << "                   or 'Pe' (Pearson dissimilarity). Default: L2.\n";
 cerr << "   vtype:          Data type for the output dissimilarity/distance matrix.\n";
 cerr << "                   It must be one of the strings 'float' or 'double'. Default: float.\n";
 cerr << "   numthreads:     Requested number of threads.\n";
 cerr << "                   Setting it to 0 will make the program to choose according to the number of processors/cores\n";
 cerr << "                   of your machine (default value).\n";
 cerr << "                   Setting to -1 forces serial implementation (no threads)\n";
 cerr << "   comment         Comment to be attached to the dissimilarity matrix. Default: no comment will be added.\n";
 cerr << "   out_file_name:  Name of the file contaning the dissimilarity matrix as a binary jmatrix.\n";
 cerr << "                   If the input matrix has row names, these names will be copied to the dissimilarity matrix as row names, too.\n";
 cerr << "                   This argument is compulsory and must be the last one.\n\n";
 cerr << "   Calling this program as pardisd turns on debugging; calling it as pardisdd turns on the jmatrix library debugging, too.\n";
 cerr << "   The distance/dissimilarity matrix in the output file will be a SymmetricMatrix of the requested data type and size (n x n).\n";
 cerr << "   The used memory is quadratic with n (concretely, n*(n+1)/2) so it can be very big.\n";
 cerr << "   The program refuses to create it if not enough RAM is available, and shows a warning\n";
 cerr << "   if the required amount of memory is above 75% of the available RAM.\n\n";

 if (error.length()>0)
  cerr << "Error was: " << error << "\n\n";

 exit(1);
}

// Looks if the name of the input exists as a file and contains a valid input matrix.
// If so, finds out its type (full, sparse, symmetric) and its value type
// To be valid it must be either full or sparse and store either floats or doubles
void VerifyInputMatrix(string inpname,unsigned char &imattype,unsigned char &imatvaltype)
{
 unsigned char e,md;
 indextype nr,nc;
 MatrixType(inpname,imattype,imatvaltype,e,md,nr,nc);

 if (DEB & DEBJM)
 {
  std::cout << "Input matrix is ";
  switch (imattype)
  {
   case MTYPEFULL: if (DEB & DEBJM)
                    std::cout << "a full matrix ";
                   break;
   case MTYPESPARSE: if (DEB & DEBJM)
                      std::cout << "a sparse matrix ";
                     break;
   case MTYPESYMMETRIC: if (DEB & DEBJM)
                        std::cout << "a symmetric matrix. This is not allowed; it must be full or sparse.\n";
                        ParallelpamStop("Invalid matrix type.\n");
                        break;
   default: if (DEB & DEBJM)
             std::cout << "of unknown type (neither full, sparse of symmetric). Was it created with jmatrix?\n";
            ParallelpamStop("Unknown matrix type.\n");
            break;
  }
  switch (imatvaltype)
  {
   case FTYPE: if (DEB & DEBJM)
               std::cout << " with elements of type 'float' and size (" << nr << "," << nc << ")\n";
               break;
   case DTYPE: if (DEB & DEBJM)
                std::cout << " with elements of type 'double' and size (" << nr << "," << nc << ")\n";
               break;
   default: if (DEB & DEBJM)
             std::cout << " with elements which are neither 'float' nor 'double'. This is not allowed to calculate dissimilarity matrix. Sorry.\n";
            ParallelpamStop("Data type of input matrix not allowed.\n");
            break;
  }
 }
}

void VerifyDistanceType(vector<string> args,unsigned char &dtype)
{
 vector<string>::iterator it=find(args.begin(),args.end(),"-dis");
 string distype;
 if (it!=args.end())
 {
  distype=*(it+1);
  if ((distype!="L1") && (distype!="L2") && (distype!="Pe"))
   ParallelpamStop("Distance/dissimilarity type (value following -dis argument) must be L1, L2 or Pe.");
 }
 else
  distype="L2";

 if (distype=="L1")
  dtype=DL1;
 if (distype=="L2")
  dtype=DL2;
 if (distype=="Pearson")
  dtype=DPe;

 if (DEB & DEBPP)
 {
  std::cout << "Used distance is ";
  switch (dtype)
  {
    case DL1: std::cout << "L1 (Manhattan).\n"; break;
    case DL2: std::cout << "L2 (Euclidean).\n"; break;
    case DPe: std::cout << "Pearson dissimilarity.\n"; break;
    default: std::cout << "unknown?\n"; break;
  }
 }
}

void VerifyOutputValueType(vector<string> args,unsigned char &vrestype)
{
 string vtype;
 vector<string>::iterator it=find(args.begin(),args.end(),"-vtype");
 if (it!=args.end())
 {
  vtype=*(it+1);
  if ((vtype!="float") && (vtype!="double"))
   ParallelpamStop("Value type of output (value following -vtype aregument) must be float or double.");
 }
 else
  vtype="float";

 vrestype = (vtype=="float") ? FTYPE : DTYPE;
 if (DEB & DEBPP)
  std::cout << "Output distance/dissimilarity matrix will contain values of type " << ((vrestype==FTYPE) ? "float.\n" : "double.\n");
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

void VerifyComment(vector<string> args,string &comment)
{
 vector<string>::iterator it=find(args.begin(),args.end(),"-com");
 if (it!=args.end())
  comment=*(it+1);
 else
  comment="";
 if (DEB & DEBPP)
 {
  if (comment=="")
   std::cout << "No comment will be attached to output matrix.\n";
  else
   std::cout << "The comment '" << comment << "' will be attached to output matrix.\n";
 }
}

void ParseArguments(int argc,char *argv[],string &inpname,unsigned char &imattype,unsigned char &imatvaltype,
                    string &outname,unsigned char &dtype,unsigned char &vrestype,unsigned int &nt,
                    string &comment)
{
 if (argc==1)
  Usage(argv[0],"");
 if ((argc<4) || (argc>11))
  Usage(argv[0],"Incorrect number of arguments.");

 inpname=string(argv[1]);

 if (string(argv[argc-2])!="-o")
  Usage(argv[0],"Last but one argument must be -o.");

 outname=string(argv[argc-1]);

 VerifyInputMatrix(inpname,imattype,imatvaltype);

 vector<string> args;
 for (int i=2;i<argc-2;i++)
  args.push_back(string(argv[i]));

 VerifyDistanceType(args,dtype);

 VerifyOutputValueType(args,vrestype);

 VerifyNThreads(args,nt);

 VerifyComment(args,comment);
}

template<typename ivaltype,typename ovaltype>
SymmetricMatrix<ovaltype> &CalcDist(bool input_is_full,string iname,unsigned char disttype,unsigned int nt)
{
 if (input_is_full)
 {
  FullMatrix<ivaltype> M(iname);
  if (DEB & DEBPP)
  {
   std::cout << "Read full matrix from file " << iname << ". ";
   std::cout << "Its size is [" << M.GetNRows() << " x " << M.GetNCols() << "] and it uses " << M.GetUsedMemoryMB() << " MBytes.\n";
  }
  return CalcDistFromFull<ivaltype,ovaltype>(M,disttype,nt);
 }
 else
 {
  SparseMatrix<ivaltype> M(iname);
  if (DEB & DEBPP)
  {
   std::cout << "Read sparse matrix from file " << iname << ". ";
   std::cout << "Its size is [" << M.GetNRows() << " x " << M.GetNCols() << "] and it uses " << M.GetUsedMemoryMB() << " MBytes.\n";
  }
  return CalcDistFromSparse<ivaltype,ovaltype>(M,disttype,nt);
 }
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
 * <h2>pardis</h2>
 * A program to calculate the distance/dissimilarity matrix between the rows of an input matrix considering
 * each row as a vector/individual and each column as a dimension/feature.
 *
 * The program must be called as
 *
 * pardis input_file [-dis distype] [-vtype valuetype] [-nt numthreads] [-com comment] -o out_file_name
 *
 * where\n
 * \n
 * <b>input_file</b>:     File with the input matrix in jmatrix format.\n
 *                 It must be a matrix of float or double with dimension (n x p) where the individuals (points/vectors,\n
 *                 which are n) must be the rows and components/dimensions (which are p) must be the columns.\n
 *                 Remember that you can use the program 'jmatrix csvread ...' to create this file from a .csv table\n
 *                 This argument is compulsory and must be immediately after the program name.\n
 * \n
 * <b>dis</b>:            Type of metrics/dissimilarity, which must be one of the strings 'L1' (Manhattan), 'L2' (Euclidean)\n
 *                 or 'Pe' (Pearson dissimilarity). Default: L2.\n
 * \n
 * <b>vtype</b>:          Data type for the output dissimilarity/distance matrix.\n
 *                 It must be one of the strings 'float' or 'double'. Default: float.\n
 * \n
 * <b>numthreads</b>:     Requested number of threads.\n
 *                 Setting it to 0 will make the program to choose according to the number of processors/cores of your machine (default value).\n
 *                 Setting to -1 forces serial implementation (no threads)\n
 * \n
 * <b>comment</b>:            Comment to be attached to the dissimilarity matrix. Default: no comment will be added.\n
 * \n
 * <b>out_file_name</b>:  Name of the file contaning the dissimilarity matrix as a binary jmatrix.\n
 *                 If the input matrix has row names, these names will be copied to the dissimilarity matrix as row names, too.\n
 *                 This argument is compulsory and must be the last one.\n
 * \n
 * Calling this program as <b>pardisd</b> turns on debugging; calling it as <b>pardisdd</b> turns on the jmatrix library debugging, too.\n
 * The distance/dissimilarity matrix in the output file will be a SymmetricMatrix of the requested data type and size (n x n).\n
 * The used memory is quadratic with n (concretely, n*(n+1)/2) so it can be very big.\n
 * The program refuses to create it if not enough RAM is available, and shows a warning\n
 * if the required amount of memory is above 75\% of the available RAM.\n
 *
 */
int main(int argc,char *argv[])
{
 int call=CheckProgName(string(argv[0]),{"pardis","pardisd","pardisdd"});
 // if call is 0 (pardis) debug is off by default.
 if (call==1)
  ParallelpamSetDebug(true,false);
 if (call==2)
  ParallelpamSetDebug(true,true);

 string iname;
 string oname;
 unsigned char imattype,imatvaltype;
 unsigned char disttype;
 unsigned char omatvaltype;
 unsigned int nt;
 string comment;

 ParseArguments(argc,argv,iname,imattype,imatvaltype,oname,disttype,omatvaltype,nt,comment);

 if (omatvaltype==FTYPE)
 {
  SymmetricMatrix<float> &D =
    ((imatvaltype==FTYPE) ? CalcDist<float,float>((imatvaltype==MTYPEFULL),iname,disttype,nt) : CalcDist<double,float>((imatvaltype==MTYPEFULL),iname,disttype,nt));
  if (comment!="")
   D.SetComment(comment);
  D.WriteBin(oname);
 }
 else
 {
  SymmetricMatrix<double> &D =
    ((imatvaltype==FTYPE) ? CalcDist<float,double>((imatvaltype==MTYPEFULL),iname,disttype,nt) : CalcDist<double,double>((imatvaltype==MTYPEFULL),iname,disttype,nt));
  if (comment!="")
   D.SetComment(comment);
  D.WriteBin(oname);
 }
 return 0;
}
