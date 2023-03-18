#include <jmatrixlib/symmetricmatrix.h>
#include "../headers/diftimehelper.h"
#include "../headers/debugpar_ppam.h"
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include <iostream>
#include <cstdlib>

extern unsigned char DEB;

using namespace std;

int main(int argc,char *argv[])
{
    string dfile=string(argv[1]);
    unsigned char mtype,ctype,e,md;
    indextype nr,nc;

    ParallelpamSetDebug(true,true);
    MatrixType(dfile,mtype,ctype,e,md,nr,nc);
    if (mtype!=MTYPESYMMETRIC)
    {
     cerr << "Error: not a symmetric matrix.\n";
     exit(1);
    }
    if (ctype!=FTYPE)
    {
     cerr << "Error: not a matrix of floats.\n";
     exit(1);
    }
    SymmetricMatrix<float> D(dfile,true);

    DifftimeHelper tfmed;
    tfmed.StartClock("\n");
    float s=0.0;

    size_t x;
    for (indextype r=0;r<D.GetNRows();r++)
     for (indextype c=r;c<D.GetNCols();c++)
      s+=D.Get(r,c);
    /*
    unsigned long long nrows=(unsigned long long)D.GetNRows();
    unsigned long long ne=(nrows*(nrows+1))/2;
    for (unsigned long long x=0;x<ne;x++)
      s+=D.data[x];
    */
    tfmed.EndClock(true);
    cout << "s=" << s << endl;
    cout << "x=" << x << endl;
    return 0;
}

#endif
