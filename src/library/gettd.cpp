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

#include "../headers/gettd.h"

//' GetTD
//'
//' Function that takes a PAM classification (as returned by ApplyPAM) and the dissimilarity matrix and returns the value of the TD function
//' (sum of dissimilarities between each point and its closest medoid, divided by the number of points).
//' This function is mainly for debugging/internal use.
//'
//' @param Lmed         The vector Lmed as returned by ApplyPAM (please, consult the help of ApplyPAM for details)
//' @param Lclasif      The vector Lclasif as returned by ApplyPAM (please, consult the help of ApplyPAM for details)
//' @param D            A reference to a symmetric matrix which is the distance/dissimilarity matrix.
//' @return TD          The value of the TD function.
template <typename disttype>
double GetTD(std::vector<indextype> Lmed,std::vector<indextype> Lclasif,SymmetricMatrix<disttype> &D)
{
 double TD=0.0;
 for (indextype k=0;k<Lclasif.size();k++)
  TD += double(D.Get(k,Lmed[Lclasif[k]]));

 return(TD/double(Lclasif.size()));
}

template double GetTD<float>(std::vector<indextype> Lmed,std::vector<indextype> Lclasif,SymmetricMatrix<float> &D);
template double GetTD<double>(std::vector<indextype> Lmed,std::vector<indextype> Lclasif,SymmetricMatrix<double> &D);
