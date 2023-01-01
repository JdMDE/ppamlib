parallelpamlib: a library to implement the Partitioning Around Medoids (PAM)
algorithm in parallel.
It uses the data in jmatrix format, a specific library of matrix manipulation
that allows extremely big matrices (as long as the RAM of the machine allows).
Apart from the PAM itself the library also implements in parallel the calculation
of the distance/dissimilarity matrix (metrics L1 and L2 and Pearson dissimilarity)
and the silhouette of the resulting clustering.
It includes four test programs: pardis, parpam, parsil and tdvalue.
