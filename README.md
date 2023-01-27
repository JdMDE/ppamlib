ppamlib: a library to implement the Partitioning Around Medoids (PAM)
algorithm in parallel.
It uses the data in jmatrix format, a specific library of matrix manipulation
that allows extremely big matrices (as long as the RAM of the machine allows).
Apart from the PAM itself the library also implements in parallel the calculation
of the distance/dissimilarity matrix (metrics L1 and L2 and Pearson dissimilarity)
and the silhouette of the resulting clustering.

It includes four test programs:

pardis: Parallel calculation of distance/dissimilarity matrix from a jmatrix with data.

parpam: Parallel implementation of the Partitioning Around Medoids (PAM) algorithm from a distance matrix.

parsil: Parallel calculation of the silhouette of each points after the clustering has been applied.

tdvalue: Calculation of the value of the optimization function of the PAM algorithm for a given clusterization result.

These library uses the library jmatlib (see https://github.com/JdMDE/jmatlib) which therefore needs to be
installed before compilation and use of ppamlib.

The code of this library with interface modifications is also used inside the parallelpam R package
(https://CRAN.R-project.org/package=parallelpam) and inside the scellpam package (https://CRAN.R-project.org/package=scellpam)

See documentation at subdirectory html, file index.html
