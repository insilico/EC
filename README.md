Evaporative Cooling (EC)
========================

#### A feature selection tool for GWAS and other biological data ####

### Description ###
EC is a free, open-source command-line tool for analysis of GWAS (SNP) and
other types of biological data.  Several modes are available for various types
of analysis, including:

Evaporative Cooling (EC) is C++ library that provides a flexible feature 
selection algorithm for SNPs and quantitative data, using ReliefF and Random 
Jungle for interactions and main effects, respectively.  EC is also available 
as a standalone [tool](http://insilico.utulsa.edu/evaporative-cooling).

EC is being developed by the In Silico Research Group at the Tandy School
of Computer Science of the [University of Tulsa](http://www.utulsa.edu).  Our
research is sponsored by the NIH and William K. Warren foundation.  For more
details, visit our research [website](http://insilico.utulsa.edu).

### Dependencies ###
* EC library, available as a source release on the EC project page, and its
dependencies:

  * [Random Jungle](http://github.com/insilico/randomjungle)

  * gfortran, sometimes installed alongside compiler tools

  * GNU Scientific library (libgsl)

  * libxml2

* [Boost](http://www.boost.org) system, filesystem, and program-options libraries 
 
* The libz/zlib compression library is required, but this is installed by default
on most Unix systems.  In MinGW libz is installed via mingw-get.

* OpenMP is required to take advantage of the parallelized tree growth in 
Random Jungle and distance matrix calculations for ReliefF.  This is another 
library typically installed alongside the compiler toolchain.

### Compilation Environment and Instructions ###
To compile this code, a GNU toolchain and suitable environment are required.
GNU g++ has been used to successfully compile the code.

We have successfully built and run EC on:

 * Linux (64-bit Ubuntu) (gcc-4.6)
 * Mac (10.6 - 10.7) (gcc-4.2.1)
 * Windows 7 (32-bit) using the [MinGW](http://www.mingw.org) compiler system
  (gcc-4.6)

To build EC, first run the bootstrap script

    ./bootstrap.sh

Ignore any extraneous warnings. This calls autoreconf and generates the 
configure script.  From this point, a standard

    ./configure && make && sudo make install

will generate the `Makefile`, compile and link the code, and copy the objects to
the installation directory (default of `/usr/local`).  As is convention, headers
are installed in `$PREFIX/include`, binary in `$PREFIX/bin`, and the library in
`$PREFIX/lib`.

The resulting binary src/ec_static.exe will run as a command-line tool.

### Usage ###

Allowed options:
  --help                                produce help message
  --verbose                             verbose output
  --convert                             convert data set to data set - no ec
  -T [ --optimize-temp ]                optimize coupling constant T
  -c [ --config-file ] arg              read configuration options from file - 
                                        command line overrides these
  -s [ --snp-data ] arg                 read SNP attributes from genotype 
                                        filename: txt, ARFF, plink (map/ped, 
                                        binary, raw)
  --snp-file-type arg                   Ignore file extension and use type: 
                                        textwhitesp, wekaarff, plinkped, 
                                        plinkbed, plinkraw, mayogeo, birdseed
  -n [ --numeric-data ] arg             read continuous attributes from 
                                        PLINK-style covar file
  -X [ --numeric-transform ] arg        perform numeric transformation: 
                                        normalize, standardize, zscore, log, 
                                        sqrt
  -a [ --alternate-pheno-file ] arg     specifies an alternative 
                                        phenotype/class label file; one value 
                                        per line
  -g [ --ec-algorithm-steps ] arg (=all)
                                        EC steps to run (all|rj|rf)
  -t [ --ec-num-target ] arg (=0)       EC N_target - target number of 
                                        attributes to keep
  -r [ --ec-iter-remove-n ] arg (=0)    Evaporative Cooling number of 
                                        attributes to remove per iteration
  -p [ --ec-iter-remove-percent ] arg   Evaporative Cooling precentage of 
                                        attributes to remove per iteration
  -O [ --out-dataset-filename ] arg     write a new tab-delimited data set with
                                        EC filtered attributes
  -o [ --out-files-prefix ] arg (=ec_run)
                                        use prefix for all output files
  --snp-metric arg (=gm)                metric for determining the difference 
                                        between subjects (gm|am|nca|nca6)
  -B [ --snp-metric-nn ] arg (=gm)      metric for determining the difference 
                                        between subjects (gm|am|nca|nca6|km)
  -W [ --snp-metric-weights ] arg (=gm) metric for determining the difference 
                                        between SNPs (gm|am|nca|nca6)
  -N [ --numeric-metric ] arg (=manhattan)
                                        metric for determining the difference 
                                        between numeric attributes 
                                        (manhattan=|euclidean)
  -R [ --rj-run-mode ] arg (=1)         Random Jungle run mode: 1 
                                        (default=library call) / 2 (system 
                                        call)
  -j [ --rj-num-trees ] arg (=1000)     Random Jungle number of trees to grow
  --rj-mtry arg (=0)                    Random Jungle size of randomly chosen 
                                        variable sets, DEFAULT: sqrt(ncol)
  --rj-nimpvar arg (=1)                 Random Jungle only necessary if 
                                        backsel>0. SIZE=[1-...] how many 
                                        variable should remain
  --rj-impmeasure arg (=1)              Random Jungle importance method (see RJ
                                        docs)
  --rj-backsel arg (=0)                 Random Jungle backward elimination (see
                                        RJ docs)
  -Y [ --rj-tree-type ] arg (=1)        Random Jungle tree type: 1 (default)-5 
                                        (see RJ docs)
  -M [ --rj-memory-mode ] arg (=0)      Random Jungle memory mode: 0 
                                        (default=double) / 1 (float) / 2 (char)
  -x [ --snp-exclusion-file ] arg       file of SNP names to be excluded
  -k [ --k-nearest-neighbors ] arg (=10)
                                        set k nearest neighbors
  -m [ --number-random-samples ] arg (=0)
                                        number of random samples (0=all|1 <= n 
                                        <= number of samples)
  -b [ --weight-by-distance-method ] arg (=equal)
                                        weight-by-distance method 
                                        (equal|one_over_k|exponential)
  --weight-by-distance-sigma arg (=2)   weight by distance sigma
  -d [ --diagnostic-tests ] arg         performs diagnostic tests and sends 
                                        output to filename without running EC
  -D [ --diagnostic-levels-file ] arg   write diagnostic attribute level counts
                                        to filename
  --dge-counts-data arg                 read digital gene expression counts 
                                        from text file
  --dge-norm-factors arg                read digital gene expression 
                                        normalization factors from text file
  --birdseed-snps-data arg              read SNP data from a birdseed formatted
                                        file
  --birdseed-phenos-data arg            read birdseed subjects phenotypes from 
                                        a text file
  --birdseed-subjects-labels arg        read subject labels from filename to 
                                        override names from data file
  --birdseed-include-snps arg           include the SNP IDs listed in the text 
                                        file
  --birdseed-exclude-snps arg           exclude the SNP IDs listed the text 
                                        file
  --distance-matrix arg                 create a distance matrix for the loaded
                                        samples and exit
  --gain-matrix arg                     create a GAIN matrix for the loaded 
                                        samples and exit
  --dump-titv-file arg                  file for dumping SNP 
                                        transition/transversion information

All commands will include an input file (`-s/--snp-data`), and, optionally, 
an output file prefix (`-o/--output-files-prefix`).

To perform a standard, all-default-parameters analysis,

    ./ec_static -s snpdata.ped -o result

This will use genotype/phenotype information from `snpdata.ped`, a PLINK
plaintext GWAS file, in the feature selection.  All of the output files 
produced will be prepended with 'result'.

This produces a file called `result.ec`, in which the SNPs are ranked 
in descending order.

For additional examples, see the [EC](http://insilico.utulsa.edu/evaporative-cooling)
page on our research website.

### Contributors ###
See [AUTHORS](https://github.com/insilico/EC/blob/master/AUTHORS) file.

### References ###
B.A. McKinney, J.E.  Crowe, Jr.,  J. Guo, and D. Tian,  ÒCapturing the
 spectrum of interaction effects in genetic  association  studies  by
simulated evaporative cooling network  analysis,Ó  PLoS Genetics.
5(3):  e1000432. doi:10.1371/journal.pgen.1000432; 2009.

McKinney, B.A., Reif, D.M., White, B.C., Crowe, J.E., Moore, J.H. Evaporative 
cooling feature selection for genotypic data involving interactions. 
Bioinformatics 23, 2113-2120 (2007). [PubMed]
