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

Ignore any extranneous warnings. This calls autoreconf and generates the 
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
  -c [ --config-file ] arg              read configuration options from file - 
                                        command line overrides these
  -s [ --snp-data ] arg                 read SNP attributes from genotype 
                                        filename: txt, ARFF, plink (map/ped, 
                                        binary, raw)
  -n [ --numeric-data ] arg             read continuous attributes from 
                                        PLINK-style covar file
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
  -S [ --snp-metric ] arg (=gm)         metric for determining the difference 
                                        between SNPs (gm|am|nca)
  -N [ --numeric-metric ] arg (=manhattan)
                                        metric for determining the difference 
                                        between numeric attributes 
                                        (manhattan=|euclidean)
  -j [ --rj-num-trees ] arg (=1000)     Random Jungle number of trees to grow
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
  --dge-phenos-data arg                 read digital gene expression phenotypes
                                        from text file
  --dge-norm-factors arg                read digital gene expression 
                                        normalization factors from text file
  --birdseed-snps-data arg              read SNP data from a birdseed formatted
                                        file
  --birdseed-phenos-data arg            read birdseed subjects phenotypes from 
                                        text file
  --birdseed-subjects-labels arg        read subject labels from filename to 
                                        override names from data file
  --birdseed-include-snps arg           read data SNPs data only for the 
                                        subject IDs in file
  --birdseed-exclude-snps arg           read data SNPs data only for the 
                                        subject IDs in file
  --distance-matrix arg                 create a distance matrix for the loaded
                                        samples and exit
  --gain-matrix arg                     create a GAIN matrix for the loaded 
                                        samples and exit

All commands will include an input file (`-s/--snp-data`), and, optionally, 
an output file prefix (`-o/--output-files-prefix`).

To perform a standard, all-default-parameters analysis,

    ./ec_static -s snpdata.ped -o result

This will use genotype/phenotype information from `snpdata.ped`, a PLINK
plaintext GWAS file, in the feature selection.  All of the output files 
produced will be prepended with 'result'.

This produces a file called `ec_run.ec`, in which the SNPs are ranked 
in descending order.

For additional examples, see the [EC](http://insilico.utulsa.edu/evaporative-cooling)
page on our research website.

### Contributors ###
See AUTHORS file.

### References ###
Moore, J.H., White, B.C. Tuning ReliefF for genome-wide genetic analysis. 
Lecture Notes in Computer Science 4447, 166-175 (2007). [Springer]

McKinney, B.A., Reif, D.M., White, B.C., Crowe, J.E., Moore, J.H. Evaporative 
cooling feature selection for genotypic data involving interactions. 
Bioinformatics 23, 2113-2120 (2007). [PubMed]
