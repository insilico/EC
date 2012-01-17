#!/bin/bash

time ec --snp-data ~/analysis/bd/after_review2/chrX_and_rs17480050_removed_no-LD.bed \
--ec-iter-remove-percent 50 \
--ec-num-target 1000 \
--verbose 1 \
--out-files-prefix wtccc-all-benchmark_thor_matrix_d30 > wtccc-all-benchmark_thor_matrix_s4.stdout

