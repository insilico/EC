#!/bin/bash

ec --snp-data ~/analysis/bd/after_review2/chrX_and_rs17480050_removed.no-LD.bed \
--ec-iter-remove-percent 50 \
--ec-num-target 1000 \
--verbose 1 \
--out-files-prefix wtccc-all-benchmark_mac

