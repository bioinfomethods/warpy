#!/usr/bin/bash

set -eo pipefail

module load R

Rscript get_vcf_qc_dist.R $1 $2
