#!/usr/bin/bash

set -eo pipefail

module load R

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

Rscript $SCRIPT_DIR/get_vcf_qc_dist.R $1 $2 $3
