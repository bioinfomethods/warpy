#!/usr/bin/env bash

set -eo pipefail

module load htslib

zcat "$1" | awk 'BEGIN{OFS="\t"} $1 ~ /^chr([1-9]|1[0-9]|2[0-2]|X|Y|M)$/' | bgzip -c > "$2"

tabix -p bed "$2"
