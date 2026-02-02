#!/bin/bash

# Determine the script's directory and warpy base directory
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
WARPY_BASE="$( cd "$SCRIPT_DIR/.." && pwd )"

# Detect environment
if [ -n "$BPIPE_DEFAULT_ENVIRONMENT" ]; then
    ENV="$BPIPE_DEFAULT_ENVIRONMENT"
elif [[ "$OSTYPE" == "darwin"* ]]; then
    ENV="mac_dev"
elif command -v sbatch &> /dev/null; then
    ENV="vcgs_mk4"
else
    ENV="ont_server"
fi

function usage() {
      echo "Usage: $0 [options]"
      echo
      echo "-s <sample id> :   sample id"
      echo "-d <dir>       :   directory containing pod5 files" 
      echo ""
      exit 0
}

while getopts "s:d:h" arg; do
  case $arg in
    s) # sample id
      export sample_id=$OPTARG
      ;;
    d) # data
      # Convert relative path to absolute
      if [[ "$OPTARG" = /* ]]; then
        export data_path=$OPTARG
      else
        export data_path="$(cd "$(dirname "$OPTARG")" && pwd)/$(basename "$OPTARG")"
      fi
      ;;
     h | *) # Display help.
         usage
      ;;
    :)
      echo "Parameter -$OPTARG requires a positional argument"
      exit 1
      ;;
  esac
done

if [ -z "$sample_id" ] || [ -z "$data_path" ];
then
    echo "Please specify sample and data path"

    usage
fi

# Set batch directory at warpy base level
BATCH_DIR=$WARPY_BASE/batches/${sample_id}

# Create batches directory if it doesn't exist
mkdir -p "$WARPY_BASE/batches"

mkdir $BATCH_DIR

echo '
samples:
  - identifier: $sample_id
    familyId: $sample_id
    sex: female
    inputs:
        - $data_path
    parents: null

' | envsubst | tee $BATCH_DIR/samples.yaml


echo
echo "Done."
echo
echo "
To run this:

cd $BATCH_DIR

# TODO: ADJUST SEX OF SAMPLES in samples.yaml

bpipe run -p ANALYSIS_QUAL=hac --env $ENV $WARPY_BASE/src/pipeline.groovy -targets $WARPY_BASE/designs/WGS_REFFLAT_10_NOALT_38/WGS_REFFLAT_10_NOALT_38.bed -samples samples.yaml
"
