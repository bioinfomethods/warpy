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
      echo "Single-sample mode:"
      echo "  -s <sample id> :   sample id"
      echo "  -d <dir>       :   directory containing pod5 files"
      echo ""
      echo "Multi-sample mode:"
      echo "  -D <parent dir> :  parent directory where each subdirectory is a sample"
      echo ""
      echo "Common options:"
      echo "  -b <batch name> :  batch name (default: sample_id or parent dir basename)"
      echo ""
      exit 0
}

while getopts "s:d:D:b:h" arg; do
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
    D) # parent directory (multi-sample mode)
      if [[ "$OPTARG" = /* ]]; then
        export parent_dir=$OPTARG
      else
        export parent_dir="$(cd "$(dirname "$OPTARG")" && pwd)/$(basename "$OPTARG")"
      fi
      ;;
    b) # batch name override
      export batch_name=$OPTARG
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

# Validate arguments: must use either single-sample mode or multi-sample mode
if [ -n "$parent_dir" ] && ([ -n "$sample_id" ] || [ -n "$data_path" ]); then
    echo "Error: Cannot combine -D with -s/-d. Use one mode or the other."
    usage
fi

if [ -z "$parent_dir" ] && ([ -z "$sample_id" ] || [ -z "$data_path" ]); then
    echo "Please specify either -s and -d (single sample) or -D (multi-sample)"
    usage
fi

if [ -n "$parent_dir" ] && [ ! -d "$parent_dir" ]; then
    echo "Error: Parent directory does not exist: $parent_dir"
    exit 1
fi

# Determine batch name
if [ -z "$batch_name" ]; then
    if [ -n "$parent_dir" ]; then
        batch_name=$(basename "$parent_dir")
    else
        batch_name=$sample_id
    fi
fi

# Set batch directory at warpy base level
BATCH_DIR=$WARPY_BASE/batches/${batch_name}

# Create batches directory if it doesn't exist
mkdir -p "$WARPY_BASE/batches"

mkdir $BATCH_DIR

if [ -n "$parent_dir" ]; then
    # Multi-sample mode: each subdirectory is a sample
    echo "samples:" > $BATCH_DIR/samples.yaml
    found_samples=0
    for subdir in "$parent_dir"/*/; do
        [ -d "$subdir" ] || continue
        sub_sample_id=$(basename "$subdir")
        # Convert to absolute path (remove trailing slash)
        sub_data_path="${subdir%/}"
        cat >> $BATCH_DIR/samples.yaml <<EOF
  - identifier: ${sub_sample_id}
    familyId: ${batch_name}
    sex: female
    inputs:
        - ${sub_data_path}
    parents: null

EOF
        found_samples=$((found_samples + 1))
    done
    if [ "$found_samples" -eq 0 ]; then
        echo "Error: No subdirectories found in $parent_dir"
        rm -rf "$BATCH_DIR"
        exit 1
    fi
    echo "Generated samples.yaml with $found_samples sample(s)"
    cat $BATCH_DIR/samples.yaml
else
    # Single-sample mode (original behavior)
    echo '
samples:
  - identifier: $sample_id
    familyId: $sample_id
    sex: female
    inputs:
        - $data_path
    parents: null

' | envsubst | tee $BATCH_DIR/samples.yaml
fi


echo
echo "Done."
echo
echo "
To run this:

cd $BATCH_DIR

# TODO: ADJUST SEX OF SAMPLES in samples.yaml

bpipe run -p ANALYSIS_QUAL=hac --env $ENV $WARPY_BASE/src/pipeline.groovy -targets $WARPY_BASE/designs/WGS_REFFLAT_10_NOALT_38/WGS_REFFLAT_10_NOALT_38.bed -samples samples.yaml
"
