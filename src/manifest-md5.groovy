#!/hpc/bpipeLibrary/shared/cpipe-wgs/bpipe

import java.time.LocalDate
import java.time.format.DateTimeFormatter

TMP_DIR="/ibm/hpcfs1/tmp"
HPC_SHARED="/hpc/bpipeLibrary/shared"

BASE="$HPC_SHARED/cpipe-wgs-test"
ARCHIE_HOME="$HPC_SHARED/archie-test"

println "************************************************************"
println "This script creates a manifest.md5 file for the input files."
println "************************************************************"

pipeline_dir = new File('.').canonicalPath
if (pipeline_dir.contains('-test') || pipeline_dir.contains('-dev')) {
    println "Detected test environment"

    BASE="$HPC_SHARED/cpipe-wgs"
    ARCHIE_HOME="$HPC_SHARED/archie-test"
}


md5_manifest = {
    produce("manifest.md5") {
        exec """
            source $ARCHIE_HOME/archie-cli/scripts/activate-archie.sh

            python $BASE/pipeline/scripts/md5_progress.py $inputs > manifest.md5
        """
    }
}


run {
    md5_manifest
}
