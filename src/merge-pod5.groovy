#!/hpc/bpipeLibrary/shared/cpipe-wgs/bpipe

import java.time.LocalDate
import java.time.format.DateTimeFormatter

TMP_DIR="/ibm/hpcfs1/tmp"
HPC_SHARED="/hpc/bpipeLibrary/shared"

BASE="$HPC_SHARED/cpipe-wgs-test"
ARCHIE_HOME="$HPC_SHARED/archie-test"
DEFAULT_GROUP_COUNT=8

options {
    sample 'Sample identifier', args: 1, required: true
    flowcell 'Flowcell identifier', args: 1, required: true
    groups 'Number of merged pod5 files to procude', args: 1, required: false
}

optGroups = opts.groups
groupCount = optGroups && optGroups.toInteger() > 0 ? opts.groups.toInteger() : DEFAULT_GROUP_COUNT

println "*********************************************************************************"
println "This script merges all given ${args.size()} pod5 files into ${groupCount} groups."
println "*********************************************************************************"

pipeline_dir = new File('.').canonicalPath
if (pipeline_dir.contains('-test') || pipeline_dir.contains('-dev')) {
    println "Detected test environment"

    BASE="$HPC_SHARED/cpipe-wgs"
    ARCHIE_HOME="$HPC_SHARED/archie-test"
}


def argsPerGroup = args.size().intdiv(groupCount)
def groups = args.collate(argsPerGroup)
  .collect()
  .withIndex()
  .collectEntries {
    groupArgs, groupNum ->
    [(groupNum):groupArgs]
  }


set_group_info = {
    branch.groupNum = branch.name
    branch.sample = "${opts.sample}"
    branch.flowcell = "${opts.flowcell}"
    forward groups[branch.groupNum]
}


merge_pod5s = {
    output.dir = "merged"

    produce("${sample}.${flowcell}.${groupNum}.pod5") {
        println "Merging pod5 files: $inputs"

        exec """
            pod5 merge $inputs -o $output.pod5 --force-overwrite
        """, "pod5_file_format_sl"
    }
}


run {
    groups * [ set_group_info + merge_pod5s ]
}
