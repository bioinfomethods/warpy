bam2bedmethyl = {
    var bam_ext = 'bam'

    def modkit_args = calling.modkit_args ?: '' 
    
    def regionFlag = region.isEmpty() ? "" : "--region $region"
    
    output.dir = "methyl"

    produce("${sample}.methyl.cpg.bed.gz", "${sample}.methyl.cpg.log.txt") {
        exec """
            modkit pileup
                --preset traditional
                --ref $REF
                --only-tabs
                --log-filepath $output.txt
                --threads $threads $regionFlag
                --prefix $output.dir/${sample}.methyl ${modkit_args} ${input[bam_ext]} -
                | bgzip -c > $output.bed.gz
        """
    }
}

modbam2bed  = {

//    label "wf_human_methyl"
//    cpus params.threads
//    input:
//        tuple path(alignment), path(alignment_index), val(alignment_meta)
//        tuple path(reference), path(reference_index), path(reference_cache)
//    output:
//        path "${params.sample_name}.methyl.*", emit: methyl_outputs

//    script:
//    def ref_path = "${reference_cache}/%2s/%2s/%s:" + System.getenv("REF_PATH")

    def modbam2bed_args = calling.modbam2bed_args ?: ''
    
    def regionFlag = region.isEmpty() ? "" : "-r $region"
    
    output.dir = "methyl"

    produce("${sample}.methyl.cpg.bed.gz", "${sample}.methyl.cpg.acc.bed.gz") {
        exec """
            modbam2bed
                -e -m 5mC --cpg -t $threads $regionFlag
                $REF
                $input.bam
                --aggregate
                --prefix $output.dir/${sample}.methyl ${modbam2bed_args}
                | bgzip -c > $output.bed.gz

            bgzip -f $output2.bed.gz.prefix
        """
    }
}
