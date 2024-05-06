bam2bedmethyl = {
    def modkit_args = calling.modkit_args ?: '' 
    
    def regionFlag = region.isEmpty() ? "" : "--region $region"
    
    output.dir = "methyl"

    produce("${opts.sample}.methyl.cpg.bed.gz", "${opts.sample}.methyl.cpg.log.txt") {
        exec """
            modkit pileup
                --preset traditional
                --ref $REF
                --only-tabs
                --log-filepath $output.txt
                --threads $threads $regionFlag
                --prefix $output.dir/${opts.sample}.methyl ${modkit_args} $input.bam -
                | bgzip -c > $output.bed.gz
        """
    }
}
