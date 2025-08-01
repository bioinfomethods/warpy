import groovy.json.JsonBuilder

call_str = {
    var bam_ext: 'bam'

    branch.dir = "str/$sample/${chr}"

    produce("${sample}.${chr}.straglr.tsv", "${sample}.${chr}.straglr.vcf.gz") {
        def sampleSex = meta.get(sample).sex

        exec """

            { grep $chr -Fw $calling.repeats_bed || true; } > $output.dir/repeats_subset.bed

            if [[ -s $output.dir/repeats_subset.bed ]]; then
                straglr-genotype --loci $output.dir/repeats_subset.bed
                    --sample ${sample}
                    --tsv $output.tsv
                    -v $output.vcf.gz.prefix
                    --sex ${sampleSex.toLowerCase()} ${input[bam_ext]} $REF
                    --min_support 1
                    --min_cluster_size 1

                set -o pipefail

                cat $output.vcf.gz.prefix | vcfstreamsort | bgziptabix $output.vcf.gz;
            else
                echo "blank subset BED for $chr";
            fi
        """
    }
}

annotate_repeat_expansions = {

    exec """
        set -o pipefail

        stranger -f ${calling.variant_catalogue_hg38} $input.vcf.gz
            | sed 's/\\ /_/g'
            | bgziptabix $output.vcf.gz

        SnpSift extractFields $output.vcf.gz CHROM POS ALT FILTER REF RL RU REPID VARID STR_STATUS > $output.stranger.tsv

        SnpSift extractFields $output.vcf.gz
            CHROM POS DisplayRU STR_NORMAL_MAX STR_PATHOLOGIC_MIN VARID Disease > $output.plot.tsv
    """
}

merge_str_tsv = {

    output.dir = 'str'

    produce("${sample}.plot.tsv", "${sample}.stranger.tsv", "${sample}.straglr.tsv") {
        exec """
            awk 'NR == 1 || FNR > 1' $input.plot.tsv > $output.plot.tsv

            awk 'NR == 1 || FNR > 1' $input.stranger.tsv > $output.stranger.tsv

            sed '1d' $input.straglr.tsv | awk 'NR == 1 || FNR > 1' > $output.straglr.tsv
        """
    }
}

merge_str_vcf = {

    output.dir = 'str'

    produce("${sample}.wf_str.vcf.gz") {
        exec """
            bcftools concat $inputs.vcf.gz > $output.vcf.gz.prefix

            bgzip -f $output.vcf.gz.prefix

            tabix $output.vcf.gz
        """
    }
}

// This part was too difficult to port: it depends on custom scripts 
// specific to the epi2me labs environment
// would like to review possibility of integrating <a href='https://stripy.org'>STRipy</a>.
//
//make_str_report = {
//   
//    output.dir = "str"
//
//    new File("str/params.json").text = new JsonBuilder(calling).toPrettyString()
//
//    produce("${sample}.wf-human-str-report.html") {
//        exec """
//            $BASE/scripts/workflow-glue report_str
//                -o $output.html
//                --params str/params.json
//                --sample_name ${sample}
//                --version versions.txt
//                --vcf $input.vcf.gz
//                --straglr $input.straglr.tsv
//                --stranger $input.plot.tsv
//                --stranger_annotation $input.stranger.tsv
//                --read_stats $input.readstats.tsv.gz
//        """
//    }
//}
//


