somalier_extract = {
    doc "Extract sites for use by somalier relate"

    var bam_ext : 'bam'

    output.dir="qc/somalier"

    produce("${sample}.somalier") {
        exec """
            set -eo pipefail

            export SOMALIER_SAMPLE_NAME=${sample}

            $tools.SOMALIER extract
                -d $output.dir
                -f $REF
                -s $SOMALIER_SITES
                ${input[bam_ext]}
        """
    }

    sample_somaliers.add(output.toString())
}

somalier_relate = {
    doc "Evaluate the relatedness between all samples in the batch from extracted sites and ped file"

    var pedigree_path : false

    def pedigree_source = pedigree_path?:input.ped

    output.dir="qc"

    println "Inputs are: " + sample_somaliers

    from([pedigree_source, *sample_somaliers]) produce("somalier.html") {
        exec """
            $tools.SOMALIER relate
                --ped $input.ped
                --output-prefix $output.prefix
                $inputs.somalier
        """
    }
}

check_variant_fraction = {
    var MAX_VAF_OUTLIER_FRAC : 1.0

    output.dir = "qc"

    produce(sample + '.variant_fraction.txt') {
        // Note: we set the fraction to 1.0 below because the actual check
        // of the correct fraction is done below now
        exec """
             set -eo pipefail

             unset GROOVY_HOME

             JAVA_OPTS="-Xmx3g" $tools.GROOVY -cp $tools.GNGS_JAR $BASE/scripts/qc/CalculateVAFOutliers.groovy
                -q
                -s $sample
                -f 1.0
                $input.vcf.gz
             | tee $output.txt
        """
    }

    def vaf_outliers = 1.0
    try { vaf_outliers = file(output.txt).text.find(/0\.\d{1,10}/) } catch(Exception e) { println "WARNING: could not read vaf outlier fraction from $output.txt" }

    check('variant fraction for sample ' + sample) {
        groovy """
            System.exit($vaf_outliers > $MAX_VAF_OUTLIER_FRAC ? 1 : 0)
        """
    } otherwise {
        send text {"Sample $sample had excessive fraction of VAF outliers ($vaf_outliers) outside expected range "} to channel: 'cpipe_operator'
    }
}

calc_sample_read_len_dist = {
    output.dir = "qc/read_length_dist"

    def sample_bams = sample_processed_bams.values().flatten().collect { "${it}" }.join(" ")

    produce("sample_read_length_dist.txt") {
        exec """
            set -eo pipefail

            echo -e \"sample_id\\tlength\\tcount\" > $output.txt

            for bam in ${sample_bams}

            do
                fname=\$(basename $bam)

                sample_id="\${fname%%.*}\"

                $tools.SAMTOOLS stats -@ $threads $bam \
                    | awk -v id=\"$sample_id\" -F'\\t' '\$1==\"RL\" {print id\"\\t\"\$2\"\\t\"\$3}' \
                >> $output.txt

            done
        """
    }
}

plot_read_length_dist = {
    output.dir = "qc/read_length_dist"

    def target_chrs = targets_by_chr*.chr.join(' ')
    def sample_bam = sample_processed_bams.values().flatten().first()

    transform(".txt") to(".pdf", ".summary.txt") {
        exec """
            set -eo pipefail

            total_len=\$(awk 'NR==FNR{chroms[\$1]; next} \$1 in chroms{sum+=\$2} END{print sum}' \
                <(echo "${target_chrs}" | tr ' ' '\\n') \
                <(samtools view -H $sample_bam | \
                    awk -F'\t' '\$1=="@SQ"{split(\$2,chr,":"); split(\$3,len,":"); print chr[2]"\t"len[2]}'))

            $tools.RSCRIPT $BASE/scripts/plot_sample_coverage_vs_read_len_dist.R \
                $input.txt \
                $output.pdf \
                $output.summary.txt \
                \$total_len
        """
    }
}
