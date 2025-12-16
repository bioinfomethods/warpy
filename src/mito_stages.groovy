register_mito_bam = {
    mito_sample_bams.get(sample, []).add(input.toString())
}

forward_mito_sample_bam = {
    forward mito_sample_bams[sample]
}

run_clairs_to = {
    var bam_ext : 'bam'

    def sample_output_dir = "mito/${sample}"

    output.dir = "${sample_output_dir}/variants"

    transform("(.*).$bam_ext") to('$1.mito.dp.vcf') {
        def output_prefix = "${output.prefix.prefix}".replaceFirst("^${sample_output_dir}", "")

        exec """
            mkdir -p ${output.dir.prefix}

            /opt/bin/run_clairs_to \
                --tumor_bam_fn ${input[bam_ext]} \
                --ref_fn /misc/bioinf-ops/jimmy.chiu/garvan_LRS_chrM_coverage/chrM_ref/Homo_sapiens_assembly38.no_alt.chrM_no_tp_circular.fasta \
                --sample_name $sample \
                --ctg_name chrM \
                --threads $threads \
                --platform ${clairs_to_model.clairs_to_model_name} \
                --output_dir $sample_output_dir \
                --snv_output_prefix "${output_prefix}.snv" \
                --indel_output_prefix "${output_prefix}.indel" \
                --conda_prefix /opt/micromamba/envs/clairs-to

            bcftools concat \
                "${output.prefix.prefix}.snv.vcf.gz" \
                "${output.prefix.prefix}.indel.vcf.gz" | \
            bcftools sort \
                -Oz \
                -T "$output.dir/XXXXXX" \
                -o ${output.prefix.prefix}.vcf.gz

            tabix -p vcf ${output.prefix.prefix}.vcf.gz

            rm -rf "${output.prefix.prefix}.snv.vcf.gz*" \
                "${output.prefix.prefix}.indel.vcf.gz*"

            bcftools view \
                -H ${output.prefix.prefix}.vcf.gz | \
            awk '!/^#/ { print \$1 "\t" (\$2 - 1) "\t" \$2 }' > ${output.prefix.prefix}.variants.bed

            samtools depth \
                -b ${output.prefix.prefix}.variants.bed \
                $input.bam | \
            bgzip -c > ${output.prefix.prefix}.variants.dp.tab.gz

            tabix -s 1 -b 2 -e 2 ${output.prefix.prefix}.variants.dp.tab.gz

            bcftools annotate \
                -a ${output.prefix.prefix}.variants.dp.tab.gz \
                -c "CHROM,POS,INFO/DP" \
                -H "##INFO=<ID=DP,Number=1,Type=Integer,Description=\\"Site read depth; no filtering applied\\">" \
                ${output.prefix.prefix}.vcf.gz > $output.mito.dp.vcf
        """
    }
}

annotate_mito_variants = {
    output.dir = "mito/${sample}/variants"

    transform('.mito.dp.vcf') to('.mito.dp.vep.vcf.gz') {
        exec """
            vep --vcf_info_field ANN \
                --assembly GRCh38 \
                --protein \
                --pubmed \
                --hgvs \
                --no_stats \
                --dir $VEP_CACHE \
                --offline \
                --vcf \
                -i $input.vcf \
                --everything \
                -o STDOUT > \
                $output.prefix

            bgzip -c $output.prefix > $output.vcf.gz

            tabix -p vcf $output.vcf.gz

            rm $output.prefix
        """

        mito_sample_vcfs[sample] = output.toString()
    }
}

run_mitoreport = {
    var maternal_ids : [:]

    def proband_mito_vcf = mito_sample_vcfs[sample]

    output.dir = "mito/mitoreports/${sample}"

    produce('index.html') {
        if (sample in maternal_ids) {
            def maternal_sample = maternal_ids[sample]
            def maternal_mito_vcf = mito_sample_vcfs[maternal_sample]

            exec """
                $tools.JAVA -jar $MITOREPORT_JAR \
                    mito-report \
                    --sample $sample \
                    -mann $MITOMAP_ANNOTATIONS \
                    -gnomad $GNOMAD_MITO_REF_VCF \
                    -vcf $proband_mito_vcf \
                    --maternal-vcf $maternal_mito_vcf \
                    --maternal-sample $maternal_sample \
                    --output-dir $output.dir \
                    $input.bam \
                    "${MITO_CONTROL_DIR}/control_1.bam" \
                    "${MITO_CONTROL_DIR}/control_2.bam" \
                    "${MITO_CONTROL_DIR}/control_3.bam"
            """
        }
        else {
            exec """
                $tools.JAVA -jar $MITOREPORT_JAR \
                    mito-report \
                    --sample $sample \
                    -mann $MITOMAP_ANNOTATIONS \
                    -gnomad $GNOMAD_MITO_REF_VCF \
                    -vcf $proband_mito_vcf \
                    --output-dir $output.dir \
                    $input.bam \
                    "${MITO_CONTROL_DIR}/control_1.bam" \
                    "${MITO_CONTROL_DIR}/control_2.bam" \
                    "${MITO_CONTROL_DIR}/control_3.bam"
            """
        }
    }
}
