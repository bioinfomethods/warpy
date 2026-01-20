

ext = { File file ->
    file.name.tokenize('.')[-1]
}

scanInputDirectories = { args ->

if(args.size()==0) 
    throw new bpipe.PipelineError(
        """
        No data directories were provided to analyse. 

        Please provide one or more data directories containing pod5, fast5 or blow5 files as arguments.
        """
    )
    
    
nonDirectories = args.grep { ! new File(it).directory }    
if(nonDirectories) {
       throw new bpipe.PipelineError(
        """
        One or more arguments given was not a valid directory. Please pass directories containing files for analysis as arguments.

        Non-directory arguments: ${nonDirectories.join(", ")}
        """
    )
 
}

List input_files = args.collect { println(it); new File(it) }*.listFiles().flatten().grep { println(it); ext(it) in ['fast5','blow5','pod5','bam']  }

}

check_tools = {
    assert tools.POD5 != null : "pod5 utility is not configured"
    
    
}

convert_fast5_to_pod5 = {
    output.dir = 'pod5'

    transform('*.fast5') to('.pod5') {
        exec """
            set -o pipefail

            $tools.POD5 convert fast5 $input.fast5 --output $output.pod5
        """
    }
}

dorado = {

    output.dir='dorado/' + branch.name

    uses(dorados: 1) {
        exec """
            set -o pipefail

            echo "Processing $inputs.x5 with dorado"

            ${inputs.x5.collect { file(it) }*.absoluteFile.collect { "ln -sf $it $output.dir/$it.name;"}.join("\n") }

            $tools.DORADO basecaller $DRD_MODELS_PATH/$model.params.drd_model $output.dir --modified-bases 5mCG_5hmCG | 
                $tools.SAMTOOLS view -b -o $output.ubam -
        """
    }
}

make_mmi = {
    output.dir = new File(REF).parentFile.absolutePath

    def map_platform = (lrs_platform == 'hifi') ? 'map-hifi' : 'map-ont'

    // nb: this could write into the reference directory
    // which is a little naughty
    produce(REF_MMI) {
        exec """
            $tools.MINIMAP2 -t ${threads} -x $map_platform -d $output ${REF}
        """
    }
}

demultiplex = {

    doc "Runs dorado demultiplexing"

    try {
        //demultiplex_limit.acquire()

        output.dir="demux/${branch.name}"

        produce("summary_counts.txt", "*.fastq.gz") {

            exec """
               mkdir $output.dir/input_tmp

               cat $inputs.fastq.gz > $output.dir/input_tmp/combined_fastq_tmp.fastq.gz

               $tools.DORADO demux 
                    --barcode-both-ends
                    --barcode-sequences $DESIGN/barcodes.fa
                    --barcode-arrangement $DESIGN/arrangement.toml  
                    --emit-fastq  
                    --output-dir $output.dir  
                    $output.dir/input_tmp/combined_fastq_tmp.fastq.gz

                 TOTAL_LINES=\$(zcat $inputs.fastq.gz  | wc --lines)

                 wc $output.dir/*.fastq | awk '{ print (\$1 / '\$TOTAL_LINES' ) "\\t" \$4 }' > $output.txt

                 rm -v $output.dir/input_tmp/combined_fastq_tmp.fastq.gz

                 pigz --best -p ${threads} $output.dir/*.fastq
            """
        }
    }
    finally {
        //demultiplex_limit.release()
    }
}



rename_and_merge_demux_output = {

    branch.barcode = sample_info.find { it.Sample_ID == sample }.barcodeno

    output.dir = "fastq"

    if(sample.startsWith('NTC') || sample.endsWith('NTC')) {
        succeed "Not merging fastq for NTC"
    }

    from("*_${barcode}.fastq.gz") produce("${sample}.fastq.gz") {
        exec """
            echo "Renaming data for $sample with barcode no. $barcode"

            zcat $inputs.fastq.gz > $output.fastq.gz
        """
    }
}

add_sample_read_group = {
    def SAMTOOLS = tools.SAMTOOLS

    output.dir = "align"

    filter("sample_rg") {
        exec """
            set -eo pipefail

            STR=\$($SAMTOOLS view -@ $threads $input.bam | 
                awk '{for(i=12;i<=NF;i++){split(\$i,a,":"); if(a[1]!="RG") print a[1]}}' | 
                sort -u | 
                paste -sd',')

            $SAMTOOLS view -h --keep-tag \$STR $input.bam | 
                $SAMTOOLS addreplacerg -w -r '@RG\tID:${sample}\tSM:${sample}' -@ $threads -o $output.bam -

            $SAMTOOLS index -@ $threads $output.bam
        """, "replace_read_group"
    }
}

unmap_bam = {
    requires sample : 'the sample id being processed'
    
    def SAMTOOLS = tools.SAMTOOLS

    output.dir = 'align'

    def tmp_bam = "$output.dir/${sample}.sample_RG.bam"

    transform('.bam') to('.unmap.ubam') {
        exec """
            set -eo pipefail

            STR=\$($SAMTOOLS view -@ $threads $input.bam | 
                awk '{for(i=12;i<=NF;i++){split(\$i,a,":"); if(a[1]!="RG") print a[1]}}' | 
                sort -u | 
                paste -sd',')

            $SAMTOOLS view -h --keep-tag \$STR $input.bam | 
                $SAMTOOLS addreplacerg -w -r '@RG\tID:${sample}\tSM:${sample}' -@ $threads -o $tmp_bam -
        """, "replace_read_group"

        exec """
            set -eo pipefail

            gatk RevertSam -I $tmp_bam
                -O $output.bam
                --KEEP_FIRST_DUPLICATE
                --SANITIZE
                --OQ
                --VALIDATION_STRINGENCY LENIENT
                --MAX_RECORDS_IN_RAM 20000
                --TMP_DIR $TMPDIR

            rm $tmp_bam
        """, "revert_bam"
    }
}

minimap2_align = {

    doc "Align FASTQ stored in unaligned BAM file (ubam)"

    requires sample : 'the sample id being processed'

    var bam_ext : 'bam'

    def cram_ext = bam_ext.replaceFirst(/bam$/,'cram')

    def map_platform = (lrs_platform == 'hifi') ? 'map-hifi' : 'map-ont'
    def rg_platform = (lrs_platform == 'hifi') ? 'PacBio' : 'ONT'
    def rg_lib = (lrs_platform == 'hifi') ? 'LIB_HIFI' : 'LIB_ONT'

    def SAMTOOLS = tools.SAMTOOLS
    
    output.dir = 'align'

    def output_pass_cram = "${sample}.pass." + cram_ext
    def output_fail_cram = "${sample}.fail." + cram_ext

    produce(output_pass_cram, output_fail_cram) {
        exec """
            mkdir -p $output.dir/tmp

            $SAMTOOLS bam2fq -@ $threads -T 1 $input.ubam
                | $tools.MINIMAP2 -y -t $threads -ax $map_platform -R "@RG\\tID:${sample}\\tPL:${rg_platform}\\tPU:1\\tLB:${rg_lib}\\tSM:${sample}" $REF_MMI - 
                | $SAMTOOLS sort -@ $threads -T $output.dir/tmp
                | tee >($SAMTOOLS view -e '[qs] < $calling.qscore_filter' -C -T $REF -o ${output.fail[cram_ext]} - )
                | $SAMTOOLS view -e '(![qs] && [qs] != 0) || [qs] >= $calling.qscore_filter' -C -T $REF -o ${output.pass[cram_ext]} -

            $SAMTOOLS index ${output.pass[cram_ext]}
        """
    }

    forward(output.pass[cram_ext])
}

minimap2_align_fastq = {
    
    def SAMTOOLS = tools.SAMTOOLS
    
    output.dir = 'align'

    exec """
            $tools.MINIMAP2 -y -t $threads -ax map-ont $REF_MMI $input.fastq.gz
            | $SAMTOOLS sort -@ $threads > $output.bam

        $SAMTOOLS index $output.bam

    """
}


merge_bams = {

    doc "Merge arbitrary BAMs"

    output.dir = 'merged'

    produce(opts.sample + '.merged.bam') {
	    exec """
		    $tools.SAMTOOLS merge $output.bam $inputs.bam 
		    -f 
		    --no-PG 
		    --write-index 
		    --reference $REF 
		    --threads 4
	    """
    }
}


merge_pass_calls = {
    var bam_ext : 'bam'

    def cram_ext = bam_ext.replaceFirst(/bam$/,'cram')

    output.dir = 'align'

    def output_pass_bam = "${sample}.merged.pass." + bam_ext

    produce(output_pass_bam) {
        exec """
            $tools.SAMTOOLS merge ${output[bam_ext]} ${inputs.pass[cram_ext]}
            -f 
            -c 
            -p 
            --no-PG 
            --reference $REF 
            --threads $threads

            $tools.SAMTOOLS index -@ $threads ${output[bam_ext]}
        """
    }
}

zip_ref = {
    output.dir = 'align/ref'

    def ref_fasta = new File(REF).name

    produce(ref_fasta + '.gz') {
        exec """
            bgzip -c $REF > $output.gz
        """
    }
}

mosdepth = {
    var bam_ext : 'bam'
   
    output.dir = 'qc/mosdepth'
    
    produce("${sample}.regions.bed.gz",
            "${sample}.mosdepth.bed",
            "${sample}.mosdepth.global.dist.txt",
            "${sample}.mosdepth.summary.txt",
            "${sample}.thresholds.bed.gz") {

        exec """
            export REF_PATH=$REF

            export MOSDEPTH_PRECISION=3

            cut -f 1,2,3 $opts.targets > $output.mosdepth.bed

            mosdepth 
            -x 
            -t $threads
            -b $output.mosdepth.bed
            --thresholds 1,10,20,30
            --no-per-base
            $output.dir/${sample}
            ${input[bam_ext]}
        """
    }
}

read_stats = {
    output.dir = "qc/bamstats"

    produce("${sample}.readstats.tsv.gz") {
        exec """

            set -o pipefail

            $tools.BAMSTATS --threads $threads $input.bam | gzip > $output.gz
        """
    }
}

call_short_variants = {
    var bam_ext : 'bam'

    output.dir = new File("variants/${sample}").absolutePath

    produce("${sample}.wf_snp.g.vcf.gz") {
        def gvcfFlags = ""

        if(calling.enable_gvcf) {
            gvcfFlags = "--gvcf" 
        }

        exec """
            set -uo pipefail

            /opt/bin/run_clair3.sh --bam_fn=${input[bam_ext]}
                                   --bed_fn=$opts.targets
                                   --ref_fn=$REF
                                   --threads=$threads
                                   --platform=$lrs_platform
                                   --model_path=$calling.CLAIR3_MODELS_PATH/${clair3_model.clair3_model_name}
                                   --output=$output.dir
                                   --sample_name=$sample
                                   --snp_min_af=$calling.snp_min_af
                                   --indel_min_af=$calling.indel_min_af
                                   --min_mq=$calling.min_mq
                                   --min_coverage=$calling.min_cov
                                   --var_pct_phasing=$calling.phasing_pct
                                   --chunk_size=$calling.chunk_size
                                   $gvcfFlags

            $BASE/scripts/fix_clair3_gvcf_v2.py
                -i ${output.dir}/merge_output.gvcf.gz
                -o $output.g.vcf.gz.prefix
                -r $REF

            bgzip -c $output.g.vcf.gz.prefix > $output.g.vcf.gz

            tabix -p vcf $output.g.vcf.gz

            rm ${output.dir}/merge_output.gvcf.gz* $output.g.vcf.gz.prefix
        """
    }
}

normalize_gvcf = {

    doc "split VCF lines so that each line contains one and only one variant, and left-normalize all VCF lines"

    output.dir="variants/${sample}"
    transform('g.vcf.gz') to ("norm.g.vcf.gz") {
        exec """
            set -o pipefail

            gunzip -c $input.gz | bcftools norm -m -both - | bcftools norm -f $REF - | bgzip -c > $output

            tabix -p vcf $output
        """
    }
}

gvcf_to_vcf = {
    output.dir = "variants"

    transform('g.vcf.gz') to ('vcf.gz') {
        exec """
            set -o pipefail

            bcftools view -i 'TYPE="snp" || TYPE="indel"' $input.g.vcf.gz | bgzip -c > $output.vcf.gz

            tabix -p vcf $output.vcf.gz
        """
    }
}

phase_variants = {
    var bam_ext : 'bam'

    output.dir = "variants/${sample}"

    transform("norm.g.vcf.gz") to("norm.phased.g.vcf.gz") {
        def tmp_vcf = "${input.vcf.gz.prefix.prefix.prefix}.tmp.g.vcf"

        def platform_opt = (lrs_platform == 'hifi') ? '--pb' : '--ont'

        exec """
            set -uo pipefail

            echo "Using longphase for final variant phasing"

            bgzip -@ $threads -dc $input.vcf.gz > $tmp_vcf  

            $tools.LONGPHASE phase $platform_opt
                -o ${output.prefix.prefix}
                -s $tmp_vcf
                -b ${input[bam_ext]}
                -r $REF
                -t $threads

            rm -f $tmp_vcf

            bgzip -c $output.prefix > $output

            tabix -f -p vcf $output

            rm -f $output.prefix
         """

         sample_vcfs.get(sample, []).add(output.toString())
    }
}

haplotag_bam = {
    var bam_ext : 'bam'

    def cram_ext = bam_ext.replaceFirst(/bam$/,'cram')

    output.dir = "align"

    def output_cram_ext = 'haplotagged.' + cram_ext

    transform(bam_ext) to(output_cram_ext) {
        exec """
            set -o pipefail

            $tools.LONGPHASE haplotag 
            -r $REF
            -s $input.norm.phased.vcf.gz
            -b ${input[bam_ext]}
            --cram
            -t $threads 
            -o $output.dir/${file(output[cram_ext].prefix).name}

            samtools index ${output[cram_ext]}
        """
    }
}

combine_family_vcfs = {
    
    requires family : 'family to process'

    var XIMMER_GNGS_JAR : "$tools.XIMMER/tools/groovy-ngs-utils/1.0.9/groovy-ngs-utils.jar"

    def family_samples = meta*.value.grep { println(it);  it.family_id == family }
    
    println "Samples to process for $family are ${family_samples*.identifier}"
    
    def family_vcfs = family_samples*.identifier.collectEntries { [ it,  sample_vcfs[it]] }
    
    def samples_missing_vcfs = family_vcfs.grep { it.value == null }*.key
    if(samples_missing_vcfs)
        fail "Sample samples $samples_missing_vcfs from family $family are missing gVCFs. Please check appropriate inputs have been provided"
    
    output.dir = "variants/${family}"
    
    println "Inputs are: " + family_vcfs*.value.flatten()
    
    from(family_vcfs*.value.flatten()) produce("${family}.family.vcf.gz") {
        exec """
            set -o pipefail

            export JAVA=$tools.JAVA

            $tools.GROOVY -cp $XIMMER_GNGS_JAR $BASE/src/MergeFamilyVCF.groovy $inputs.g.vcf.gz | 
                bcftools view -i 'TYPE="snp" || TYPE="indel"' - | 
                bgzip -c > $output.vcf.gz

            tabix -p vcf $output.vcf.gz
        """
    }
}


dorado_demux = {

    output.dir='dorado'

    // ${inputs.x5.collect { file(it) }*.absoluteFile.collect { "ln -sf $it $output.dir/inputs/$it.name;"}.join("\n") }

    uses(dorados: 1) {
		exec """
			$tools.DORADO basecaller
				$DRD_MODELS_PATH/$model.params.drd_model
				$opts.dir
				--barcode-both-ends
				--barcode-sequences $DESIGN/barcodes.fa
				--barcode-arrangement $DESIGN/arrangement.toml  
				> $output.ubam

		"""
	}
}

dorado_demux_classify = {

    output.dir = "dorado/demultiplexed"

    def files = sample_info*.Sample_ID.collect { "${it}.ubam" }

    produce(files) {
        exec """
            $tools.DORADO demux --output-dir $output.dir --no-classify $input.ubam

            ${sample_info.collect { "mv -v twist_${it.barcodeno}.bam ${it.Sample_ID}.ubam;"}.join('\n')}
        """
    }
}

rename_and_merge_demux_output = {

    branch.barcode = sample_info.find { it.Sample_ID == sample }.barcodeno

    output.dir = "fastq"

    if(sample.startsWith('NTC') || sample.endsWith('NTC')) {
        succeed "Not merging fastq for NTC"
    }

    from("*_${barcode}.fastq.gz") produce("${sample}.fastq.gz") {
        exec """
            echo "Renaming data for $sample with barcode no. $barcode"

            zcat $inputs.fastq.gz > $output.fastq.gz
        """
    }
}
