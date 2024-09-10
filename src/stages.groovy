

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

    // nb: this could write into the reference directory
    // which is a little naughty
    produce(REF_MMI) {
        exec """
            $tools.MINIMAP2 -t ${threads} -x map-ont -d $output ${REF}
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


minimap2_align = {

    doc "Align FASTQ stored in unaligned BAM file (ubam)"

    requires sample : 'the sample id being processed'
    
    def SAMTOOLS = tools.SAMTOOLS
    
    output.dir = 'align'

    produce("${sample}.pass.bam", "${sample}.fail.bam") {
        exec """
            $SAMTOOLS bam2fq -@ $threads -T 1 $input.ubam
                | $tools.MINIMAP2 -y -t $threads -ax map-ont -R "@RG\\tID:${sample}\\tPL:ONT\\tPU:1\\tLB:ONT_LIB\\tSM:${sample}" $REF_MMI - 
                | $SAMTOOLS sort -@ $threads
                | tee >($SAMTOOLS view -e '[qs] < $calling.qscore_filter' -o $output.fail.bam - )
                | $SAMTOOLS view -e '[qs] >= $calling.qscore_filter' -o $output.pass.bam -

            $SAMTOOLS index $output.pass.bam
        """
    }

    forward(output.pass.bam)
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

    output.dir = 'align'

    produce(sample + '.merged.pass.bam') {
        exec """
            $tools.SAMTOOLS merge $output.bam $inputs.pass.bam 
            -f 
            -c 
            -p 
            --no-PG 
            --write-index 
            --reference $REF 
            --threads $threads
        """
    }
}

mosdepth = {
   
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
            $input.bam
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

pileup_variants = {
    
    output.dir="clair3_output/pileup"
   
    println("Clair chunk: " + region)

    def chunk = new gngs.Region(region.toString())
    
    // gngs.Region region = new gngs.Region(clair_chunk.toString())
    
    produce("${sample}_${chunk.chr}_${chunk.from}.vcf", "${sample}_${chunk.chr}_${chunk.from}.txt") {


        uses(clair3: 1) {
            exec """
                set -uo pipefail

                export REF_PATH=cram_cache/%2s/%2s/%s

                python $tools.CLAIR3/clair3.py CallVariantsFromCffi
                    --chkpnt_fn $calling.CLAIR3_MODELS_PATH/${clair3_model.clair3_model_name}/pileup
                    --bam_fn $input.bam
                    --bed_fn $opts.targets
                    --call_fn $output.vcf.optional
                    --ref_fn $REF
                    --ctgName $chunk.chr
                    --ctgStart $chunk.from
                    --ctgEnd $chunk.to
                    --platform ont 
                    --fast_mode False
                    --snp_min_af $calling.snp_min_af
                    --indel_min_af $calling.indel_min_af
                    --minMQ $calling.min_mq
                    --minCoverage $calling.min_cov
                    --call_snp_only False
                    --gvcf ${calling.enable_gvcf ? "True" : "False"}
                    --enable_long_indel False
                    --temp_file_dir $output.dir/${sample}_${chunk.chr}_${chunk.from}_gvcf_tmp_path
                    --pileup

                touch -a $output.vcf.optional


                echo "`date` : succesfully called variants" > $output.txt
            """
        }
    }
}

aggregate_pileup_variants = {
    
    output.dir="variants/$sample"
    
    from('CONTIGS') produce("${sample}.aggregate.pileup.vcf.gz") {
        exec """
            set -uo pipefail

            $tools.PYPY $tools.CLAIR3/clair3.py SortVcf
                --contigs_fn $input
                --input_dir ${file(input1.vcf).parentFile.path}
                --vcf_fn_prefix $sample
                --output_fn $output.vcf.gz.prefix
                --sampleName $sample
                --print_ref_calls ${calling.enable_gvcf ? "True" : "False"}
                --ref_fn $REF
                --cmd_fn $output.dir/tmp/CMD

            if [ "\$( bgzip -@ $threads -fdc $output.vcf.gz | grep -v '#' | wc -l )" -eq 0 ]; 
            then echo "[INFO] Exit in pileup variant calling"; exit 1; fi

            bgzip -@ $threads -fdc $output.vcf.gz |
                $tools.PYPY $tools.CLAIR3/clair3.py SelectQual --var_pct_phasing $calling.phasing_pct --phase --output_fn .
        """
    }
}


select_het_snps = {

    branch.dir="split_folder/$sample"
    
    transform('.vcf.gz') to('.het_snps.vcf.gz') {
        exec """
            set -uo pipefail

            $tools.PYPY $tools.CLAIR3/clair3.py SelectHetSnp
                --vcf_fn $input.vcf.gz
                --split_folder $output.dir
                --ctgName $chr

            bgzip -c $output.dir/${chr}.vcf > $output.vcf.gz

            tabix $output.vcf.gz
        """
    }
}



phase_contig = {

    transform('.vcf.gz') to('.phased.vcf.gz') {

        def tmp_vcf = "${input.vcf.gz.prefix.prefix}.tmp.vcf"
        exec """
                set -uo pipefail

                echo "Using longphase for phasing"

                bgzip -@ $threads -dc $input.vcf.gz > $tmp_vcf  

                $tools.LONGPHASE phase 
                    --ont 
                    -o ${output.prefix.prefix}
                    -s $tmp_vcf
                    -b $input.bam 
                    -r $REF
                    -t $threads

                rm -f $output

                bgzip $output.prefix

                tabix -f -p vcf $output
         """
    }
}


get_qual_filter = {
    
    doc """
        Full alignment calling is expensive, so only a configurable proportion of variants are selected for it.
        The proportion is fixed, but to make it easier to implement downstream, this step determines the appropriate
        qual score threshold to use to get the top X% of low quality variants in subsequent steps.
        """
    
    output.dir =  "variants/$sample"
    
    exec """
        set -uo pipefail

        echo "[INFO] 5/7 Determine qual score cutoffs to to select fixed proportion of candidates for full-alignment calling"

        bgzip -fdc $input.vcf.gz |
        $tools.PYPY $tools.CLAIR3/clair3.py SelectQual
                --output_fn $output.dir
                --var_pct_full $calling.var_pct_full
                --ref_pct_full $calling.ref_pct_full
                --platform ont 

        mv $output.dir/qual $output.qual.txt
    """
}

create_candidates = {
    
    doc """
        This creates BED files as candidate_bed/<ctg>.0_14 with candidates
        along with a file the FULL_ALN_FILE_<ctg> listing all of the BED
        files.  All we really want are the BEDs, the file of filenames is
        used for the purposes of parallel in the original workflow.

        https://github.com/HKU-BAL/Clair3/blob/329d09b39c12b6d8d9097aeb1fe9ec740b9334f6/scripts/clair3.sh#L218
    """
    
    output.dir = branch.dir + "/clair3_output/full_aln_candidates/$sample/" + chr

    from('vcf.gz') produce('*.bed') {
        exec """
            set -uo pipefail

            $tools.PYPY $tools.CLAIR3/clair3.py SelectCandidates
                --pileup_vcf_fn $input.vcf.gz
                --split_folder $output.dir
                --ref_fn $REF
                --var_pct_full $calling.var_pct_full
                --ref_pct_full $calling.ref_pct_full
                --platform ont
                --ctgName $chr

            rm -f $output.dir/*.bed

            src_bed_file_name_pattern="^$output.dir/$chr.+\$"

            for i in $output.dir/*.*;
            do 
                if [[ \$i =~ \$src_bed_file_name_pattern ]]; then
                    mv -v \$i \${i}.bed;
                fi
            done
        """
    }
}

evaluate_candidates = {

    def bed_prefix = file(input.bed).name.replaceAll('.bed$','')

    output.dir = "variants/full_alignments"
    
    produce("${sample}.${bed_prefix}.full_alignment.vcf") {
        exec """
            set -uo pipefail

            echo "[INFO] 6/7 Call low-quality variants using full-alignment model"

            python $tools.CLAIR3/clair3.py CallVariantsFromCffi
                --chkpnt_fn $calling.CLAIR3_MODELS_PATH/${clair3_model.clair3_model_name}/full_alignment
                --bam_fn $input.bam 
                --call_fn $output.vcf
                --sampleName ${sample}
                --ref_fn $REF
                --full_aln_regions $input.bed
                --ctgName $chr
                --add_indel_length
                --gvcf ${calling.enable_gvcf?"True":"False"}
                --minMQ $calling.min_mq
                --minCoverage $calling.min_cov
                --snp_min_af $calling.snp_min_af
                --indel_min_af $calling.indel_min_af
                --platform ont 
                --phased_vcf_fn $input.vcf.gz
        """
    }
}

aggregate_full_align_variants = {
    
    output.dir="variants/$sample"
    
    from('CONTIGS', '*.full_alignment.vcf') produce("${sample}.full_alignment.vcf.gz"){
        exec """

            echo "First input VCF is $input1.vcf"

            $tools.PYPY $tools.CLAIR3/clair3.py SortVcf
                --input_dir ${file(input1.vcf).parentFile.path}
                --output_fn $output.vcf.gz.prefix
                --vcf_fn_prefix $sample
                --sampleName $sample
                --ref_fn $REF
                --contigs_fn $input1
                --print_ref_calls True
                --cmd_fn $output.dir/tmp/CMD

            if [ "\$( bgzip -fdc $output.vcf.gz | grep -v '#' | wc -l )" -eq 0 ]; then
                echo "[INFO] Exit in full-alignment variant calling";
                exit 0;
            fi
        """
    }
}

merge_pileup_and_full_vars = {
    
    output.dir = "variants"
    
    List output_files = ["${sample}.merged_methods.${chr}.vcf"]
    if(calling.enable_gvcf) {
        output_files.add "$output.dir/gvcf/${sample}.merged_methods.${chr}.gvcf"
    }

    produce(output_files) {
        
        def gvcfFlags = ""
        if(calling.enable_gvcf) {
            gvcfFlags = "--print_ref_calls True --gvcf True --gvcf_fn $output.gvcf" 
        }
        
        exec """
            echo "[INFO] 7/7 Merge pileup VCF and full-alignment VCF"

            mkdir -p $output.dir/gvcf

            $tools.PYPY $tools.CLAIR3/clair3.py MergeVcf
                --pileup_vcf_fn $input.aggregate.pileup.vcf.gz
                --bed_fn_prefix candidates/$chr
                --full_alignment_vcf_fn $input.full_alignment.vcf.gz
                --output_fn $output.vcf
                --platform ont $gvcfFlags
                --haploid_precise False
                --haploid_sensitive False
                --non_var_gvcf_fn ${sample}-non_var.gvcf
                --ref_fn $REF
                --ctgName $chr
        """
    }
}

aggregate_all_variants = {
    
    output.dir = "variants/${sample}"

    from("CONTIGS") produce("${sample}.wf_snp.vcf.gz") {
        exec """

            rm -rf $output.dir/merge_output

            mkdir $output.dir/merge_output

            for i in $inputs.vcf.gz ; do cp \$i $output.dir/merge_output ; done 

            $tools.PYPY $tools.CLAIR3/clair3.py SortVcf
                --input_dir variants
                --vcf_fn_prefix $sample
                --output_fn $output.vcf.gz.prefix
                --sampleName $sample
                --ref_fn $REF
                --print_ref_calls True
                --contigs_fn $input1
                --cmd_fn $output.dir/tmp/CMD

            if [ "\$( bgzip -fdc $output.vcf.gz | grep -v '#' | wc -l )" -eq 0 ]; then
                echo "[INFO] Exit in all contigs variant merging : no variants?";
                exit 0;
            fi

            echo "[INFO] Finish calling, output file: $output.vcf.gz"
        """
    }
}


normalize_vcf = {

    doc "split VCF lines so that each line contains one and only one variant, and left-normalize all VCF lines"

    output.dir="variants"
    transform('vcf.gz') to ("norm.vcf.gz") {
        exec """
            set -o pipefail

            gunzip -c $input.gz | bcftools norm -m -both - | bcftools norm -f $REF - | bgzip -c > $output

            tabix -p vcf $output
        """
        
        sample_vcfs.get(sample, []).add(output.toString())
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

            $tools.GROOVY -cp $XIMMER_GNGS_JAR $BASE/src/MergeFamilyVCF.groovy $inputs.vcf.gz | bgzip -c >  $output.vcf.gz

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
