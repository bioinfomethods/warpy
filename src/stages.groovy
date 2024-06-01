

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

    transform('.x5') to ('.ubam') {
        uses(dorados: 1) {
            exec """
                set -o pipefail

                ${inputs.x5.collect { file(it) }*.absoluteFile.collect { "ln -sf $it $output.dir/$it.name;"}.join("\n") }

                $tools.DORADO basecaller $DRD_MODELS_PATH/$model.params.drd_model $output.dir/${file(input.x5).name} --modified-bases 5mCG_5hmCG | 
                    $tools.SAMTOOLS view -b -o $output.ubam -
            """
        }
    }
}

make_mmi = {
    output.dir = BASE + "/data/hg38"

    produce(REF_MMI) {
        exec """
            $tools.MINIMAP2 -t ${threads} -x map-ont -d $output ${REF}
        """
    }
}

minimap2_align = {
    
    def SAMTOOLS = tools.SAMTOOLS
    
    output.dir = 'align'

    exec """
        $SAMTOOLS bam2fq -@ $threads -T 1 $input.ubam
            | $tools.MINIMAP2 -y -t $threads -ax map-ont -R "@RG\\tID:${opts.sample}\\tPL:ONT\\tPU:1\\tLB:ONT_LIB\\tSM:${opts.sample}" $REF_MMI - 
            | $SAMTOOLS sort -@ $threads
            | tee >($SAMTOOLS view -e '[qs] < $calling.qscore_filter' -o $output.fail.bam - )
            | $SAMTOOLS view -e '[qs] >= $calling.qscore_filter' -o $output.pass.bam -

        $SAMTOOLS index $output.pass.bam

    """
}

merge_pass_calls = {

    output.dir = 'align'

    produce(opts.sample + '.pass.bam') {
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
    
    produce("${opts.sample}.regions.bed.gz",
            "${opts.sample}.mosdepth.bed",
            "${opts.sample}.mosdepth.global.dist.txt",
            "${opts.sample}.mosdepth.summary.txt",
            "${opts.sample}.thresholds.bed.gz") {

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
            $output.dir/${opts.sample}
            $input.bam
        """
    }
}

read_stats = {
    output.dir = "qc/bamstats"

    produce("${opts.sample}.readstats.tsv.gz") {
        exec """
            $tools.BAMSTATS --threads $threads $input.bam | gzip > $output.gz
        """
    }
}

pileup_variants = {
    
    output.dir="clair_output/pileup"

    var GVCF : false
    
    
    println("Clair chunk: " + clair_chunk)
    
    gngs.Region region = new gngs.Region(clair_chunk.toString())
    
    produce("${opts.sample}_${region.chr}_${region.from}.vcf", "${opts.sample}_${region.chr}_${region.from}.txt") {


        uses(clair3: 1) {
            exec """
                set -uo pipefail

                export REF_PATH=cram_cache/%2s/%2s/%s

                python $tools.CLAIR3/clair3.py CallVariantsFromCffi
                    --chkpnt_fn $CLAIR3_MODELS_PATH/${clair3_model.clair3_model_name}/pileup
                    --bam_fn $input.bam
                    --bed_fn $opts.targets
                    --call_fn $output.vcf.optional
                    --ref_fn $REF
                    --ctgName $region.chr
                    --ctgStart $region.from
                    --ctgEnd $region.to
                    --platform ont 
                    --fast_mode False
                    --snp_min_af $calling.snp_min_af
                    --indel_min_af $calling.indel_min_af
                    --minMQ $calling.min_mq
                    --minCoverage $calling.min_cov
                    --call_snp_only False
                    --gvcf $GVCF
                    --enable_long_indel False
                    --temp_file_dir gvcf_tmp_path
                    --pileup

                touch -a $output.vcf.optional


                echo "`date` : succesfully called variants" > $output.txt
            """
        }
    }
}

aggregate_pileup_variants = {
    
    output.dir="variants/$opts.sample"
    
    from('CONTIGS') produce("${opts.sample}.aggregate.pileup.vcf.gz") {
        exec """
            set -uo pipefail

            $tools.PYPY $tools.CLAIR3/clair3.py SortVcf
                --contigs_fn $input
                --input_dir ${file(input1.vcf).parentFile.path}
                --vcf_fn_prefix $opts.sample
                --output_fn $output.vcf.gz.prefix
                --sampleName $opts.sample
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

    branch.dir="split_folder/$opts.sample"
    
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

                conda deactivate

                unset __conda_setup
         """
    }
}


get_qual_filter = {
    
    output.dir =  "variants/$opts.sample"
    
    exec """
        set -uo pipefail

        echo "[INFO] 5/7 Select candidates for full-alignment calling"

        bgzip -fdc $input.vcf.gz |
        $tools.PYPY $tools.CLAIR3/clair3.py SelectQual
                --output_fn $output.dir
                --var_pct_full $calling.var_pct_full
                --ref_pct_full $calling.ref_pct_full
                --platform ont 

        mv $output.dir/qual $output.vcf
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
    
    output.dir = branch.dir + '/candidates/' + chr

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
    
    produce("${opts.sample}.${bed_prefix}.full_alignment.vcf") {
        exec """
            set -uo pipefail

            echo "[INFO] 6/7 Call low-quality variants using full-alignment model"

            python $tools.CLAIR3/clair3.py CallVariantsFromCffi
                --chkpnt_fn $CLAIR3_MODELS_PATH/${clair3_model.clair3_model_name}/full_alignment
                --bam_fn $input.bam 
                --call_fn $output.vcf
                --sampleName ${opts.sample}
                --ref_fn $REF
                --full_aln_regions $input.bed
                --ctgName $chr
                --add_indel_length
                --gvcf false
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
    
    output.dir="variants/$opts.sample"
    
    from('CONTIGS', '*.full_alignment.vcf') produce("${opts.sample}.full_alignment.vcf.gz"){
        exec """

            echo "First input VCF is $input1.vcf"

            $tools.PYPY $tools.CLAIR3/clair3.py SortVcf
                --input_dir ${file(input1.vcf).parentFile.path}
                --output_fn $output.vcf.gz.prefix
                --sampleName $opts.sample
                --ref_fn $REF
                --contigs_fn $input1
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

    produce("${opts.sample}.merged_methods.${chr}.vcf.gz") {
        exec """

            echo "[INFO] 7/7 Merge pileup VCF and full-alignment VCF"

            $tools.PYPY $tools.CLAIR3/clair3.py MergeVcf
                --pileup_vcf_fn $input.aggregate.pileup.vcf.gz
                --bed_fn_prefix candidates/$chr
                --full_alignment_vcf_fn $input.full_alignment.vcf.gz
                --output_fn $output.vcf.gz.prefix
                --platform ont
                --print_ref_calls False
                --gvcf false
                --haploid_precise False
                --haploid_sensitive False
                --non_var_gvcf_fn non_var.gvcf
                --ref_fn $REF
                --ctgName $chr

            bgzip -c $output.vcf.gz.prefix > $output.vcf.gz

            tabix $output.vcf.gz
        """
    }
}

aggregate_all_variants = {
    
    output.dir = "variants/${opts.sample}"

    from("CONTIGS") produce("${opts.sample}.wf_snp.vcf.gz") {
        exec """

            rm -rf $output.dir/merge_output

            mkdir $output.dir/merge_output

            for i in $inputs.vcf.gz ; do cp \$i $output.dir/merge_output ; done 

            $tools.PYPY $tools.CLAIR3/clair3.py SortVcf
                --input_dir variants
                --vcf_fn_prefix $opts.sample
                --output_fn $output.vcf.gz.prefix
                --sampleName $opts.sample
                --ref_fn $REF
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
