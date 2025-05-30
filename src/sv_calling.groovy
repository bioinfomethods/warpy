

filterBam = {
    var bam_ext: 'bam'

    output.dir = "align"
    
    def output_bam_ext = 'filtered.' + bam_ext

    transform(bam_ext) to(output_bam_ext) {
        exec """
            samtools view -@ $threads 
                ${input[bam_ext]}
                -F 260 
                -o ${output[bam_ext]}
                --write-index 
                --reference $REF
        """

        sample_bams.get(sample, []).add(output[bam_ext].toString())
    }
}


// NOTE VCF entries for alleles with no support are removed to prevent them from
//      breaking downstream parsers that do not expect them
sniffles2 = {
    
    var sniffles_args : ''
    
    branch.dir = "sv/$sample"

    produce("${sample}.sniffles.vcf") {
        exec """
            export REF_PATH=$REF

            sniffles
                --threads $threads
                --sample-id ${sample}
                --reference $REF
                --output-rnames
                --cluster-merge-pos $calling.cluster_merge_pos
                --input $input.bam
                --allow-overwrite
                --tandem-repeats ${calling.tr_bed} $sniffles_args
                --vcf ${output.vcf.prefix}.tmp.vcf

            sed '/.:0:0:0:NULL/d' ${output.vcf.prefix}.tmp.vcf | 
            awk -f $BASE/scripts/fix_allele_seq.awk > $output.vcf
        """
    }
}

sniffles2_for_trios = {
    var sniffles_args : ''
    
    branch.dir = "sv/snf"

    produce("${sample}.sniffles.snf") {
        exec """
            export REF_PATH=$REF

            sniffles
                --threads $threads
                --sample-id ${sample}
                --reference $REF
                --output-rnames
                --allow-overwrite
                --cluster-merge-pos $calling.cluster_merge_pos
                --input $input.bam
                --tandem-repeats ${calling.tr_bed} $sniffles_args
                --snf $output.snf
        """

        sample_snfs.get(sample, []).add(output.snf.toString())
    }
}

sniffles2_joint_call = {
    requires family : 'family to process'

    var sniffles_args : ''

    def family_samples = meta*.value.grep { println(it);  it.family_id == family }

    println "Samples to process for $family are ${family_samples*.identifier}"

    def family_snfs = family_samples*.identifier.collectEntries { [ it,  sample_snfs[it]] }
    
    def samples_missing_snfs = family_snfs.grep { it.value == null }*.key
    if(samples_missing_snfs)
        fail "Sample samples $samples_missing_snfs from family $family are missing SNFs (sample SV output file). Please check appropriate inputs have been provided"

    output.dir = "sv/$family"

    println "Inputs are: " + family_snfs*.value.flatten()
    
    from(family_snfs*.value.flatten()) produce("${family}.sniffles.family.sv.vcf.gz") {
        exec """
            sniffles
                --threads $threads
                --input $inputs
                --vcf $output.prefix

            $BASE/scripts/vcfsort -T $TMPDIR -N $threads $output.prefix | 
            awk -f $BASE/scripts/fix_allele_seq.awk | bgzip -c - > $output.vcf.gz

            bcftools index -t $output.vcf.gz
        """
    }
}

init_jasmine_bams = {
    def family_samples = meta*.value.grep { println(it); it.family_id == family }
    def family_sample_identifiers = family_samples.collect { it.identifier }

    def bamsListings = sample_bams
      .findAll { k, v -> k in family_sample_identifiers }
      .sort()
      .collect { k, v -> v }
      .flatten()
      .unique()
      .join('\\n')

    produce("${family}.bam.listings.txt") {
        output.dir = "sv/$family"

        groovy """
            new File('$output.bam.listings.txt').text = '''$bamsListings'''
        """
    }
}

init_jasmine_vcfs = {
    var sv_tool: 'sniffles'

    def family_samples = meta*.value.grep { println(it); it.family_id == family }
    def family_sample_identifiers = family_samples.collect { it.identifier }

    def vcfsListings = sample_sv_vcfs
      .findAll { k, v -> k.split("#")[0] == sv_tool && k.split("#")[1] in family_sample_identifiers }
      .sort()
      .collect { k, v -> v }
      .flatten()
      .unique()
      .join('\\n')

    def vcf_list = (sv_tool == 'sniffles'? "${family}.vcf.listings.txt":"${family}.${sv_tool}.vcf.listings.txt")

    output.dir = "sv/$family"

    produce(vcf_list) {
        groovy """
            new File('$output.vcf.listings.txt').text = '''$vcfsListings'''
        """
    }
}

jasmine_merge = {
    doc 'Joint call SVs for a family using Jasmine'

    requires family : 'family to process'

    var sv_tool: 'sniffles'

    def out_vcf = (sv_tool == 'sniffles'? "${family}.jasmine.family.sv.vcf.gz":"${family}.${sv_tool}.jasmine.family.sv.vcf.gz")

    output.dir = "sv/$family"

    produce(out_vcf) {
        exec """
            set -o pipefail

            jasmine 
                threads=$threads 
                out_dir=$output.dir 
                genome_file=$REF 
                file_list=$input.vcf.listings.txt 
                bam_list=$input.bam.listings.txt 
                out_file=$output.prefix
                min_support=$jasmine_sv.min_support 
                --mark_specific 
                spec_reads=$jasmine_sv.spec_reads 
                spec_len=$jasmine_sv.spec_len 
                --pre_normalize 
                --output_genotypes 
                --clique_merging 
                --dup_to_ins 
                --normalize_type 
                --require_first_sample 
                --default_zero_genotype
                --run_iris 
                    iris_args=
                        min_ins_length=$jasmine_sv.iris.min_ins_length,--rerunracon,--keep_long_variants

            $BASE/scripts/vcfsort -T $TMPDIR -N $threads $output.prefix | 
            awk -f $BASE/scripts/fix_allele_seq.awk | bgzip -c - > $output.vcf.gz

            tabix -p vcf $output.vcf.gz
        """
    }
}

//Assume ONT data here
cutesv = {
    var cutesv_args : ''

    def max_cluster_bias_ins = (lrs_platform == 'hifi') ? cutesv_params.hifi_max_cluster_bias_ins : cutesv_params.ont_max_cluster_bias_ins
    def max_cluster_bias_del = (lrs_platform == 'hifi') ? cutesv_params.hifi_max_cluster_bias_del : cutesv_params.ont_max_cluster_bias_del
    def diff_ratio_merging_ins = (lrs_platform == 'hifi') ? cutesv_params.hifi_diff_ratio_merging_ins : cutesv_params.ont_diff_ratio_merging_ins
    def diff_ratio_merging_del = (lrs_platform == 'hifi') ? cutesv_params.hifi_diff_ratio_merging_del : cutesv_params.ont_diff_ratio_merging_del
    
    branch.dir = "sv/$sample"

    def cutesv_tmp_dir = new File(TMPDIR, UUID.randomUUID().toString()).path

    produce("${sample}.cutesv.raw.vcf") {
        exec """
            set -o pipefail

            mkdir -p $cutesv_tmp_dir

            cuteSV
                --threads $threads
                --sample ${sample}
                --report_readid
                --max_cluster_bias_INS $max_cluster_bias_ins
	            --diff_ratio_merging_INS $diff_ratio_merging_ins
	            --max_cluster_bias_DEL $max_cluster_bias_del
	            --diff_ratio_merging_DEL $diff_ratio_merging_del
                --min_support $cutesv_params.min_support
                --genotype
                $input.bam
                $REF
                ${branch.dir}/${sample}.cutesv.tmp.vcf
                $cutesv_tmp_dir

            bcftools view ${branch.dir}/${sample}.cutesv.tmp.vcf |
                awk -f $BASE/scripts/count_sv_support.awk |
                awk -f $BASE/scripts/fix_allele_seq.awk > ${branch.dir}/${sample}.cutesv.tmp2.vcf

            echo '##INFO=<ID=SUPPORT,Number=1,Type=Integer,Description="Number of reads supporting the structural variation">' > ${branch.dir}/vcf_header.txt

            bcftools annotate -h ${branch.dir}/vcf_header.txt -o $output.vcf ${branch.dir}/${sample}.cutesv.tmp2.vcf

            rm ${branch.dir}/${sample}.cutesv.tmp*.vcf ${branch.dir}/vcf_header.txt

            rm -rf $cutesv_tmp_dir
        """
    }
}

filter_sv_calls = {
//    input:
//        file vcf
//        tuple path(mosdepth_bed), path(mosdepth_dist), path(mosdepth_threshold) // MOSDEPTH_TUPLE
//        file target_bed
//    output:
//        path "*.filtered.vcf", emit: vcf
//    script:
//        def sv_types_joined = params.sv_types.split(',').join(" ")
    var sv_tool: 'sniffles'

    def in_file = ''
    def out_files = []
    
    if (sv_tool == 'sniffles') {
        in_file = "${sample}.sniffles.vcf"
        out_files = ["${sample}.wf_sv.vcf.gz", "${sample}_filter.sh"]
    }
    else {
        in_file = "${sample}.${sv_tool}.raw.vcf"
        out_files = ["${sample}.${sv_tool}.vcf.gz", "${sample}_${sv_tool}_filter.sh"]
    }

    from(in_file) produce(out_files) {
        exec """
            set -o pipefail

            $BASE/scripts/get_filter_calls_command.py 
                --target_bedfile $opts.targets
                --vcf $input.vcf
                --depth_bedfile $input.bed.gz
                --min_sv_length $calling.min_sv_length
                --max_sv_length $calling.max_sv_length
                --sv_types ${calling.sv_types.split(',').join(' ')}
                --min_read_support $calling.min_read_support
                --min_read_support_limit $calling.min_read_support_limit > $output.sh

            bash $output.sh > $output.vcf.gz.prefix

            $BASE/scripts/vcfsort -T $TMPDIR -N $threads $output.vcf.gz.prefix | bgzip -c > $output.vcf.gz

            tabix -p vcf $output.vcf.gz
        """
        
        sample_sv_vcfs.get("${sv_tool}#${sample}", []).add(output.vcf.gz.prefix.toString())
    }
}

// NOTE This is the last touch the VCF has as part of the workflow,
//  we'll rename it with its desired output name here
sort_sv_vcf = {
//    label "wf_human_sv"
//    cpus 1
//    input:
//        file vcf
//    output:
//        path "*.wf_sv.vcf", emit: vcf
    
    produce("${sample}.wf_sv.vcf.gz") {
        exec """
            vcfsort $input.vcf | bgzip -c > $output.vcf.gz

            tabix -p vcf $output.vcf.gz
        """
    }
}

ximmer_summarize = {
    
    if(!file("$REF_BASE/decipher_population_cnvs.txt.gz").exists())
        fail "Please install the decipher population CNVs in $REF_BASE/decipher_population_cnvs.txt.gz"
    
    if(!file("$REF_BASE/dgvMerged.txt.gz").exists())
        fail "Please install the DGV population CNVs in $REF_BASE/dgvMerged.txt.gz"
 
    output.dir = sample
    
    var EXCLUDE_CNV_REGIONS : "$REF_BASE/cnv_tel_cen_excluded_regions_v2.bed",
        target_bed : opts.targets
    
    if(!file(EXCLUDE_CNV_REGIONS))
        fail "Please set the excluded CNV regions in $EXCLUDE_CNV_REGIONS or override the  EXCLUDE_CNV_REGIONS variable"
     
    if(!file(GNOMAD_SV_VCF).exists())
        fail "Please set the excluded CNV regions in $EXCLUDE_CNV_REGIONS or override the  EXCLUDE_CNV_REGIONS variable"
 
    requires GNOMAD_SV_VCF : "Please configure the location of the gnomAD-sv VCF file"
    
    produce('local_combined_cnvs.json', 'local_cnv_report.tsv') {
        exec """

            export JAVA_OPTS="-Xmx${memory}g"

            $tools.GROOVY -cp $XIMMER_GNGS_JAR:$tools.XIMMER/src/main/groovy:$tools.XIMMER/src/main/resources:$tools.XIMMER/src/main/js $tools.XIMMER/src/main/groovy/SummarizeCNVs.groovy  
                    -ddd $REF_BASE/decipher_population_cnvs.txt.gz
                    -dgv $REF_BASE/dgvMerged.txt.gz  
                    -refgene $REF_BASE/refGene.txt.gz  
                    -target $target_bed  
                    -sniffle $input.vcf.gz
                    -o $output.dir/cnv_report.html  
                    -gmd $GNOMAD_SV_VCF
                    -o $output.dir/cnv_report.html   
                    -x50 $EXCLUDE_CNV_REGIONS
                    -json $output.dir/local_combined_cnvs.json  
                    -tsv $output.dir/local_cnv_report.tsv 



        """, "wgs_summarize_cnvs"
    }
}

format_qc = {
    
    doc "Format the coverage depth output from mosdepth as Ximmer compatible JSON / JS file"
    
    output.dir="qc"
    
    produce("${sample}.cov.js") {
    
        groovy """

            def cov_data = new graxxia.TSV("$input.mosdepth.summary.txt").toListMap().find { it.chrom == 'total_region'}

            new File("$output.js").text = "covs = // NOJSON\\n" + groovy.json.JsonOutput.toJson([
                "means": [
                    "$sample": cov_data.mean
                ],
                "medians": [
                    "$sample": cov_data.mean
                ]
            ])

        """, "local"
    } 
}

zip_summary = {

    doc "Create zip files ready to load into CXP"

    output.dir = sample

    produce("${sample}.zip", "${sample}.qc.zip") {
        exec """
            zip $output.zip $input.json $input.tsv $sample/cnv_calls.js 

            mkdir -p $sample/qc

            cp -v $input.cov.js $sample/qc

            cd $sample; zip -r ${file(output2.zip).name} qc

        """     
    }
}


post_to_cxp = {
    
    var POST_RESULTS : false,
        cnv_batch: new File('.').absoluteFile.parentFile.parentFile.name,
        CXP_PROJECT : 'R0001_residual_project',
        XIMMER_GNGS_JAR : "$tools.XIMMER/tools/groovy-ngs-utils/1.0.9/groovy-ngs-utils.jar",
        STAGE_CNV_RESULTS_TARGET : null
    
    requires CXP_URL : 'URL of CXP server to POST to'

    output.dir = sample
    
    if(!POST_RESULTS) 
        succeed "Skipping posting results for ${sample} because the POST_RESULTS flag is false. Set this to true to send results to CXGo."

    uses(cxp_posts:1) {
        from(input.bam, "${sample}.zip", "${sample}.qc.zip") {
            produce("${sample}.post.log.txt") {
                
                def stageCommand = "";
                def IMPORT_ANALYSIS_ZIP="$input1.zip"
                def IMPORT_QC_ZIP="$input.qc.zip"
                
                if(STAGE_CNV_RESULTS_TARGET) {
                    stageCommand = 
                    """
                        scp $input1.zip $STAGE_CNV_RESULTS_TARGET

                        scp $input.qc.zip $STAGE_CNV_RESULTS_TARGET
                    """
                    
                    // Strip content prior to the : to allow staging with scp to remote host,
                    // but plain path to be passed in the POST command
                    IMPORT_ANALYSIS_ZIP=STAGE_CNV_RESULTS_TARGET.replaceAll('^.*:','') + "/" + file(input1.zip).name
                    IMPORT_QC_ZIP=STAGE_CNV_RESULTS_TARGET.replaceAll('^.*:','') + "/" + file(input.qc.zip).name
                }

                var sampleSex : sex // meta.get(sample).sex
                 
                exec """
                    $stageCommand

                    set -o pipefail

                    $tools.GROOVY -cp $XIMMER_GNGS_JAR:$tools.XIMMER/src/main/groovy:$tools.XIMMER/src/main/resources:$tools.XIMMER/src/main/js $tools.XIMMER/src/main/groovy/PostToCXPWGS.groovy
                        -project $CXP_PROJECT
                        -analysis $IMPORT_ANALYSIS_ZIP
                        -sex $sample:$sampleSex
                        -cxp $CXP_URL 
                        -target $target_bed
                        -batch $cnv_batch
                        -bam $sample:$input.bam
                        -qc $IMPORT_QC_ZIP  | tee $output.txt


                """
            }
        }
    }
}

symbolic_alt = {
    
    var XIMMER_GNGS_JAR : "$tools.XIMMER/tools/groovy-ngs-utils/1.0.9/groovy-ngs-utils.jar"

    branch.sv_out = branch.family_branch ? branch.family_id : branch.sample
    
    output.dir = "sv/${branch.sv_out}"

    transform(".vcf.gz") to(".vcf") {
        exec """
            gunzip -c $input.vcf.gz |
            $tools.GROOVY -cp $XIMMER_GNGS_JAR -e 'gngs.VCF.filter() { it.update { v -> { if(v.info.SVTYPE in ["INS", "DEL"]) { v.alt = "<" + v.info.SVTYPE + ">" } } } }'
            > $output.vcf
        """
    }
}

sv_annotate = {

    output.dir = "sv/${branch.sv_out}"

    exec """
        gatk SVAnnotate 
            -V $input.vcf
            --protein-coding-gtf $GENCODE_GTF
            --lenient
            -O $output.vcf.gz
    """
}

strvctvre_annotate = {

    doc "Run StrVCTVRE"
    
    requires PHYLOP100WAY : "Please set the location of PHYLOP100WAY"

    output.dir = "sv/${branch.sv_out}"

    transform('vcf.gz') to('strvctvre.vcf.bgz') {
        exec """
            python $tools.STRVCTVRE_HOME/StrVCTVRE.py -p $PHYLOP100WAY -i "$input.vcf.gz" -o $output.vcf.bgz.prefix

            bgzip -c $output.vcf.bgz.prefix > $output.vcf.bgz

            tabix -p vcf $output.vcf.bgz
        """
    }
}


prepare_sv_alignment_plots = {

    output.dir = "results/plots"

    produce("${sample}.zip") {
        exec """
            set -uo pipefail

            python $tools.LRS_PLOTTING_HOME/prepare-alignments.py
                $input.bam
                $input.vcf.gz
                $output.zip
        """
    }
}
