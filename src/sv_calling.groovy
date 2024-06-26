

filterBam = {

    output.dir = "align"

    filter('filtered') {
        exec """
            samtools view -@ $threads 
                $input.bam
                -F 2308 
                -o $output.bam
                --write-index 
                --reference $REF
        """
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
                --output-rnames
                --cluster-merge-pos $calling.cluster_merge_pos
                --input $input.bam
                --tandem-repeats ${calling.tr_bed} $sniffles_args
                --vcf ${output.vcf.prefix}.tmp.vcf

            sed '/.:0:0:0:NULL/d' ${output.vcf.prefix}.tmp.vcf > $output.vcf
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
    
    from(family_snfs*.value.flatten()) produce("${family}.family.sv.vcf.gz") {
        exec """
            sniffles
                --threads $threads
                --input $inputs
                --vcf $output.prefix

            bgzip -c $output.prefix > $output.vcf.gz

            bcftools index -t $output.vcf.gz
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
    

    from('*.regions.bed.gz') produce("${sample}.wf_sv.vcf.gz", "${sample}_filter.sh") {
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

            $BASE/scripts/vcfsort -T $TMPDIR -N $threads $output.vcf.gz.prefix | bgziptabix $output.vcf.gz
        """
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
            vcfsort $input.vcf | bgziptabix $output.vcf.gz
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

            $tools.GROOVY -cp $XIMMER_GNGS_JAR:$XIMMER/src/main/groovy:$XIMMER/src/main/resources:$XIMMER/src/main/js $XIMMER/src/main/groovy/SummarizeCNVs.groovy  
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
        XIMMER_GNGS_JAR : "$XIMMER/tools/groovy-ngs-utils/1.0.9/groovy-ngs-utils.jar",
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

                var sampleSex = meta.get(sample).sex
                 
                exec """
                    $stageCommand

                    set -o pipefail

                    $tools.GROOVY -cp $XIMMER_GNGS_JAR:$XIMMER/src/main/groovy:$XIMMER/src/main/resources:$XIMMER/src/main/js $XIMMER/src/main/groovy/PostToCXPWGS.groovy
                        -project $CXP_PROJECT
                        -analysis $IMPORT_ANALYSIS_ZIP
                        -sex $sample:$sampleSex
                        -cxp $CXP_URL 
                        -target $target_bed
                        -batch $cnv_batch
                        -bam $sample:$input.bam
                        -qc $IMPORT_QC_ZIP 
                """
            }
        }
    }
}

symbolic_alt = {
    
    var XIMMER_GNGS_JAR : "$tools.XIMMER/tools/groovy-ngs-utils/1.0.9/groovy-ngs-utils.jar"

    branch.sv_out = branch.family_branch ? branch.name : branch.sample
    
    output.dir = "sv/${branch.sv_out}"

    transform(".vcf.gz") to(".vcf") {
        exec """
            gunzip -c $input.vcf.gz | 
            $tools.GROOVY -cp $XIMMER_GNGS_JAR -e 'gngs.VCF.filter() { it.update { v -> if(v.info.SVTYPE=="INS") v.alt = "<INS>" }}'
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
