extract_snps = {
    output.dir = "cnv/snp"

    transform('wf_snp.norm.phased.vcf.gz') to('wf_snp.norm.phased.snp.vcf.gz') {
        exec """
            set -eo pipefail

            bcftools view -v snps -f PASS $input.vcf.gz | bgzip -c > $output.vcf.gz
            
            tabix -p vcf $output.vcf.gz
        """
    }
}

//mosdepth runs using SAM flag filter 1796, which is a superset of the stage filterBAM using filter 260

spectre_mosdepth = {
    var bin_size : 1000
   
    output.dir = "cnv/mosdepth/$sample"
    
    produce("${sample}.regions.bed.gz",
            "${sample}.mosdepth.global.dist.txt",
            "${sample}.mosdepth.summary.txt",
            "${sample}.thresholds.bed.gz") {

        exec """
            set -eo pipefail

            export REF_PATH=$REF

            export MOSDEPTH_PRECISION=3

            mosdepth
                -x 
                -t $threads
                -b $bin_size
                -f $REF
                --thresholds 1,10,20,30
                --no-per-base
                $output.dir/${sample}
                $input.cram
        """
    }
}

spectre = {
    var min_cnv_len : 2000
    var bin_size : 1000

    def cnv_target_chrs = targets_by_chr*.chr.join(',')

    def ref_gz = "align/ref/" + new File(REF).name + '.gz'
    
    output.dir = "cnv/spectre/${sample}"
    
    transform('wf_snp.norm.phased.snp.vcf.gz') to('spectre.vcf.gz') {
        exec """
            spectre CNVCaller
                --bin-size $bin_size
                --coverage cnv/mosdepth/$sample
                --snv $input.snp.vcf.gz
                --only-chr $cnv_target_chrs
                --sample-id ${sample}.spectre
                --output-dir $output.dir
                --reference $ref_gz
                --min-cnv-len $min_cnv_len
                --threads $threads
        """, "spectre"

        exec """
            set -eo pipefail

            bgzip -c $output.vcf.gz.prefix > $output.vcf.gz

            tabix -p vcf $output.vcf.gz
        """, "vcf_utils"
    }
}

cnvpytor = {
    var bin_size : 1000

    def cnv_target_chrs = targets_by_chr*.chr.join(' ')

    def ref_gz = "align/ref/" + new File(REF).name + '.gz'
    
    output.dir = 'cnv/cnvpytor'

    def pytor_file = "$output.dir/${sample}.pytor"

    from(".cram") produce("${sample}.cnvpytor.tsv", "${sample}.cnvpytor.vcf.gz") {
        exec """
            set -eo pipefail

            cnvpytor -root $pytor_file -chrom $cnv_target_chrs -T $ref_gz -rd $input.cram
            
            cnvpytor -root $pytor_file -ls
            
            cnvpytor -root $pytor_file -his $bin_size
            
            cnvpytor -root $pytor_file -partition $bin_size
            
            cnvpytor -root $pytor_file -genotype $bin_size

            cnvpytor -root $pytor_file -call $bin_size > $output.tsv

            echo -e "set Q0_range -1 0.5\\nset p_range 0 0.0001\\nset p_N 0 0.5\\nset size_range $bin_size inf\\nset dG_range 100000 inf\\nset print_filename $output.vcf.gz.prefix\\nprint calls\\nquit\\n" | cnvpytor -root $pytor_file -view $bin_size
        """, "cnvpytor"

        exec """
            set -eo pipefail

            bgzip -c $output.vcf.gz.prefix > $output.vcf.gz

            tabix -p vcf $output.vcf.gz
        """, "vcf_utils"
    }
}

ximmer_summarize_cnv = {
    output.dir = "cnv/ximmer/${sample}"
    
    var EXCLUDE_CNV_REGIONS : "$REF_BASE/centromeres.hg38.bed"

    requires GNOMAD_SV_VCF : "Please configure the location of the gnomAD-sv VCF file"
    
    produce('local_combined_cnvs.json', 'local_cnv_report.tsv') {
        exec """
            export JAVA_OPTS="-Xmx${memory}g"

            $tools.GROOVY -cp $XIMMER_GNGS_JAR:$tools.XIMMER/src/main/groovy:$tools.XIMMER/src/main/resources:$tools.XIMMER/src/main/js $tools.XIMMER/src/main/groovy/SummarizeCNVs.groovy
                    -ddd $REF_BASE/decipher_population_cnvs.txt.gz
                    -dgv $REF_BASE/dgvMerged.txt.gz  
                    -refgene $REF_BASE/refGene.txt.gz  
                    -target $opts.targets ${input.wf_sv.vcf.gz.optional.flag('-sniffle')} ${input.cutesv.vcf.gz.optional.flag('-cutesv')} ${input.spectre.vcf.gz.optional.flag('-spectre')} ${input.cnvpytor.vcf.gz.optional.flag('-cnvpytor')}
                    -o $output.dir/cnv_report.html
                    -x50 $EXCLUDE_CNV_REGIONS
                    -json $output.dir/local_combined_cnvs.json  
                    -tsv $output.dir/local_cnv_report.tsv 
        """
    }
}
