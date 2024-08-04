

XIMMER_GNGS_JAR="$tools.XIMMER/tools/groovy-ngs-utils/1.0.9/groovy-ngs-utils.jar"
XIMMER_CLASSPATH="$XIMMER_GNGS_JAR:$tools.XIMMER/src/main/groovy:$tools.XIMMER/src/main/resources:$tools.XIMMER/src/main/js"

load 'sv_calling.groovy'

options  {
    targets 'Target regions to call variants in', args:1, type: File, required: true
    sex 'Sex of sample', args:1, required: true
}

inputs 'vcf.gz' : 'Sniffles VCF file to load'

target_bed = opts.targets

sample = gngs.VCF.parse(args[0]).samples[0]

sex = opts.sex


run {
    ximmer_summarize + format_qc + zip_summary + post_to_cxp
}
