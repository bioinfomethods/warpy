
scriptFile = new File(bpipe.Config.config.script)

if(scriptFile.name == 'qc_reports.groovy') {
    options {
        design 'name of the design containing target regions to calculate coverage for', args:1, required: true
    }

    DESIGN=opts.design
    TARGET_BED="$BASE/designs/$DESIGN/${DESIGN}.bed"
}

calc_coverage = {

    requires TARGET_BED : 'Target regions to calculate QC for'

    var GAP_THRESHOLD : 6


    def GAP_TARGET = null
    def customGapTarget = TARGET_BED.absolutePath.replaceAll('.bed', '_GAP.bed')
    if(file(customGapTarget).exists()) {
        GAP_TARGET = customGapTarget
    }
    else {
        GAP_TARGET = TARGET_BED
    }

    output.dir = "qc"

    produce("${sample}.cov.stats.tsv", "${sample}.intervalsummary.tsv", "${sample}.gaps.tsv.gz") {
    exec """
         $tools.GNGS/bin/gngstool Cov  -L $TARGET_BED
                           -samplesummary $output.stats.tsv 
                           -intervalsummary $output.intervalsummary.tsv
                           -gaps $output.gaps.tsv.gz
                           -gaptarget $GAP_TARGET
                           -gt $GAP_THRESHOLD
                           -refgene $REF_BASE/refGene.txt.gz
                           $input.bam
        """
    }
}

samtools_stats = {
    output.dir = "qc"

    produce(sample + '.samtools.stats.tsv') {
        exec """
        $tools.SAMTOOLS stats --threads $threads $input.bam > $output.stats.tsv
        """
    }
}

read_lengths = {

    output.dir = "qc"

    var sample : file(input.bam).name.tokenize('.')[0]

    produce("${sample}.lengths.png", "${sample}.lengthsgb.png", "${sample}.n50.txt") {
		exec """
			$tools.GROOVY -cp $tools.GNGS/groovy-ngs-utils.jar $BASE/scripts/read_lengths.groovy
				-bam $input.bam
				-prefix $output.prefix.prefix
		"""
    }
}

send_report = {
    def sample = file(input.bam).name.tokenize('.')[0]

    def n50 = new File(input.n50.txt).text.trim()
    def stats = new graxxia.TSV(input.stats.tsv.toString()).toListMap()
    stats[0].n50 = n50

    def statsTable = gngs.Utils.table(stats[0].collect { [Metric: it.key, Value: it.value] }, out:new StringWriter())    

    send issue(
	    title: [ search: sample, match: sample ],
        channel: "qc_report_to_gitlab",
        description: "## Coverage Statistics for $sample\n\n$statsTable"
    ) to channel: 'qc_report_to_gitlab'

    send issue(
	    title: [ search: sample, match: sample ],
        channel: "qc_report_to_gitlab",
        description: "## Read Length Distribution\n![QC Plot]($input.lengths.png)\n"
    ) to channel: 'qc_report_to_gitlab', file: input.lengths.png

    send issue(
	    title: [ search: sample, match: sample ],
        channel: "qc_report_to_gitlab",
        description: "## Length-Weighted Read Length Distribution\n![QC Plot]($input.lengthsgb.png)\n"
    ) to channel: 'qc_report_to_gitlab', file: input.lengthsgb.png
}

if(scriptFile.name == 'qc_reports.groovy') {
    run {
        '%.bam' * [ 
            calc_coverage,
            read_lengths,
            samtools_stats
         ] + send_report
    }
}
