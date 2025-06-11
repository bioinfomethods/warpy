options {
    design 'name of the design containing target regions to calculate coverage for', args:1, required: true
}

calc_coverage = {
    output.dir = "qc"
    
   transform('.bam') to('.cov.stats.tsv', '.intervalsummary.tsv') {
		exec """
			 $tools.GNGS/bin/gngstool Cov  -L $BASE/designs/$opts.design/${opts.design}.bed 
						   -samplesummary $output.stats.tsv 
						   -intervalsummary $output.intervalsummary.tsv $input.bam
		"""
    }
}

read_lengths = {

    output.dir = "qc"

    def sample = file(input.bam).name.tokenize('.')[0]

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

run {
    '%.bam' * [ 
		calc_coverage,
		read_lengths
	 ] + send_report
}
