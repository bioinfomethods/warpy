/**
 * Read length QC plots
 *
 * Creates some simple plots of the read length distribution and calculates an estimate
 * of the N50 value
 */
import gngs.*
import graxxia.*
import gngs.plot.*
import gngs.plot.bx.*

cli = new CliBuilder(usage: 'read_lengths -bam <bam file>').tap {
    bam 'BAM file to calculate statistics from', args:1, required: true
    prefix 'Prefix for output files', args:1, required: true
    bed 'Regions to calculate read lengths for', args:1, required: false
}

opts = cli.parse(args)
if(!opts)  {
	System.exit(1)
}

log = java.util.logging.Logger.getLogger('read_lengths')
Utils.configureSimpleLogging()

bam = new SAM(opts.bam)

if(opts.bed) {
    regions = new BED(opts.bed).load()
}
else {
    // Everything is computed from a subset of the genome to make it fast
    // Pick completely arbitrary set of regions to select read lengths from
    regions = [
        new Region('chr1:20000000-24000000'),
        new Region('chr2:40000000-44000000'),
        new Region('chr9:40000000-44000000'),
        new Region('chr15:70000000-74000000'),
    ] as Regions
}

lens = new ArrayList(100000)

regions.each { region ->
	bam.withIterator(region) { i ->
		i.each { r ->
			if(r.secondaryOrSupplementary)
				return
			lens.add(r.readLength)        
		}
	}
}

// Raw read length distribution
p = new Plot(title: 'Read Length Distribution', xLabel: 'Read Length', yLabel: 'Frequency') << 
	new Density.Area(data: lens.grep { it < 30000 } )

readLengthOutputFileName = opts.prefix + '.lengths.png'
p.save(readLengthOutputFileName)
log.info "Wrote $readLengthOutputFileName"

// Bin the reads by length so we can calculate weighted distribution
binner = new Binner(200, 0, 50000) // 200 bins from 0 to 50kb
len_counts = lens.countBy { binner.bin(it)}.collect {
    [
        bin: it.key,
        length: binner.midPoints[it.key],
        count: it.value
    ]
}
.grep { it.length } // fell into a bin
.each {
    it.bases = it.length * it.count
}
.sort { it.length }


// Read lengths weighed by bases
p = new gngs.plot.Plot(title: 'Read Length by Sequenced Mb (subset)', xLabel: 'Read Length', yLabel: 'Sequenced Bases') << \
new gngs.plot.Bars(x: len_counts*.length, y: len_counts*.bases.collect { it / 1000000 }, width: 1000)

readLengthByGBOutputFileName = opts.prefix + '.lengthsgb.png'
p.save(readLengthByGBOutputFileName)
log.info "Wrote $readLengthByGBOutputFileName"

// Calculate the n50
total_gb = len_counts.sum { it.bases }
cum_n50 = 0
n50_table = len_counts.sort { -it.length }.takeWhile { (cum_n50 += it.bases) < (total_gb/2) }
n50 = n50_table[-1].length

n50OutputFileName = opts.prefix + '.n50.txt'
new File(n50OutputFileName).text = String.valueOf(n50) + '\n'
log.info "Wrote $n50OutputFileName"

log.info "Done"

