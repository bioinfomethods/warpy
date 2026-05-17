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
    tsv 'Output file for weighted read length distribution (TSV, default: <prefix>.lengths.tsv)', args:1, required: false
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

log.info "Sampled ${lens.size()} length values from ${regions.numberOfRanges} regions"

// Raw read length distribution
p = new Plot(title: 'Read Length Distribution', xLabel: 'Read Length', yLabel: 'Frequency') << 
	new Density.Area(data: lens.grep { it < 30000 } )

readLengthOutputFileName = opts.prefix + '.lengths.png'
p.save(readLengthOutputFileName)
log.info "Wrote $readLengthOutputFileName"

def binLengths(lens, numBins, maxLen) {
    def binner = new Binner(numBins, 0, maxLen)
    def binWidth = maxLen / numBins
    return lens.countBy { binner.bin(it) }.collect {
        [
            bin: it.key,
            length: binner.midPoints[it.key],
            upperBound: (it.key + 1) * binWidth,
            count: it.value
        ]
    }
    .grep { it.length }
    .each {
        it.bases = it.length * it.count
    }
    .sort { it.length }
}

// Bin the reads by length so we can calculate weighted distribution
len_counts = binLengths(lens, 200, 50000)


// Read lengths weighed by bases
p = new gngs.plot.Plot(title: 'Read Length by Sequenced Mb (subset)', xLabel: 'Read Length', yLabel: 'Sequenced Bases') << \
new gngs.plot.Bars(x: len_counts*.length, y: len_counts*.bases.collect { it / 1000000 }, width: 1000)

readLengthByGBOutputFileName = opts.prefix + '.lengthsgb.png'
p.save(readLengthByGBOutputFileName)
log.info "Wrote $readLengthByGBOutputFileName"

def tsvOutputFileName = opts.tsv ?: opts.prefix + '.lengths.tsv'
def tsv_counts = binLengths(lens, 40, 20000)
def totalCount = tsv_counts.sum { it.count }
def totalBases = tsv_counts.sum { it.bases }
double cumulativeFraction = 0.0
double cumulativeFractionBases = 0.0
new File(tsvOutputFileName).withWriter { w ->
    w.writeLine("length\tcount\tfraction\tcumulative_fraction\tbases\tfraction_bases\tcumulative_fraction_bases")
    tsv_counts.each { entry ->
        def fraction = (entry.count / (double) totalCount).round(4)
        def fractionBases = (entry.bases / (double) totalBases).round(4)
        cumulativeFraction += fraction
        cumulativeFractionBases += fractionBases
        w.writeLine("${(int) entry.upperBound}\t${entry.count}\t${fraction}\t${cumulativeFraction.round(4)}\t${String.format('%d', (long) entry.bases)}\t${fractionBases}\t${cumulativeFractionBases.round(4)}")
    }
}
log.info "Wrote $tsvOutputFileName"


// Calculate the n50
total_gb = len_counts.sum { it.bases }
cum_n50 = 0
n50_table = len_counts.sort { -it.length }.takeWhile { (cum_n50 += it.bases) < (total_gb/2) }
n50 = n50_table[-1].length

n50OutputFileName = opts.prefix + '.n50.txt'
new File(n50OutputFileName).text = String.valueOf(n50) + '\n'
log.info "Wrote $n50OutputFileName"

log.info "Done"

