import gngs.*
import groovy.transform.CompileStatic

/**
 * Takes a VCF file and calculates per sample the fraction of heterozygous SNPs
 * that have an allele frequency outside the bounds of the expected frequency. 
 * Will then throw a warning if the fraction exceeds the given threshold. 
 * 
 * @author Simon Sadedin & Macgregor Todd
 */

class CalculateVAFOutliers extends ToolBase {

    double maxVAFOutlierFrac = 0.05
    double filterThreshold = 0.25
    double lowerBound = 0.75
    double upperBound = 0.975

    @Override
	@CompileStatic
	public void run() {
		
		if(!opts.arguments()) 
			throw new IllegalArgumentException('Please provide a VCF to calculate VAF outliers on as an argument')
		
		String vcfPath = opts.arguments()[0]
		List<String> samples = (List<String>)(opts['s'] ? [opts['s']] : new VCF(vcfPath).samples)
	
		if (opts['f']) {
			maxVAFOutlierFrac = opts['f'] as Double
		}
		if (opts['l']) {
			lowerBound = opts['l'] as Double
		}
		if (opts['u']) {
			upperBound = opts['u'] as Double
		}
        
        boolean fail = false
	        
        for(String sample : samples) {
            double outlierFrac = calculateOutlierFraction(vcfPath, sample)

            if(opts['q']) {
                println outlierFrac
            } 
            else {
                println "-~" * 40
                println("| " + ("VAF outlier fraction for $sample is " + outlierFrac).center(76) + " |")
                println "-~" * 40
            }

            if(outlierFrac > maxVAFOutlierFrac) {
                fail = true
                if(!opts['q']) {
                    println "ERROR: Sample $sample outlier fraction exceeds " + maxVAFOutlierFrac
                }
            }

        }
        
        if(fail) {
            System.exit(1)
        }
	}
    
    double calculateOutlierFraction(String vcfPath, String sample) {
	
		final int sampleIndex = new VCF(vcfPath).samples.indexOf(sample) 
		int total = 0
		int outliers = 0
        int filtered = 0
		
		VCF.parse(vcfPath) { Variant v ->
		
			double vaf = v.getVaf(1,sampleIndex)
            
            if(v.type != 'SNP')
                return false

            if (vaf < filterThreshold) {
                ++filtered
                return false
            }
            
			if ( vaf > lowerBound && vaf < upperBound ) {
				++outliers
			}
			
			++total
			
			return false
		}
		
		return total > 0 ? (outliers / total) : 0
    }
	
	static void main(String[] args) {
		cli('CalculateVAFOutlier -s <sample> -f <fraction> -l <lower bound> -u <upper bound> <vcf>', args) {
			s 'Sample to calculate VAF outliers for (if unspecified, calculate for all)', args:1, required: false
			f 'Fraction of outliers that will throw a warning', args:1, required: false
            q 'Quiet mode - print VAF outlier fraction and exit'
			l 'Lower bound of expected VAF', args:1, required: false
			u 'Upper bound of expected VAF', args:1, required: false
		}
	}
}
