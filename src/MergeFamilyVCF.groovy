import gngs.SAM
import gngs.ToolBase
import gngs.VCFSiteWalker
import groovy.transform.CompileStatic
import groovy.util.logging.Log
import htsjdk.variant.variantcontext.VariantContext
import htsjdk.variant.variantcontext.VariantContextBuilder
import htsjdk.variant.vcf.VCFEncoder
import htsjdk.variant.vcf.VCFFileReader
import htsjdk.variant.vcf.VCFHeader
import htsjdk.variant.variantcontext.Allele
import htsjdk.variant.variantcontext.Genotype
import htsjdk.variant.variantcontext.GenotypeBuilder
import javax.annotation.Nonnull
import javax.annotation.Nullable

/**
 * Pseudo Joint Genotyping tool that combines samples from multiple VCFs, inferring
 * hom-ref genotypes with heuristic based values for DP, AD and GQ for inferred variants.
 */
@Log
class MergeFamilyVCF extends ToolBase {
    
    Map<String, Integer> genotypeIndexes
    
    List<String> samples
    
    List<String> vcfs

    List<SAM> bams
    
    private VCFHeader header

    private Writer output

    private VCFEncoder encoder

    List<VariantContext> lastGenotypes

    @Override
    public void run() {
        
        if(!opts.arguments()) {
            throw new IllegalArgumentException("Please provide one or more VCFs to combine")
        }
        
        this.vcfs = opts.arguments().grep { it.endsWith('.vcf.gz') }
        if(!this.vcfs)
            throw new IllegalArgumentException("Please provide one or more VCFs to combine")

//        this.bams = opts.arguments().grep { it.endsWith('.bam') }.collect { new SAM(it) }
//        if(!this.bams || this.bams.size() != this.vcfs.size())
//            throw new IllegalArgumentException("Please provide a BAM file for each VCF")

        log.info "Merging samples from VCFs ${vcfs}"
        
        this.samples = opts.arguments().collect { new gngs.VCF(it).samples[0] }
        
        log.info "Processing samples $samples"
        
        initializeOutput()
        createHeader()
        writeHeader()
        this.encoder = new VCFEncoder(header, true, true)
        
        this.lastGenotypes = this.samples.collect { null }
        
        def walker = new VCFSiteWalker(opts.arguments())
        walker.walkByAllele { List<VariantContext> ctxs ->
            VariantContext outputVariant = synthesizeOutputVariant(ctxs)
            output.write(this.encoder.encode(outputVariant))
            output.write('\n')
        }
        
        output.flush()
        output.close()
    }
    
    
    
    /**
     * Construct a variant context containing all the genotypes of the input samples, which are
     * 
     * - derived from actual genotype at site if it is present
     * - inferred as hom-ref if it is not
     * 
     * Hom-ref inferred variants are assigned GQ, AD and DP values derived from the
     * last variant for the sample at the site. This is assuming that the VCF is written
     * with gVCF-like properties where blocked values are in use and the last value observed
     * is representative until a different value is included.
     * 
     * @param   presentCtxs     List of contexts that are actually present at the site for the given allele
     * @return  VariantContext  representing combined genotypes for all samples
     */
    @CompileStatic
    @Nonnull
    VariantContext synthesizeOutputVariant(List<VariantContext> presentCtxs) {
        
        List<VariantContext>  ctxs = this.samples.collect { String sample ->
            presentCtxs.find {
                 sample in it.sampleNames
             }
        }
        
        // Update the observed genotype for each sample that has a genotype at this allele
        ctxs.eachWithIndex { VariantContext ctx, int index ->
            if(ctx != null) {
                lastGenotypes[index] = ctx
            }
        }

        // Use the first non-null variant context with an alternative allele
        VariantContext base = ctxs.find  { it != null && it.isVariant() }
        
        // If all samples are hom-ref, just choose the first
        if(!base)
            base = ctxs.find  { it != null }
            
        Allele ref = base.reference
        Allele alt = base.alternateAlleles[0] ?: ref

        List<Allele> baseAlleles = (List<Allele>)presentCtxs*.alleles.flatten().unique { a -> ((Allele)a).baseString }

        int i = 0
        List<Genotype> genotypes = ctxs.collect { VariantContext ctx ->
            Genotype result = createSampleGenotype(baseAlleles, samples[i], ctx?.genotypes?.getAt(0), lastGenotypes[i])
            ++i
            return result
        }
        
        // log.info "Process $ref/$alt with genotypes : " + genotypes.join(", ")
        double maxQual = 1.0d
   
        VariantContext combined =
                new VariantContextBuilder(base)
                .alleles([ref,alt] as Set) // groovy calls wrong method if List used directly
                .genotypes(genotypes)
                .log10PError(maxQual/-10d)
                .make()
        return combined
    }
    
    /**
     * Return an output genotype for the given sample input genotype. If the input genotype
     * is complete then it returns a genotype comprised of the primary alternate allele at the site. If it is not,
     * the genotype is inferred as hom-ref and realistic values are heuristically assigned to the DP, AD and GQ fields.
     * 
     * @param preferredAlleles  the allele(s) to use as the primary alternate allele, if one of them is present
     * @param sampleId          the sample id to create the genotype for
     * @param sampleGenotype
     * @param lastVariant
     * @return  output genotype
     */
    @CompileStatic
    Genotype createSampleGenotype(final List<Allele> preferredAlleles, final String sampleId, final Genotype sampleGenotype, final VariantContext lastVariant) {
        if(sampleGenotype == null) {
            Allele ref = preferredAlleles.find { Allele a -> a.reference }
            def gtb = new GenotypeBuilder(sampleId)
                .alleles([ref, ref])
                
            if(lastVariant) {
                gtb.DP(lastVariant.genotypes.get(0).DP)
                gtb.GQ(lastVariant.genotypes.get(0).GQ)
            }

            return gtb.make()
        }
        
        List<Allele> alleles = sampleGenotype.alleles.collect { Allele a ->
            Allele baseAllele = preferredAlleles.find { a.baseString == it.baseString }
            if(baseAllele != null)
                return baseAllele
            else
                return a
        }
        return new GenotypeBuilder(sampleGenotype)
            .alleles(alleles)
            .make()
    }
     
    private void initializeOutput() {
        if(opts.o) {
            output = new File(opts.o).newWriter()
        }
        else
            output = new PrintWriter(System.out)
    }

    /**
     * Create the header for the output VCF from the first input VCF
     */
    private void createHeader() {
        VCFFileReader vcfReader = new VCFFileReader(new File(opts.arguments()[0]))
        header = new VCFHeader(vcfReader.fileHeader.metaDataInInputOrder, samples)
    }

    @CompileStatic
    private void writeHeader() {
        output.write(header.metaDataInInputOrder.collect { '##' + it }.join('\n') + '\n')
        output.write((['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT'] + samples).join('\t') + '\n')
    }
    
    static void main(String [] args) {
        cli('MergeFamilyVCF <vcfs> <bams>', 'Merges VCF for a family including ref-filling for variants not present in parents', args) {
            o 'Output file (default: stdout)', args:1, required: false
        }
    }
}
