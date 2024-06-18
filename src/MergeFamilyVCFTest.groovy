import static org.junit.Assert.*

import org.junit.Test

import htsjdk.variant.variantcontext.Allele
import htsjdk.variant.variantcontext.Genotype
import htsjdk.variant.variantcontext.GenotypeBuilder
import htsjdk.variant.variantcontext.GenotypesContext
import htsjdk.variant.variantcontext.SimpleAllele
import htsjdk.variant.variantcontext.VariantContext
import htsjdk.variant.variantcontext.VariantContextBuilder
import MergeFamilyVCF

class MergeFamilyVCFTest {
    
    MergeFamilyVCF mfv = new MergeFamilyVCF(samples: ['sample1','sample2'], lastGenotypes: [null,null])
    
    /**
     * Test allele representing reference allele
     */
    Allele a = new SimpleAllele('A', true)
    
    /**
     * Test allele representing alternate allele
     */
    Allele c = new SimpleAllele('C', false)
        
    @Test
    public void 'identical het preserved in output'() {
        
        def variants = [
            createVariant('sample1', [a, c]),
            createVariant('sample2', [a, c]),
        ]

        VariantContext v = mfv.synthesizeOutputVariant(variants)
        assert v.alleles.find { it.reference }.baseString == 'A' 
        assert v.alleles.find { !it.reference }.baseString == 'C' 
        assert v.alleles.size() == 2
        
        def gt0 = v.genotypes[0]
        println(gt0)

        assert v.genotypes[0].alleles[0] == a
        assert v.genotypes[0].alleles[1] == c
        
        assert v.genotypes[1].alleles[0] == a
        assert v.genotypes[1].alleles[1] == c
        
        assert v.genotypes[1].isHet()
    }
    
    @Test
    public void 'missing sample2 homref in output'() {
        
        def variants = [
            createVariant('sample1', [a, c]),
        ]

        VariantContext v = mfv.synthesizeOutputVariant(variants)
        assert v.alleles.find { it.reference }.baseString == 'A' 
        assert v.alleles.find { !it.reference }.baseString == 'C' 
        assert v.alleles.size() == 2
        
        def gt0 = v.genotypes[0]

        assert v.genotypes[0].alleles[0] == a
        assert v.genotypes[0].alleles[1] == c
        
        assert v.genotypes[1].alleles[0] == a
        assert v.genotypes[1].alleles[1] == a
        
        assert v.genotypes[1].isHom()
    }    
    
    @Test
    public void 'missing sample1 homref in output'() {
        
        def variants = [
            createVariant('sample2', [a, c]),
        ]

        VariantContext v = mfv.synthesizeOutputVariant(variants)
        assert v.alleles.find { it.reference }.baseString == 'A' 
        assert v.alleles.find { !it.reference }.baseString == 'C' 
        assert v.alleles.size() == 2
        
        def gt0 = v.genotypes[0]

        assert v.genotypes[0].alleles[0] == a
        assert v.genotypes[0].alleles[1] == a
        
        assert v.genotypes[1].alleles[0] == a
        assert v.genotypes[1].alleles[1] == c
        
        assert v.genotypes[0].isHom()
        assert v.genotypes[1].isHet()
    }    
    
   
    private VariantContext createVariant(id,List<Allele> alleles, qual=100d) {
        Genotype resultGt = new GenotypeBuilder(id)
                .alleles(alleles)
                .make()

        VariantContext result = new VariantContextBuilder('T', 'chr1', 100, 100, [a,c,*alleles].unique())
        .genotypes(resultGt)
        .log10PError(qual/-10d)
        .make()
        return result
    }


}
