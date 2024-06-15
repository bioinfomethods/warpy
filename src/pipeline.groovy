
import gngs.*
import org.yaml.snakeyaml.*

title 'Warpy - Oxford Nanopore Super High Accuracy Pipeline'

options {
    samples 'Sample metadata file', args:1, type: File, required: true
    targets 'Target regions to call variants in', args:1, type: File, required: true
    sex 'Sample sex (required for STR analysis)', args:1, type: String, required: true
}

load 'stages.groovy'
load 'sv_calling.groovy'
load 'str_calling.groovy'
load 'methylation.groovy'

meta =
 file(opts.samples)
    .withReader { r -> Collections.synchronizedMap(new Yaml().load(r)) }
        .samples
        .collectEntries { [ it.identifier,  it ] } as Map

println "Analysis samples are: " + meta*.key
    
// input_files = scanInputDirectories(args)
input_files = meta.collectEntries { sample, sampleMeta ->

   def sampleInputs = sampleMeta.inputs.collect { i ->
       if(new File(i).directory) {
           new File(i).listFiles().grep { f -> ext(f) in ['pod5', 'blow5','fast5'] }*.path
       }
       else
           return i
   }
           
   [ sample, sampleInputs.flatten()]
}

println(new groovy.json.JsonBuilder(input_files).toPrettyString())

println("Flattened: " + input_files*.value.flatten())

by_extension = input_files*.value.flatten().groupBy { ext(new File(it)) }
if(by_extension.size() > 1) 
    throw new bpipe.PipelineError("Warpy can only accept a single file type for analysis in one execution")
    
input_pattern = '%' + by_extension*.key[0]

// to make pipeline generic to work for either fast5 or blow5,
// define virtual file extentions 'x5' that can map to either
filetype x5 : ['pod5', 'blow5','fast5']

input_data_type = (by_extension.bam ? 'bam' : 'x5')

Map params = model.params

DRD_MODELS_PATH="$BASE/models/dorado"
CLAIR3_MODELS_PATH="$BASE/models/Clair3"

VERSION="1.0"

// Convert basecaller model to Clair3 model
model_map_file = "$BASE/data/clair3_models.tsv"
model_map = new graxxia.TSV(model_map_file).toListMap()
clair3_model = model_map.find { it.basecaller == 'dorado' && it.basecall_model_name == params.drd_model }
assert clair3_model != null : 'No dorado model for base caller model ' + params.drd_model + ' could be found in ' + model_map_file

if(clair3_model.clair3_model_name == '-') 
    throw new bpipe.PipelineError("No suitable clair3 model could be found: $clair3_model.clair3_nomodel_reason")

targets = bed(opts.targets)

str_chrs = new File(calling.repeats_bed).readLines()*.tokenize()*.getAt(0).unique()

println "The chromosomes for STR calling are: $str_chrs"

dorado_group_size = 10

input_groups = input_files.collectEntries { [ it.key, it.value.collate(dorado_group_size).indexed().collectEntries { [ "dorado_group_" + it.key, it.value] } ] }

println "The input groups are: \n\n"  + input_groups

calling_chunk_size = 10000000

targets_by_chr = new gngs.BED(opts.targets).load().groupBy { it.chr }.collect { chr, List<Region> regions ->
   new gngs.Region(chr, regions*.from.min(), regions*.to.max() )
}

genome 'hg38'


contigs = channel(targets_by_chr*.chr.unique()).named('chr')

// target_channel = channel(targets_by_chr).named('clair_chunk')

sample_channel = channel(input_files).named('sample')

// Set the family ids equal to sample id if no family id is provided
meta.each { if(!it.value.family_id) it.value.family_id = it.value.identifier }

family_channel = channel(meta*.value*.family_id).named('family')
   
init = {
    println "\nProcessing ${input_files.size()} input files ...\n"
    println "\nUsing base calling model: $params.drd_model"
    println "\nUsing clair3 model: $clair3_model"
    
    produce('versions.txt') {
        exec """
            echo $VERSION > versions.txt
        """
    }
    
    produce("CONTIGS") {
        groovy """
            new File('CONTIGS').text = new gngs.BED('$opts.targets').load()*.chr.unique().join('\\n') + '\\n'
        """
    }
}

basecall_align_reads = segment {
    make_mmi.when { ! new File(REF_MMI).exists() } + input_groups * [ convert_fast5_to_pod5.when { input.x5.endsWith('.fast5') } +
        dorado + minimap2_align ] + merge_pass_calls
}

forward_sample_bam = {
    
    // println "Forwarding bam file $input.bam for sample $sample"
    
    forward input.bam
}

sample_gvcfs = Collections.synchronizedMap([:])
sample_snfs = Collections.synchronizedMap([:])

run(input_files*.value.flatten()) {
    init + check_tools + 
    
    // Phase 1: resolve or create BAM files
    sample_channel * [
        basecall_align_reads.when { input_data_type == 'x5' } + forward_sample_bam.when { input_data_type == 'bam' } + read_stats 
    ] + 
    
    // Phase 2: single sample variant calling
    [
         snp_calling : sample_channel * [ 
             hg38.partition(10000000)  * [ pileup_variants ] + aggregate_pileup_variants +
             [ 
                    get_qual_filter,
                    contigs * [
                        select_het_snps + phase_contig,
                        create_candidates + '%.bed' * [ evaluate_candidates ] 
                    ] 
                    + aggregate_full_align_variants
             ] +
                contigs * [ merge_pileup_and_full_vars ] + aggregate_all_variants,
         ],

         sv_calling: sample_channel * [ mosdepth + filterBam + [
            sniffles2_for_trios,
            sniffles2 + filter_sv_calls
         ] ],

         methylation: sample_channel * [ bam2bedmethyl ],
         
         str_calling: sample_channel * [ chr(*str_chrs) * [ call_str + annotate_repeat_expansions ] + merge_str_tsv + merge_str_vcf ]
    ] +

    // Phase 3: family merging
    family_channel * [ sniffles2_joint_call, (combine_family_gvcfs + genotype_gvcfs).when { calling.enable_gvcf } ] +
    
    // annotate
    symbolic_alt + sv_annotate + strvctvre_annotate
}
