cwlVersion: v1.2
class: CommandLineTool
id: deepsomatic_postprocess_variants
doc: "Deepsomatic Postprocess Variants"
requirements:
  - class: InlineJavascriptRequirement
  - class: ShellCommandRequirement
  - class: ResourceRequirement
    ramMin: $(inputs.ram*1000)
    coresMin: $(inputs.cpu)
  - class: DockerRequirement
    dockerPull: 'google/deepsomatic:1.7.0-gpu'
  - class: InitialWorkDirRequirement
    listing: [$(inputs.infile), $(inputs.nonvariant_site_tfrecords), $(inputs.small_model_cnv_records)]
baseCommand: []
arguments:
  - position: 0
    shellQuote: false
    valueFrom: >-
      postprocess_variants
inputs:
  gvcf_outfile: { type: 'string?', inputBinding: { position: 2, prefix: "--gvcf_outfile"}, doc: "Optional. Destination path where we will write the Genomic VCF output." }
  infile: { type: 'File[]', doc: "Required. Path(s) to CallVariantOutput protos in TFRecord format to postprocess. These should be the complete set of outputs for call_variants.py." }
  infile_path: { type: 'string', inputBinding: { position: 2, prefix: "--infile" }, doc: "Smart path through which all CallVariantOutput protos in TFRecord format to postprocess." }
  nonvariant_site_tfrecords: { type: 'File[]?', doc: "Optional. Path(s) to the non-variant sites protos in TFRecord format to convert to gVCF file. This should be the complete set of outputs from the --gvcf flag of make_examples.py." }
  nonvariant_site_tfrecord_path: { type: 'string?', inputBinding: { position: 2, prefix: "--nonvariant_site_tfrecord_path", shellQuote: false }, doc: "Smart path through which all nonvariant_site_tfrecords files can be found." } 
  outfile: { type: 'string', inputBinding: { position: 2, prefix: "--outfile"}, doc: "Required. Destination path where we will write output variant calls in VCF format." }
  par_regions_bed: { type: 'File?', inputBinding: { position: 2, prefix: "--par_regions_bed"}, doc: "Optional BED file containing Human Pseudoautosomal Region (PAR) regions.Variants within this region are unaffected by genotype reallocation applied on regions supplied by --haploid_contigs flag." }
  pon_filtering: { type: 'string?', inputBinding: { position: 2, prefix: "--pon_filtering"}, doc: "Optional. Path to VCF file (in the docker) with Panel of Normals (PON) data.If set, the output VCF will be filtered: any variants that appear in PON will be marked with a PON filter, and PASS filter value will be removed." }
  pon_filtering_custom: { type: 'File?', inputBinding: { position: 2, prefix: "--pon_filtering"}, doc: "Optional. Custom VCF with Panel of Normals (PON) data.If set, the output VCF will be filtered: any variants that appear in PON will be marked with a PON filter, and PASS filter value will be removed." }
  ref: { type: 'File', secondaryFiles: [{ pattern: '.fai', required: true}], inputBinding: { position: 2, prefix: "--ref"}, doc: "Required. Genome reference in FAI-indexed FASTA format. Used to determine the sort order for the emitted variants and the VCF header." }
  regions_bed: { type: 'File[]?', inputBinding: { position: 2, prefix: "--regions", itemSeparator: "," }, doc: "Optional. Space-separated list of regions we want to process. Elements can be region literals (e.g., chr20:10-20) or paths to BED/BEDPE files. This should match the flag passed to make_examples.py." }
  sample_name: { type: 'string?', inputBinding: { position: 2, prefix: "--sample_name" }, doc: "Optional. If set, this will only be used if the sample name cannot be determined from the CallVariantsOutput or non-variant sites protos." }
  small_model_cnv_records: { type: 'File[]?', doc: "Optional. Path(s) to CallVariantOutput protos in TFRecord format that were called by the small model to include in postprocess." }
  small_model_cnv_records_path: { type: 'string?', inputBinding: { position: 2, prefix: "--small_model_cnv_records" }, doc: "Smart path through which all small_model_cnv_records files can be found." }

  # Catchall for Remaining Options
  extra_args: { type: 'string?', inputBinding: { position: 2, shellQuote: false }, doc: "Any additional args for deepsomatic make_call_variants" }

  cpu:
    type: 'int?'
    default: 1
    doc: "Number of CPUs to allocate to this task."
  ram:
    type: 'int?'
    default: 16
    doc: "Maximum GB of RAM to allocate for this tool."
outputs:
  output: 
    type: File
    secondaryFiles: [{pattern: ".tbi", required: true}]
    outputBinding:
      glob: $(inputs.outfile)
  gvcf:
    type: File?
    secondaryFiles: [{pattern: ".tbi", required: true}]
    outputBinding:
      glob: $(inputs.gvcf_outfile)
