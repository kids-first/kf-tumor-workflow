cwlVersion: v1.2
class: CommandLineTool
id: deepsomatic_make_examples 
doc: "Deepsomatic make_examples"
requirements:
  - class: InlineJavascriptRequirement
  - class: ShellCommandRequirement
  - class: ResourceRequirement
    ramMin: $(inputs.max_memory*1000)
    coresMin: $(inputs.cpu)
  - class: DockerRequirement
    dockerPull: 'google/deepsomatic:1.7.0'
baseCommand: []
arguments:
  - position: 0
    shellQuote: false
    valueFrom: >-
      make_examples_somatic
  - position: 2
    shellQuote: false
    valueFrom: >-
      --examples make_examples_somatic.tfrecord@$(inputs.task_total).gz
  - position: 2
    shellQuote: false
    valueFrom: >-
      $(inputs.output_gvcf ? "--gvcf gvcf.tfrecord@" + inputs.task_total + ".gz" : "")
  - position: 2
    shellQuote: false
    valueFrom: >-
      $(inputs.output_candidate_positions ? "--candidate_positions candidate_positions@" + inputs.task_total : "")
inputs:
  # Main Args
  ref: { type: 'File', secondaryFiles: [{ pattern: '.fai', required: true }], inputBinding: { position: 2, prefix: "--ref"}, doc: "Required. Genome reference to use. Must have an associated FAI index as well. Supports text or gzipped references. Should match the reference used to align the BAM file provided to --reads." }
  reads_normal: { type: 'File[]?', secondaryFiles: [{ pattern: '.bai', required: false }, { pattern: '^.bai', required: false }, { pattern: '.crai', required: false }, { pattern: '^.crai', required: false }], inputBinding: { position: 2, prefix: "--reads_normal", itemSeparator: ","}, doc: "Required. Reads from the normal matched sample. Aligned, sorted, indexed BAM/CRAM file. Should be aligned to a reference genome compatible with --ref. Can provide multiple BAMs/CRAMs (comma-separated)." }
  reads_tumor: { type: 'File[]', secondaryFiles: [{ pattern: '.bai', required: false }, { pattern: '^.bai', required: false }, { pattern: '.crai', required: false }, { pattern: '^.crai', required: false }], inputBinding: { position: 2, prefix: "--reads_tumor", itemSeparator: ","}, doc: "Required. Reads from the tumor sample. Aligned, sorted, indexed BAM/CRAM file. Should be aligned to a reference genome compatible with --ref. Can provide multiple BAMs/CRAMs (comma-separated)." }
  checkpoint: { type: 'string?', inputBinding: { position: 2, prefix: "--checkpoint"}, doc: "Path to the TensorFlow model checkpoint." }
  checkpoint_custom: { type: 'File?', inputBinding: { position: 2, prefix: "--checkpoint"}, doc: "Custom TensorFlow model checkpoint." }
  task: { type: 'int?', inputBinding: { position: 2, prefix: "--task"}, doc: "Task ID of this task" }
  task_total: { type: 'int?', doc: "Total number of task shards being run." }
  mode:
    type:
      - 'null'
      - type: enum
        name: mode
        symbols: ["calling", "training", "candidate_sweep"]
    inputBinding:
      prefix: "--mode"
      position: 2
    doc: |
      Mode to run. Must be one of calling, training or candidate_sweep.
      calling - examples are prepared for inference only.
      training - examples are prepared with labels for training. 
      candidate_sweep - (advanced pre-step) - candidate positions are prepared for the subsequent run of make_examples with intervals created with equal number of candidates. NOTE: When this option is used, make_examples must be run again with the mode set to calling.
  candidate_positions: { type: 'File?', inputBinding: { position: 2, prefix: "--candidate_positions"}, doc: "Path to the binary file containing candidate positions." }
  output_candidate_positions: { type: 'boolean?', doc: "Should we write the binary file containing candidate positions?" }
  runtime_by_region: { type: 'string?', inputBinding: { position: 2, prefix: "--runtime_by_region"}, doc: "[optional] Output filename for a TSV file of runtimes and other stats by region. If examples are sharded, this should be sharded into the same number of shards as the examples." }
  # examples_outname: { type: 'string', inputBinding: { position: 2, prefix: "--examples"}, doc: "Required. Path to write tf.Example protos in TFRecord format." }
  output_gvcf: { type: 'boolean?', doc: "Should we write gVCF records in TFRecord of Variant proto format?" }
  regions_file: { type: 'File[]?', inputBinding: { position: 2, prefix: "--regions"}, doc: "Optional. Space-separated list of regions we want to process." }
  exclude_regions_file:  { type: 'File[]?', inputBinding: { position: 2, prefix: "--exclude_regions"}, doc: "Optional. List of regions we want to exclude from processing in BED/BEDPE files. Region exclusion happens after processing the --regions argument, so --region 20 --exclude_regions 20   100 does everything on chromosome 20 excluding base 100" }
  population_vcf_files: { type: 'string[]?', inputBinding: { position: 2, prefix: "--population_vcfs"}, doc: "Path to built-in tabix-indexed VCF file (or list of VCFs broken by chromosome), separated by comma or space, containing population allele frequencies." } 
  population_vcf_files_custom: { type: 'File[]?', secondaryFiles: [{pattern: ".tbi", required: true}], inputBinding: { position: 2, prefix: "--population_vcfs"}, doc: "Optional. Tabix-indexed VCF file (or list of VCFs broken by chromosome), separated by comma or space, containing population allele frequencies. Each of the item can be a file path, or a wildcard pattern." }
  sample_name_normal: { type: 'string?', inputBinding: { position: 2, prefix: "--sample_name_normal"}, doc: "Sample name for normal match to use for our sample_name in the output Variant/DeepVariantCall protos. If not specified, will be inferred from the header information from --reads_normal." }
  sample_name_tumor: { type: 'string?', inputBinding: { position: 2, prefix: "--sample_name_tumor"}, doc: "Sample name for tumor to use for our sample_name in the output Variant/DeepVariantCall protos. If not specified, will be inferred from the header information from --reads_tumor." }

  # Other File Args
  exclude_variants_vcf_file: { type: 'File?', secondaryFiles: [{pattern: ".tbi", required: false}], inputBinding: { position: 2, prefix: "--exclude_variants_vcf_filename"}, doc: "Optional. Population VCF (with AF) to exclude variants from. In our use case, this is currently only used for DeepSomatic tumor-only training examples creation." }  
  confident_regions_file: { type: 'File?', inputBinding: { position: 2, prefix: "--confident_regions"}, doc: "Regions that we are confident are hom-ref or a variant in BED format. In BED or other equivalent format, sorted or unsorted. Contig names must match those of the reference genome." }
  denovo_regions_file: { type: 'File?', inputBinding: { position: 2, prefix: "--denovo_regions"}, doc: "Regions where variants are de novo. Used to label variants as de novo." }
  
  # Catchall for Remaining Options
  extra_args: { type: 'string?', inputBinding: { position: 2, shellQuote: false }, doc: "Any additional args for deepsomatic make_somatic_examples" }

  cpu:
    type: 'int?'
    default: 1
    doc: "Number of CPUs to allocate to this task."
  ram:
    type: 'int?'
    default: 1
    doc: "Maximum GB of RAM to allocate for this tool."
outputs:
  examples:
    type: File
    outputBinding:
      glob: "*somatic.tfrecord*" 
  candidates:
    type: 'File?'
    outputBinding:
      glob: "candidate_positions*"
  gvcf:
    type: 'File?'
    outputBinding:
      glob: "*.gvcf.tfrecord*"
