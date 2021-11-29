cwlVersion: v1.1
class: CommandLineTool
id: gatk_mutect2
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: DockerRequirement
    dockerPull: 'broadinstitute/gatk:4.2.2.0'
  - class: ResourceRequirement
    ramMin: $(inputs.max_memory * 1000)
    coresMin: $(inputs.cores)
baseCommand: []
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      gatk --java-options "-Xmx${return Math.floor(inputs.max_memory*1000/1.074-1)}m" Mutect2

inputs:
  reference: { type: 'File', secondaryFiles: [^.dict, .fai], inputBinding: { prefix: "--reference", position: 2 }, doc: "Refernce fasta, with fai and dict indices" }
  input_tumor_aligned: { type: 'File', secondaryFiles: [{ pattern: ".bai", required: false },{ pattern: "^.bai", required: false },{ pattern: ".crai", required: false },{ pattern: "^.crai", required: false }], inputBinding: { prefix: "--input", position: 2 }, doc: "BAM/SAM/CRAM file containing tumor reads" }
  input_normal_aligned: { type: 'File?', secondaryFiles: [{ pattern: ".bai", required: false },{ pattern: "^.bai", required: false },{ pattern: ".crai", required: false },{ pattern: "^.crai", required: false }], inputBinding: { prefix: "--input", position: 2 }, doc: "BAM/SAM/CRAM file containing normal reads" }
  interval_list: { type: 'File?', inputBinding: { prefix: "--intervals", position: 2 }, doc: "One or more genomic intervals over which to operate in a GATK-interval compatible file format" }
  germline_resource_vcf: { type: 'File?', secondaryFiles: ['.tbi'], inputBinding: { prefix: "--germline-resource", position: 2 }, doc: "Population vcf of germline sequencing containing allele fractions in VCF format." }
  panel_of_normals: { type: 'File?', secondaryFiles: ['.tbi'], inputBinding: { prefix: "--panel-of-normals", position: 2 }, doc: "VCF file (and index) of sites observed in normal. A panel of normals can be a useful (optional) input to help filter out commonly seen sequencing noise that may appear as low allele-fraction somatic variants." }
  alleles: { type: 'File?', secondaryFiles: ['.tbi'], inputBinding: { prefix: "--alleles", position: 2 }, doc: "VCF file (and index) containing set of alleles for which to force genotyping regardless of evidence" }

  input_tumor_name: { type: 'string', inputBinding: { prefix: "--tumor-sample", position: 2 }, doc: "BAM sample name of tumor. May be URL-encoded as output by GetSampleName with -encode argument." }
  input_normal_name: { type: 'string?', inputBinding: { prefix: "--normal-sample", position: 2 }, doc: "BAM sample name of normal(s), if any. May be URL-encoded as output by GetSampleName with -encode argument." }
  output_vcf_name: { type: 'string', inputBinding: { prefix: "--output", position: 2 }, doc: "Name of vcf.gz file into which variants should be written" }
  output_f1r2_name: { type: 'string?', inputBinding: { prefix: "--f1r2-tar-gz", position: 2 }, doc: "Name of tar.gz file to collect F1R2 counts and output files. Optional but required to run LearnReadOrientationModel." }
  output_bam_name: { type: 'string?', inputBinding: { prefix: "--bam-output", position: 2 }, doc: "File to which assembled haplotypes should be written" }

  disable_read_filter: { type: ['null', { type: array, items: string, inputBinding: { prefix: "--disable-read-filter" }}], inputBinding: { position: 2 }, doc: "Read filters to be disabled before analysis" }
  disable_adaptive_pruning: { type: 'boolean?', inputBinding: { prefix: "--disable-adaptive-pruning", position: 2 }, doc: "Disable the adaptive algorithm for pruning paths in the graph. A single edge multiplicity cutoff for pruning doesn't work in samples with variable depths, for example exomes and RNA." }
  extra_args: { type: 'string?', inputBinding: { position: 3, shellQuote: false }, doc: "Any additional arguments for this tool. See GATK Documentation for complete list of options. Example input: --interval-merging-rule OVERLAPPING_ONLY" }

  cores: { type: 'int?', default: 3, doc: "Minimum reserved number of CPU cores for the task." }
  max_memory: { type: 'int?', default: 6, doc: "GB of RAM to allocate to the task." }

outputs:
  mutect2_vcf:
    type: File
    outputBinding:
      glob: '*.vcf.gz'
    secondaryFiles: [.tbi]
  mutect_stats:
    type: File
    outputBinding:
      glob: '*.stats'
  f1r2_counts:
    type: 'File?'
    outputBinding:
      glob: '*.f1r2_counts.tar.gz'
  mutect2_bam:
    type: 'File?'
    outputBinding:
      glob: '*.bam'
    secondaryFiles: [^.bai]
