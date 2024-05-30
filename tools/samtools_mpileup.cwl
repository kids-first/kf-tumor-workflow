cwlVersion: v1.2
class: CommandLineTool
id: samtools_mpileup
doc: |
  produces "pileup" textual format from an alignment
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'pgc-images.sbgenomics.com/danmiller/samtools_parallel:1.20'
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    ramMin: 32000
    coresMin: $(inputs.threads)

baseCommand: []
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      zcat $(inputs.snp_vcf.path) | grep -v "#" | awk {'printf ("%s\t%s\t%s\t%s\t%s\n", $1,$2-1,$2,$4,$5)'} > snps.bed
  - position: 10
    shellQuote: false
    prefix: "&&"
    valueFrom: >-
      parallel --jobs $(inputs.threads) --keep-order --colsep '\t' samtools mpileup -f $(inputs.reference.path) -r {1} -l snps.bed -d 8000 -Q 0 -q 1 $(inputs.input_reads.path) :::: $(inputs.calling_regions ? inputs.calling_regions.path : inputs.reference.secondaryFiles[0].path) > $(inputs.output_basename).miniPileup

inputs:
  input_reads: {type: File, secondaryFiles: ['^.bai']}
  threads:
    type: int? 
    default: 16
  reference: {type: File, secondaryFiles: [.fai]}
  calling_regions: {type: 'File?'}
  snp_vcf: {type: File}
  output_basename: {type: 'string?', default: "delme_test", doc: "String to use as basename for outputs"}
outputs:
  pileup:
    type: File
    outputBinding:
      glob: '*.miniPileup'
