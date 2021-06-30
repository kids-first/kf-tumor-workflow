cwlVersion: v1.2
class: CommandLineTool
id: gatk_mutect2_create_pon
label: GATK Mutect2 Create Panel of Normals

requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: DockerRequirement
    dockerPull: 'pgc-images.sbgenomics.com/d3b-bixu/gatk:4.1.1.0'
  - class: ResourceRequirement
    ramMin: 4000
    coresMin: 2

baseCommand: [/gatk, GenomicsDBImport]

arguments:
  - position: 0
    shellQuote: false
    valueFrom: >-
      --java-options ${var str = "-Xmx"; return str.concat(inputs.ram, "g")}
      --genomicsdb-workspace-path pon_db -R ${inputs.reference.path}
      -L ${inputs.interval_list.path}
      ${
        var vcf_string;
        for (var i=0; i<inputs.input_vcfs.length; i++){
          vcf_string += " -V" + inputs.input_vcfs[i].path;
        }
        return(vcf_string);
      }
      && gatk CreateSomaticPanelOfNormals -R ${inputs.reference.path}
      -V genedb://pon_db -O ${inputs.output_basename}.vcf.gz

inputs:
  reference: {type: "File", secondaryFiles: [^.dict, .fai]}
  input_vcfs: {type: "File[]", doc: "List of input vcf files"}
  interval_list: {type: "File", doc: "Interval list file"}
  output_basename: {type: string, doc: "Output panel of normal vcf basename"}
  ram: {type: "int?", default: 2, doc: "In GB"}
  cpu: {type: "int?", default: 4, doc: "Number of CPUs to request"}

outputs:
  panel_of_normal_vcf:
    type: File
    outputBinding:
      glob: '*.vcf.gz'
