# Kids First Mutect2 Tumor-Only Filtering

The Kids First Mutect2 tumor-only filtering approach began with an investigation of strategies that performed
cutoff-based filters to reduce the false discovery rate. These approaches were:
- FDA HIVE DRAGEN: https://www.fda.gov/science-research/fda-science-forum/hive-dragen-pipeline-enables-somatic-tumor-only-filtration
- TOSCA: https://doi.org/10.1093/bioadv/vbac070

These cutoff-based approaches primarily rely on population frequencies from large genomic sequencing projects (ESP, ExAC,
1000G, dbSNP, and gnomAD). Additionally, these approaches apply filters based on:
- Clinvar's clinical severity
- VEP's variant consequence
- Read Depth
- Variant Allele Frequency

In addition to these approaches, the investigation also looked at the effect of providing a panel of normals and/or
calling blacklists to Mutect2.

## Benchmarking Dataset

In order to evaluate these approaches, the investigation made use of a robust variant calling dataset defined here:
https://pubmed.ncbi.nlm.nih.gov/34504347/. The data from that paper can be found here: https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/data/.
This dataset contains two sets of variant calling truth standards. The first of these is a synthetic truth. In this set,
the normal sample reads have been spiked with variants. The second of these is a consensus set. This set was created
through aligning the sample reads with multiple alignment software then calling variants in each of those alignment with
a variety of variant calling software. Calls were then assigned confidence based on the preponderance of support from the
matrix of calls.

The truth sets for all samples can be found here:
- Synthetic Truth
   - Per-sample: https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/data/DeepLearning_bams/truth_vcf/
- Tumor-Normal Truth
   - Per-sample, Per-software: https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/analysis/SNVs/vcfs/WGS/
   - All-samples, All-software (Superset): https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/release/v1.2/

## Benchmarking Approach

First Mutect2 tumor-only variant calling was performed on two inputs samples:
- Synthetic Tumor
- Cell Line Tumor
Both of these samples has variants called multiple times with accompanying inputs:
- With/Without Panel of Normals (both Genomic Data Commons (GDC) and Mutect2 1000G)
- With/Without Calling Blacklist

The resultant variant files were then passed through the filters described above and compared against their respective
truth sets. The synthetic tumor variants were compared against the per-sample synthetic truth VCF. The cell line tumor-only
variants were compared against the sample's BWA-Mutect2 tumor-normal results. The comparison was done on a by-site basis.
Variants were considered a match if they matched all: chromosome, position, reference allele, and alternative allele.

## Recommendations from Benchmarking

Tumor-only calling with Mutect2 is inherently very noisy. Therefore, most of our recommendations revolve around reducing
the noise where possible. These recommendations are:
- Restrict the callable regions with a blacklist and GDC Panel of Normals (PON)
- Filter low support reads:
   - Allele Depth (AD) > 0: WGS uninformative reads
   - Variant Allele Frequency (VAF) > 1%: WXS noise
- Filter potential germline variants with gnomAD AF < 0.00003
- Only keep variants that have no FILTERS/are PASS
- Rescue any variants that fall in hotspot regions/genes

### Calling Parameters

#### Panel of Normals

Recommendation: Use GDC PON

In our investigation, we looked at the difference in performance when providing a panel of normals. We tested two panels.
The first was the public Panel of Normals created from 1000G samples by the Broad Institute here: https://storage.cloud.google.com/gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz
The second was the controlled-access Panel of Normals created by the GDC team here: https://gdc.cancer.gov/about-data/gdc-data-processing/gdc-reference-files
When compared to providing no panel at all, both panels had better precision. Between the panels, the GDC panel outperformed
the Broad panel. The difference is significant enough that, if you can obtain it, we recommend using the GDC panel over
the Broad panel. If you cannot access the GDC panel and cannot make your own panel, we recommend using the Broad panel.
Any panel is better than no panel for tumor-only.

#### Blacklist

Recommendation: Use GENCODE Blacklist

Another good way to cut down on noise is to restrict calling in difficult-to-map and repetitive regions. These regions
tend to produce a lot of calls for short reads technology. In our investigation, we looked at the effect of providing a
blacklist to Mutect2. The blacklist we used was sourced from here: https://github.com/Boyle-Lab/Blacklist/blob/master/lists/hg38-blacklist.v2.bed.gz
In all cases, providing the blacklist resulting in an increase in precision with almost no cost to sensitivity. As such,
we recommend removing difficult to call regions from your input intervals.

### Remove Low Support Variants

Recommended Filter:
- AD > 0
- VAF > 0.01/1%

In our investigation we came across a couple of unique scenarios for Mutect2 tumor-only calls. The first of these was the
presence of variants with no read support. These would be variants where the alternative allele depth was zero. After a
bit of digging we discovered the cause to be what GATK refers to as "uninformative reads". You can read about that here:
https://gatk.broadinstitute.org/hc/en-us/articles/360035532252-Allele-Depth-AD-is-lower-than-expected. Long story short,
the reads that support the variants in question have been deemed "uninformative" which is defined as a read that "passes
the quality filters, but the likelihood of the most likely allele given the read is not significantly larger than the
likelihood of the second most likely allele given the read." These reads tend to show up for whole genome experiments
in intronic regions, particularly around repetitive areas. Ultimately, we decided that it would be best to remove any
reads that do not have any "informative" support. As such, we recommend only collecting the reads that have an allele
depth greater than zero.

The second of these scenarios was variants with extremely low allele frequencies. This scenario was found in whole exome
and other targeted sequencing approaches. It seems like Mutect2 tumor-only is not as efficient as somatic with removing
variants with extremely low AF. With somatic we typically would see a minimum VAF around 1%. For tumor-only, we would
sometimes see samples with thousands of variants below 0.5%. This observation was not universal but it was happening
enough that we recommend setting a minimum VAF of 1%.

Overall, these values are only meant to act as a floor to exclude highly suspect variants. We found that placing these
filters on the data had no adverse effect on our sensitivity.

### Population Allele Frequency Filter

Recommended filter:
- gnomAD AF < 0.00003

As mentioned above, several different population allele frequency databases are available and can be used for annotating
and filtering variants. In our investigation, however, we found gnomAD to be the most impactful when it came to filtering.
gnomAD filtering already exists as part of our germline filtering approach for the somatic workflow so it makes sense to
also use it here.

A key difference from our somatic pipeline is that the gnomAD filter is much more stringent. For somatic filtration,
variants with a gnomAD AF over 0.001 are removed for being potentially germline. For tumor-only filtration, this threshold
was found to be overly permissive. FDA HIVE DRAGEN's poster discusses this threshold in a bit more detail and ultimately
decided that a much more stringent gnomAD AF value of 0.00003 minimized the false discovery rate without negatively
impacting the sensitivity of Mutect2. Our findings agreed with theirs. Therefore, we recommend removing variants with
an associated gnomAD AF value over 0.00003/0.003%.

### The Public File Filter: PASS or HotSpot

Recommended filter:
- FILTER = PASS | HotSpotAllele = 1

The final filter we recommend is our public file filter from our somatic and germline workflows. This filter is not so
much a filter as it is a rescue for genes of interest. In short, we recommend keeping variants that pass all the filters
mentioned above OR are present in a hotspot region/gene.
