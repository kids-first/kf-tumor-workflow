# Tumor-Only SNV Benchmarking Methods and Results
When it comes to somatic variant calling, having a matched normal sample from the patient is always preferred.
However, this is not always possible.
Mutect2, in our opinion as of this writing, is the best caller to use when a matched normal is not available in tumor-only calling mode.
Using the same gold standard dataset that we used to benchmark (see the first two sections of [this doc](https://github.com/kids-first/kf-somatic-workflow/blob/master/docs/SOMATIC_SNV_BENCHMARK.md)) our somatic with matched tumor-normal variant calling methods, we reviewed filtering strategies from other publications ([TOSCA](https://academic.oup.com/bioinformaticsadvances/article/2/1/vbac070/6717791?login=false) and [HIVE Dragen Pipeline](https://www.fda.gov/science-research/fda-science-forum/hive-dragen-pipeline-enables-somatic-tumor-only-filtration)), and chose the best combination that maximized F1 score while minimizing true positive (TP) loss as aggressively pursuing false positive (FP) reduction and results in a massive loss of TP calls.
Our resultant filtering strategy can be found [here](https://github.com/kids-first/kf-tumor-workflow?tab=readme-ov-file#kids-first-drc-tumor-only-pipeline) with benchmark results below.

## Initial Assessment
Using the same gold standard dataset for somatic, we compared the tumor-only calls.
The focus was on Whole Genome Sequencing (WGS); at the end we address/added some extra criteria after examining proposed filtering on whole exome sequencing (WXS) samples
This initial assessment was used to gauge our current filtering strategy, which generally matched that of the somatic workflow of:
 - Variant `FILTER` = `PASS` or `INFO HotSpotAllele=1`
 - `FILTER` GNOMAD_AF < 0.001 applied
 - GDC Panel of Normals (PON) used - this is different from somatic workflow

### Cell Line Tumor Only Benchmark, High + Medium Confidence Calls
| Sample ID | True Positive | False Positive | False Negative | F1 Score | Precision | Sensitivity |
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| FD_T_2|  36812 |  31817 |  78222 |  0.4 | 0.54 | 0.32 |
| FD_T_3|  37045 |  31991 |  78371 |  0.4 | 0.54 | 0.32 |
| NV_T_2|  38596 |  39479 |  79858 |  0.39 | 0.49 | 0.33 |
| NV_T_3|  38780 |  37809 |  79629 |  0.4 | 0.51 | 0.33 |
| IL_T_1|  37817 |  36644 |  78340 |  0.4 | 0.51 | 0.33 |
| IL_T_2|  37445 |  34935 |  77361 |  0.4 | 0.52 | 0.33 |
| NS_T_1|  37504 |  40680 |  76893 |  0.39 | 0.48 | 0.33 |
| NS_T_2|  35418 |  31957 |  75577 |  0.4 | 0.53 | 0.32 |
### Synthetic Tumor Only Benchmark, High + Medium Confidence Calls
| Sample ID | True Positive | False Positive | False Negative | F1 Score | Precision | Sensitivity |
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| syn_WHS_FD_N_2|  87309 |  48421 |  27757 |  0.7 | 0.64 | 0.76 |
| syn_WGS_FD_N_3|  86609 |  48848 |  28215 |  0.69 | 0.64 | 0.75 |
| syn_WGS_NV_N_2|  90573 |  59350 |  24607 |  0.68 | 0.6 | 0.79 |
| syn_WGS_NV_N_3|  89412 |  60020 |  25588 |  0.68 | 0.6 | 0.78 |
| syn_WGS_IL_N_1|  88403 |  58522 |  25918 |  0.68 | 0.6 | 0.77 |
| syn_WGS_IL_N_2|  87446 |  56718 |  26852 |  0.68 | 0.61 | 0.77 |
| syn_WGS_NS_N_1|  86835 |  75729 |  27403 |  0.63 | 0.53 | 0.76 |
| syn_WGS_NS_N_2|  84797 |  55039 |  29212 |  0.67 | 0.61 | 0.74 |

## Variant Filter fine-tuning
Given that the focus is on reducing FPs from tumor only calling and that all synthetic and cell line samples have very high overlap in calls, only one samples was used to benchmark various filtering strategies and combinations.
Many combinations of filters were performed from the two publications.
Below we attempt to concisely summarize the story and results.  

### Test synthetic sample `syn_WGS_IL_N_1` HIVE filter components and final combination
- RAW: the protected file as output by our pipeline
- PASS: Stripped the public file soft filters from the RAW and grabbed only the PASS records
- PASS+GNO 0.1%: Applied the KF gnomAD filter (exclude gnomad >= 0.001) directly to the PASS file. These results should be a pretty close approximation of what our public files look like
- PASS+GNO 0.003%: Applied the HIVE gnomAD filter (exclude gnomad >= 0.00003) directly to the PASS file. I used this as a comparison to see just how much value we were getting from the other filtration steps
- PASS+PON: Removed any records from the PASS file that were also found in the GATK PON
- PASS+PON+dbSNP: Removed any records from PASS+PON that were labeled COMMON in dbSNP
- PASS+PON+dbSNP+ClinVar: Removed any records from PASS+PON+dbSNP that had "benign" in their worst VEP assigned CLN_SIG
- HIVE: Removed any records from PASS+PON+dbSNP+ClinVar that had a gnomad AF >= 0.00003

| Mode | True Positive | False Positive | False Negative | F1 Score | Precision | sensitivity |
:-|-:|-:|-:|:-|:-|:-
| RAW |  101252 |  1586166 |  13069 |  0.11 | 0.06 | 0.89 |
| PASS |  87868 |  85432 |  26453 |  0.61 | 0.51 | 0.77 |
| PASS+GNO 0.1% |  87815 |  46161 |  26506 |  0.71 | 0.66 | 0.77 |
| PASS+GNO 0.003% |  87259 |  41473 |  27062 |  0.72 | 0.68 | 0.76 |
| PASS+PON |  87868 |  85185 |  26453 |  0.61 | 0.51 | 0.77 |
| PASS+PON+dbSNP |  87789 |  78862 |  26532 |  0.62 | 0.53 | 0.77 |
| PASS+PON+dbSNP+ClinVar |  87789 |  78855 |  26532 |  0.62 | 0.53 | 0.77 |
| HIVE |  87252 |  41083 |  27069 |  0.72 | 0.68 | 0.76 |

The `HIVE` method is not much better than `PASS+GNO 0.003%` and is more complicated to implement.

### Test call line sample `WGS_IL_T_1` HIVE filter components and final combination

| Mode | True Positive | False Positive | False Negative | F1 Score | Precision | sensitivity |
:-|-:|-:|-:|:-|:-|:-
| RAW|  49474 |  1590429 |  60438 |  0.06 | 0.03 | 0.45 |
| PASS|  40116 |  85168 |  69796 |  0.34 | 0.32 | 0.36 |
| PASS+PON|  40113 |  84922 |  69799 |  0.34 | 0.32 | 0.36 |
| PASS+PON+dbSNP|  40010 |  76907 |  69902 |  0.35 | 0.34 | 0.36 |
| PASS+PON+dbSNP+ClinVar|  40010 |  76898 |  69902 |  0.35 | 0.34 | 0.36 |
| HIVE|  38728 |  33237 |  71184 |  0.43 | 0.54 | 0.35 |
| PASS+GNO 0.1%|  39659 |  41303 |  70253 |  0.42 | 0.49 | 0.36 |
| PASS+GNO 0.003%|  38734 |  33608 |  71178 |  0.43 | 0.54 | 0.35 |

### A note on the TOSCA methodology
This method's main difference was adding variant allele frequency (VAF > 5%) and read depth (DP > 250) cutoffs, as well as consequence (CSQ).
Consequence ended up being extremely aggressive, ~99% of true positives, while the depth filter is obviously too high for WGS, but as a concept can be helpful.
### Final summary
The table below summarized efforts of trying PON, blacklists, and various filters, with `GNOMAD_HARD` representing gnomad >= 0.00003:
| Sample | Mode | True Positive | False Positive | False Negative | F1 Score | Precision | sensitivity |
:-|:-|-:|-:|-:|:-|:-|:-
| No PoN | PASS |  41382 |  194272 |  68530 |  0.24 | 0.18 | 0.38 |
| No PoN | Blacklist PASS|  41120 |  174740 |  68792 |  0.25 | 0.19 | 0.37 |
| Broad PoN | PASS |  41106 |  161252 |  68806 |  0.26 | 0.2 | 0.37 |
| Broad PoN | Blacklist PASS|  40872 |  144031 |  69040 |  0.28 | 0.22 | 0.37 |
| GDC PoN | PASS |  40116 |  85168 |  69796 |  0.34 | 0.32 | 0.36 |
| GDC PoN | Blacklist PASS|  39899 |  72728 |  70013 |  0.36 | 0.35 | 0.36 |
| No PoN | HIVE (PASS+GNOMAD_HARD) |  38878 |  53869 |  71034 |  0.38 | 0.42 | 0.35 |
| No PoN | Blacklist HIVE (PASS+GNOMAD_HARD) |  38726 |  38762 |  71186 |  0.41 | 0.5 | 0.35 |
| Broad PoN | HIVE (PASS+GNOMAD_HARD) |  38867 |  45983 |  71045 |  0.4 | 0.46 | 0.35 |
| Broad PoN | Blacklist HIVE (PASS+GNOMAD_HARD) | 38734 |  32378 |  71178 |  0.43 | 0.54 | 0.35 |
| GDC PoN | HIVE (PASS+GNOMAD_HARD) |  38734 |  33608 |  71178 |  0.43 | 0.54 | 0.35 |
| GDC PoN | Blacklist HIVE (PASS+GNOMAD_HARD) | 38568 |  22674 |  71344 |  0.45 | 0.63 | 0.35 |
| No PoN | TOSCA (PASS+GNOMAD_HARD+DP30)|  37140 |  32411 |  72772 |  0.41 | 0.53 | 0.34 |
| No PoN | Blacklist TOSCA (PASS+GNOMAD_HARD+DP30)| 36996 |  27029 |  72916 |  0.43 | 0.58 | 0.34 |
| Broad PoN | TOSCA (PASS+GNOMAD_HARD+DP30)|  37129 |  25406 |  72783 |  0.43 | 0.59 | 0.34 |
| Broad PoN | Blacklist TOSCA (PASS+GNOMAD_HARD+DP30)| 37004 |  21426 |  72908 |  0.44 | 0.63 | 0.34 |
| GDC PoN | TOSCA (PASS+GNOMAD_HARD+DP30)|  36995 |  15999 |  72917 |  0.45 | 0.7 | 0.34 |
| GDC PoN | Blacklist TOSCA (PASS+GNOMAD_HARD+DP30)| 36837 |  13898 |  73075 |  0.46 | 0.73 | 0.34 |

Lastly, for WGS filtering, Mutect2 results may have 0 allele depths for some calls. This is explained [here](https://gatk.broadinstitute.org/hc/en-us/articles/360035532252-Allele-Depth-AD-is-lower-than-expected) and has spurred a common sense addition of `AD > 0`

## WXS Filtering
WGS filtering is generalizable and safe to use for WXS, but WXS just needs a little more. After examining a well-characterized matched tumor-normal WXS set from the OpenPedCan project, we learned the following: Many high depth regions have VAF < 1%. An example was a sample having 6495 of the 7357 are present at an VAF below 1%. If compared to the corresponding matched-normal analysis, only one of those calls were deemed valid with VAF < 1%. Therefore, setting 1% as a floor for WXS is more than reasonable
