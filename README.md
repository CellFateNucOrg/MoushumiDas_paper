# Mouhsumi Das et al. paper figures

Code and data objects for RNAseq figures from Das _et al._ (2024) **Condensin I folds the Caenorhabditis elegans genome** 
The pipeline scripts for processing the raw RNAseq reads with Salmon and DESeq2 and other custom scripts can be found at the following repository:
[SMC_RNAseq](https://github.com/CellFateNucOrg/SMC_RNAseq/tree/v0.2) repository.

Old version: https://www.biorxiv.org/content/10.1101/2022.06.14.495661v2
Published version: 

## Main figures

### Figure 2
Panels C & D created with _`Fig_2c&d_eigenVectors.R`_ script
Chromatin states were taken from Evans et al. (2016, PMID: 27791097) and lifted over to Ce11

### Figure 3
Panel I created using _`Fig_3i_compartmentSizes.R`_ script

### Figure 4
Panel G is based on plot created using _`Fig_4g_ChIPseq.R`_ script and publically available ChIP seq data lifted over to Ce11:

```
DPY27_N2_L3_GSE67650_ce11.bw, # Kramer et al (2015)
GSE87741_DPY30_N2_MxEmb_MACS141_e5_average_ce11.bw, # Albritton et al (2017)
GSE87741_SDC2_N2_MxEmb_MACS141_e5_average_ce11.bw, # Albritton et al (2017)
GSE87741_SDC3_N2_MxEmb_MACS141_e5_average_ce11.bw, # Albritton et al (2017)
KLE2_N2_L3_GSE45678_ce11.bw, # Kranz et al. (2013)
GSM4293382_SCC1_Q0835_mE16_E_N2_L3_rep1_ce11.bw # Huang et al. (2022)
```

### Figure 5
Figure of kleisin subunit cleavage RNAseq and SDC-3AID RNAseq were created with _`Fig_5_RNAseq_TEV.R`_ script.

## Supplementary figures

### Figure S1
Panel A quantification was carried out in imageJ using the Analyze>Gels tools to select lanes and plot intensity profile. Lines were drawn at the base of peaks to subtract local background and the area of the peak was recorded. Plots of quantification created with _`Fig_S1_SMCcleavageQuantification.R`_ script.

### Figure S2
Panel B created with _`Fig_S2b_eigenVectors.R`_ scripts
Chromatin states were taken from Evans et al. (2016, PMID: 27791097) and lifted over to Ce11

### Figure S4
Created using _`Fig_S4_compartmentSwitch.R`_ script
- Chromatin states were taken from Evans et al. (2016) and lifted over to Ce11
- The same Chromatin IP tracks that were used in Evans et al. (2016) were used here after lifting over to Ce11
- Genome features were obtained by downloading WS280 genome annotation GFF file and extracting different types into separate BED files (see _getGFFdata.R_ script)

### Figure S6
Created using _`Fig_S6_RNAseqTEV.R`_ script

### Figure S7
Created using _`Fig_S7_RNAseqTEV_2.R`_ script

### Figure S9
Panels F & G Created using _`Fig_S9f&g_RNAseq_sdc3.R`_ script

# DESeq2 results tables
These can be found in the folder _`p0.05_lfc0.5_filtCycChrAX`_ and are used by the scripts.

| **Pertubation**                            | **Description** | **Filename**                                                              | 
|--------------------------------------------|-----------------|-------------------------------------------------------------------------------|
| SCC-1 cleavage         | Cohesin<sup>SCC-1</sup> cleavage by TEV     | _`filtCycChrAX_X.wt.wt.0mM_scc16cs_vs_wt_DESeq2_fullResults_p0.05.rds`_    |
| COH-1 cleavage         | Cohesin<sup>COH-1</sup> cleavage by TEV  | _`filtCycChrAX_X.wt.wt.0mM_coh1cs_vs_wt_DESeq2_fullResults_p0.05.rds`_ |
| SCC-1 & COH-1 cleavage | Cohesin<sup>SCC-1</sup> and cohesin<sup>COH-1</sup> cleavage by TEV   | _`filtCycChrAX_X.wt.wt.0mM_scc1coh1cs_vs_wt_DESeq2_fullResults_p0.05.rds`_    |
| DPY-26 cleavage        | Condensin I/I<sup>DC</sup> cleavage by TEV      | _`filtCycChrAX_X.wt.wt.0mM_dpy26cs_vs_wt_DESeq2_fullResults_p0.05.rds`_  |
| KLE-2 cleavage         | Condensin II cleavage by TEV  | _`filtCycChrAX_X.wt.wt.0mM_kle2cs_vs_wt_DESeq2_fullResults_p0.05.rds`_   |
| SDC-3 degron           | SDC-3 degradation by TIR1+auxin  | _`filtCycChrAX_wt.TIR1.sdc3deg.X_1mM_vs_0mM_DESeq2_fullResults_p0.05.rds`_   |

