# Mouhsumi Das et al. paper figures

Partial code and data objects for figures from Das et al. 

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
"GSE87741_DPY30_N2_MxEmb_MACS141_e5_average_ce11.bw", # Albritton et al (2017)
"GSE87741_SDC2_N2_MxEmb_MACS141_e5_average_ce11.bw", # Albritton et al (2017)
"GSE87741_SDC3_N2_MxEmb_MACS141_e5_average_ce11.bw", # Albritton et al (2017)
"KLE2_N2_L3_GSE45678_ce11.bw", # Kranz et al. (2013)
"GSM4293382_SCC1_Q0835_mE16_E_N2_L3_rep1_ce11.bw" # Huang et al. (2022)
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

