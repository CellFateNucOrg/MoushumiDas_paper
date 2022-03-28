# Mouhsumi Das et al. paper figures

Partial code and data objects for figures from Das et al. 

## Main figures

### Figure 2
Panels C, D & H created with _HiCfeatures.R_ script
Chromatin states were taken from Evans et al. (2016) and lifted over to Ce11

### Figure 3
Panel G created using _ChIPseqPlots.R_ script and publically available ChIP seq data lifted over to Ce11:

'''
"DPY27_N2_L3_GSE67650_ce11.bw", # Kramer et al (2015)
"GSE87741_DPY30_N2_MxEmb_MACS141_e5_average_ce11.bw", # Albritton et al (2017)
"GSE87741_SDC2_N2_MxEmb_MACS141_e5_average_ce11.bw", # Albritton et al (2017)
"GSE87741_SDC3_N2_MxEmb_MACS141_e5_average_ce11.bw", # Albritton et al (2017)
"KLE2_N2_L3_GSE45678_ce11.bw", # Kranz et al. (2013)
"GSM4293382_SCC1_Q0835_mE16_E_N2_L3_rep1_ce11.bw" # Huang et al. (2022)
'''

### Figure 4 
Panels for SDC-3AID RNAseq created with _RNAseqFig_deg.R_ script

### Figure 5
Figure of kleisin subunit cleavage RNAseq created with _RNAseqFig_TEV.R_ script

## Supplementary figures

### Figure S2
Panel B created with _HICfeatures.R_ scripts
Chromatin states were taken from Evans et al. (2016) and lifted over to Ce11

### Figure S4
Created using _compartmentSwitch_complexHeatmap.R_ script
- Chromatin states were taken from Evans et al. (2016) and lifted over to Ce11
- The same Chromatin IP tracks that were used in Evans et al. (2016) were used here after lifting over to Ce11
- Genome features were obtained by downloading WS280 genome annotation GFF file and extracting different types into separate BED files (see _getGFFdata.R_ script)

### Figure S6
Created using _RNAseqFigSupl_TEV.R_ script

### Figure S7
Created using _RNAseqFigSupl_TEV.R_ script

### Figure S9
Created using _RNAseqFigSupl_TEV.R_ script

### Figure S11
Created using _RNAseqFigSupl_deg.R_ script
