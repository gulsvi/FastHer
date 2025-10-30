# Data Files Description

## UK Biobank Datasets

### 200k Samples (BPIFC gene)
- 'ukb_200k_chr22_Z_BPIFC.txt' - Simulated phenotype Z-scores
- 'ukb_200k_chr22_SNPorder_BPIFC.bim' - SNP info 
- 'ukb_200k_chr22_LD_BPIFC.npz' - LD matrix for region from chromosome 22

### 200k Samples (HIRA gene)
- 'ukb_200k_chr22_Z_HIRA.txt' - Simulated phenotype Z-statistics  
- 'ukb_200k_chr22_SNPorder_HIRA.bim' - SNP info
- 'ukb_200k_chr22_LD_HIRA.npz' - LD matrix for region from chromosome 22

## File Formats
- **.txt**: Tab-separated, columns: Predictor, Z, n
- **.bim**: PLINK format, columns: CHR, SNP, cM, POS, A1, A2  
- **.npz**: Numpy compressed LD matrix