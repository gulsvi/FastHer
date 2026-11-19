# Data Files Description

## UK Biobank Datasets

### 200k Samples (BPIFC gene)
- 'ukb_200k_chr22_Z_BPIFC.txt' - Simulated phenotype Z-scores
- 'ukb_200k_chr22_SNPorder_BPIFC.bim' - SNP info 
- 'ukb_200k_chr22_LD_BPIFC.npz' - LD matrix for region from chromosome 22


## File Formats
- **.txt**: Tab-separated, columns: Predictor, Z, n
- **.bim**: PLINK format, columns: CHR, SNP, cM, POS, A1, A2  
- **.npz**: Numpy compressed LD matrix 
### Simulated locus with allele frequencies (`sim5k_chr22_1KG_*`)
Example for the SNP-weighted model (`FastHer(..., use_weights = TRUE)`).
- 617 SNPs of chromosome 22 (1000 Genomes, European); genotypes expanded to N = 5000 individuals.
- The phenotype is simulated: LDAK model (power -0.25, all SNPs causal), true local h2 = 0.1.
- 'sim5k_chr22_1KG_Z.txt' - z-statistics (BETA/SE from GCTA fastGWA-lr), n, p-value and the allele
  frequency of A1 (column AF1, used to compute the SNP weights).
- 'sim5k_chr22_1KG_SNPorder.bim' - SNP order, positions and alleles.
- 'sim5k_chr22_1KG_LD.npz' - in-sample LD matrix (Pearson correlations).
