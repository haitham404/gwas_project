# GWAS QC Experiments

## 1. Baseline

Parameters
- MAF: 0.05
- SNP Missingness (`--geno`): 0.05
- Sample Missingness (`--mind`): 0.02
- HWE: Default threshold

Summary
- Minimum MAF: 0.0512821
- Maximum MAF: 0.5000000
- Minimum SNP Missingness: 0
- Maximum SNP Missingness: 0.00641026
- Minimum Sample Missingness: 5.90537e-05
- Maximum Sample Missingness: 0.0108806

---

## 2. geno10 Experiment

Changed Parameter
- SNP Missingness (`--geno`): 0.10

Summary
- Minimum MAF: 0.0512821
- Maximum MAF: 0.5000000
- Minimum SNP Missingness: 0
- Maximum SNP Missingness: 0.00641026
- Minimum Sample Missingness: 5.90537e-05
- Maximum Sample Missingness: 0.0108806

---

## 3. relaxMAF Experiment

Changed Parameter
- MAF: 0.10

Summary
- Minimum MAF: 0.1000000
- Maximum MAF: 0.5000000
- Minimum SNP Missingness: 0
- Maximum SNP Missingness: 0.00641026
- Minimum Sample Missingness: 5.86717e-05
- Maximum Sample Missingness: 0.0116170

---

## 4. strictHWE Experiment

Changed Parameter
- HWE: More stringent threshold

Summary
- Minimum MAF: 0.0512821
- Maximum MAF: 0.5000000
- Minimum SNP Missingness: 0
- Maximum SNP Missingness: 0.00641026
- Minimum Sample Missingness: 5.90537e-05
- Maximum Sample Missingness: 0.0108806