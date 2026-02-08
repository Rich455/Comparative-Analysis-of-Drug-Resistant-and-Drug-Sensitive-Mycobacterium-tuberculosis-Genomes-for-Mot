

#  Hi, I’m Richard

**M.D. | M.S. Computer Science (Bioinformatics & HPC)**  
  

---

##  Comparative Analysis of Drug-Sensitive and Drug-Resistant *Mycobacterium tuberculosis*

###  Project Overview

This project compares two strains of *Mycobacterium tuberculosis* to understand the genetic basis of drug resistance.
Source of the genomes data is National Library of the Medicine.

- **Drug-sensitive strain:** *Mycobacterium tuberculosis* H37Rv  
- **Drug-resistant strain:** *Mycobacterium tuberculosis* SCAID 320.0  

H37Rv is a globally used laboratory reference strain, while SCAID 320.0 is a clinical strain associated with drug resistance.

---

##  Strain Comparison

| Feature | Drug-Sensitive Strain | Drug-Resistant Strain |
|-------|----------------------|----------------------|
| Strain | H37Rv | SCAID 320.0 |
| Resistance Status | Drug-sensitive | Drug-resistant |
| Genome Type | Complete genome | Complete genome |
| Research Use | Global reference | Clinical resistance model |

---

##  Resistance Genes Analyzed

| Gene | Associated Drug | Biological Function |
|-----|---------------|-------------------|
| **rpoB** | Rifampicin | RNA polymerase β-subunit |
| **katG** | Isoniazid | Catalase-peroxidase (drug activation) |
| **inhA** | Isoniazid | Mycolic acid synthesis |
| **gyrA** | Fluoroquinolones | DNA gyrase |
| **embB** | Ethambutol | Cell wall arabinan synthesis |

---

##  Project Objectives

- Align resistance genes between two TB genomes  
- Identify true nucleotide mutations (SNPs)  
- Calculate mutation counts and percentages per gene  
- Interpret biological significance of detected mutations  

---

##  Method Summary

1. Load complete genomes from FASTA files  
2. Extract resistance genes using H37Rv reference coordinates  
3. Align each gene **gene-to-gene** using pairwise alignment  
4. Detect nucleotide-level differences  
5. Calculate mutation percentages per gene  

---

##  Results

### Genome Sizes

- **Drug-sensitive genome (H37Rv):** 4,411,532 bp  
- **Drug-resistant genome (SCAID 320.0):** 4,406,628 bp  

---

###  Gene-Level Mutation Results

#### **rpoB (Rifampicin Resistance)**
- Gene length: 3,518 bp  
- Mutations detected: 2  
- Mutation percentage: **0.06%**

 Image 01
 
 ![Alt text describing the image](https://github.com/Rich455/Comparative-Analysis-of-Drug-Resistant-and-Drug-Sensitive-Mycobacterium-tuberculosis-Genomes-for-Mot/blob/main/MTB%20analysis%2001.png)
---




**Interpretation:**  
Low mutation frequency aligns with known rifampicin resistance hotspots.

---

#### **katG (Isoniazid Resistance)**
- Gene length: 2,223 bp  
- Mutations detected: 845  
- Mutation percentage: **38.01%**

**Interpretation:**  
High variability suggests strong evolutionary pressure.  
This gene is essential for isoniazid activation.

---

#### **inhA (Isoniazid Resistance)**
- Gene length: 1,347 bp  
- Mutations detected: 460  
- Mutation percentage: **34.15%**
- 

**Interpretation:**  
Extensive mutations may reflect drug-target alteration mechanisms.

---

#### **gyrA (Fluoroquinolone Resistance)**
- Gene length: 1,517 bp  
- Mutations detected: 518  
- Mutation percentage: **34.15%**

**Interpretation:**  
Mutations in *gyrA* are well known to confer fluoroquinolone resistance.

---

Image 02

![Alt text describing the image](https://github.com/Rich455/Comparative-Analysis-of-Drug-Resistant-and-Drug-Sensitive-Mycobacterium-tuberculosis-Genomes-for-Mot/blob/main/MTB%20analysis%2002.png)
---


#### **embB (Ethambutol Resistance)**
- Gene length: 3,291 bp  
- Mutations detected: 1,074  
- Mutation percentage: **32.63%**

**Interpretation:**  
High mutation frequency suggests adaptive changes in cell wall synthesis.

---

## Biological Interpretation

- Very low mutation rate in **rpoB** reflects known resistance biology  
- High mutation rates in other genes may indicate:
  - Structural gene variation  
  - Long-term drug pressure  
  - Insertions/deletions affecting alignment  

This highlights the importance of **correct gene-level alignment** in comparative genomics.

---

##  Why This Project Matters

- Demonstrates proper comparative genomics workflow  
- Avoids false mutation inflation caused by misalignment  
- Applicable to:
  - TB drug-resistance research  
  - Bioinformatics education  
  - Clinical genomics  
  - Genomic surveillance  

---

## Future Improvements

- Translate nucleotide mutations to amino acid changes  
- Focus on WHO-validated resistance hotspots  
- Add mutation density visualizations  
- Compare multiple resistant strains  

---

##  Tools Used

- Python  
- Biopython  
- Pairwise sequence alignment  
- FASTA genome data from NCBI  

---

