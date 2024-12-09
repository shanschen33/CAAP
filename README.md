# CAAP Pipeline
CAAP (Convergence of Amino Acid Property) detects convergence by grouping amino acids based on their physicochemical properties. CAAP can be implemented in different convergence detection methods. Here, we recommend using the R method (`Zou & Zhang in 2015`) as the primary approach, with the CCS method (`Xu et al. in 2017`) serving as an alternative for detecting CAAP.

# Pipeline Overview
The CAAP pipeline is organized into sequential steps, including data preparation, ancestral parameter inference, and convergence calculation using the preferred R method or the alternative CCS method. Below is a detailed breakdown of each stage in the pipeline:
```
CAAP Pipeline
│
├── 01_data
│   ├── LG substitution matrix
│   ├── Cleaned amino acid FASTA file
│   └── Respective tree file
│
├── 02_ancestal
│   └── Inferred ancestral parameters
│
├── 03_output
│   └── Convergence calculation results (R method)
│
├── 04_scripts
│   ├── 01_gnrt_dic_ctf.py
│   └── 02_calc_R_group_aas.py
│
└── 05_CCS
    ├── 01_data
    │   └── Cleaned amino acid FASTA file
    ├── 02_output
    │   └── Convergence calculation results (CCS method)
    └── 03_scripts
        └── 03_CCS_group_aas.py
```
# Step-by-Step Pipeline
## 1. Data and Dependencies Preparation
### Data Preparation
See example files in `01_data`.
- Amino acid substitution matrix: We used LG here, other amino acid substitution matrices are also acceptable. LG is an empirical amino acid substitution matrix estimated using 3,912 alignments from Pfam. See details in `Le & Gascuel, 2008`.
- Cleaned amino acid fasta files: You should prepare an amino acid alignment without gaps and ambiguous sites.
- Tree files: You should prepare a gene tree pruned from the species tree corresponding to the alignment. Branch lengths of trees are not needed as we will infer them later.

### PAML Installation
Please refer to this page: https://github.com/abacus-gene/paml/wiki/Installation.

After installing PAML, type `codeml` in the command line to see if this function works. If you installed it successfully, the output will be:
```
error when opening file codeml.ctl
tell me the full path-name of the file?
```
Don't worry. We will generate the `codeml.ctl` later in this pipeline.

## 2. Evolutionary Parameter Inference
- Process:
Generate control file for codeml and infer evolutionary parameters using PAML.

- Execution:
```
# This will create a dictionary named by the gene. The dictionary contains the codeml.ctl, which is the control file to run codeml in PAML.
python 04_scripts/01_gnrt_dic_ctf.py 01_data/ 02_ancestral/ 01_data/lg.dat
cd 02_ancestral/ENSG00000000003.AA.cleaned
# Use codeml in PAML. This will allow you to have rst and rates files, which are essential for following convergence detection.
codeml codeml.ctl
```
- Output:
Inferred ancestral parameters are stored in the `02_ancestral/ENSG00000000003.AA.cleaned` directory.

## 3. Convergence Calculation (R)
- Process:
Calculate CAAP based on R using evolutionary parameters. The choice of group schemes (GS) of amino acids is flexible (GS1-4), we recommend trying all the group schemes to take different aspects of physiochemical properties into consideration.
- Execution:
```
# We use GS1 here as an example. This script counts observed CAAP sites (O) and calculates the expected CAAP sites (E).
python 04_scripts/02_calc_R_group_aas.py 02_ancestral 03_output GS1
# This script summarizes the statistics, calculates R, and tests whether O significantly deviated from E. We select R<0.05 as the threshold for candidate genes (See details in CAAP paper).
python 04_scripts/03_sig_gene.py 03_output/ 03_output/
```
- Output:

Results are saved in the `03_output` directory. `02_calc_R_group_aas.py` will generate two files with names formatted as `{gene_name}.{group_scheme}.expconv.tsv` and `{gene_name}.{group_scheme}.obsconv.tsv` for each gene. 

The file with the `expconv.tsv` suffix contains the expected number of CAAP sites calculated under neutral evolution. The first, second and third column indicates tested branch pairs, expected parallel CAAP sites and expected total CAAP sites (including parallel substitutions and convergent substitutions), respectively.

The file with the `obsconv.tsv` suffix contains the expected number of CAAP sites calculated under neutral evolution. The first, second, third, fourth and fifth column indicates tested branch pairs, parallel CAAP sites, convergent CAAP sites, list of parallel CAAP sites and list of convergent CAAP sites, respectively.

`03_sig_gene.py` will generate a file named `R_value_significance.tsv`.

Please see the details in the CAAP paper.

## 4. Based on CCS
```
cd 05_CCS
```

To detect genes with CAAP in CCS method, you only need to prepare the amino acid alignment without gaps and ambiguous sites for each gene. Please see examples in `01_data`.

Then detect the genes with CAAP:
```
python 03_scripts/01_CCS_group_aas.py 01_data/ 02_output/ GS1
```
The output file named `{Group_scheme}.output.tsv` will be saved in `02_output`. The first, second, third, fourth and fifth column indicates names of alignments, parallel CAAP sites, convergent CAAP sites, list of parallel CAAP sites and list of convergent CAAP sites, respectively.

# References
1. Zhengting Zou, Jianzhi Zhang, Are Convergent and Parallel Amino Acid Substitutions in Protein Evolution More Prevalent Than Neutral Expectations?, Molecular Biology and Evolution, Volume 32, Issue 8, August 2015, Pages 2085–2096, https://doi.org/10.1093/molbev/msv091
2. Shaohua Xu, Ziwen He, Zixiao Guo, Zhang Zhang, Gerald J. Wyckoff, Anthony Greenberg, Chung-I. Wu, Suhua Shi, Genome-Wide Convergence during Evolution of Mangroves from Woody Plants, Molecular Biology and Evolution, Volume 34, Issue 4, April 2017, Pages 1008–1015, https://doi.org/10.1093/molbev/msw277
3. Si Quang Le, Olivier Gascuel, An Improved General Amino Acid Replacement Matrix, Molecular Biology and Evolution, Volume 25, Issue 7, July 2008, Pages 1307–1320, https://doi.org/10.1093/molbev/msn067


