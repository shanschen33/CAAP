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
Don't worry. we will generate the `codeml.ctl` later in this pipeline.

## 2. Evolutionary Parameter Inference

Directory: 02_ancestal
- Process:
  - Utilize `codeml` to infer evolutionary parameters.
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
# We use GS1 here as an example.
python 04_scripts/02_calc_R_group_aas.py 02_ancestral 03_output GS1
```
- Output:
Convergence results are saved in the `03_output` directory. Each gene will have two files with names formatted as `{gene_name}.{group_scheme}.expconv.tsv` and `{gene_name}.{group_scheme}.obsconv.tsv`

The file with the `expconv.tsv` suffix contains the expected number of CAAP sites calculated under neutral evolution. The first, second and third column indicates tested branch pairs, expected parallel CAAP sites and expected total CAAP sites (including parallel substitutions and convergent substitutions), respectively.

The file with the `obsconv.tsv` suffix contains the expected number of CAAP sites calculated under neutral evolution. The first, second, third, fourth and fifth column indicates tested branch pairs, parallel CAAP sites, convergent CAAP sites, list of parallel CAAP sites and list of convergent CAAP sites, respectively.

Please see the details in the CAAP paper.

# citation
1. R
2. CCS
3. Si Quang Le, Olivier Gascuel, An Improved General Amino Acid Replacement Matrix, Molecular Biology and Evolution, Volume 25, Issue 7, July 2008, Pages 1307–1320, https://doi.org/10.1093/molbev/msn067


