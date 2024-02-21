# calcon-aas
Calculated convergence after grouping amino acids by physicochemical properties, modified origin method in Zou &amp; Zhang, 2015 and Xu et al., 2017.

# modified R
*01_data* contained the LG substitution matrix used to infer ancestral, 6 cleaned amino acid fasta file and respective tree file.

*02_ancestal* contained the inferred parameter of ancestors, attained by running *04_scripts/01_run_codeml*.

*03_output* contained the output of sample data, attained by run `python 04_scripts/02_calc_R_group_aas.py 02_ancestral 03_output GS1` (Here take GS1 as an example).

*04_scripts* contained scripts.

# modified CCS
*01_data* contained 10 cleaned amino acid fasta file.

*02_output* contained the output of the data, attained by running *03_scripts/01_CCS_group_aas.py*.
