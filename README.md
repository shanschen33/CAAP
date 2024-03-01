# CAAP
CAAP (Convergence of Amino Acid Property) calculates convergence by grouping amino acids based on their physicochemical properties, using the modified origin method proposed by Zou & Zhang in 2015 and Xu et al. in 2017.

# modified R
*01_data* contained the LG substitution matrix used to infer ancestral, 6 cleaned amino acid fasta file and respective tree file.

*02_ancestal* contained the inferred parameter of ancestors, attained by running *04_scripts/01_run_codeml*. 

*03_output* contained the output of sample data, attained by run `python 04_scripts/02_calc_R_group_aas.py 02_ancestral 03_output GS1` (Here take GS1 as an example).

*04_scripts* contained scripts. 
- *01_run_codeml/01_gnrt_dic_ctf.py* was used to generate control file and directory for codeml.
- *01_run_codeml/02_mv_codeml.sh* and *01_run_codeml/03_run_codeml.sh* were used to ensure the parallel running of multiple codeml.

# modified CCS
*01_data* contained 10 cleaned amino acid fasta file.

*02_output* contained the output of the data, attained by running *03_scripts/01_CCS_group_aas.py*.

*03_scripts* contained scripts.
