import os
import pandas as pd
import sys
import scipy.stats as stats

CAAP_results_dir = sys.argv[1]
output_dir = sys.argv[2]

def main():
    file_list_O = [x for x in os.listdir(CAAP_results_dir) if x.endswith('obsconv.tsv')]
    file_list_E = [x for x in os.listdir(CAAP_results_dir) if x.endswith('expconv.tsv')]

    result_O = pd.DataFrame(columns=["gene", "O_p", "O_c", "brp"])
    result_E = pd.DataFrame(columns=["gene", "E_pc", "brp"])

    for file in file_list_O:
        df = pd.read_csv(f'{CAAP_results_dir}/{file}', sep="\t", header=None)
        gene = file.split(".")[0]
        temp_df = pd.DataFrame({
            "gene": [gene] * len(df),
            "O_p": df[1],
            "O_c": df[2]
        })
        result_O = pd.concat([result_O, temp_df], ignore_index=True)
    result_O["brp"] = [1, 2] * (len(result_O) // 2)

    for file in file_list_E:
        df = pd.read_csv(f'{CAAP_results_dir}/{file}', sep="\t", header=None)
        gene = file.split(".")[0]
        temp_df = pd.DataFrame({
            "gene": [gene] * len(df),
            "E_pc": df[2]
        })
        result_E = pd.concat([result_E, temp_df], ignore_index=True)
    result_E["brp"] = [1, 2] * (len(result_E) // 2)

    combine_OE = pd.merge(result_O, result_E, on=["gene", "brp"])

    combine_OE["O_pc"] = combine_OE["O_p"] + combine_OE["O_c"]
    combine_OE["R"] = combine_OE["O_pc"] / combine_OE["E_pc"]

    combine_OE['poisson_p_value'] = combine_OE.apply(
    lambda row: poisson_test(row['O_pc'], row['E_pc']), axis=1)

    output_file = os.path.join(output_dir, "R_value_significance.tsv")
    combine_OE.to_csv(output_file, sep='\t', index=False)

    print(f"Results saved to: {output_file}")

def poisson_test(O, E):
    if O < E:
        p = stats.poisson.cdf(O, E)
    else:
        p = 1 - stats.poisson.cdf(O - 1, E)
    return p

if __name__ == "__main__":
    main()
