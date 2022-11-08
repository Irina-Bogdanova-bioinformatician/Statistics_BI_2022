import argparse
from statsmodels.stats.weightstats import ztest
from pathlib import Path
import scipy.stats as st
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests


def check_intervals_intersect(first_ci: tuple, second_ci: tuple) -> bool:
    """ Check if intervals intersect (True or False)
    :param first_ci: confidence intervals 1
    :param second_ci: confidence intervals 2 to compare with 1
    :return: True or False
    """
    are_intersect = False
    intersection = min(first_ci[1], second_ci[1]) - max(first_ci[0], second_ci[0])
    if intersection >= 0:
        are_intersect = True
    return are_intersect


def check_dge_with_ci(first_table: pd.DataFrame, second_table: pd.DataFrame, genes: list) -> list:
    """ Check differential gene expression using confidence intervals.
        The function returns list with expression comparison results for each gene (order as in 'genes' list)
        (True if there is a difference)
    :param first_table: data for one of the compared cell types
    :param second_table: data for another one of the compared cell types
    :param genes: list with genes
    :return: ci_test_results [True/False, ...]
    """
    ci_test_results = []

    def get_ci(data):
        ci = st.t.interval(alpha=0.95, df=len(data) - 1, loc=np.mean(data),
                           scale=st.sem(data))
        return ci

    for gene in genes:
        ci_1 = get_ci(first_table[gene])
        ci_2 = get_ci(second_table[gene])
        ci_test_results.append(not check_intervals_intersect(ci_1, ci_2))

    return ci_test_results


def check_dge_with_ztest(first_table: pd.DataFrame, second_table: pd.DataFrame, genes: list, correction_method=None) \
        -> (list, list):
    """ Check differential gene expression using z test.
        The function returns list with expression comparison results for each gene (order as in 'genes' list)
        (True if there is a difference), as well as list of p-values for each gene
        :param first_table: data for one of the compared cell types
        :param second_table: data for another one of the compared cell types
        :param genes: list with genes
        :param correction_method: statsmodels method of multiple comparisons correction - not required.
        :return: z_test_results [True/False, ...], p_values list
    """
    z_test_results = []
    p_values = []

    for gene in genes:
        p = ztest(first_table[gene], second_table[gene])[1]
        p_values.append(p)

    if correction_method:
        p_values = multipletests(p_values, alpha=0.05, method=correction_method, is_sorted=False, returnsorted=False)[1]

    z_test_results = []
    for p in p_values:
        result = True if p < 0.05 else False
        z_test_results.append(result)

    return z_test_results, p_values


def get_mean_diff(first_table: pd.DataFrame, second_table: pd.DataFrame, genes: list) -> list:
    """ Count difference in mean gene expressions.
    :param first_table: data for one of the compared cell types
    :param second_table: data for another one of the compared cell types
    :param genes: genes list
    :return: list with mean expression differences
    """
    mean_diff = []

    for gene in genes:
        mean_diff.append(round(first_table[gene].mean() - second_table[gene].mean(), 3))
    return mean_diff


def main(expression_path_1: Path, expression_path_2: Path, out_path: Path, correction_method: str):
    expression_data_1 = pd.read_csv(expression_path_1, index_col=0)
    expression_data_2 = pd.read_csv(expression_path_2, index_col=0)

    methods = {'bonferroni', 'sidak', 'holm-sidak', 'holm', 'simes-hochberg', 'hommel', 'fdr_bh',
               'fdr_by', 'fdr_tsbh', 'fdr_tsbky'}
    if correction_method and correction_method not in methods:
        print(f'There is no {correction_method} multiple comparisons correction method in statsmodels. '
              f'No correction will be carried out')

    genes_1 = list(expression_data_1.columns)
    genes_2 = list(expression_data_2.columns)
    if genes_1 != genes_2:
        print("The headers in compared columns must be identical")
    else:
        if "Cell_type" in genes_1:
            genes_1.remove("Cell_type")

        ci_test_result = check_dge_with_ci(expression_data_1, expression_data_2, genes_1)
        ztest_result, p_values = check_dge_with_ztest(expression_data_1, expression_data_2, genes_1, correction_method)
        mean_diff = get_mean_diff(expression_data_1, expression_data_2, genes_1)

        results = {
            "ci_test_results": ci_test_result,
            "z_test_results": ztest_result,
            "z_test_p_values": p_values,
            "mean_diff": mean_diff
        }

        results = pd.DataFrame(results)
        results.index = genes_1
        results = results.round({'z_test_p_values': 3})

        results.to_csv(out_path)
        print("We have finished there")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generating the expression difference table. "
                                                 "The input tables must include expressions of the same genes")
    parser.add_argument("--first_cell_type_expressions_path", required=True, type=Path,
                        help="Path to the table (csv file) with gene expressions of first cell type")
    parser.add_argument("--second_cell_type_expressions_path", required=True, type=Path,
                        help="Path to the table with gene expressions of second cell type")
    parser.add_argument("--save_results_table", type=Path, default="expression_comparison_results.csv",
                        help="Path to the output table with gene expressions comparison results")
    parser.add_argument('--correction_method', help='The name of the method for correcting for multiple comparisons '
                                                    'implemented in statsmodels', type=str, default=None)
    args = parser.parse_args()

    main(args.first_cell_type_expressions_path, args.second_cell_type_expressions_path,
         args.save_results_table, args.correction_method)

