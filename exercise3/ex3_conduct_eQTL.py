import pandas as pd
import numpy as np
from scipy.stats import f
import multiprocessing
from concurrent.futures import ProcessPoolExecutor, as_completed


def SXX(X):
    """Calculate the sum of squares of X"""
    return np.sum(X ** 2) - len(X) * X.mean() ** 2


def SXY(X, Y):
    """Calculate the sum of squares of X and Y"""
    return np.sum(X * Y) - len(X) * X.mean() * Y.mean()


def calc_regression(X, Y):
    """Calculate the regression coefficients"""
    b1 = SXY(X, Y) / SXX(X)
    b0 = Y.mean() - b1 * X.mean()

    return b0, b1


def SSE(X, Y):
    """Calculate the sum of squared errors"""
    b0, b1 = calc_regression(X, Y)
    y_hat = lambda x: b0 + x * b1
    sse = np.sum((Y - y_hat(X)) ** 2)

    return sse, sse / (len(X) - 2)


def SSR(X, Y):
    """Calculate the sum of squared regression"""
    b0, b1 = calc_regression(X, Y)
    y_hat = lambda x: b0 + x * b1

    ssr = np.sum((Y.mean() - y_hat(X)) ** 2)

    return ssr, ssr


def SST(X, Y):
    """Calculate the total sum of squares"""

    sst = np.sum((Y.mean() - Y) ** 2)
    return sst, sst / (len(X) - 1)


def calc_corr(X, Y):
    """Calculate the correlation coefficient"""
    ssr, msr = SSR(X, Y)
    sst = SST(X, Y)[0]
    mse = SSE(X, Y)[1]

    R = ssr / sst
    F = msr / mse

    return R, F, f.sf(F, 1, len(X) - 2)


def convert_to_int(x, mode=0):
    """
    Converts genotype data from string to int
    :param x: data in specific cell in the data frame
    :param mode: 0- for Linear regression without Hetrozygotes, 1- for Linear regression with Hetrozygotes, 2- for ANOVA without Hetrozygotes
    """

    if mode == 0:
        d = {"B": 0, "D": 1, "H": np.NaN, "U": np.NaN, 'b': 0, 'd': 1, 'h': np.NaN, 'u': np.nan}
    elif mode == 1:
        d = {"B": 0, "D": 2, "H": 1, "U": np.NaN, 'b': 0, 'd': 2, 'h': 1, 'u': np.nan}

    return d[x]


def snp_gene_association(snp, gene, input_matrix, genos, mode=0):
    """Calculate the association between a SNP and a gene"""

    snp_genotypes = genos.loc[snp]
    gene_expression = input_matrix.loc[gene]

    if mode == 0:
        X = snp_genotypes[~snp_genotypes.isin(["U", "H", "u", "h"])].map(lambda x: convert_to_int(x, mode)).values
        Y = gene_expression[~snp_genotypes.isin(["U", "H", "u", "h"])].values

    elif mode == 1:
        X = snp_genotypes[~snp_genotypes.isin(["U", "u"])].map(lambda x: convert_to_int(x, mode)).values
        Y = gene_expression[~snp_genotypes.isin(["U", "u"])].values

    p = calc_corr(X, Y)[2]
    return p


def run_association_test_on_gene(gene, input_matrix, genos, par, mode=0):
    """Run the association test for a single gene"""

    print(f"Processing gene: {gene}, Progress: {par * 100:.2f}%")

    p_values = []
    for snp in genos.index:
        p = snp_gene_association(snp, gene, input_matrix, genos, mode)
        if p is not None:
            p_values.append(p)
        else:
            p_values.append(np.nan)
    return p_values, gene


def run_association_test_all_genes(matrix, genos, name, mode=0):
    """main function to run the association test for all genes"""

    # get only BXD strains that exist in the expression matrix
    common_strains = genos.columns.intersection(matrix.columns).tolist()
    input_genos = genos[common_strains].copy()

    input_matrix = np.zeros((len(matrix.index), len(genos.index)))

    gene_list = []
    futures = []

    # Run the association test for all genes in parallel
    # Use all available CPUs, but leave 2 for the system
    cpus = max(1, multiprocessing.cpu_count() - 2)
    with ProcessPoolExecutor(cpus) as executor:
        for s_ in range(len(matrix)):
            g = matrix.index[s_]
            kwargs = dict(gene=g, input_matrix=matrix, genos=input_genos, par=s_ / len(matrix), mode=mode)
            future = executor.submit(run_association_test_on_gene, **kwargs)
            futures.append(future)

    for i, future in enumerate(as_completed(futures)):
        p_vals, gene = future.result()
        input_matrix[i] = p_vals
        gene_list.append(gene)
    df = pd.DataFrame(input_matrix, index=gene_list, columns=input_genos.index)
    df.to_csv(f"output/step4_{name}_association_test_results.csv")

    return df


if __name__ == "__main__":

    # Load the processed expression matrices and unique genotypes
    processed_expression_mtxs = {
        "Erythroid": pd.read_csv("output/step3_erythroid_preprocessed_expression_data.csv", index_col=0),
        "Stem": pd.read_csv("output/step3_stem_preprocessed_expression_data.csv", index_col=0)
    }
    unique_genotypes = pd.read_csv("output/step3_unique_genotypes.csv", index_col=0)

    # Run the association tests for the erythroid and stem cells datasets
    erythroid_results = run_association_test_all_genes(
        processed_expression_mtxs['Erythroid'],
        unique_genotypes,
        "erythroid_eqtls", mode=1
    )
    stem_results = run_association_test_all_genes(
        processed_expression_mtxs['Stem'],
        unique_genotypes,
        "stem_eqtls", mode=1
    )