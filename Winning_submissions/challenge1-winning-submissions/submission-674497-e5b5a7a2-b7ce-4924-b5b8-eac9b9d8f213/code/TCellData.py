'''
--------------------------------------------------------------------------------
TCells.src.TCellData T Cell Data Specific Class
--------------------------------------------------------------------------------
Loading and processing raw T cell and gene expression data
--------------------------------------------------------------------------------
Brody Langille
Omega Funds
2023
--------------------------------------------------------------------------------
'''

import logging, traceback
from typing import Tuple, List, Dict, Union
from scanpy import read_h5ad, AnnData
from pandas import DataFrame, Index, concat, Series, set_option
from scipy.sparse._csr import csr_matrix
from numpy import array, mean, median, argmin
from numpy.linalg import norm as euclideanDistance

def closestVector(x: array, useMedian: bool = True) -> array:
    '''
    determine the vector (row) from an input (2D) matrix
    that is closest to the (element-wise columnular) mean or median of all vectors in the matrix
    i.e. the mean or median vector is a single vector/row that contains the mean or median of each "column"
    '''
    assert(len(x.shape) == 2)
    logging.debug(f'Computing the closest vector to the mean/median (vector) from input matrix:\n\n{x}')

    # calculate the mean or median of each column in the input matrix,
    # store into a single vector/row
    if useMedian == True:
        metric: str = "Median"
        referenceVector: array = array([median(x[:, i]) for i in range(0, x.shape[1])]) # x.shape[1] is number of "columns" in the input matrix
    else:
        metric: str = "Mean"
        logging.info(f'Received argument of useMedian == False. Computing mean vector now...')
        referenceVector: array = array([mean(x[:, i]) for i in range(0, x.shape[1])]) # x.shape[1] is number of "columns" in the input matrix
    
    logging.debug(f'{metric} vector: {referenceVector}')

    # calculate the Euclidean distance between each vector/row and the calculated mean vector/row
    distancesFromReferenceVector: array = array([euclideanDistance(vector - referenceVector) for vector in x])
    logging.debug(f'Vector (Euclidean) distances from the {metric} vector: {distancesFromReferenceVector}')

    i: int = argmin(distancesFromReferenceVector) # find the min distance (index)
    r: array = x[i] # get the min distance vector/row based on index
    logging.debug(f'Vector (i = {i}) is the closest to the {metric} (vector): {r}')
    return r

class TCellGeneExpression:
    def __init__(self, df: DataFrame, mean: bool):
        self.data: DataFrame = self._getTCellGeneExpression(df=df, mean=mean)

    def getTCellGeneExpressionVector(self, gene: str) -> array:
        '''
        return a vector of gene expression values
        measured across all experiments
        i.e. against all (66) knocked out genes
        '''
        try:
            r: array = self.data[gene].to_numpy()
            logging.debug(f'Expression vector for gene {gene}: {r}')
            return r

        except Exception as e:
            msg: str = f'Failed to find expression vector for gene {gene}. Due to exception: {e}. Traceback: {traceback.format_exc()}'
            logging.error(msg)
            raise Exception(msg)

    @staticmethod
    def _getTCellGeneExpression(df: DataFrame, mean: bool = True, write2File: bool = False) -> DataFrame:
        '''
        get the gene (RNA) expression values for all measured genes
        either:
            1. per condition/knockout gene (mean expression) - if enabled (default) - results in vectors of length 66, or
            2. per sample (T cell) - results in vectors of length ~29k

        example:
                      condition    Mrpl15    Lypla1     Tcea1   Atp6v1h  ...  9130016M20Rik      Htr7   Col17a1     Awat2      Amot
        0        Arid4b  0.380150  0.242047  0.583878  0.124692  ...       0.000000  0.000858  0.000000  0.000865  0.000000
        1        Arid5b  0.433201  0.249483  0.582998  0.131191  ...       0.000000  0.000000  0.000000  0.000000  0.000000
        2          Atf2  0.329280  0.235661  0.541153  0.145771  ...       0.000000  0.000000  0.000000  0.000000  0.000000
        3          Batf  0.357733  0.155210  0.515990  0.066413  ...       0.000000  0.000000  0.000000  0.000000  0.000000
        4          Crem  0.284578  0.224179  0.515363  0.114598  ...       0.000000  0.000000  0.000996  0.000000  0.000000
        5        Ctnnb1  0.299226  0.183942  0.604708  0.185781  ...       0.000000  0.000000  0.000000  0.000000  0.000000
        6          Dkk3  0.342953  0.256452  0.586786  0.169227  ...       0.000000  0.000000  0.000000  0.000000  0.000000
        7          Dvl1  0.348138  0.236937  0.557947  0.123965  ...       0.000000  0.000000  0.000000  0.000000  0.000000
        8          Dvl2  0.418156  0.268286  0.596147  0.144550  ...       0.000470  0.000283  0.000000  0.000425  0.000000
        9          Dvl3  0.365085  0.219149  0.566919  0.092630  ...       0.000000  0.000000  0.000000  0.000000  0.000000
        10         Eef2  0.277175  0.000000  1.014028  0.000000  ...       0.000000  0.000000  0.000000  0.000000  0.000000
        11         Egr1  0.346773  0.236045  0.526344  0.135105  ...       0.001441  0.000000  0.000000  0.000000  0.000000
        12         Elf1  0.277793  0.197251  0.434838  0.114430  ...       0.000000  0.000000  0.000000  0.000000  0.000000
        '''
        logging.info(f'Generating T cell gene expression...')

        x: DataFrame = df.copy() # dataframe is passed by reference, so make a copy, change it and return it

        # remove the (T cell end) state column
        # leaving only the condition and gene columns
        x.drop(columns=["state"], inplace=True)

        # calculate the mean expression value for each condition/knockout gene as a row
        if mean == True:
            x: DataFrame = x.groupby(by="condition", as_index=False).mean()

        set_option('display.max_rows', None)
        logging.info(f'T cell gene expression (mean = {mean}):\n\n{x}')

        if write2File == True:
            path: str = "data/t-cell-gene-expression-per-knockout-gene.csv"
            x.to_csv(path, index=False)
            logging.info(f'Wrote T cell gene expression per knockout gene to file: {path}')

        return x

class TCellEndStates:
    def __init__(self, df: DataFrame):
        self.data: DataFrame = self._getTCellEndStateDistributionsPerKnockoutGene(df=df)

    def getTCellEndStateVectors(self) -> array:
        # convert to labels (numpy arrays)
        r: array = self.data.copy().drop(columns=["condition", "nTCells"]).to_numpy()
        logging.debug(f'Converted T cell end states dataframe -> numpy array: {r}')
        return r

    def getConditions(self) -> List[str]:
        r: List[str] = self.data["condition"].tolist()
        logging.debug(f'{len(r)} conditions: {r}')
        return r

    @staticmethod
    def _getTCellEndStateDistributionsPerKnockoutGene(df: DataFrame, write2File: bool = False) -> DataFrame:
        '''
        get the distribution of T cell end states
        per knockout gene
        (not including unperturbed state and genes missing from expression data due to QC filtering)
        i.e. a dataframe with 66 rows and 7 columns [gene, nTCells, progenitor, effector, terminal exhausted, cycling, other]
        example:
                        condition  nTCells  progenitor  effector  terminal exhausted   cycling     other
        0       Rps6        2    1.000000  0.000000            0.000000  0.000000  0.000000
        1      Ep300       24    1.000000  0.000000            0.000000  0.000000  0.000000
        2      Tbx21        2    1.000000  0.000000            0.000000  0.000000  0.000000
        3       Eef2        1    1.000000  0.000000            0.000000  0.000000  0.000000
        4       Klf2        5    1.000000  0.000000            0.000000  0.000000  0.000000
        5      Runx3       10    1.000000  0.000000            0.000000  0.000000  0.000000
        6       Batf        6    0.666667  0.000000            0.166667  0.000000  0.166667
        7      Foxm1       56    0.642857  0.071429            0.142857  0.142857  0.000000
        8        Yy1       25    0.600000  0.200000            0.080000  0.120000  0.000000
        9     Ctnnb1       43    0.534884  0.209302            0.139535  0.116279  0.000000
        10     Sp140       33    0.515152  0.121212            0.121212  0.242424  0.000000
        11      Tpt1       25    0.440000  0.160000            0.120000  0.280000  0.000000
        12     Eomes       16    0.437500  0.312500            0.062500  0.187500  0.000000
        ...
        '''
        logging.info(f'Generating T cell end state distributions per knockout gene...')

        # drop all of the (extra) gene expression data that we don't need for this method
        # only retain (T cell end) state and condition (knockout genes)
        x: DataFrame = df.copy()
        x: DataFrame = x[["state", "condition"]]

        # get a unique list of knockout genes (incl. the Unperturbed label)
        conditions: List[str] = x["condition"].unique().tolist()
        logging.debug(f'{len(conditions)} unique conditions (knockout genes + unperturbed): {conditions}')

        excludeConditions: List[str] = ["Unperturbed", "Fzd1", "P2rx7"] # latter 2 genes have no expression data (QC filter)
        states: List[str] = ["progenitor", "effector", "terminal exhausted", "cycling", "other"] # T cell measured states NOTE: in order relative to the probability distribution P defined in the challenge requirements
        columns: List[str] = ["condition", "nTCells"] + states # columns to reconstruct a final dataframe
        probabilities: List[List[Union[str, float]]] = [] # tracking of probability distributions across states, per knockout gene

        for condition in conditions:
            if condition not in excludeConditions:
                # filter main dataframe to only the current condition
                # i.e. knockout gene
                r: DataFrame = x[x["condition"] == condition]
                n: int = len(r)
                logging.debug(f'Condition {condition} ({n} samples):\n\n{r}')

                # get the distribution across states for this specific condition
                # i.e. knockout gene
                # name the values column as "probability" and the index of the dataframe as the (T cell end) state
                # for easier processing downstream
                dist: DataFrame = r[["state"]].value_counts(normalize=True).to_frame("probability").reset_index("state")
                logging.debug(f'Condition {condition} T cell end state distribution:\n\n{dist}')

                # map the probability values for each state to a list
                # (where the indices match the states list, in order)
                probs: List[float] = [dist[dist["state"] == state]["probability"].item() for state in states]
                probs: List[Union[str, float]] = [condition, n] + probs # add the knockout gene and count of samples to the list/row as well
                probabilities.append(probs)

            else:
                logging.debug(f'Excluding condition {condition} from output.')

        # convert the results into a single, consolidated dataframe
        r: DataFrame = DataFrame(probabilities, columns=columns)
        # sort the results such that the conditions/knockout genes with
        #  the highest % of T cells in the end state "progenitor" appear
        # at the top and the list is descending on this metric
        r.sort_values(by=["progenitor"], ascending=False, inplace=True)
        r.reset_index(drop=True, inplace=True) # remove index post-sort
        set_option('display.max_rows', None) # show all rows of the dataframe when printed
        logging.info(f'T cell end state distributions per knockout gene:\n\n{r}')

        if write2File == True:
            path: str = "data/t-cell-end-state-distributions-per-knockout-gene.csv"
            r.to_csv(path, index=False)
            logging.info(f'Wrote T cell end state distributions per knockout gene to file: {path}')

        return r

class TCellData:
    def __init__(self, path: str, labelMap: List[dict]):
        self.path: str = path
        self.labelMap: List[dict] = labelMap

    def _getRawData(self) -> Tuple[csr_matrix, DataFrame, Index, Index]:
        '''
        get the provided, raw data from a single .h5ad file
        and return the consituent components, separately, for further processing
        '''
        logging.info(f'Getting raw data from .h5ad file: {self.path}...')

        adata: AnnData = read_h5ad(self.path)
        logging.debug(f'adata ({type(adata)}):\n\n{adata}')

        # spare, csr_matrix, data matrix
        X: csr_matrix = adata.X
        logging.debug(f'adata.X ({type(X)}, shape: {X.shape}):\n\n{X}')

        # metadata, dataframe, related to data matrix
        obs: DataFrame = adata.obs
        logging.debug(f'adata.obs ({type(obs)}):\n\n{obs}')

        # column names/index of X csr_matrix
        varNames: Index = adata.var_names
        logging.debug(f'adata.var_name ({type(varNames)}):\n\n{varNames}')

        # row names/index of X csr_matrix
        obsNames: Index = adata.obs_names
        logging.debug(f'adata.obsNames ({type(obsNames)}):\n\n{obsNames}')

        logging.info(f'X, obs, varNames, obsNames data retrieved from .h5ad file: {self.path}.')
        return X, obs, varNames, obsNames

    def getRawDataFrame(self) -> DataFrame:
        '''
        combine the raw data from the .h5ad data file
        into a single dataframe that can be used by downstream processes
        '''
        logging.info(f'Generating a single dataframe from raw data...')

        # get raw data components, read and unpacked from .h5ad file
        r: Tuple[csr_matrix, DataFrame, Index, Index] = self._getRawData()

        # reconstruct into a single dataframe
        # with rows as T cells/samples gene expression and columns as genes (that expression is measured for)
        df: DataFrame = DataFrame(r[0].toarray(), columns=r[2], index=r[3])

        # also add the knocked out genes (condition) and the T cell end state (ex. progenitor) as columns as well
        x: DataFrame = r[1][["state", "condition"]]

        # confirm that the 2 dataframes natch up
        # row-by-row in their indices before merging/joining and returning
        for i in range(0, len(x.index)):
            assert(x.index[i] == df.index[i])

        r: DataFrame = concat([df, x], axis="columns")
        logging.info(f'Generated a single dataframe from the raw data:\n\n{r}')
        return r
