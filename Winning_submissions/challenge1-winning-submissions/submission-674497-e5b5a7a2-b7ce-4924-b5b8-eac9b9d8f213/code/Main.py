'''
--------------------------------------------------------------------------------
TCells.src.Main Main module
--------------------------------------------------------------------------------
Cancer Immunotherapy Data Science Grand Challenge - Challenge 1
--------------------------------------------------------------------------------
Brody Langille, Jordan Trajkovski
Omega Funds
2023
--------------------------------------------------------------------------------
'''

import logging, os
from typing import Tuple, List, Dict, Union
from pandas import DataFrame, concat
from numpy import array
from tensorflow.keras.models import Sequential
from sklearn.model_selection import train_test_split
from TCellData import TCellData, TCellGeneExpression, TCellEndStates, closestVector
from TCellModel import build, train, predict, plotModelTrainingHistory

config: dict = {
    "logLevel": "INFO",
    "data": {
        "training": {
            "path": "sc_training.h5ad",
            "testPercent": 0.05,
            "shuffle": True,
            "averageExpressionPerCondition": False
        },
        "inference": {
            "validation": {
                "path": "../solution/validation_output.csv",
                "knockoutGenes": [
                    "Aqr",
                    "Bach2",
                    "Bhlhe40"
                ]
            },
            "test": {
                "path": "../solution/test_output.csv",
                "knockoutGenes": [
                    "Ets1",
                    "Fosb",
                    "Mafk",
                    "Stat3"
                ]
            }
        },
        "labelMap": [
            {
                "state": "progenitor",
                "letter": "a"
            },
            {
                "state": "effector",
                "letter": "b"
            },
            {
                "state": "terminal exhausted",
                "letter": "c"
            },
            {
                "state": "cycling",
                "letter": "d"
            },
            {
                "state": "other",
                "letter": "e"
            }
        ]
    },
    "model": {
        "vectorLabels": True,
        "epochs": 1000,
        "batchSize": 20,
        "learningRate": 0.000005,
        "regularizationRate": 0.001,
        "dropoutRate": 0.001,
        "activation": "tanh",
        "inputLayerUnits": 14,
        "metrics": [
            "accuracy"
        ],
        "showPlot": False
    },
    "modelTrainingIterations": 25,
    "useMedianClosestVector": True
}

'''
MAIN
'''
def main(config: dict) -> array:
    logging.getLogger().setLevel(config["logLevel"])

    # define general objects used in both load/build (model) cases below
    tCellData: TCellData = TCellData(
        path=config["data"]["training"]["path"],
        labelMap=config["data"]["labelMap"]
    )

    # get (all) raw data from .h5ad file
    # as a single, joined dataframe
    df: DataFrame = tCellData.getRawDataFrame()

    # create a (T cell) gene expression object
    # which loads and calculates the expression vectors for all genes
    # and can return a vector for any gene (in the future)
    tCellGeneExpression: TCellGeneExpression = TCellGeneExpression(df=df, mean=config["data"]["training"]["averageExpressionPerCondition"])

    # create a (T cell) end state object
    # which loads and calculates the end state probability distributions for
    # the knocked out genes in the input data only
    tCellEndStates: TCellEndStates = TCellEndStates(df=df)

    # create object to store outputs from each iteration
    # for each output type [test, validation]
    predictGenes: List[str] = config["data"]["inference"]["test"]["knockoutGenes"] + config["data"]["inference"]["validation"]["knockoutGenes"]
    predictions: Dict[str, List[array]] = {gene: [] for gene in predictGenes}
    outputHeader: List[str] = ["gene"] + [f'{x["letter"]}_i' for x in config["data"]["labelMap"]] # output .csv file header i.e. "gene,a_i,b_i,c_i,d_i,e_i"

    '''
    compute the label (y) data
    '''
    y: array = tCellEndStates.getTCellEndStateVectors()
    logging.info(f'y (shape: {y.shape}):\n\n{y}')

    '''
    compute the input/feature (x) data
    '''
    # get the conditions (in the same order as they appear in the label/y data)
    # to use to get the corresponding gene expression vectors, per gene, for the input/feature data
    conditions: List[str] = tCellEndStates.getConditions()

    # get the input/feature vectors
    # that match up line-by-line to the conditions in the labels/y data
    x: List[array] = []

    for condition in conditions:
        expression: array = tCellGeneExpression.getTCellGeneExpressionVector(gene=condition)
        x.append(expression)

    x: array = array(x)
    logging.info(f'x (shape: {x.shape}):\n\n{x}')

    '''
    main loop

        - re-build, train ML model on each iteration (using the same input/output data and model config/hyperparameters)
        - compute and store/aggregate predictions for each gene
    '''

    # train a model multiple times. count defined by the loopCount config value
    for i in range(0, config["modelTrainingIterations"]):

        '''
        split the x and y sets (input vectors + labels)
        into training and test sets
        do this split for each modelTrainingIterations
        '''
        xTrain, xTest = train_test_split(
            x,
            test_size=config["data"]["training"]["testPercent"],
            shuffle=config["data"]["training"]["shuffle"]
        )
        logging.info(f'xTrain (shape: {xTrain.shape}):\n\n{xTrain}')
        logging.info(f'xTest (shape: {xTest.shape}):\n\n{xTest}')
        
        yTrain, yTest = train_test_split(
            y,
            test_size=config["data"]["training"]["testPercent"],
            shuffle=config["data"]["training"]["shuffle"]
        )
        logging.info(f'yTrain (shape: {yTrain.shape}):\n\n{yTrain}')
        logging.info(f'yTest (shape: {yTest.shape}):\n\n{yTest}')

        '''
        build the ML model
        '''
        model: Sequential = build(
            inputLayerUnits=config["model"]["inputLayerUnits"],
            inputSize=xTrain.shape[1],
            outputSize=yTrain.shape[1],
            metrics=config["model"]["metrics"],
            learningRate=config["model"]["learningRate"],
            vectorLabels=config["model"]["vectorLabels"],
            regularizationRate=config["model"]["regularizationRate"],
            dropoutRate=config["model"]["dropoutRate"],
            activation=config["model"]["activation"]
        )

        '''
        train the ML model
        '''
        model, history = train(
            model=model,
            x=xTrain,
            y=yTrain,
            epochs=config["model"]["epochs"],
            batchSize=config["model"]["batchSize"]
        )

        # plot the ML model training (metrics) curve(s)
        if config["model"]["showPlot"] == True:
            plotModelTrainingHistory(
                history=history,
                metricsList=config["model"]["metrics"],
            )

        '''
        predict the effect of each gene knocked out,
        on the T cell end state distribution
        '''
        for gene in predictGenes:
            # get the inference data
            xInference: array = tCellGeneExpression.getTCellGeneExpressionVector(gene=gene)
            logging.debug(f'Knockout gene: {gene} inference input: {xInference}')

            # make prediction using the trained model
            yInference: array = predict(
                model=model,
                x=xInference
            )

            logging.debug(f'Knockout gene: {gene} prediction:\n\n{yInference}')
            predictions[gene].append(yInference)

    '''
    after all (model) iterations (and predictions for the knockout genes) are completed,
    perform the final results computations

    and generate the (2) output (files):

        1. solution/validation_output.csv
            * 3 genes knocked out
        2. solution/test_output.csv
            * 4 genes knocked out
    '''
    logging.info(f'Final processing after {config["modelTrainingIterations"]} iterations, on results:\n\n{predictions}...')

    # ech of test, validation (sets of knockout genes, and separate output files)
    for typ, data in config["data"]["inference"].items():
        path: str = data["path"]
        folder: str = os.path.dirname(path)

        if not os.path.exists(folder):
            os.makedirs(folder)

        knockoutGenes: List[str] = data["knockoutGenes"]
        logging.info(f'Processing final results for {typ} ({len(knockoutGenes)} knockout genes: {knockoutGenes})...')
        results: List[Union[str, float]] = []

        for gene in knockoutGenes:
            r: array = closestVector(x=array(predictions[gene]), useMedian=config["useMedianClosestVector"])
            results.append([gene] + r.tolist())

        df: DataFrame = DataFrame(results, columns=outputHeader)
        logging.info(f'Final {typ} results:\n\n{df}')
        df.to_csv(path, index=False)

if __name__ == "__main__":
    main(config=config)
