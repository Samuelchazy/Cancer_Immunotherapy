'''
--------------------------------------------------------------------------------
TCells.src.TCellModel T Cell ML Model Class
--------------------------------------------------------------------------------
Build, train, evaluate and predict with deep ANN w/ Softmax output layer
--------------------------------------------------------------------------------
Brody Langille, Jordan Trajkovski
Omega Funds
2023
--------------------------------------------------------------------------------
'''

import logging
from typing import Tuple, List, Dict, Union
from numpy import array
from pandas import DataFrame
from matplotlib import pyplot as plt
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, Dropout
from tensorflow.keras.optimizers import Adam
from tensorflow.keras.regularizers import L1, L2
from tensorflow.keras.callbacks import History

def build(inputLayerUnits: int, inputSize: int, outputSize: int, metrics: List[str], learningRate: float = 0.001, vectorLabels: bool = True, dropoutRate: float = 0.1, regularizationRate: float = 0.0001, activation: str = "relu") -> Sequential:
    '''
    build a deep neural network
    with (multiclass) softmax output layer
    '''
    logging.info(f'Building model...')

    model: Sequential = Sequential(
        [
            Dense(units=inputLayerUnits, input_shape=(inputSize,), activation=activation, kernel_regularizer=L2(regularizationRate)),
            Dropout(rate=dropoutRate),
            Dense(units=outputSize, activation="softmax")
        ]
    )

    if vectorLabels == True:
        loss: str = "categorical_crossentropy"

    else:
        loss: str = "sparse_categorical_crossentropy"

    optimizer: Adam = Adam(
        learning_rate=learningRate
    )

    model.compile(
        optimizer=optimizer,
        loss=loss,
        metrics=metrics
    )

    logging.info(f'Built model ({type(model)}): {model}.')
    return model

def train(model: Sequential, x: array, y: array, epochs: int, batchSize: int) -> Tuple[Sequential, DataFrame, DataFrame]:
    '''
    train the model using "training" set of data x and y
    where x is input vectors and y is labels (vectors)
    '''
    logging.debug(f'Training model...')
    history: History = model.fit(x=x, y=y, epochs=epochs, batch_size=batchSize)
    hist: DataFrame = DataFrame(history.history)
    epochs: DataFrame = history.epoch
    logging.debug(f'History:\n\n{hist}')
    logging.debug(f'Epochs:\n\n{epochs}.')
    logging.debug(f'Model trained.')
    return model, history

def evaluate(model: Sequential, x: array, y: array, batchSize: int = 32, pause: bool = False) -> Tuple[float, float]:
    '''
    evaluate the accuracy and loss of a given trained model
    using the reserved "test" set of data x and y
    where x is input vectors and y is labels (vectors)
    '''
    logging.debug(f'Evaluating model (loss and accuracy)...')
    testLoss, testAccuracy = model.evaluate(x=x, y=y, batch_size=batchSize)
    logging.info(f'Test loss: {testLoss}.')
    logging.info(f'Test accuracy: {testAccuracy}.')

    if pause == True:
        input('Inspect evaluation/model test results ^^')

    return testLoss, testAccuracy

def predict(model: Sequential, x: array) -> array:
    '''
    run an prediction/inference using a trained model
    '''
    logging.debug(f'Running prediction on {x}...')

    # check for input shape
    if len(x.shape) == 1:
        # 1D numpy array passed, need to wrap in 2D format for inference,
        # AND unpack result back to 1D before returning
        logging.debug(f'1D input detected, processing...')
        r: array = model.predict(array([x]))[0]

    else:
        r: array = model.predict(x)

    logging.info(f'Prediction for {x} = {r}.')
    return r

def plotModelTrainingHistory(history: DataFrame, metricsList: List[str]) -> None:
    # plot the curve of model training results i.e. plot the loss of the model per epoch
    # list_of_metrics should be one of the names shown in:
    # https://www.tensorflow.org/tutorials/structured_data/imbalanced_data#define_the_model_and_metrics
    plt.figure()
    plt.xlabel("Epoch")
    plt.ylabel("Value")
    logging.info(f'History:\n\n{history.history}')

    for metric in metricsList:
        if metric in history.history:
            x: DataFrame = history.history[metric]
            plt.plot(history.epoch[1:], x[1:], label=metric) # excludes index 0 since it's the first iteration (no real learning has been done yet)

        else:
            msg: str = f'Could not generate plot. {metric} not provided in history.'
            logging.error(msg)
            raise Exception(msg)

    plt.legend()
    plt.show()
    return None
