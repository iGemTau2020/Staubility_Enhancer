from sklearn.ensemble import RandomForestRegressor
from sklearn import metrics
import pickle
from os import getcwd
from os.path import join, sep
import pandas as pd
import numpy as np


def prediction(input_features):
    """
           Precdiction
                    this function gets genes with his features and predict the best gene
                    based on RandomForest trained algorithm

           ...

           Input
           ----------
                input_features - numpy array
                           normilzed array with shape of NxM (N -number of genes, M number of featurs (201)

           Output
           -------
           top_names - list of strings
                    list of the 10 best genes names
            top_vals - list of strings of floats
                    list of the 10 best genes normlized scores
           """
    path_to_model = join(getcwd(), r'In\\models_and_features\\finalized_model.sav')
    names = input_features.index.values.tolist()
    test_set = np.array(input_features.values)
    loaded_model = pickle.load(open(path_to_model, 'rb'))
    y_pred = loaded_model.predict(test_set)
    y_top = np.array(y_pred).argsort()[::-1][:10]
    top_names = [names[x] for x in y_top]
    best_val = y_pred[y_top[0]]
    top_vals = [y_pred[x]/best_val for x in y_top]
    return top_names , top_vals
