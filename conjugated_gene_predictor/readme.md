 Prediction
 -----------
 This script is the conjugated gene prediction function. This function gets a table of genes and features, and predicts the ten best candidate genes most likely to get higher fluorescent protein expression. This is predicted using a RandomForest algorithm, with saved weights file.


 Input
 ----------
 - input_features - numpy array
       normalized array with shape of MxN (M -number of genes, N number of features)

 Output
 -------
 - top_names - list of strings
               list of the 10 best genes names
 - top_vals - list of strings of floats
               list of the 10 best genes normalized scores

 Requirements
 -------
 - pickle
 - numpy
 - pandas
 - sklearn.ensemble.RandomForestRegressor
  
