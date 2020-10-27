 Precdiction
 -----------
      This function gets a table of genes and features and predict the best 10 genes that are most likely genes to get the higher fluorescent protein expression
      based on RandomForest algorithm, with saved weights file


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

 Requirements
 -------
     - pickle
     - numpy
     - pandas
     - sklearn.ensemble.RandomForestRegressor
  
