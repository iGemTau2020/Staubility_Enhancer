In the feature generation process, we calculated ~2K features. 

However, when using this many features, we encountered multiple problems.
First, each feature we add creates another dimension and therefore increases the computation and memory required for our algorithm exponentially. 
A bigger problem is that each redundant feature also inserts noise, and increases the chance of overfitting. In overfitting, the model that we are developing begins to describe the random error in the data rather than the actual relationships between variables. 

Following our conversation with Prof. Zvika Shinar, we realized that our model probably performed overfitting. In order to reduce the risk of overfitting and computational complexity, we decided we should filter these features. We used two filters for this purpose:
Decrease the number of missing values for each instance (gene). 
Choose approximately 300 features from our entire set of 1000+ features. The feature selector we implemented filters features that do not contribute to the prediction, based on a correlation-based metric.

The model, code, and data for this analysis are presented within this directory.

**input**: The first input file is a csv, file containing all the features. It is too big to be uploaded to this respiratory. 

The other inputs are the fluorescence values that serve as "labels" for our model: flourescence_table_NATIVEpr_GFP, flourescence_table_NOP1pr_GFP. 

**code**: The main code is **PreProcFeatures.m**, that filters the features according to the CFS criteria (using the function called calculate_CFS.m).

**output**: name of the 300 chosen features. 
