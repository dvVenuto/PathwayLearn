# PathwayLearn
Machine learning applied to RNA-seq pathway analysis \n
Authors: David Venuto, Wanling Yang

**About**
__ __
The PathwayLearn function applies L1 regularized logistic regression to RNA-seq count data sets to identify significantly deregulated KEGG pathways.  An AUC is generated for each pathway after training and testing the model and pathways are ranked in by AUC.  We propose that the AUC is roughly proportional to the level of deregulation.  Coefficients for each gene within a pathway are also reported which are proportional to there predictive power.

The pathway learn function takes arguments of:

1.	A training data frame representing the data you wish to train the algorithm on
2.	A testing or holdout data frame
3.	The column number on which your training data is separated between case and control
4.	The column number which separates your testing data between case and control
5.	Alpha value
6.	The name of the first condition (ex. Case)
7.	The name of your 2nd condition (ex. Control)
8.	A boolean value which dictates if the data will be quantile normalized.

Returned is an object with parameters:

1.	Predication accuracy rate.
2.	List of Genes in each pathway
3.	A list of predictions for each sample within each pathway.
4.	The coefficient values of each gene within each pathway.
