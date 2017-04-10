# PathwayLearn
Machine learning applied to RNA-seq pathway analysis (Publication Pending)

Full package release after publication.

Authors: David Venuto, Wanling Yang

The University Of Hong Kong

**About**
__ __
The PathwayLearn function applies L1 regularized logistic regression to RNA-seq count data sets to identify significantly deregulated KEGG pathways.  An AUC is generated for each pathway after training and testing the model and pathways are ranked in by AUC.  We propose that the AUC is roughly proportional to the level of deregulation.  Coefficients for each gene within a pathway are also reported, which are proportional to there predictive power of the gene in classifying the sample's expression signature.  Full package implementation coming pendin publication.

The PathwayLearn function takes arguments of:

1.	A training data frame representing the data you wish to train the algorithm on.
2.	The column number on which your training data is separated between case and control.
3.	Alpha value.
4.	The name of the first condition (ex. Case).
5.	The name of your 2nd condition (ex. Control).
6.	A boolean value which dictates if the data will be quantile normalized.
7.	A boolean value which dictates if the data will be limma voom normalized.
8.	The number of folds in the cross validation.

Returned is an object with parameters:

1.	Predication accuracy rate.
2.	List of Genes in each pathway.
3.	A list of predictions for each sample within each pathway.
4.  A data frame with the AUC values for each pathway.
5.	The coefficient values of each gene within each pathway.
