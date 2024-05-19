## Optimizing the EM Algorithm with Initial Clustering

This project script is designed to enhance the Expectation-Maximization (EM) algorithm by utilizing initial clustering methods to provide starting values for the EM steps. This approach aims to improve the accuracy and convergence of the EM algorithm when estimating the parameters of mixed distributions in statistical data.

# Description

The script uses various clustering techniques (K-means, hierarchical median, hierarchical clustering, DBSCAN) to segment the data into clusters. These clusters serve as initial guesses for the parameters of the EM algorithm, which is applied subsequently to refine these estimates and model the underlying distributions effectively.

# Libraries Used
janitor: For cleaning dataset names.
dplyr: For data manipulation.
MASS: For traditional statistical models.
cluster: For clustering algorithms.
dbscan: For density-based clustering.
mixtools: For mixture model analysis.
factoextra: For extracting and visualizing the results of multivariate data analyses.
ggplot2, ggpubr: For creating advanced graphics.
tidyr: For data tidying.
insight: For easy access to model information.
Data

The script processes NBA player statistics, focusing on player performances across various seasons. The data is cleaned, and NA values are removed to ensure the quality and reliability of the analysis.

# Functions

clus_init_EM
This function takes a clustering method, cluster size, and dataset as input. It applies the specified clustering algorithm to the data to form initial clusters, then uses these clusters to initialize the parameters for the EM algorithm. The function checks for sufficient sample sizes in each cluster and handles potential issues by returning NULL if a cluster is inadequately populated.

EM Algorithm Execution
The EM algorithm (mvnormalmixEM from the mixtools library) is applied to the dataset with parameters initialized from the clusters. The script calculates the mixture model parameters such as means, variances, and mixture weights.

# Plots and Summaries

Several plots are generated to visualize the distribution of data points and the effectiveness of the clustering and EM algorithms:

Histograms of player statistics.
Scatter plots comparing the actual values with estimates derived from various priming methods.
Bar plots showing the number of iterations for each priming method.