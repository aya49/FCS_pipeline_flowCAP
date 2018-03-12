# Feature-based Comparison of FCM Files: A Data Processing Pipeline

- **Input**: FlowType files derived from Flow Cytometry (FCM) FCS files
- **Output**: FCM file features, their distance matrices, classifications/clusterings based on those distance matrices, and plots & tables containing scores for each feature.

Flow Cytometry (FCM) bioinformatics is a growing sub-field of computation biology -- aimed at developing effective and efficient computational tools to store, organize, and analyze high dimensional FCM data. FCM is a high-throughput biological apparatus capable of analyzing thousands of cells per second on up to 50 features or biological markers. Hence contributing large amounts of data towards the big biological data paradigm.
	
The International Mouse Phenotyping Consortium (IMPC), on the other hand, is a collaboration between international institutions to decipher the function of 20,000 target mice genes, using mice as a model organism. IMPC is doing so by breeding mice with a certain gene knocked out -- which theoretically cancels the function of that gene. In turn, FCM is used to measure the immunological changes caused by this cancellation of gene.
	
Many tools exist to cluster or classify the cells in each FCM file into its respective cell population. However, there is a lack of tools to compare FCM files amongst each other downstream. Given that each FCM file comes from a control mouse or a mouse with a specific gene knocked out, IMPC becomes a prime motivation for this problem.
	
Therefore this repository contains a FCM data processing pipeline used to isolate features for each FCM file. It then tests the different types of features extracted on a benchmark data set, FlowCAP-II AML.


## Input Data
- The original FlowCAP-II AML Flow Cytometry FCS files can be found here: https://flowrepository.org/id/FR-FCM-ZZYA
- Each FCS file need to be first converted to a FlowType Rdata file. This can be done using the 'FlowType' function, on parameter: method='Thresholds'. Find a tutorial of the FlowType package here: https://www.bioconductor.org/packages/devel/bioc/html/flowType.html

## Pre-requisites

- Everything is written in R. All the libraries needed are listed at the top of every script -- please install as needed.
- Paths and Script options are listed at the top of every script -- please change as needed.

## Running the Code

The code is numbered based on dependency (e.g. 06 is depentant on 05, 04, etc.) and can be run in the following order:

- [00_collateFT.R](00_collateFT.R) collates the FlowType and meta-data files from a user specified repository. This script should be changed according to how the user saves the FlowType files.
-	[01_childparent_indices.R](01_childparent_indices.R) lists out the edges in the cell hierarchy representation of an FCM (in the form of an acyclic graphical ) based on the cell populations it contains.
-	[01a_normalize.R](01a_normalize.R) **normalizes** the cell counts.
-	[02_count_stats.R](02_count_stats.R) plots statistics regarding cell counts in different files.
-	[02_pvalue_single.R](02_pvalue_single.R) creates the node-based phenodeviance FCM file **features**.
-	[02_pvalue_single_half.R](02_pvalue_single_half.R) creates the phenodeviance features for if the user wishes to use half the control files as experiment files.
-	[03_childparent_matrix.R](03_childparent_matrix.R) creates the rest of the FCM file **features**.
-	[04_checkmatrix.R](04_checkmatrix.R) checks the features made for extreme values.
-	[04_dist.R](04_dist.R) creates **distance matrices** from the features.
-	[04_dist_w.R](04_dist_w.R) creates weighted manhattan **distance matrices** from the features
-	[05_dist_lin.R](05_dist_lin.R) creates distance matrices from linear combinations of existing distance matrices.
-	[06a_dist_score.R](06a_dist_score.R) creates Gaussian kernel **distance matrices**, and clusters original features using spectral clustering and DensityCut.
-	[06b_dist_score.R](06b_dist_score.R) **clusters/classifies** distance matrices. Also **scores** and plots those clusterings/classifications.
-	[06c_dist_score.R](06c_dist_score.R) creates summary plots and tables for distances & clustering/classification scores.

