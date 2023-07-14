# Source Code
The main part of the code is written in Julia or R. Splitting the data, training the NN, and labeling the second data set is done in Julia for the first two experiments. Everything else, the creation and testing of the RBM and the analysis of misclassified objects was carried out in R.
"scDeepSort" is a python package and therefore called in python, the further analysis of this use case is also done in R. The necessary packages are imported at the beginning of each file, but not installed. scDeepSort was installed according to its installation guide (https://github.com/ZJUFanLab/scDeepSort/wiki). <br>

## Julia
DataPreparation.jl: Prepares data for use in NN (Amino Acid Encoding, Cross-Validation) <br>
DataSplit.jl: Split the data into the two training and the test set <br>
*EvaluateNetwork.jl: Script to train and evaluate the NNs <br>
Evaluation.jl: Compute quality measures for model <br>
ModelDefinitions.jl: Defines architecture of NNs <br>
*NewLabels.jl: Assigns labels by NN to data objects <br>
Parameters.jl: Defines parameters for NN <br>
PrepareTbData.jl: Concatenate the sequences of the used proteins <br>
TrainNetwork.jl: Functions to train the NN <br>
( *: Main scripts that call the other scripts)

## R
AlignAASeq2MSA.R: Align sequence from PDB to MSA for HPAIV analysis <br>
CheckHPAIVSerotypes.R: Check serotypes of strains to identify H5 and H7 strains <br>
*DeepSortAnalysis.R: Analysis of the results from scDeepSort (Experiment 3) <br>
*HPAIV_RosettaOnNNLabels.R: Train RBM for Experiment 1 and analyse misclassifications <br>
MCFS_Tuberculosis.R: Run Monte Carlo Feature Selection for Tuberculosis data <br>
*MTB_RosettaOnNNLabels.R: Train RBM for Experiment 2 and analyse misclassifications <br>
ProteinStructureAnalysis.R: Compare found rules of RBM with secondary structure <br>
RecalculatePositions.R: Compute original sequence positions of MSA positions <br>
( *: Main analysis of the three use cases)
