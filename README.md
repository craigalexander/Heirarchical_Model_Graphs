# Hierarchical_Model_Graphs
Code for the construction of hierarchical models for multiple responses with graphical model visualisation 

## File Description

Here is a brief description of all the code files in the folder

ChordalGraph.R - This file has code to compute the model density values for the undirected graph between responses for the chordal case
DoubleMHMultReffs.R - This is the model selection code for the undirected graph between responses (when 4 responses or greater)
GWishGen.cpp - C++ code for the G-Wishart sampler used for the graphical model selection
MixModelFunctions.R - Contains functions used within the sampler for the Bayesian hierarchical model
ModelCode.R - Demonstration code used for running the Bayesian hierarcical model on sample from the Sounds of the City corpus.
PlotGM.R - The plotting function used to create a graphical model from model output
PlotGraph.R - Wrapper function to produce multiple graphical model plots for model output.
SimplifiedSampler.R - The main sampler code used to compute the Bayesian hierarcical model.
Summary.R - The summary function used to summarise the model output.
