main_sample.cpp creates synthetic networks. proportion of nodes, mean degree, gamma, beta, among other parameters can be set there.
main_filter_final_synt.cpp specifically for benchmark. reads the synthetic networks, computes clustering, and generates the clustering data for the rewired network. Differenciates type of nodes 1 or 2 according if it is geometric or not.
main_filter_final_real.cpp reads networks, computes clustering, and generates the clustering data for the rewired network.
filtering_synt determines the percentiles for the different synthetic samples> plots histograms and the resolution
filtering_data was used to determine the percentiles for the different datasets> plots histograms and the resolution

FOlders>
2. synthetic - fig2 has the data used later for fig2
3.DATASETS has the datasets and the clustering data used for figure 3 (see SI and references for explanation of each of them). it also has the datasets with features
3. fig3 
4. fig4 swapp sets generates fig4 for plant pollinator
5. fig5 code and data to generate fig5

features filter contains the filters used to feature selection techniques used to compare with the performance of the filter in 
