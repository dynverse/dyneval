function SCUBA(dataFolder, cluster_mode)
%Main function.
%dataset      -- name of the dataset
%cluster_mode -- metric used for clustering data. 'original': Euclidean
%distance of original data; 'pca': Euclidean distance for PCA transformed
%data; 'pca2': Euclidean distance for PCA transformed data, with PCA
%applied the samples selected from the last cell-stage.

if 0,
    HD = helpdlg(sprintf(['SCUBA will now analyze the dataset']));
    waitfor(HD)
    pause(0.1)
end

if ~exist('cluster_mode'),
    cluster_mode = 'pca2';
end

% initialization of file names
[dataFile processDataMat processDataTxt PCAdataFile dataFolder resultsDir intermediate_filesDir figuresDir] = initialization(dataFolder);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% additional preprocessing %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pca_analysis(dataFolder);

%%%%%%%%%%%%%%%%%%%%%%%%%
% 1) Dynamic Clustering %
%%%%%%%%%%%%%%%%%%%%%%%%%

% Initial estimation of the lineage tree.
initialize_tree(dataFolder, cluster_mode);

% Refinement of tree structure by maximizing penalized likelihood function.
refine_tree(dataFolder, cluster_mode);


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2) Bifurcation Analysis %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Project data along bifurcation directions.
bifurcation_direction(dataFolder, cluster_mode);

% Model the dynamic change of gene expression pattern along the bifurcation
% direction by fitting the Fokker-Planck equation.
bifurcation_analysis(dataFolder)

end




