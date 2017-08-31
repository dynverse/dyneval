function [dataFile processDataMat processDataTxt PCAdataFile dataFolder resultsDir intermediate_filesDir figuresDir] = initialization(dataset);

% set path of data files
dataFolder = dataset; %fullfile(pwd,'sample_data',dataset);
% save('dataFolder.mat','dataFolder');

% Alternatively, uncomment the following line and select interactively
% [dataFile dataFolder] = setDataFile;

% Set the results folder variable
resultsDir = fullfile(dataFolder,'results');
intermediate_filesDir = fullfile(dataFolder,'intermediate_files');
figuresDir = fullfile(dataFolder,'figures');

% Create results folder in data folder if not previously present
createSubfolders(dataFolder)
createSubfolders(resultsDir)
createSubfolders(intermediate_filesDir)
createSubfolders(figuresDir)

dataFile = fullfile(dataFolder, '/Data.txt');
processDataMat = fullfile(intermediate_filesDir, '/PData.mat');
processDataTxt = fullfile(intermediate_filesDir, '/PData.txt');
PCAdataFile = fullfile(intermediate_filesDir, '/PDataPCA.mat');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function createSubfolders(dataFolder)

if ~exist(fullfile(dataFolder),'dir')
    mkdir(dataFolder)
end

end
