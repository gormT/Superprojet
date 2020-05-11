myDir = uigetdir; %gets directory
myFiles = dir(fullfile(myDir,'*.msh')); %gets all msh files in struct
for k = 1:length(myFiles)
  baseFileName = myFiles(k).name;
%   fullFileName = fullfile(myFolder, baseFileName);
  fprintf(1, 'Now reading %s\n', baseFileName);
  msh2mesh(baseFileName);
end