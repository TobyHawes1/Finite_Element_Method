function [] = beginlabfem2d (myPath)

[folder,subFolders] = feval(myPath);

addpath(folder);
for k=1:length(subFolders)
    subFolders{k} = strcat(folder,subFolders{k});
    addpath(genpath(subFolders{k}));
end

end