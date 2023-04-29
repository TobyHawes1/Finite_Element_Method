function [] = endlabfem2d (myPath)

[folder,subFolders] = feval(myPath);

for k=1:length(subFolders)
    subFolders{k} = strcat(folder,subFolders{k});
    rmpath(genpath(subFolders{k}));
end
rmpath(folder);

end