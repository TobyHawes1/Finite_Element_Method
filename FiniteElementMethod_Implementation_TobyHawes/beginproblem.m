function [] = beginproblem (myPath,problem)

[folder,subFolders] = feval(myPath);

problem = strcat(folder,'data/',problem);

addpath(problem);

end