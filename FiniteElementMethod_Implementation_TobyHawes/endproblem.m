function [] = endproblem (myPath,problem)

[Folder,SubFolders] = feval(myPath);

problem = strcat(Folder,'data/',problem);

rmpath(problem);

end