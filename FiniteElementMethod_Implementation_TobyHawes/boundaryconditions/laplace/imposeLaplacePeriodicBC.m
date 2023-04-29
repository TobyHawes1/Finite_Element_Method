function [A,rhs,dof] = imposeLaplacePeriodicBC (A,matrixData,rhs,dof,domain,variable)

domainVariable  = strcat(domain,'_',variable);
periodicBCnodes = strcat(domainVariable,'_periodicBC');
internalNodes   = strcat(domainVariable,'_1_internal');
nRows = size(A.(matrixData.C),1);

nPBC = size(dof.(periodicBCnodes),2);
% Process matrix
for i = 1:nRows
    for j = 2:nPBC-1
        index1 = dof.(periodicBCnodes)(1,j);
        index2 = dof.(periodicBCnodes)(2,j);
        A.(matrixData.C)(i,index1) = A.(matrixData.C)(i,index1) + A.(matrixData.C)(i,index2);
        A.(matrixData.C)(index1,i) = A.(matrixData.C)(index1,i) + A.(matrixData.C)(index2,i);
    end
end
index1 = dof.(periodicBCnodes)(1,1);
index2 = dof.(periodicBCnodes)(2,1);
A.(matrixData.C)(index1,index1) = A.(matrixData.C)(index1,index1) + A.(matrixData.C)(index2,index2);
index1 = dof.(periodicBCnodes)(1,nPBC);
index2 = dof.(periodicBCnodes)(2,nPBC);
A.(matrixData.C)(index1,index1) = A.(matrixData.C)(index1,index1) + A.(matrixData.C)(index2,index2);
% Process rhs
for j = 1:nPBC
    index1 = dof.(periodicBCnodes)(1,j);
    index2 = dof.(periodicBCnodes)(2,j);
    rhs.(domainVariable)(index1,1) = rhs.(domainVariable)(index1,1) + rhs.(domainVariable)(index2,1);
end
% Process internal nodes
dof.(internalNodes) = setdiff(dof.(internalNodes),dof.(periodicBCnodes)(2,:));

end