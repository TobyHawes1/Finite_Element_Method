function [A,matrixData] = selectStokesMatrices (A,domain,matrixData,settings)

if strcmp(settings.stokesForm,'Reduced')
    dName = strcat('d',domain);
    if matrixData.stiffness==1
        matrixData.C22 = strcat(dName,'_C11');
    elseif matrixData.stiffness==2
        matrixData.C22 = strcat(dName,'_C22');
        A.(matrixData.C22) = A.(matrixData.C11);
    end
elseif strcmp(settings.stokesForm,'Complete')
    matrixData.stiffness = 2;
end

end