function [error] = laplaceError (u,mesh,domain,variable,data,time,varargin)

domainVariable = strcat(domain,'_',variable);
errL2     = strcat(domainVariable,'_L2');
errH1     = strcat(domainVariable,'_H1');

if ~isempty(varargin)
    error = varargin{1};
end

error.(errL2) = [];
error.(errH1) = [];
[error.(errL2),error.(errH1)] = error2D(u.(domainVariable),data,1,mesh,domain,variable,30);
%[error.(errL2)] = error2D(u.(domainVariable),data,1,mesh,domain,variable,30);

end