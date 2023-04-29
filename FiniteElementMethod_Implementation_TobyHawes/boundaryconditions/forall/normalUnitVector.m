function [nv] = normalUnitVector(x,y,edgeLength,domains,domain)

if (domains(1)~=0 && domains(2)==0)
    % Domain on the left
    nv = [y(2)-y(1), -(x(2)-x(1))];
elseif (domains(1)==0 && domains(2)~=0)
    % Domain on the right
    nv = [-(y(2)-y(1)), x(2)-x(1)];
elseif (domains(1)~=0 && domains(2)~=0)
    % Interface
    if domains(1)==domain
        nv = [y(2)-y(1), -(x(2)-x(1))];
    else%if dominio(1)==2
        nv = [y(2)-y(1), x(2)-x(1)];
    end
end

nv = nv./edgeLength;

end