function [] = plotLaplace (mesh,domain,u,plotExact,t,tIndex)

fem = strcat(domain,'_v1_fem');
femPQ = mesh.(fem)(1,1);
switch femPQ
    case 'P'
       [mesh,U] = preparePlotMeshP(mesh,domain,u,'v1');
        nodes    = strcat(domain,'_v1_p');
        elements = strcat(domain,'_v1_t');
        figure;
        pdeplot(mesh.(nodes),[],mesh.(elements),...
            'XYData',U.d1_v1,'ZData',U.d1_v1,'ColorMap','default','Mesh','on');
        grid on; title('Computed solution');
        if plotExact
            U = laplaceData(30,mesh.(nodes)(1,:),mesh.(nodes)(2,:),t);
            figure;
            pdeplot(mesh.(nodes),[],mesh.(elements),...
                'XYData',U,'ZData',U,'ColorMap','default','Mesh','on');
            grid on; title('Exact solution');
        end
        %
    case 'Q'
        [X,Y,U] = preparePlotMeshQ (mesh,domain,u,'v1',[],tIndex);
        figure;
        sup = surf(X,Y,U);
        set(sup,'EdgeColor','black','FaceColor','interp');
        grid on; title('Computed solution');
        if plotExact
            U = laplaceData(30,X,Y,t);
            figure;
            sup = surf(X,Y,U);
            set(sup,'EdgeColor','black','FaceColor','interp');
            grid on; title('Exact solution');
        end
end

end