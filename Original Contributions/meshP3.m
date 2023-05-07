function [mesh3] = meshP3(mesh2, p1)

%Input p1 -> number of nodes in the domain in the linear case. This is
%required to calculate the number of nodes that will be in the cubic case

p = mesh2.p;
e = mesh2.e;
t = mesh2.t;

num_nodes = size(p,2);
num_tri = size(t,2);

%Number of nodes in cubic mesh -> number of nodes in P2 + number of edges +
%number triangles (trivial)

%Remember, the edges of triangles that lie on the edge of the domain have
%nodes shared ONLY by the triangle that edge belongs to.

%Eulers formula derivation - calculate new number of nodes
num_edges = (3 * p1) - 3 - size(e,2);
new_num_nodes = num_nodes + num_edges + num_tri;

%Redefine sizes of matrices
p = [p,zeros(2,new_num_nodes - num_nodes)];
e = [e;zeros(1,size(e,2))];
t_new = [t(1:4, :); zeros(1,num_tri); t(5:end,:)];
t_new = [t_new(1:6, :); zeros(1,num_tri); t_new(7:end,:)];
t_new = [t_new(1:8, :); zeros(2,num_tri)];

t = t_new;

nodes_added = 0;
%iterate through all triangles in mesh
for i = 1:num_tri
    for k = 1:3

        %find difference in x and y values of outside nodes
        x1y1 = p(:,t(k,i));
        if k == 3
            k2 = 1;
        else
            k2 = k+1;
        end
        x2y2 = p(:,t(k2,i));
        dxy = x2y2 - x1y1;
        centre_node1 = p(:,t(k,i)) + (1/3)*dxy;
        centre_node2 = p(:,t(k,i)) + (2/3)*dxy;

        %check whether new node position has already been calculated 
        found = "False";
        temp_position = 0;
        %Check if nodes coordinates are the same
        for j = 1:num_nodes+nodes_added
            %Check to 3 d.p. to remove inaccuracies
            if round(p(:,j),3) == round(centre_node1,3)
                found = "True";
                temp_position = j;
            end
            if round(p(:,j),3) == round(centre_node2,3)
                found = "True";
            end
        end

        %Different cases:
        %1. Found = False -> this edge of the triangle has not had
        %the two nodes added to it -> still in original state of one
        % central midpoint

        %2. Found = True -> Nodes have already been added, will need
        %to rearrange matrix t to ensure anticlockwise trend still followed
        %by nodes

        %Change coordinates of old midpoint
        if found == "False"
            %Change coordinates
            p(:,t(2+(2*k),i)) = centre_node1;
            %Add new node
            nodes_added = nodes_added + 1;
            p(:,num_nodes + nodes_added) = centre_node2;
            t(3+(2*k),i) = num_nodes + nodes_added;
            
            %Add second centre node to matrix e if necessary.
            for m = 1:size(e,2)
                if round(p(:,e(4,m)),3) == round(x2y2,3) & round(p(:,e(5,m)),3) == round(x1y1,3)
                    e(6,m) = num_nodes + nodes_added;
                    e(7,m) = t(2+(2*k),i);
                elseif round(p(:,e(4,m)),3) == round(x1y1,3) & round(p(:,e(5,m)),3) == round(x2y2,3)
                    e(6,m) = num_nodes + nodes_added;
                    e(7,m) = t(2+(2*k),i);
                end
            end
        end
        %If node has already been calculated, we will need to account for
        %this in the matrix t in the triangle that shares the edge.
        if found == "True"
            t(3+(2*k),i) = t(2+(2*k),i);
            t(2+(2*k),i) = temp_position;
        end
    end
    %Find position of midpoint of triangle -> average of the coordinates of
    %the vertices
    vertice1 = p(:,t(1,i));
    vertice2 = p(:,t(2,i));
    vertice3 = p(:,t(3,i));
    midxy = (vertice1+vertice2+vertice3)/3;
    nodes_added = nodes_added + 1;
    p(:,num_nodes + nodes_added) = midxy;
    t(10,i) = num_nodes + nodes_added;
end

%Create mesh
mesh3 = struct;

mesh3.p = p;
mesh3.e = e;
mesh3.t = t;

end