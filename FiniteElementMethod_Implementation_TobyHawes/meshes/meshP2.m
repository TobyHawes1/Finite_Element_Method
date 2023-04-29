function [mesh2] = meshP2(mesh1)

p = mesh1.p;
e = mesh1.e;
t = mesh1.t;

%Redefine the sizes of matrices e,t
e = [e;zeros(1,size(e,2))];
t = [t;zeros(3,size(t,2))];

nodes_added = 0;
num_nodes = size(p,2);

%Eulers formula Derivation -> Number of midpoints added = number of edges
num_edges = (3 * size(p,2)) - 3 - size(e,2);
new_nodes = num_edges;
%redefine size of matrix p
p = [p,zeros(2,new_nodes)];

%Change sizes of matrices

%Calculate and create nodes at midpoint of each node
for i=1:size(t,2)
    for k = 1:3

        %Get x and y values of nodes at vertices of triangles
        x1y1 = p(:,t(k,i));

        if k == 3
            k2 = 1;
        else
            k2 = k + 1;
        end
        x2y2 = p(:,t(k2,i));
        dxy = x2y2 - x1y1;
        mid_node = p(:,t(k,i)) + (1/2)*dxy;
        found = "False";
        
        %Check if node already in matrix "p" to avoid duplicate nodes.
        %Check to 3d.p. to remove inaccuracies
        for j = 1:num_nodes+nodes_added
            if round(p(:,j),3) == round(mid_node,3)
                found = "True";
                %Store what index the node is saved under, so we can label
                %the node in the matrix with this position.
                temp_position = j;
            end
        end

        %Different cases:
        %1. found = True -> midpoint has already been put into matrix p, it
        %now needs to have its index added to the matrix t

        %2. found = False -> midpoint has not yet been calculated. Add the
        %midpoint to the matrix p and its index into t

        if found == "False"
            nodes_added = nodes_added + 1;
            p(:,num_nodes + nodes_added) = mid_node;
            t(k+3,i) = num_nodes + nodes_added;

            %Check if nodes we are claculating the midpoint between lie on
            %the edge of the domain. In this case, add the midpoint index
            %to the matrix "e"
            for y = 1:size(e,2)
                if e(4,y) == t(k,i) && e(5,y) == t(k2,i)
                    e(6,y) = num_nodes + nodes_added;
                elseif e(4,y) == t(k2,i) && e(5,y) == t(k,i)
                    e(6,y) = num_nodes + nodes_added;
                end
            end
        end

        %Add the index to matrix "t"
        if found == "True"
            t(3+k,i) = temp_position;
        end
    end
end

%Create the mesh
mesh2 = struct;

mesh2.p = p;
mesh2.e = e;
mesh2.t = t;
