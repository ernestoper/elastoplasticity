function [coordinates,elements3,dirichlet,neumann]=refinement_uniform(coordinates,elements3,dirichlet,neumann) 
%function: [coordinates,elements3,dirichlet]=refinement_uniform(coordinates,elements3,dirichlet) 
%requires: getEdges, symetrizeMatrix  
%uniform refinement of a 2D triangulation 

%uniform refinement   
[element2edges, edge2nodes]=getEdges(elements3);    
nodes2edge=sparse(edge2nodes(:,1),edge2nodes(:,2),1:size(edge2nodes,1),size(coordinates,1),size(coordinates,1));
nodes2edge=symetrizeMatrix(nodes2edge); 
    
%elements on uniformly refined mesh
elements3_internal=element2edges+size(coordinates,1);
elements3_refin1= [elements3(:,1) elements3_internal(:,3) elements3_internal(:,2)];
elements3_refin2= [elements3(:,2) elements3_internal(:,1) elements3_internal(:,3)];
elements3_refin3= [elements3(:,3) elements3_internal(:,2) elements3_internal(:,1)];    
elements3=[elements3_internal; elements3_refin1; elements3_refin2; elements3_refin3];  

%dirichlet edges of uniformly refined mesh
dirichlet_edges=diag(nodes2edge(dirichlet(:,1),dirichlet(:,2)));
dirichlet=[dirichlet(:,1) dirichlet_edges+size(coordinates,1); dirichlet_edges+size(coordinates,1) dirichlet(:,2)];

%neumann edges of uniformly refined mesh
if ~isempty(neumann)
    neumann_edges=diag(nodes2edge(neumann(:,1),neumann(:,2)));
    neumann=[neumann(:,1) neumann_edges+size(coordinates,1); neumann_edges+size(coordinates,1) neumann(:,2)];
end

%coordinates of uniformly refined mesh
coordinates_internal=(coordinates(edge2nodes(:,1),:)+coordinates(edge2nodes(:,2),:))/2;
coordinates=[coordinates; coordinates_internal];   

    function [element2edges, edge2nodes]=getEdges(elements)
    %function: [element2edges, edge2nodes]=edge_numbering(elements)
    %requires: deleterepeatedrows
    %generates edges of (triangular) triangulation defined in elements
    %elements is matrix, whose rows contain numbers of its element nodes 
    %element2edges returns edges numbers of each triangular element
    %edge2nodes returns two node numbers of each edge
    %example in 2D: [element2edges, edge2nodes]=getEdges([1 2 3; 2 4 3])
    %example in 3D: [element2edges, edge2nodes]=getEdges([1 2 3 4; 1 2 3 5; 1 2 4 6])

    %2D case
    if (size(elements,2)==3)
        %extracts sets of edges 
        edges1=elements(:,[2 3]);
        edges2=elements(:,[3 1]);
        edges3=elements(:,[1 2]);

        %as sets of their nodes (vertices)
        vertices=zeros(size(elements,1)*3,2);
        vertices(1:3:end,:)=edges1;
        vertices(2:3:end,:)=edges2;
        vertices(3:3:end,:)=edges3;

        %repeated sets of nodes (joint edges) are eliminated 
        [edge2nodes,element2edges]=deleterepeatedrows(vertices);
        element2edges=reshape(element2edges,3,size(elements,1))';
    end

    %3D case
    if (size(elements,2)==4)
        %extracts sets of edges 
        edges1=elements(:,[1 2]);
        edges2=elements(:,[1 3]);
        edges3=elements(:,[1 4]);
        edges4=elements(:,[2 3]);
        edges5=elements(:,[2 4]);
        edges6=elements(:,[3 4]);

        %as sets of their nodes (vertices)
        vertices=zeros(size(elements,1)*6,2);
        vertices(1:6:end,:)=edges1;
        vertices(2:6:end,:)=edges2;
        vertices(3:6:end,:)=edges3;
        vertices(4:6:end,:)=edges4;
        vertices(5:6:end,:)=edges5;
        vertices(6:6:end,:)=edges6;

        %repeated sets of nodes (joint edges) are eliminated 
        [edge2nodes,element2edges]=deleterepeatedrows(vertices);
        element2edges=reshape(element2edges,6,size(elements,1))';
    end
    end

    function A_sym = symetrizeMatrix(A)
    [i,j,k]=find(A);
    W=sparse([i; j], [j; i], ones(size(k,1)*2,1));
    A_help=sparse([i; j], [j; i], [k; k]);
    [i,j,k]=find(A_help);
    [i,j,kk]=find(W);
    A_sym=sparse(i,j,(kk.^(-1)).*k); %Now Kantennr_sym is a symetric form of Kantennr
    end
   
    function [matrix,I]=deleterepeatedrows(matrix)
    %function: [element2edges, edge2nodes]=edge_numbering(elements)
    %requires: deleterepeatedrows
    %generates edges of (triangular) triangulation defined in elements
    %elements is matrix, whose rows contain numbers of its element nodes 
    %element2edges returns edges numbers of each triangular element
    %edge2nodes returns two node numbers of each edge
    %example: [element2edges, edge2nodes]=edge_numbering([1 2 3; 2 4 3])


    %fast and short way suggested by John D'Ericco working in both 2D and 3D
    matrixs=sort(matrix,2);
    [dummy,J,I] = unique(matrixs,'rows');
    %I=reshape(I,size(matrixs,2),size(I,1)/size(matrixs,2));
    matrix=matrix(J,:);
    end


end


