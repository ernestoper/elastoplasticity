function mesh_preparation
global Koordinaten Elemente Dirichlet
global STEMAelastic maske Areaglobal DirichletKnoten Rglobal

NE=size(Elemente,1); %number of triangles
NLB=6; %number of local basic functions, it must be known!
Elemente_elasticity=2*Elemente(:,[1 1 2 2 3 3])-kron(ones(NE,1),[1,0,1,0,1,0]);

Rglobal=zeros(3,6,NE);
Areaglobal=zeros(NE,1); 

Y_3Dmatrix=reshape(repmat(Elemente_elasticity,1,NLB)',NLB,NLB,NE);
X_3Dmatrix=permute(Y_3Dmatrix,[2 1 3]);
Z_3Dmatrix=0*X_3Dmatrix;
for j=1: NE;
   Elementelocal=Elemente(j,:);
   Koordinatenlocal=Koordinaten(Elementelocal,:);
   Rglobal(:,:,j)=R(Koordinatenlocal); 
   Areaglobal(j)=det([1 1 1;Koordinatenlocal']);
   Z_3Dmatrix(:,:,j)=STEMA3(Koordinatenlocal);
end
STEMAelastic=sparse(X_3Dmatrix(:),Y_3Dmatrix(:),Z_3Dmatrix(:));

%Preparation - extracting Dirichlet nodes
maske=zeros(size(Koordinaten,1),1);
maske(Dirichlet)=ones(size(Dirichlet));
DirichletKnoten=find(maske);
incorporate_Dirichlet;

    function R=R(Knoten)
    PhiGrad=[1 1 1;Knoten']\[zeros(1,2);eye(2)];
    R([1,3],[1,3,5])=PhiGrad';    
    R([3,2],[2,4,6])=PhiGrad';
    end

    function STEMA3=STEMA3(Knoten)
    %calculates \int_elemet C epsilon_i:\epsilon_j, 
    %i,j=1..6 for tracefree and symmetic Plasticstrain
    global C
    Rlocal=R(Knoten);
    STEMA3=det([1 1 1;Knoten'])/2*Rlocal'*C*Rlocal;
    end

    function incorporate_Dirichlet
    %global DirichletKnoten Koordinaten
    global W B 
    %Incomporating Dirichlet conditions
      [W,M] = u_D(Koordinaten(DirichletKnoten,:));
       B= sparse(size(W,1),2*size(Koordinaten,1));
       for j=1:size(DirichletKnoten,1) 
         B(2*j-1,2*DirichletKnoten(j)-1:2*DirichletKnoten(j)) = M(2*j-1,1:2);
         B(2*j,2*DirichletKnoten(j)-1:2*DirichletKnoten(j)) = M(2*j,1:2);
       end
       maske = find(sum(abs(B)'));
       B = B(maske,:);
       W = W(maske,:);
    end

    function [W,M] = u_D(locKoordinaten)
        global problem
        M = zeros(2*size(locKoordinaten,1),2);
        W = zeros(2*size(locKoordinaten,1),1);

        switch problem
              case 'beam1D',
                 temp = find(locKoordinaten(:,1)==0);
                 M(2*temp-1,1:2) = ones(size(temp,1),1)*[1 0];
                 temp = find(locKoordinaten(:,2)==0); 
                 M(2*temp,1:2) = ones(size(temp,1),1)*[0 1];
              case {'beam2D','cook'}
                 temp = find(locKoordinaten(:,1)==0);
                 M(2*temp-1,1:2) = ones(size(temp,1),1)*[1 0];
                 M(2*temp,1:2) = ones(size(temp,1),1)*[0 1];
              case {'ring','platehole'},
                 temp = find(locKoordinaten(:,1)==0);
                 M(2*temp-1,1:2) = ones(size(temp,1),1)*[1 0];
                 temp = find(locKoordinaten(:,2)==0); 
                 M(2*temp,1:2) = ones(size(temp,1),1)*[0 1];
              case 'Lshape',
                 temp1 = find(abs(locKoordinaten(:,1)-locKoordinaten(:,2))==2);
                 temp2 = find(abs(locKoordinaten(:,2)+locKoordinaten(:,1))==2);
                 temp3a = find(locKoordinaten(:,1)-locKoordinaten(:,2)==0);
                 temp3b=  find(locKoordinaten(:,1)<=0);
                 temp3=intersect(temp3a,temp3b);
                 temp4a = find(locKoordinaten(:,1)+locKoordinaten(:,2)==0);
                 temp4b=  find(locKoordinaten(:,1)<=0);
                 temp4=intersect(temp4a,temp4b);
                 temp=union(temp1,temp2);temp=union(temp,temp3);temp=union(temp,temp4);
                 M(2*temp-1,1:2) = ones(size(temp,1),1)*[1 0];
                 M(2*temp,1:2) = ones(size(temp,1),1)*[0 1];
                 alpha=0.544483737;
                 C1=-(cos((alpha+1)*3*pi/4))/(cos((alpha-1)*3*pi/4));
                 global lambda mu
                 C2=(2*(lambda+2*mu))/(lambda+mu);
                 for i=1:size(temp,1);
                    x=locKoordinaten(temp(i),1);
                    y=locKoordinaten(temp(i),2);
                 end
        end
    end
end
