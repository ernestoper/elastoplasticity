function [Phi,plasticelements,P1,P2]=evaluate_Phi(U,P1prev,P2prev)
global fv gv STEMAelastic Rglobal Areaglobal Elemente
%global yield_type
global mu

%calculating Phi - right side for the Newton matrix - elementwise
%plasticelements=[];
Phi=zeros(1,2*size(U,1));
%Calculation C\epsilon(U):\epsilon(V), V=V_i ==> vector
Phi=(vectorform(U)'*STEMAelastic);
%-int_\Omega fv dx
Phi=Phi-fv-gv; 

[P1,P2,plasticelements]=evaluate_P_global(U,P1prev,P2prev);

for j=1: size(Elemente,1);
      %assignment of real indices
      I=2*Elemente(j,[1 1 2 2 3 3])-[1,0,1,0,1,0];       
      
      P1local=P1(j,:);
      P2local=P2(j,:);
      if plasticelements(j)~=0
         Rlocal=Rglobal(:,:,j);
         Arealocal=Areaglobal(j);
         %Calculating CP0:\epsilon(V), V=V_i ==> vector
         STEMA3Plasticlocal=integral_plastic(Rlocal,Arealocal,P1local+P2local,mu);
         Phi(I)=Phi(I)-STEMA3Plasticlocal;   
      end
end
Phi=Phi';

    function integral=integral_plastic(Rlocal,Arealocal,Plasticstrain,mu)
    %calculates \int_element CPlasticstrainn:\epsilon_i, i=1..6 
    %for trace-free and symmetric Plasticstrain
    %global mu 
    Plasticvector=[Plasticstrain(1) -Plasticstrain(1) Plasticstrain(2)];
    integral=mu*Arealocal*(Plasticvector*Rlocal); %Area/2 * R'*2*mu*Plastic
    end

    function vectorform=vectorform(matrix)
    if size(matrix,1)>=size(matrix,2)
       matrix=matrix';
    end
    vectorform=reshape(matrix,prod(size(matrix)),1);
    vectorform=vectorform;
    end

end


