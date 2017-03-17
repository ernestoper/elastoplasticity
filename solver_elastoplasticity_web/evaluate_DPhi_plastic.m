function DPhi=evaluate_DPhi_plastic(U,P1prev,P2prev,plasticelements)
global STEMAelastic Elemente Rglobal Areaglobal
global C mu sigmay1 sigmay2 h1 h2 tolerance

DPhi=STEMAelastic;
for j=1:size(Elemente,1);
   if plasticelements(j)>0
      %assignment of real indices
      I=2*Elemente(j,[1 1 2 2 3 3])-[1,0,1,0,1,0];
      P1prevlocal=P1prev(j,:);
      P2prevlocal=P2prev(j,:);
      Ulocal=U(Elemente(j,:),:);
      Rlocal=Rglobal(:,:,j);
      Arealocal=Areaglobal(j);

      STEMA3Plastic=zeros(6);
      for i=1:6
         Uvectorlocal=vectorform(Ulocal);
         %setting up the increment for the approximation of DPhi
         epsilon=sqrt(eps)*max(1,abs(Uvectorlocal(i)));
         Upluslocal=Uvectorlocal; 
         Uminuslocal=Uvectorlocal;
         Upluslocal(i)=Upluslocal(i)+epsilon;
         Uminuslocal(i)=Uminuslocal(i)-epsilon;
         Upluslocal=matrix2form(Upluslocal);
         Uminuslocal=matrix2form(Uminuslocal);
   
         [P1pluslocal,P2pluslocal]=evaluate_P_on_element(Rlocal,Upluslocal,P1prevlocal,P2prevlocal,C,mu,sigmay1,sigmay2,h1,h2,tolerance);
         [P1minuslocal,P2minuslocal]=evaluate_P_on_element(Rlocal,Uminuslocal,P1prevlocal,P2prevlocal,C,mu,sigmay1,sigmay2,h1,h2,tolerance);

         %approximation of the derivation by the difference
         STEMA3Plastic(i,:)=integral_plastic(Rlocal,Arealocal,(P1pluslocal+P2pluslocal-P1minuslocal-P2minuslocal)/2/epsilon,mu);
      end
      DPhi(I,I)=DPhi(I,I)-STEMA3Plastic';
   end
end

    function integral=integral_plastic(Rlocal,Arealocal,Plasticstrain,mu)
    %calculates \int_element CPlasticstrainn:\epsilon_i, i=1..6 
    %for trace-free and symmetric Plasticstrain
    %global mu 
    Plasticvector=[Plasticstrain(1) -Plasticstrain(1) Plasticstrain(2)];
    integral=mu*Arealocal*(Plasticvector*Rlocal); %Area/2 * R'*2*mu*Plastic
    end

    function matrix2form=matrix2form(vector)
    matrix2form=reshape(vector,2,prod(size(vector))/2)';
    end

    function vectorform=vectorform(matrix)
    if size(matrix,1)>=size(matrix,2)
       matrix=matrix';
    end
    vectorform=reshape(matrix,prod(size(matrix)),1);
    vectorform=vectorform;
    end

end
