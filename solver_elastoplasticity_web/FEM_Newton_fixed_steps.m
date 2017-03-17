function [U,P1,P2]=FEM_Newton_fixed_steps(P1prev,P2prev,U,number_of_steps)
global STEMAelastic B W 
[dummy,feste2DKnoten]=find(B);

%fprintf('Evaluating Phi, ');
[Phi,plasticelements,P1,P2]=evaluate_Phi(U,P1prev,P2prev);
          
residualvector=Phi; residualvector(feste2DKnoten)=[];
residual=sqrt(residualvector'*(residualvector));
fprintf('residual = %f ', residual);
if (norm(P1)+norm(P2)==0)
         fprintf('(elasticity) \n'); 
   else
      if (norm(P2)>0)
         fprintf('(two-yield plasticity) \n');   
      else fprintf('(single-yield plasticity) \n'); 
      end      
end

for i=1:number_of_steps
   if isempty(find(plasticelements, 1))	%elasticity only
      %fprintf('Substituting elastic DPhi, ');
      %global matrix - set up
      stabilization=max(max(STEMAelastic));
      A=[STEMAelastic, stabilization*B';stabilization*B,...
                   sparse(size(B,1),size(B,1))];
   else %plasticity!!
      %fprintf('Evaluating plastic DPHi, ');
      DPhi=evaluate_DPhi_plastic(U,P1prev,P2prev,plasticelements);
      %global matrix - set up
      stabilization=max(max(DPhi));
      %stabilization=1;
      A=[DPhi,stabilization*B';stabilization*B,sparse(size(B,1),...
                   size(B,1))];
   end 
   %condestA=[condestA condest(A)];
   %N=[N size(STEMAelastic,2)];
   b=[Phi; stabilization*W];
   % one Newton iteration
   solution=A\b;
   lambda=solution(size(STEMAelastic,1)+1:end);
   Udeltavector=solution(1:size(STEMAelastic,1));
   Udelta=matrix2form(Udeltavector);
   U=U-Udelta;
   %[Phi,plasticelements,P1,P2]=evaluate_Phi2(U,P1prev,P2prev);
   [Phi,plasticelements,P1,P2]=evaluate_Phi(U,P1prev,P2prev);
   
   
   %residualvector=Phi;residualvector(feste2DKnoten)=[];
   %residual=sqrt(residualvector'*(residualvector));
   residual=norm(stabilization*B'*lambda-Phi);
   fprintf('residual = %f ', residual);
   
   %which state - elastic, single-yield or two-yield?
   if (norm(P1)+norm(P2)==0)
         fprintf('(elasticity) \n'); 
   else
      if (norm(P2)>0)
         fprintf('(two-yield plasticity) \n');   
      else fprintf('(single-yield plasticity) \n'); 
      end
   end
end

    function matrix2form=matrix2form(vector)
    matrix2form=reshape(vector,2,prod(size(vector))/2)';
    end


end
