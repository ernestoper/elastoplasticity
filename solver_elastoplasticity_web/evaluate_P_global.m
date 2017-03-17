function [P1,P2,plasticelements]=evaluate_P_global(U,P1prev,P2prev)
global Elemente Rglobal
global C mu sigmay1 sigmay2 h1 h2 tolerance
plasticelements=zeros(size(Elemente,1),1);

P1=P1prev; P2=P2prev;
for j=1: size(Elemente,1);
      Ulocal=U(Elemente(j,:),:);
      Rlocal=Rglobal(:,:,j);
      [P1local,P2local]=evaluate_P_on_element(Rlocal,Ulocal,P1prev(j,:),P2prev(j,:),C, mu, sigmay1, sigmay2, h1, h2, tolerance);
      P1(j,:)=P1local;
      P2(j,:)=P2local;
      if (norm(P1local)+norm(P2local)==0)
         plasticelements(j)=0; %elastic element!!
      else
         if norm(P2local)>0
            plasticelements(j)=2; %multi yield element !!   
         else
            plasticelements(j)=1; %sigle yield element!!
         end
      end   
end
   

       