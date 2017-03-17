function evaluate_fv(time)
global Koordinaten Elemente
global fv
%Calculating fv
   fv=zeros(1,2*size(Koordinaten,1));
   for j=1: size(Elemente,1);
      I=2*Elemente(j,[1 1 2 2 3 3])-[1,0,1,0,1,0]; %assignment of real indices
      fs=f(sum(Koordinaten(Elemente(j,:),:))/3,time); % f value at the middle point
      fv(I)=fv(I)+(det([1 1 1; Koordinaten(Elemente(j,:),:)'])*[fs fs fs]'/6)';   
  end
  fv=sparse(fv);   

  function f=f(x,t);
    f=[0 0];
    %f=1000*[norm(x)-1.5 norm(x)-1.5]*t;
    %if t==1
    %   f=[120 120];
    %else
    %  f=[140 140];
    %end
  end




end


  
