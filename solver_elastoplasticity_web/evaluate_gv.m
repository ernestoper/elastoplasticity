function evaluate_gv(t)
global Koordinaten Neumann
global gv
%Calculating gv
  gv=zeros(1,2*size(Koordinaten,1));
  if ~isempty(Neumann)
         n=(Koordinaten(Neumann(:,2),:)-Koordinaten(Neumann(:,1),:))*[0 -1; 1 0];
         for j=1:size(Neumann,1)
            I=2*Neumann(j,[1 1 2 2])-[1 0 1 0];
            gm=g(sum(Koordinaten(Neumann(j,:),:))/2,n(j,:)/norm(n(j,:)),t);
            gv(I)=gv(I)+norm(n(j,:))*[gm gm]/2;
         end
  end
  gv=sparse(gv);
end
