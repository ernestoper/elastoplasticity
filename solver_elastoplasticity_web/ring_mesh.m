function [Koordinaten,Elemente,Dirichlet,Neumann]=ring_mesh(number_of_refinements)
%A square beam fixed on the left edge in both directions and pulled on the right edge (Neumann conditions)

load model_ring

%Uniforme Verfeinerung
netz=1;
parameter = 0;
while netz <= number_of_refinements
   %[Koordinaten,Elemente,Dirichlet,Neumann]=Rotverfeinerung(Koordinaten,Elemente,Dirichlet,Neumann);
   [Koordinaten,Elemente,Dirichlet,Neumann]=refinement_uniform(Koordinaten,Elemente,Dirichlet,Neumann);
   Koordinaten=correct_ring_coordinates(Koordinaten,Neumann);
	netz = netz + 1;
end

    function Koordinaten_new=correct_ring_coordinates(~,~)
    index_Neumann=unique(Neumann);

    r=sqrt(Koordinaten(index_Neumann,1).^2+Koordinaten(index_Neumann,2).^2);
    r_ring=round(r);

    index_change=index_Neumann(find(r-r_ring));

    Koordinaten_new=Koordinaten;

    if ~isempty(index_change)
       r_to_change=r_ring(find(r-r_ring));
       Koordinaten_to_change=Koordinaten(index_change,:);
       for i=1:size(index_change,1)
          Koordinaten_changed(i,:)=Koordinaten_to_change(i,:)/norm(Koordinaten_to_change(i,:));
          Koordinaten_changed(i,:)=Koordinaten_changed(i,:)*r_to_change(i);
       end
       Koordinaten_new(index_change,:)=Koordinaten_changed;
    end
end
end
