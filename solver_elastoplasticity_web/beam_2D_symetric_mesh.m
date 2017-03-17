function [Koordinaten,Elemente,Dirichlet,Neumann]=beam_2D_symetric_mesh(number_of_refinements)
%A square beam fixed on the left edge in both directions and pulled on the right edge (Neumann conditions)

load model_beam

%Uniforme Verfeinerung
netz=1;
parameter = 0;
while netz <= number_of_refinements
   [Koordinaten,Elemente,Dirichlet,Neumann]=refinement_uniform(Koordinaten,Elemente,Dirichlet,Neumann);
  netz = netz + 1;
end
