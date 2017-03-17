function [Koordinaten,Elemente,Dirichlet,Neumann]=Lshape_mesh(number_of_refinements)

load model_Lshape

%Uniforme Verfeinerung
netz=1;
parameter = 0;
while netz <= number_of_refinements
   [Koordinaten,Elemente,Dirichlet,Neumann]=refinement_uniform(Koordinaten,Elemente,Dirichlet,Neumann);
  netz = netz + 1;
end
