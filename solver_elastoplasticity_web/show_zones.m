function show_zones(U,P1,P2)
global Koordinaten Elemente
colors_zones='color'; %'gray' or 'black_and_white'

for j=1:size(Elemente,1)
   if norm(P2(j,:))>0
      value(j)=2;
   else
      if norm(P1(j,:))>0
         value(j)=1;
      else value(j)=0;   
      end
   end
end

show_piecewise_constant(U,value,Elemente,Koordinaten);

if ismember(0,value)
   mymap=[0.3 0.3 0.3];  %elastic color
else
   mymap=[];
end
switch colors_zones
   case 'gray'
      if ismember(1,value)
         mymap=[mymap; 0.7 0.7 0.7];
      end
      if ismember(2,value)
         mymap=[mymap;  0.95 0.95 0.95];
      end
   case 'color'
      if ismember(1,value)
         mymap=[mymap; 0.6250 0.3906 0.2487];
      end
      if ismember(2,value)
         mymap=[mymap; 1 0.7812 0.4975];
      end
   end

colormap(mymap);

    function show_piecewise_constant(U,RefSpannung,Elemente,Koordinaten)
    % shows the deformed mesh and piecewise constant 
    % magnification factor for the displacement
    if norm(U)==0
       factor=1;
    else
       factor =10^(-round(log10(max(max(U)))));
    end
    factor=2;
    fprintf('Displacement magnified by factor %d \n',factor);
    fprintf(' \n',factor);

     Koordinatendeformed = Koordinaten + factor * U;
     X = reshape(Koordinatendeformed(Elemente',1),size(Elemente,2), ...
                 size(Elemente,1));
     Y = reshape(Koordinatendeformed(Elemente',2),size(Elemente,2), ...
         size(Elemente,1));
     fill(X,Y,RefSpannung)
     %axis([0 1.25 0 1])
     hold off

     shading flat %to see finite elemens or not
     %colorbar
     %drawnow
    end
end

