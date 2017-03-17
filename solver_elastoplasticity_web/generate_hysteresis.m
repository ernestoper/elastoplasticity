if strcmp(problem,'cook')
   gvalue=g([48 60],[0 0],t(counter));
   hysteresis_g=[hysteresis_g gvalue(2)];
   hysteresis_u=[hysteresis_u U(find(Koordinaten(:,2)==60),2)];
end
if strcmp(problem,'beam1D')
   gvalue=g([1 0.5],[1 0],t(counter));
   hysteresis_g=[hysteresis_g gvalue(1)];
   index=intersect(find(Koordinaten(:,2)==1),find(Koordinaten(:,1)==1));
   hysteresis_u=[hysteresis_u U(index,1)];
end
if strcmp(problem,'beam2D')
   gvalue=g([1 0.5],[1 0],t(counter));
   
   hysteresis_g=[hysteresis_g gvalue(1)];
   index=intersect(find(Koordinaten(:,2)==1),find(Koordinaten(:,1)==1));
   hysteresis_u=[hysteresis_u U(index,1)];
end


