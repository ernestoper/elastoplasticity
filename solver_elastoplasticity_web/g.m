function g=g(x,normal,t)
global problem
switch problem
   case {'beam1D'}
         if (x(1)==1)
            g=12*[1 0]*sin(t*(pi/20));
         else
            g=[0 0];
         end
   case {'beam2D'}
         if (x(1)==1)&(normal==[1 0])
            g=12*[1 0]*sin(t*(pi/20)); %first experiment - hysteresis
            %g=1.2*[1 0]*t; %second experiment
            
         else
            g=[0 0];
         end      
   case 'ring'
         if norm(x,'fro')<1.5
            %g=-normal*t;
            g=t*x/norm(x);
         else
            %g=-(1/4)*normal*t;
            g=-0.25*t*x/norm(x);
         end
   case 'cook'
      if (x(1)==48)&(norm(normal-[1 0])<eps)
         g=12*[0 1]*sin(t*(pi/20));
      else
          g=[0 0];
      end
    case 'platehole'
       if (normal==[0 1])&(x(2)==100)
          g=[0 600]*t;
       else
          g=[0 0];
       end       
end
  end