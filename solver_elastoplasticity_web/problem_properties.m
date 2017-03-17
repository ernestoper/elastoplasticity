switch problem
case 'beam1D',
   mu=1000;lambda=1000; sigmay1=5; h1=100; sigmay2=7; h2=50;
   t=0:1:10*5; %hysteresis
   %t=8.5;       %one time step   
   [Koordinaten,Elemente,Dirichlet,Neumann]=beam_2D_symetric_mesh(3); 
case 'beam2D',
   mu=1000;lambda=1000; sigmay1=5; h1=100; sigmay2=7; h2=50;  
   t=0:1:10*5; %hysteresis
   %t=7;        %one time step
   [Koordinaten,Elemente,Dirichlet,Neumann]=beam_2D_symetric_mesh(0); 
case 'cook'
   mu=1000; lambda=1000; sigmay1=5; h1=100; sigmay2=7; h2=50; % - hysteresis
   sigmay2=6;  % - one time step 
   t=0:0.5:10*5; %hysteresis
   t=1.7; %one time step
   [Koordinaten,Elemente,Dirichlet,Neumann]=cook_mesh(5);      
case 'ring',
   t=0:10:430; %evolution
   t=250;	% one time step
   E=70000; nu=0.33; lambda=E*nu/((1+nu)*(1-2*nu)); mu=E/(2+2*nu); %glass
   sigmay1=243*sqrt(2/3); h1=1; sigmay2=250*sqrt(2/3); h2=1;  %only for comparison with JA
   [Koordinaten,Elemente,Dirichlet,Neumann]=ring_mesh(6); 
case 'platehole'
   t=0.4:0.1:0.8; 
   %t=[0.7]; %one time step 
   E=206900; nu=0.29; lambda=E*nu/((1+nu)*(1-2*nu)); mu=E/(2+2*nu); %glass
   mu=1000; lambda=1000;
   sigmay1=450*sqrt(2/3); h1=1; sigmay2=500*sqrt(2/3); h2=1;
   [Koordinaten,Elemente,Dirichlet,Neumann]=platehole_mesh(8);     
case 'Lshape'
   t=0.5; %one time step 
   E=100000; nu=0.3; lambda=E*nu/((1+nu)*(1-2*nu)); mu=E/(2+2*nu); %glass
   sigmay1=1; h1=2; sigmay2=1.41; h2=0.02; 
   %h1=1;
   [Koordinaten,Elemente,Dirichlet,Neumann]=Lshape_mesh(4); 
end  

mu_times_2=mu*2; C=mu*[2 0 0;0 2 0;0 0 1] + lambda*[1 1 0;1 1 0;0 0 0];

