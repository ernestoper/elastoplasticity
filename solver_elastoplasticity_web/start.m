%Date 14.4. 2010, Reykjavik
%Author: Jan Valdman, email: Jan.Valdman@gmail.com
%details are provided at http://sites.google.com/site/janvaldman/software

clear all; close all;
declaration_of_variables % specifies global variables

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
problem='beam1D'; %'beam2D' or 'beam1D' or 'cook' or 'ring' or 'platehole' or 'Lshape' ;
Newtonsteps=4; %number of Newtonsteps 
tolerance=1e-6; %for the plastic-dependence scheme
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

problem_properties; %Initial mesh, material and time properties
mesh_preparation; %gets mesh ready - stifness matrix, dirichlet nodes etc. 

%initial conditions
Uprev=zeros(size(Koordinaten,1),2); %defined on every node
P1prev=zeros(size(Elemente,1),2); %defined on every element - constant
P2prev=zeros(size(Elemente,1),2); %defined on every element - constant

%set up
hysteresis_u=[]; hysteresis_g=[];
counter=1; numberoftimes=size(t,2);

further=1;
%loop over all discrete times
while further
   %improving the space error - mesh-refining
   [U,P1,P2]=solve_one_time_step(P1prev,P2prev,Uprev,t(counter),Newtonsteps);
   Uprev=U; P1prev=P1; P2prev=P2;
   
   %generating figures
   figure(1);
   subplot(1,2,1)
   show_zones(U,P1,P2);
   axis([0 1.25 -0.1 1.1])
   title('elastoplastic zones')
   
   %output for hysteresis behavior
   generate_hysteresis;
   subplot(1,2,2)
   plot(hysteresis_u,hysteresis_g,'x-');
   axis([-0.05 0.05 -15 15])
   title('hysteresis: displacement versus surface force')
   
   figure(1);
   M(counter)=getframe(gcf);
   
   
   if (counter==numberoftimes)
      further=0;
   else 
      numberoftimes=size(t,2);
      counter=counter+1;
   end
end

%figure(2);
%show_sigma(U,P1,P2);
%title('stress')

%animation
figure(2)
axis off;
movie(M,1,3);
%movie2avi(M,'cook_video.mp4','FPS',2,'compression', 'Cinepak');
movie2avi(M,'cook_video.avi','compression', 'None');


