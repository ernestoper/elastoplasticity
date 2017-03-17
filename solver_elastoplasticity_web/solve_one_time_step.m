function [U,P1,P2]=solve_one_time_step(P1prev,P2prev,U,time,Newtonsteps)
global Koordinaten Elemente Dirichlet %mesh invariants
global problem

fprintf('Number of elements: %d \n',size(Elemente,1));

if strcmp(problem,'Lshape')
   for j=1:size(Dirichlet,1)
      U(Dirichlet(j,1),:)=Lshape_Dirichlet(Koordinaten(Dirichlet(j,1),:));
      U(Dirichlet(j,2),:)=Lshape_Dirichlet(Koordinaten(Dirichlet(j,2),:));
   end
end

evaluate_fv(time);
evaluate_gv(time);
[U,P1,P2]=FEM_Newton_fixed_steps(P1prev,P2prev,U,Newtonsteps); 

    function uu=Lshape_Dirichlet(x)
    global lambda mu
    [t r]=cart2pol(x(1),x(2));
    alph = .544483737;
    omega = pi * 3 / 4;
    C_1 = 1;
    C_2 = -C_1*cos((alph+1)*omega)/cos((alph-1)*omega);
    C_3 = 2*(lambda+2*mu)/(lambda+mu);
    ut = (1/(2*mu)) * r.^alph.*((alph+1)*C_1*sin((alph+1)*t)+(C_3+alph-1)*...
                  C_2*sin((alph-1)*t));
    ur = (1/(2*mu))*r.^alph .*(-(alph+1)*C_1*cos((alph+1)*t)+(C_3-(alph+1))*...
                  C_2*cos((alph-1)*t));
    uu(:,1) = ur .* cos(t) - ut .* sin(t);
    uu(:,2) = ur .* sin(t) + ut .* cos(t);
    end

end
