function [P1local,P2local]=evaluate_P_on_element(Rlocal,Ulocal,P1prevlocal,P2prevlocal, C, mu, sigmay1, sigmay2, h1, h2, tolerance)

P1prevlocalmatrix=[P1prevlocal;P1prevlocal(2) -P1prevlocal(1)];
P2prevlocalmatrix=[P2prevlocal;P2prevlocal(2) -P2prevlocal(1)];
Ulocalvector([1 3 5])=Ulocal(:,1);
Ulocalvector([2 4 6])=Ulocal(:,2);

%gamma (epsilon11,epsilon22,epsilon12) on the element
gamma=Rlocal*Ulocalvector'; 
CepsU=matrixform(C*gamma);
      
%A=C epsilon(U)-(C+H)Pprev
A1=CepsU-((2*mu+h1)*P1prevlocalmatrix+2*mu*P2prevlocalmatrix);   
A2=CepsU-(2*mu*P1prevlocalmatrix+(2*mu+h2)*P2prevlocalmatrix);
devA1=dev(A1); 
devA2=dev(A2);
      
[deltaP1localmatrix,deltaP2localmatrix]=dependence(devA1,devA2, 2*mu, sigmay1, sigmay2, h1, h2, tolerance);

P1local=deltaP1localmatrix(1,:)+P1prevlocal;
P2local=deltaP2localmatrix(1,:)+P2prevlocal;

    function [P1,P2]=dependence(devA1,devA2, mu_times_2, sigmay1, sigmay2, h1, h2, tolerance)
    P2=dependence_single_general(devA2,mu_times_2,sigmay2,h2); 
    P1=dependence_single_general(devA1-mu_times_2*P2,mu_times_2,sigmay1,h1);

    normold=norm(P1,'fro')+norm(P2,'fro');

    while 1
     P2=dependence_single_general(devA2-mu_times_2*P1,mu_times_2,sigmay2,h2);           
     P1=dependence_single_general(devA1-mu_times_2*P2,mu_times_2,sigmay1,h1);       

     normnew=norm(P1,'fro')+norm(P2,'fro');
     difference=abs(normnew-normold);  
     if difference>0
        difference=difference/(normnew+normold);
     end

     if difference<tolerance 
        break
     end
     normold=normnew;
    end

        function P=dependence_single_general(devA,mu_times_2,sigmay,h)
            normdevA=norm(devA,'fro');
            if normdevA<sigmay
               P=zeros(2);
            else
               P=devA*(normdevA-sigmay)/normdevA/(mu_times_2+h);
            end

        end
    end

    function matrixform=matrixform(vector)
    matrixform=[vector(1) vector(3);
    vector(3)	vector(2)];
    end

    function A=dev(B)
    tr_over_2=(B(1,1)+B(2,2))/2;
    A=B-[tr_over_2 0; 0 tr_over_2];
    end
end


