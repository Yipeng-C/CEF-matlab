function[O20,O22,O40,O42,O43,O44,O60,O62,O63,O64,O66,Jx,Jy,Jz,Jplus,Jminus,Jsquare,Unit] = OperatorTotalmomentum(J)
%function[O20,O40,O43,O60,O63,O66,Jx,Jy,Jz,Jplus,Jminus] = OperatorTotalmomentum(J)
% basis states are the L,S,Lz,Sz states written in such a way that
% the first 1st line : (-3,-S) 2nd line : (-3,-S+1)...(-3,S),(-2,-S),....etc

%set the size of the operators that will be use for the Hamiltonian
Jplus  = zeros((2.*J+1),(2.*J+1));
Jminus = zeros((2.*J+1),(2.*J+1));
Jz     = zeros((2.*J+1),(2.*J+1));
Unit = eye(2.*J+1);
%set down the matrix Jz
Jz = zeros((2.*J)+1,(2.*J)+1);
for i = 1:((2.*J)+1)
    Jz(i,i) = ((-J) + (i-1));
end

%set down the matrix J+ and J-
for i = 1:(2.*J)
    Jplus(i+1,i) = sqrt( ((J) - (-(J  )+(i-1))).*((J)+((-J)+(i-1))+1) );
    Jminus(i,i+1)  = sqrt( ((J) + (-(J-1)+(i-1))).*((J)-((-(J-1))+(i-1))+1) );
end

%definition of the Jx and Jy operator
imag = sqrt(-1);
Jx = (0.5).*(Jplus + Jminus);
%Jy = (0.5).*imag.*(Jminus - Jplus);
Jy = (0.5)*imag.*(Jminus - Jplus);

%set the L square matrix
Jsquare = (Jx*ctranspose(Jx)) + (Jy*ctranspose(Jy)) + (Jz*ctranspose(Jz));
Jplussquare = Jplus*Jplus;
Jminussquare = Jminus*Jminus;

%set the CEF hamiltonian
O20 = 3.*(Jz^2) - Jsquare;
O22 = 0.5.*(Jplus^2+Jminus^2);
O40 = 35.*(Jz^4) - 30.*(Jsquare*(Jz^2)) + 25.*(Jz^2) -6.*(Jsquare) + 3.*(Jsquare^2);
O42 = 0.25.*((7*(Jz^2)-Jsquare-5*Unit)*(Jplus^2+Jminus^2)+(Jplus^2+Jminus^2)*(7*(Jz^2) - Jsquare - 5*Unit));
O44 = 0.5.*(Jplus^4+Jminus^4);
O60 = 231.*(Jz^6) - 315.*(Jsquare*(Jz^4)) + 735.*(Jz^4) + 105.*((Jsquare^2)*(Jz^2)) - 525.*(Jsquare*(Jz^2)) + 294.*(Jz^2) - 5.*(Jsquare^3) + 40.*(Jsquare^2) - 60.*(Jsquare);
O43 = 0.25.*(Jz*(Jplus^3+Jminus^3) + (Jplus^3+Jminus^3)*Jz);
O62 = 0.25.*((33.*(Jz^4) - (18.*Jsquare + 123*Unit)*Jz^2 + Jsquare^2 + 10.*Jsquare + 102*Unit) * (Jplus^2 + Jminus^2) + (Jplus^2 + Jminus^2) * (33.*(Jz^4) - (18.*Jsquare + 123*Unit)*Jz^2 + Jsquare^2 + 10.*Jsquare + 102.*Unit));
O63 = 0.25.*( ((11.*Jz^3 -3.*(Jsquare*Jz) - 59.*Jz)*(Jplus^3+Jminus^3)) +  ((Jplus^3+Jminus^3))*(11.*Jz^3 -3.*(Jsquare*Jz) - 59.*Jz));
O64 = 0.25.*((11.*Jz^2 - Jsquare - 38*Unit)*(Jplus^4+Jminus^4) + (Jplus^4+Jminus^4)*(11.*Jz^2 - Jsquare - 38*Unit));
O66 = 0.5.*(Jplus^6 + Jminus^6);

end