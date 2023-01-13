%function [BCoeff,Echeck,Vcheck]=YbBCoeff()
clear all
clc
% for Yb2Ti2O7
a=10.091e-10;
oneover4piepsilon0=9e9; % (1/(4*pi*epsilon) unit is N*m^2*C^-2)
q=-1.602176e-19;
J=7/2;
%hbar=6.582119514e-15;
%from N*m to eV is 1N*m= 6.24150913e18eV
unittrans=6.24150913e18;
bohr_radius=5.29177e-11;

Atom_Yb_position= a*[0.25,0.75,0.5];

%since Er is the magnetic ion set as origin need to define a (x,y,z) close
%to this magnetic ion and let r=sqrt(x^2+y^2+z^2) << R(i) as defined
%below, lets random set to like delta x&y*z= 0.01
pointxyz_position_r= a*[0.255,0.745,0.495];  %randomly chosen as pertubation.
BondYbtoxyz=pointxyz_position_r-Atom_Yb_position(1,:);
Atom_Erxyz_r=sqrt(sum(BondYbtoxyz.^2));



Atom_O_position = a*[0.375,0.875,0.375; 0.375,0.875,0.6694;
    0.125,0.91940,0.625; 0.08060,0.875,0.375;
    0.375,0.5806,0.375; 0.125,0.625,0.33060;
    0.41940,0.625,0.625; 0.125,0.625,0.625];
%Zvector is the high symmetry direction, local [111].
Zvector=Atom_O_position(1,:) - Atom_O_position(8,:);
%Zunit=Zvector/sqrt(sum(Zvector.^2)); %(which is atcually -1-11) then the defined symmetry Y-aix y= 010
Zunit=[1,1,-1]/sqrt(3);
Yunit=[1,1,2]/sqrt(1);



%Atom_Er_r = sqrt(sum(Atom_Er_position(1,:).^2));


%[phi,theta,Atom_O_R]=cart2sph(Atom_O_position(:,1),Atom_O_position(:,2),Atom_O_position(:,3));
for i =1:8
    Atom_YbObond(i,:) = Atom_Yb_position(1,:)-Atom_O_position(i,:);
    Atom_YbO_R(i) = sqrt(sum(Atom_YbObond(i,:).^2));
    theta(i) = atan2(norm(cross(Atom_YbObond(i,:),Zunit)),dot(Atom_YbObond(i,:),Zunit));
    phi(i) = atan2(norm(cross(Atom_YbObond(i,:)* sin(theta(i)),Yunit)),...
        dot(Atom_YbObond(i,:)* sin(theta(i)),Yunit));
end
A_CF_para=[];

% A_CF_para= [A20,A22,A40,A42,A44,A60,A62,A64,A66]
%calculate the theta and phi
%  for i=1:8
%      Atom_ErO_R(i)=sqrt(sum(Atom_O_position(i,:).^2));
%      theta1(i)= atan2(norm(cross(Zunit,Atom_O_position(i,:))), ...
%      dot(Atom_O_position(i,:),Zunit));
%      phi1(i) = atan2(norm( cross(Yunit,Atom_O_position(i,:) * sin(theta(i)))), ...
%         dot(Atom_O_position(i,:) * sin(theta(i)),Yunit));
% end


for j=1 %B20
    for i = 1:8
        temp1(i)= q*(4*pi/(2*2+1)) *2*q/(Atom_YbO_R(i)^(3)) *ylm(2,0,theta(i),phi(i));
    end
    A_CF_para(j)=sum(temp1);
end


for j=2  %B40
    for i = 1:8
        temp2(i)=q*(4*pi/(2*4+1)) *2*q/((Atom_YbO_R(i))^(5)) * ylm(4,0,theta(i),phi(i));
    end
    A_CF_para(j)=sum(temp2);
end

for j=3 %B43
    for i = 1:8
        temp3(i)=-q*(4*pi/(2*4+1)) *2*q/((Atom_YbO_R(i))^(5)) *1/sqrt(2) * (ylm(4,-j,theta(i),phi(i)) + (-1)^(j) * ylm(4,j,theta(i),phi(i)));
    end
    A_CF_para(j)=sum(temp3);
end


for j=4 %B60
    for i =1:8
        temp4(i)=q*(4*pi/(2*6+1))*2*q/((Atom_YbO_R(i))^(7))  * ylm(6,(j-4),theta(i),phi(i));
    end
    A_CF_para(j)=sum(temp4);
end

for j=5 %B63
    for i =1:8
        temp5(i)=-q*(4*pi/(2*6+1)) *2*q/((Atom_YbO_R(i))^(7)) *(1/sqrt(2)) * (ylm(6,-(j-2),theta(i),phi(i))+ (-1)^(j-2)*ylm(6,j-2,theta(i),phi(i)));
    end
    A_CF_para(j)=sum(temp5);
end

for j=6 %B66
    for i =1:8
        temp6(i)=q*(4*pi/(2*6+1)) *2*q/((Atom_YbO_R(i))^(7)) *1/sqrt(2) * (ylm(6,-j,theta(i),phi(i))+ (-1)^(j)*ylm(6,j,theta(i),phi(i)));
    end
    A_CF_para(j)=sum(temp6);
end

alphaEr=4/(9*25*7);
betaEr = 2/(9*5*7*11*13);
gammaEr = 8/(27*7*11^2*13^2);
alphaYb = 2/(9*7);
betaYb = -2/(3*5*7*11);
gammaYb = 4/(27*7*11*13);
%sigma=[0.460,0.0190,-0.0283]; % for Er
sigma= [0.450,0.0188,-0.0277]; % for Yb
%meanr=[0.773,1.677,8.431]; % for Er
meanr=[0.710,1.448,7.003]; % for Yb
%meanr=[0.6382,1.036,3.627]; % for Yb from Freeman Jonathans; cite.

for i=1;
    %B(i)=alphaEr*(1-sigma(1))*A_CF_para(i)*meanr(1);
    B(i)=alphaYb*A_CF_para(i)*meanr(1)*bohr_radius^2;%*(1-sigma(1));
    
end
for i=2:3
    % B(i)=betaEr*(1-sigma(2))*A_CF_para(i)*meanr(2);
    B(i)=betaYb*A_CF_para(i)*meanr(2)*bohr_radius^4;%*(1-sigma(2));
    
end
for i=4:6
    % B(i)=gammaEr*(1-sigma(3))*A_CF_para(i)*meanr(3);
    B(i)=gammaYb*A_CF_para(i)*meanr(3)*bohr_radius^6;%*(1-sigma(3));
    
end

% for garnet
% Coeff_Zna=[0.25*sqrt(5/pi),0.25*sqrt(15/pi),3/(sqrt(pi)*16),sqrt(5/pi)*3/8,sqrt(35/pi)*3/16,...
 %   sqrt(13/pi)/32,sqrt(2730/pi)/64,sqrt(13/7*pi)*21/32,sqrt(26/231*pi)*231/64];
 
 %For Yb2Ti2O7
Coeff_Zna=[0.25*sqrt(5/pi),3/(16*sqrt(pi)),sqrt(70/pi)*3/8,sqrt(13/pi)/32,sqrt(2730/pi)/32,sqrt(26/(231*pi))*231/64];

BCoeff=1000*B.*Coeff_Zna*unittrans*oneover4piepsilon0 % 1000 for meV 

[O20,O22,O40,O42,O43,O44,O60,O62,O63,O64,O66,Jx,Jy,Jz,Jplus,Jminus,Jsquare,Unit] = OperatorTotalmomentum(J);
Hcefcheck = BCoeff(1)*O20 + BCoeff(2)*O40 + BCoeff(3)*O43 + BCoeff(4)*O60 +BCoeff(5)*O63 +BCoeff(6)*O66;
Hcefcheck=round(Hcefcheck,5); 
[Vcheck,Echeck]=eig(Hcefcheck,'Vector');
 [Echeck,index]=sort(Echeck);
 Vcheck=Vcheck(index,:);
 Echeck = Echeck + abs(min(Echeck(:,1)));


 %end
