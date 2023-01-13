function BCoeff=Point_charge_cal_NN()

clear all;
%point charge calculation
J=15/2;
%kb=0.0862; from K to meV

%  location of the Er and Garnet. Anst
%  a is lattice parameter
a= 12.352e-10;
%a = 12.352;
q=-1.602176e-19;
%q=1;
% hbar=1.0545718e-34; ( unit is J*s)
%hbar = 4.135667662e-15; % unit is eV*s/Rad
oneover4piepsilon0=9e9; % (1/(4*pi*epsilon) unit is N*m^2*C^-2)
epsilon0=8.854187817e-12; % with unit C^2 * N-1 * m^-2
%from N*m to eV is 1N*m= 6.24150913e18eV
unittrans=6.24150913e18;
% bohr magneton
bohr_magneton=5.7883818012e-5;   % with unit eV/T
g_J=6/5;
bohr_radius=5.29177e-11; %   this is the unit for the meanr^2 and meanr^4!!!!!!
t=pi;
Ry = [cos(t) 0 sin(t); 0 1 0; -sin(t) 0 cos(t)];

%to generate the Zvector, use three Er atom in the triangluar lattice.
Atom_Er_position = a*[0.25, 0.125, 0]; 
%Atom_Er_position = a*[0, 0.25, 0.125]; 

%O position for 0, 0.25, 0.125 Er with 001 ising direction
% Atom_O_position = a*[0.05670  0.15110 -0.02660; -0.05670  0.34890 -0.02660;
%     0.19330  0.27660  0.09890; 0.02660  0.44330  0.15110;
%     0.09890  0.19330  0.27660; -0.09890  0.30670  0.27660;
%     -0.02660  0.05670  0.15110; -0.19330  0.22340  0.09890];  % 8 Nearest and then 8 NN


% O position for 0.25, 0.125, 0 Er with 010 ising direction
Atom_O_position = a*[0.3492,-0.02636,-0.05657; 0.22364,0.0992,-0.19343;
    0.44343,0.15080,0.02636; 0.30657,0.27636,-0.0992;
    0.15080,-0.02636,0.05657; 0.27636,0.09920,0.19343;
    0.05657,0.15080,-0.02636; 0.19343,0.27636,0.09920;  % 8 Nearest and then 8 NN
    0.47340,-0.05670,0.15110; 0.40110,0.30670,0.22340;
    0.09890,0.19330,0.27660; -0.02660,0.05670,0.15110;
    0.52660,0.05670,-0.15110; 0.40110,0.19330,-0.27660;
    0.09890,0.30670,-0.22340; 0.02660,-0.05670,-0.15110; % 4 Neaest Er ions with positive charge
    0.125,0.0, 0.25; 0, 0.25, 0.125; 0.375, 0.0 -0.25; 0.5,0.25,-0.125];
%Zvector is the 110 direction, triangle center.
%Zvector=Atom_Er_position(2,:) + Atom_Er_position(3,:)-2*Atom_Er_position(1,:);
%Zunit=-1*Zvector/sqrt(sum(Zvector.^2)); %(which is atcually 10-1) then the defined symmetry Y-aix y= 010
Xunit =[1,0,1]/sqrt(2);
Yunit=[1,0,-1]/sqrt(2);
Zunit=[0,1,0]/sqrt(1);

% Xunit =[1,0,1]/sqrt(2);
% Yunit=[1,0,-1]/sqrt(2);
% Zunit=[0,1,0]/sqrt(1);

%Atom_Er_r = sqrt(sum(Atom_Er_position(1,:).^2));

%considering 16 O around Er atom (with extra 4 Next nearest bneighbour)
n=8;

%[phi,theta,Atom_O_R]=cart2sph(Atom_O_position(:,1),Atom_O_position(:,2),Atom_O_position(:,3));
for i =1:n
    Atom_ErObond(i,:) = Atom_O_position(i,:)-Atom_Er_position(1,:);
    Atom_ErObond_ProjPlane(i,:) = Atom_ErObond(i,:)-dot(Atom_ErObond(i,:),Zunit)*Zunit;
    Atom_ErO_R(i) = sqrt(sum(Atom_ErObond(i,:).^2));
    theta(i) = acos(dot(Atom_ErObond(i,:)/norm(Atom_ErObond(i,:)),Zunit));
        if(Atom_ErObond_ProjPlane(i,2) < 0)
            phi(i) = 2*pi-acos(dot(Atom_ErObond_ProjPlane(i,:)/norm(Atom_ErObond_ProjPlane(i,:)),Xunit));
        else
            phi(i) = acos(dot(Atom_ErObond_ProjPlane(i,:)/norm(Atom_ErObond_ProjPlane(i,:)),Xunit));
        end 
end

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
    for i = 1:n
        if i<=16
            temp(i)= q*(4*pi/(2*2+1)) *2*q/(Atom_ErO_R(i)^(2*2-1))  *ylm(2,-2*(j-1),theta(i),phi(i));
        else
            temp(i)= -q*(4*pi/(2*2+1)) * 3*q/(Atom_ErO_R(i)^(2*2-1))  *ylm(2,-2*(j-1),theta(i),phi(i));
        end
    end
    A_CF_para(j)=sum(temp);
end

for j=2 %B22
    for i = 1:n
        if i <=16
            temp(i)= q*(4*pi/(2*2+1)) *2*q/(Atom_ErO_R(i)^(3)) *1/sqrt(2)*(ylm(2,-2*(j-1),theta(i),phi(i))+ylm(2,2*(j-1),theta(i),phi(i)));
        else
            temp(i)= -q*(4*pi/(2*2+1)) *3*q/(Atom_ErO_R(i)^(3)) *1/sqrt(2)*(ylm(2,-2*(j-1),theta(i),phi(i))+ylm(2,2*(j-1),theta(i),phi(i)));
        end
    end
    A_CF_para(j)=sum(temp);
end

for j=3  %B40
    for i = 1:n
        if i<=16
            temp(i)=q*(4*pi/(2*4+1)) *2*q/((Atom_ErO_R(i))^(5)) * ylm(4,-2*(j-3),theta(i),phi(i));
        else
            temp(i)=-q*(4*pi/(2*4+1)) *3*q/((Atom_ErO_R(i))^(5)) * ylm(4,-2*(j-3),theta(i),phi(i));
        end
    end
    A_CF_para(j)=sum(temp);
end

for j=4:5 %B42 B44
    for i = 1:n
        if i<=16
            temp(i)=q*(4*pi/(2*4+1)) *2*q/((Atom_ErO_R(i))^(5)) *1/sqrt(2) * (ylm(4,-2*(j-3),theta(i),phi(i)) + (-1)^(2*(j-3)) * ylm(4,2*(j-3),theta(i),phi(i)));
        else
            temp(i)=-q*(4*pi/(2*4+1)) *3*q/((Atom_ErO_R(i))^(5)) *1/sqrt(2) * (ylm(4,-2*(j-3),theta(i),phi(i)) + (-1)^(2*(j-3)) * ylm(4,2*(j-3),theta(i),phi(i)));
        end
    end
    A_CF_para(j)=sum(temp);
end


for j=6 %B60
    for i =1:n
        if i<=16
            temp2(i)=q*(4*pi/(2*6+1)) *2*q/((Atom_ErO_R(i))^(7))  * ylm(6,-2*(j-6),theta(i),phi(i));
        else
            temp2(i)=-q*(4*pi/(2*6+1)) *3*q/((Atom_ErO_R(i))^(7))  * ylm(6,-2*(j-6),theta(i),phi(i));
        end
    end
    A_CF_para(j)=sum(temp2);
end

for j=7:9 %B62 64 66
    for i =1:n
        if i<=16
            temp2(i)=q*(4*pi/(2*6+1)) *2*q/((Atom_ErO_R(i))^(7)) *1/sqrt(2) * (ylm(6,-2*(j-6),theta(i),phi(i))+ (-1)^(2*(j-6))*ylm(6,2*(j-6),theta(i),phi(i)));
        else
            temp2(i)=-q*(4*pi/(2*6+1)) *3*q/((Atom_ErO_R(i))^(7)) *1/sqrt(2) * (ylm(6,-2*(j-6),theta(i),phi(i))+ (-1)^(2*(j-6))*ylm(6,2*(j-6),theta(i),phi(i)));
        end
    end
    A_CF_para(j)=sum(temp2);
end

alphaEr=4/(9*25*7);
betaEr = 2/(9*5*7*11*13);
gammaEr = 8/(27*7*11^2*13^2);
sigma=[0.460,0.0190,-0.0283];
%meanr=[0.773,1.677,8.431];
meanr = [0.750,1.49,6.52];
%meanr=[0.666,1.126,3.978]; % from https://onlinelibrary.wiley.com/doi/pdf/10.1002/qua.560080816
%meanr=[0.698,1.218,4.506]
for i=1:2;
    %B(i)=alphaEr*(1-sigma(1))*A_CF_para(i)*meanr(1);
    B(i)=alphaEr*A_CF_para(i)*meanr(1)*bohr_radius^2;
    
end
for i=3:5
    % B(i)=betaEr*(1-sigma(2))*A_CF_para(i)*meanr(2);
    B(i)=betaEr*A_CF_para(i)*meanr(2)*bohr_radius^4;
    
end
for i=6:9
    % B(i)=gammaEr*(1-sigma(3))*A_CF_para(i)*meanr(3);
    B(i)=gammaEr*A_CF_para(i)*meanr(3)*bohr_radius^6;
    
end
Coeff_Zna=[0.25*sqrt(5/pi),0.25*sqrt(15/pi),3/(sqrt(pi)*16),sqrt(5/pi)*3/8,sqrt(35/pi)*3/16,...
    sqrt(13/pi)/32,sqrt(2730/pi)/64,sqrt(13/(7*pi))*21/32,sqrt(26/(231*pi))*231/64];




BCoeff=1000*B.*Coeff_Zna*oneover4piepsilon0*unittrans;

[O20,O22,O40,O42,O43,O44,O60,O62,O63,O64,O66,Jx,Jy,Jz,Jplus,Jminus,Jsquare,Unit] = OperatorTotalmomentum(J);
Hcef_pointcharge_cal=BCoeff(1)*O20+BCoeff(2)*O22+BCoeff(3)*O40+BCoeff(4)*O42+BCoeff(5)*O44+BCoeff(6)*O60+BCoeff(7)*O62+BCoeff(8)*O64+BCoeff(9)*O66;% - 1000*g_J*bohr_magneton*(0.1e-3)*Jz;
%Hcef_pointcharge_cal = round(Hcef_pointcharge_cal,6);
[V_pointcharge_cal,E_pointcharge_cal] = eig(Hcef_pointcharge_cal,'Vector');
[Energy_pointcharge_cal,index]=sort(E_pointcharge_cal);
V_pointcharge_cal=V_pointcharge_cal(index,:);
Energy_pointcharge_cal = (Energy_pointcharge_cal + abs(min(E_pointcharge_cal(:,1))));

gz=abs(g_J*2*V_pointcharge_cal(:,1)'*Jz*V_pointcharge_cal(:,1));
gx=abs(g_J*2*V_pointcharge_cal(:,1)'*Jx*V_pointcharge_cal(:,2));
gy=abs(g_J*2*V_pointcharge_cal(:,1)'*Jy*V_pointcharge_cal(:,2));
if gz<0
   g_xy=V_pointcharge_cal(:,2)'*Jplus*V_pointcharge_cal(:,1);
else
    g_xy=V_pointcharge_cal(:,1)'*Jplus*V_pointcharge_cal(:,2);
end  
G_tensor_xy=[gz,g_xy];
G_tensor_abc=[gx,gy,gz];
[V_pointcharge_cal(:,1),V_pointcharge_cal(:,2)];
%Phi_ is V'*Jz*V when it is negative
%Phi+ is V'*Jz*V when it is positive
% then  g//=2g_J * V'(:,1)*Jz*V(:,1)   and g_|_ = g_J * (Phi+*Jplus*Phi_)
end

