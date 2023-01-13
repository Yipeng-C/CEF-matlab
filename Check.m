%function [BCoeff,chi2_energy_lvs,chi2_intensiy,chi2_initial,PeakIntensity,Peakposition]=Initialchi2()
bohr_magneton=5.7883818012e-5;   % with unit eV/T
g_J=6/5;

BCoeff=[0.1271	-0.4371	0.00066574	-0.0017	0.0033	1.03E-05	9.01E-05	5.03E-05	8.51E-06];%Point_charge_cal_NN();


Peakposition=[0;0;5.7204;5.7204;9.7955;9.7955;13.0862;13.0862;16.0332;16.0332; ...
    53.7546;53.7546;61.4186;61.4186;65.1326;65.1326];

PeakIntensity=[0.013658,0.0052646,0.00042553,0.00038729, ...
    0.00035254,0.00058745,0.00080699];

J=15/2;
[O20,O22,O40,O42,O43,O44,O60,O62,O63,O64,O66,Jx,Jy,Jz,Jplus,Jminus,Jsquare,Unit] = OperatorTotalmomentum(J);
Hcef_pointcharge_cal=BCoeff(1)*O20+BCoeff(2)*O22+BCoeff(3)*O40+BCoeff(4)*O42+BCoeff(5)*O44+BCoeff(6)*O60+BCoeff(7)*O62+BCoeff(8)*O64+BCoeff(9)*O66;% - 1000*g_J*bohr_magneton*1e-3*Jz;
%Hcef_pointcharge_cal=round(Hcef_pointcharge_cal,6);
[V_pointcharge_cal,E_pointcharge_cal] = eig(Hcef_pointcharge_cal,'Vector');
[E_pointcharge_cal,index]=sort(E_pointcharge_cal);
V_pointcharge_cal=V_pointcharge_cal(index,:);
E_pointcharge_cal = (E_pointcharge_cal + abs(min(E_pointcharge_cal(:,1))))
gz=abs(g_J*2*V_pointcharge_cal(:,1)'*Jz*V_pointcharge_cal(:,1));
gx=abs(g_J*2*V_pointcharge_cal(:,1)'*Jx*V_pointcharge_cal(:,2));
gy=-(g_J*2*V_pointcharge_cal(:,1)'*i*Jy*V_pointcharge_cal(:,2));
gplus=abs(g_J*V_pointcharge_cal(:,1)'*Jplus*V_pointcharge_cal(:,2));
G_tensor=[gx,gy,gz]
[V_pointcharge_cal(:,1),V_pointcharge_cal(:,2)]
