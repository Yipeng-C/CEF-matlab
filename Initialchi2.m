function [BCoeff,chi2_energy_lvs,chi2_intensiy,chi2_initial,PeakIntensity,Peakposition]=Initialchi2()
bohr_magneton=5.7883818012e-5;   % with unit eV/T
g_J=6/5;
BCoeff = Point_charge_cal_NN();
%answer for PRB
%BCoeff=[0.127147143	-0.437091444	0.000665743	-0.001723674	0.003272632	1.03E-05	9.01E-05	5.03E-05	-8.51E-06];%Point_charge_cal_NN();

%for; area viogt fit 
%BCoeff=[0.137143691264374,-0.427632091141611,0.00013674036647264,-0.00219787927145609,0.00392591840989928,1.32004880328859e-05,7.12493329897568e-05,5.05035684760091e-05,3.08595909122684e-05];
%BCoeff=[0.116456100258849,-0.444954482767922,-8.42450453889297e-05,-0.00157925764555956,0.00541519539905552,1.44605715041516e-05,7.21856560576377e-05,4.06425704994395e-05,-1.54755700124284e-05];
%BCoeff=[0.135760319776703,-0.447708390368447,-3.96480556419359e-05,-0.00153611227099010,0.00527703229575165,1.44079153894841e-05,6.27479416022996e-05,4.88648020170999e-05,-3.24896073696224e-05];

%Peak Intensity solution
%BCoeff=[0.190688104505208,-0.406166645965357,0.000508134376225400,-0.00138930116767481,0.00513316694598045,1.50086925262041e-05,5.29499995817002e-05,3.50169775564476e-05,9.96616473165506e-07];

%Area solutionfor G fit result (answer)
%BCoeff=[0.135521741178622,-0.471112787614797,0.000166662471166580,-0.00176779261302444,0.00386196226999526,1.25048666863994e-05,6.51590971016271e-05,5.53729266104967e-05,2.58991839871017e-5];

% Point center peak intensity
%BCoeff=[0.139688100166598,-0.412138050355336,9.90961197823570e-05,-0.00170893166707160,0.00222987737711345,1.14493167786201e-05,9.39305547011383e-05,7.49356715117257e-05,5.63481593317837e-05];
%BCoeff=[0.110325373217127,-0.448601243394356,-5.39710360760268e-05,-0.000934876730535080,0.00520669639147618,1.03922711771856e-05,0.000102490950782354,1.82275153546413e-05,-1.51643289903760e-05];

Peakposition=[0;0;5.7204;5.7204;9.7955;9.7955;13.0862;13.0862;16.0332;16.0332; ...
    53.7546;53.7546;61.4186;61.4186;65.1326;65.1326];

 %  Peakposition=[0;0;5.616;5.616;9.6326;9.6326;12.801;12.801;15.748;15.748; ...
 %     54.708;54.708;62.372;62.372;66.086;66.086];
% 120meV zero point is 0.95343, for low energy lvs are 6.7105 and 10.912
% with 54.708 62.372 66.086 
% 30meV energy lvs zero point is -0.038219  5.6822, 9.7573,  13.048  15.995
% 
% Following is the area not the intensity. from fit(Answer)
PeakIntensity=[0.013658,0.0052646,0.00042553,0.00038729, ...
  0.00035254,0.00058745,0.00080699];

% add shoulder intensity(area) into (s1: 0.0045925, s2: 0.0018698)
%PeakIntensity=[0.0182505,0.0071344,0.00042553,0.00038729, ...
 %  0.00035254,0.00058745,0.00080699];


%following from Voigt area
%PeakIntensity=[0.012604,0.005074,0.0003014,0.0003052, ...
%   0.00024912,0.00068112,0.00081941];
%from Voigt peakintensity
% PeakIntensity =[0.0099116,0.0047067,0.00017727,0.00024064,...
%     3.8184e-5,0.0001358,0.00024978]; 


 %Peakposition=[0;0;5.6596;5.6596;9.6722;9.6722;12.8122;12.8122;15.7982;15.7982; ...
 %      53.7546;53.7546;61.4186;61.4186;65.1326;65.1326];


%Now using Intensity to fit the data.
%  PeakIntensity=[0.010041,0.0047477,0.00022422,0.00026768...
%      3.185e-5,0.00012839,0.00023648];



J=15/2;
[O20,O22,O40,O42,O43,O44,O60,O62,O63,O64,O66,Jx,Jy,Jz,Jplus,Jminus,Jsquare,Unit] = OperatorTotalmomentum(J);
Hcef_pointcharge_cal=BCoeff(1)*O20+BCoeff(2)*O22+BCoeff(3)*O40+BCoeff(4)*O42+BCoeff(5)*O44+BCoeff(6)*O60+BCoeff(7)*O62+BCoeff(8)*O64+BCoeff(9)*O66 ; %- 1000*g_J*bohr_magneton*1*Jz;
Hcef_pointcharge_cal=round(Hcef_pointcharge_cal,6);
[V_pointcharge_cal,E_pointcharge_cal] = eig(Hcef_pointcharge_cal,'Vector');
[E_pointcharge_cal,index]=sort(E_pointcharge_cal);
V_pointcharge_cal=V_pointcharge_cal(index,:);
E_pointcharge_cal = (E_pointcharge_cal + abs(min(E_pointcharge_cal(:,1))));
gz=abs(g_J*2*V_pointcharge_cal(:,1)'*Jz*V_pointcharge_cal(:,1));
gx=abs(g_J*2*V_pointcharge_cal(:,1)'*Jx*V_pointcharge_cal(:,2));
gy=-(g_J*2*V_pointcharge_cal(:,1)'*i*Jy*V_pointcharge_cal(:,2));
gplus=abs(g_J*V_pointcharge_cal(:,1)'*Jplus*V_pointcharge_cal(:,2));
G_tensor=[gx,gy,gz];
[V_pointcharge_cal(:,1),V_pointcharge_cal(:,2)];
scattering1 = scattering_CEF(V_pointcharge_cal(:,1),V_pointcharge_cal(:,2),Jx,Jy,Jz)+scattering_CEF(V_pointcharge_cal(:,2),V_pointcharge_cal(:,1),Jx,Jy,Jz);  %3meV

scattering2 = scattering_CEF(V_pointcharge_cal(:,1),V_pointcharge_cal(:,3),Jx,Jy,Jz) + scattering_CEF(V_pointcharge_cal(:,1),V_pointcharge_cal(:,4),Jx,Jy,Jz)...
    + scattering_CEF(V_pointcharge_cal(:,2),V_pointcharge_cal(:,3),Jx,Jy,Jz) + scattering_CEF(V_pointcharge_cal(:,2),V_pointcharge_cal(:,4),Jx,Jy,Jz); %3meV

scattering3 = scattering_CEF(V_pointcharge_cal(:,1),V_pointcharge_cal(:,5),Jx,Jy,Jz) + scattering_CEF(V_pointcharge_cal(:,1),V_pointcharge_cal(:,6),Jx,Jy,Jz)...
    + scattering_CEF(V_pointcharge_cal(:,2),V_pointcharge_cal(:,5),Jx,Jy,Jz) + scattering_CEF(V_pointcharge_cal(:,2),V_pointcharge_cal(:,6),Jx,Jy,Jz);

scattering4 = scattering_CEF(V_pointcharge_cal(:,1),V_pointcharge_cal(:,7),Jx,Jy,Jz) + scattering_CEF(V_pointcharge_cal(:,1),V_pointcharge_cal(:,8),Jx,Jy,Jz)...
    + scattering_CEF(V_pointcharge_cal(:,2),V_pointcharge_cal(:,7),Jx,Jy,Jz) + scattering_CEF(V_pointcharge_cal(:,2),V_pointcharge_cal(:,8),Jx,Jy,Jz);


scattering5 = scattering_CEF(V_pointcharge_cal(:,1),V_pointcharge_cal(:,9),Jx,Jy,Jz) + scattering_CEF(V_pointcharge_cal(:,1),V_pointcharge_cal(:,10),Jx,Jy,Jz)...
    + scattering_CEF(V_pointcharge_cal(:,2),V_pointcharge_cal(:,9),Jx,Jy,Jz) + scattering_CEF(V_pointcharge_cal(:,2),V_pointcharge_cal(:,10),Jx,Jy,Jz);


scattering6 = scattering_CEF(V_pointcharge_cal(:,1),V_pointcharge_cal(:,11),Jx,Jy,Jz) + scattering_CEF(V_pointcharge_cal(:,1),V_pointcharge_cal(:,12),Jx,Jy,Jz)...
    + scattering_CEF(V_pointcharge_cal(:,2),V_pointcharge_cal(:,11),Jx,Jy,Jz) + scattering_CEF(V_pointcharge_cal(:,2),V_pointcharge_cal(:,12),Jx,Jy,Jz);

scattering7 = scattering_CEF(V_pointcharge_cal(:,1),V_pointcharge_cal(:,13),Jx,Jy,Jz) + scattering_CEF(V_pointcharge_cal(:,1),V_pointcharge_cal(:,14),Jx,Jy,Jz)...
    + scattering_CEF(V_pointcharge_cal(:,2),V_pointcharge_cal(:,13),Jx,Jy,Jz) + scattering_CEF(V_pointcharge_cal(:,2),V_pointcharge_cal(:,14),Jx,Jy,Jz);

scattering8 = scattering_CEF(V_pointcharge_cal(:,1),V_pointcharge_cal(:,15),Jx,Jy,Jz) + scattering_CEF(V_pointcharge_cal(:,1),V_pointcharge_cal(:,16),Jx,Jy,Jz)...
    + scattering_CEF(V_pointcharge_cal(:,2),V_pointcharge_cal(:,15),Jx,Jy,Jz) + scattering_CEF(V_pointcharge_cal(:,2),V_pointcharge_cal(:,16),Jx,Jy,Jz);


scattering9 = scattering_CEF(V_pointcharge_cal(:,3),V_pointcharge_cal(:,5),Jx,Jy,Jz) + scattering_CEF(V_pointcharge_cal(:,3),V_pointcharge_cal(:,6),Jx,Jy,Jz)...
    + scattering_CEF(V_pointcharge_cal(:,4),V_pointcharge_cal(:,5),Jx,Jy,Jz) + scattering_CEF(V_pointcharge_cal(:,4),V_pointcharge_cal(:,6),Jx,Jy,Jz);

s1=scattering1/(scattering2);
s2=scattering2/(scattering2);
s3=scattering3/(scattering2);
s4=scattering4/(scattering2);
s5=scattering5/(scattering2);
s6=scattering6/(scattering8);
s7=scattering7/(scattering8);
s8=scattering8/(scattering8);
s9=scattering9/(scattering2);
calscattering=[s2,s3,s4,s5,s6,s7,s8];
expintensityratio(1:4)=PeakIntensity(1:4)./PeakIntensity(1);
expintensityratio(5:7)=PeakIntensity(5:7)./PeakIntensity(7);
%expintensityratio(3)=0.355;
expintensityratio;
pointcharge_cal_differ = (E_pointcharge_cal - Peakposition)./Peakposition;
chi2_energy_lvs = 100*sum(pointcharge_cal_differ(3:16).^2);

chi2_intensiy = 100*  sum(((expintensityratio(1:7) - calscattering(1:7))./expintensityratio(1:7)).^2);
chi2_initial=chi2_energy_lvs+chi2_intensiy;
data=[chi2_energy_lvs,chi2_intensiy,chi2_initial];
Data=[E_pointcharge_cal,Peakposition];
A=[calscattering;expintensityratio];
[V_pointcharge_cal(:,1),V_pointcharge_cal(:,2)];
% fileID = fopen('solution.txt','w');
% fprintf(fileID,'% 6.5e\n',BCoeff');
% fprintf(fileID,'% 6.5e',data');
% fclose(fileID);
end