function     [Solution,Total_chi_to_see,gtensor,wave_function_temp]=Main_PointCharge_fitted_answer();


clc;
clear;
bohr_magneton=5.7883818012e-5;   % with unit eV/T
g_J=6/5;
[BCoeff,chi2_energy_lvs,chi2_intensiy,chi2_initial,PeakIntensity,Peakposition]=Initialchi2();
%BCoeff(9)=1e-7;

%Answer
%BCoeff=[0.135521741178622,-0.471112787614797,0.000166662471166580,-0.00176779261302444,0.00386196226999526,1.25048666863994e-05,6.51590971016271e-05,5.53729266104967e-05,2.58991839871017e-5];

%BCoeff = [-0.1389268882366492,0.047,0.0003368948619330090,1.173314220590409e-4,-0.002454432841458166,2.16391211615225e-06,2.69785000736794e-6,3.93724711468979e-05,00];

J=15/2;
[O20,O22,O40,O42,O43,O44,O60,O62,O63,O64,O66,Jx,Jy,Jz,Jplus,Jminus,Jsquare,Unit] = OperatorTotalmomentum(J);

for i=1:1000000
    if i<=500000
        Atemp=zeros(1,9);
        nonzerotemp= randi([1 9],1,1);
        Atemp(nonzerotemp) = 3*(-1 + 2*rand(1));
        Btemp=BCoeff + 0.5*(-1 + 2*rand(1,9)).*BCoeff;
    else if i>500000
            Btemp = BCoeff + 0.01*(-1 + 2*rand(1,9)).*BCoeff;
        end
    end
    Hcef = Btemp(1).*O20 + Btemp(2).*O22 + Btemp(3).*O40 + Btemp(4).*O42 + Btemp(5).*O44 + Btemp(6).*O60 + Btemp(7).*O62 + Btemp(8).*O64 + Btemp(9).*O66;% - 1000*g_J*bohr_magneton*(0.1e-6)*Jz;
    Hcef = round(Hcef,6);
    [V,E] = eig(Hcef,'Vector');
    [E,index]=sort(E);
    V=V(index,:);
    E = E + abs(min(E(:,1)));
    gz=abs(2*g_J*V(:,1)'*Jz*V(:,1));
    gx=abs(2*g_J*V(:,1)'*Jx*V(:,2));
    gy=abs(2*g_J*V(:,1)'*Jy*V(:,2));
if gz<0
   g_xy=2*g_J*V(:,2)'*Jplus*V(:,1);
else
    g_xy=2*g_J*V(:,1)'*Jplus*V(:,2);
end  
G_tensor=[gx,gy,gz,g_xy];
    
    
    scattering1=scattering_CEF(V(:,1),V(:,2),Jx,Jy,Jz)+scattering_CEF(V(:,2),V(:,1),Jx,Jy,Jz);  %3meV
    
    scattering2 = scattering_CEF(V(:,1),V(:,3),Jx,Jy,Jz) + scattering_CEF(V(:,1),V(:,4),Jx,Jy,Jz)...
        + scattering_CEF(V(:,2),V(:,3),Jx,Jy,Jz) + scattering_CEF(V(:,2),V(:,4),Jx,Jy,Jz); %3meV
    
    scattering3 = scattering_CEF(V(:,1),V(:,5),Jx,Jy,Jz) + scattering_CEF(V(:,1),V(:,6),Jx,Jy,Jz)...
        + scattering_CEF(V(:,2),V(:,5),Jx,Jy,Jz) + scattering_CEF(V(:,2),V(:,6),Jx,Jy,Jz);
    
    scattering4 = scattering_CEF(V(:,1),V(:,7),Jx,Jy,Jz) + scattering_CEF(V(:,1),V(:,8),Jx,Jy,Jz)...
        + scattering_CEF(V(:,2),V(:,7),Jx,Jy,Jz) + scattering_CEF(V(:,2),V(:,8),Jx,Jy,Jz);
    
    
    scattering5 = scattering_CEF(V(:,1),V(:,9),Jx,Jy,Jz) + scattering_CEF(V(:,1),V(:,10),Jx,Jy,Jz)...
        + scattering_CEF(V(:,2),V(:,9),Jx,Jy,Jz) + scattering_CEF(V(:,2),V(:,10),Jx,Jy,Jz);
    
    
    scattering6 = scattering_CEF(V(:,1),V(:,11),Jx,Jy,Jz) + scattering_CEF(V(:,1),V(:,12),Jx,Jy,Jz)...
        + scattering_CEF(V(:,2),V(:,11),Jx,Jy,Jz) + scattering_CEF(V(:,2),V(:,12),Jx,Jy,Jz);
    
    scattering7 = scattering_CEF(V(:,1),V(:,13),Jx,Jy,Jz) + scattering_CEF(V(:,1),V(:,14),Jx,Jy,Jz)...
        + scattering_CEF(V(:,2),V(:,13),Jx,Jy,Jz) + scattering_CEF(V(:,2),V(:,14),Jx,Jy,Jz);
    
    scattering8 = scattering_CEF(V(:,1),V(:,15),Jx,Jy,Jz) + scattering_CEF(V(:,1),V(:,16),Jx,Jy,Jz)...
        + scattering_CEF(V(:,2),V(:,15),Jx,Jy,Jz) + scattering_CEF(V(:,2),V(:,16),Jx,Jy,Jz);
    
    scattering9 = scattering_CEF(V(:,3),V(:,5),Jx,Jy,Jz) + scattering_CEF(V(:,3),V(:,6),Jx,Jy,Jz)...
        + scattering_CEF(V(:,4),V(:,5),Jx,Jy,Jz) + scattering_CEF(V(:,4),V(:,6),Jx,Jy,Jz);
    
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
    Energy_difference = (E - Peakposition)./Peakposition;
    chi2_E_temp =100 * sum(Energy_difference(3:16).^2);
    
    chi2_I_temp = 100*  sum(((expintensityratio(1:7) - calscattering(1:7))./expintensityratio(1:7)).^2);
    %chi2_I_HighmeV_temp= 100* (expintensityratio(5:6) - calscattering(5:6)).^2 / (expintensityratio(5:6)).^2;
    chi2_temp=chi2_E_temp+chi2_I_temp;%+chi2_I_HighmeV_temp;
    
       while chi2_temp < chi2_initial % for total chi2
            chi2_initial = chi2_temp;
            BCoeff=Btemp;
    %  while (chi2_I_temp < 2.5 && chi2_E_temp < chi2_energy_temp)
    %          chi2_energy_temp=chi2_E_temp;
    %          chi2_initial=chi2_temp;
    
%     while chi2_E_temp < chi2_energy_lvs; % for energy lvs only
%         chi2_energy_lvs=chi2_E_temp;
%         BCoeff=Btemp;
        clf;
        plot(E, 'b--o')
        E_fitted=E;
        V_fitted=V;
        Intensity_fitted=calscattering;
        Total_chi_to_see=[chi2_E_temp,chi2_I_temp,chi2_initial];
        gtensor=G_tensor;
        wave_function_temp=V_fitted(:,1);
        chi2_E_temp;
        data=[E_fitted,Peakposition];
        Solution=Btemp;
    end
end
%end
hold on
plot(Peakposition,'b')

end