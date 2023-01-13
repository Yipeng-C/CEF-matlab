clc

clear

n=50;

Solution_temp=zeros(n,9);

chi_sqaure_temp=zeros(n,1);

wave_function=zeros(16,n);

gxyz=zeros(n,4);

for i=1:n
    
    [Solution,Total_chi_to_see,gtensor,wave_function_temp]=Main_PointCharge_fitted_answer();

    Solution_temp(i,:)=Solution;

    chi_sqaure_temp(i)= Total_chi_to_see(3)
     
    wave_function(:,i) = wave_function_temp;
    
     gxyz(i,:) = gtensor
     

end



[chi_sqaure,index]=sort(chi_sqaure_temp);
Solution_temp=Solution_temp(index,:);
gxyz=gxyz(index,:);
wave_function=wave_function(:,index);


Coeff = fopen('Largestepruns.txt','w');
Chisquare = fopen('chi2.txt','w');
gtensors = fopen('g_tensor.txt','w');

fprintf(Coeff,'% 6.5e % 6.5e % 6.5e % 6.5e % 6.5e % 6.5e % 6.5e % 6.5e % 6.5e \n',Solution_temp);
fprintf(Chisquare,'% 6.5e\n',chi_sqaure);
fprintf(gtensors,'% 6.5e % 6.5e % 6.5e % 6.5e\n', gxyz);

fclose(Coeff);
fclose(Chisquare);
fclose(gtensors);
