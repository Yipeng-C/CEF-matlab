function S=scattering_CEF(i,j,Jx,Jy,Jz)
k=8.6173324*10^(-2);
x=ctranspose(i)*Jx*j;
y=ctranspose(i)*Jy*j;
z=ctranspose(i)*Jz*j;
S=conj(x)*x+conj(y)*y+conj(z)*z;
