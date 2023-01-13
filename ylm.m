function ylm = ylm (l,m,theta,phi)
%theta sit in [0,pi];
% Phi and theta are in radians
theta = theta(:);
phi = phi(:)';
Pn = legendre(l,cos(theta),'norm');
if m >= 0
    Pn = (-1)^m *Pn(abs(m)+1,:);
else
    Pn = Pn(abs(m)+1,:);
end
ylm = (1/sqrt(2*pi)) * Pn' * exp(i*m*phi);