dtheta = 0.01;
theta = [0:dtheta:2*pi];
dphi = 0.01;
phi = [0:dphi:pi];

[THETA,PHI] = meshgrid(theta,phi);
Y22 = 0.25*sqrt(15/(2*pi))*(sin(THETA).^2).*exp(2*i*PHI);
Y11 =-0.50*sqrt(3/(2*pi))*sin(THETA).*exp(i*PHI);
Y10 = 0.50*sqrt(3/pi)*cos(THETA);

integral = sum(sum(conj(Y22).*Y11.*Y10))*dtheta*dphi;


disp(integral);