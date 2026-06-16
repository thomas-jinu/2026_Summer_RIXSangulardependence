
% General outline of the code


clc;
clear all;

%Index for 3d orbitals
d3z2r2=1;
dx2y2=2;
dxy=3;
dzy=4;
dxz=5;

%Index for 2p orbitals {Strong LS coupling implies the J-mj states are used}
p3232=1;
p3212=2;
p32n12=3;
p32n32=4;

%Index for spins of 3d orbital electrons
s_up=1;
s_down=2;

%Polarization Index
p_sigma=1;
p_pi=2;
p_x=3;
p_y=4;
p_z=5;

%Spherical Harmonics index for d orbirtals
sY22=1;
sY21=2;
sY20=3;
sY2n1=4;
sY2n2=5;

%Spherical Harmonics index for the p orbitals
sY11=1;
sY10=2;
sY1n1=3;

%{
Experimental angles
delta - refers to the tilting angle of the substrate
theta - Half angle of the incoming and outgoing orientation. Either 45 ( 2theta=90 ) or 65 for 2theta=130;
%}
delta=0;
theta=0;  

%{
Incident and reflected angles 
Theta  gives information on the momentum that is passed into the 2d substrate. 
The phi angles gives information on crystal axis along which it is passed.
%}
theta_i=0;
phi_i=0;
theta_o=0;
phi_o=0;

%2-d space of the angles used for the spherical harmonics 

dphi_sh=.001;
dtheta_sh = .001;
Theta_sh=0:dtheta_sh:pi;
Phi_sh =0:dphi_sh:2*pi;

[phi_sh,theta_sh] = meshgrid(Phi_sh,Theta_sh);

dtheta = sin(theta_sh);   % sin(theta) is the jacobian of angle theta


%Defining the spherical harmonics used for orbital representations for integration theta_sh and phi_sh
Y22=(.25)*sqrt(15/(2*pi))*exp(2*1i*phi_sh).*(sin(theta_sh).*sin(theta_sh));
Y21=(-.5)*sqrt(15/(2*pi))*exp(1i*phi_sh).*(sin(theta_sh).*cos(theta_sh));
Y20=(.25)*sqrt(5/(pi))*(3*cos(theta_sh).*cos(theta_sh)-1);
Y2n1=(.5)*sqrt(15/(2*pi))*exp(-1i*phi_sh).*(sin(theta_sh).*cos(theta_sh));
Y2n2=(.25)*sqrt(15/(2*pi))*exp(-2*1i*phi_sh).*(sin(theta_sh).*sin(theta_sh));
Y11=(-.5)*sqrt(3/(2*pi))*sin(theta_sh).*exp(1i*phi_sh) ;
Y10= (.5)*sqrt(3/pi)*cos(theta_sh);
Y1n1=(.5)*sqrt(3/(2*pi))*exp(-1i*phi_sh).*sin(theta_sh);


%{
Define the polarization operators in each of the linear directions.
Not the T-matrix here which is entirely dependent on the intial and final
orbitals is the expectation value of 
polarization of light*position of the electron (e,r)

Here the postion vector r can be defined in terms of spherical harmonics.
%}

op_z=sqrt(4*pi/3).*Y10;
op_y=(1i)*sqrt(.5)*sqrt(4*pi/3).*(Y1n1+Y11);
op_x=sqrt(.5)*sqrt(4*pi/3).*(Y1n1-Y11);
op_pi=0;
op_sigma=0;

% In order to calculate the matrix elements for the dipole operator we have
% to have incorporate both the spin and the orbital part. The spin part
%  s represented in the Sz basis with |+> and |->  as the basis states. 


% The spins can have arbitrary orientations represented using the 
% Spin Orientation angles . 
theta_s = 0;
phi_s =0;

%Spinors which are the basis vectors |+> and |-> aligned along z axis
spinor_up = [1;0];
spinor_down = [0;1];

%Operators

os_up=  cos(theta_s/2)*exp(-1i*phi_s).*spinor_up + sin(theta_s/2).*spinor_down;
os_down = sin(theta_s/2)*exp(-1i*phi_s).*spinor_up - cos(theta_s/2).*spinor_down;


%The orbitals part include spherical harmonics. With the 3d orbital having
%Y2m and 2p having the Y1m harmonics. Because of the mixing of spin and
%orbital angular momentum of 2p states the following algorithm is defined
%for calculating the scattering cross-section
%The matrix is defined with the following index:
%Number of Y1m spherical harmonics (3)
%Number of Y2m spherical harmonics (5)
%Polarization of Light (x,y,z, sigma and pi) - (5)

SH = zeros(3,5,5);

% From Y11 to Y2m
% S = sum(A) returns the sum of the elements of A along the first array
% dimension whose size does not equal 1.Hence the double sum sums it along
% both rows and columns. 
% You could also use S = sum(A,"all") to sum over all matrix elements.

% A.*B does element by element multiplication necessary for doing numerical
% integration.

% Computing Simpson's weight for integration x = theta and y = phi
Nx = 3142;
Ny = 6284;
wx = ones(1,Nx); wx(2:2:Nx-1) = 4;  wx(3:2:Nx-2) = 2;
wy = ones(1,Ny); wy(2:2:Ny-1) = 4;  wy(3:2:Ny-2) = 2;

%Combining to create a single matrix
w = wy.'*wx; 
w = w';



SH(sY11,sY22,p_x)=sum(sum(conj(Y11).*op_x.*Y22.*dtheta.*w))*dtheta_sh*dphi_sh/9;
SH(sY11,sY21,p_x)=sum(sum(conj(Y11).*op_x.*Y21.*dtheta.*w))*dtheta_sh*dphi_sh/9;
SH(sY11,sY20,p_x)=sum(sum(conj(Y11).*op_x.*Y20.*dtheta.*w))*dtheta_sh*dphi_sh/9;
SH(sY11,sY2n1,p_x)=sum(sum(conj(Y11).*op_x.*Y2n1.*dtheta.*w))*dtheta_sh*dphi_sh/9;
SH(sY11,sY2n2,p_x)=sum(sum(conj(Y11).*op_x.*Y2n2.*dtheta.*w))*dtheta_sh*dphi_sh/9;

SH(sY11,sY22,p_y)=sum(sum(conj(Y11).*op_y.*Y22.*dtheta.*w))*dtheta_sh*dphi_sh/9;
SH(sY11,sY21,p_y)=sum(sum(conj(Y11).*op_y.*Y21.*dtheta.*w))*dtheta_sh*dphi_sh/9;
SH(sY11,sY20,p_y)=sum(sum(conj(Y11).*op_y.*Y20.*dtheta.*w))*dtheta_sh*dphi_sh/9;
SH(sY11,sY2n1,p_y)=sum(sum(conj(Y11).*op_y.*Y2n1.*dtheta.*w))*dtheta_sh*dphi_sh/9;
SH(sY11,sY2n2,p_y)=sum(sum(conj(Y11).*op_y.*Y2n2.*dtheta.*w))*dtheta_sh*dphi_sh/9;

SH(sY11,sY22,p_z)=sum(sum(conj(Y11).*op_z.*Y22.*dtheta.*w))*dtheta_sh*dphi_sh/9;
SH(sY11,sY21,p_z)=sum(sum(conj(Y11).*op_z.*Y21.*dtheta.*w))*dtheta_sh*dphi_sh/9;
SH(sY11,sY20,p_z)=sum(sum(conj(Y11).*op_z.*Y20.*dtheta.*w))*dtheta_sh*dphi_sh/9;
SH(sY11,sY2n1,p_z)=sum(sum(conj(Y11).*op_z.*Y2n1.*dtheta.*w))*dtheta_sh*dphi_sh/9;
SH(sY11,sY2n2,p_z)=sum(sum(conj(Y11).*op_z.*Y2n2.*dtheta.*w))*dtheta_sh*dphi_sh/9;

SH(sY11,sY22,p_pi)=sum(sum(conj(Y11).*op_pi.*Y22.*dtheta.*w))*dtheta_sh*dphi_sh/9;
SH(sY11,sY21,p_pi)=sum(sum(conj(Y11).*op_pi.*Y21.*dtheta.*w))*dtheta_sh*dphi_sh/9;
SH(sY11,sY20,p_pi)=sum(sum(conj(Y11).*op_pi.*Y20.*dtheta.*w))*dtheta_sh*dphi_sh/9;
SH(sY11,sY2n1,p_pi)=sum(sum(conj(Y11).*op_pi.*Y2n1.*dtheta.*w))*dtheta_sh*dphi_sh/9;
SH(sY11,sY2n2,p_pi)=sum(sum(conj(Y11).*op_pi.*Y2n2.*dtheta.*w))*dtheta_sh*dphi_sh/9;

SH(sY11,sY22,p_sigma)=sum(sum(conj(Y11).*op_sigma.*Y22.*dtheta.*w))*dtheta_sh*dphi_sh/9;
SH(sY11,sY21,p_sigma)=sum(sum(conj(Y11).*op_sigma.*Y21.*dtheta.*w))*dtheta_sh*dphi_sh/9;
SH(sY11,sY20,p_sigma)=sum(sum(conj(Y11).*op_sigma.*Y20.*dtheta.*w))*dtheta_sh*dphi_sh/9;
SH(sY11,sY2n1,p_sigma)=sum(sum(conj(Y11).*op_sigma.*Y2n1.*dtheta.*w))*dtheta_sh*dphi_sh/9;
SH(sY11,sY2n2,p_sigma)=sum(sum(conj(Y11).*op_sigma.*Y2n2.*dtheta.*w))*dtheta_sh*dphi_sh/9;

% From Y10 to Y2m

SH(sY10,sY22,p_x)=sum(sum(conj(Y10).*op_x.*Y22.*dtheta.*w))*dtheta_sh*dphi_sh/9;
SH(sY10,sY21,p_x)=sum(sum(conj(Y10).*op_x.*Y21.*dtheta.*w))*dtheta_sh*dphi_sh/9;
SH(sY10,sY20,p_x)=sum(sum(conj(Y10).*op_x.*Y20.*dtheta.*w))*dtheta_sh*dphi_sh/9;
SH(sY10,sY2n1,p_x)=sum(sum(conj(Y10).*op_x.*Y2n1.*dtheta.*w))*dtheta_sh*dphi_sh/9;
SH(sY10,sY2n2,p_x)=sum(sum(conj(Y10).*op_x.*Y2n2.*dtheta.*w))*dtheta_sh*dphi_sh/9;

SH(sY10,sY22,p_y)=sum(sum(conj(Y10).*op_y.*Y22.*dtheta.*w))*dtheta_sh*dphi_sh/9;
SH(sY10,sY21,p_y)=sum(sum(conj(Y10).*op_y.*Y21.*dtheta.*w))*dtheta_sh*dphi_sh/9;
SH(sY10,sY20,p_y)=sum(sum(conj(Y10).*op_y.*Y20.*dtheta.*w))*dtheta_sh*dphi_sh/9;
SH(sY10,sY2n1,p_y)=sum(sum(conj(Y10).*op_y.*Y2n1.*dtheta.*w))*dtheta_sh*dphi_sh/9;
SH(sY10,sY2n2,p_y)=sum(sum(conj(Y10).*op_y.*Y2n2.*dtheta.*w))*dtheta_sh*dphi_sh/9;

SH(sY10,sY22,p_z)=sum(sum(conj(Y10).*op_z.*Y22.*dtheta.*w))*dtheta_sh*dphi_sh/9;
SH(sY10,sY21,p_z)=sum(sum(conj(Y10).*op_z.*Y21.*dtheta.*w))*dtheta_sh*dphi_sh/9;
SH(sY10,sY20,p_z)=sum(sum(conj(Y10).*op_z.*Y20.*dtheta.*w))*dtheta_sh*dphi_sh/9;
SH(sY10,sY2n1,p_z)=sum(sum(conj(Y10).*op_z.*Y2n1.*dtheta.*w))*dtheta_sh*dphi_sh/9;
SH(sY10,sY2n2,p_z)=sum(sum(conj(Y10).*op_z.*Y2n2.*dtheta.*w))*dtheta_sh*dphi_sh/9;

SH(sY10,sY22,p_pi)=sum(sum(conj(Y10).*op_pi.*Y22.*dtheta.*w))*dtheta_sh*dphi_sh/9;
SH(sY10,sY21,p_pi)=sum(sum(conj(Y10).*op_pi.*Y21.*dtheta.*w))*dtheta_sh*dphi_sh/9;
SH(sY10,sY20,p_pi)=sum(sum(conj(Y10).*op_pi.*Y20.*dtheta.*w))*dtheta_sh*dphi_sh/9;
SH(sY10,sY2n1,p_pi)=sum(sum(conj(Y10).*op_pi.*Y2n1.*dtheta.*w))*dtheta_sh*dphi_sh/9;
SH(sY10,sY2n2,p_pi)=sum(sum(conj(Y10).*op_pi.*Y2n2.*dtheta.*w))*dtheta_sh*dphi_sh/9;

SH(sY10,sY22,p_sigma)=sum(sum(conj(Y10).*op_sigma.*Y22.*dtheta.*w))*dtheta_sh*dphi_sh/9;
SH(sY10,sY21,p_sigma)=sum(sum(conj(Y10).*op_sigma.*Y21.*dtheta.*w))*dtheta_sh*dphi_sh/9;
SH(sY10,sY20,p_sigma)=sum(sum(conj(Y10).*op_sigma.*Y20.*dtheta.*w))*dtheta_sh*dphi_sh/9;
SH(sY10,sY2n1,p_sigma)=sum(sum(conj(Y10).*op_sigma.*Y2n1.*dtheta.*w))*dtheta_sh*dphi_sh/9;
SH(sY10,sY2n2,p_sigma)=sum(sum(conj(Y10).*op_sigma.*Y2n2.*dtheta.*w))*dtheta_sh*dphi_sh/9;

% From Y1n1 to Y2m

SH(sY1n1,sY22,p_x)=sum(sum(conj(Y1n1).*op_x.*Y22.*dtheta.*w))*dtheta_sh*dphi_sh/9;
SH(sY1n1,sY21,p_x)=sum(sum(conj(Y1n1).*op_x.*Y21.*dtheta.*w))*dtheta_sh*dphi_sh/9;
SH(sY1n1,sY20,p_x)=sum(sum(conj(Y1n1).*op_x.*Y20.*dtheta.*w))*dtheta_sh*dphi_sh/9;
SH(sY1n1,sY2n1,p_x)=sum(sum(conj(Y1n1).*op_x.*Y2n1.*dtheta.*w))*dtheta_sh*dphi_sh/9;
SH(sY1n1,sY2n2,p_x)=sum(sum(conj(Y1n1).*op_x.*Y2n2.*dtheta.*w))*dtheta_sh*dphi_sh/9;

SH(sY1n1,sY22,p_y)=sum(sum(conj(Y1n1).*op_y.*Y22.*dtheta.*w))*dtheta_sh*dphi_sh/9;
SH(sY1n1,sY21,p_y)=sum(sum(conj(Y1n1).*op_y.*Y21.*dtheta.*w))*dtheta_sh*dphi_sh/9;
SH(sY1n1,sY20,p_y)=sum(sum(conj(Y1n1).*op_y.*Y20.*dtheta.*w))*dtheta_sh*dphi_sh/9;
SH(sY1n1,sY2n1,p_y)=sum(sum(conj(Y1n1).*op_y.*Y2n1.*dtheta.*w))*dtheta_sh*dphi_sh/9;
SH(sY1n1,sY2n2,p_y)=sum(sum(conj(Y1n1).*op_y.*Y2n2.*dtheta.*w))*dtheta_sh*dphi_sh/9;

SH(sY1n1,sY22,p_z)=sum(sum(conj(Y1n1).*op_z.*Y22.*dtheta.*w))*dtheta_sh*dphi_sh/9;
SH(sY1n1,sY21,p_z)=sum(sum(conj(Y1n1).*op_z.*Y21.*dtheta.*w))*dtheta_sh*dphi_sh/9;
SH(sY1n1,sY20,p_z)=sum(sum(conj(Y1n1).*op_z.*Y20.*dtheta.*w))*dtheta_sh*dphi_sh/9;
SH(sY1n1,sY2n1,p_z)=sum(sum(conj(Y1n1).*op_z.*Y2n1.*dtheta.*w))*dtheta_sh*dphi_sh/9;
SH(sY1n1,sY2n2,p_z)=sum(sum(conj(Y1n1).*op_z.*Y2n2.*dtheta.*w))*dtheta_sh*dphi_sh/9;

SH(sY1n1,sY22,p_pi)=sum(sum(conj(Y1n1).*op_pi.*Y22.*dtheta.*w))*dtheta_sh*dphi_sh/9;
SH(sY1n1,sY21,p_pi)=sum(sum(conj(Y1n1).*op_pi.*Y21.*dtheta.*w))*dtheta_sh*dphi_sh/9;
SH(sY1n1,sY20,p_pi)=sum(sum(conj(Y1n1).*op_pi.*Y20.*dtheta.*w))*dtheta_sh*dphi_sh/9;
SH(sY1n1,sY2n1,p_pi)=sum(sum(conj(Y1n1).*op_pi.*Y2n1.*dtheta.*w))*dtheta_sh*dphi_sh/9;
SH(sY1n1,sY2n2,p_pi)=sum(sum(conj(Y1n1).*op_pi.*Y2n2.*dtheta.*w))*dtheta_sh*dphi_sh/9;

SH(sY1n1,sY22,p_sigma)=sum(sum(conj(Y1n1).*op_sigma.*Y22.*dtheta.*w))*dtheta_sh*dphi_sh/9;
SH(sY1n1,sY21,p_sigma)=sum(sum(conj(Y1n1).*op_sigma.*Y21.*dtheta.*w))*dtheta_sh*dphi_sh/9;
SH(sY1n1,sY20,p_sigma)=sum(sum(conj(Y1n1).*op_sigma.*Y20.*dtheta.*w))*dtheta_sh*dphi_sh/9;
SH(sY1n1,sY2n1,p_sigma)=sum(sum(conj(Y1n1).*op_sigma.*Y2n1.*dtheta.*w))*dtheta_sh*dphi_sh/9;
SH(sY1n1,sY2n2,p_sigma)=sum(conj(Y1n1).*op_sigma.*Y2n2.*dtheta.*w,'all')*dtheta_sh*dphi_sh/9;


%{
Now we will define the orbitals; 
        
2p orbitals with strong LS coupling will have the following orbital and spin parts attached:

op32n32 = Y1n1 * |spin_down >

op32n12 = sqrt(1/3)*Y1n1*|spin_up> + sqrt(2/3)*Y10*|spin_down>;

op3212 = sqrt(2/3)*Y10*|spin_up> + sqrt(1/3)*Y11*|spin_down>;

op3232 = Y11*|spin_up>;


3d orbitals will have the L and S decoupled. So we will currently only
worry about the orbital part

o3dx2y2 = sqrt(1/2)*(Y22 + Y2n2);

o3dz2r2= (Y20);

o3dxy = (-1i/sqrt(2))*(Y22-Y2n2);

o3dzy = (-1i/sqrt(2))*(Y21+Y2n1);

o3dxz = sqrt(.5)*(Y21-Y2n1);
%}


%{
Defining the M matrix:
Since the matrix elements in the cross-section depend on the initial 
orbital, final orbital and the polarization we will define a 4-d matrix 
to store each value.
 
1)  First index gives the p orbitals from whcih transition occurs (No spin
accounted here since 2p orbitals have strong LS coupling) (4 values)
2)  Second index is the allowed 3d orbitals (5 values)
3)  Third index is the spins for electrons in 3d orbitals (up or down,2 values)
4)  Fourth index represents the polarization(5 values, sigma, pi and linearly
polarized light in x y and z directions)
%}

M =  zeros(4,5,2,5);

%Transition from op32n32 = op3232 = Y11*|spin_up> to 3d orbitals

disp("2p3232 to 3d orbitals, polarization p_x");

M(p3232,dx2y2,s_up,p_x)= sqrt(.5)* (SH(sY11,sY22,p_x)+SH(sY11,sY2n2,p_x))*(os_up'*os_up);
M(p3232,d3z2r2,s_up,p_x)= (SH(sY11,sY20,p_x))*(os_up'*os_up);
M(p3232,dxy,s_up,p_x)= (-1i*sqrt(.5))*(SH(sY11,sY22,p_x)-SH(sY11,sY2n2,p_x))*(os_up'*os_up);
M(p3232,dzy,s_up,p_x)= (-1i*sqrt(.5))*(SH(sY11,sY21,p_x)+SH(sY11,sY2n1,p_x))*(os_up'*os_up);
M(p3232,dxz,s_up,p_x)= sqrt(.5)*(SH(sY11,sY21,p_x)-SH(sY11,sY2n1,p_x))*(os_up'*os_up);

M(p3232,dx2y2,s_down,p_x)= sqrt(.5)* (SH(sY11,sY22,p_x)+SH(sY11,sY2n2,p_x))*(os_up'*os_down);
M(p3232,d3z2r2,s_down,p_x)= (SH(sY11,sY20,p_x))*(os_up'*os_down);
M(p3232,dxy,s_down,p_x)= (-1i*sqrt(.5))*(SH(sY11,sY22,p_x)-SH(sY11,sY2n2,p_x))*(os_up'*os_down);
M(p3232,dzy,s_down,p_x)= (-1i*sqrt(.5))*(SH(sY11,sY21,p_x)+SH(sY11,sY2n1,p_x))*(os_up'*os_down);
M(p3232,dxz,s_down,p_x)= sqrt(.5)*(SH(sY11,sY21,p_x)-SH(sY11,sY2n1,p_x))*(os_up'*os_down);


disp("2p3232 3dx2y2_up");
disp (M(p3232,dx2y2,s_up,p_x)); 
disp("2p3232 3dx2y2_down");
disp (M(p3232,dx2y2,s_down,p_x)); 
disp("2p3232 3d3z2r2_up");
disp (M(p3232,d3z2r2,s_up,p_x)); 
disp("2p3232 3d3z2r2_down");
disp (M(p3232,d3z2r2,s_down,p_x)); 
disp("2p3232 3dxy_up");
disp (M(p3232,dxy,s_up,p_x)); 
disp("2p3232 3dxy_down");
disp (M(p3232,dxy,s_down,p_x));
disp("2p3232 3dzy_up");
disp (M(p3232,dzy,s_up,p_x));
disp("2p3232 3dzy_down");
disp (M(p3232,dzy,s_down,p_x));
disp("2p3232 3dxz_up");
disp (M(p3232,dxz,s_up,p_x)); 
disp("2p3232 3dxz_down");
disp (M(p3232,dxz,s_down,p_x)); 



disp("2p3232 to 3d orbitals, polarization p_y");

M(p3232,dx2y2,s_up,p_y)= sqrt(.5)* (SH(sY11,sY22,p_y)+SH(sY11,sY2n2,p_y))*(os_up'*os_up);
M(p3232,d3z2r2,s_up,p_y)= (SH(sY11,sY20,p_y))*(os_up'*os_up);
M(p3232,dxy,s_up,p_y)= (-1i*sqrt(.5))*(SH(sY11,sY22,p_y)-SH(sY11,sY2n2,p_y))*(os_up'*os_up);
M(p3232,dzy,s_up,p_y)= (-1i*sqrt(.5))*(SH(sY11,sY21,p_y)+SH(sY11,sY2n1,p_y))*(os_up'*os_up);
M(p3232,dxz,s_up,p_y)= sqrt(.5)*(SH(sY11,sY21,p_y)-SH(sY11,sY2n1,p_y))*(os_up'*os_up);

M(p3232,dx2y2,s_down,p_y)= sqrt(.5)* (SH(sY11,sY22,p_y)+SH(sY11,sY2n2,p_y))*(os_up'*os_down);
M(p3232,d3z2r2,s_down,p_y)= (SH(sY11,sY20,p_y))*(os_up'*os_down);
M(p3232,dxy,s_down,p_y)= (-1i*sqrt(.5))*(SH(sY11,sY22,p_y)-SH(sY11,sY2n2,p_y))*(os_up'*os_down);
M(p3232,dzy,s_down,p_y)= (-1i*sqrt(.5))*(SH(sY11,sY21,p_y)+SH(sY11,sY2n1,p_y))*(os_up'*os_down);
M(p3232,dxz,s_down,p_y)= sqrt(.5)*(SH(sY11,sY21,p_y)-SH(sY11,sY2n1,p_y))*(os_up'*os_down);


disp("2p3232 3dx2y2_up");
disp (M(p3232,dx2y2,s_up,p_y)); 
disp("2p3232 3dx2y2_down");
disp (M(p3232,dx2y2,s_down,p_y)); 
disp("2p3232 3d3z2r2_up");
disp (M(p3232,d3z2r2,s_up,p_y)); 
disp("2p3232 3d3z2r2_down");
disp (M(p3232,d3z2r2,s_down,p_y)); 
disp("2p3232 3dxy_up");
disp (M(p3232,dxy,s_up,p_y)); 
disp("2p3232 3dxy_down");
disp (M(p3232,dxy,s_down,p_y));
disp("2p3232 3dzy_up");
disp (M(p3232,dzy,s_up,p_y));
disp("2p3232 3dzy_down");
disp (M(p3232,dzy,s_down,p_y));
disp("2p3232 3dxz_up");
disp (M(p3232,dxz,s_up,p_y)); 
disp("2p3232 3dxz_down");
disp (M(p3232,dxz,s_down,p_y)); 



disp("2p3232 to 3d orbitals, polarization p_z");

M(p3232,dx2y2,s_up,p_z)= sqrt(.5)* (SH(sY11,sY22,p_z)+SH(sY11,sY2n2,p_z))*(os_up'*os_up);
M(p3232,d3z2r2,s_up,p_z)= (SH(sY11,sY20,p_z))*(os_up'*os_up);
M(p3232,dxy,s_up,p_z)= (-1i*sqrt(.5))*(SH(sY11,sY22,p_z)-SH(sY11,sY2n2,p_z))*(os_up'*os_up);
M(p3232,dzy,s_up,p_z)= (-1i*sqrt(.5))*(SH(sY11,sY21,p_z)+SH(sY11,sY2n1,p_z))*(os_up'*os_up);
M(p3232,dxz,s_up,p_z)= sqrt(.5)*(SH(sY11,sY21,p_z)-SH(sY11,sY2n1,p_z))*(os_up'*os_up);

M(p3232,dx2y2,s_down,p_z)= sqrt(.5)* (SH(sY11,sY22,p_z)+SH(sY11,sY2n2,p_z))*(os_up'*os_down);
M(p3232,d3z2r2,s_down,p_z)= (SH(sY11,sY20,p_z))*(os_up'*os_down);
M(p3232,dxy,s_down,p_z)= (-1i*sqrt(.5))*(SH(sY11,sY22,p_z)-SH(sY11,sY2n2,p_z))*(os_up'*os_down);
M(p3232,dzy,s_down,p_z)= (-1i*sqrt(.5))*(SH(sY11,sY21,p_z)+SH(sY11,sY2n1,p_z))*(os_up'*os_down);
M(p3232,dxz,s_down,p_z)= sqrt(.5)*(SH(sY11,sY21,p_z)-SH(sY11,sY2n1,p_z))*(os_up'*os_down);


disp("2p3232 3dx2y2_up");
disp (M(p3232,dx2y2,s_up,p_z)); 
disp("2p3232 3dx2y2_down");
disp (M(p3232,dx2y2,s_down,p_z)); 
disp("2p3232 3d3z2r2_up");
disp (M(p3232,d3z2r2,s_up,p_z)); 
disp("2p3232 3d3z2r2_down");
disp (M(p3232,d3z2r2,s_down,p_z)); 
disp("2p3232 3dxy_up");
disp (M(p3232,dxy,s_up,p_z)); 
disp("2p3232 3dxy_down");
disp (M(p3232,dxy,s_down,p_z));
disp("2p3232 3dzy_up");
disp (M(p3232,dzy,s_up,p_z));
disp("2p3232 3dzy_down");
disp (M(p3232,dzy,s_down,p_z));
disp("2p3232 3dxz_up");
disp (M(p3232,dxz,s_up,p_z)); 
disp("2p3232 3dxz_down");
disp (M(p3232,dxz,s_down,p_z)); 


% 2p32n32 rto 3d orbitals

disp("2p32n32 to 3d orbitals, polarization p_x");

M(p32n32,dx2y2,s_up,p_x)= sqrt(.5)* (SH(sY1n1,sY22,p_x)+SH(sY1n1,sY2n2,p_x))*(os_down'*os_up);
M(p32n32,d3z2r2,s_up,p_x)= (SH(sY1n1,sY20,p_x))*(os_down'*os_up);
M(p32n32,dxy,s_up,p_x)= (-1i*sqrt(.5))*(SH(sY1n1,sY22,p_x)-SH(sY1n1,sY2n2,p_x))*(os_down'*os_up);
M(p32n32,dzy,s_up,p_x)= (-1i*sqrt(.5))*(SH(sY1n1,sY21,p_x)+SH(sY1n1,sY2n1,p_x))*(os_down'*os_up);
M(p32n32,dxz,s_up,p_x)= sqrt(.5)*(SH(sY1n1,sY21,p_x)-SH(sY1n1,sY2n1,p_x))*(os_down'*os_up);

M(p32n32,dx2y2,s_down,p_x)= sqrt(.5)* (SH(sY1n1,sY22,p_x)+SH(sY1n1,sY2n2,p_x))*(os_down'*os_down);
M(p32n32,d3z2r2,s_down,p_x)= (SH(sY1n1,sY20,p_x))*(os_down'*os_down);
M(p32n32,dxy,s_down,p_x)= (-1i*sqrt(.5))*(SH(sY1n1,sY22,p_x)-SH(sY1n1,sY2n2,p_x))*(os_down'*os_down);
M(p32n32,dzy,s_down,p_x)= (-1i*sqrt(.5))*(SH(sY1n1,sY21,p_x)+SH(sY1n1,sY2n1,p_x))*(os_down'*os_down);
M(p32n32,dxz,s_down,p_x)= sqrt(.5)*(SH(sY1n1,sY21,p_x)-SH(sY1n1,sY2n1,p_x))*(os_down'*os_down);


disp("2p32n32 3dx2y2_up");
disp (M(p32n32,dx2y2,s_up,p_x)); 
disp("2p32n32 3dx2y2_down");
disp (M(p32n32,dx2y2,s_down,p_x)); 
disp("2p32n32 3d3z2r2_up");
disp (M(p32n32,d3z2r2,s_up,p_x)); 
disp("2p32n32 3d3z2r2_down");
disp (M(p32n32,d3z2r2,s_down,p_x)); 
disp("2p32n32 3dxy_up");
disp (M(p32n32,dxy,s_up,p_x)); 
disp("2p32n32 3dxy_down");
disp (M(p32n32,dxy,s_down,p_x));
disp("2p32n32 3dzy_up");
disp (M(p32n32,dzy,s_up,p_x));
disp("2p32n32 3dzy_down");
disp (M(p32n32,dzy,s_down,p_x));
disp("2p32n32 3dxz_up");
disp (M(p32n32,dxz,s_up,p_x)); 
disp("2p32n32 3dxz_down");
disp (M(p32n32,dxz,s_down,p_x)); 



disp("2p32n32 to 3d orbitals, polarization p_y");

M(p32n32,dx2y2,s_up,p_y)= sqrt(.5)* (SH(sY1n1,sY22,p_y)+SH(sY1n1,sY2n2,p_y))*(os_down'*os_up);
M(p32n32,d3z2r2,s_up,p_y)= (SH(sY1n1,sY20,p_y))*(os_down'*os_up);
M(p32n32,dxy,s_up,p_y)= (-1i*sqrt(.5))*(SH(sY1n1,sY22,p_y)-SH(sY1n1,sY2n2,p_y))*(os_down'*os_up);
M(p32n32,dzy,s_up,p_y)= (-1i*sqrt(.5))*(SH(sY1n1,sY21,p_y)+SH(sY1n1,sY2n1,p_y))*(os_down'*os_up);
M(p32n32,dxz,s_up,p_y)= sqrt(.5)*(SH(sY1n1,sY21,p_y)-SH(sY1n1,sY2n1,p_y))*(os_down'*os_up);

M(p32n32,dx2y2,s_down,p_y)= sqrt(.5)* (SH(sY1n1,sY22,p_y)+SH(sY1n1,sY2n2,p_y))*(os_down'*os_down);
M(p32n32,d3z2r2,s_down,p_y)= (SH(sY1n1,sY20,p_y))*(os_down'*os_down);
M(p32n32,dxy,s_down,p_y)= (-1i*sqrt(.5))*(SH(sY1n1,sY22,p_y)-SH(sY1n1,sY2n2,p_y))*(os_down'*os_down);
M(p32n32,dzy,s_down,p_y)= (-1i*sqrt(.5))*(SH(sY1n1,sY21,p_y)+SH(sY1n1,sY2n1,p_y))*(os_down'*os_down);
M(p32n32,dxz,s_down,p_y)= sqrt(.5)*(SH(sY1n1,sY21,p_y)-SH(sY1n1,sY2n1,p_y))*(os_down'*os_down);


disp("2p32n32 3dx2y2_up");
disp (M(p32n32,dx2y2,s_up,p_y)); 
disp("2p32n32 3dx2y2_down");
disp (M(p32n32,dx2y2,s_down,p_y)); 
disp("2p32n32 3d3z2r2_up");
disp (M(p32n32,d3z2r2,s_up,p_y)); 
disp("2p32n32 3d3z2r2_down");
disp (M(p32n32,d3z2r2,s_down,p_y)); 
disp("2p32n32 3dxy_up");
disp (M(p32n32,dxy,s_up,p_y)); 
disp("2p32n32 3dxy_down");
disp (M(p32n32,dxy,s_down,p_y));
disp("2p32n32 3dzy_up");
disp (M(p32n32,dzy,s_up,p_y));
disp("2p32n32 3dzy_down");
disp (M(p32n32,dzy,s_down,p_y));
disp("2p32n32 3dxz_up");
disp (M(p32n32,dxz,s_up,p_y)); 
disp("2p32n32 3dxz_down");
disp (M(p32n32,dxz,s_down,p_y)); 



disp("2p32n32 to 3d orbitals, polarization p_z");

M(p32n32,dx2y2,s_up,p_z)= sqrt(.5)* (SH(sY1n1,sY22,p_z)+SH(sY1n1,sY2n2,p_z))*(os_down'*os_up);
M(p32n32,d3z2r2,s_up,p_z)= (SH(sY1n1,sY20,p_z))*(os_down'*os_up);
M(p32n32,dxy,s_up,p_z)= (-1i*sqrt(.5))*(SH(sY1n1,sY22,p_z)-SH(sY1n1,sY2n2,p_z))*(os_down'*os_up);
M(p32n32,dzy,s_up,p_z)= (-1i*sqrt(.5))*(SH(sY1n1,sY21,p_z)+SH(sY1n1,sY2n1,p_z))*(os_down'*os_up);
M(p32n32,dxz,s_up,p_z)= sqrt(.5)*(SH(sY1n1,sY21,p_z)-SH(sY1n1,sY2n1,p_z))*(os_down'*os_up);

M(p32n32,dx2y2,s_down,p_z)= sqrt(.5)* (SH(sY1n1,sY22,p_z)+SH(sY1n1,sY2n2,p_z))*(os_down'*os_down);
M(p32n32,d3z2r2,s_down,p_z)= (SH(sY1n1,sY20,p_z))*(os_down'*os_down);
M(p32n32,dxy,s_down,p_z)= (-1i*sqrt(.5))*(SH(sY1n1,sY22,p_z)-SH(sY1n1,sY2n2,p_z))*(os_down'*os_down);
M(p32n32,dzy,s_down,p_z)= (-1i*sqrt(.5))*(SH(sY1n1,sY21,p_z)+SH(sY1n1,sY2n1,p_z))*(os_down'*os_down);
M(p32n32,dxz,s_down,p_z)= sqrt(.5)*(SH(sY1n1,sY21,p_z)-SH(sY1n1,sY2n1,p_z))*(os_down'*os_down);


disp("2p32n32 3dx2y2_up");
disp (M(p32n32,dx2y2,s_up,p_z)); 
disp("2p32n32 3dx2y2_down");
disp (M(p32n32,dx2y2,s_down,p_z)); 
disp("2p32n32 3d3z2r2_up");
disp (M(p32n32,d3z2r2,s_up,p_z)); 
disp("2p32n32 3d3z2r2_down");
disp (M(p32n32,d3z2r2,s_down,p_z)); 
disp("2p32n32 3dxy_up");
disp (M(p32n32,dxy,s_up,p_z)); 
disp("2p32n32 3dxy_down");
disp (M(p32n32,dxy,s_down,p_z));
disp("2p32n32 3dzy_up");
disp (M(p32n32,dzy,s_up,p_z));
disp("2p32n32 3dzy_down");
disp (M(p32n32,dzy,s_down,p_z));
disp("2p32n32 3dxz_up");
disp (M(p32n32,dxz,s_up,p_z)); 
disp("2p32n32 3dxz_down");
disp (M(p32n32,dxz,s_down,p_z));


% 2p3212 to 3d orbitals

disp("2p3212 to 3d orbitals, polarization p_x");

M(p3212,dx2y2,s_up,p_x)= sqrt(.5)* sqrt(2/3)*(SH(sY10,sY22,p_x)+SH(sY10,sY2n2,p_x))*(os_up'*os_up) + sqrt(.5)* sqrt(1/3)*(SH(sY11,sY22,p_x)+SH(sY11,sY2n2,p_x))*(os_down'*os_up);
M(p3212,d3z2r2,s_up,p_x)= sqrt(2/3)*SH(sY10,sY20,p_x)*(os_up'*os_up) + sqrt(1/3)*SH(sY11,sY20,p_x)*(os_down'*os_up);
M(p3212,dxy,s_up,p_x)= (-1i*sqrt(.5))*sqrt(2/3)*(SH(sY10,sY22,p_x)-SH(sY10,sY2n2,p_x))*(os_up'*os_up) + (-1i*sqrt(.5))*sqrt(1/3)*(SH(sY11,sY22,p_x)-SH(sY11,sY2n2,p_x))*(os_down'*os_up);
M(p3212,dzy,s_up,p_x)= (-1i*sqrt(.5))*sqrt(2/3)*(SH(sY10,sY21,p_x)+SH(sY10,sY2n1,p_x))*(os_up'*os_up) + (-1i*sqrt(.5))*sqrt(1/3)*(SH(sY11,sY21,p_x)+SH(sY11,sY2n1,p_x))*(os_down'*os_up);
M(p3212,dxz,s_up,p_x)= sqrt(.5)*sqrt(2/3)*(SH(sY10,sY21,p_x)-SH(sY10,sY2n1,p_x))*(os_up'*os_up) + sqrt(.5)*sqrt(1/3)*(SH(sY11,sY21,p_x)-SH(sY11,sY2n1,p_x))*(os_down'*os_up);

M(p3212,dx2y2,s_down,p_x)= sqrt(.5)* sqrt(2/3)*(SH(sY10,sY22,p_x)+SH(sY10,sY2n2,p_x))*(os_up'*os_down) + sqrt(.5)* sqrt(1/3)*(SH(sY11,sY22,p_x)+SH(sY11,sY2n2,p_x))*(os_down'*os_down);
M(p3212,d3z2r2,s_down,p_x)= sqrt(2/3)*SH(sY10,sY20,p_x)*(os_up'*os_down) + sqrt(1/3)*SH(sY11,sY20,p_x)*(os_down'*os_down);
M(p3212,dxy,s_down,p_x)= (-1i*sqrt(.5))*sqrt(2/3)*(SH(sY10,sY22,p_x)-SH(sY10,sY2n2,p_x))*(os_up'*os_down) + (-1i*sqrt(.5))*sqrt(1/3)*(SH(sY11,sY22,p_x)-SH(sY11,sY2n2,p_x))*(os_down'*os_down);
M(p3212,dzy,s_down,p_x)= (-1i*sqrt(.5))*sqrt(2/3)*(SH(sY10,sY21,p_x)+SH(sY10,sY2n1,p_x))*(os_up'*os_down) + (-1i*sqrt(.5))*sqrt(1/3)*(SH(sY11,sY21,p_x)+SH(sY11,sY2n1,p_x))*(os_down'*os_down);
M(p3212,dxz,s_down,p_x)= sqrt(.5)*sqrt(2/3)*(SH(sY10,sY21,p_x)-SH(sY10,sY2n1,p_x))*(os_up'*os_down) + sqrt(.5)*sqrt(1/3)*(SH(sY11,sY21,p_x)-SH(sY11,sY2n1,p_x))*(os_down'*os_down);


disp("2p3212 3dx2y2_up");
disp (M(p3212,dx2y2,s_up,p_x)); 
disp("2p3212 3dx2y2_down");
disp (M(p3212,dx2y2,s_down,p_x)); 
disp("2p3212 3d3z2r2_up");
disp (M(p3212,d3z2r2,s_up,p_x)); 
disp("2p3212 3d3z2r2_down");
disp (M(p3212,d3z2r2,s_down,p_x)); 
disp("2p3212 3dxy_up");
disp (M(p3212,dxy,s_up,p_x)); 
disp("2p3212 3dxy_down");
disp (M(p3212,dxy,s_down,p_x));
disp("2p3212 3dzy_up");
disp (M(p3212,dzy,s_up,p_x));
disp("2p3212 3dzy_down");
disp (M(p3212,dzy,s_down,p_x));
disp("2p3212 3dxz_up");
disp (M(p3212,dxz,s_up,p_x)); 
disp("2p3212 3dxz_down");
disp (M(p3212,dxz,s_down,p_x)); 


disp("2p3212 to 3d orbitals, polarization p_y");

M(p3212,dx2y2,s_up,p_y)= sqrt(.5)* sqrt(2/3)*(SH(sY10,sY22,p_y)+SH(sY10,sY2n2,p_y))*(os_up'*os_up) + sqrt(.5)* sqrt(1/3)*(SH(sY11,sY22,p_y)+SH(sY11,sY2n2,p_y))*(os_down'*os_up);
M(p3212,d3z2r2,s_up,p_y)= sqrt(2/3)*SH(sY10,sY20,p_y)*(os_up'*os_up) + sqrt(1/3)*SH(sY11,sY20,p_y)*(os_down'*os_up);
M(p3212,dxy,s_up,p_y)= (-1i*sqrt(.5))*sqrt(2/3)*(SH(sY10,sY22,p_y)-SH(sY10,sY2n2,p_y))*(os_up'*os_up) + (-1i*sqrt(.5))*sqrt(1/3)*(SH(sY11,sY22,p_y)-SH(sY11,sY2n2,p_y))*(os_down'*os_up);
M(p3212,dzy,s_up,p_y)= (-1i*sqrt(.5))*sqrt(2/3)*(SH(sY10,sY21,p_y)+SH(sY10,sY2n1,p_y))*(os_up'*os_up) + (-1i*sqrt(.5))*sqrt(1/3)*(SH(sY11,sY21,p_y)+SH(sY11,sY2n1,p_y))*(os_down'*os_up);
M(p3212,dxz,s_up,p_y)= sqrt(.5)*sqrt(2/3)*(SH(sY10,sY21,p_y)-SH(sY10,sY2n1,p_y))*(os_up'*os_up) + sqrt(.5)*sqrt(1/3)*(SH(sY11,sY21,p_y)-SH(sY11,sY2n1,p_y))*(os_down'*os_up);

M(p3212,dx2y2,s_down,p_y)= sqrt(.5)* sqrt(2/3)*(SH(sY10,sY22,p_y)+SH(sY10,sY2n2,p_y))*(os_up'*os_down) + sqrt(.5)* sqrt(1/3)*(SH(sY11,sY22,p_y)+SH(sY11,sY2n2,p_y))*(os_down'*os_down);
M(p3212,d3z2r2,s_down,p_y)= sqrt(2/3)*SH(sY10,sY20,p_y)*(os_up'*os_down) + sqrt(1/3)*SH(sY11,sY20,p_y)*(os_down'*os_down);
M(p3212,dxy,s_down,p_y)= (-1i*sqrt(.5))*sqrt(2/3)*(SH(sY10,sY22,p_y)-SH(sY10,sY2n2,p_y))*(os_up'*os_down) + (-1i*sqrt(.5))*sqrt(1/3)*(SH(sY11,sY22,p_y)-SH(sY11,sY2n2,p_y))*(os_down'*os_down);
M(p3212,dzy,s_down,p_y)= (-1i*sqrt(.5))*sqrt(2/3)*(SH(sY10,sY21,p_y)+SH(sY10,sY2n1,p_y))*(os_up'*os_down) + (-1i*sqrt(.5))*sqrt(1/3)*(SH(sY11,sY21,p_y)+SH(sY11,sY2n1,p_y))*(os_down'*os_down);
M(p3212,dxz,s_down,p_y)= sqrt(.5)*sqrt(2/3)*(SH(sY10,sY21,p_y)-SH(sY10,sY2n1,p_y))*(os_up'*os_down) + sqrt(.5)*sqrt(1/3)*(SH(sY11,sY21,p_y)-SH(sY11,sY2n1,p_y))*(os_down'*os_down);


disp("2p3212 3dx2y2_up");
disp (M(p3212,dx2y2,s_up,p_y)); 
disp("2p3212 3dx2y2_down");
disp (M(p3212,dx2y2,s_down,p_y)); 
disp("2p3212 3d3z2r2_up");
disp (M(p3212,d3z2r2,s_up,p_y)); 
disp("2p3212 3d3z2r2_down");
disp (M(p3212,d3z2r2,s_down,p_y)); 
disp("2p3212 3dxy_up");
disp (M(p3212,dxy,s_up,p_y)); 
disp("2p3212 3dxy_down");
disp (M(p3212,dxy,s_down,p_y));
disp("2p3212 3dzy_up");
disp (M(p3212,dzy,s_up,p_y));
disp("2p3212 3dzy_down");
disp (M(p3212,dzy,s_down,p_y));
disp("2p3212 3dxz_up");
disp (M(p3212,dxz,s_up,p_y)); 
disp("2p3212 3dxz_down");
disp (M(p3212,dxz,s_down,p_y)); 


disp("2p3212 to 3d orbitals, polarization p_z");

M(p3212,dx2y2,s_up,p_z)= sqrt(.5)* sqrt(2/3)*(SH(sY10,sY22,p_z)+SH(sY10,sY2n2,p_z))*(os_up'*os_up) + sqrt(.5)* sqrt(1/3)*(SH(sY11,sY22,p_z)+SH(sY11,sY2n2,p_z))*(os_down'*os_up);
M(p3212,d3z2r2,s_up,p_z)= sqrt(2/3)*SH(sY10,sY20,p_z)*(os_up'*os_up) + sqrt(1/3)*SH(sY11,sY20,p_z)*(os_down'*os_up);
M(p3212,dxy,s_up,p_z)= (-1i*sqrt(.5))*sqrt(2/3)*(SH(sY10,sY22,p_z)-SH(sY10,sY2n2,p_z))*(os_up'*os_up) + (-1i*sqrt(.5))*sqrt(1/3)*(SH(sY11,sY22,p_z)-SH(sY11,sY2n2,p_z))*(os_down'*os_up);
M(p3212,dzy,s_up,p_z)= (-1i*sqrt(.5))*sqrt(2/3)*(SH(sY10,sY21,p_z)+SH(sY10,sY2n1,p_z))*(os_up'*os_up) + (-1i*sqrt(.5))*sqrt(1/3)*(SH(sY11,sY21,p_z)+SH(sY11,sY2n1,p_z))*(os_down'*os_up);
M(p3212,dxz,s_up,p_z)= sqrt(.5)*sqrt(2/3)*(SH(sY10,sY21,p_z)-SH(sY10,sY2n1,p_z))*(os_up'*os_up) + sqrt(.5)*sqrt(1/3)*(SH(sY11,sY21,p_z)-SH(sY11,sY2n1,p_z))*(os_down'*os_up);

M(p3212,dx2y2,s_down,p_z)= sqrt(.5)* sqrt(2/3)*(SH(sY10,sY22,p_z)+SH(sY10,sY2n2,p_z))*(os_up'*os_down) + sqrt(.5)* sqrt(1/3)*(SH(sY11,sY22,p_z)+SH(sY11,sY2n2,p_z))*(os_down'*os_down);
M(p3212,d3z2r2,s_down,p_z)= sqrt(2/3)*SH(sY10,sY20,p_z)*(os_up'*os_down) + sqrt(1/3)*SH(sY11,sY20,p_z)*(os_down'*os_down);
M(p3212,dxy,s_down,p_z)= (-1i*sqrt(.5))*sqrt(2/3)*(SH(sY10,sY22,p_z)-SH(sY10,sY2n2,p_z))*(os_up'*os_down) + (-1i*sqrt(.5))*sqrt(1/3)*(SH(sY11,sY22,p_z)-SH(sY11,sY2n2,p_z))*(os_down'*os_down);
M(p3212,dzy,s_down,p_z)= (-1i*sqrt(.5))*sqrt(2/3)*(SH(sY10,sY21,p_z)+SH(sY10,sY2n1,p_z))*(os_up'*os_down) + (-1i*sqrt(.5))*sqrt(1/3)*(SH(sY11,sY21,p_z)+SH(sY11,sY2n1,p_z))*(os_down'*os_down);
M(p3212,dxz,s_down,p_z)= sqrt(.5)*sqrt(2/3)*(SH(sY10,sY21,p_z)-SH(sY10,sY2n1,p_z))*(os_up'*os_down) + sqrt(.5)*sqrt(1/3)*(SH(sY11,sY21,p_z)-SH(sY11,sY2n1,p_z))*(os_down'*os_down);


disp("2p3212 3dx2y2_up");
disp (M(p3212,dx2y2,s_up,p_z)); 
disp("2p3212 3dx2y2_down");
disp (M(p3212,dx2y2,s_down,p_z)); 
disp("2p3212 3d3z2r2_up");
disp (M(p3212,d3z2r2,s_up,p_z)); 
disp("2p3212 3d3z2r2_down");
disp (M(p3212,d3z2r2,s_down,p_z)); 
disp("2p3212 3dxy_up");
disp (M(p3212,dxy,s_up,p_z)); 
disp("2p3212 3dxy_down");
disp (M(p3212,dxy,s_down,p_z));
disp("2p3212 3dzy_up");
disp (M(p3212,dzy,s_up,p_z));
disp("2p3212 3dzy_down");
disp (M(p3212,dzy,s_down,p_z));
disp("2p3212 3dxz_up");
disp (M(p3212,dxz,s_up,p_z)); 
disp("2p3212 3dxz_down");
disp (M(p3212,dxz,s_down,p_z)); 

disp(M)



M(p3232,dx2y2,s_up,p_y)=sum(sum(conj(op3232).*op_y.*o3dx2y2_up.*dtheta)*dtheta_sh*dphi_sh);
M(p3232,dx2y2,s_up,p_z)=sum(sum(conj(op3232).*op_z.*o3dx2y2_up.*dtheta)*dtheta_sh*dphi_sh);


























%{








 
%Now we will define the orbitals; 2p orbitals with strong LS coupling 
op32n32 = Y1n1;
op32n12 = sqrt(1/3)*Y1n1 + sqrt(2/3)*Y10;
op3212 = sqrt(2/3)*Y10 + sqrt(1/3)*Y11;
op3232 = Y11;
 
 
 
%Define 3d orbitals with the spin
 
o3dx2y2_down = sqrt(1/2)*(Y22 + Y2n2);
o3dx2y2_up = sqrt(.5)*(Y22 + Y2n2);
o3dz2r2_down = (Y20);
o3dz2r2_up =  (Y20);
o3dxy_down = (-1i/sqrt(2))*(Y22-Y2n2);
o3dxy_up = (-1i/sqrt(2))*(Y22-Y2n2);
o3dzy_down = (-1i/sqrt(2))*(Y21+Y2n1);
o3dzy_up = (-1i/sqrt(2))*(Y21+Y2n1);
o3dxz_down = sqrt(.5)*(Y21-Y2n1);
o3dxz_up = sqrt(.5)*(Y21-Y2n1);





%{
Calculating the matrix elements for linearly polarized light. 
There will (4*5*2 = 40 ) matrix elements.
We will calculate the matrix elements in terms of holes, where
hole transitions from 3d orbital to 2p orbital
%}

%Transition from 2p orbitals to 3dx2y2 up

disp("2p3232 to 3dx2y2 spin up");
s=sum((conj(spinor_up).*os_up));
M(p3232,dx2y2,s_up,p_x)=sum(sum(conj(op3232).*op_x.*o3dx2y2_up.*dtheta)*dtheta_sh*dphi_sh);
M(p3232,dx2y2,s_up,p_y)=sum(sum(conj(op3232).*op_y.*o3dx2y2_up.*dtheta)*dtheta_sh*dphi_sh);
M(p3232,dx2y2,s_up,p_z)=sum(sum(conj(op3232).*op_z.*o3dx2y2_up.*dtheta)*dtheta_sh*dphi_sh);

%M(p3232,dx2y2,s_up,p_y)=sum(sum(conj(op3232).*op_y.*o3dx2y2_up))*(conj(spinor_down).*os_up)*dtheta_sh*dphi_sh;
%M(p3232,dx2y2,s_up,p_z)=sum(sum(conj(op3232).*op_z.*o3dx2y2_up))*(conj(spinor_down).*os_up)*dtheta_sh*dphi_sh;
disp(M(p3232,dx2y2,s_up,p_x));
disp(M(p3232,dx2y2,s_up,p_y));
disp(M(p3232,dx2y2,s_up,p_z));

disp("2p32n32 to 3dx2y2 spin down");

M(p32n32,dx2y2,s_down,p_x)=sum(sum(conj(op32n32).*op_x.*o3dx2y2_down.*dtheta))*dtheta_sh*dphi_sh;
M(p32n32,dx2y2,s_down,p_y)=sum(sum(conj(op32n32).*op_y.*o3dx2y2_down.*dtheta))*dtheta_sh*dphi_sh;
M(p32n32,dx2y2,s_down,p_z)=sum(sum(conj(op32n32).*op_z.*o3dx2y2_down.*dtheta))*dtheta_sh*dphi_sh;

%M(p3232,dx2y2,s_up,p_y)=sum(sum(conj(op3232).*op_y.*o3dx2y2_up))*(conj(spinor_down).*os_up)*dtheta_sh*dphi_sh;
%M(p3232,dx2y2,s_up,p_z)=sum(sum(conj(op3232).*op_z.*o3dx2y2_up))*(conj(spinor_down).*os_up)*dtheta_sh*dphi_sh;
disp(M(p32n32,dx2y2,s_down,p_x));
disp(M(p32n32,dx2y2,s_down,p_y));
disp(M(p32n32,dx2y2,s_down,p_z));




disp("2p3232 to 3dz2r2 spin up");
s=sum((conj(spinor_up).*os_up));
M(p3232,d3z2r2,s_up,p_x)=sum(sum(conj(op3232).*op_x.*o3dz2r2_up.*dtheta)*dtheta_sh*dphi_sh);
M(p3232,d3z2r2,s_up,p_y)=sum(sum(conj(op3232).*op_y.*o3dz2r2_up.*dtheta)*dtheta_sh*dphi_sh);
M(p3232,d3z2r2,s_up,p_z)=sum(sum(conj(op3232).*op_z.*o3dz2r2_up.*dtheta)*dtheta_sh*dphi_sh);

%M(p3232,dx2y2,s_up,p_y)=sum(sum(conj(op3232).*op_y.*o3dx2y2_up))*(conj(spinor_down).*os_up)*dtheta_sh*dphi_sh;
%M(p3232,dx2y2,s_up,p_z)=sum(sum(conj(op3232).*op_z.*o3dx2y2_up))*(conj(spinor_down).*os_up)*dtheta_sh*dphi_sh;
disp(M(p3232,d3z2r2,s_up,p_x));
disp(M(p3232,d3z2r2,s_up,p_y));
disp(M(p3232,d3z2r2,s_up,p_z));

disp("2p32n32 to 3dz2r2 spin down");

M(p32n32,d3z2r2,s_down,p_x)=sum(sum(conj(op32n32).*op_x.*o3dz2r2_down.*dtheta))*dtheta_sh*dphi_sh;
M(p32n32,d3z2r2,s_down,p_y)=sum(sum(conj(op32n32).*op_y.*o3dz2r2_down.*dtheta))*dtheta_sh*dphi_sh;
M(p32n32,d3z2r2,s_down,p_z)=sum(sum(conj(op32n32).*op_z.*o3dz2r2_down.*dtheta))*dtheta_sh*dphi_sh;

%M(p3232,dx2y2,s_up,p_y)=sum(sum(conj(op3232).*op_y.*o3dx2y2_up))*(conj(spinor_down).*os_up)*dtheta_sh*dphi_sh;
%M(p3232,dx2y2,s_up,p_z)=sum(sum(conj(op3232).*op_z.*o3dx2y2_up))*(conj(spinor_down).*os_up)*dtheta_sh*dphi_sh;
disp(M(p32n32,d3z2r2,s_down,p_x));
disp(M(p32n32,d3z2r2,s_down,p_y));
disp(M(p32n32,d3z2r2,s_down,p_z));


disp("2p3232 to 3dxy spin up");
s=sum((conj(spinor_up).*os_up));
M(p3232,dxy,s_up,p_x)=sum(sum(conj(op3232).*op_x.*o3dxy_up.*dtheta)*dtheta_sh*dphi_sh);
M(p3232,dxy,s_up,p_y)=sum(sum(conj(op3232).*op_y.*o3dxy_up.*dtheta)*dtheta_sh*dphi_sh);
M(p3232,dxy,s_up,p_z)=sum(sum(conj(op3232).*op_z.*o3dxy_up.*dtheta)*dtheta_sh*dphi_sh);

%M(p3232,dx2y2,s_up,p_y)=sum(sum(conj(op3232).*op_y.*o3dx2y2_up))*(conj(spinor_down).*os_up)*dtheta_sh*dphi_sh;
%M(p3232,dx2y2,s_up,p_z)=sum(sum(conj(op3232).*op_z.*o3dx2y2_up))*(conj(spinor_down).*os_up)*dtheta_sh*dphi_sh;
disp(M(p3232,dxy,s_up,p_x));
disp(M(p3232,dxy,s_up,p_y));
disp(M(p3232,dxy,s_up,p_z));

disp("2p32n32 to 3dxy spin down");

M(p32n32,dxy,s_down,p_x)=sum(sum(conj(op32n32).*op_x.*o3dxy_down.*dtheta))*dtheta_sh*dphi_sh;
M(p32n32,dxy,s_down,p_y)=sum(sum(conj(op32n32).*op_y.*o3dxy_down.*dtheta))*dtheta_sh*dphi_sh;
M(p32n32,dxy,s_down,p_z)=sum(sum(conj(op32n32).*op_z.*o3dxy_down.*dtheta))*dtheta_sh*dphi_sh;

%M(p3232,dx2y2,s_up,p_y)=sum(sum(conj(op3232).*op_y.*o3dx2y2_up))*(conj(spinor_down).*os_up)*dtheta_sh*dphi_sh;
%M(p3232,dx2y2,s_up,p_z)=sum(sum(conj(op3232).*op_z.*o3dx2y2_up))*(conj(spinor_down).*os_up)*dtheta_sh*dphi_sh;
disp(M(p32n32,dxy,s_down,p_x));
disp(M(p32n32,dxy,s_down,p_y));
disp(M(p32n32,dxy,s_down,p_z));

disp("2p3232 to 3dzy spin up");
s=sum((conj(spinor_up).*os_up));
M(p3232,dzy,s_up,p_x)=sum(sum(conj(op3232).*op_x.*o3dzy_up.*dtheta)*dtheta_sh*dphi_sh);
M(p3232,dzy,s_up,p_y)=sum(sum(conj(op3232).*op_y.*o3dzy_up.*dtheta)*dtheta_sh*dphi_sh);
M(p3232,dzy,s_up,p_z)=sum(sum(conj(op3232).*op_z.*o3dzy_up.*dtheta)*dtheta_sh*dphi_sh);

%M(p3232,dx2y2,s_up,p_y)=sum(sum(conj(op3232).*op_y.*o3dx2y2_up))*(conj(spinor_down).*os_up)*dtheta_sh*dphi_sh;
%M(p3232,dx2y2,s_up,p_z)=sum(sum(conj(op3232).*op_z.*o3dx2y2_up))*(conj(spinor_down).*os_up)*dtheta_sh*dphi_sh;
disp(M(p3232,dzy,s_up,p_x));
disp(M(p3232,dzy,s_up,p_y));
disp(M(p3232,dzy,s_up,p_z));

disp("2p32n32 to 3dzy spin down");

M(p32n32,dzy,s_down,p_x)=sum(sum(conj(op32n32).*op_x.*o3dzy_down.*dtheta))*dtheta_sh*dphi_sh;
M(p32n32,dzy,s_down,p_y)=sum(sum(conj(op32n32).*op_y.*o3dzy_down.*dtheta))*dtheta_sh*dphi_sh;
M(p32n32,dzy,s_down,p_z)=sum(sum(conj(op32n32).*op_z.*o3dzy_down.*dtheta))*dtheta_sh*dphi_sh;

%M(p3232,dx2y2,s_up,p_y)=sum(sum(conj(op3232).*op_y.*o3dx2y2_up))*(conj(spinor_down).*os_up)*dtheta_sh*dphi_sh;
%M(p3232,dx2y2,s_up,p_z)=sum(sum(conj(op3232).*op_z.*o3dx2y2_up))*(conj(spinor_down).*os_up)*dtheta_sh*dphi_sh;
disp(M(p32n32,dzy,s_down,p_x));
disp(M(p32n32,dzy,s_down,p_y));
disp(M(p32n32,dzy,s_down,p_z));

disp("2p3232 to 3dxz spin up");
s=sum((conj(spinor_up).*os_up));
M(p3232,dxz,s_up,p_x)=sum(sum(conj(op3232).*op_x.*o3dxz_up.*dtheta)*dtheta_sh*dphi_sh);
M(p3232,dxz,s_up,p_y)=sum(sum(conj(op3232).*op_y.*o3dxz_up.*dtheta)*dtheta_sh*dphi_sh);
M(p3232,dxz,s_up,p_z)=sum(sum(conj(op3232).*op_z.*o3dxz_up.*dtheta)*dtheta_sh*dphi_sh);

%M(p3232,dx2y2,s_up,p_y)=sum(sum(conj(op3232).*op_y.*o3dx2y2_up))*(conj(spinor_down).*os_up)*dtheta_sh*dphi_sh;
%M(p3232,dx2y2,s_up,p_z)=sum(sum(conj(op3232).*op_z.*o3dx2y2_up))*(conj(spinor_down).*os_up)*dtheta_sh*dphi_sh;
disp(M(p3232,dxz,s_up,p_x));
disp(M(p3232,dxz,s_up,p_y));
disp(M(p3232,dxz,s_up,p_z));

disp("2p32n32 to 3dxy spin down");

M(p32n32,dxz,s_down,p_x)=sum(sum(conj(op32n32).*op_x.*o3dxz_down.*dtheta))*dtheta_sh*dphi_sh;
M(p32n32,dxz,s_down,p_y)=sum(sum(conj(op32n32).*op_y.*o3dxz_down.*dtheta))*dtheta_sh*dphi_sh;
M(p32n32,dxz,s_down,p_z)=sum(sum(conj(op32n32).*op_z.*o3dxz_down.*dtheta))*dtheta_sh*dphi_sh;

%M(p3232,dx2y2,s_up,p_y)=sum(sum(conj(op3232).*op_y.*o3dx2y2_up))*(conj(spinor_down).*os_up)*dtheta_sh*dphi_sh;
%M(p3232,dx2y2,s_up,p_z)=sum(sum(conj(op3232).*op_z.*o3dx2y2_up))*(conj(spinor_down).*os_up)*dtheta_sh*dphi_sh;
disp(M(p32n32,dxz,s_down,p_x));
disp(M(p32n32,dxz,s_down,p_y));
disp(M(p32n32,dxz,s_down,p_z));

%{
disp("p3212 to 3dx2y2 spin up");
M(p3212,dx2y2,s_up,p_x)=sum(sum(conj(op3212).*op_x.*o3dx2y2_up))**dtheta_sh*dphi_sh;
M(p3212,dx2y2,s_up,p_y)=sum(sum(conj(op3212).*op_y.*o3dx2y2_up))*dtheta_sh*dphi_sh;
M(p3212,dx2y2,s_up,p_z)=sum(sum(conj(op3212).*op_z.*o3dx2y2_up))*dtheta_sh*dphi_sh;
disp(M(p3212,dx2y2,s_up,p_x));
disp(M(p3212,dx2y2,s_up,p_y));
disp(M(p3212,dx2y2,s_up,p_z));

disp("p32n12 to 3dx2y2 spin up");
M(p32n12,dx2y2,s_up,p_x)=sum(sum(conj(op32n12).*op_x.*o3dx2y2_up))*dtheta_sh*dphi_sh;
M(p32n12,dx2y2,s_up,p_y)=sum(sum(conj(op32n12).*op_y.*o3dx2y2_up))*dtheta_sh*dphi_sh;
M(p32n12,dx2y2,s_up,p_z)=sum(sum(conj(op32n12).*op_z.*o3dx2y2_up))*dtheta_sh*dphi_sh;
disp(M(p32n12,dx2y2,s_up,p_x));
disp(M(p32n12,dx2y2,s_up,p_y));
disp(M(p32n12,dx2y2,s_up,p_z));

disp("p32n32 to 3dx2y2 spin up");
M(p32n32,dx2y2,s_up,p_x)=sum(sum(conj(op32n32).*op_x.*o3dx2y2_up))*dtheta_sh*dphi_sh;
M(p32n32,dx2y2,s_up,p_y)=sum(sum(conj(op32n32).*op_y.*o3dx2y2_up))*dtheta_sh*dphi_sh;
M(p32n32,dx2y2,s_up,p_z)=sum(sum(conj(op32n32).*op_z.*o3dx2y2_up))*dtheta_sh*dphi_sh;
disp(M(p32n32,dx2y2,s_up,p_x));
disp(M(p32n32,dx2y2,s_up,p_y));
disp(M(p32n32,dx2y2,s_up,p_z));

% For pi polarization
M(p3232,dx2y2,s_up,p_pi)=sum(sum(conj(op3232).*op_pi.*o3dx2y2_up))*dtheta_sh*dphi_sh;
M(p3212,dx2y2,s_up,p_pi)=sum(sum(conj(op3212).*op_pi.*o3dx2y2_up))*dtheta_sh*dphi_sh;
M(p32n12,dx2y2,s_up,p_pi)=sum(sum(conj(op32n12).*op_pi.*o3dx2y2_up))*dtheta_sh*dphi_sh;
M(p32n32,dx2y2,s_up,p_pi)=sum(sum(conj(op32n32).*op_pi.*o3dx2y2_up))*dtheta_sh*dphi_sh;







%Transition from 2p orbitals to 3dx2y2 down

disp("3p3232 to 3dx2y2 spin down");
M(p3232,dx2y2,s_down,p_x)=sum(sum(conj(op3232).*op_x.*o3dx2y2_down))*dtheta_sh*dphi_sh;
M(p3232,dx2y2,s_down,p_y)=sum(sum(conj(op3232).*op_y.*o3dx2y2_down))*dtheta_sh*dphi_sh;
M(p3232,dx2y2,s_down,p_z)=sum(sum(conj(op3232).*op_z.*o3dx2y2_down))*dtheta_sh*dphi_sh;
disp(M(p3232,dx2y2,s_down,p_x));
disp(M(p3232,dx2y2,s_down,p_y));
disp(M(p3232,dx2y2,s_down,p_z));

disp("p3212 to 3dx2y2 spin down");
M(p3212,dx2y2,s_down,p_x)=sum(sum(conj(op3212).*op_x.*o3dx2y2_down))*dtheta_sh*dphi_sh;
M(p3212,dx2y2,s_down,p_y)=sum(sum(conj(op3212).*op_y.*o3dx2y2_down))*dtheta_sh*dphi_sh;
M(p3212,dx2y2,s_down,p_z)=sum(sum(conj(op3212).*op_z.*o3dx2y2_down))*dtheta_sh*dphi_sh;
disp(M(p3212,dx2y2,s_down,p_x));
disp(M(p3212,dx2y2,s_down,p_y));
disp(M(p3212,dx2y2,s_down,p_z));

disp("p32n12 to 3dx2y2 spin down");
M(p32n12,dx2y2,s_down,p_x)=sum(sum(conj(op32n12).*op_x.*o3dx2y2_down))*dtheta_sh*dphi_sh;
M(p32n12,dx2y2,s_down,p_y)=sum(sum(conj(op32n12).*op_y.*o3dx2y2_down))*dtheta_sh*dphi_sh;
M(p32n12,dx2y2,s_down,p_z)=sum(sum(conj(op32n12).*op_z.*o3dx2y2_down))*dtheta_sh*dphi_sh;
disp(M(p32n12,dx2y2,s_down,p_x));
disp(M(p32n12,dx2y2,s_down,p_y));
disp(M(p32n12,dx2y2,s_down,p_z));

disp("p32n32 to 3dx2y2 spin down");
M(p32n32,dx2y2,s_down,p_x)=sum(sum(conj(op32n32).*op_x.*o3dx2y2_down))*dtheta_sh*dphi_sh;
M(p32n32,dx2y2,s_down,p_y)=sum(sum(conj(op32n32).*op_y.*o3dx2y2_down))*dtheta_sh*dphi_sh;
M(p32n32,dx2y2,s_down,p_z)=sum(sum(conj(op32n32).*op_z.*o3dx2y2_down))*dtheta_sh*dphi_sh;
disp(M(p32n32,dx2y2,s_down,p_x));
disp(M(p32n32,dx2y2,s_down,p_y));
disp(M(p32n32,dx2y2,s_down,p_z));


%}









%{

%Alternatively you can use the "for loop" implementation if needed

m11=0;


for theta_sh=0:dtheta_sh:pi
    for phi_sh=0:dphi_sh:2*pi
        Y22=(.25)*sqrt(15/(2*pi))*exp(2i*phi_sh)*(sin(theta_sh)*sin(theta_sh));
        Y21=(-.5)*sqrt(15/(2*pi))*exp(1i*phi_sh)*(sin(theta_sh)*cos(theta_sh));
        Y20=(.25)*sqrt(5/(pi))*(3*cos(theta_sh)*cos(theta_sh)-1);
        Y2n1=(.5)*sqrt(15/(2*pi))*exp(-1i*phi_sh)*(sin(theta_sh)*cos(theta_sh));
        Y2n2=(.25)*sqrt(15/(2*pi))*exp(-2i*phi_sh)*(sin(theta_sh)*sin(theta_sh));
        Y11=(-.5)*sqrt(3/(2*pi))*sin(theta_sh)*exp(1i*phi_sh) ;        
        Y10= (.5)*sqrt(3/pi)*cos(theta_sh);
        Y1n1=(.5)*sqrt(3/(2*pi))*exp(-1i*phi_sh)*sin(theta_sh);
        
       
    
        o3dx2y2_up = os_up*sqrt(.5)*(Y22 + Y2n2);
    

        %Definie 2p orbitals in the J Jz basis
        op32n32 = Y1n1;
        op32n12 = sqrt(1/3)*Y1n1 + sqrt(2/3)*Y10;
        op3212 = sqrt(2/3)*Y10 + sqrt(1/3)*Y11;
        op3232 = Y11;

        
        op_x=(1/sqrt(2))*sqrt(4*pi/3)*(Y1n1-Y11);
        
        
        m11 = m11 + conj(op3232).*op_x.*o3dx2y2_up*sin(theta_sh)*dtheta_sh*dphi_sh;
    end
end

disp(m11);

%}

 















%}



