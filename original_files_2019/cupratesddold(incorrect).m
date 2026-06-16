clc;

%Index for 3d orbitals
d3z2r2=1;
dx2y2=2;
dxy=3;
dyz=4;
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
p_x=1;
p_y=2;
p_z=3;


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

% All the matrix elements not defined here are zero;
M(p3232,dx2y2,s_up,p_x)=(-sqrt(6)/2);
M(p3232,dx2y2,s_up,p_y)=(-1i*sqrt(6)/2);

M(p32n12,dx2y2,s_up,p_x)=sqrt(.5);
M(p32n12,dx2y2,s_up,p_y)=(-1i*sqrt(.5));

M(p3212,dx2y2,s_down,p_x)=-sqrt(.5);
M(p3212,dx2y2,s_down,p_y)=(-1i*sqrt(.5));

M(p32n32,dx2y2,s_down,p_x)=(sqrt(6)/2);
M(p32n32,dx2y2,s_down,p_y)=(-1i*sqrt(6)/2);


M(p3232,d3z2r2,s_up,p_x)=sqrt(.5);
M(p3232,d3z2r2,s_up,p_y)=(-1i*sqrt(.5));

M(p3212,d3z2r2,s_up,p_z)=4/sqrt(6);

M(p32n12,d3z2r2,s_up,p_x)=-1/sqrt(6);
M(p32n12,d3z2r2,s_up,p_y)=(-1i/sqrt(6));

M(p3212,d3z2r2,s_down,p_x)=1/sqrt(6);
M(p3212,d3z2r2,s_down,p_y)=(-1i/sqrt(6));

M(p32n12,d3z2r2,s_down,p_z)=4/sqrt(6);

M(p32n32,d3z2r2,s_down,p_x)=-sqrt(.5);
M(p32n32,d3z2r2,s_down,p_y)=(-1i*sqrt(.5));


M(p3232,dxy,s_up,p_x)=(1i*sqrt(6)/2);
M(p3232,dxy,s_up,p_y)=-sqrt(6)/2;

M(p32n12,dxy,s_up,p_x)=1i*sqrt(.5);
M(p32n12,dxy,s_up,p_y)=sqrt(.5);

M(p3212,dxy,s_down,p_x)=1i*sqrt(.5);
M(p3212,dxy,s_down,p_y)=-sqrt(.5);

M(p32n32,dxy,s_down,p_x)=1i*sqrt(6)/2;
M(p32n32,dxy,s_down,p_y)=sqrt(6)/2;


M(p3232,dyz,s_up,p_z)=-1i*sqrt(6)/2;

M(p3212,dyz,s_up,p_y)=-sqrt(2);

M(p32n12,dyz,s_up,p_z)=-1i*sqrt(.5);

M(p3212,dyz,s_down,p_z)=-1i*sqrt(.5);

M(p32n12,dyz,s_down,p_y)=-sqrt(2);

M(p32n32,dyz,s_down,p_z)=-1i*sqrt(6)/2;


M(p3232,dxz,s_up,p_z)=sqrt(6)/2;

M(p3212,dxz,s_up,p_x)=-sqrt(2);

M(p32n12,dxz,s_up,p_z)=-sqrt(.5);

M(p3212,dxz,s_down,p_z)=sqrt(.5);

M(p32n12,dxz,s_down,p_x)=-sqrt(2);

M(p32n32,dxz,s_down,p_z)=-sqrt(6)/2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
Incident and reflected angles 
Theta  gives information on the momentum that is passed into the 2d substrate. 
The phi angles gives information on crystal axis along which it is passed.
%}
theta_i=0;
phi=0;          % incoming and outgoing phi angles are same
theta_o=0;


%{
Experimental angles
delta - refers to the tilting angle of the substrate
theta - Half angle of the incoming and outgoing orientation. Either 45 ( 2theta=90 ) or 65 for 2theta=130;
%}
delta=-pi/4:.01:pi/4;
theta=pi/4;  

% The crossection dependence is calculated with respect to q_par which is
% the momentum transferred parallel into the surface

k=.5;   % Place holder for incoming photon momentum since we plol arbitraty units of intensity we are free to choose its value.

q_par=2*k*sin(theta)*sin(delta);

%Plot the q-par dependence w.r.t to delta

%plot(delta,q_par);


%{ 
So far we have defined all of the necessary paramters associated with the
experiment. However our cross-section values are all w.r.t x, y z
polarizations alone. With delta and phi being our probing parameters we
will from now on we will have sigma and pi polarizations. These can then be
expressed in the x-y-z polarization basis with appropriate prefactors
dependent on probing parameters there will be seperate parameteris for
incoming and outgoing. However pi polarization parameters are exactly the
same
%}

p_pi_x = cos(phi-pi/2);
p_pi_y = sin(phi-pi/2);
p_pi_z = 0;

p_sigma_x_in = -sin(phi-pi/2)*cos(pi/2-theta);
p_sigma_y_in = cos(delta).*cos(phi-pi/2).*cos(pi/2-theta) - sin(delta).*sin(pi/2-theta);
p_sigma_z_in = sin(delta).*cos(phi-pi/2).*cos(pi/2-theta) + cos(delta).*sin(pi/2-theta);

p_sigma_x_out = sin(phi-pi/2)*cos(pi/2-theta);
p_sigma_y_out = -cos(delta).*cos(phi-pi/2).*cos(pi/2-theta) - sin(delta).*sin(pi/2-theta);
p_sigma_z_out = -sin(delta).*cos(phi-pi/2).*cos(pi/2-theta) + cos(delta).*sin(pi/2-theta);

% We will start calculation for pi-polarizad light that excites a spin-up
% electron.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The ground state of the Cuprate ion is dominated by the dx2y2 orbit. The
% ground state spin can have arbitrary orientations represented using the 
% Spin Orientation angles .

theta_s = pi/2;
phi_s =pi/4;

%The spin part is represented in the Sz basis with |+> and |->  as the basis states. 
%Spinors which are the basis vectors |+> and |-> aligned along z axis are
%defined below:
spinor_up = [1;0];
spinor_down = [0;1];

%Ground State Spins are then given as follows:

os_up=  cos(theta_s/2)*exp(-1i*phi_s).*spinor_up + sin(theta_s/2).*spinor_down;
os_down = sin(theta_s/2)*exp(-1i*phi_s).*spinor_up - cos(theta_s/2).*spinor_down;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{ 
Our transitions occur in a 2 step process. In hole language:
Step 1
    3dx2y2 hole is de-excited into any of the 2p orbitals.
Step 2
    2p hole is then excited to any of the other 3d orbitals

The orbitals do have corresponding spins.
The allowed channels for each of these transitions correspond to the
non-zero matrix elements of the Dipole operator.
%}




% Elastic line - Ground state to dx2y2|spin_down>. There are two possible
% transitions for this channel thru 2p3232 and 2p32n12. Incoming light could
% be sigma or pi polarized. The outgoing light polarization is not measured
% so it is summed over both sigma and pi polarization


% Incoming light pi polarized \\Same as outgoing light polarization

dx2y2up3232 = M(p3232,dx2y2,s_up,p_x)*p_pi_x + M(p3232,dx2y2,s_up,p_y)*p_pi_y;  % We know that p_pi_z=0
dx2y2up32n12= M(p32n12,dx2y2,s_up,p_x)*p_pi_x + M(p32n12,dx2y2,s_up,p_y)*p_pi_y;  % We know that p_pi_z=0

% Incoming light sigma polarized

dx2y2up3232_sigma = M(p3232,dx2y2,s_up,p_x)*p_sigma_x_in + M(p3232,dx2y2,s_up,p_y)*p_sigma_y_in + M(p3232,dx2y2,s_up,p_z)*p_sigma_z_in;
dx2y2up32n12_sigma = M(p32n12,dx2y2,s_up,p_x)*p_sigma_x_in + M(p32n12,dx2y2,s_up,p_y)*p_sigma_y_in + M(p32n12,dx2y2,s_up,p_z)*p_sigma_z_in;


% Outgoing light sigma polarized

dx2y2up3232_sigma_out = M(p3232,dx2y2,s_up,p_x)*p_sigma_x_out + M(p3232,dx2y2,s_up,p_y)*p_sigma_y_out + M(p3232,dx2y2,s_up,p_z)*p_sigma_z_out;
dx2y2up32n12_sigma_out = M(p32n12,dx2y2,s_up,p_x)*p_sigma_x_out + M(p32n12,dx2y2,s_up,p_y)*p_sigma_y_out + M(p32n12,dx2y2,s_up,p_z)*p_sigma_z_out;



C1 = conj(conj(dx2y2up3232_sigma_out).*dx2y2up3232).*(conj(dx2y2up3232_sigma_out).*dx2y2up3232)*(conj(os_up'*spinor_up).*(os_up'*spinor_up));
C2 = conj(conj(dx2y2up3232).*dx2y2up3232).*(conj(dx2y2up3232).*dx2y2up3232)*(conj(os_up'*spinor_up).*(os_up'*spinor_up));
C3 = conj(conj(dx2y2up32n12_sigma_out).*dx2y2up3232).*(conj(dx2y2up32n12_sigma_out).*dx2y2up3232)*(conj(os_up'*spinor_up).*(os_up'*spinor_up));
C4 = conj(conj(dx2y2up32n12).*dx2y2up3232).*(conj(dx2y2up32n12).*dx2y2up3232)*(conj(os_up'*spinor_up).*(os_up'*spinor_up));

cs = C1+C2+C3+C4;

plot(q_par,cs,'g');
































%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
 we assume that the dx2y2_down is the ground state of the hole from which
% all transitions occur.

% csd3z2r2_up = dx2y2-> 2p orbital then 2p->d3z2r2_up


%For incoming pi-polarized light these are the orbital parts
dx2y2down3212 = M(p3212,dx2y2,s_down,p_x)*p_pi_x + M(p3212,dx2y2,s_down,p_y)*p_pi_y;  % We know that p_pi_z=0

d3z2r2upp3212_pi_out = M(p3212,d3z2r2,s_up,p_x)*p_pi_x + M(p3212,d3z2r2,s_up,p_y)*p_pi_y;

d3z2r2upp3212_sigma_out = M(p3212,d3z2r2,s_up,p_x)*p_sigma_x_out+ M(p3212,d3z2r2,s_up,p_y)*p_sigma_y_out +M(p3212,d3z2r2,s_up,p_z)*p_sigma_z_out;

% Here we have to account for the spin parts as well. The ground state is
% in a mixed state of up and down spinors while the final state is either
% in spin up or down along the quantization axis. 

p=conj((conj(d3z2r2upp3212_pi_out).*dx2y2down3212)).*conj(d3z2r2upp3212_pi_out).*dx2y2down3212;
disp(p);

csd3zr2up_pi =(conj((conj(d3z2r2upp3212_sigma_out).*dx2y2down3212)).*conj(d3z2r2upp3212_sigma_out).*dx2y2down3212 + conj((conj(d3z2r2upp3212_pi_out).*dx2y2down3212)).*conj(d3z2r2upp3212_pi_out).*dx2y2down3212) *(conj(os_down'*spinor_up)*(os_down'*spinor_up))  ;

plot(q_par,csd3zr2up_pi,'g'); hold on;

%For incoming sigma polarized light
dx2y2down3212_sigma = M(p3212,dx2y2,s_down,p_x)*p_sigma_x_in + M(p3212,dx2y2,s_down,p_y)*p_sigma_y_in + M(p3212,dx2y2,s_down,p_z)*p_sigma_z_in;

csd3z2r2up_sigma = (conj(d3z2r2upp3212_sigma_out).*dx2y2down3212_sigma).*conj(conj(d3z2r2upp3212_sigma_out).*dx2y2down3212_sigma) + (conj(d3z2r2upp3212_pi_out).*dx2y2down3212_sigma).*conj(conj(d3z2r2upp3212_pi_out).*dx2y2down3212_sigma);

plot(q_par,csd3z2r2up_sigma,'b');




%{
%From p3232 the following transitions are possible. (Spin also accounted for)

%p3232

dx2y2upp3232_pi = M(p3232,dx2y2,s_up,p_x)*p_pi_x + M(p3232,dx2y2,s_up,p_y)*p_pi_y;
dx2y2upp3232_sigma = M(p3232,dx2y2,s_up,p_x)*p_sigma_x_out+ M(p3232,dx2y2,s_up,p_y)*p_sigma_y_out +M(p3232,dx2y2,s_up,p_z)*p_sigma_z_out;
d3z2r2upp3232_pi = M(p3232,d3z2r2,s_up,p_x)*p_pi_x + M(p3232,d3z2r2,s_up,p_y)*p_pi_y;
d3z2r2upp3232_sigma = M(p3232,d3z2r2,s_up,p_x)*p_sigma_x_out+ M(p3232,d3z2r2,s_up,p_y)*p_sigma_y_out +M(p3232,d3z2r2,s_up,p_z)*p_sigma_z_out;
dxyupp3232_pi = M(p3232,dxy,s_up,p_x)*p_pi_x + M(p3232,dxy,s_up,p_y)*p_pi_y;
dxyupp3232_sigma = M(p3232,dxy,s_up,p_x)*p_sigma_x_out+ M(p3232,dxy,s_up,p_y)*p_sigma_y_out +M(p3232,dxy,s_up,p_z)*p_sigma_z_out;
dyzupp3232_pi = M(p3232,dyz,s_up,p_x)*p_pi_x + M(p3232,dyz,s_up,p_y)*p_pi_y;
dyzupp3232_sigma = M(p3232,dyz,s_up,p_x)*p_sigma_x_out+ M(p3232,dyz,s_up,p_y)*p_sigma_y_out +M(p3232,dyz,s_up,p_z)*p_sigma_z_out;
dxzupp3232_pi = M(p3232,dxz,s_up,p_x)*p_pi_x + M(p3232,dxz,s_up,p_y)*p_pi_y;
dxzupp3232_sigma = M(p3232,dxz,s_up,p_x)*p_sigma_x_out+ M(p3232,dxz,s_up,p_y)*p_sigma_y_out +M(p3232,dxz,s_up,p_z)*p_sigma_z_out;

% The second p oribtal - p3212
dx2y2up3212 = M(p3212,dx2y2,s_up,p_x)*p_pi_x + M(p3212,dx2y2,s_up,p_y)*p_pi_y;

dx2y2downp3212_pi = M(p3212,dx2y2,s_down,p_x)*p_pi_x + M(p3212,dx2y2,s_down,p_y)*p_pi_y;
dx2y2downp3212_sigma = M(p3212,dx2y2,s_down,p_x)*p_sigma_x_out+ M(p3212,dx2y2,s_down,p_y)*p_sigma_y_out +M(p3212,dx2y2,s_down,p_z)*p_sigma_z_out;
d3z2r2upp3212_pi = M(p3212,d3z2r2,s_up,p_x)*p_pi_x + M(p3212,d3z2r2,s_up,p_y)*p_pi_y;
d3z2r2upp3212_sigma = M(p3212,d3z2r2,s_up,p_x)*p_sigma_x_out+ M(p3212,d3z2r2,s_up,p_y)*p_sigma_y_out +M(p3212,d3z2r2,s_up,p_z)*p_sigma_z_out;
d3z2r2downp3212_pi = M(p3212,d3z2r2,s_down,p_x)*p_pi_x + M(p3212,d3z2r2,s_down,p_y)*p_pi_y;
d3z2r2downp3212_sigma = M(p3212,d3z2r2,s_down,p_x)*p_sigma_x_out+ M(p3212,d3z2r2,s_down,p_y)*p_sigma_y_out +M(p3212,d3z2r2,s_down,p_z)*p_sigma_z_out;
dyzupp3212_pi = M(p3212,dyz,s_up,p_x)*p_pi_x + M(p3212,dyz,s_up,p_y)*p_pi_y;
dyzupp3212_sigma = M(p3212,dyz,s_up,p_x)*p_sigma_x_out+ M(p3212,dyz,s_up,p_y)*p_sigma_y_out +M(p3212,dyz,s_up,p_z)*p_sigma_z_out;
dyzdownp3212_pi = M(p3212,dyz,s_down,p_x)*p_pi_x + M(p3212,dyz,s_down,p_y)*p_pi_y;
dyzdownp3212_sigma = M(p3212,dyz,s_down,p_x)*p_sigma_x_out+ M(p3212,dyz,s_down,p_y)*p_sigma_y_out +M(p3212,dyz,s_down,p_z)*p_sigma_z_out;
dxzupp3212_pi = M(p3212,dxz,s_up,p_x)*p_pi_x + M(p3212,dxz,s_up,p_y)*p_pi_y;
dxzupp3212_sigma = M(p3212,dxz,s_up,p_x)*p_sigma_x_out+ M(p3212,dxz,s_up,p_y)*p_sigma_y_out +M(p3212,dxz,s_up,p_z)*p_sigma_z_out;
dxzdownp3212_pi = M(p3212,dxz,s_down,p_x)*p_pi_x + M(p3212,dxz,s_down,p_y)*p_pi_y;
dxzdownp3212_sigma = M(p3212,dxz,s_down,p_x)*p_sigma_x_out+ M(p3212,dxz,s_down,p_y)*p_sigma_y_out +M(p3212,dxz,s_down,p_z)*p_sigma_z_out;

%The third p-orbital - p32n12
dx2y2up32n12 = M(p32n12,dx2y2,s_up,p_x)*p_pi_x + M(p32n12,dx2y2,s_up,p_y)*p_pi_y;

dx2y2upp32n12_pi = M(p32n12,dx2y2,s_up,p_x)*p_pi_x + M(p32n12,dx2y2,s_up,p_y)*p_pi_y;
dx2y2upp32n12_sigma = M(p32n12,dx2y2,s_up,p_x)*p_sigma_x_out+ M(p32n12,dx2y2,s_up,p_y)*p_sigma_y_out +M(p32n12,dx2y2,s_up,p_z)*p_sigma_z_out;
d3z2r2upp32n12_pi = M(p32n12,d3z2r2,s_up,p_x)*p_pi_x + M(p32n12,d3z2r2,s_up,p_y)*p_pi_y;
d3z2r2upp32n12_sigma = M(p32n12,d3z2r2,s_up,p_x)*p_sigma_x_out+ M(p32n12,d3z2r2,s_up,p_y)*p_sigma_y_out +M(p32n12,d3z2r2,s_up,p_z)*p_sigma_z_out;
d3z2r2downp32n12_pi = M(p32n12,d3z2r2,s_down,p_x)*p_pi_x + M(p32n12,d3z2r2,s_down,p_y)*p_pi_y;
d3z2r2downp32n12_sigma = M(p32n12,d3z2r2,s_down,p_x)*p_sigma_x_out+ M(p32n12,d3z2r2,s_down,p_y)*p_sigma_y_out +M(p32n12,d3z2r2,s_down,p_z)*p_sigma_z_out;
dxyupp32n12_pi = M(p32n12,dxy,s_up,p_x)*p_pi_x + M(p32n12,dxy,s_up,p_y)*p_pi_y;
dxyupp32n12_sigma = M(p32n12,dxy,s_up,p_x)*p_sigma_x_out+ M(p32n12,dxy,s_up,p_y)*p_sigma_y_out +M(p32n12,dxy,s_up,p_z)*p_sigma_z_out;
dyzupp32n12_pi = M(p32n12,dyz,s_up,p_x)*p_pi_x + M(p32n12,dyz,s_up,p_y)*p_pi_y;
dyzupp32n12_sigma = M(p32n12,dyz,s_up,p_x)*p_sigma_x_out+ M(p32n12,dyz,s_up,p_y)*p_sigma_y_out +M(p32n12,dyz,s_up,p_z)*p_sigma_z_out;
dyzdownp32n12_pi = M(p32n12,dyz,s_down,p_x)*p_pi_x + M(p32n12,dyz,s_down,p_y)*p_pi_y;
dyzdownp32n12_sigma = M(p32n12,dyz,s_down,p_x)*p_sigma_x_out+ M(p32n12,dyz,s_down,p_y)*p_sigma_y_out +M(p32n12,dyz,s_down,p_z)*p_sigma_z_out;
dxzupp32n12_pi = M(p32n12,dxz,s_up,p_x)*p_pi_x + M(p32n12,dxz,s_up,p_y)*p_pi_y;
dxzupp32n12_sigma = M(p32n12,dxz,s_up,p_x)*p_sigma_x_out+ M(p32n12,dxz,s_up,p_y)*p_sigma_y_out +M(p32n12,dxz,s_up,p_z)*p_sigma_z_out;
dxzdownp32n12_pi = M(p32n12,dxz,s_down,p_x)*p_pi_x + M(p32n12,dxz,s_down,p_y)*p_pi_y;
dxzdownp32n12_sigma = M(p32n12,dxz,s_down,p_x)*p_sigma_x_out+ M(p32n12,dxz,s_down,p_y)*p_sigma_y_out +M(p32n12,dxz,s_down,p_z)*p_sigma_z_out;

dyz_sum_p_pi = dyzup3232*(dx2y2upp3232_pi+d3z2r2upp3232_pi+dxyupp3232_pi+dyzupp3232_pi+dxzupp3232_pi) + dyzupp3212*(dx2y2downp3212_pi+d3z2r2upp3212_pi+d3z2r2downp3212_pi+dyzupp3212_pi +dyzdownp3212_pi+dxzupp3212_pi+dxzdownp3212_pi) + dyzupp32n12*(dx2y2upp32n12_pi + d3z2r2upp32n12_pi + d3z2r2downp32n12_pi + dxyupp32n12_pi + dyzupp32n12_pi + dyzdownp32n12_pi + dxzupp32n12_pi + dxzdownp32n12_pi); 
dyz_sum_p_sigma = dyzup3232*(dx2y2upp3232_sigma+d3z2r2upp3232_sigma+dxyupp3232_sigma+dyzupp3232_sigma+dxzupp3232_sigma) + dyzupp3212*(dx2y2downp3212_sigma+d3z2r2upp3212_sigma+d3z2r2downp3212_sigma+dyzupp3212_sigma +dyzdownp3212_sigma+dxzupp3212_sigma+dxzdownp3212_sigma) + dyzupp32n12*(dx2y2upp32n12_sigma + d3z2r2upp32n12_sigma + d3z2r2downp32n12_sigma + dxyupp32n12_sigma + dyzupp32n12_sigma + dyzdownp32n12_sigma + dxzupp32n12_sigma + dxzdownp32n12_sigma); 

csdd = conj(dyz_sum_p_pi).*dyz_sum_p_pi + conj(dyz_sum_p_sigma).*dyz_sum_p_sigma;

%dyz pi spin_up
plot(delta,csdd);

for p2=1:1:4
    for d3=1:1:5
        for s=1:1:2
            csdd = (M(p2,dxy,s_up,p_x)*p_pi_x + M(p2,dxy,s_up,p_y)*p_pi_y) * (M(p2,d3,s,p_x)*p_pi_x + M(p2,d3,s,p_x)*p_pi_x
     

%}















%}
