    %{
Code Updates
12/2/19 -   Corrected the spin representations in the Sz basis
            Prefix of o represents the operator
12/3/19      _p_ represents pi polarization
            _s_ represents sigma polarization
12/4/19     rechecked dipole matrix elements
            Plots for dxy and dx2y2 are correct. Not the others
            All plots corrected.
            Note dxy/dyz implies dxy+dyz responses. 
12/5/19     Used Umesh's code for formatting plots


%}
clc;
clear;

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

M =  zeros(4,5,2,3);

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

M(p32n12,d3z2r2,s_up,p_x)= -1/sqrt(6);
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
Experimental angles
delta - refers to the tilting angle of the substrate
theta - Half angle of the incoming and outgoing orientation. Either 
45 (2theta=90) or 65 for 2theta=130;
%}

theta= pi/4; 
theta_deg=45;
delta_deg=-theta_deg:5:theta_deg;
delta= pi*delta_deg/180;


%{
Incident and reflected angles 
Theta  gives information on the momentum that is passed into the 2d 
substrate. The phi angles gives information on crystal axis along which 
it is passed.
%}
theta_i=pi/2-theta;
phi=0 ;          % incoming and outgoing phi angles are same
theta_o=pi/2-theta;
 

% The crossection dependence is calculated with respect to q_par or delta
%which is the momentum transferred parallel into the surface

k=.5;% Place holder for incoming photon momentum (arbitrarily chosen)
q_par=2*k*sin(theta)*sin(delta);



%{
Plot the q-par dependence w.r.t to delta
plot(delta,q_par);
%}


% Polarization vectors and their angular dependence currently phi-0 case

op_pi_x_in = sin(theta+delta);
op_pi_y_in = 0;
op_pi_z_in = cos(theta+delta);

op_pi_x_out = -sin(theta-delta);
op_pi_y_out = 0;
op_pi_z_out = cos(theta-delta);

op_sigma_x_in = 0;
op_sigma_y_in = -1;
op_sigma_z_in = 0;

op_sigma_x_out = 0;
op_sigma_y_out = -1;
op_sigma_z_out = 0;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

os_u =  cos(theta_s/2).*spinor_up + sin(theta_s/2)*exp(1i*phi_s).*spinor_down;
os_d = -sin(theta_s/2)*exp(-1i*phi_s).*spinor_up + cos(theta_s/2).*spinor_down;


% Since the arbitrary spin is has specific values in spinor representation
% we will now create M matrices for up and down spins with appropriate
% co-factors



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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




%Ground state to dx2y2|spin_down>. However as spins are arbitrary the 
%the dx2y2 state might be a mixture of spin up and down's along the z axis

%Incoming light could be sigma or pi polarized.
%The outgoing light polarization is not measured
%so it is summed over both sigma and pi polarization


%The cross-section response to a particular orbital(d_out) with a
%particular spin(spin_out) is obtained by summing over all 2p orbitals and
%both sigma and pi outgoing polarizations for a given incoming
%polarization.

d_in=dx2y2; %incoming d-orbital
d_out=0; %outgoing d-orbital
s_out=0; % outgoing spin
pol_out = 0; %outgoing polarization
p_orb =0; %transitioning p orbital
cs=0;
cs_t=0;
a_in = os_d; % Spinor Coefficients
cs_f = zeros(5,2,2,19);  %(dout,sout,pol_in,:)

for pol_in=1:1:2
    for d_out=1:1:5
        for s_out=1:1:2
            for pol_out=1:1:2
                for p_orb=1:1:4
                    
                      if s_out == 1
                             a_out = os_d;
                      elseif s_out == 2
                             a_out = os_u;   
                      end
                      
                      if pol_in == 1
                        cs_in_up = M(p_orb,d_in,s_up,p_x).*op_pi_x_in.*a_in(1,1) +M(p_orb,d_in,s_up,p_y).*op_pi_y_in.*a_in(1,1)...
                         + M(p_orb,d_in,s_up,p_z).*op_pi_z_in.*a_in(1,1);
                        cs_in_down = M(p_orb,d_in,s_down,p_x).*op_pi_x_in.*a_in(2,1)...
                         +M(p_orb,d_in,s_down,p_y).*op_pi_y_in.*a_in(2,1) + M(p_orb,d_in,s_down,p_z).*op_pi_z_in.*a_in(2,1);
                      elseif pol_in == 2
                        cs_in_up = M(p_orb,d_in,s_up,p_x)*op_sigma_x_in.*a_in(1,1) +M(p_orb,d_in,s_up,p_y)*op_sigma_y_in.*a_in(1,1)...
                         + M(p_orb,d_in,s_up,p_z)*op_sigma_z_in.*a_in(1,1);
                        cs_in_down = M(p_orb,d_in,s_down,p_x)*op_sigma_x_in.*a_in(2,1)...
                         +M(p_orb,d_in,s_down,p_y)*op_sigma_y_in.*a_in(2,1) + M(p_orb,d_in,s_down,p_z)*op_sigma_z_in.*a_in(2,1);
                      end     
                
                     if pol_out == 1    
                        cs_out_up = M(p_orb,d_out,s_up,p_x).*op_pi_x_out.*a_out(1,1) +M(p_orb,d_out,s_up,p_y).*op_pi_y_out.*a_out(1,1)...
                         + M(p_orb,d_out,s_up,p_z).*op_pi_z_out.*a_out(1,1);
                        cs_out_down = M(p_orb,d_out,s_down,p_x).*op_pi_x_out.*a_out(2,1)...
                         +M(p_orb,d_out,s_down,p_y).*op_pi_y_out.*a_out(2,1) + M(p_orb,d_out,s_down,p_z).*op_pi_z_out.*a_out(2,1);
                     elseif pol_out == 2
                       cs_out_up = M(p_orb,d_out,s_up,p_x)*op_sigma_x_out.*a_out(1,1) +M(p_orb,d_out,s_up,p_y)*op_sigma_y_out.*a_out(1,1)...
                         + M(p_orb,d_out,s_up,p_z)*op_sigma_z_out.*a_out(1,1);
                       cs_out_down = M(p_orb,d_out,s_down,p_x)*op_sigma_x_out.*a_out(2,1)...
                         +M(p_orb,d_out,s_down,p_y)*op_sigma_y_out.*a_out(2,1) + M(p_orb,d_out,s_down,p_z)*op_sigma_z_out.*a_out(2,1);
                     end  
                cs_t = conj(cs_out_up+cs_out_down).*cs_in_up + conj(cs_out_up+cs_out_down).*cs_in_down + cs_t;
                end
                cs_sq = conj(cs_t).*cs_t;  
                cs = cs_sq + cs;
                cs_t=0;
            end
            cs_f(d_out,s_out,pol_in,:)=cs;
            cs=0;
        end
    end
end



%Plotting the cross-sections - Formatting of the plots taken from Umesh's
%code

figure;
box on;

set(0,'DefaultAxesFontName', 'times');
set(0,'DefaultAxesFontSize', 14);
slength=0.8; sheight=0.105;
sxstart=0.125; systart=0.115; 

% For the first 3 orbitals are plotted here
for d_out=1:1:3
    
 pos=[sxstart,systart+(d_out-1)*0.2,slength,sheight];    
 splot(d_out) = subplot('Position', pos);
 
if d_out==1;   title('3z^2-r^2'); xlabel('\delta(Degrees)'); end
if d_out==2;   title('x^2-y^2');  end
if d_out==3;   title('xy'); end

hold on;

%cs_f (d_out,s_out,pol_in,:)
%d-out1:5 same designation as index for M orbitals
%s_out; 1 is down, 2 is up
%pol_in; 1 is pi, 2 is sigma
plot(delta_deg, squeeze(cs_f(d_out,2,2,:)),'b--','linewidth', 1.5,'DisplayName','cos(2x)');
plot(delta_deg, squeeze(cs_f(d_out,2,1,:)),'g-','linewidth', 1.5);
plot(delta_deg, squeeze(cs_f(d_out,1,2,:)),'r--','linewidth', 1.5);
plot(delta_deg, squeeze(cs_f(d_out,1,1,:)),'k-','linewidth', 1.5);

box on;
xlim([-45 45]);
ylim([-0.1 4.5])
hold off;
ylabel('I (a.u.)');

end

%Since dxz and dyz are plotted together they are in a seperate section
d_out=4;
pos=[sxstart,systart+(d_out-1)*0.2,slength,sheight];    
splot(d_out) = subplot('Position', pos);    
title('yz/zx'); 
hold on;
plot(delta_deg, squeeze(cs_f(d_out,2,2,:))+squeeze(cs_f(d_out+1,2,2,:)),'b--','linewidth', 1.5);
plot(delta_deg, squeeze(cs_f(d_out,2,1,:))+squeeze(cs_f(d_out+1,2,1,:)),'g-','linewidth', 1.5);
plot(delta_deg, squeeze(cs_f(d_out,1,2,:))+squeeze(cs_f(d_out+1,1,2,:)),'r--','linewidth', 1.5);
plot(delta_deg, squeeze(cs_f(d_out,1,1,:))+squeeze(cs_f(d_out+1,1,1,:)),'k-','linewidth', 1.5);
box on;
xlim([-45 45]);
ylim([-0.1 4.5])
hold off;
ylabel('I (a.u.)');

box on;
str = '$$ \Theta = 90^\circ,~\phi = 0^\circ $$';
text(0.0,6.5,str,'Interpreter','latex', 'FontSize', 14)
saveas(gcf,'Plotpsi90','epsc');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Obtaining the intensity now as a function of energy lossis done using eqn 18 ( Moretti
%Sala et al 2011 New J of Physics). Intensity responses to d orbitals are
%treated as delta functions which are broadened using the two free
%parameters E_f and Gamma. 

%Gamma represents the lifetime broadening of the states. Equal for orbitals
%with same symmetry
%E_f is the energy of the final state and given spin
%J is the superexchange coupling  (Energy difference between the two spin
%orientations of a given d orbital. J is fixed to 130meV

% All these values are listed on Table 1 in the paper


J=.13;
E_f=zeros(5,1);  % First index corresponds to the d-orbital. Same convention as defined for M matrix
J_f=zeros(5,1);
Gamma_f=zeros(5,1);

J_f(d3z2r2,1)=J/6;
J_f(dx2y2,1)=J;

E_f(d3z2r2,1)=1.97;
E_f(dx2y2,1)=0;
E_f(dxy,1)=1.5;
E_f(dxz,1)=1.84;
E_f(dyz,1)=1.84;

Gamma(d3z2r2,1)=.10;
Gamma(dx2y2,1)=.10;
Gamma(dxy,1)=.08;
Gamma(dxz,1)=.10;
Gamma(dyz,1)=.10;

% Energy loss is identified as w
w_start=1;
w_end=3;
n=1000;    %Number of divisions
w(numel(delta_deg),n)=0;

%To store intensities 
I_s(numel(delta_deg),n)=0;
I_p(numel(delta_deg),n)=0;

w_in=1;

    
for i=1:1:1000
    for d_out=1:1:5
        
        %w-out =0; Therefore w =-w_in
        w(:,i) = -w_in;
        
        f1=repmat(squeeze(cs_f(d_out,2,2,:)),1,n);
        f2=repmat(squeeze(cs_f(d_out,1,2,:)),1,n);
        f3=repmat(squeeze(cs_f(d_out,2,1,:)),1,n);
        f4=repmat(squeeze(cs_f(d_out,1,1,:)),1,n);
        
        I_s(:,i)= I_s(:,i) + Gamma(d_out,1) .* ( f1(:,i)./((w(:,i) + E_f(d_out,1)).^2 + Gamma(d_out,1).^2) ...
        + f2(:,i)./((w(:,i) + E_f(d_out,1) + 2*J_f(d_out)).^2 + Gamma(d_out,1)^2));
    
        I_p(:,i)= I_p(:,i) + Gamma(d_out,1) .* (f3(:,i)./(( w(:,i) + E_f(d_out,1)).^2 + Gamma(d_out,1)^2) ...
        + f4(:,i)./((w(:,i) + E_f(d_out,1) + 2*J_f(d_out)).^2 + Gamma(d_out,1)^2));
    end
    w_in = 1 + (w_end-w_start)*(i)/n;
end



% Again formatting of plots is borrowed from Umesh's code

norm_s(numel(delta_deg),1)=0;

figure;
for ang=1:1:numel(delta)
norm_s(ang,:)=sum(I_s(ang,:),2);
plot(w(ang,:), I_s(ang,:)/(pi*norm_s(ang,:))+.0001*ang,'r-', 'linewidth', 1.0);
hold on;
end

xlabel('Energy Loss (eV)');
ylabel('Intensity (a.u.)');
title('\sigma - Pol: SCOC');
xlabel('Energy Loss (eV)');
ylabel('Intensity (a.u.)');
saveas(gcf,'Plotsigmapolpsi90','epsc');
hold off

norm_p(numel(delta_deg),1)=0;

figure;
for ang=1:1:numel(delta)
norm_p(ang,:)=sum(I_p(ang,:),2);
plot(w(ang,:), I_p(ang,:)/(pi*norm_p(ang,:))+.0001*ang,'r-', 'linewidth', 1.0);
hold on;
end
xlabel('Energy Loss (eV)');
ylabel('Intensity (a.u.)');

title('\pi - Pol: SCOC');
saveas(gcf,'Plotpipolpsi90','epsc');
hold off


%{
k=squeeze(cs_f(2,1,1,:));
plot(q_par,k);
                
%}                
                
                
                




                
                
                
                
                
                
 %{
%cs = zeros(2,5,2,dx2y2,2,


% Incoming light sigma polarized

        dx2y2up3232_sigma = (M(p3232,dx2y2,s_up,p_x)*p_sigma_x_in + M(p3232,dx2y2,s_up,p_y)*p_sigma_y_in + M(p3232,dx2y2,s_up,p_z)*p_sigma_z_in);
        dx2y2up32n12_sigma = (M(p32n12,dx2y2,s_up,p_x)*p_sigma_x_in + M(p32n12,dx2y2,s_up,p_y)*p_sigma_y_in + M(p32n12,dx2y2,s_up,p_z)*p_sigma_z_in);
        dx2y2down3212_sigma = (M(p3212,dx2y2,s_down,p_x)*p_sigma_x_in + M(p3212,dx2y2,s_down,p_y)*p_sigma_y_in + M(p3212,dx2y2,s_down,p_z)*p_sigma_z_in);
        dx2y2down32n32_sigma = (M(p32n32,dx2y2,s_down,p_x)*p_sigma_x_in + M(p32n32,dx2y2,s_down,p_y)*p_sigma_y_in + M(p32n32,dx2y2,s_down,p_z)*p_sigma_z_in);

% Incoming light pi polarized 

        dx2y2up3232_pi = (M(p3232,dx2y2,s_up,p_x)*p_pi_x_in + M(p3232,dx2y2,s_up,p_y)*p_pi_y_in + M(p3232,dx2y2,s_up,p_z)*p_pi_z_in)*(os_down'*spinor_up);
        dx2y2up32n12_pi = (M(p32n12,dx2y2,s_up,p_x)*p_pi_x_in + M(p32n12,dx2y2,s_up,p_y)*p_pi_y_in + M(p32n12,dx2y2,s_up,p_z)*p_pi_z_in)*(os_down'*spinor_up);
        dx2y2down3212_pi = (M(p3212,dx2y2,s_down,p_x)*p_pi_x_in + M(p3212,dx2y2,s_down,p_y)*p_pi_y_in + M(p3212,dx2y2,s_down,p_z)*p_pi_z_in)*(os_down'*spinor_down);
        dx2y2down32n32_pi = (M(p32n32,dx2y2,s_down,p_x)*p_pi_x_in + M(p32n32,dx2y2,s_down,p_y)*p_pi_y_in + M(p32n32,dx2y2,s_down,p_z)*p_pi_z_in)*(os_down'*spinor_down);
            
%}               
                
                
                
      %{          
                
                
            % Outgoing light pi polarized
               dout3232_out    = M(p3232,d_out,s_out,p_x)*p_pi_x_out + M(p3232,d_out,s_out,p_y)*p_pi_y_out + M(p3232,d_out,s_out,p_z)*p_pi_z_out ;      
               dout3212_out    = M(p3212,d_out,s_out,p_x)*p_pi_x_out + M(p3212,d_out,s_out,p_y)*p_pi_y_out + M(p3212,d_out,s_out,p_z)*p_pi_z_out ;  
               dout32n12_out   = M(p32n12,d_out,s_out,p_x)*p_pi_x_out + M(p32n12,d_out,s_out,p_y)*p_pi_y_out + M(p32n12,d_out,s_out,p_z)*p_pi_z_out;  
               dout32n32_out   = M(p32n32,d_out,s_out,p_x)*p_pi_x_out + M(p32n32,d_out,s_out,p_y)*p_pi_y_out + M(p32n32,d_out,s_out,p_z)*p_pi_z_out;  
       
            % Outgoing light sigma polarized
               dout3232_out_s    = M(p3232,d_out,s_out,p_x)*p_sigma_x_out + M(p3232,d_out,s_out,p_y)*p_sigma_y_out + M(p3232,d_out,s_out,p_z)*p_sigma_z_out ;      
               dout3212_out_s    = M(p3212,d_out,s_out,p_x)*p_sigma_x_out + M(p3212,d_out,s_out,p_y)*p_sigma_y_out + M(p3212,d_out,s_out,p_z)*p_sigma_z_out ;  
               dout32n12_out_s   = M(p32n12,d_out,s_out,p_x)*p_sigma_x_out + M(p32n12,d_out,s_out,p_y)*p_sigma_y_out + M(p32n12,d_out,s_out,p_z)*p_sigma_z_out;  
               dout32n32_out_s   = M(p32n32,d_out,s_out,p_x)*p_sigma_x_out + M(p32n32,d_out,s_out,p_y)*p_sigma_y_out + M(p32n32,d_out,s_out,p_z)*p_sigma_z_out; 
               
             %Define the elements for each channel/incoming sigma/outgoing pi
             c_pi = conj(dout3232_out).*dx2y2up3232_sigma.*(spinor_down'*spinor_up) + conj(dout32n12_out).*dx2y2up32n12_sigma.*(spinor_down'*spinor_up)  + conj(dout3212_out).*dx2y2down3212_sigma.*(spinor_down'*spinor_down) + conj(dout32n32_out).*dx2y2down32n32_sigma.*(spinor_down'*spinor_down);
            %Define the elements for each channel/incoming sigma/outgoing sigma
             c_sigma = conj(dout3232_out_s).*dx2y2up3232_sigma.*(spinor_down'*spinor_up)  +  conj(dout32n12_out_s).*dx2y2up32n12_sigma.*(spinor_down'*spinor_up)  + conj(dout3212_out_s).*dx2y2down3212_sigma.*(spinor_down'*spinor_down) + conj(dout32n32_out_s).*dx2y2down32n32_sigma.*(spinor_down'*spinor_down);         
       elseif spin==1
            %Define the elements for each channel/incoming sigma/outgoing pi
             c_pi = conj(dout3232_out).*dx2y2up3232_sigma.*(spinor_up'*spinor_up) + conj(dout32n12_out).*dx2y2up32n12_sigma.*(spinor_up'*spinor_up)  + conj(dout3212_out).*dx2y2down3212_sigma.*(spinor_up'*spinor_down) + conj(dout32n32_out).*dx2y2down32n32_sigma.*(spinor_up'*spinor_down);
            %Define the elements for each channel/incoming sigma/outgoing sigma
             c_sigma = conj(dout3232_out_s).*dx2y2up3232_sigma.*(spinor_up'*spinor_up)  + conj(dout32n12_out_s).*dx2y2up32n12_sigma.*(spinor_up'*spinor_up)  + conj(dout3212_out_s).*dx2y2down3212_sigma.*(spinor_up'*spinor_down) + conj(dout32n32_out_s).*dx2y2down32n32_sigma.*(spinor_up'*spinor_down); 
           
           
       end
       %{
       c5_p = conj(dout3232_out_s).*dx2y2up3232_pi + conj(dout3232_out).*dx2y2up32n12_pi + conj(dout3232_out).*dx2y2down3212_pi + conj(dout3232_out).*dx2y2down32n32_pi; 
       c6_p = conj(dout3212_out_s).*dx2y2up3232_pi + conj(dout3212_out).*dx2y2up32n12_pi + conj(dout3212_out).*dx2y2down3212_pi + conj(dout3212_out).*dx2y2down32n32_pi;
       c7_p = conj(dout32n12_out_s).*dx2y2up3232_pi + conj(dout32n12_out).*dx2y2up32n12_pi + conj(dout32n12_out).*dx2y2down3212_pi + conj(dout32n12_out).*dx2y2down32n32_pi;
       c8_p = conj(dout32n32_out_s).*dx2y2up3232_pi + conj(dout32n32_out).*dx2y2up32n12_pi + conj(dout32n32_out).*dx2y2down3212_pi + conj(dout32n32_out).*dx2y2down32n32_pi;
       %}
     
       
       cs = conj(c_sigma).*c_sigma + conj(c_pi).*c_pi;
       plot(q_par,cs,'g');
       
   %     end
%end

           
%{

c1=(conj(dx2y2up3232_sigma_out).*dx2y2up3232_in);
c2 = (conj(dx2y2up3232_out).*dx2y2up3232_in);
c3 = (conj(dx2y2up32n12_sigma_out).*dx2y2up3232_in);
c4 = (conj(dx2y2up32n12_out).*dx2y2up3232_in);
c5 = c1+c3;
c6 = c2+c4;
cs = (conj(c5).*c5 + conj(c6).*c6)/5.6;

plot(q_par,cs,'g');

%}
%}
%{
% Incoming light pi polarized \\Same as outgoing light polarization

dx2y2up3232 = M(p3232,dx2y2,s_up,p_x)*p_pi_x_in + M(p3232,dx2y2,s_up,p_y)*p_pi_y_in + M(p3232,dx2y2,s_up,p_z)*p_pi_z_in; 
dx2y2up32n12= M(p32n12,dx2y2,s_up,p_x)*p_pi_x_in + M(p32n12,dx2y2,s_up,p_y)*p_pi_y_in+ M(p32n12,dx2y2,s_up,p_z)*p_pi_z_in;  % We know that p_pi_z=0

% Incoming light sigma polarized

dx2y2up3232_sigma = M(p3232,dx2y2,s_up,p_x)*p_sigma_x_in + M(p3232,dx2y2,s_up,p_y)*p_sigma_y_in + M(p3232,dx2y2,s_up,p_z)*p_sigma_z_in;
dx2y2up32n12_sigma = M(p32n12,dx2y2,s_up,p_x)*p_sigma_x_in + M(p32n12,dx2y2,s_up,p_y)*p_sigma_y_in + M(p32n12,dx2y2,s_up,p_z)*p_sigma_z_in;


% Outgoing light sigma polarized

dx2y2up3232_sigma_out = M(p3232,dx2y2,s_up,p_x)*p_sigma_x_out + M(p3232,dx2y2,s_up,p_y)*p_sigma_y_out + M(p3232,dx2y2,s_up,p_z)*p_sigma_z_out;
dx2y2up32n12_sigma_out = M(p32n12,dx2y2,s_up,p_x)*p_sigma_x_out + M(p32n12,dx2y2,s_up,p_y)*p_sigma_y_out + M(p32n12,dx2y2,s_up,p_z)*p_sigma_z_out;


c1=(conj(dx2y2up3232_sigma_out).*dx2y2up3232_sigma);
c2 = (conj(dx2y2up3232).*dx2y2up3232_sigma)
c3 = (conj(dx2y2up32n12_sigma_out).*dx2y2up3232_sigma);
c4 = (conj(dx2y2up32n12).*dx2y2up3232_sigma);
c5 = c1+c3;
c6 = c2+c4;
cs = ((conj(c5).*c5 + conj(c6).*c6))/5;

plot(q_par,cs,'g');



%}



%Previous version of the polarization angles
%{ 
So far we have defined all of the necessary paramters associated with the
experiment. However our cross-section values are all w.r.t x, y z
polarizations alone. With delta and phi being our probing parameters we
will from now on use sigma and pi polarizations. These can then be
expressed in the x-y-z polarization basis with appropriate prefactors.



p_pi_x_in = cos(phi)*cos(theta_i);
p_pi_y_in = -sin(delta).*sin(theta_i) + cos(delta).*sin(phi).*cos(theta_i);
p_pi_z_in = cos(delta).*sin(theta_i)+ sin(delta).*sin(phi).*cos(theta_i);

p_pi_x_out = -cos(phi).*cos(theta_o);
p_pi_y_out = -sin(delta).*sin(theta_o) - cos(delta).*sin(phi).*cos(theta_o);
p_pi_z_out = cos(delta).*sin(theta_o) - sin(delta).*sin(phi).*cos(theta_o);

p_sigma_x_in = -sin(phi);
p_sigma_y_in = cos(phi).*cos(delta);
p_sigma_z_in = sin(delta).*cos(phi);

p_sigma_x_out = -sin(phi);
p_sigma_y_out = cos(phi).*cos(delta);
p_sigma_z_out = sin(delta).*sin(phi);

%}














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
