%% sample_simulation_file.m
%
% DESCRIPTION: just an example simulating a GI
%
%
%
%
%
%
%
%
%
%% Initialization parameters
clear;
% Constants
c = 299792458;
h = 4.135*10^(-15);
delta_PMMA = 0.45e-6;
beta_PMMA = 0.6828e-10;


% Spectral parameters
E_central = 25000;
BW = 0.00001;
N_bands = 3;
DE = E_central*BW;
Ei = linspace(E_central-DE/2,E_central+DE/2,N_bands);
I_E = exp(-((2.35*(Ei-E_central)).^2)/2/(DE^2));
I_E = I_E/sum(I_E);
lambda_central = h*c/E_central;
lambdai = h*c./Ei;
ki = 2*pi./lambdai;
k_central = 2*pi/lambda_central;


% Imaging parameters
FOV = 8000*1e-6;
pxs = 8*1e-6;
Nph = 5;

% Gi parameters
g1 = 4*1e-6; % Period of phase grating
g2 = g1/2;
dc = 0.5;
z = (1-1/2)*g1^2/4/lambda_central; % Intergrating distance pi

% numerical parameters
N = 1e+6; %total number of points for the FOV
x = linspace(0,FOV,N);
dx = diff(x(1:2));
dz = dx;


%define source
source_size = 124.6e-6;
proj_source_size = source_size./2.355*z/22;
sconv = exp(-(x-FOV/2).^2/2./proj_source_size.^2);
sconv = sconv./sum(sconv); % Source Kernel







%% sample part! Define your sample here

Sph = real(k_central*2*sqrt((3000e-6)^2-(x-FOV/2).^2));
Ph = Sph*delta_PMMA+1i*Sph*beta_PMMA;
Ds = exp(1i*Ph);


Ds = repmat(Ds,3,1);






%% define G1

disp('Create gratings')
tic
G1 = create_grating('G1_pi','Si',Ei,E_central,x,g1,dc);
toc

%% propagate wave

disp('Propagate wave')
tic
D_flat = fresnel_propagation_poly_1D(G1,FOV,lambdai,z);
D_samp = fresnel_propagation_poly_1D(G1.*Ds,FOV,lambdai,z);
toc

%% phase stepping and detection

DQE = ones(size(Ei));

disp('Phase stepping')
tic
[PSC_flat,PSC_samp] = phase_stepping_1D(D_flat,D_samp,Nph,Ei,E_central,x,g2,dc,I_E,DQE,pxs,sconv,14,1);
toc


%% FCA

disp('Get signals')
tic
[A,B,P] = FCA(PSC_flat,PSC_samp,1);
toc

%% plot
subplot(311)
plot(A);title('Abs')
subplot(312)
plot(B);title('Drk')
subplot(313)
plot(P);title('DPC')





 


