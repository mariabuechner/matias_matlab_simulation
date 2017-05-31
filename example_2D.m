%% example_2D.m
%
% DESCRIPTION: just an example simulating a 2D GI
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


% Spectral parameters
E_central = 25000;
lambda_central = h*c/E_central;
k_central = 2*pi/lambda_central;


% Imaging parameters
dx = 1e-7;
dy = 1e-7;
FOV = [dx*2^13 dy*2^13];
pxs = 8*1e-6;
Nph = 9;

% Gi parameters
g1 = 4e-6; % Period of phase grating
g2 = g1/2;
dc = 0.5;
z = (2-1/2)*g1^2/4/lambda_central; % Intergrating distance pi

% numerical parameters
N = FOV./[dx dy]; %total number of points for the FOV
[x,y] = meshgrid(linspace(0,FOV(1),N(1)),linspace(0,FOV(2),N(2)));
x = x';
y = y';
dz = dx;


%define source
source_size = [124.6e-6 40e-6];
proj_source_size = source_size./2.355*z/22;
sconv = exp(-(x-FOV(1)/2).^2/2./proj_source_size(1).^2-(y-FOV(2)/2).^2/2./proj_source_size(2).^2);
sconv = sconv./sum(sum(sconv)); % Source Kernel







%% sample part! Define your sample here


delta = 0.45e-6;
beta = 0.6828e-9;

Sph = real(k_central*2*sqrt((240e-6)^2-(x-FOV(1)/2).^2-(y-FOV(2)/2).^2));
Ph = Sph*delta+1i*Sph*beta;
Ds = exp(1i*Ph);






%% define G1

disp('Create gratings')
tic
G1 = create_grating_2D('G1_pi','Si',E_central,E_central,x,y,g1,dc,0*pi/180);
toc

%% propagate wave

disp('Propagate wave')
tic
D_flat = fresnel_propagation_mono_2D(G1,FOV,lambda_central,z);
D_samp = fresnel_propagation_mono_2D(Ds.*G1,FOV,lambda_central,z);
toc

%% phase stepping and detection

DQE = 1;

disp('Phase stepping')
tic
[PSC_flat,PSC_samp] = phase_stepping_2D(D_flat,D_samp,Nph,E_central,E_central,x,y,g2,dc,0,DQE,pxs,sconv,14,1);
toc


%% FCA

disp('Get signals')
tic
[A,B,P] = FCA_2D(PSC_flat,PSC_samp,1);
toc

%% plot
subplot(221)
imagesc(A);title('Abs');
subplot(222)
imagesc(B);title('Drk')
subplot(223)
imagesc(P);title('DPC')





 


