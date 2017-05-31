%% fresnel_propagation_mono_2D.m
%
% DESCRIPTION: monochromatic fresnel propagation for 2D waves
% 
% CALL: Do = fresnel_propagation_mono_2D(Di,L,lambda,z)
%   - Di: intput wave, DIM 1 x and DIM 2 y   
%   - L : vector containing field of view for both dimention. DIM 1 x and
%   DIM 2 y
%   - lambda: wavelength
%   - z : propagation distance
%   - Do: output wave
%
%
% UPDATES:
%   02.10.2013 (Matias) : first version
%
%%
function Do = fresnel_propagation_mono_2D(Di,L,lambda,z) 
     
     M = size(Di); %get input field array size
     dXY = L./M; %sample interval
     fx = -1/(2*dXY(1)):1/L(1):1/(2*dXY(1))-1/L(1); %freq coords
     fy = -1/(2*dXY(2)):1/L(2):1/(2*dXY(2))-1/L(2); %freq coords
     
     [fx,fy] = meshgrid(fx,fy);
     fx = fx';
     fy = fy';

     H = exp(-1i*pi*lambda*z.*(fx.^2+fy.^2)); %trans func
     H = ifftshift(H); %shift trans func
     di = fft2(fftshift(Di)); %shift, fft src fields
     do = H.*di; %multiply
     Do = ifftshift(ifft2(do)); %inv fft, center obs field
     
 end