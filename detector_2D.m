%% detector_2D.m
%
% DESCRIPTION:
%   simulates a simple 1D pixel detector. Poisson noise included  
%
%
% CALL : [S_flat,S_samp] = detector(Di,DQE,FOV,pxs,sconv,G2)
%   - Di_flat: incident flat wave; DIM 1 x, DIM 2 y
%   - Di_samp: incident sample wave; DIM 1 x, DIM 2 y
%   - DQE: detector quantum efficiency at each energy
%   - pxs: pixel size
%   - sconv: source distribution (due to finite source size)
%   - G2: transmition function of G2 (fix this)
%   - S_flat: flat detected signal
%   - S_samp: samp detected signal
%
%
% UPDATES:
%   04.10.2013 (Matias): first version
%
%
%%
function [S_flat,S_samp] = detector_2D(Di_flat,Di_samp,DQE,FOV,pxs,sconv,G2)
    

  

    npx_x = length(0:pxs:FOV(1));
    npx_y = length(0:pxs:FOV(2));

    
    
    % calculate intensities
    I_flat = abs(Di_flat).^2; 
    I_samp = abs(Di_samp).^2;
    
%     S_flat_tmp = I_flat;
%     S_samp_tmp = I_samp;
    
    
    % introduce finite source size
    S_flat_tmp = abs(fftshift(ifft2(fft2(I_flat).*fft2(sconv))));
    S_samp_tmp = abs(fftshift(ifft2(fft2(I_samp).*fft2(sconv))));
    
    
    % find a more elegant solution for this
    
    
    if ~isempty(G2)
        S_flat_tmp = S_flat_tmp.*G2;
        S_samp_tmp = S_samp_tmp.*G2;
    end
    
    % take into consideration DQE  
    S_flat_tmp = DQE.*S_flat_tmp;
    S_samp_tmp = DQE.*S_samp_tmp;
    
   
    
    
    %pixel indices
    pxl_x=floor(size(Di_flat,1)/npx_x);
    pxl_y=floor(size(Di_flat,2)/npx_y);
    ipx_x=1:npx_x;
    ipx_y=1:npx_y;
    I1p_x=floor(1+(ipx_x-1)*pxl_x);
    I2p_x=ceil(pxl_x*(ipx_x+0));
    I1p_y=floor(1+(ipx_y-1)*pxl_y);
    I2p_y=ceil(pxl_y*(ipx_y+0));
    
    if I2p_x(end) < size(Di_flat,1)
        S_flat_tmp(I2p_x(end)+1:end,:) = [];
        S_samp_tmp(I2p_x(end)+1:end,:) = [];      
    else
        S_flat_tmp(end:end+I2p_x(1)-I1p_x(2)-(I2p_x(end)-I1p_x(end)),:) = 0;
        S_samp_tmp(end:end+I2p_x(1)-I1p_x(2)-(I2p_x(end)-I1p_x(end)),:) = 0;
    end
    
    
    if I2p_y(end) < size(Di_flat,2)
        S_flat_tmp(:,I2p_y(end)+1:end) = [];
        S_samp_tmp(:,I2p_y(end)+1:end) = [];      
    else
        S_flat_tmp(:,end:end+I2p_y(1)-I1p_y(2)-(I2p_y(end)-I1p_y(end))) = 0;
        S_samp_tmp(:,end:end+I2p_y(1)-I1p_y(2)-(I2p_y(end)-I1p_y(end))) = 0;
    end
    
    
    % calculate signal at each pixel 
    
    S_flat = reshape(S_flat_tmp,pxl_x,npx_x*pxl_y*npx_y);
    S_samp = reshape(S_samp_tmp,pxl_x,npx_x*pxl_y*npx_y);
     
   
    S_flat = sum(S_flat,1);
    S_samp = sum(S_samp,1);
    
    S_flat = reshape(S_flat,npx_x,pxl_y*npx_y);
    S_samp = reshape(S_samp,npx_x,pxl_y*npx_y);
    
    S_flat = S_flat';
    S_samp = S_samp';
    
    S_flat = reshape(S_flat,pxl_y,npx_x*npx_y);
    S_samp = reshape(S_samp,pxl_y,npx_x*npx_y);
    
    
    S_flat = sum(S_flat,1);
    S_samp = sum(S_samp,1);
    
        
    S_flat = reshape(S_flat,npx_y,npx_x);
    S_samp = reshape(S_samp,npx_y,npx_x);
    
        
    S_flat = S_flat';
    S_samp = S_samp';

    

end


