%% detector.m
%
% DESCRIPTION:
%   simulates a simple 1D pixel detector. Poisson noise included  
%
%
% CALL : [S_flat,S_samp] = detector(Di,I_E,DQE,FOV,pxs,sconv)
%   - Di_flat: incident flat wave; DIM 1 energy dimension, DIM 2 spatial dimension
%   - Di_SAMP: incident sample wave; DIM 1 energy dimension, DIM 2 spatial dimension
%   - I_E: vector containing the intensities for each energy
%   - DQE: vectror containing the detector quantum efficiency at each energy
%   - pxs: pixel size
%   - sconv: source distribution (due to finite source size)
%   - nbits: number of bits 
%   - chi: proportionality factor depending on detector properties
%   - S: detected signal
%   - G2: transmition function of G2 (fix this)
%
%
% UPDATES:
%   01.10.2013 (Matias): removed noise
%   01.10.2013 (Matias): fixed memory issue
%   30.09.2013 (Matias): first version
%
%
%%
function [S_flat,S_samp] = detector(Di_flat,Di_samp,I_E,DQE,FOV,pxs,sconv,G2)
    

  

    npx = length(0:pxs:FOV);

    
    % 
    DQE = repmat(DQE',1,size(Di_flat,2));
    I_E = repmat(I_E',1,size(Di_flat,2));
    
    % calculate intensities
    I_flat = abs(Di_flat).^2; 
    I_samp = abs(Di_samp).^2;
    
    sconv = repmat(sconv,size(I_flat,1),1);
    
    % introduce finite source size
    S_flat_tmp = abs(squeeze(fftshift(ifft(fft(I_flat,[],2).*fft(sconv,[],2),[],2),2)));
    S_samp_tmp = abs(squeeze(fftshift(ifft(fft(I_samp,[],2).*fft(sconv,[],2),[],2),2)));
    
    % find a more elegant solution for this
    
    
    if ~isempty(G2)
        S_flat_tmp = S_flat_tmp.*G2;
        S_samp_tmp = S_samp_tmp.*G2;
    end
    
    % take into consideration DQE and intensity at every energy band 
    S_flat_tmp = squeeze(sum(DQE.*I_E.*S_flat_tmp,1));
    S_samp_tmp = squeeze(sum(DQE.*I_E.*S_samp_tmp,1));
    
   
    
    
    %pixel indices
    pxl=floor(size(Di_flat,2)/npx);
    ipx=1:npx;
    I1p=floor(1+(ipx-1)*pxl);
    I2p=ceil(pxl*(ipx+0));
    
    if I2p(end) < size(Di_flat,2)
        S_flat_tmp(I2p(end)+1:end) = [];
        S_samp_tmp(I2p(end)+1:end) = [];
        
    else
        S_flat_tmp(end:end+I2p(1)-I1p(2)-(I2p(end)-I1p(end))) = 0;
        S_samp_tmp(end:end+I2p(1)-I1p(2)-(I2p(end)-I1p(end))) = 0;
    end
    
    
    
    S_flat = reshape(S_flat_tmp,pxl,npx);
    S_samp = reshape(S_samp_tmp,pxl,npx);
    

    
    % calculate signal at each pixel 
    S_flat = sum(S_flat,1);
    S_samp = sum(S_samp,1);
    
   
    

end


