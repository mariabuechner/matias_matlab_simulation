%% phase_stepping_1D.m
%
% DESCRIPTION: this function performs the phase stepping and the detection
% of the signal. Poisson noise is also added at this step.
%
%
%
% CALL: [PSC_flat,PSC_samp] = phase_stepping_1D(Di_flat,Di_samp,Nph,E,E_design,x,g2,dc,I_E,DQE,pxs,sconv,nbits,chi)
%   - Di_flat: incident flat wave; DIM 1 energy dimension, DIM 2 spatial dimension
%   - Di_samp: incident sample wave; DIM 1 energy dimension, DIM 2 spatial dimension
%   - Nph: number of phase steps
%   - E: vector containing the different energy bands
%   - E_design: nominal energy
%   - x: spatial coordinates
%   - g2: period of absorption grating
%   - dc: duty cycle of absorption grating
%   - I_E: vector containing the intensities for each energy
%   - DQE: vectror containing the detector quantum efficiency at each energy
%   - pxs: pixel size
%   - sconv: source distribution (due to finite source size)
%   - nbits: number of bits 
%   - chi: proportionality factor depending on detector properties
%   - PSC_flat: flat phase stepping curve
%   - PSC_samp: sample phase stepping curve
%
%
%
%
% UPDATES:
%   01.10.2013 (Matias) : first version 
%
%
%
%
%% 
function [PSC_flat,PSC_samp] = phase_stepping_1D(Di_flat,Di_samp,Nph,E,E_design,x,g2,dc,I_E,DQE,pxs,sconv,nbits,chi)

    FOV = max(x);
    % define phase steps 
    s = linspace(0,g2,Nph+1);
    s = s(1:end-1);
    
    
    PSC_flat = zeros(length(0:pxs:FOV),Nph);
    PSC_samp = zeros(length(0:pxs:FOV),Nph);
    
    
    for i=1:length(s)
        
        G2 = create_grating('G2',[],E,E_design,x-s(i),g2,dc);
    
        % readout signal
        [PSC_flat(:,i),PSC_samp(:,i)] = detector(Di_flat,Di_samp,I_E,DQE,FOV,pxs,sconv,G2);
        
        
       
    end
    
    
    % add noise 
    Sc = repmat(squeeze(max(PSC_flat,[],2)),1,size(PSC_flat,2));
    PSC_flat = (2^nbits)*PSC_flat./Sc;
    PSC_samp = (2^nbits)*PSC_samp./Sc;
    PSC_flat = countstat(PSC_flat,chi);
    PSC_samp = countstat(PSC_samp,chi);
    

%     x_values = [1:67108864];
%     tic
%     y_values = sin(x_values);
%     toc
%     
%     tic
%     fft_values = fft(y_values);
%     toc
%   
    
   


end