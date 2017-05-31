%% phase_stepping_2D.m
%
% DESCRIPTION: this function performs the phase stepping and the detection
% of the signal. Poisson noise is also added at this step.
%
%
%
% CALL: [PSC_flat,PSC_samp] = phase_stepping_2D(Di_flat,Di_samp,Nph,E_op,E_design,x,y,g2,dc,angle,DQE,pxs,sconv,nbits,chi)
%   - Di_flat: incident flat wave; DIM 1 energy dimension, DIM 2 spatial dimension
%   - Di_samp: incident sample wave; DIM 1 energy dimension, DIM 2 spatial dimension
%   - Nph: number of phase steps
%   - E_op: operational energy
%   - E_design: nominal energy
%   - x: spatial coordinates in x direction
%   - y: spatial coordinates in y direction
%   - g2: period of absorption grating
%   - dc: duty cycle of absorption grating
%   - angle: angle between absorbtion grating and vetical direction 
%   - DQE: detector quantum efficiency
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
%   04.10.2013 (Matias) : first version 
%
%
%
%
%% 
function [PSC_flat,PSC_samp] = phase_stepping_2D(Di_flat,Di_samp,Nph,E_op,E_design,x,y,g2,dc,angle,DQE,pxs,sconv,nbits,chi)

    FOV(1) = max(max(x));
    FOV(2) = max(max(y));
    % define phase steps 
    s = linspace(0,g2,Nph+1);
    s = s(1:end-1);
    
    
    PSC_flat = zeros(length(0:pxs:FOV(1)),length(0:pxs:FOV(2)),Nph);
    PSC_samp = zeros(length(0:pxs:FOV(1)),length(0:pxs:FOV(2)),Nph);
    
    
    for i=1:length(s)

        G2 = create_grating_2D('G2',[],E_op,E_design,x-s(i),y,g2,dc,angle);
    
        % readout signal
        [PSC_flat(:,:,i),PSC_samp(:,:,i)] = detector_2D(Di_flat,Di_samp,DQE,FOV,pxs,sconv,G2);
        
        
       
    end
    
    
    % add noise 
    Sc = repmat(squeeze(max(PSC_flat(:))),size(PSC_flat));
    PSC_flat = (2^nbits)*PSC_flat./Sc;
    PSC_samp = (2^nbits)*PSC_samp./Sc;
    PSC_flat = countstat(PSC_flat,chi);
    PSC_samp = countstat(PSC_samp,chi);
    

        
    
  
    
   


end