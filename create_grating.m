%% create_grating.m
%
% DESCRIPTION:
%   generates grating of different kinds (1D),v  (pi phase shift, absorption etc.)
%   from different materials. A design energy is provided in order to
%   calculate the height of the structure. From the specified design
%   values the dehaviour of the grating is calculated for a given energy
%   range. At the moment the absorption grating and source gratings are
%   considered to be ideal.
%
%
% CALL: G = create_grating(type,material,E,E_design,x,period,dc)
%   - type: specify type of grating; G1_pi, G1_pi_half, G0, G2
%   - material: for G1 'Si' more options will be added, for G0 and G2 set to
%   []
%   - E: vector containing energy range
%   - E_design: design energy
%   - x: spatial range for which grating is designed
%   - period: pitch of grating
%   - dc: duty cycle
%   - G: grating as a transmission function,m DIM 1 energy dimension, DIM 2 spatial dimension 
%
%
%
% UPDATES:
%   30.09.2013 (Matias): first version
%
%%
function G = create_grating(type,material,E,E_design,x,period,dc)

    % constants
    c = 299792458;
    h = 4.135*10^(-15);


    % 
    k = sin((pi-2*pi*dc)/2);
    mask = (sin(2*pi*x/period)>k);
    mask = repmat(mask,length(E),1);


    switch type
        case 'G1_pi'
            switch material
                case 'Si'
                    H = c*h/2/Si_delta_LUT(E_design)/E_design;
                    Phase = 2*pi*Si_delta_LUT(E)*H.*E/c/h;
                    Phase = repmat(Phase',1,length(x));
                    G = exp(1i*Phase.*mask);
            end

        case 'G1_pi_half'
            switch material
                case 'Si'
                    H = c*h/2/2/Si_delta_LUT(E_design)/E_design;
                    Phase = 2*pi*Si_delta_LUT(E)*H.*E/c/h;
                    Phase = repmat(Phase',1,length(x));
                    G = exp(1i*Phase.*mask);
            end


        case 'G2'

            G = mask;

        case 'G0'

            G = mask;




    end

end