%% create_grating_2D.m
%
% DESCRIPTION:
%   generates grating of different kinds (2D) (pi phase shift, absorption etc.)
%   from different materials. A design energy is provided in order to
%   calculate the height of the structure. From the specified design
%   values the dehaviour of the grating is calculated for an opertion energy
%   . At the moment the absorption grating and source gratings are
%   considered to be ideal. The grating lines are allong the x direction (DIM 1)
%
%
% CALL: G = create_grating(type,material,E_design,E_op,L,period,dc,angle)
%   - type: specify type of grating; G1_pi, G1_pi_half, G0, G2
%   - material: for G1 'Si' more options will be added, for G0 and G2 set to
%   []
%   - E_design: design energy
%   - E_op: operation energy
%   - L: vector containing the field of view in both directions, DIM 1 x,
%   DIM 2 y
%   - period: pitch of grating
%   - dc: duty cycle
%   - angle: anlge between the grating lines and the vertical direction
%   - G: grating as a transmission function,m DIM 1 energy dimension, DIM 2 spatial dimension
%
%
%
% UPDATES:
%   02.10.2013 (Matias): first version
%
%%
function G = create_grating_2D(type,material,E_design,E_op,x,y,period,dc,angle)

    % constants
    c = 299792458;
    h = 4.135*10^(-15);


    % define mask 
    
    k = sin((pi-2*pi*dc)/2);
    mask = (sin(2*pi*(x*cos(angle)-y*sin(angle))/period)>k);


    switch type
        case 'G1_pi'
            switch material
                case 'Si'
                    H = c*h/2/Si_delta_LUT(E_design)/E_design;
                    Phase = 2*pi*Si_delta_LUT(E_op)*H.*E_op/c/h;
                    Phase = pi*ones(size(mask));
                    G = exp(1i*Phase.*mask);
            end

        case 'G1_pi_half'
            switch material
                case 'Si'
                    H = c*h/2/2/Si_delta_LUT(E_design)/E_design;
                    Phase = 2*pi*Si_delta_LUT(E_op)*H.*E_op/c/h;
                    Phase = Phase*ones(size(mask));
                    G = exp(1i*Phase.*mask);
            end


        case 'G2'

            G = mask;

        case 'G0'

            G = mask;



    end

end