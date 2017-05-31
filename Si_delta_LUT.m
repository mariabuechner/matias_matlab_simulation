%% Si_delta_LUT
% 
% DESCRIPTION:
%   LUT for Si
%
%
%%
function delta=Si_delta_LUT(E)
    E_tab=[1 5 10 15 20 25 30 35 40]*(1e+3);% 45 50 55 60 65 70 75 80 85 90 95 100 105 110 115 120];
    delta_tab=[0.88413e-03 0.39495e-04 0.97631e-05 0.43171e-05 0.24223e-05 0.15481e-05 0.10741e-05 0.78871e-06 ...
         0.60360e-06];
    delta=interp1(E_tab,delta_tab,E,'cubic');
end
     