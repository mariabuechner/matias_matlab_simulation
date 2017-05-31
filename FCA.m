%% FCA.m
%
% DESCRIPTION: simple function performing fourier component analysis based
%   on Peters fba.m function 
%
% CALL: [A,B,P] = FCA(flat,samp,periods)
%   - flat: flat PSCs
%   - samp: samp PSCs
%   - periods: number of periods of PSCs
%   - A: absorption image
%   - B: dark field image
%   - P: DPC image 
%
%
% UPDATES:
%   01.10.2013 (Matias) : fixed typo
%   30.09.2013 (Matias) : first version
%
%%
function [A,B,P] = FCA(flat,samp,periods)
    
    
    f_flat = fft(flat,[],2);
    f_samp = fft(samp,[],2);
    A_f = abs(f_flat(:,1));
    A_s = abs(f_samp(:,1));
    B_f = abs(f_flat(:,periods+1));
    B_s = abs(f_samp(:,periods+1));
    P_f = angle(f_flat(:,periods+1));
    P_s = angle(f_samp(:,periods+1));
    
    A = A_s./A_f;
    B = B_s./B_f./A;
    P = angle(exp(1i*(P_f-P_s)));
    
end