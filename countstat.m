% [data,pdata] = countstat(data,fac)
%
% transforms a given data into countstatistcs with poisson statistics
% each new value is mynrand with mu = value; sigma = fac*sqrt(value).
% Quantization noise is included by rounding the data.

function [data,pdata] = countstat(data,fac)

if nargin < 2,
    fac = 1;
end

pdata = round(fac*sqrt(data).*randn(size(data)));
data = pdata+data;