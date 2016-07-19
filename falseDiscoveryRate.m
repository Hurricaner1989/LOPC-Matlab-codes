function p_adjusted = falseDiscoveryRate(n, p)
% Compute adjusted p-value matrix based on FDR from p-value matrix
% 
% Depend on mafdr function from bioinformatics toolbox
% 
% Inputs:
% 'N'    number of variables
% 'P'    N /times N p-value matrix
% Outputs:
% 'r_2'  second order partial correlation matrix
% 'p_2'  second order partial correlation p-value matrix 
% 
% References:
% [1] Zuo, Yiming, Guoqiang Yu, Mahlet G. Tadesse, and Habtom W. Ressom. 
%     Biological network inference using low order partial correlation. 
%     Methods 69, no. 3 (2014): 266-273.
%
% Copyright 2014-2016, Yiming Zuo.

% calculate FDR
p_tri = zeros(n*(n-1)/2, 1);
num = 1;
for i = 1:n-1
    for j = i+1:n
        p_tri(num) = p(i, j);
        num = num+1;
    end
end
[~, q] = mafdr(p_tri);

% create adjusted p-value matrix
p_adjusted = ones(n, n);
num = 1;
for i = 1:n-1
    for j = i+1:n
        p_adjusted(i, j) = q(num);
        num = num+1;
    end
end