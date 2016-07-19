function [r_1, p_1] = firstParCorr(n, r, m)
% Compute first order partial correlation
%
% Depend on normcdf function from statistics and machine learning toolbox
%
% Inputs:
% 'n'    number of variables
% 'r'    Pearson's correlation matrix
% 'm'    number of samples
% Outputs:
% 'r_1'  first order partial correlation matrix
% 'p_1'  first order partial correlation p-value matrix 
%
% References:
% [1] Zuo, Yiming, Guoqiang Yu, Mahlet G. Tadesse, and Habtom W. Ressom. 
%     Biological network inference using low order partial correlation. 
%     Methods 69, no. 3 (2014): 266-273.
%
% Copyright 2014-2016, Yiming Zuo.

% first-order partial correlation
r_1 = ones(n,n); % r_1 is by default set to 1
for i = 1:n-1
    for j = i+1:n
        for k = 1:n
            if (k~=i) && (k~=j)
                r1 = (r(i,j)-r(i,k)*r(j,k)) / sqrt((1-r(i,k)^2)*(1-r(j,k)^2));
                if abs(r1) < abs(r_1(i, j)) % to get the smallest first order partial correlation
                    r_1(i, j) = r1;
                end
            end
        end
        r_1(j, i) = r_1(i, j);
    end
end

% Fisher transformation for first-order partial correlation
z_1 = 0.5 * log((1+r_1)./(1-r_1));
p_1 = ones(n, n);
for i = 1:n-1
    for j = i+1:n
        p_1(i,j) = (1- normcdf(abs(z_1(i,j)), 0, sqrt(1/(m-3-1)))) * 2;
    end
end
