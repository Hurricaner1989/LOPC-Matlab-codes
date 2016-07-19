function [r_2, p_2] = lopc(data, thres)
% Implement LOPC algorithm to build networks based on low order partial
% correlation up to the second order
%
% Depend on mafdr function from bioinformatics toolbox
% Depend on normcdf function from statistics and machine learning toolbox
%
% Input parameters:
%     'data'    An m \times n numeric matrix, where m is the number of
%               samples and n is the number of features.
%     'thres'   A numeric value for picking out statistically significant
%               pairwise association for 0-th order partial correlation 
%               (i.e., Pearson's correlation) and first order partial 
%               correlation after correcting multiple testing problem 
%               using False Discovery Rate (FDR). Recommended values 
%               are 0.01, 0.05 and 0.1. 
% Output parameters:
%     'r_2'     An n \times n numeric and symmetric matrix, whose 
%               diagonal elements are all 1s. r_2(i,j) stores the minimal
%               second order partial correlation between feature i and 
%               feature j given all possible combinations of the remaining 
%               features.
%     'p_2'     An n \times n numeric matrix, storing the p-values
%               associated with corresponding elements in r_2. Since r_2 
%               is symmetric, p_2 only calculates the p-values for the 
%               upper triangle of r_2, all other elements are 1s.
%    
% Note, lopc returns p-values for upper-triangular eleemnts in the data
% matrix. Adjusted p-values based FDR can be easily obtained using the
% falseDiscoveryRate function provided by me.
%
% References:
% [1] Zuo, Yiming, Guoqiang Yu, Mahlet G. Tadesse, and Habtom W. Ressom. 
%     Biological network inference using low order partial correlation. 
%     Methods 69, no. 3 (2014): 266-273.
%
% Copyright 2014-2016, Yiming Zuo.

%% 0th order partial correlation (i.e., Pearson's correlation)
m = size(data, 1); % sample size
n = size(data, 2); % variable number
[r, p] = corrcoef(data);
p_a = falseDiscoveryRate(n, p); % adjusted p-value based on FDR

%% First-order partial correlation
r_1 = ones(n, n); % initialize first order partial correlation matrix 
r1 = ones(n, n, n);
for i = 1:n-1
    for j = i+1:n
        if (p_a(i, j) < thres) % skip the computation if 0-th order partial 
% correlation is already insignificant
            for k = 1:n
                if (k ~= i) && (k ~= j)
                    r1(i,j,k) = (r(i,j)-r(i,k)*r(j,k)) / sqrt((1-r(i,k)^2)*(1-r(j,k)^2));
                    r1(j, i, k) = r1(i, j, k);
                    if abs(r1(i, j, k)) < abs(r_1(i, j)) % only store the minimal value
                        r_1(i, j) = r1(i, j, k);
                    end
                end
            end
        else
            r_1(i, j) = 0;
        end
        r_1(j, i) = r_1(i, j);
    end
end

% Fisher transformation for first-order partial correlation
z_1 = 0.5 * log((1+r_1)./(1-r_1));
p_1 = ones(n,n); % initialize p-value matrix
for i = 1:n-1
    for j=i+1:n % only store the upper-triangular elements
        p_1(i,j) = (1 - normcdf(abs(z_1(i,j)),0,sqrt(1/(m-3-1)))) * 2;
    end
end
p_1_a = falseDiscoveryRate(n, p_1); % adjusted p-value based on FDR

%% second-order correlation
r_2 = ones(n, n);
r2 = 1;
for i = 1:n-1
    for j = i+1:n
        if (p_a(i,j)<thres && p_1_a(i,j)<thres) % skip the computation if 0-th order partial 
% correlation or first-order partial correlation is already insignificant
            for k = 1:n-1
                if (k~=i) && (k~=j)
                    for q = k+1:n
                        if (q~=i) && (q~=j)
                            r2 = (r1(i,j,k)-r1(i,q,k)*r1(j,q,k))/sqrt((1-r1(i,q,k)^2)*(1-r1(j,q,k)^2));
                            if abs(r2) < abs(r_2(i,j))
                                r_2(i,j) = r2;
                            end
                        end
                    end
                end
            end
        else
            r_2(i, j) = 0;
        end
        r_2(j, i) = r_2(i, j);
    end
end

% Fisher transformation for second-order partial correlation
z_2 = 0.5 * log((1+r_2)./(1-r_2));
p_2 = ones(n, n);
for i = 1:n-1
    for j = i+1:n
        p_2(i,j) = (1 - normcdf(abs(z_2(i,j)), 0, sqrt(1/(m-3-2)))) * 2;
    end
end

