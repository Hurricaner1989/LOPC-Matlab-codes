function [r_2, p_2] = secParCorr(n, r, m)
% Compute second order partial correlation
% Inputs:
% 'n'    number of variables
% 'r'    Pearson's correlation matrix
% 'm'    number of samples
% Outputs:
% 'r_2'  second order partial correlation matrix
% 'p_2'  second order partial correlation p-value matrix 

% first-order partial correlation
r1 = ones(n, n, n); % r1 is by default set to 1
for i = 1:n-1
    for j = i+1:n
        for k = 1:n
            if (k~=i) && (k~=j)
                r1(i,j,k) = (r(i,j)-r(i,k)*r(j,k)) / sqrt((1-r(i,k)^2)*(1-r(j,k)^2));
                r1(j, i, k) = r1(i, j, k);
            end
        end
    end
end

% second-order partial correlation
r_2 = ones(n,n);
for i = 1:n-1
    for j = i+1:n
        for k = 1:n-1
            if (k~=i) && (k~=j)
                for q = k+1:n
                    if (q~=i) && (q~=j)
                        r2 = (r1(i,j,k)-r1(i,q,k)*r1(j,q,k))/sqrt((1-r1(i,q,k)^2)*(1-r1(j,q,k)^2));
                        if abs(r2) < abs(r_2(i, j))
                            r_2(i, j) = r2;
                        end
                    end
                end
            end
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
