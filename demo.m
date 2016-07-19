% A demo to illustrate the use of firstParCorr, secParCorr and lopc
% functions using a real proteomic dataset

%% Start clean
clc; clear; close all;

%% Read in proteomic data 
% make sure pwd is LOPC matlab package
pwd 
% set the working directory to where the data is
cd 'Data'
% prepare proteomic data
[peakProtein, ~] = readProtein();
% normalization
peakProtein_n = [];
peakProtein_n = (peakProtein - repmat(mean(peakProtein,2),1,size(peakProtein,2))) ...
    ./repmat(std(peakProtein,0,2),1,size(peakProtein,2));
% go back one level up to the package folder
cd '..'
% make sure pwd is LOPC matlab package
pwd

%% 0-th partial correlation (Pearson's correlation)
n = size(peakProtein_n, 1); % variable number
m = size(peakProtein_n, 2); % sample size
[r, p] = corrcoef(peakProtein_n'); % Pearson's correlation
p_a = falseDiscoveryRate(n, p); % obtain adjusted p-value based on FDR

%% low order partial correlation
% this section lists traditional methods to calculate low order partial
% correlation up to first and second order

% first-order partial correlation
[r_1, p_1] = firstParCorr(n, r, m);
p_1_a = falseDiscoveryRate(n, p_1);
% second-order partial correlation
[r_2, p_2] = secParCorr(n, r, m);
p_2_a = falseDiscoveryRate(n, p_2); 
% some warning comes out but we can safely proceed without considering them 

%% LOPC
% a more efficient algorithm to calculate low order partial correlation up
% to the second order
thres = 0.1; % the cutoff for adjusted p-value, 0.01, 0.05 or 0.1 are recommended
[re_2, pe_2] = lopc(peakProtein_n', thres);
pe_2_a = falseDiscoveryRate(n, pe_2);

%% Find significant protein pairs
% LOPC is used here, but one can easily replace LOPC with other methods
[i, j] = find(pe_2_a < thres);
[i, j] % display their (row,col) indices.

% show the network
coords = [cos(2*pi*(1:n)/n); sin(2*pi*(1:n)/n)]';
gplot(pe_2_a < thres, coords, '-*')
text(coords(:,1), coords(:,2), num2str((1:n)'), 'FontSize', 10)

% % save significant pairs to excel to import into Cytostape for better visulization 
% xlswrite('C:\Users\yz335\Desktop\ProteinNet.xlsx',i,1,'A1')
% xlswrite('C:\Users\yz335\Desktop\ProteinNet.xlsx',j,1,'B1')
% xlswrite('C:\Users\yz335\Desktop\ProteinNet.xlsx',re_2(find(pe_2_a<0.1)),1,'C1')

