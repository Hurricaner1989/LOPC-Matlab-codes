%% Start clean
clc; clear; close all;
%% Read in proteomic data
% make sure pwd is LOPC matlab package
pwd 
% set the working directory to where the data is
cd 'Data'
% read in protein names
filename = 'ProteinRowNames.csv';
[info, title] = xlsread(filename);
% preprocess for proteomic data
[peakProtein, idProtein] = mergeProtein();
% normalization
peakProtein_n = [];
peakProtein_n = (peakProtein - repmat(mean(peakProtein),size(peakProtein,1),1)) ...
    ./repmat(std(peakProtein),size(peakProtein,1),1);
% go back one level to the package folder
cd '..'
pwd
%% 0-th partial correlation
n = size(peakProtein_n,1); % variable number
m = size(peakProtein_n,2);% sample size
[r, p] = corrcoef(peakProtein_n');
p_a = falseDiscoveryRate(n, p); % obtain adjusted p-value based on FDR
%% low order partial correlation
% first-order partial correlation
[r_1, p_1] = firstParCorr(n, r, m);
p_1_a = falseDiscoveryRate(n, p_1);
% second-order partial correlation
[r_2, p_2] = secParCorr(n, r, m);
p_2_a = falseDiscoveryRate(n, p_2); % some warning comes out but we can safely 
% proceed without considering them right now
%% LOPC
thres = 0.1; % the cutoff for adjusted p-value, 0.05, 0.1 are recommended
[re_2, pe_2] = lopc(peakProtein_n', thres);
pe_2_a = falseDiscoveryRate(n, pe_2);
%% Find significant LOPC pairs after correcting multiple testing problem.
[i, j] = find(pe_2_a < thres);
[i, j] % display their (row,col) indices.

% show the network
coords = [cos(2*pi*(1:n)/n); sin(2*pi*(1:n)/n)]';
gplot(pe_2_a < thres, coords, '-*')
text(coords(:,1), coords(:,2), num2str((1:n)'), 'FontSize', 10)

% % save connections to excel and later import to Cytostape 
% xlswrite('C:\Users\yz335\Desktop\ProteinNet.xlsx',i,1,'A1')
% xlswrite('C:\Users\yz335\Desktop\ProteinNet.xlsx',j,1,'B1')
% xlswrite('C:\Users\yz335\Desktop\ProteinNet.xlsx',re_2(find(pe_2_a<0.1)),1,'C1')






