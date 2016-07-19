function [commonpeak_intensity, commonpeak_id] = readProtein()
filename = 'Proteomics_LCMS_TU_refined140518.xlsx';
[info, title] = xlsread(filename);
commonpeak_intensity = info(:, 1:89); % read in protein intensity
commonpeak_id = title(1, 2:90);