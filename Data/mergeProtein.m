function [commonpeak_intensity, commonpeak_id] = mergeProtein()
filename = 'Proteomics_LCMS_TU_refined140518.xlsx';
[info, title] = xlsread(filename);
commonpeak_intensity = info(:, 1:89);
commonpeak_id = title(1, 1:89);