% Yongbin Wei: use brainsmash to generate 5000 surrogate maps
% https://brainsmash.readthedocs.io/en/latest/
% run_brainsmash.py

%% Spin test for significant correlation

clear; clc; close all;

%%%%% Load surrogated maps from Yongbin %%%%%

surrogate_vol_cort_sz = readmatrix('R22_GAMBA_Antipsychotic_Yongbin_surrogate_map_5000_2022.06.28\surrogate_vol_sz.txt');
surrogate_vol_subcort_sz = readmatrix('R22_GAMBA_Antipsychotic_Yongbin_surrogate_map_5000_2022.06.28\surrogate_vol_subcort_SZ.txt');
surrogate_thick_sz = readmatrix('R22_GAMBA_Antipsychotic_Yongbin_surrogate_map_5000_2022.06.28\surrogate_thick_sz.txt');
surrogate_area_sz = readmatrix('R22_GAMBA_Antipsychotic_Yongbin_surrogate_map_5000_2022.06.28\surrogate_area_sz.txt');

surrogate_vol_sz = [surrogate_vol_subcort_sz, surrogate_vol_cort_sz]; % 1-14 subcortical, 15-82 cortical


%%%%% Load gene expression of interests %%%%%

load('Literature\Genetic\GOI\R1_GOI.mat','*_WBC*')
GE_GOI_mean = GE_2011_plos_WBC(:,4)/nanstd(GE_2011_plos_WBC(:,4));   % mean expression of 2011_plos_WBC Monocyte 
% GE_GOI_mean = GE_2014_HMG_WBC(:,3)/nanstd(GE_2014_HMG_WBC(:,3));   % mean expression of 2014_HMG_WBC Monocyte



%%%%% gene expression & delta_thick %%%%%
load('R1_Subject_Info.mat')
load('R9_freesurfer_regional_prepost.mat','delta_thick_*','region_label')

delta_thick_SZ(excl_age(1),:)=nan;
dataIMG_thick_SZ = nanmean(delta_thick_SZ)';
dataIMG_thick_SZ = dataIMG_thick_SZ./nanstd(dataIMG_thick_SZ);

reg = regstats(dataIMG_thick_SZ, GE_GOI_mean(15:end), 'linear', 'tstat');
beta_thick_sz = reg.tstat.beta(2);
pval_thick_sz = reg.tstat.pval(2); % same with corr p_GG_delta_vol

% get surrogated maps
Nmaps = size(surrogate_thick_sz,1);
beta_thick_sz_surrogate = zeros(Nmaps,1);

for i=1:Nmaps
    
    tmp_thick_sz = surrogate_thick_sz(i,:)';
    reg = regstats(tmp_thick_sz, GE_GOI_mean(15:end), 'linear', 'tstat');
    beta_thick_sz_surrogate(i) = reg.tstat.beta(2);

end

% p spin test value
p_spin_thick_sz = nnz(beta_thick_sz_surrogate> beta_thick_sz)/Nmaps        




