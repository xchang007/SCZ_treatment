
%% 1.Ratio before and after treatment
clear;clc; close all;

load('R1_Subject_Info.mat','SubInfo_merge_FYF*')
load('R9_freesurfer_regional_fitlme_updateAgeSex.mat');

SZ_BL = find((SubInfo_merge_FYF.Group=='SZ').*(SubInfo_merge_FYF.Wave=='BL'));
SZ_FU1 = find((SubInfo_merge_FYF.Group=='SZ').*(SubInfo_merge_FYF.Wave=='FU1'));
HC_BL = find((SubInfo_merge_FYF.Group=='HC').*(SubInfo_merge_FYF.Wave=='BL'));
HC_FU1 = find((SubInfo_merge_FYF.Group=='HC').*(SubInfo_merge_FYF.Wave=='FU1'));

delta_vol_SZ = (table2array(measure_vol(SZ_FU1,:))-table2array(measure_vol(SZ_BL,:)))./table2array(measure_vol(SZ_BL,:));
delta_vol_HC = (table2array(measure_vol(HC_FU1,:))-table2array(measure_vol(HC_BL,:)))./table2array(measure_vol(HC_BL,:));

delta_thick_SZ = (table2array(measure_thick(SZ_FU1,:))-table2array(measure_thick(SZ_BL,:)))./table2array(measure_thick(SZ_BL,:));
delta_thick_HC = (table2array(measure_thick(HC_FU1,:))-table2array(measure_thick(HC_BL,:)))./table2array(measure_thick(HC_BL,:));

delta_area_SZ = (table2array(measure_area(SZ_FU1,:))-table2array(measure_area(SZ_BL,:)))./table2array(measure_area(SZ_BL,:));
delta_area_HC = (table2array(measure_area(HC_FU1,:))-table2array(measure_area(HC_BL,:)))./table2array(measure_area(HC_BL,:));


% save('R9_freesurfer_regional_prepost.mat','delta_*')

%% 2.Brain changes correlated with gene set of interests (GOI)
% S22_gene_expression_GAMBA_GOI.m  %2

clear;clc; close all;

addpath(genpath('Project_2023_Cui_treatment/Toolbox_Fudan'))
load('Literature/Genetic/GOI/R1_GOI.mat','GE_*', 'genes_*','region*')

load('R1_Subject_Info.mat','SubInfo_merge_FYF*','excl*')
load('R9_freesurfer_regional_prepost.mat');

GE_lists = [GE_2018_lancet_Antipsychotic, GE_2011_plos_WBC, GE_2014_HMG_WBC];
GE_lists_name = {'18Antipsychotic','11WBC','11Neutrophil','11Basophil',...
    '11Monocyte','11Lymphocytes','14WBC','14Neutrophil','14Monocyte'};
GG = GE_lists./nanstd(GE_lists);



%%%%% gene expression & delta_vol %%%%%
% delta_vol_SZ(excl_age(1),:)=nan;
% dataIMG_SZ = nanmean(delta_vol_SZ);
% dataIMG_HC = nanmean(delta_vol_HC);
% dataIMG_SZ = [dataIMG_SZ(69:end),dataIMG_SZ(1:68)]';
% dataIMG_HC = [dataIMG_HC(69:end),dataIMG_HC(1:68)]';
% dataIMG_SZ = dataIMG_SZ./nanstd(dataIMG_SZ);
% dataIMG_HC = dataIMG_HC./nanstd(dataIMG_HC);

% [r_GG_delta_vol, p_GG_delta_vol] = corr(dataIMG_SZ,GG,'rows','complete')


%%%%% gene expression & delta_thick %%%%%
delta_thick_SZ(excl_age(1),:)=nan;
dataIMG_SZ = nanmean(delta_thick_SZ)';
dataIMG_HC = nanmean(delta_thick_HC)';
dataIMG_SZ = dataIMG_SZ./nanstd(dataIMG_SZ);
dataIMG_HC = dataIMG_HC./nanstd(dataIMG_HC);

[r_GG_delta_thick, p_GG_delta_thick] = corr([dataIMG_SZ,dataIMG_HC],GG(15:end,:),'rows','complete');


%%%%% gene expression & delta_area %%%%%
delta_area_SZ(excl_age(1),:)=nan;
dataIMG_SZ = nanmean(delta_area_SZ)';
dataIMG_HC = nanmean(delta_area_HC)';
dataIMG_SZ = dataIMG_SZ./nanstd(dataIMG_SZ);
dataIMG_HC = dataIMG_HC./nanstd(dataIMG_HC);

[r_GG_delta_area, p_GG_delta_area] = corr([dataIMG_SZ,dataIMG_HC],GG(15:end,:),'rows','complete');


%%%%% FDR for p_GG_delta_thick p_GG_delta_area %%%%%
p_all = [p_GG_delta_thick,p_GG_delta_area];
[h, crit_p]=fdr_bh(p_all,0.05)          


%% 3.Brain changes correlated with gene set of interests (GOI), subgroup analysis
% S22_gene_expression_GAMBA_GOI.m  %2

clear;clc; close all;

addpath(genpath('Project_2023_Cui_treatment/Toolbox_Fudan'))
load('Literature/Genetic/GOI/R1_GOI.mat','GE_*', 'genes_*','region*')


load('R1_Subject_Info.mat','SubInfo_merge_FYF*','excl*')
load('R9_freesurfer_regional_prepost.mat');
load('R9_freesurfer_regional_fitlme_updateAgeSex_subgroup.mat', 'sub_*') % for subgroup analysis

GE_lists = [GE_2018_lancet_Antipsychotic, GE_2011_plos_WBC, GE_2014_HMG_WBC];
GE_lists_name = {'18Antipsychotic','11WBC','11Neutrophil','11Basophil',...
    '11Monocyte','11Lymphocytes','14WBC','14Neutrophil','14Monocyte'};
GG = GE_lists./nanstd(GE_lists);


%%%%% gene expression & delta_thick %%%%%
delta_thick_SZ(excl_age(1),:)=nan;

dataIMG_SZ_medicated   = nanmean(delta_thick_SZ(sub_medicated(1:39),:))';
dataIMG_SZ_unmedicated = nanmean(delta_thick_SZ(sub_unmedicated(1:39),:))';
dataIMG_SZ_medicated   = dataIMG_SZ_medicated./nanstd(dataIMG_SZ_medicated);
dataIMG_SZ_unmedicated = dataIMG_SZ_unmedicated./nanstd(dataIMG_SZ_unmedicated);

[r_GG_delta_thick_subgroup, p_GG_delta_thick_subgroup] = corr([dataIMG_SZ_medicated,dataIMG_SZ_unmedicated],GG(15:end,:),'rows','complete');


%%%%% gene expression & delta_area %%%%%
delta_area_SZ(excl_age(1),:)=nan;

dataIMG_SZ_medicated   = nanmean(delta_area_SZ(sub_medicated(1:39),:))';
dataIMG_SZ_unmedicated = nanmean(delta_area_SZ(sub_unmedicated(1:39),:))';
dataIMG_SZ_medicated   = dataIMG_SZ_medicated./nanstd(dataIMG_SZ_medicated);
dataIMG_SZ_unmedicated = dataIMG_SZ_unmedicated./nanstd(dataIMG_SZ_unmedicated);

[r_GG_delta_area_subgroup, p_GG_delta_area_subgroup] = corr([dataIMG_SZ_medicated,dataIMG_SZ_unmedicated],GG(15:end,:),'rows','complete');


