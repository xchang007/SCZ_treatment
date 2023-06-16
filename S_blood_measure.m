
%% 1. white blood cell count 
% check blood cell count distribution
% swtest normalitytest package 
% https://www.mathworks.com/matlabcentral/fileexchange/13964-shapiro-wilk-and-shapiro-francia-normality-tests

blood_name = {'Leukocyte (10⁹/L)','Neutrophils %','Lymphocyte %','Monocyte %','Eosinophils %','Basophils %'};
tmp_blood = table2array(blood_measure_0628(:,3:8));
figure;
for i=1:6

    [H(i), pValue(i), SWstatistic(i)] = swtest(tmp_blood(:,i), 0.05);
   
    subplot(2,3,i);histogram(tmp_blood(:,i),8,'FaceColor',[33,113,181]/255,'EdgeColor','w');
    title(blood_name(i));
end

varfun(@nanmedian,blood_measure_0628(:,3:8))
%         5.645                   0.546                       0.362                       0.066                   0.017                        0.002          
varfun(@(X) prctile(X,25),blood_measure_0628(:,3:8))
%         4.79                0.469                   0.298                 0.055                 0.011                  0.002       
varfun(@(X) prctile(X,75),blood_measure_0628(:,3:8))
%        6.58                 0.63                   0.441                 0.072                 0.026                  0.004       


%% 2.Correlation between blood cell with PANSS, WAIS
clc;close all;
load('R1_Subject_Info.mat')

%%%%% Corr with PANSS %%%%%
[r_blood_PANSS_BL_Spearman,p_blood_PANSS_BL_Spearman]=corr([SCZ_BL.PANSS_pos,SCZ_BL.PANSS_neg,SCZ_BL.PANSS_gen,...
    table2array(blood_measure_0628(:,3:8))],'rows','complete','type','Spearman');

[r_blood_PANSS_FU_Spearman,p_blood_PANSS_FU_Spearman]=corr([SCZ_FU.PANSS_pos,SCZ_FU.PANSS_neg,SCZ_FU.PANSS_gen,...
    table2array(blood_measure_0628(:,3:8))],'rows','complete','type','Spearman');

%%%%% Corr with WAIS %%%%%
[r_blood_WAIS_BL_Spearman,p_blood_WAIS_BL_Spearman]=corr([SCZ_BL.WAIS_DigitSymbol,SCZ_BL.WAIS_DigitSpan_forward,SCZ_BL.WAIS_DigitSpan_backward,...
    table2array(blood_measure_0628(:,3:8))],'rows','complete','type','Spearman');

[r_blood_WAIS_FU_Spearman,p_blood_WAIS_FU_Spearman]=corr([SCZ_FU.WAIS_DigitSymbol,SCZ_FU.WAIS_DigitSpan_forward,SCZ_FU.WAIS_DigitSpan_backward,...
    table2array(blood_measure_0628(:,3:8))],'rows','complete','type','Spearman');

[h, crit_p]=fdr_bh([p_blood_PANSS_BL_Spearman,p_blood_PANSS_FU_Spearman,p_blood_WAIS_BL_Spearman,p_blood_WAIS_FU_Spearman],0.05);   


%% 3.Correlation between blood cell with FS
clc;clear;close all;
load('R9_freesurfer_regional_prepost.mat','delta_*')

index = 6; % Monocyte 
[r_region_vol_inflammation_spearman,p_region_vol_inflammation_spearman] = corr(delta_vol_SZ,table2array(blood_measure_0628(:,index)),'Rows','complete','Type','Spearman');
[r_region_thick_inflammation_spearman,p_region_thick_inflammation_spearman] = corr(delta_thick_SZ,table2array(blood_measure_0628(:,index)), 'Rows','complete','Type','Spearman');
[r_region_area_inflammation_spearman,p_region_area_inflammation_spearman] = corr(delta_area_SZ, table2array(blood_measure_0628(:,index)), 'Rows','complete','Type','Spearman');


load('R9_freesurfer_regional_fitlme_updateAgeSex.mat');
h = or(h_grp,h_grpwave);

p_region_inflammation_spearman_sig = [p_region_vol_inflammation_spearman(h,:);p_region_thick_inflammation_spearman(h,:); p_region_area_inflammation_spearman(h,:)];
r_region_inflammation_spearman_sig = [r_region_vol_inflammation_spearman(h,:);r_region_thick_inflammation_spearman(h,:); r_region_area_inflammation_spearman(h,:)];

[h, crit_p]=fdr_bh(p_region_inflammation_spearman_sig,0.05); 

