% group*wave interaction: regional_vol~Group*Wave+Sex+Age+global_vol+(1|SubID)


%% 1.Volume

clear;clc; close all;

load('R1_Subject_Info.mat','SubInfo_merge_FYF*','excl*')
load('R7_freesurfer_measure.mat')

sub = [1:180];
tbl = SubInfo_merge_FYF_updateAgeSex;

tbl.SubID = categorical(tbl.SubID);

%%%% excl one sub with age>med+3std %%%%%
global_measure = global_vol(sub);
global_measure(excl_age) = nan;

%%%% freesurfer output %%%%%
ROI= 2:35;
aseg_roi=[6:9,13:14,16,24:30];  
measure_vol = [aparc_volume_lh_match(sub,ROI),aparc_volume_rh_match(sub,ROI),aseg_stats_match(sub,aseg_roi)];
region_label_vol=[aparc_volume_lh_match.Properties.VariableNames(ROI),...
    aparc_volume_rh_match.Properties.VariableNames(ROI),...
    aseg_stats_match.Properties.VariableNames(aseg_roi)]'; 

[p_fitlme_grpwave_vol,t_fitlme_grpwave_vol, coef_fitlme_grpwave_vol, measure_reg_vol]=...
    S_func_fitlme_grpwave(tbl, measure_vol, global_measure);


% save('R9_freesurfer_regional_fitlme_updateAgeSex.mat','*_vol','region_label*');

%% 2.Cortical thickness
clear;clc; close all;
load('R1_Subject_Info.mat','SubInfo_merge_FYF*','excl*')
load('R7_freesurfer_measure.mat')

sub = [1:180];
tbl = SubInfo_merge_FYF_updateAgeSex;

tbl.SubID = categorical(tbl.SubID);

%%%% excl one sub with age>med+3std %%%%%
global_measure = global_thick(sub);
global_measure(excl_age) = nan;

%%%% freesurfer output %%%%%
ROI= 2:35;
measure_thick = [aparc_thickness_lh_match(sub,ROI),aparc_thickness_rh_match(sub,ROI)];
region_label_thick=[aparc_thickness_lh_match.Properties.VariableNames(ROI),...
    aparc_thickness_rh_match.Properties.VariableNames(ROI)]'; 


[p_fitlme_grpwave_thick,t_fitlme_grpwave_thick, coef_fitlme_grpwave_thick, measure_reg_thick]=...
    S_func_fitlme_grpwave(tbl, measure_thick, global_measure);


% save('R9_freesurfer_regional_fitlme_updateAgeSex.mat','*_thick','-append');


%% 3.Surface area
clear;clc; close all;
load('R1_Subject_Info.mat','SubInfo_merge_FYF*','excl*')
load('R7_freesurfer_measure.mat')

sub = [1:180];
tbl = SubInfo_merge_FYF_updateAgeSex;

tbl.SubID = categorical(tbl.SubID);

%%%% excl one sub with age>med+3std %%%%%
global_measure = global_area(sub);
global_measure(excl_age) = nan;

%%%% freesurfer output %%%%%
ROI= 2:35;
measure_area = [aparc_area_lh_match(sub,ROI),aparc_area_rh_match(sub,ROI)];
region_label_area=[aparc_area_lh_match.Properties.VariableNames(ROI),...
    aparc_area_rh_match.Properties.VariableNames(ROI)]'; 


[p_fitlme_grpwave_area,t_fitlme_grpwave_area, coef_fitlme_grpwave_area, measure_reg_area]=...
    S_func_fitlme_grpwave(tbl, measure_area, global_measure);



% save('R9_freesurfer_regional_fitlme_updateAgeSex.mat','*_area','-append');


%% 4.FDR for all
clear;clc; close all;
load('R9_freesurfer_regional_fitlme_updateAgeSex.mat');

p_all = [p_fitlme_grpwave_vol(69:end,:); p_fitlme_grpwave_thick; p_fitlme_grpwave_area];
t_all = [t_fitlme_grpwave_vol(69:end,:); t_fitlme_grpwave_thick; t_fitlme_grpwave_area];
region_label_all = [region_label_vol(69:end); region_label_thick; region_label_area];
measure_all = [measure_vol(:,69:end), measure_thick, measure_area];


h_grp =fdr_bh(p_all.Group,0.05);        
h_grpwave =fdr_bh(p_all.GroupWave,0.05);

region_sig = region_label_all(h_grp)
tbl = [t_all.Group(h_grp),p_all.Group(h_grp), t_all.Wave(h_grp),p_all.Wave(h_grp), t_all.GroupWave(h_grp),p_all.GroupWave(h_grp)];

% region_sig = region_label_all(h_grpwave)
% tbl = [t_all.Group(h_grpwave),p_all.Group(h_grpwave), t_all.Wave(h_grpwave),p_all.Wave(h_grpwave), t_all.GroupWave(h_grpwave),p_all.GroupWave(h_grpwave)];

% save('R9_freesurfer_regional_fitlme_updateAgeSex.mat','*_all','h_grp*','-append');

