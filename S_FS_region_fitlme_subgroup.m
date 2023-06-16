% group*wave interaction: regional_vol~Group*Wave+Sex+Age+global_vol+(1|SubID)
% separate for medicated and unmedicated patients at recruitment

%% 1.Volume

clear;clc; close all;

load('R1_Subject_Info.mat','SubInfo_merge_FYF*','excl*')
load('R7_freesurfer_measure.mat')

tbl = SubInfo_merge_FYF_updateAgeSex;
tbl.SubID = categorical(tbl.SubID);

sub=1:180;
sub_medicated = ones(length(sub),1);
sub_unmedicated = ones(length(sub),1);

sub_medicated(1:78) = SubInfo_merge_FYF_updateAgeSex.Days_at_recruitment(1:78)~=0;
sub_unmedicated(1:78) = SubInfo_merge_FYF_updateAgeSex.Days_at_recruitment(1:78)==0;
[nnz(sub_medicated(1:78)),nnz(sub_unmedicated(1:78))]/2 % 32     7
sub_medicated = logical(sub_medicated); 
sub_unmedicated = logical(sub_unmedicated); 

%%%% excl one sub with age>med+3std %%%%%
global_measure = global_vol(sub);
global_measure(excl_age) = nan;

%%%% freesurfer output %%%%%
ROI= 2:35;
aseg_roi=[6:9,13:14,16,24:30];  
measure_vol = [aparc_volume_lh_match(sub,ROI),aparc_volume_rh_match(sub,ROI),aseg_stats_match(sub,aseg_roi)];

[p_fitlme_grpwave_vol_medicated,t_fitlme_grpwave_vol_medicated, coef_fitlme_grpwave_vol_medicated, measure_reg_vol_medicated]=...
    S_func_fitlme_grpwave(tbl(sub_medicated,:), measure_vol(sub_medicated,:), global_measure(sub_medicated));

[p_fitlme_grpwave_vol_unmedicated,t_fitlme_grpwave_vol_unmedicated, coef_fitlme_grpwave_vol_unmedicated, measure_reg_vol_unmedicated]=...
    S_func_fitlme_grpwave(tbl(sub_unmedicated,:), measure_vol(sub_unmedicated,:), global_measure(sub_unmedicated));


% save('R9_freesurfer_regional_fitlme_updateAgeSex_subgroup.mat','*medicated');

%% 2.Cortical thickness
clear;clc; close all;
load('R1_Subject_Info.mat','SubInfo_merge_FYF*','excl*')
load('R7_freesurfer_measure.mat')

tbl = SubInfo_merge_FYF_updateAgeSex;
tbl.SubID = categorical(tbl.SubID);

sub=1:180;
sub_medicated = ones(length(sub),1);
sub_unmedicated = ones(length(sub),1);

sub_medicated(1:78) = SubInfo_merge_FYF_updateAgeSex.Days_at_recruitment(1:78)~=0;
sub_unmedicated(1:78) = SubInfo_merge_FYF_updateAgeSex.Days_at_recruitment(1:78)==0;
[nnz(sub_medicated(1:78)),nnz(sub_unmedicated(1:78))]/2 % 32     7
sub_medicated = logical(sub_medicated); 
sub_unmedicated = logical(sub_unmedicated); 

%%%% excl one sub with age>med+3std %%%%%
global_measure = global_thick(sub);
global_measure(excl_age) = nan;

%%%% freesurfer output %%%%%
ROI= 2:35;
measure_thick = [aparc_thickness_lh_match(sub,ROI),aparc_thickness_rh_match(sub,ROI)];


[p_fitlme_grpwave_thick_medicated,t_fitlme_grpwave_thick_medicated, coef_fitlme_grpwave_thick_medicated, measure_reg_thick_medicated]=...
    S_func_fitlme_grpwave(tbl(sub_medicated,:), measure_thick(sub_medicated,:), global_measure(sub_medicated));

[p_fitlme_grpwave_thick_unmedicated,t_fitlme_grpwave_thick_unmedicated, coef_fitlme_grpwave_thick_unmedicated, measure_reg_thick_unmedicated]=...
    S_func_fitlme_grpwave(tbl(sub_unmedicated,:), measure_thick(sub_unmedicated,:), global_measure(sub_unmedicated));


% save('R9_freesurfer_regional_fitlme_updateAgeSex_subgroup.mat','*medicated','-append');


%% 3.Surface area
clear;clc; close all;
load('R1_Subject_Info.mat','SubInfo_merge_FYF*','excl*')
load('R7_freesurfer_measure.mat')

tbl = SubInfo_merge_FYF_updateAgeSex;
tbl.SubID = categorical(tbl.SubID);

sub=1:180;
sub_medicated = ones(length(sub),1);
sub_unmedicated = ones(length(sub),1);

sub_medicated(1:78) = SubInfo_merge_FYF_updateAgeSex.Days_at_recruitment(1:78)~=0;
sub_unmedicated(1:78) = SubInfo_merge_FYF_updateAgeSex.Days_at_recruitment(1:78)==0;
[nnz(sub_medicated(1:78)),nnz(sub_unmedicated(1:78))]/2 % 32     7
sub_medicated = logical(sub_medicated); 
sub_unmedicated = logical(sub_unmedicated); 

%%%% excl one sub with age>med+3std %%%%%
global_measure = global_area(sub);
global_measure(excl_age) = nan;

%%%% freesurfer output %%%%%
ROI= 2:35;
measure_area = [aparc_area_lh_match(sub,ROI),aparc_area_rh_match(sub,ROI)];

[p_fitlme_grpwave_area_medicated,t_fitlme_grpwave_area_medicated, coef_fitlme_grpwave_area_medicated, measure_reg_area_medicated]=...
    S_func_fitlme_grpwave(tbl(sub_medicated,:), measure_area(sub_medicated,:), global_measure(sub_medicated));

[p_fitlme_grpwave_area_unmedicated,t_fitlme_grpwave_area_unmedicated, coef_fitlme_grpwave_area_unmedicated, measure_reg_area_unmedicated]=...
    S_func_fitlme_grpwave(tbl(sub_unmedicated,:), measure_area(sub_unmedicated,:), global_measure(sub_unmedicated));

% save('R9_freesurfer_regional_fitlme_updateAgeSex_subgroup.mat','*medicated','-append');


%% 4.FDR for all
clear;clc; close all;

load('R9_freesurfer_regional_fitlme_updateAgeSex_subgroup.mat');
load('R9_freesurfer_regional_fitlme_updateAgeSex.mat','*_all','h_grp*');

p_all_medicated = [p_fitlme_grpwave_vol_medicated(69:end,:); p_fitlme_grpwave_thick_medicated; p_fitlme_grpwave_area_medicated];
t_all_medicated = [t_fitlme_grpwave_vol_medicated(69:end,:); t_fitlme_grpwave_thick_medicated; t_fitlme_grpwave_area_medicated];

p_all_unmedicated = [p_fitlme_grpwave_vol_unmedicated(69:end,:); p_fitlme_grpwave_thick_unmedicated; p_fitlme_grpwave_area_unmedicated];
t_all_unmedicated = [t_fitlme_grpwave_vol_unmedicated(69:end,:); t_fitlme_grpwave_thick_unmedicated; t_fitlme_grpwave_area_unmedicated];

h_grp =fdr_bh(p_all.Group,0.05);        
h_grpwave =fdr_bh(p_all.GroupWave,0.05);

tbl_grp = [t_all_medicated.Group(h_grp),p_all_medicated.Group(h_grp),...
    t_all_medicated.Wave(h_grp),p_all_medicated.Wave(h_grp),...
    t_all_medicated.GroupWave(h_grp),p_all_medicated.GroupWave(h_grp)];

tbl_grpwave = [t_all_medicated.Group(h_grpwave),p_all_medicated.Group(h_grpwave),...
    t_all_medicated.Wave(h_grpwave),p_all_medicated.Wave(h_grpwave),...
    t_all_medicated.GroupWave(h_grpwave),p_all_medicated.GroupWave(h_grpwave)];

% tbl_grp = [t_all_unmedicated.Group(h_grp),p_all_unmedicated.Group(h_grp),...
%     t_all_unmedicated.Wave(h_grp),p_all_unmedicated.Wave(h_grp),...
%     t_all_unmedicated.GroupWave(h_grp),p_all_unmedicated.GroupWave(h_grp)];
% 
% tbl_grpwave = [t_all_unmedicated.Group(h_grpwave),p_all_unmedicated.Group(h_grpwave),...
%     t_all_unmedicated.Wave(h_grpwave),p_all_unmedicated.Wave(h_grpwave),...
%     t_all_unmedicated.GroupWave(h_grpwave),p_all_unmedicated.GroupWave(h_grpwave)];


