function [p_fitlme_grpwave,t_fitlme_grpwave, coef_fitlme_grpwave, measure_reg]=...
    S_func_fitlme_grpwave(tbl, measure, global_measure)

% tbl              variables: Group,Wave,Sex,Age,SubID
% measure          freesurfer regional measurement, Nsub*Nregion, table
% global_measure   global measurement to be regressed

p_fitlme_grpwave = zeros(size(measure,2),6);                   % p: group, wave, sex, age, global, group*wave
t_fitlme_grpwave = zeros(size(measure,2),6);                   % t: group, wave, sex, age, global, group*wave
coef_fitlme_grpwave = zeros(size(measure,2),6);
measure_reg = zeros(size(measure));                                 % regressed by age,sex

for i= 1:size(measure,2)
    i
    tbl.global = global_measure;
    tbl.regional = table2array(measure(:,i));
    tmp_lme = fitlme(tbl,'regional~Group*Wave+Sex+Age+global+(1|SubID)');
    p_fitlme_grpwave(i,:) = double(tmp_lme.Coefficients(2:end,6));   % group, wave, sex, age, global_vol, group*wave
    t_fitlme_grpwave(i,:) = double(tmp_lme.Coefficients(2:end,4));   % group, wave, sex, age, global_vol, group*wave
    coef_fitlme_grpwave(i,:) = double(tmp_lme.Coefficients(2:end,2));% group, wave, sex, age, global_vol, group*wave

    % regress age,sex
    measure_reg(:,i) = tbl.regional-(coef_fitlme_grpwave(i,3)*(tbl.Sex=='Male')+coef_fitlme_grpwave(i,4)*tbl.Age);
end

p_fitlme_grpwave=array2table(p_fitlme_grpwave);
p_fitlme_grpwave.Properties.VariableNames = {'Group','Wave','Sex','Age','Global','GroupWave'};
t_fitlme_grpwave=array2table(t_fitlme_grpwave);
t_fitlme_grpwave.Properties.VariableNames = {'Group','Wave','Sex','Age','Global','GroupWave'};
coef_fitlme_grpwave=array2table(coef_fitlme_grpwave);
coef_fitlme_grpwave.Properties.VariableNames = {'Group','Wave','Sex','Age','Global','GroupWave'};

measure_reg = array2table(measure_reg); 
measure_reg.Properties.VariableNames = measure.Properties.VariableNames;

end