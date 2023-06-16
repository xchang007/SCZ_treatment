%% 1.Load GAMBA toolbox gene expression 
% https://github.com/yongbin-wei/GAMBA
% GE_aparc_2mm_denoised.mat courtesy of Yongbin Wei

clear; clc; close all;
GAMBAPath = 'E:\Matlab\Toolbox_Fudan\Genetic\GAMBA-main\';
addpath(genpath(GAMBAPath))

% load gene expression
ge = load(fullfile(GAMBAPath, 'data', 'GE_aparc_2mm_denoised.mat'));
gene_symbol = ge.gene_symbol;
GG = nanmean(ge.gene_expression_region_FS_z,3); % 82 regions * 20949 genes
regionGE = ge.regionDescriptions;

% save('R1_GOI.mat', 'GG','regionGE')

%% 2.2018_lancet_Antipsychotic, doi: 10.1016/S2215-0366(18)30049-X

study_2018_lancet_Antipsychotic = readtable('2018-lancet-Antipsychotic_Table2.xlsx');
genes_2018_lancet_Antipsychotic = study_2018_lancet_Antipsychotic.NearbyGene;

[tf,idx_gs] = ismember(gene_symbol, genes_2018_lancet_Antipsychotic);
nnz(tf)   % 5/5 genes found in gene_symbol
GE_2018_lancet_Antipsychotic = nanmean(GG(:, tf), 2);              

% save('R1_GOI.mat', 'GE_*','genes_*','-append')


%% 3.1 Inflammation: 2011-plos-WBC_Table2, doi: 10.1371/journal.pgen.1002113

genes_2011_plos_WBC = readtable('2011-plos-WBC_Table2.xlsx');
genes = genes_2011_plos_WBC;

GE = zeros(size(GG,1),size(genes,2));
num_of_genes =  zeros(2,size(genes,2));
for i=1:size(genes,2)
    tf = ismember(gene_symbol, table2cell(genes(:,i)));
    GE(:,i) = nanmean(GG(:, tf), 2);
    num_of_genes(:,i)  = [nnz(~cellfun(@isempty,table2cell(genes(:,i)))); nnz(tf)];    % num_of_genes found in gene_symbol
end

GE_2011_plos_WBC = GE;

% save('R1_GOI.mat', 'GE_*','genes_*','-append')


%% 3.2 Inflammation: 2014_HMG_WBC, doi: 10.1093/hmg/ddu401.

genes_2014_HMG_WBC = readtable('2014-HMG-WBC_Table2.xlsx');
genes = genes_2014_HMG_WBC;

GE = zeros(size(GG,1),size(genes,2));
num_of_genes =  zeros(2,size(genes,2));
for i=1:size(genes,2)
    tf = ismember(gene_symbol, table2cell(genes(:,i)));
    GE(:,i) = nanmean(GG(:, tf), 2);
    num_of_genes(:,i)  = [nnz(~cellfun(@isempty,table2cell(genes(:,i)))); nnz(tf)];    % num_of_genes found in gene_symbol
end

GE_2014_HMG_WBC = GE;

% save('R1_GOI.mat', 'GE_*','genes_*','-append')






