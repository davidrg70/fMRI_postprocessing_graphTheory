% dg_calculate_ModuleDegree_ParticipationCoefficient.m
% David Garnica, UNC Psychiatry, February 2025
% using BCT 2019_03_03
clear all; close all; clc;

BMAP_dir = '/users/d/g/dga/BrainMAP/';
conn_work_dir = fullfile(BMAP_dir,'conn_work');
gtmetrics_dir = fullfile(BMAP_dir,'EF_data3/GT_metrics/');
EFclusters_dir = fullfile(BMAP_dir,'EF_data3/EF_clustering');
tools_dir = '/users/d/g/dga/tools/';
atlas_dir = fullfile(tools_dir,'atlases/Seitzman');
bct_dir = fullfile(tools_dir,'BCT_2019_03_03');
boxplots_dir = fullfile(tools_dir,'DataViz/daboxplot/');
netviewer_dir = fullfile(tools_dir,'BrainNetViewer/');
addpath(genpath(bct_dir));
addpath(BMAP_dir);
addpath(tools_dir);
addpath(boxplots_dir);
addpath(netviewer_dir);
cd(gtmetrics_dir);
fprintf('Please, first select the FC matrices file to work with! \n');
fprintf('Within: %s \n', gtmetrics_dir);
fcfile = uigetfile(gtmetrics_dir);
load(fcfile);

%% determine the proper community affiliation vector (Ci)
% based on Seitzman atlas!!
sourcenames = string(fc_data.sources); % fc_data.sources has the ROIs ANALYZED!
sourcenames = erase(sourcenames, "S300."); % remove the S300 suffix that CONN added
sourcenames = regexprep(sourcenames, '^\d+', ''); % remove the numbers

S_MNI_allInfo = readtable(fullfile(atlas_dir, 'ROIs_300inVol_MNI_allInfo.txt'),...
    "FileType", "text", 'Delimiter', 'space');
S_AnatLabels = readtable(fullfile(atlas_dir, 'ROIs_anatomicalLabels.txt'),...
    "FileType", "text", 'Delimiter', 'space');
[uniqueNetworks,~,ix] = unique(S_MNI_allInfo.netName);
[uniqueLabels, ~, Ci] = unique(sourcenames, 'stable'); % get unique labels and assign numerical values

%% apply BCT Module Degree (without Z score)
% 1 Z (WMD) array per subject!
% 1 sum result (Z) value per network!

flag = 0; % default, to specify an undirected graph
Z = cell(1,numel(fc_data.matrices));
for i = 1:numel(fc_data.matrices)
    W = fc_data.matrices{i,1};
    W(isnan(W)) = 0; %  Replace NaNs with Zeros in the FC Matrix!!!!
    W = (W + W') / 2; % Ensure symmetry
    W(W < 0) = 0; % Remove negative weights (keep positive only)
    Z{i} = module_degree_NOTzscore(W,Ci,flag);
end

% Mean within-module degree per module: an average of all node values within each module/network
% informing the typical within-module degree for each network
mean_WD_per_module = struct;
for k = 1:numel(uniqueLabels) % N number of modules/networks!
    field = strcat(uniqueLabels(k));
    mean_WD_per_module.(field) = [];

    for j = 1:numel(fc_data.matrices) % loop across subjects!
        zvals = Z{1,j}(Ci==k); % indexes all Z/WMD values of that module/network
        average = mean(zvals);
        mean_WD_per_module.(field){j} = average;
    end
end

%% apply BCT Participation Coeff
% 1 P (PC) array per subject

flag = 0; % default, to specify an undirected graph
P = cell(1,numel(fc_data.matrices));
for i = 1:numel(fc_data.matrices)
    W = fc_data.matrices{i,1};
    W(isnan(W)) = 0;  %  Replace NaNs with Zeros in the FC Matrix!!!!
    W = (W + W') / 2; % Ensure symmetry
    [Ppos,Pneg]=participation_coef_sign(W, Ci);
    P{i} = Ppos; % after Jessica's suggestion June 4th 2025, take only positive values!
end

% Jessica considers relevant the PC per network! April 11th 2025
% so I calculate averages per network and per subject too!
mean_PC_per_module = struct;
for k = 1:numel(uniqueLabels) % N number of modules/networks!
    field = strcat(uniqueLabels(k));
    mean_PC_per_module.(field) = [];

    for j = 1:numel(fc_data.matrices) % loop across subjects!
        pcvals = P{1,j}(Ci==k); % indexes all Z/WMD values of that module/network
        average = mean(pcvals);
        mean_PC_per_module.(field){j} = average;
    end
end

%% prepare the EF clusters comparison
load(fullfile(EFclusters_dir, 'Three_clusters_membership_ADHDandTD.mat'));
% it loads 2 variables: membership, and names

load(fullfile(conn_work_dir, 'rest_EF3.mat'));
% I need to load the CONN_x struct to determine the actual subjects analyzed with FC
% which are not the same that were analyzed with npsych and clustering-Walktrap!!!!!
subjectIDs = [];
for i = 1:CONN_x.Setup.nsubjects
    [~,fileName,~] = fileparts(CONN_x.Setup.structural{1,i}{1,1}{1,1});
    % Use regular expression to extract the number after 'Subject0'
    tokens = regexp(fileName, 'sub-(\d+)', 'tokens');
    if ~isempty(tokens)
        % Convert extracted token to a number and add to array
        subjectIDs(i) = str2double(tokens{1}{1});
    end
end
subjectIDs = subjectIDs';

% remove those ids in names and membership that are not in subjectIDs
numeric_from_names = cellfun(@(x) str2double(regexp(x, '\d+', 'match', 'once')), names);
[C,ia] = setdiff(numeric_from_names,subjectIDs);
names(ia) = [];
membership(ia) = [];

% now, get WD and PartCoeff values per cluster!
idx = (membership==1);
WDcluster1_CinguloOpercular = cell2mat(mean_WD_per_module.CinguloOpercular(idx))';
WDcluster1_DefaultMode = cell2mat(mean_WD_per_module.DefaultMode(idx))';
WDcluster1_FrontoParietal = cell2mat(mean_WD_per_module.FrontoParietal(idx))';
WDcluster1 = Z(idx);
PCcluster1_CinguloOpercular = cell2mat(mean_PC_per_module.CinguloOpercular(idx))';
PCcluster1_DefaultMode = cell2mat(mean_PC_per_module.DefaultMode(idx))';
PCcluster1_FrontoParietal = cell2mat(mean_PC_per_module.FrontoParietal(idx))';
PCcluster1 = P(idx);

idx = (membership==2);
WDcluster2_CinguloOpercular = cell2mat(mean_WD_per_module.CinguloOpercular(idx))';
WDcluster2_DefaultMode = cell2mat(mean_WD_per_module.DefaultMode(idx))';
WDcluster2_FrontoParietal = cell2mat(mean_WD_per_module.FrontoParietal(idx))';
WDcluster2 = Z(idx);
PCcluster2_CinguloOpercular = cell2mat(mean_PC_per_module.CinguloOpercular(idx))';
PCcluster2_DefaultMode = cell2mat(mean_PC_per_module.DefaultMode(idx))';
PCcluster2_FrontoParietal = cell2mat(mean_PC_per_module.FrontoParietal(idx))';
PCcluster2 = P(idx);

idx = (membership==3);
WDcluster3_CinguloOpercular = cell2mat(mean_WD_per_module.CinguloOpercular(idx))';
WDcluster3_DefaultMode = cell2mat(mean_WD_per_module.DefaultMode(idx))';
WDcluster3_FrontoParietal = cell2mat(mean_WD_per_module.FrontoParietal(idx))';
WDcluster3 = Z(idx);
PCcluster3_CinguloOpercular = cell2mat(mean_PC_per_module.CinguloOpercular(idx))';
PCcluster3_DefaultMode = cell2mat(mean_PC_per_module.DefaultMode(idx))';
PCcluster3_FrontoParietal = cell2mat(mean_PC_per_module.FrontoParietal(idx))';
PCcluster3 = P(idx);

% bring all WD values, of the 3 clusters, to a matrix, as the kruskalwallis function requires
WD_CinguloOpercular_3clusters = [WDcluster1_CinguloOpercular;...
                                 WDcluster2_CinguloOpercular;...
                                 WDcluster3_CinguloOpercular];
WD_DefaultMode_3clusters = [WDcluster1_DefaultMode;...
                            WDcluster2_DefaultMode;...
                            WDcluster3_DefaultMode];
WD_FrontoParietal_3clusters = [WDcluster1_FrontoParietal;...
                               WDcluster2_FrontoParietal;...
                               WDcluster3_FrontoParietal];

% bring all PC values, of the 3 clusters, to a matrix, as the kruskalwallis function requires
PC_CinguloOpercular_3clusters = [PCcluster1_CinguloOpercular;...
                                 PCcluster2_CinguloOpercular;...
                                 PCcluster3_CinguloOpercular];
PC_DefaultMode_3clusters = [PCcluster1_DefaultMode;...
                            PCcluster2_DefaultMode;...
                            PCcluster3_DefaultMode];
PC_FrontoParietal_3clusters = [PCcluster1_FrontoParietal;...
                               PCcluster2_FrontoParietal;...
                               PCcluster3_FrontoParietal];

% Create grouping variable (applies to the 3 brain networks)
% AND IT IS STRICTLY NECESSARY FOR KRUSKAL-WALLIS: it tells MATLAB which data point belongs to which group
group = [repmat(1, length(WDcluster1_DefaultMode), 1); ...
         repmat(2, length(WDcluster2_DefaultMode), 1); ...
         repmat(3, length(WDcluster3_DefaultMode), 1)];

% add-on BUT IMPORTANT: determine #, ADHDs and TDs per cluster!
demo_data = readtable('/users/d/g/dga/BrainMAP/Demographics Form.csv');
for i = 1:size(subjectIDs,1)
    [row,~] = find(demo_data.sub_id == subjectIDs(i));
    dx = demo_data.adhd_diag{row,1};
    if contains(dx,'ADHD')
        clusterNum = membership(i);
        switch clusterNum
            case 1
                ADHDs_cluster1(i) = subjectIDs(i);
            case 2
                ADHDs_cluster2(i) = subjectIDs(i);
            case 3
                ADHDs_cluster3(i) = subjectIDs(i);
        end
    elseif contains(dx,'TD')
        clusterNum = membership(i);
        switch clusterNum
            case 1
                TDs_cluster1(i) = subjectIDs(i);
            case 2
                TDs_cluster2(i) = subjectIDs(i);
            case 3
                TDs_cluster3(i) = subjectIDs(i);
        end
    end
end
ADHDs_cluster1 = ADHDs_cluster1(ADHDs_cluster1~=0);
ADHDs_cluster2 = ADHDs_cluster2(ADHDs_cluster2~=0);
ADHDs_cluster3 = ADHDs_cluster3(ADHDs_cluster3~=0);
TDs_cluster1 = TDs_cluster1(TDs_cluster1~=0);
TDs_cluster2 = TDs_cluster2(TDs_cluster2~=0);
TDs_cluster3 = TDs_cluster3(TDs_cluster3~=0);
ClusterCounts.NCluster1 = length(ADHDs_cluster1)+length(TDs_cluster1);
ClusterCounts.NCluster2 = length(ADHDs_cluster2)+length(TDs_cluster2);
ClusterCounts.NCluster3 = length(ADHDs_cluster3)+length(TDs_cluster3);
ClusterCounts.ADHDs_Cluster1 = length(ADHDs_cluster1);
ClusterCounts.ADHDs_Cluster2 = length(ADHDs_cluster2);
ClusterCounts.ADHDs_Cluster3 = length(ADHDs_cluster3);
ClusterCounts.TDs_Cluster1 = length(TDs_cluster1);
ClusterCounts.TDs_Cluster2 = length(TDs_cluster2);
ClusterCounts.TDs_Cluster3 = length(TDs_cluster3);
fprintf('Cluster 1 has %d participants \n', ClusterCounts.NCluster1);
fprintf('... %d participants with ADHD \n', ClusterCounts.ADHDs_Cluster1);
fprintf('... %d participants as TD \n', ClusterCounts.TDs_Cluster1);
fprintf('Cluster 2 has %d participants \n', ClusterCounts.NCluster2);
fprintf('... %d participants with ADHD \n', ClusterCounts.ADHDs_Cluster2);
fprintf('... %d participants as TD \n', ClusterCounts.TDs_Cluster2);
fprintf('Cluster 3 has %d participants \n', ClusterCounts.NCluster3);
fprintf('... %d participants with ADHD \n', ClusterCounts.ADHDs_Cluster3);
fprintf('... %d participants as TD \n', ClusterCounts.TDs_Cluster3);

% % save data needed to plot in RStudio!
% writematrix(ADHDs_cluster1, fullfile(gtmetrics_dir,'ADHDs_cluster1.txt'));
% writematrix(ADHDs_cluster2, fullfile(gtmetrics_dir,'ADHDs_cluster2.txt'));
% writematrix(ADHDs_cluster3, fullfile(gtmetrics_dir,'ADHDs_cluster3.txt'));
% writematrix(TDs_cluster1, fullfile(gtmetrics_dir,'TDs_cluster1.txt'));
% writematrix(TDs_cluster2, fullfile(gtmetrics_dir,'TDs_cluster2.txt'));
% writematrix(TDs_cluster3, fullfile(gtmetrics_dir,'TDs_cluster3.txt'));

% find Dx of each ID, and add it to a TOTAL table
Dx = [];
for i = 1:length(subjectIDs)
    sub = subjectIDs(i);
    [row,~] = find(demo_data.sub_id == sub);
    if ~isempty(row)
        Dx{i} = string(demo_data.adhd_diag(row));
    end
end

GT_table = table(subjectIDs,...
    membership, Dx',...
    cell2mat(mean_WD_per_module.CinguloOpercular'), ...
    cell2mat(mean_WD_per_module.DefaultMode'), ...
    cell2mat(mean_WD_per_module.FrontoParietal'), ...
    cell2mat(mean_PC_per_module.CinguloOpercular'), ...
    cell2mat(mean_PC_per_module.DefaultMode'),...
    cell2mat(mean_PC_per_module.FrontoParietal'));
GT_table.Properties.VariableNames = {'IDs','membership','Dx',...
    'WD_CO','WD_DM','WD_FP','PC_CO','PC_DM','PC_FP'};
writetable(GT_table, 'WD_PC_fMRI_analyzed_IDs.csv', 'WriteVariableNames', true);

% DETERMINE MEDIANS OF EACH CLUSTER PER GT METRIC, AND PASTE THIS IN THE EXCEL TABLE!!
Median_Cluster1_WD_CO = median(GT_table.WD_CO(GT_table.membership == 1));
Median_Cluster2_WD_CO = median(GT_table.WD_CO(GT_table.membership == 2));
Median_Cluster3_WD_CO = median(GT_table.WD_CO(GT_table.membership == 3));
Median_Cluster1_WD_DM = median(GT_table.WD_DM(GT_table.membership == 1));
Median_Cluster2_WD_DM = median(GT_table.WD_DM(GT_table.membership == 2));
Median_Cluster3_WD_DM = median(GT_table.WD_DM(GT_table.membership == 3));
Median_Cluster1_WD_FP = median(GT_table.WD_FP(GT_table.membership == 1));
Median_Cluster2_WD_FP = median(GT_table.WD_FP(GT_table.membership == 2));
Median_Cluster3_WD_FP = median(GT_table.WD_FP(GT_table.membership == 3));

Median_Cluster1_PC_CO = median(GT_table.PC_CO(GT_table.membership == 1));
Median_Cluster2_PC_CO = median(GT_table.PC_CO(GT_table.membership == 2));
Median_Cluster3_PC_CO = median(GT_table.PC_CO(GT_table.membership == 3));
Median_Cluster1_PC_DM = median(GT_table.PC_DM(GT_table.membership == 1));
Median_Cluster2_PC_DM = median(GT_table.PC_DM(GT_table.membership == 2));
Median_Cluster3_PC_DM = median(GT_table.PC_DM(GT_table.membership == 3));
Median_Cluster1_PC_FP = median(GT_table.PC_FP(GT_table.membership == 1));
Median_Cluster2_PC_FP = median(GT_table.PC_FP(GT_table.membership == 2));
Median_Cluster3_PC_FP = median(GT_table.PC_FP(GT_table.membership == 3));
Medians_table = table(Median_Cluster1_WD_CO,Median_Cluster2_WD_CO,Median_Cluster3_WD_CO,...
                      Median_Cluster1_WD_DM,Median_Cluster2_WD_DM,Median_Cluster3_WD_DM,...
                      Median_Cluster1_WD_FP,Median_Cluster2_WD_FP,Median_Cluster3_WD_FP,...
                      Median_Cluster1_PC_CO,Median_Cluster2_PC_CO,Median_Cluster3_PC_CO,...
                      Median_Cluster1_PC_DM,Median_Cluster2_PC_DM,Median_Cluster3_PC_DM,...
                      Median_Cluster1_PC_FP,Median_Cluster2_PC_FP,Median_Cluster1_PC_FP);
Medians_table = rows2vars(Medians_table);

%% PLOT WD and PC DATA DISTRIBUTION, USING THE DABOXPLOT FUNCTION!
close all;
[figures] = WD_PC_distribution(GT_table);
warning('Y-axis of text within plots (referring to network names) is hard coded. Do check it if changes are done!');

%% PLOT EDGES AND NODES DISTRIBUTION WITH BRAIN SPACE GRAPHS (BrainNetViewer)

% 1) Build cluster means and write BrainNet files
% ========== INPUTS ==========
N   = length(fc_data.sources);
S   = numel(fc_data.matrices);   % should be the number of subjects
coords = fc_data.MNIcoordinates; % [N nodes x 3]
Ci_mod = Ci(:);                  % [N nodes x 1] module/network per ROI

% Choose a sparsity for display (tune as needed)
prop_keep = 0.10;   % 10% strongest edges
pos_only  = true;   % drop negative weights for visualization

% Output folders
outDir = fullfile(gtmetrics_dir,'BNV_outputs');
edgeDir = fullfile(outDir,'edges');
if ~exist(outDir,'dir'), mkdir(outDir); end
if ~exist(edgeDir,'dir'), mkdir(edgeDir); end

% ========== SUBJECT INDICES PER CLUSTER ==========
idx1 = find(membership == 1);
idx2 = find(membership == 2);
idx3 = find(membership == 3);

fprintf('Subjects per cluster: C1=%d, C2=%d, C3=%d\n', numel(idx1), numel(idx2), numel(idx3));

% ========== CLUSTER MEAN MATRICES ==========
mean_stack = @(C) mean(cat(3, C{:}), 3, 'omitnan');  % C is a cell array of N×N matrices
M1 = mean_stack(fc_data.matrices(idx1));
M2 = mean_stack(fc_data.matrices(idx2));
M3 = mean_stack(fc_data.matrices(idx3));

% Clean / standardize (sym, zero diag)
clean = @(W) (W + W.')/2;
M1 = clean(M1); M1(1:N+1:end) = 0;
M2 = clean(M2); M2(1:N+1:end) = 0;
M3 = clean(M3); M3(1:N+1:end) = 0;

% ========== NODE SIZES ==========
% Use a *global* size reference so node sizes are comparable across clusters.
M_all = clean(mean(cat(3, M1, M2, M3), 3, 'omitnan'));
M_all(1:size(M_all,1)+1:end) = 0;
str_all = node_strength(M_all, 'pos'); % node_strength(W, 'pos'|'abs'|'raw')

% Robust scaling 2..8
minSz = 2; maxSz = 8;
rngVal = max(str_all) - min(str_all);
if rngVal < eps
    sz = repmat((minSz+maxSz)/2, N, 1);
else
    sz = minSz + (maxSz - minSz) * (str_all - min(str_all)) / (rngVal + eps);
end

% ========== .NODE FILE (shared across clusters) ==========
% Node columns: x y z  color size label
% color := Ci_mod (module ID), size := sz, label := simple placeholders
nodeFile = fullfile(outDir, 'Seitzman131_modules.node');
fid = fopen(nodeFile, 'w');
for i = 1:N
    fprintf(fid, '%.2f %.2f %.2f %d %.2f ROI%03d\n', ...
        coords(i,1), coords(i,2), coords(i,3), Ci_mod(i), sz(i), i);
end
fclose(fid);
fprintf('Wrote node file: %s\n', nodeFile);

% ========== .EDGE FILES (one per cluster, thresholded) ==========
E1 = threshold_proportional(M1, prop_keep, pos_only);
E2 = threshold_proportional(M2, prop_keep, pos_only);
E3 = threshold_proportional(M3, prop_keep, pos_only);

edgeFile1 = fullfile(edgeDir, 'Cluster1.edge');
edgeFile2 = fullfile(edgeDir, 'Cluster2.edge');
edgeFile3 = fullfile(edgeDir, 'Cluster3.edge');

writematrix(E1, edgeFile1, 'FileType','text', 'Delimiter',' ');
writematrix(E2, edgeFile2, 'FileType','text', 'Delimiter',' ');
writematrix(E3, edgeFile3, 'FileType','text', 'Delimiter',' ');

fprintf('Wrote edge files:\n  %s\n  %s\n  %s\n', edgeFile1, edgeFile2, edgeFile3);

% 2) View/export with BrainNet Viewer
% OPTION A — GUI
% Open BrainNet Viewer
% Load File → Surface: BrainMesh_ICBM152.nv (BNV’s template: BrainNetViewer/Data/SurfTemplate/BrainMesh_ICBM152.nv)
% Node: Seitzman131_modules.node
% Edge: e.g., Cluster1.edge
% Click Draw (camera icon), adjust styles, repeat for Cluster2/3.

% OPTION B — Batch PNGs from MATLAB
surfFile = fullfile(fileparts(which('BrainNet_MapCfg')), 'Data', 'SurfTemplate', 'BrainMesh_ICBM152.nv');
outImgDir = fullfile(outDir, 'png');
if ~exist(outImgDir,'dir'), mkdir(outImgDir); end

% default config (you can save a .mat cfg from BNV GUI and pass it as 5th arg)
BrainNet_MapCfg(surfFile, nodeFile, edgeFile1, fullfile(outImgDir,'Cluster1.png'));
BrainNet_MapCfg(surfFile, nodeFile, edgeFile2, fullfile(outImgDir,'Cluster2.png'));
BrainNet_MapCfg(surfFile, nodeFile, edgeFile3, fullfile(outImgDir,'Cluster3.png'));
fprintf('Exported PNGs to %s\n', outImgDir);
warning('USE BNV GUI TO SET COLORS FOR NODES+EDGES!');

%% perform comparisons of WD!
close all;
[p1_wd,tbl1_wd,stats1_wd] = kruskalwallis(WD_CinguloOpercular_3clusters, group);
title('Kruskal-Wallis test for WD of 3 EF clusters: CinguloOpercular');
figure;
c = multcompare(stats1_wd, 'CType', 'dunn-sidak');
tbl1_WD_CO_multComp = array2table(c,"VariableNames",["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);

close all;
[p2_wd,tbl2_wd,stats2_wd] = kruskalwallis(WD_DefaultMode_3clusters, group);
title('Kruskal-Wallis test for WD of 3 EF clusters: DefaultMode');
figure;
c = multcompare(stats2_wd, 'CType', 'dunn-sidak');
tbl2_WD_DM_multComp = array2table(c,"VariableNames",["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);

close all;
[p3_wd,tbl3_wd,stats3_wd] = kruskalwallis(WD_FrontoParietal_3clusters, group);
title('Kruskal-Wallis test for WD of 3 EF clusters: FrontoParietal');
figure;
c = multcompare(stats3_wd, 'CType', 'dunn-sidak');
tbl3_WD_FP_multComp = array2table(c,"VariableNames",["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);

% Note: 'dunn-sidak' is is more conservative than no correction, 
% but it less conservative than Bonferroni. The other options are
% 'bonferroni' — very conservative, 'lsd' — no correction (not recommended
% unless exploratory), 'tukey-kramer' — for parametric (ANOVA)

%% perform comparisons of PC!
close all;
[p1_pc,tbl1_pc,stats1_pc] = kruskalwallis(PC_CinguloOpercular_3clusters, group);
title('Kruskal-Wallis test for WD of 3 EF clusters: CinguloOpercular');
figure;
c = multcompare(stats1_pc, 'CType', 'dunn-sidak');
tbl1_PC_CO_multComp = array2table(c,"VariableNames",["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);

close all;
[p2_pc,tbl2_pc,stats2_pc] = kruskalwallis(PC_DefaultMode_3clusters, group);
title('Kruskal-Wallis test for WD of 3 EF clusters: DefaultMode');
figure;
c = multcompare(stats2_pc, 'CType', 'dunn-sidak');
tbl2_PC_DM_multComp = array2table(c,"VariableNames",["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);

close all;
[p3_pc,tbl3_pc,stats3_pc] = kruskalwallis(PC_FrontoParietal_3clusters, group);
title('Kruskal-Wallis test for WD of 3 EF clusters: FrontoParietal');
figure;
c = multcompare(stats3_pc, 'CType', 'dunn-sidak');
tbl3_PC_FP_multComp = array2table(c,"VariableNames",["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);

%% Again kruskalwallis to compare the WHOLE DISTRIBUTION (not averaged per network and subject)
% % between clusters in WD and PartCoeff (all networks/values!)
% Z = Z';
% idx = (membership==1);
% c1 = cell2mat(Z(idx));
% idx = (membership==2);
% c2 = cell2mat(Z(idx));
% idx = (membership==3);
% c3 = cell2mat(Z(idx));
% WD_3clusters_Cumulative = [c1;c2;c3];
% 
% % update the group variable!!
% % now to allow the comparison of the all subjects data and their cumulative distribution
% group = [repmat(1, length(c1), 1); ...
%          repmat(2, length(c2), 1); ...
%          repmat(3, length(c3), 1)];
% 
% P = P';
% idx = (membership==1);
% c1 = cell2mat(P(idx));
% idx = (membership==2);
% c2 = cell2mat(P(idx));
% idx = (membership==3);
% c3 = cell2mat(P(idx));
% PC_3clusters_Cumulative = [c1;c2;c3];

%% About comparisons of WD CUMULATIVE or PC CUMULATIVE (WHOLE DISTRIBUTION)!
% close all;
% [p4_wd,tbl4_wd,stats4_wd] = kruskalwallis(WD_3clusters_Cumulative, group);
% title('KW test for cumulative WD values of 3 EF clusters: 3 networks');
% if p4_wd < 0.05
%     figure;
%     c = multcompare(stats4_wd, 'CType', 'dunn-sidak');
%     tbl4_wd_multComp = array2table(c,"VariableNames", ...
%     ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
% end
% 
% close all;
% [p5_pc,tbl5_pc,stats5_pc] = kruskalwallis(PC_3clusters_Cumulative, group);
% title('KW test for cumulative WD values of 3 EF clusters: 3 networks');
% if p5_pc < 0.05
%     figure;
%     c = multcompare(stats5_pc, 'CType', 'dunn-sidak');
%     tbl4_pc_multComp = array2table(c,"VariableNames", ...
%     ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
% end

%%%%% I INQUIRED CHATGPT ON THIS (Model 5, 08.20.2025):
% After analysis network-level fMRI, I compared numeric values of 3 groups in 3 sets variables: 
% Cluster 1 PC in CinguloOpercular network 
% Cluster 2 PC in CinguloOpercular network 
% Cluster 3 PC in CinguloOpercular network 
% Cluster 1 PC in Default mode network 
% Cluster 2 PC in Default mode network 
% Cluster 3 PC in Default mode network 
% Cluster 1 PC in FrontoParietal network 
% Cluster 2 PC in FrontoParietal network 
% Cluster 3 PC in FrontoParietal network 
% For each network, I ran a Kruskal-Wallis test with posthoc analysis. That is, I compared Clusters 1,2,3 for PC in CinguloOpercular network (1 test), another test for PC in Default mode network, and another test for PC in FrontoParietal network. Any of the 3 tests showed significant differences, neither the posthoc pairwise compairsons. After these comparisons, does it make sense to compare the clusters cumulative distribution of PC? That is, I concatenate the PC values of each cluster (CinguloOpercular, Default mode, and FrontoParietal), and I compare the cumulative WD between clusters.
% Ultimately, I want to know if the clusters integration (participation coefficient -PC-) differs significantly across networks but also across clusters.

%%%%% CHATGPt'S ANSWER
% Short answer: don’t concatenate across networks to compare cumulative distributions. 
% That destroys the network structure and can manufacture significance. 
% It’s fine as a descriptive plot (ECDFs by cluster), but not for inference.

% Either WD or PC, What to do instead (two good options):

% Option A — Mixed model (beta family; best if PC is in [0,1])
% Model PC with Cluster (between-subjects), Network (within-subjects: CO, DMN, FP), plus Subject random intercept.
% Tests: main effects and Cluster×Network interaction; post-hoc EMMs within each network.

% Option B — Repeated-measures MANOVA (fast, interpretable)
% Wide data with columns PC_CO, PC_DMN, PC_FP per subject + cluster.
% Tests the within-subject factor (Network), between-subjects factor (Cluster), and their interaction.

% Extras you might consider
% Within-cluster network effect: Friedman test per cluster (nonparam RM).
% Effect sizes & power: report η² (MANOVA) or semi-partial R²/Δ pseudo-R² (mixed model); non-significant per-network KW could be power-limited.
% Visualization: ECDF or ridge plots of PC by cluster, faceted by network (descriptive only).
% Multiple testing: control FDR for families of post-hoc contrasts.

%% correct for multiple comparisons!
% INCLUDING THE PAIRWISE COMPARISONS P-VALUES! AND ACROSS PC AND WD!
p_vals_toCorrect = [p1_wd;p2_wd;p3_wd;p1_pc;p2_pc;p3_pc;...
                    tbl1_WD_CO_multComp.("P-value");...
                    tbl2_WD_DM_multComp.("P-value");...
                    tbl3_WD_FP_multComp.("P-value");...
                    tbl1_PC_CO_multComp.("P-value");...
                    tbl2_PC_DM_multComp.("P-value");...
                    tbl3_PC_FP_multComp.("P-value")];
[h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(p_vals_toCorrect,.05,'pdep','yes');
