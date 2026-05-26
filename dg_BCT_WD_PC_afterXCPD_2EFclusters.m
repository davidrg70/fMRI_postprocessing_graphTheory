% dg_BCT_WD_PC_afterXCPD_2EFclusters.m
% David Garnica, UNC Psych & Neuro, April 2026
% using BCT 2019_03_03
clear all; close all; clc;

% WORKING WITH WAVE1 AND WAVE2 DATA HOPE POSTPROCESSED WITH XCP-D
% FD threshold already set at 0.1 mm!! ..........................
preproc_dir = '/proj/cohenlab/projects/ADHDBrainMAP/BrainMap_Proc_clpipe191_Nov2024/data_fmriprep/';
postproc_dir = '/proj/cohenlab/projects/ADHDBrainMAP/BrainMap_Proc_clpipe191_Nov2024/data_postproc/XCPD_.1FD_GSR/';

proj_dir = '/users/d/g/dga/BrainMAP/EF_data3/EF_clustering/';
gtmetrics_dir = '/users/d/g/dga/BrainMAP/EF_data3/GT_metrics/';
atlas_dir = '/users/d/g/dga/tools/atlases/Seitzman';
tools_dir = '/users/d/g/dga/tools/';
bct_dir = '/users/d/g/dga/tools/BCT_2019_03_03';
vplots_dir = '/users/d/g/dga/tools/bastibe-Violinplot';
addpath(tools_dir);
addpath(genpath(bct_dir));
addpath(genpath(vplots_dir));
cd(gtmetrics_dir);

% informative count of preprocessed data!
list = dir(preproc_dir);
names = {list.name};
idx = find(contains(names,'sub-')&~contains(names,'.html')&~contains(names,'.sh')...
           &~contains(names,'.json')&~contains(names,'logs')&~contains(names,'HPS')...
           &~contains(names,'atlases')&~contains(names,'working'));
subjects_ids = names(idx);
nW1preproc = sum(contains(subjects_ids,'W01'));
nW2preproc = sum(contains(subjects_ids,'W02'));
nW3preproc = sum(contains(subjects_ids,'W03'));
fprintf('Subjects with preprocessed Wave 1 data (any condition and Nruns): %s \n',string(nW1preproc));
fprintf('Subjects with preprocessed Wave 2 data (any condition and Nruns): %s \n',string(nW2preproc));
fprintf('Subjects with preprocessed Wave 3 data (any condition and Nruns): %s \n',string(nW3preproc));

%% get postproc directory and fMRI condition data
list = dir(postproc_dir);
names = {list.name};
idx = find(contains(names,'sub-')&~contains(names,'.html')&~contains(names,'.sh')...
           &~contains(names,'.json')&~contains(names,'logs')&~contains(names,'HPS')...
           &~contains(names,'atlases')&~contains(names,'working'));
subjects_ids = names(idx)';

% informative count of postprocessed data!
nW1postproc = sum(contains(subjects_ids,'W01'));
nW2postproc = sum(contains(subjects_ids,'W02'));
nW3postproc = sum(contains(subjects_ids,'W03'));
fprintf('Subjects with postprocessed Wave 1 data (any condition and Nruns): %s \n',string(nW1postproc));
fprintf('Subjects with postprocessed Wave 2 data (any condition and Nruns): %s \n',string(nW2postproc));
fprintf('Subjects with postprocessed Wave 3 data (any condition and Nruns): %s \n',string(nW3postproc));

% make user to select WAVE
list = {'1','2','3'};
[indx,~] = listdlg('PromptString', {'Select BrainMAP wave data.', ...
                                     'BL/W1, W2, W3'}, ...
                    'SelectionMode','single', ...
                    'ListString', list, ...
                    'ListSize', [250 50], ...
                    'Name','Postprocessing Selection');
if indx == 1
    waveInd = 'W01';
elseif indx == 2
    waveInd = 'W02';
elseif indx == 3
    waveInd = 'W03';
end
idx = find(~contains(subjects_ids,waveInd));
subjects_ids(idx,:) = []; % REMOVES HERE WAVES NOT INTENDED FOR ANALYSIS

% make user to select either Unscrubbed or Scrubbed data!
% that is, between "_stat-mean_timeseries.tsv" or
% "_stat-mean_timeseries_scrubbed_contigXn.tsv" files
list = {'Not scrubbed','Scrubbed'};
[indx,~] = listdlg('PromptString', {'Select postproc data type.', ...
                                     'After Hope postprocessing Apr 2026.'}, ...
                    'SelectionMode','single', ...
                    'ListString', list, ...
                    'ListSize', [250 50], ...
                    'Name','Postprocessing Selection');
if indx == 1
    scrInd = '_';
elseif indx == 2 % if scrubbed, select CONTIGUITY REQUIREMENTS
    list = {'0','5','9','13'};
    [indx,~] = listdlg('PromptString', {'Select postproc data type.', ...
                                     'After Hope postprocessing Apr 2026.'}, ...
                    'SelectionMode','single', ...
                    'ListString', list, ...
                    'ListSize', [250 80], ...
                    'Name','Postprocessing Selection');
    if indx == 1
        scrInd = 'scrubbed_contigX0';
    elseif indx == 2
        scrInd = 'scrubbed_contigX5';
    elseif indx == 3
        scrInd = 'scrubbed_contigX9';
    elseif indx == 4
        scrInd = 'scrubbed_contigX13';
    end
end

% make user to select fMRI condition (bart,gngreg1/2,gngrew1/2/3/4,rest1/2)
list = {'bart','gngregular','gngreward','rest'};
    [indx,~] = listdlg('PromptString', {'Select fMRI condition.'},...
                    'SelectionMode','single', ...
                    'ListString', list, ...
                    'ListSize', [250 80], ...
                    'Name','Postprocessing Selection');
warning('Number of runs determines number of times WD+PC calculations will be performed!');
warning('WD+PC calculations will be averaged per participant later!')
if indx == 1
    condInd = 'bart';
elseif indx == 2
    condInd = 'gngregular';
elseif indx == 3
    condInd = 'gngreward';
elseif indx == 4
    condInd = 'rest';
end

% report selection of features!
fprintf('------------------------------------------------------ \n');
warning('DATA SELECTED: %s, %s, %s', waveInd, condInd, scrInd);
fprintf('------------------------------------------------------ \n');

% loop over participants and get FC matrices data !!!!!!!!!!!!!!!!!!!!!!!!
% AND GET HOW MANY RUNS THERE ARE, AND STORE .TSV FILENAMES PER PARTICIPANT, TO BE OPENED LATER FOR WD+PC CALCULATIONS
fc_filenames = []; fc_unscrubbed_filenames = [];
for i = 1:numel(subjects_ids)
    subj_dir = fullfile(postproc_dir,subjects_ids{i},'ses-01/func');
    list = dir(subj_dir);
    names = {list.name}';
    idx = find(contains(names,scrInd)&contains(names,condInd)&contains(names,'.tsv') ...
           &~contains(names,'.html')&~contains(names,'.sh')&~contains(names,'.json')...
           &~contains(names,'logs')&~contains(names,'HPS')&~contains(names,'atlases')...
           &~contains(names,'working')&~contains(names,'.gz')&~contains(names,'.nii'));
    if isempty(idx) & contains(scrInd,'scrubbed') % scrubbing selected before
        warning('No scrubbed .tsv/FC matrices of %s',string(subjects_ids{i}));
    elseif isempty(idx) & ~contains(scrInd,'_contigX') % no scrubbing selected before
        warning('No .tsv/FC matrices of %s',string(subjects_ids{i}));
    end
    fc_filenames{i} = names(idx)'; % get files of FC matrices

    % STILL, store the timeseries.tsv files without scrubbing that XCP-D outputs
    % I do this to later know how many volumes were removed with scrubbing with 5/9/13 contiguity requirements
    idx = find(~contains(names,scrInd)&contains(names,condInd)&contains(names,'mean_timeseries.tsv') ...
        &~contains(names,'.html')&~contains(names,'.sh')&~contains(names,'.json')...
        &~contains(names,'logs')&~contains(names,'HPS')&~contains(names,'atlases')...
        &~contains(names,'working')&~contains(names,'.gz')&~contains(names,'.nii'));
    fc_unscrubbed_filenames{i} = names(idx)';
end

% filter after Proportion of Volumes surviving! --> 60% of volumes should be available after scrubbing!
volsThreshold = 40;
FC_matrices = cell(numel(subjects_ids),1);
FCz_matrices = cell(numel(subjects_ids),1);
scrubbed_runsToRemove = [];
for i = 1:numel(subjects_ids) % loop over number of participants
    subj_data = fc_filenames{i};
    subj_unscrubbed_data = fc_unscrubbed_filenames{i};

    if isempty(subj_data)
        FC_matrices{i} = [];
        FCz_matrices{i} = [];
    else
        FC_matrices{i} = cell(1, size(subj_data,2));
        FCz_matrices{i} = cell(1, size(subj_data,2));

        for j = 1:size(subj_data,2)  % loop over number of runs
            run_file = fullfile(postproc_dir, subjects_ids{i}, 'ses-01', 'func', subj_data{j});
            unscrubbed_run_file = fullfile(postproc_dir, subjects_ids{i}, 'ses-01', 'func', subj_unscrubbed_data{j});

            % open each run file... read numeric matrix directly
            X = readmatrix(run_file, 'FileType', 'text', 'Delimiter', '\t');
            Xu = readmatrix(unscrubbed_run_file, 'FileType', 'text', 'Delimiter', '\t');

            % determine proportion of volumes scrubbed
            nvols_scrubbed = size(X,1);
            nvols_unscrubbed = size(Xu,1);
            scrubbing_proportion = 100 - ((nvols_scrubbed * 100) / nvols_unscrubbed);
            if scrubbing_proportion > volsThreshold
                scrubbed_runsToRemove(i,j) = 1;
                % report how many participants and which ones did not survive threshold ...
                warning('Consider excluding %s due to >40p vols scrubbed, run %s of %s', string(subjects_ids{i}), string(j), condInd)
            else
                scrubbed_runsToRemove(i,j) = 0;
            end

            % FC matrix
            FC = corr(X, 'Rows', 'pairwise');
            FC_matrices{i}{j} = FC;
            % Fisher z
            FCz = atanh(FC);  % XCP-D does NOT do Fisher transform (https://neurostars.org/t/does-xcp-d-complete-a-fishers-z-transformation/27760)
            FCz(1:size(FCz,1)+1:end) = 0; % optional: reset diagonal after atanh(1) ... Inf values converted to zeros!
            FCz_matrices{i}{j} = FCz;
        end
    end
end

%% determine the proper community affiliation vector (Ci)
% based on Seitzman atlas!!
atlas_file = '/users/d/g/dga/BrainMAP/xcpd_work/allTasks/Seitzman/atlas-Seitzman/atlas-Seitzman_dseg.tsv';
atlas_tbl  = readtable(atlas_file, 'FileType','text', 'Delimiter','\t');
% Second column contains ROI labels
roi_labels = string(atlas_tbl{:,2});
roi_labels = roi_labels(:);   % column vector
% Extract network names
network_labels = regexprep(roi_labels, '_\d+$', '');
network_labels = strrep(network_labels, 'VentralAttetion', 'VentralAttention');
% Assign one integer per network
[unique_networks, ~, Ci] = unique(network_labels, 'stable');
Ci = Ci(:);
% Summary table
Ci_lookup = table((1:numel(unique_networks))', unique_networks, 'VariableNames', {'CommunityID','Network'});
uniqueLabels = unique_networks;

%% apply BCT Module Degree (without Z score)
% 1 Z (WMD) array per run!
flag = 0; % default, to specify an undirected graph
Z = cell(1,numel(FCz_matrices));
for i = 1:numel(FCz_matrices) % loop over number of participants
    subj_data = FCz_matrices{i};
    if isempty(subj_data)
        Z{i} = [];
    else
        for j = 1:size(subj_data,2)  % loop over number of runs
            W = FCz_matrices{i,1}{1,j};
            W(isnan(W)) = 0; %  Replace NaNs with Zeros in the FC Matrix!!!!
            Z{i}{j} = module_degree_NOTzscore(W,Ci,flag);
        end
    end
end
Z = Z';

% average WD values between 2+ runs... each participant will have a 1x300 double array!
Z_avg = cell(1, numel(Z));
for i = 1:numel(Z)
    subj_data = Z{i};
    if isempty(subj_data) % if both runs are empty
        Z_avg{i} = [];
    else % if 1,2,3 runs, convert to array and average!
        Z_data = cell2mat(subj_data); 
        Z_avg{i} = mean(Z_data,2); % average rows / runs!
    end
end
Z_avg = Z_avg';

% Mean within-module degree per module: an average of all node values within each module/network
% informing the typical within-module degree for each network
mean_WD_per_module = struct;
for k = 1:numel(uniqueLabels) % N number of modules/networks!
    field = strcat(uniqueLabels(k));
    mean_WD_per_module.(field) = [];

    for j = 1:numel(Z_avg) % loop across subjects!
        zvals = Z_avg{j,1};
        if isempty(zvals)
            mean_WD_per_module.(field){j} = [];
        else
            zvals = Z_avg{j,1}(Ci==k); % indexes all Z/WMD values of that module/network
            average = mean(zvals);
            mean_WD_per_module.(field){j} = average;
        end
    end
end

%% apply BCT Participation Coeff
% 1 PC array per run
flag = 0; % default, to specify an undirected graph
P = cell(1,numel(FCz_matrices));
W = []; Z = [];
for i = 1:numel(FCz_matrices)
    subj_data = FCz_matrices{i};
    if isempty(subj_data)
        P{i} = [];
    else
        for j = 1:size(subj_data,2)  % loop over number of runs
            W = FCz_matrices{i,1}{1,j};
            W(isnan(W)) = 0;  %  Replace NaNs with Zeros in the FC Matrix!!!!
            W = (W + W') / 2; % Ensure symmetry
            [Ppos,Pneg]=participation_coef_sign(W, Ci);
            P{i}{j} = Ppos; % after Jessica's suggestion June 4th 2025, take only positive values!
        end
    end
end
P = P';

% average WD values between 2+ runs... each participant will have a 1x300 double array!
P_avg = cell(1, numel(P));
for i = 1:numel(P)
    subj_data = P{i};
    if isempty(subj_data) % if both runs are empty
        P_avg{i} = [];
    else % if 1,2,3 runs, convert to array and average!
        P_data = cell2mat(subj_data); 
        P_avg{i} = mean(P_data,2); % average rows / runs!
    end
end
P_avg = P_avg';

% Jessica considers relevant the PC per network! April 11th 2025
% so I calculate averages per network and per subject too!
mean_PC_per_module = struct;
for k = 1:numel(uniqueLabels) % N number of modules/networks!
    field = strcat(uniqueLabels(k));
    mean_PC_per_module.(field) = [];

    for j = 1:numel(P_avg) % loop across subjects!
        pcvals = P_avg{j,1};
        if isempty(pcvals)
            mean_PC_per_module.(field){j} = [];
        else
            pcvals = P_avg{j,1}(Ci==k); % indexes all Z/WMD values of that module/network
            average = mean(pcvals);
            mean_PC_per_module.(field){j} = average;
        end
    end
end

%% compare EF clusters!
load(fullfile(proj_dir, 'clusters_membership_ADHDandTD.mat')); % data saved with dg_Walktrap.R!!
cd(proj_dir);

subjects_ids = erase(subjects_ids,'sub-');
subjects_ids = erase(subjects_ids,'W01');

uniqueClusters = unique(membership);
% get PC and WD data, and FCz matrices of participant of each cluster!
cluster1_idx = find(membership == 1);
cluster2_idx = find(membership == 2);
cluster1_names = names(cluster1_idx);
cluster2_names = names(cluster2_idx);
subjects_ids = str2double(subjects_ids);

cluster1_WD_FrontoParietal = []; cluster1_WD_DefaultMode = []; cluster1_WD_CinguloOpercular = [];
cluster1_PC_FrontoParietal = []; cluster1_PC_DefaultMode = []; cluster1_PC_CinguloOpercular = [];
cluster1_averaged_FCz_matrices = [];
for i = 1:numel(cluster1_names)
    id = cluster1_names{i,1};
    if contains(id,'a')
        id = erase(id,'a'); id = str2double(id);
    elseif contains(id,'c')
        id = erase(id,'c'); id = str2double(id);
    end
    idx = find(id == subjects_ids);
    if isempty(idx)
        fprintf('No data posprocessed for %s (cluster 1). No PC/WD values indexed \n', num2str(id));
    else
        cluster1_WD_FrontoParietal{i} = cell2mat(mean_WD_per_module.FrontoParietal(idx));
        cluster1_WD_DefaultMode{i} = cell2mat(mean_WD_per_module.DefaultMode(idx));
        cluster1_WD_CinguloOpercular{i} = cell2mat(mean_WD_per_module.CinguloOpercular(idx));
        cluster1_PC_FrontoParietal{i} = cell2mat(mean_PC_per_module.FrontoParietal(idx));
        cluster1_PC_DefaultMode{i} = cell2mat(mean_PC_per_module.DefaultMode(idx));
        cluster1_PC_CinguloOpercular{i} = cell2mat(mean_PC_per_module.CinguloOpercular(idx));

        % also get FCz matrices!
        runs = FCz_matrices{idx,1};
        if isempty(runs)
            cluster1_averaged_FCz_matrices{i} = [];
        else
            valid_runs = runs(~cellfun(@isempty, runs));
            if numel(valid_runs) == 1
                cluster1_averaged_FCz_matrices{i} = valid_runs{1};
            else            
                M = cat(3, valid_runs{:});
                cluster1_averaged_FCz_matrices{i} = mean(M, 3, 'omitnan');
            end
        end
    end
end

cluster2_WD_FrontoParietal = []; cluster2_WD_DefaultMode = []; cluster2_WD_CinguloOpercular = [];
cluster2_PC_FrontoParietal = []; cluster2_PC_DefaultMode = []; cluster2_PC_CinguloOpercular = [];
cluster2_averaged_FCz_matrices = [];
for i = 1:numel(cluster2_names)
    id = cluster2_names{i,1};
    if contains(id,'a')
        id = erase(id,'a'); id = str2double(id);
    elseif contains(id,'c')
        id = erase(id,'c'); id = str2double(id);
    end
    idx = find(id == subjects_ids);
    if isempty(idx)
        fprintf('No data posprocessed for %s (cluster 2). No PC/WD values indexed \n', num2str(id));
    else
        cluster2_WD_FrontoParietal{i} = cell2mat(mean_WD_per_module.FrontoParietal(idx));
        cluster2_WD_DefaultMode{i} = cell2mat(mean_WD_per_module.DefaultMode(idx));
        cluster2_WD_CinguloOpercular{i} = cell2mat(mean_WD_per_module.CinguloOpercular(idx));
        cluster2_PC_FrontoParietal{i} = cell2mat(mean_PC_per_module.FrontoParietal(idx));
        cluster2_PC_DefaultMode{i} = cell2mat(mean_PC_per_module.DefaultMode(idx));
        cluster2_PC_CinguloOpercular{i} = cell2mat(mean_PC_per_module.CinguloOpercular(idx));

        % also get FCz matrices!
        runs = FCz_matrices{idx,1};
        if isempty(runs)
            cluster2_averaged_FCz_matrices{i} = [];
        else
            valid_runs = runs(~cellfun(@isempty, runs));
            if numel(valid_runs) == 1
                cluster2_averaged_FCz_matrices{i} = valid_runs{1};
            else            
                M = cat(3, valid_runs{:});
                cluster2_averaged_FCz_matrices{i} = mean(M, 3, 'omitnan');
            end
        end
    end
end

% remove empty elements of WD and PC clusters cells AND CONVERT TO DOUBLE
cluster1_WD_FrontoParietal = cell2mat(cluster1_WD_FrontoParietal);
cluster1_WD_DefaultMode = cell2mat(cluster1_WD_DefaultMode);
cluster1_WD_CinguloOpercular = cell2mat(cluster1_WD_CinguloOpercular);
cluster1_PC_FrontoParietal = cell2mat(cluster1_PC_FrontoParietal);
cluster1_PC_DefaultMode = cell2mat(cluster1_PC_DefaultMode);
cluster1_PC_CinguloOpercular = cell2mat(cluster1_PC_CinguloOpercular);

cluster2_WD_FrontoParietal = cell2mat(cluster2_WD_FrontoParietal);
cluster2_WD_DefaultMode = cell2mat(cluster2_WD_DefaultMode);
cluster2_WD_CinguloOpercular = cell2mat(cluster2_WD_CinguloOpercular);
cluster2_PC_FrontoParietal = cell2mat(cluster2_PC_FrontoParietal);
cluster2_PC_DefaultMode = cell2mat(cluster2_PC_DefaultMode);
cluster2_PC_CinguloOpercular = cell2mat(cluster2_PC_CinguloOpercular);

% use Mann-Whitney test (Wilcoxon rank sum test) to compare PC and WD between clusters
[p_WD_FP, h_WD_FP, stats] = ranksum(cluster1_WD_FrontoParietal, cluster2_WD_FrontoParietal); zval_WD_FP = stats.zval;
[p_WD_DM, h_WD_DM, stats] = ranksum(cluster1_WD_DefaultMode, cluster2_WD_DefaultMode); zval_WD_DM = stats.zval;
[p_WD_CO, h_WD_CO, stats] = ranksum(cluster1_WD_CinguloOpercular, cluster2_WD_CinguloOpercular); zval_WD_CO = stats.zval;

[p_PC_FP, h_PC_FP, stats] = ranksum(cluster1_PC_FrontoParietal, cluster2_PC_FrontoParietal); zval_PC_FP = stats.zval;
[p_PC_DM, h_PC_DM, stats] = ranksum(cluster1_PC_DefaultMode, cluster2_PC_DefaultMode); zval_PC_DM = stats.zval;
[p_PC_CO, h_PC_CO, stats] = ranksum(cluster1_PC_CinguloOpercular, cluster2_PC_CinguloOpercular); zval_PC_CO = stats.zval;

%% FDR correction and other stats to report!
% SETTINGS
q = 0.05;
method = 'pdep';
report = 'yes';

% sample sizes
n1 = numel(cluster1_WD_FrontoParietal);   % assuming same subjects across WD/PC
n2 = numel(cluster2_WD_FrontoParietal);

% EFFECT SIZES (r = z / sqrt(n1+n2))
r_WD_FP = zval_WD_FP / sqrt(n1 + n2);
r_WD_DM = zval_WD_DM / sqrt(n1 + n2);
r_WD_CO = zval_WD_CO / sqrt(n1 + n2);

r_PC_FP = zval_PC_FP / sqrt(n1 + n2);
r_PC_DM = zval_PC_DM / sqrt(n1 + n2);
r_PC_CO = zval_PC_CO / sqrt(n1 + n2);

% FDR CORRECTION SEPARATELY FOR WD AND PC
pvals_WD = [p_WD_FP; p_WD_DM; p_WD_CO];
[h_WD_fdr, crit_p_WD, adj_ci_cvrg_WD, adj_p_WD] = fdr_bh(pvals_WD, q, method, report);

pvals_PC = [p_PC_FP; p_PC_DM; p_PC_CO];
[h_PC_fdr, crit_p_PC, adj_ci_cvrg_PC, adj_p_PC] = fdr_bh(pvals_PC, q, method, report);

% DESCRIPTIVE STATS
% (median/IQR recommended for ranksum; mean/SD also included)
% WD - FrontoParietal
med1_WD_FP = median(cluster1_WD_FrontoParietal, 'omitnan');
med2_WD_FP = median(cluster2_WD_FrontoParietal, 'omitnan');
iqr1_WD_FP = iqr(cluster1_WD_FrontoParietal);
iqr2_WD_FP = iqr(cluster2_WD_FrontoParietal);
mean1_WD_FP = mean(cluster1_WD_FrontoParietal, 'omitnan');
mean2_WD_FP = mean(cluster2_WD_FrontoParietal, 'omitnan');
sd1_WD_FP = std(cluster1_WD_FrontoParietal, 'omitnan');
sd2_WD_FP = std(cluster2_WD_FrontoParietal, 'omitnan');

% WD - DefaultMode
med1_WD_DM = median(cluster1_WD_DefaultMode, 'omitnan');
med2_WD_DM = median(cluster2_WD_DefaultMode, 'omitnan');
iqr1_WD_DM = iqr(cluster1_WD_DefaultMode);
iqr2_WD_DM = iqr(cluster2_WD_DefaultMode);
mean1_WD_DM = mean(cluster1_WD_DefaultMode, 'omitnan');
mean2_WD_DM = mean(cluster2_WD_DefaultMode, 'omitnan');
sd1_WD_DM = std(cluster1_WD_DefaultMode, 'omitnan');
sd2_WD_DM = std(cluster2_WD_DefaultMode, 'omitnan');

% WD - CinguloOpercular
med1_WD_CO = median(cluster1_WD_CinguloOpercular, 'omitnan');
med2_WD_CO = median(cluster2_WD_CinguloOpercular, 'omitnan');
iqr1_WD_CO = iqr(cluster1_WD_CinguloOpercular);
iqr2_WD_CO = iqr(cluster2_WD_CinguloOpercular);
mean1_WD_CO = mean(cluster1_WD_CinguloOpercular, 'omitnan');
mean2_WD_CO = mean(cluster2_WD_CinguloOpercular, 'omitnan');
sd1_WD_CO = std(cluster1_WD_CinguloOpercular, 'omitnan');
sd2_WD_CO = std(cluster2_WD_CinguloOpercular, 'omitnan');

% PC - FrontoParietal
med1_PC_FP = median(cluster1_PC_FrontoParietal, 'omitnan');
med2_PC_FP = median(cluster2_PC_FrontoParietal, 'omitnan');
iqr1_PC_FP = iqr(cluster1_PC_FrontoParietal);
iqr2_PC_FP = iqr(cluster2_PC_FrontoParietal);
mean1_PC_FP = mean(cluster1_PC_FrontoParietal, 'omitnan');
mean2_PC_FP = mean(cluster2_PC_FrontoParietal, 'omitnan');
sd1_PC_FP = std(cluster1_PC_FrontoParietal, 'omitnan');
sd2_PC_FP = std(cluster2_PC_FrontoParietal, 'omitnan');

% PC - DefaultMode
med1_PC_DM = median(cluster1_PC_DefaultMode, 'omitnan');
med2_PC_DM = median(cluster2_PC_DefaultMode, 'omitnan');
iqr1_PC_DM = iqr(cluster1_PC_DefaultMode);
iqr2_PC_DM = iqr(cluster2_PC_DefaultMode);
mean1_PC_DM = mean(cluster1_PC_DefaultMode, 'omitnan');
mean2_PC_DM = mean(cluster2_PC_DefaultMode, 'omitnan');
sd1_PC_DM = std(cluster1_PC_DefaultMode, 'omitnan');
sd2_PC_DM = std(cluster2_PC_DefaultMode, 'omitnan');

% PC - CinguloOpercular
med1_PC_CO = median(cluster1_PC_CinguloOpercular, 'omitnan');
med2_PC_CO = median(cluster2_PC_CinguloOpercular, 'omitnan');
iqr1_PC_CO = iqr(cluster1_PC_CinguloOpercular);
iqr2_PC_CO = iqr(cluster2_PC_CinguloOpercular);
mean1_PC_CO = mean(cluster1_PC_CinguloOpercular, 'omitnan');
mean2_PC_CO = mean(cluster2_PC_CinguloOpercular, 'omitnan');
sd1_PC_CO = std(cluster1_PC_CinguloOpercular, 'omitnan');
sd2_PC_CO = std(cluster2_PC_CinguloOpercular, 'omitnan');

% RESULTS TABLE
Measure = {'WD'; 'WD'; 'WD'; 'PC'; 'PC'; 'PC'};
Network = {'FrontoParietal'; 'DefaultMode'; 'CinguloOpercular'; ...
           'FrontoParietal'; 'DefaultMode'; 'CinguloOpercular'};
Cluster1_N = repmat(n1, 6, 1);
Cluster2_N = repmat(n2, 6, 1);

Cluster1_Median = [med1_WD_FP; med1_WD_DM; med1_WD_CO; med1_PC_FP; med1_PC_DM; med1_PC_CO];
Cluster1_IQR    = [iqr1_WD_FP; iqr1_WD_DM; iqr1_WD_CO; iqr1_PC_FP; iqr1_PC_DM; iqr1_PC_CO];
Cluster2_Median = [med2_WD_FP; med2_WD_DM; med2_WD_CO; med2_PC_FP; med2_PC_DM; med2_PC_CO];
Cluster2_IQR    = [iqr2_WD_FP; iqr2_WD_DM; iqr2_WD_CO; iqr2_PC_FP; iqr2_PC_DM; iqr2_PC_CO];

Cluster1_Mean = [mean1_WD_FP; mean1_WD_DM; mean1_WD_CO; mean1_PC_FP; mean1_PC_DM; mean1_PC_CO];
Cluster1_SD   = [sd1_WD_FP; sd1_WD_DM; sd1_WD_CO; sd1_PC_FP; sd1_PC_DM; sd1_PC_CO];
Cluster2_Mean = [mean2_WD_FP; mean2_WD_DM; mean2_WD_CO; mean2_PC_FP; mean2_PC_DM; mean2_PC_CO];
Cluster2_SD   = [sd2_WD_FP; sd2_WD_DM; sd2_WD_CO; sd2_PC_FP; sd2_PC_DM; sd2_PC_CO];

Z_value = [zval_WD_FP; zval_WD_DM; zval_WD_CO; zval_PC_FP; zval_PC_DM; zval_PC_CO];
P_raw   = [p_WD_FP; p_WD_DM; p_WD_CO; p_PC_FP; p_PC_DM; p_PC_CO];
Effect_r = [r_WD_FP; r_WD_DM; r_WD_CO; r_PC_FP; r_PC_DM; r_PC_CO];

P_FDR = [adj_p_WD(1); adj_p_WD(2); adj_p_WD(3); ...
         adj_p_PC(1); adj_p_PC(2); adj_p_PC(3)];

H_FDR = [h_WD_fdr(1); h_WD_fdr(2); h_WD_fdr(3); ...
         h_PC_fdr(1); h_PC_fdr(2); h_PC_fdr(3)];

RESULTS_TABLE = table(Measure, Network, ...
    Cluster1_N, Cluster2_N, ...
    Cluster1_Median, Cluster1_IQR, ...
    Cluster2_Median, Cluster2_IQR, ...
    Cluster1_Mean, Cluster1_SD, ...
    Cluster2_Mean, Cluster2_SD, ...
    Z_value, P_raw, P_FDR, Effect_r, H_FDR);

fprintf('---> See table RESULTS_TABLE !! \n');

%% Violinplots for WD and PC
close all;
% Colors for clusters
color = [0.2 0.4 0.9; 
         0.93 0.62, 0.13];

% FIGURE 1: WD
figure('Units','normalized','Position',[0.1 0.1 0.8 0.45]);
wd_titles = {'FrontoParietal', 'DefaultMode', 'CinguloOpercular'};
wd_cluster1 = {cluster1_WD_FrontoParietal, ...
               cluster1_WD_DefaultMode, ...
               cluster1_WD_CinguloOpercular};
wd_cluster2 = {cluster2_WD_FrontoParietal, ...
               cluster2_WD_DefaultMode, ...
               cluster2_WD_CinguloOpercular};
for k = 1:3
    subplot(1,3,k)
    x = wd_cluster1{k};
    y = wd_cluster2{k};
    % combine raw values
    vplotsdata = vertcat(x(:), y(:));
    category1 = repmat("Cluster 1", length(x), 1);
    category2 = repmat("Cluster 2", length(y), 1);
    categories = cellstr(vertcat(category1, category2));
    violinplot(vplotsdata, categories, ...
        'ViolinColor', color, ...
        'ShowWhiskers', true, ...
        'ShowNotches', false, ...
        'ShowMean', false, ...
        'ShowMedian', true);
    title(sprintf('WD - %s', wd_titles{k}));
    ylabel('WD values');
    xlim([0.5, 2.5]);
    set(gca, 'FontSize', 12);
end
sgtitle('Within-module degree by cluster');

% FIGURE 2: PC
figure('Units','normalized','Position',[0.1 0.1 0.8 0.45]);
pc_titles = {'FrontoParietal', 'DefaultMode', 'CinguloOpercular'};
pc_cluster1 = {cluster1_PC_FrontoParietal, ...
               cluster1_PC_DefaultMode, ...
               cluster1_PC_CinguloOpercular};
pc_cluster2 = {cluster2_PC_FrontoParietal, ...
               cluster2_PC_DefaultMode, ...
               cluster2_PC_CinguloOpercular};
for k = 1:3
    subplot(1,3,k)
    x = pc_cluster1{k};
    y = pc_cluster2{k};
    % combine raw values
    vplotsdata = vertcat(x(:), y(:));
    category1 = repmat("Cluster 1", length(x), 1);
    category2 = repmat("Cluster 2", length(y), 1);
    categories = cellstr(vertcat(category1, category2));
    violinplot(vplotsdata, categories, ...
        'ViolinColor', color, ...
        'ShowWhiskers', true, ...
        'ShowNotches', false, ...
        'ShowMean', true, ...
        'ShowMedian', true);
    title(sprintf('PC - %s', pc_titles{k}));
    ylabel('PC values');
    xlim([0.5, 2.5]);
    set(gca, 'FontSize', 12);
end
sgtitle('Participation coefficient by cluster');

%% Boxplots for WD and PC
close all;
addpath('/users/d/g/dga/tools/DataViz/daboxplot');
color = [0.2 0.4 0.9;
         0.93 0.62 0.13];

figure('Units','normalized','Position',[0.1 0.1 0.8 0.45]);
wd_titles = {'FrontoParietal', 'DefaultMode', 'CinguloOpercular'};
wd_cluster1 = {cluster1_WD_FrontoParietal, ...
               cluster1_WD_DefaultMode, ...
               cluster1_WD_CinguloOpercular};
wd_cluster2 = {cluster2_WD_FrontoParietal, ...
               cluster2_WD_DefaultMode, ...
               cluster2_WD_CinguloOpercular};
for k = 1:3
    subplot(1,3,k)
    x = wd_cluster1{k};
    y = wd_cluster2{k};
    daboxplot({x(:), y(:)}, ...
        'fill', 1, ...
        'colors', color, ...
        'whiskers', 0, ...
        'scatter', 1, ...
        'scattersize', 12, ...
        'scatteralpha', 0.6, ...
        'jitter', 1, ...
        'mean', 0, ...
        'outliers', 1, ...
        'boxalpha', 0.7, ...
        'xtlabels', {'Cluster 1','Cluster 2'});
    hold on
    add_manual_whiskers(x(:),1);
    add_manual_whiskers(y(:),2);
    title(sprintf('WD - %s', wd_titles{k}));
    ylabel('WD values');
    set(gca, 'FontSize', 12);
end
sgtitle('WD by cluster');

figure('Units','normalized','Position',[0.1 0.1 0.8 0.45]);
pc_titles = {'FrontoParietal', 'DefaultMode', 'CinguloOpercular'};
pc_cluster1 = {cluster1_PC_FrontoParietal, ...
               cluster1_PC_DefaultMode, ...
               cluster1_PC_CinguloOpercular};
pc_cluster2 = {cluster2_PC_FrontoParietal, ...
               cluster2_PC_DefaultMode, ...
               cluster2_PC_CinguloOpercular};
for k = 1:3
    subplot(1,3,k)
    x = pc_cluster1{k};
    y = pc_cluster2{k};
    daboxplot({x(:), y(:)}, ...
        'fill', 1, ...
        'colors', color, ...
        'whiskers', 0, ...
        'scatter', 1, ...
        'scattersize', 12, ...
        'scatteralpha', 0.6, ...
        'jitter', 1, ...
        'mean', 0, ...
        'outliers', 1, ...
        'boxalpha', 0.7, ...
        'xtlabels', {'Cluster 1','Cluster 2'});
    hold on
    add_manual_whiskers(x(:),1);
    add_manual_whiskers(y(:),2);
    title(sprintf('PC - %s', pc_titles{k}));
    ylabel('PC values');
    set(gca, 'FontSize', 12);
end
sgtitle('PC by cluster');

function add_manual_whiskers(data,xpos)
data = data(~isnan(data));
q1 = prctile(data,25);
q3 = prctile(data,75);
iqr_val = q3 - q1;
lower_lim = q1 - 1.5*iqr_val;
upper_lim = q3 + 1.5*iqr_val;
lower_whisker = min(data(data >= lower_lim));
upper_whisker = max(data(data <= upper_lim));
capw = 0.12;
plot([xpos xpos],[lower_whisker q1],'k-','LineWidth',1);
plot([xpos xpos],[q3 upper_whisker],'k-','LineWidth',1);
plot([xpos-capw xpos+capw],[lower_whisker lower_whisker],'k-','LineWidth',1);
plot([xpos-capw xpos+capw],[upper_whisker upper_whisker],'k-','LineWidth',1);
end

%% Average FCz matrices of each cluster and save them to plot with brainconn in RStudio!
% deal with some empty elements!
valid_idx = ~cellfun(@isempty, cluster1_averaged_FCz_matrices);
valid_mats = cluster1_averaged_FCz_matrices(valid_idx);
nSubj = numel(valid_mats);
all_mats = nan(300,300,nSubj);
for i = 1:nSubj
    all_mats(:,:,i) = valid_mats{i};
end
cluster1_mean_FCz_matrix = mean(all_mats, 3, 'omitnan');

valid_idx = ~cellfun(@isempty, cluster2_averaged_FCz_matrices);
valid_mats = cluster2_averaged_FCz_matrices(valid_idx);
nSubj = numel(valid_mats);
all_mats = nan(300,300,nSubj);
for i = 1:nSubj
    all_mats(:,:,i) = valid_mats{i};
end
cluster2_mean_FCz_matrix = mean(all_mats, 3, 'omitnan');

% extract FC values of FPN (9), DMN(6), and CON (4) nodes! (trim the FCz matrices!)
% AND KEEP ORDER OF NODES AS GIVEN BY SEIZTMAN ATLAS!
mask = Ci == 9 | Ci == 6 | Ci == 4;
idx_target = find(mask);
cluster1_FCz_submatrix = cluster1_mean_FCz_matrix(idx_target, idx_target);
cluster2_FCz_submatrix = cluster2_mean_FCz_matrix(idx_target, idx_target);
Ci_target = Ci(idx_target);

% SAVE CLUSTERS SUBMATRICES!
filename1 = fullfile(gtmetrics_dir,'mean_cluster1_FCz.csv');
filename2 = fullfile(gtmetrics_dir,'mean_cluster2_FCz.csv');
writematrix(cluster1_FCz_submatrix, filename1);
writematrix(cluster2_FCz_submatrix, filename2);
