% dg_GNG_events_regressors_for_XCPD_v4.m
% February to May 2026 - UNC-CH

% IMPORTANT NOTES:

% RUN THIS SCRIPT AFTER RUNNING "DG_BMAP_GNG_create_onsets.R", adaptation of Monica's "MGLMasters_step1_create_onsets.R"
% Monica's scripts saves GNG Reg/Rew events as txt files, which I save at /users/d/g/dga/BrainMAP/EF_data3/gng_psychopy/onsets_wave1/
% Data downloaded from SharePoint - cohenlabteam/Documents/Research Studies/ADHD BrainMAP/Tasks/GNG/data/
% Data saved in Longleaf at /users/d/g/dga/BrainMAP/EF_data3/gng_psychopy/

% the GNG regular/reward events timepoints data is in COHENLABTEAM:
% Documents/Research Studies/ADHD BrainMAP/Tasks/GNG/data ... most updated
% file: per-run_gng_data_08-20-2025.xlsx (------THIS IS ELLIE'S SPREADSHEET!------)
% https://adminliveunc.sharepoint.com/:x:/r/sites/cohenlabteam/_layouts/15/Doc.aspx?sourcedoc=%7BBF9AC788-CD25-46E3-B911-778A64FEA7BE%7D&file=per-run_gng_data_08-20-2025.csv&action=default&mobileredirect=true

% HOW TO SET CUSTOMIZED REGRESSORS FOR XCP-D:
% https://xcp-d.readthedocs.io/en/latest/usage.html#custom-confounds

clear all; close all; clc;
analysis_dir = '/users/d/g/dga/BrainMAP';
ef_data_dir = fullfile(analysis_dir, 'EF_data3');
cd(ef_data_dir);
addpath('/users/d/g/dga/tools/spm12_7487');
addpath(analysis_dir);

% perRunGNG = readtable(fullfile(ef_data_dir,'/per-run_gng_data_12-17-2025.csv'));
filename = fullfile(ef_data_dir,'/per-run_gng_data_04-22-2026.csv');
opts = detectImportOptions(filename,'VariableNamingRule','preserve');
opts.VariableTypes(endsWith(opts.VariableNames,"_type")) = {'categorical'}; % this way, all data is brought!
perRunGNG = readtable(filename, opts);

% DEFINE FMRI DATA DIRECTORIES
rest_data_dir = '/users/d/g/dga/BrainMAP/xcpd_work/data_rest_fmriprep/';
gngreg_data_dir = '/users/d/g/dga/BrainMAP/xcpd_work/data_gngregular_fmriprep/';
gngrew_data_dir = '/users/d/g/dga/BrainMAP/xcpd_work/data_gngreward_fmriprep/';

% list all ids per fMRI condition and make blank rows to RECORD WHO I HAVE POSTPROCESSED ...
list = dir(rest_data_dir); names = {list.name}; idx = find(contains(names,'sub-'));
all_rest_ids = names(idx)';
all_rest_ids_new = cell(2*numel(all_rest_ids)-1, 1);   % creates (nRows*2 -1) x 1 cell array
all_rest_ids_new(1:2:end) = all_rest_ids;              % puts original rows in odd positions
all_rest_ids = all_rest_ids_new; clear all_rest_ids_new;
list = dir(gngreg_data_dir); names = {list.name}; idx = find(contains(names,'sub-'));
all_gngreg_ids = names(idx)';
all_gngreg_ids_new = cell(2*numel(all_gngreg_ids)-1, 1);   % creates (nRows*2 -1) x 1 cell array
all_gngreg_ids_new(1:2:end) = all_gngreg_ids;              % puts original rows in odd positions
all_gngreg_ids = all_gngreg_ids_new; clear all_gngreg_ids_new;
list = dir(gngrew_data_dir); names = {list.name}; idx = find(contains(names,'sub-'));
all_gngrew_ids = names(idx)';
all_gngrew_ids_new = cell(4*numel(all_gngrew_ids)-1, 1);   % creates (nRows*4 -1) x 1 cell array
all_gngrew_ids_new(1:4:end) = all_gngrew_ids;              % puts original rows in odd positions
all_gngrew_ids = all_gngrew_ids_new; clear all_gngrew_ids_new;

fprintf('--> KEEP UPDATING "BrainMAP runs count, preproc and postproc tracker.xlsx" !!! \n');

%% TRIM data spreadsheet first AND REMOVE WITHDRAWALS (avoid/shorten heavy computational work ahead...)

% determine wave 2 and 3 entries/rows
[rows_of2,~] = find(perRunGNG.wave == 2);
[rows_of3,~] = find(perRunGNG.wave == 3);
perRunGNG([rows_of2; rows_of3],:) = []; % remove those rows! I only want wave 1 data

% remove columns I don't care about (interview_date,interview_age,sex,group,wave,task_end_var
% ... DO NOT remove those columns having "_response","_accuracy","_rt")
vars = string(perRunGNG.Properties.VariableNames);
[~,col1] = find(contains(vars,'interview_date') | contains(vars,'interview_age') | contains(vars,'sex'));
[~,col2] = find(contains(vars,'group') | contains(vars,'wave') | contains(vars,'task_end_var'));

cols_idx = [col1, col2];
perRunGNG(:,cols_idx) = [];

% REMOVE --SOME-- WITHDRAWALS THAT WERE SCANNED AND  (KEEP THIS SECTION THE SAME AS IN " dg_trim_EF_data_v4.m")
% WITHDRAWALS ARE NOTED WITH 'Y' IN A COLUMN OF THE DEMOGRAPHICS FORM.CSV (12.03.2025 file version)
filename = [analysis_dir,'/','Demographics Form.csv'];
demographics = readtable(filename);
[ys,~] = find(contains(demographics.withdrawn,'y'));
withdrawals = demographics.sub_id(ys)';
withdrawals = cat(2,[1,2,1259,2079,2122],withdrawals); % add 2,3,and 1259,2079,2122
% The ID "2" was a pilot, coming from shuffling data/spreadsheet
% The ID "1259" is in the updated SST spreadsheet, but I do not know why
% '2079' has incomplete GNG data: missing rew 4; has 1 incomplete reg 2 and one complete reg 2
% '2122' was never enrolled, does not have EF data either

% ----> '2008','2009','2014','2097' withdrew, have 2 RS runs in wave 1, but
% do not have complete EF data! so I removed them from this analysis!!!!!!
% '2009' also has a bad T1
% '2148' withdrew and has only 1 RS run in wave 1. It seems it doesn't have EF data either
withdrawals = sort(withdrawals);
[rows,~] = find(ismember(perRunGNG.subject_id,withdrawals));
perRunGNG(rows,:) = []; % removes the rows of withdrawn participants found

% get unique IDs and GNG tasks
uniqueIDs = unique(perRunGNG.subject_id);
uniqueGNGs = unique(perRunGNG.task);

% just save the trimmed xlsx
filename = fullfile(ef_data_dir,'/per-run_gng_data_05-06-2026_trimmed.xlsx');
writetable(perRunGNG,filename);

% ----- EDITS TO SEVERAL PARTICIPANTS (NOTE: I ERASED SOME NII FILES IN THE DIRS) -----

% NOTE ABOUT sub-2029W01:
% this participant has a sub-2029W01_ses-01_task-gngregular01_run-02 and a
% REPEATED gngregular entry in the perRunGNG table... as the NII file has a
% "run-02" I assume I should take the second line of gngregular run 1 in
% the perRunGNG table... I mean, I think the run1 was repeated and the
% technician just named the file with a "01run-02" and the actual events of
% that nii are the second gngregular line... then I remove line 115 (the first gngregular line)
% before running the whole script section! ---CONFIRM WITH ELLIE IF NII FILE CORRESPONDS ---
perRunGNG(115,:) = [];

% NOTE ABOUT sub-2056W01:
% something similar with this participant. The perRunGNG table has a double
% entry for gngreward run 4. I select the second line because it has a
% duration of 369 sec, which corresponds to the 450+ volumes that the nii image has!
perRunGNG(241,:) = [];

% NOTE ABOUT sub-2075W01:
% something similar with this participant. The perRunGNG table has a double
% entry for gngregular run 1. I select the second line because it has a
% duration of 377 sec, which corresponds to the 473 volumes that the nii image has!
perRunGNG(308,:) = [];

% NOTE ABOUT sub-2104W01:
% something similar with this participant. The perRunGNG table has a double
% entry for gngregular run 1. I select the second line because it has a
% duration of 386 sec, which corresponds to the 483 volumes that the nii image has!
perRunGNG(417,:) = [];

% NOTE ABOUT sub-2117W01:
% something similar with this participant. The perRunGNG table has a double
% entry for gngregular run 1. I select the second line because it has a
% duration of 391 sec, which corresponds to the 495 volumes that the nii image has!
perRunGNG(488,:) = [];

% NOTE ABOUT sub-2125W01:
% something similar with this participant. The perRunGNG table has a double
% entry for gngreward run 3. I select the second line because it has a
% duration of 372 sec, which corresponds to the 450+ volumes that the nii image has!
perRunGNG(527,:) = [];

% NOTE ABOUT sub-2128W01:
% something similar with this participant. The perRunGNG table has a double
% entry for gngregular run 2. The directory had "task-gngregular02_run-01"
% and a "task-gngregular02_run-02" small nii files (<20Mb) that I removed.
% I left the "task-gngregular02_run-03" nii file that corresponds to the
% second line entry in the perRunGNG table, with 534 Mb, 495 vols and 390
% sec length
perRunGNG(537,:) = []; % line with 22 sec length

% NOTE ABOUT sub-2128W01:
% something similar with this participant. The perRunGNG table has a double
% entry for gngregular run 2. The directory had "task-gngregular02_run-01"
% and a "task-gngregular02_run-02_space-MNI152NLin2009cAsym_desc-brain_mask" but
% not preproc... I removed those files and left the
% "task-gngregular02_run-03" nii file that seems to correspond to the second
% entry of gngregular run 2 in the perRunGNG table, the one with 389 sec
perRunGNG(591,:) = [];

% NOTE ABOUT sub-2144W01:
% something similar with this participant. The perRunGNG table has a double
% entry for gngreward run 4. The directory had "sub-2144W01_ses-01_task-gngreward04"
% and a "task-gngreward04_run-02_space-MNI152NLin2009cAsym_desc-brain_mask" but
% not preproc... I removed those files and left the
% "task-gngreward04_run-03" nii file that seems to correspond to the second
% entry of gngreward run 4 in the perRunGNG table, the one with 369 sec
perRunGNG(622,:) = [];

% NOTE ABOUT sub-2170W01:
% something similar with this participant. The perRunGNG table has a double
% entry for gngregular run 1. I select the first line because it has a
% 2.1% of omissions while the second line has 21%! ---CONFIRM WITH ELLIE IF NII FILE CORRESPONDS ---
perRunGNG(736,:) = [];
% ALSO about sub-2170W01:
% the perRUNGNG table shows a repeated gngreward run1 for this participant, but
% there's a run2 in the files directory...so I think that the second run1
% in the perRUNGNG table actually corresponds to the run2 nii... I fix this here
perRunGNG(737,2) = array2table(2);

% NOTE ABOUT sub-2181W01:
% something similar with this participant. The perRunGNG table has a double
% entry for gngregular run 1. I select the second line because it has a
% duration of 391 sec, which corresponds to the 495 volumes that the nii image has!
perRunGNG(770,:) = [];

% ---------------------------------------------------------------------
% other unprocessed IDs: 2112, 2153, 2162, 2164, 2167, 2168, 2173, 2178,
%                        2184, 2189, 2190, 2191, 2193, 2195, 2196, 2197,
%                        2198, 2199, 2200, 2201, 2202, 2203
% ---------------------------------------------------------------------

%% 1. Compare IDs and Runs with Ellie's "per-run_gng_data_mm-dd-yyyy.xlsx" file
% I want to know if some data is not at /users/d/g/dga/BrainMAP/EF_data3/gng_psychopy/ but should be following Ellie's spreadsheet!
events_dir = '/users/d/g/dga/BrainMAP/EF_data3/gng_psychopy/onsets_wave1/';
allRegRewFiles = dir(events_dir);
allRegRewFiles = {allRegRewFiles.name}';

perRunGNG_IDs = unique(perRunGNG.subject_id);
allRegRewFiles_IDs = unique(cellfun(@(x) str2double(regexp(x, '(?<=sub-)\d+', 'match', 'once')), allRegRewFiles));

% know common IDs... the amount
common_IDs = intersect(perRunGNG_IDs, allRegRewFiles_IDs);
[C,ia] = setdiff(allRegRewFiles_IDs, perRunGNG_IDs); % just to know ...

% check if runs of each gngreg and gngrew are the same in both data bases
% each ID should have 2 gngreg and 4 gngrew in each data base............
CheckConcurrence = table;
for i = 1:numel(common_IDs)
    current_id = common_IDs(i);

    idx1 = perRunGNG.subject_id == current_id;
    id_filenames = perRunGNG.filename(idx1);
    countReg_runs1 = sum(contains(string(id_filenames), 'regular'));
    countRew_runs1 = sum(contains(string(id_filenames), 'reward'));

    idx2 = contains(string(allRegRewFiles), string(current_id));
    id_filenames = allRegRewFiles(idx2);
    log_Reg_run1 = any(contains(string(id_filenames), 'gngreg_run-1'));
    log_Reg_run2 = any(contains(string(id_filenames), 'gngreg_run-2'));
    log_Rew_run1 = any(contains(string(id_filenames), 'gngrew_run-1'));
    log_Rew_run2 = any(contains(string(id_filenames), 'gngrew_run-2'));
    log_Rew_run3 = any(contains(string(id_filenames), 'gngrew_run-3'));
    log_Rew_run4 = any(contains(string(id_filenames), 'gngrew_run-4'));

    countReg_runs2 = log_Reg_run1 + log_Reg_run2;
    countRew_runs2 = log_Rew_run1 + log_Rew_run2 + log_Rew_run3 + log_Rew_run4;

    CheckConcurrence.subject_id(i,1) = current_id;
    CheckConcurrence.countReg_runsEllie(i,1) = countReg_runs1;
    CheckConcurrence.countRew_runsEllie(i,1) = countRew_runs1;
    CheckConcurrence.countReg_runsSPd(i,1) = countReg_runs2;
    CheckConcurrence.countRew_runsSPd(i,1) = countRew_runs2;
    CheckConcurrence.ConcurRegRuns(i,1) = double(countReg_runs1 == countReg_runs2);
    CheckConcurrence.ConcurRewRuns(i,1) = double(countRew_runs1 == countRew_runs2);
end

warning('--> HARD CODED REMOVAL OF 2168,2178,2191,2200,2202');
warning('Check later if those IDs are preprocessed!!!!!!!!!');
[rows, ~] = find(ismember(common_IDs, [2168 2178 2191 2200 2202]));
common_IDs(rows,:) = [];
CheckConcurrence(rows,:) = [];

%% 2. MAIN WORK:
% 2a. plot events before
% 2b. plot events after HRF convolution
% 2c. append fMRIPrep and task regressors and save custom_confounds files for XCP-D!

xcpd_work_dir = '/users/d/g/dga/BrainMAP/xcpd_work/';
custom_confounds_root_reg = fullfile(xcpd_work_dir, 'data_gngregular_custom_confounds');
custom_confounds_root_rew = fullfile(xcpd_work_dir, 'data_gngreward_custom_confounds');

if ~exist(custom_confounds_root_reg, 'dir')
    mkdir(custom_confounds_root_reg);
end
if ~exist(custom_confounds_root_rew, 'dir')
    mkdir(custom_confounds_root_rew);
end

for i = 1:numel(common_IDs) % loop over each participant!
    idx_reg = [];
    % for each ID determine indices of each run, based on CheckConcurrence table
    % (countReg_runsSPd and countRew_runsSPd variables, iterate across those numbers)
    current_id = common_IDs(i);

    %% GNG REGULAR!!
    countReg_runs = CheckConcurrence.countReg_runsSPd(i,1);  task = 'regular';
    for j = 1:countReg_runs
        clear events scan V matlabbatch out names possible_names indx task_idx preproc_data_dir run_volumes SPM out_dir confounds_table
        idx_reg = [];
        run_dx = ['reg_','run-',num2str(j)]; run = ['run-0',num2str(j)]; % run "name"
        [idx_reg,~] = find(contains(allRegRewFiles,string(current_id)) & contains(allRegRewFiles,run_dx));
        run_files = allRegRewFiles(idx_reg);

        % open events files of that run, that participant! and populate the events struct ................
        % these are the event files saved with Monica's script!
        events = struct; scan = struct;
        for k = 1:numel(run_files)
            current_run_file = fullfile(events_dir,run_files{k});
            % extract type of event from filename
            match = regexp(current_run_file, 'run-\d+_(.*?)\.txt', 'tokens');
            event_type = match{1}{1};
            event_data = readmatrix(current_run_file);
            switch event_type
                case 'go-hit'
                    if ~isempty(event_data) && any(event_data(:) ~= 0) % if data is not all zeros!
                        events.go_hit_onsets = event_data(:,1);
                        events.go_hit_durations = event_data(:,2);
                    end
                case 'go-incorr'
                    if ~isempty(event_data) && any(event_data(:) ~= 0)
                        events.go_incorr_onsets = event_data(:,1);
                        events.go_incorr_durations = event_data(:,2);
                    end
                case 'go-omit'
                    if ~isempty(event_data) && any(event_data(:) ~= 0)
                        events.go_omit_onsets = event_data(:,1);
                        events.go_omit_durations = event_data(:,2);
                    end
                case 'nogo-comerr'
                    if ~isempty(event_data) && any(event_data(:) ~= 0)
                        events.nogo_comerr_onsets = event_data(:,1);
                        events.nogo_comerr_durations = event_data(:,2);
                    end
                case 'nogo-corr'
                    if ~isempty(event_data) && any(event_data(:) ~= 0)
                        events.nogo_corr_onsets = event_data(:,1);
                        events.nogo_corr_durations = event_data(:,2);
                    end
                case 'prem'
                    if ~isempty(event_data) && any(event_data(:) ~= 0)
                        events.pre_resp_onsets = event_data(:,1);
                        events.pre_resp_durations = event_data(:,2);
                    end
            end
        end

        % get scan data and populate that struct too.....
        preproc_data_dir = fullfile(gngreg_data_dir,['sub-',num2str(current_id),'W01'],'ses-01','func');
        list = dir(preproc_data_dir);
        names = {list.name}';
        % filter filenames by task, and try to deal with any possible run filename
        indx = j; HRF = 'before'; task = 'regular';
        if indx == 1 && contains(task,'regular') % gngregular1
            task_idx = contains(names,'gngregular01')|...
                contains(names,'gngregular_run-01')|...
                contains(names,'gngregular01_run-01')|...
                contains(names,'gngregular01_run-02');
        elseif indx == 2 && contains(task,'regular') % gngregular2
            task_idx = contains(names,'gngregular02')|...
                contains(names,'gngregular_run-02')|...
                contains(names,'gngregular02_run-01')|...
                contains(names,'gngregular02_run-02');
        end
        possible_names = names(task_idx);
        % from those, select the one that has "desc-preproc_bold.nii"
        [row,~] = find(contains(possible_names,'desc-preproc_bold.nii') & ~contains(possible_names,'.mat'));
        if ~isempty(row) & (length(row) == 1) & ~contains(possible_names(row),'.gz') % if 1 file and compressed
            preproc_file = possible_names{row};
            V = spm_vol(fullfile(preproc_data_dir,preproc_file));
            TR = V(1).private.timing.tspace;
            nvols = size(V,1);
            t  = (0:nvols-1)' * TR;
        elseif ~isempty(row) & (length(row) > 1) % if 2 files
            possible_names = possible_names(row,1);
            [row,~] = find(contains(possible_names,'desc-preproc_bold.nii') & ~contains(possible_names,'.gz')); % find the decompressed one!
            % at this point, I expect only one element, the one decompressed before
            preproc_file = possible_names{row};
            V = spm_vol(fullfile(preproc_data_dir,preproc_file));
            TR = V(1).private.timing.tspace;
            nvols = size(V,1);
            t  = (0:nvols-1)' * TR;
        elseif ~isempty(row) & (length(row) == 1) & contains(possible_names(row),'.gz')
            % at this point, there is only a .gz file and I decompress it here
            preproc_file = possible_names{row};
            gunzip(fullfile(preproc_data_dir,preproc_file));
            unzipped_preproc_file = erase(preproc_file,".gz");
            V = spm_vol(fullfile(preproc_data_dir,preproc_file));
            TR = V(1).private.timing.tspace;
            nvols = size(V,1);
            t  = (0:nvols-1)' * TR;
            % delete(fullfile(preproc_data_dir,unzipped_preproc_file));
        elseif isempty(row)
            warning('--> No preprocessed .nii of %s Regular run %s', num2str(current_id), num2str(indx));
        end

        % save events data as a .mat file!
        if ~isempty(row)
            scan.timing = t;
            scan.TR = TR;
            scan.nvols = nvols;

            clean_name = regexprep(preproc_file, '\.nii(\.gz)?$', ''); % remove the .nii.gz
            events_filename = fullfile(preproc_data_dir,[clean_name,'_events.mat']);
            events_mat_filepath = events_filename;
            nii_filepath = preproc_file;
            save(events_filename,'events','scan','events_mat_filepath','nii_filepath');
            % plot spm-alike design
            close all;
            out = dg_plot_spm_like_design(scan, events, task, HRF);
            image_filename1 = fullfile(preproc_data_dir,[clean_name,'_eventsPlot.png']);
            saveas(gcf, image_filename1);
            clean_name = string(extractBefore(clean_name, '_space-'));
            fprintf('--Events and plot saved: %s \n',clean_name);
        end

        % now, run GLM with HRF in SPM...
        preproc_filename = fullfile(preproc_data_dir, preproc_file);
        if exist(preproc_filename)
            preproc_filename = regexprep(preproc_filename, '\.gz$', ''); % at this point, I expect all nii decompressed!
            run_volumes = spm_select('Expand',preproc_filename); % GET NUMBER OF VOLUMES AND ALL VOLUMES NAMES TO PASS IT TO SPM!
        else
            error('No decompressed .nii file for %s Reward run %s', num2str(current_id), num2str(indx));
        end
        spm1stLdir = fullfile(gngreg_data_dir, 'derivatives', 'spm');
        sub_id = regexp(preproc_data_dir, 'sub-[^/]+', 'match', 'once');
        % get correct task and run!
        sub_task_run = regexp(preproc_file, '(?<=ses-01_).*?(?=_space)', 'match', 'once');
        new_sub_dir = fullfile(spm1stLdir, sub_id, '1st_level', sub_task_run);

        % create dir to store SPM/GLM results
        if ~exist(new_sub_dir)
            mkdir(new_sub_dir);
        end

        matlabbatch{1}.spm.stats.fmri_spec.dir = {new_sub_dir};
        matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
        matlabbatch{1}.spm.stats.fmri_spec.timing.RT = 0.8;
        matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
        matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;
        matlabbatch{1}.spm.stats.fmri_spec.sess.scans = cellstr(run_volumes);

        % --- Build sess.cond dynamically and SKIP empty conditions ---
        % start clean
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond = struct('name', {}, 'onset', {},...
            'duration', {},'tmod', {},...
            'pmod', {}, 'orth', {});
        % --- Build conditions (skip if onset OR duration is empty) ---
        fn = fieldnames(events);
        % Keep only the *_onsets fields
        is_onset = endsWith(fn, '_onsets');
        onset_fields = fn(is_onset);
        k = 0;  % condition counter (only increments when condition is kept)
        for ii = 1:numel(onset_fields)
            onset_fn = onset_fields{ii};                        % e.g., 'go_hit_onsets'
            cond_name = erase(onset_fn, '_onsets');             % e.g., 'go_hit'
            dur_fn = [cond_name '_durations'];                  % e.g., 'go_hit_durations'

            % Pull and force column vectors
            on  = events.(onset_fn);  on  = on(:);
            dur = events.(dur_fn);    dur = dur(:);

            % Sanity check
            if numel(on) ~= numel(dur)
                error('Condition "%s": #onsets (%d) ~= #durations (%d).', ...
                    cond_name, numel(on), numel(dur));
            end

            % Remove NaNs if needed
            keep = ~isnan(on) & ~isnan(dur);
            on = on(keep);
            dur = dur(keep);

            % Add condition to batch
            k = k + 1;
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(k).name = cond_name;
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(k).onset = on;
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(k).duration = dur;
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(k).tmod = 0;
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(k).pmod = struct('name', {}, 'param', {}, 'poly', {});
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(k).orth = 1;
        end

        matlabbatch{1}.spm.stats.fmri_spec.sess.multi = {''};
        matlabbatch{1}.spm.stats.fmri_spec.sess.regress = struct('name', {}, 'val', {});
        matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {''};
        matlabbatch{1}.spm.stats.fmri_spec.sess.hpf = 128;
        matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
        matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0]; % no time derivative and no dispersion derivative
        matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
        matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
        matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.8;
        matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
        matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
        matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
        matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
        matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

        spm_jobman('run', matlabbatch);

        % save plot of postHRF events design matrix!
        %     but first, load SPM.mat struct that has the HRF-convolved events data/design matrix
        spm_matfilename = fullfile(spm1stLdir,sub_id,'1st_level',sub_task_run,'SPM.mat');
        clear matlabbatch;
        load(spm_matfilename); % load the SPM.mat file!
        if ~isempty(preproc_filename) && ~isempty(events) && ~isempty(scan) && ~isempty(SPM)
            close all;
            clean_name = regexprep(events_filename, '\_events.mat$', ''); % remove the .nii.gz
            out = dg_plot_hrf_design(SPM, events, scan, task); % with SPM struct!
            image_filename = [clean_name,'_HRFdesignPlot.png'];
            saveas(gcf, image_filename);
            clean_name = string(extractBefore(clean_name, '_space-'));
            fprintf('--Events and plot saved: %s \n',clean_name);
        else
            warning('--> No preproc nii nor events .mat file found for %s Regular run %s', num2str(current_id), num2str(indx));
        end

        % now, APPEND fMRIPrep confounds with task-convolved confounds/regressors
        %     first, load fMRIPrep confounds .tsv file
        list = dir(preproc_data_dir);
        names = {list.name}';
        [row2,~] = find(contains(names,'desc-confounds_timeseries.tsv') & contains(names, sub_task_run));
        if ~isempty(row2) & (length(row2) == 1)
            fmriprep_confounds_file = names{row2};
            run_confounds = readtable(fullfile(preproc_data_dir,fmriprep_confounds_file),'FileType', 'text', 'Delimiter', '\t');
        else
            fmriprep_confounds_file = '-';
        end
        %     then, get task events names in the same order they were convolved
        %     take them from the SPM.mat file :)
        events_names = {SPM.Sess.Fc.name}; % tasks events names saved by SPM
        convolved_events = SPM.xX.X; % design matrix from SPM
        % WAIT! REMOVE THE CONSTANT TERM / INTERCEPT ARRAY, XCP-D DOESN'T NEED IT!
        % SPM adds a column of ones to the GLM design because that column is the
        % intercept, it allows to estimate the mean signal level independently of the task regressors...
        % .....usually it's the last line... but I do this in a more flexible/reliable way
        [~, col_all_ones] = find(all(convolved_events == 1));
        convolved_events(:,col_all_ones) = [];
        % APPEND HRF-CONVOLVED REGRESSORS WITH FMRIPREP WM AND CSF CONFOUNDS BETWEEN THIS SECTION AND THE FOLLOWING
        % KEY NOTE: THE 36P indication for the "--nuisance_regressors" parameter includes GSR!!!
        % THUS I FORCE A 32P OPTION + TASK REGRESSORS
        % MORE INFO: https://xcp-d.readthedocs.io/en/latest/workflows.html#confound-regressor-selection
        fmriprep_names = { ...
            'rot_x'
            'rot_x_derivative1'
            'rot_x_derivative1_power2'
            'rot_x_power2'
            'rot_y'
            'rot_y_derivative1'
            'rot_y_derivative1_power2'
            'rot_y_power2'
            'rot_z'
            'rot_z_derivative1'
            'rot_z_derivative1_power2'
            'rot_z_power2'
            'trans_x'
            'trans_x_derivative1'
            'trans_x_derivative1_power2'
            'trans_x_power2'
            'trans_y'
            'trans_y_derivative1'
            'trans_y_derivative1_power2'
            'trans_y_power2'
            'trans_z'
            'trans_z_derivative1'
            'trans_z_derivative1_power2'
            'trans_z_power2'
            'csf'
            'csf_derivative1'
            'csf_derivative1_power2'
            'csf_power2'
            'white_matter'
            'white_matter_derivative1'
            'white_matter_derivative1_power2'
            'white_matter_power2'};

        missing_cols = setdiff(fmriprep_names, run_confounds.Properties.VariableNames);
        if ~isempty(missing_cols)
            error('Missing fMRIPrep confound columns:\n%s', strjoin(missing_cols, '\n'));
        end
        task_table = array2table(convolved_events);
        task_table.Properties.VariableNames = events_names;
        confounds_table = [run_confounds(:, fmriprep_names), task_table];

        % % just to check / visualize
        % mu = mean(confounds_table, 1, 'omitnan');
        % sigma = std(confounds_table, 0, 1, 'omitnan'); % normalize
        % because numeric ranges across regressors are too heterogeneous and makes plot colormap unuseful
        % Xz = (confounds_table - mu) ./ sigma;
        % close all;
        % figure; imagesc(Xz); colormap copper; colorbar; set(gca, 'Color', [0.8 0.8 0.8]);
        % SAVE COVARIATES AS .TSV FILE HERE (AND AS A .MAT FILE TOO)
        custom_confounds_dir = custom_confounds_root_reg;
        out_dir = fullfile(custom_confounds_dir, sub_id, 'ses-01', 'func'); % build BIDS-like output directory: sub-XX/ses-01/func
        if ~exist(custom_confounds_dir, 'dir')
            mkdir(custom_confounds_dir);
        end
        if ~exist(out_dir, 'dir')
            mkdir(out_dir);
        end

        % a json file that makes the TSV contain only the task regressors, and the custom dataset will be recognizable to XCP-D
        custom_desc_file = fullfile(custom_confounds_dir, 'dataset_description.json');
        fid = fopen(custom_desc_file, 'w');
        fprintf(fid, '{\n');
        fprintf(fid, '  "Name": "Custom confounds dataset",\n');
        fprintf(fid, '  "BIDSVersion": "1.9.0",\n');
        fprintf(fid, '  "DatasetType": "derivatives",\n');
        fprintf(fid, '  "GeneratedBy": [\n');
        fprintf(fid, '    {\n');
        fprintf(fid, '      "Name": "MATLAB custom confounds workflow"\n');
        fprintf(fid, '    }\n');
        fprintf(fid, '  ]\n');
        fprintf(fid, '}\n');
        fclose(fid);

        out_tsv = fullfile(out_dir, fmriprep_confounds_file); % use the SAME filename as fMRIPrep confounds
        writetable(confounds_table, out_tsv, 'FileType', 'text', 'Delimiter', '\t', 'WriteVariableNames', true); % write TSV
        fprintf('Saved XCP-D confounds TSV:\n%s\n', out_tsv); % report :)
        [filepath, out_name, ext] = fileparts(out_tsv);
        out_mat = fullfile(filepath,[out_name,'.mat']);
        save(out_mat,'confounds_table');

        % OPTIONAL FOR XCP-D: creation of a .json file with the same name as the .tsv
        % Build JSON filename from TSV
        [~, baseName, ~] = fileparts(out_tsv);
        out_json = fullfile(out_dir, [baseName, '.json']);
        % Define JSON contents
        json_struct = struct('Description', 'Confounds generated after SPM HRF.');
        % Encode JSON (pretty print if supported)
        try
            json_text = jsonencode(json_struct, 'PrettyPrint', true);
        catch
            json_text = jsonencode(json_struct); % fallback for older MATLAB
        end
        % Write JSON to disk
        fid = fopen(out_json, 'w');
        assert(fid ~= -1, 'Could not create JSON file: %s', out_json);
        fwrite(fid, json_text, 'char');
        fclose(fid);
        fprintf('Saved XCP-D confounds JSON:\n%s\n', out_json);

        % I must do scrubbing AFTER postproc. USE LATER THE ARRAYS I SAVED WITH "dg_BMAP_volumes_analysis_perTask.m"
        % as example, after running XCP-D there are several outputs, among them there's a file
        % "sub-2002W01_ses-01_task-gngregular_run-01_space-MNI152NLin2009cAsym_seg-Seitzman_stat-mean_timeseries.tsv"
        % this file has a 487x300 table, and I have to remove the ROWS/VOLS to scrub from it
        % (the actual scrubbing we'll do), to later recalculate a FC matrix!!!!!!!!!!!!!

    end

    clear events scan task V matlabbatch out names possible_names indx task_idx preproc_data_dir run_volumes SPM out_dir confounds_table
    %% GNG REWARD!!
    countRew_runs = CheckConcurrence.countRew_runsSPd(i,1);  task = 'reward';
    for j = 1:countRew_runs
        clear events scan V matlabbatch out names possible_names indx task_idx preproc_data_dir run_volumes SPM out_dir confounds_table
        idx_rew = [];
        if current_id == 2170 && j == 2 % this participant did not have a "gngrew 02"
            run_dx = 'rew_run-3';
        else
            run_dx = ['rew_','run-',num2str(j)]; run = ['run-0',num2str(j)];
        end
        [idx_rew,~] = find(contains(allRegRewFiles,string(current_id)) & ...
            contains(allRegRewFiles,run_dx) & contains(allRegRewFiles,'rew'));
        run_files = allRegRewFiles(idx_rew);

        % open events files of that run, that participant! and populate the events struct ................
        % these are the event files saved with Monica's script!
        events = struct; scan = struct;
        for k = 1:numel(run_files)
            current_run_file = fullfile(events_dir,run_files{k});
            % extract type of event from filename
            match = regexp(current_run_file, 'run-\d+_(.*?)\.txt', 'tokens');
            event_type = match{1}{1};
            event_data = readmatrix(current_run_file);
            switch event_type
                % shared events with Regular
                case 'go-hit'
                    if ~isempty(event_data) && any(event_data(:) ~= 0) % if data is not all zeros!
                        events.go_hit_onsets = event_data(:,1);
                        events.go_hit_durations = event_data(:,2);
                    end
                case 'go-incorr'
                    if ~isempty(event_data) && any(event_data(:) ~= 0)
                        events.go_incorr_onsets = event_data(:,1);
                        events.go_incorr_durations = event_data(:,2);
                    end
                case 'go-omit'
                    if ~isempty(event_data) && any(event_data(:) ~= 0)
                        events.go_omit_onsets = event_data(:,1);
                        events.go_omit_durations = event_data(:,2);
                    end
                case 'nogo-comerr'
                    if ~isempty(event_data) && any(event_data(:) ~= 0)
                        events.nogo_comerr_onsets = event_data(:,1);
                        events.nogo_comerr_durations = event_data(:,2);
                    end
                case 'nogo-corr'
                    if ~isempty(event_data) && any(event_data(:) ~= 0)
                        events.nogo_corr_onsets = event_data(:,1);
                        events.nogo_corr_durations = event_data(:,2);
                    end
                case 'prem'
                    if ~isempty(event_data) && any(event_data(:) ~= 0)
                        events.pre_resp_onsets = event_data(:,1);
                        events.pre_resp_durations = event_data(:,2);
                    end
                    % events exclusive from Reward!
                case 'go-hit_rew'
                    if ~isempty(event_data) && any(event_data(:) ~= 0) % if data is not all zeros!
                        events.go_hit_rew_onsets = event_data(:,1);
                        events.go_hit_rew_durations = event_data(:,2);
                    end
                case 'go-hit_norew'
                    if ~isempty(event_data) && any(event_data(:) ~= 0)
                        events.go_hit_norew_onsets = event_data(:,1);
                        events.go_hit_norew_durations = event_data(:,2);
                    end
                case 'go-omit_norew'
                    if ~isempty(event_data) && any(event_data(:) ~= 0)
                        events.go_omit_norew_onsets = event_data(:,1);
                        events.go_omit_norew_durations = event_data(:,2);
                    end
                case 'nogo-comerr_norew'
                    if ~isempty(event_data) && any(event_data(:) ~= 0)
                        events.nogo_comerr_norew_onsets = event_data(:,1);
                        events.nogo_comerr_norew_durations = event_data(:,2);
                    end
                case 'nogo-corr_rew'
                    if ~isempty(event_data) && any(event_data(:) ~= 0)
                        events.nogo_corr_rew_onsets = event_data(:,1);
                        events.nogo_corr_rew_durations = event_data(:,2);
                    end
                case 'prem_rew'
                    if ~isempty(event_data) && any(event_data(:) ~= 0)
                        events.prem_rew_onsets = event_data(:,1);
                        events.prem_rew_durations = event_data(:,2);
                    end
                case 'prem_norew'
                    if ~isempty(event_data) && any(event_data(:) ~= 0)
                        events.prem_norew_onsets = event_data(:,1);
                        events.prem_norew_durations = event_data(:,2);
                    end
                case 'go-incorr_norew'
                    if ~isempty(event_data) && any(event_data(:) ~= 0)
                        events.go_incorr_norew_onsets = event_data(:,1);
                        events.go_incorr_norew_durations = event_data(:,2);
                    end
            end
        end

        % get scan data and populate that struct too.....
        preproc_data_dir = fullfile(gngrew_data_dir,['sub-',num2str(current_id),'W01'],'ses-01','func');
        list = dir(preproc_data_dir);
        names = {list.name}';
        % filter filenames by task, and try to deal with any possible run filename
        if current_id == 2170 && j == 2 % this participant did not have a "gngrew 02"
            indx = 3;
        else
            indx = j;
        end
        task = 'reward'; HRF = 'before';
        if indx == 1 && contains(task,'reward') % gngreward1
            task_idx = contains(names,'gngreward01')|...
                contains(names,'gngreward_run-01')|...
                contains(names,'gngreward01_run-01')|...
                contains(names,'gngreward01_run-02');
        elseif indx == 2 && contains(task,'reward') % gngreward2
            task_idx = contains(names,'gngreward02')|...
                contains(names,'gngreward_run-02')|...
                contains(names,'gngreward02_run-01')|...
                contains(names,'gngreward02_run-02');
        elseif indx == 3 && contains(task,'reward') % gngreward3
            task_idx = contains(names,'gngreward03')|...
                contains(names,'gngreward_run-03')|...
                contains(names,'gngreward02_run-01')|...
                contains(names,'gngreward02_run-02');
        elseif indx == 4 && contains(task,'reward') % gngreward4
            task_idx = contains(names,'gngreward04')|...
                contains(names,'gngreward_run-04')|...
                contains(names,'gngreward04_run-01')|...
                contains(names,'gngreward04_run-02');
        end
        possible_names = names(task_idx);
        % from those, select the one that has "desc-preproc_bold.nii"
        [row,~] = find(contains(possible_names,'desc-preproc_bold.nii') & ~contains(possible_names,'.mat'));
        if ~isempty(row) & (length(row) == 1) & ~contains(possible_names(row),'.gz') % if 1 file and compressed
            preproc_file = possible_names{row};
            V = spm_vol(fullfile(preproc_data_dir,preproc_file));
            TR = V(1).private.timing.tspace;
            nvols = size(V,1);
            t  = (0:nvols-1)' * TR;
        elseif ~isempty(row) & (length(row) > 1) % if 2 files
            possible_names = possible_names(row,1);
            [row,~] = find(contains(possible_names,'desc-preproc_bold.nii') & ~contains(possible_names,'.gz')); % find the decompressed one!
            % at this point, I expect only one element, the one decompressed before
            preproc_file = possible_names{row};
            V = spm_vol(fullfile(preproc_data_dir,preproc_file));
            TR = V(1).private.timing.tspace;
            nvols = size(V,1);
            t  = (0:nvols-1)' * TR;
        elseif ~isempty(row) & (length(row) == 1) & contains(possible_names(row),'.gz')
            % at this point, there is only a .gz file and I decompress it here
            preproc_file = possible_names{row};
            gunzip(fullfile(preproc_data_dir,preproc_file));
            unzipped_preproc_file = erase(preproc_file,".gz");
            V = spm_vol(fullfile(preproc_data_dir,preproc_file));
            TR = V(1).private.timing.tspace;
            nvols = size(V,1);
            t  = (0:nvols-1)' * TR;
            % delete(fullfile(preproc_data_dir,unzipped_preproc_file));
        elseif isempty(row)
            warning('--> No preprocessed .nii of %s Reward run %s', num2str(current_id), num2str(indx));
        end

        % save events data as a .mat file!
        if ~isempty(row)
            scan.timing = t;
            scan.TR = TR;
            scan.nvols = nvols;

            clean_name = regexprep(preproc_file, '\.nii(\.gz)?$', ''); % remove the .nii.gz
            events_filename = fullfile(preproc_data_dir,[clean_name,'_events.mat']);
            events_mat_filepath = events_filename;
            nii_filepath = preproc_file;
            save(events_filename,'events','scan','events_mat_filepath','nii_filepath');
            % plot spm-alike design
            close all;
            out = dg_plot_spm_like_design(scan, events, task, HRF);
            image_filename1 = fullfile(preproc_data_dir,[clean_name,'_eventsPlot.png']);
            saveas(gcf, image_filename1);
            clean_name = string(extractBefore(clean_name, '_space-'));
            fprintf('--Events and plot saved: %s \n',clean_name);
        end

        % now, run GLM with HRF in SPM...
        preproc_filename = fullfile(preproc_data_dir, preproc_file);
        if exist(preproc_filename)
            preproc_filename = regexprep(preproc_filename, '\.gz$', ''); % at this point, I expect all nii decompressed!
            run_volumes = spm_select('Expand',preproc_filename); % GET NUMBER OF VOLUMES AND ALL VOLUMES NAMES TO PASS IT TO SPM!
        else
            error('No decompressed .nii file for %s Reward run %s', num2str(current_id), num2str(indx));
        end
        spm1stLdir = fullfile(gngrew_data_dir, 'derivatives', 'spm');
        sub_id = regexp(preproc_data_dir, 'sub-[^/]+', 'match', 'once');
        % get correct task and run!
        sub_task_run = regexp(preproc_file, '(?<=ses-01_).*?(?=_space)', 'match', 'once');
        new_sub_dir = fullfile(spm1stLdir, sub_id, '1st_level', sub_task_run);

        % create dir to store SPM/GLM results
        if ~exist(new_sub_dir)
            mkdir(new_sub_dir);
        end

        matlabbatch{1}.spm.stats.fmri_spec.dir = {new_sub_dir};
        matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
        matlabbatch{1}.spm.stats.fmri_spec.timing.RT = 0.8;
        matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
        matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;
        matlabbatch{1}.spm.stats.fmri_spec.sess.scans = cellstr(run_volumes);

        % --- Build sess.cond dynamically and SKIP empty conditions ---
        % start clean
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond = struct('name', {}, 'onset', {},...
            'duration', {},'tmod', {},...
            'pmod', {}, 'orth', {});
        % --- Build conditions (skip if onset OR duration is empty) ---
        fn = fieldnames(events);
        % Keep only the *_onsets fields
        is_onset = endsWith(fn, '_onsets');
        onset_fields = fn(is_onset);
        k = 0;  % condition counter (only increments when condition is kept)
        for ii = 1:numel(onset_fields)
            onset_fn = onset_fields{ii};                        % e.g., 'go_hit_onsets'
            cond_name = erase(onset_fn, '_onsets');             % e.g., 'go_hit'
            dur_fn = [cond_name '_durations'];                  % e.g., 'go_hit_durations'

            % Pull and force column vectors
            on  = events.(onset_fn);  on  = on(:);
            dur = events.(dur_fn);    dur = dur(:);

            % Sanity check
            if numel(on) ~= numel(dur)
                error('Condition "%s": #onsets (%d) ~= #durations (%d).', ...
                    cond_name, numel(on), numel(dur));
            end

            % Remove NaNs if needed
            keep = ~isnan(on) & ~isnan(dur);
            on = on(keep);
            dur = dur(keep);

            % Add condition to batch
            k = k + 1;
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(k).name = cond_name;
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(k).onset = on;
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(k).duration = dur;
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(k).tmod = 0;
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(k).pmod = struct('name', {}, 'param', {}, 'poly', {});
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(k).orth = 1;
        end

        matlabbatch{1}.spm.stats.fmri_spec.sess.multi = {''};
        matlabbatch{1}.spm.stats.fmri_spec.sess.regress = struct('name', {}, 'val', {});
        matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {''};
        matlabbatch{1}.spm.stats.fmri_spec.sess.hpf = 128;
        matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
        matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0]; % no time derivative and no dispersion derivative
        matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
        matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
        matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.8;
        matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
        matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
        matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
        matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
        matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

        spm_jobman('run', matlabbatch);

        % save plot of postHRF events design matrix!
        %     but first, load SPM.mat struct that has the HRF-convolved events data/design matrix
        spm_matfilename = fullfile(spm1stLdir,sub_id,'1st_level',sub_task_run,'SPM.mat');
        clear matlabbatch;
        load(spm_matfilename); % load the SPM.mat file!
        if ~isempty(preproc_filename) && ~isempty(events) && ~isempty(scan) && ~isempty(SPM)
            close all;
            clean_name = regexprep(events_filename, '\_events.mat$', ''); % remove the .nii.gz
            out = dg_plot_hrf_design(SPM, events, scan, task); % with SPM struct!
            image_filename = [clean_name,'_HRFdesignPlot.png'];
            saveas(gcf, image_filename);
            clean_name = string(extractBefore(clean_name, '_space-'));
            fprintf('--Events and plot saved: %s \n',clean_name);
        else
            warning('--> No preproc nii nor events .mat file found for %s Regular run %s', num2str(current_id), num2str(indx));
        end

        % now, APPEND fMRIPrep confounds with task-convolved confounds/regressors
        %     first, load fMRIPrep confounds .tsv file
        list = dir(preproc_data_dir);
        names = {list.name}';
        [row2,~] = find(contains(names,'desc-confounds_timeseries.tsv') & contains(names, sub_task_run));
        if ~isempty(row2) & (length(row2) == 1)
            fmriprep_confounds_file = names{row2};
            run_confounds = readtable(fullfile(preproc_data_dir,fmriprep_confounds_file),'FileType', 'text', 'Delimiter', '\t');
        else
            fmriprep_confounds_file = '-';
        end

        %     then, get task events names in the same order they were convolved
        %     take them from the SPM.mat file :)
        events_names = {SPM.Sess.Fc.name}; % tasks events names saved by SPM
        convolved_events = SPM.xX.X; % design matrix from SPM
        % WAIT! REMOVE THE CONSTANT TERM / INTERCEPT ARRAY, XCP-D DOESN'T NEED IT!
        % SPM adds a column of ones to the GLM design because that column is the
        % intercept, it allows to estimate the mean signal level independently of the task regressors...
        % .....usually it's the last line... but I do this in a more flexible/reliable way
        [~, col_all_ones] = find(all(convolved_events == 1));
        convolved_events(:,col_all_ones) = [];
        % APPEND HRF-CONVOLVED REGRESSORS WITH FMRIPREP WM AND CSF CONFOUNDS BETWEEN THIS SECTION AND THE FOLLOWING
        % KEY NOTE: THE 36P indication for the "--nuisance_regressors" parameter includes GSR!!!
        % THUS I FORCE A 32P OPTION + TASK REGRESSORS
        % MORE INFO: https://xcp-d.readthedocs.io/en/latest/workflows.html#confound-regressor-selection
        fmriprep_names = { ...
            'rot_x'
            'rot_x_derivative1'
            'rot_x_derivative1_power2'
            'rot_x_power2'
            'rot_y'
            'rot_y_derivative1'
            'rot_y_derivative1_power2'
            'rot_y_power2'
            'rot_z'
            'rot_z_derivative1'
            'rot_z_derivative1_power2'
            'rot_z_power2'
            'trans_x'
            'trans_x_derivative1'
            'trans_x_derivative1_power2'
            'trans_x_power2'
            'trans_y'
            'trans_y_derivative1'
            'trans_y_derivative1_power2'
            'trans_y_power2'
            'trans_z'
            'trans_z_derivative1'
            'trans_z_derivative1_power2'
            'trans_z_power2'
            'csf'
            'csf_derivative1'
            'csf_derivative1_power2'
            'csf_power2'
            'white_matter'
            'white_matter_derivative1'
            'white_matter_derivative1_power2'
            'white_matter_power2'};

        missing_cols = setdiff(fmriprep_names, run_confounds.Properties.VariableNames);
        if ~isempty(missing_cols)
            error('Missing fMRIPrep confound columns:\n%s', strjoin(missing_cols, '\n'));
        end

        task_table = array2table(convolved_events);
        task_table.Properties.VariableNames = events_names;

        confounds_table = [run_confounds(:, fmriprep_names), task_table];

        % % just to check / visualize
        % mu = mean(confounds_table, 1, 'omitnan');
        % sigma = std(confounds_table, 0, 1, 'omitnan'); % normalize
        % because numeric ranges across regressors are too heterogeneous and makes plot colormap unuseful
        % Xz = (confounds_table - mu) ./ sigma;
        % close all;
        % figure; imagesc(Xz); colormap copper; colorbar; set(gca, 'Color', [0.8 0.8 0.8]);

        % SAVE COVARIATES AS .TSV FILE HERE (AND AS A .MAT FILE TOO)
        custom_confounds_dir = custom_confounds_root_rew;
        out_dir = fullfile(custom_confounds_dir, sub_id, 'ses-01', 'func'); % build BIDS-like output directory: sub-XX/ses-01/func
        if ~exist(custom_confounds_dir, 'dir')
            mkdir(custom_confounds_dir);
        end
        if ~exist(out_dir, 'dir')
            mkdir(out_dir);
        end

        % a json file that makes the TSV contain only the task regressors, and the custom dataset will be recognizable to XCP-D
        custom_desc_file = fullfile(custom_confounds_dir, 'dataset_description.json');
        fid = fopen(custom_desc_file, 'w');
        fprintf(fid, '{\n');
        fprintf(fid, '  "Name": "Custom confounds dataset",\n');
        fprintf(fid, '  "BIDSVersion": "1.9.0",\n');
        fprintf(fid, '  "DatasetType": "derivatives",\n');
        fprintf(fid, '  "GeneratedBy": [\n');
        fprintf(fid, '    {\n');
        fprintf(fid, '      "Name": "MATLAB custom confounds workflow"\n');
        fprintf(fid, '    }\n');
        fprintf(fid, '  ]\n');
        fprintf(fid, '}\n');
        fclose(fid);

        out_tsv = fullfile(out_dir, fmriprep_confounds_file); % use the SAME filename as fMRIPrep confounds
        writetable(confounds_table, out_tsv, 'FileType', 'text', 'Delimiter', '\t', 'WriteVariableNames', true); % write TSV
        fprintf('Saved XCP-D confounds TSV:\n%s\n', out_tsv); % report :)
        [filepath, out_name, ext] = fileparts(out_tsv);
        out_mat = fullfile(filepath,[out_name,'.mat']);
        save(out_mat,'confounds_table');

        % OPTIONAL FOR XCP-D: creation of a .json file with the same name as the .tsv
        % Build JSON filename from TSV
        [~, baseName, ~] = fileparts(out_tsv);
        out_json = fullfile(out_dir, [baseName, '.json']);
        % Define JSON contents
        json_struct = struct('Description', 'Confounds generated after SPM-HRF.');
        % Encode JSON (pretty print if supported)
        try
            json_text = jsonencode(json_struct, 'PrettyPrint', true);
        catch
            json_text = jsonencode(json_struct); % fallback for older MATLAB
        end
        % Write JSON to disk
        fid = fopen(out_json, 'w');
        assert(fid ~= -1, 'Could not create JSON file: %s', out_json);
        fwrite(fid, json_text, 'char');
        fclose(fid);
        fprintf('Saved XCP-D confounds JSON:\n%s\n', out_json);

        % I must do scrubbing AFTER postproc. USE LATER THE ARRAYS I SAVED WITH "dg_BMAP_volumes_analysis_perTask.m"
        % as example, after running XCP-D there are several outputs, among them there's a file
        % "sub-2002W01_ses-01_task-gngregular_run-01_space-MNI152NLin2009cAsym_seg-Seitzman_stat-mean_timeseries.tsv"
        % this file has a 487x300 table, and I have to remove the ROWS/VOLS to scrub from it
        % (the actual scrubbing we'll do), to later recalculate a FC matrix!!!!!!!!!!!!!

    end

end

%% 3. Edit/customize a .yaml file for each run (1participant) and run XCP-D per run (1participant)
clear events scan task V matlabbatch out names possible_names indx task_idx preproc_data_dir run_volumes SPM out_dir confounds_table
close all;

xcpd_work_dir = '/users/d/g/dga/BrainMAP/xcpd_work/';
gngreg_output_dir = '/users/d/g/dga/BrainMAP/xcpd_work/data_gngregular_xcpd/';
gngrew_output_dir = '/users/d/g/dga/BrainMAP/xcpd_work/data_gngreward_xcpd/';
XCP_DIR = '/work/users/d/g/dga/tools/XCP-D';
SIMG = fullfile(XCP_DIR, 'XCP-D-0.10.7.simg');
custom_confounds_root_reg = fullfile(xcpd_work_dir, 'data_gngregular_custom_confounds');

% define Seitzman as custom atlas and a .json to make XCP-D to find it also in its custom directory
% (not in the directory XCP-D creates for default atlases)
ATLAS_DIR = '/users/d/g/dga/BrainMAP/xcpd_work/seitzman_atlas_dataset';
SEITZMAN_SRC_DIR = '/users/d/g/dga/BrainMAP/xcpd_work/data_gngregular_xcpd/atlases/atlas-Seitzman';
SEITZMAN_DST_DIR = fullfile(ATLAS_DIR, 'atlas-Seitzman');

if ~exist(ATLAS_DIR, 'dir')
    mkdir(ATLAS_DIR);
end

if ~exist(SEITZMAN_DST_DIR, 'dir')
    mkdir(SEITZMAN_DST_DIR);
end

atlas_desc = struct();
atlas_desc.Name = 'Custom Seitzman atlas dataset';
atlas_desc.BIDSVersion = '1.9.0';
atlas_desc.DatasetType = 'atlas';

atlas_desc_file = fullfile(ATLAS_DIR, 'dataset_description.json');
fid = fopen(atlas_desc_file, 'w');
fprintf(fid, '%s', jsonencode(atlas_desc, 'PrettyPrint', true));
fclose(fid);

copyfile(fullfile(SEITZMAN_SRC_DIR, 'atlas-Seitzman_dseg.tsv'), SEITZMAN_DST_DIR);
copyfile(fullfile(SEITZMAN_SRC_DIR, 'atlas-Seitzman_space-MNI152NLin2009cAsym_dseg.json'), SEITZMAN_DST_DIR);
copyfile(fullfile(SEITZMAN_SRC_DIR, 'atlas-Seitzman_space-MNI152NLin2009cAsym_dseg.nii.gz'), SEITZMAN_DST_DIR);

for i = 1 % :numel(common_IDs) % for loop per subject
    %% XCP-D for GNG REGULAR runs!!
    current_id = common_IDs(i);
    countReg_runs = CheckConcurrence.countReg_runsSPd(i,1); task = 'regular'; task_id = 'gngregular';
    for j = 1:countReg_runs
        clear task_idx preproc_data_dir possible_names row1 row2 row3 sub_id sub_task_run new_sub_dir confounds_table;

        indx = j;
        custom_confounds_dir = custom_confounds_root_reg;
        sub_id = ['sub-',num2str(current_id),'W01'];
        run_id = ['0' num2str(indx)];
        session_id = '01';
        sub_custom_confounds_dir = fullfile(custom_confounds_dir,sub_id,'ses-01','func');
        list = dir(sub_custom_confounds_dir);
        names = {list.name}';

        % load only the mat file with custom confounds
        if indx == 1 && contains(task,'regular') % gngregular1
            task_match = contains(names,'gngregular01') | ...
                contains(names,'gngregular_run-01') | ...
                contains(names,'gngregular01_run-01') | ...
                contains(names,'gngregular01_run-02');
            ext_match = endsWith(names,'.mat');
            task_idx = task_match & ext_match;
        elseif indx == 2 && contains(task,'regular') % gngregular2
            task_match = contains(names,'gngregular02') | ...
                contains(names,'gngregular_run-02') | ...
                contains(names,'gngregular02_run-01') | ...
                contains(names,'gngregular02_run-02');
            ext_match = endsWith(names,'.mat');
            task_idx = task_match & ext_match;
        end
        matfile = names{task_idx};

        if ~isempty(matfile) & (size(matfile,1) == 1)
            custom_confounds_matfile = fullfile(sub_custom_confounds_dir,matfile);
            load(custom_confounds_matfile);
        end

        % get confounds name here
        confounds_names = confounds_table.Properties.VariableNames';

        % EDIT AND SAVE THE YAML FILE HERE (32P!! )
        participant_label = erase(sub_id, 'sub-');
        task_names = setdiff(confounds_names(:), fmriprep_names(:), 'stable');
        fmriprep_names = {
            'rot_x'
            'rot_x_derivative1'
            'rot_x_derivative1_power2'
            'rot_x_power2'
            'rot_y'
            'rot_y_derivative1'
            'rot_y_derivative1_power2'
            'rot_y_power2'
            'rot_z'
            'rot_z_derivative1'
            'rot_z_derivative1_power2'
            'rot_z_power2'
            'trans_x'
            'trans_x_derivative1'
            'trans_x_derivative1_power2'
            'trans_x_power2'
            'trans_y'
            'trans_y_derivative1'
            'trans_y_derivative1_power2'
            'trans_y_power2'
            'trans_z'
            'trans_z_derivative1'
            'trans_z_derivative1_power2'
            'trans_z_power2'
            'csf'
            'csf_derivative1'
            'csf_derivative1_power2'
            'csf_power2'
            'white_matter'
            'white_matter_derivative1'
            'white_matter_derivative1_power2'
            'white_matter_power2'};

        yaml_file = fullfile(xcpd_work_dir, sprintf('custom_config_%s_%s_run-%s.yaml', participant_label, task_id, run_id));
        dataset_fmriprep = 'custom';
        dataset_task = 'custom';

        fid = fopen(yaml_file, 'w');
        fprintf(fid, 'name: confounds_custom\n');
        fprintf(fid, 'description: |\n');
        fprintf(fid, '  Nuisance regression including fMRIPrep confounds and HRF-convolved task regressors.\n');
        fprintf(fid, 'confounds:\n');

        fprintf(fid, '  fmriprep_core:\n');
        fprintf(fid, '    dataset: %s\n', dataset_fmriprep);
        fprintf(fid, '    query:\n');
        fprintf(fid, '      space: null\n');
        fprintf(fid, '      cohort: null\n');
        fprintf(fid, '      res: null\n');
        fprintf(fid, '      den: null\n');
        fprintf(fid, '      desc: confounds\n');
        fprintf(fid, '      suffix: timeseries\n');
        fprintf(fid, '      extension: .tsv\n');
        fprintf(fid, '    columns:\n');
        for k = 1:numel(fmriprep_names)
            fprintf(fid, '      - %s\n', fmriprep_names{k});
        end

        fprintf(fid, '\n');
        fprintf(fid, '  task_hrf:\n');
        fprintf(fid, '    dataset: %s\n', dataset_task);
        fprintf(fid, '    query:\n');
        fprintf(fid, '      space: null\n');
        fprintf(fid, '      cohort: null\n');
        fprintf(fid, '      res: null\n');
        fprintf(fid, '      den: null\n');
        fprintf(fid, '      desc: confounds\n');
        fprintf(fid, '      suffix: timeseries\n');
        fprintf(fid, '      extension: .tsv\n');
        fprintf(fid, '    columns:\n');
        for k = 1:numel(task_names)
            fprintf(fid, '      - %s\n', task_names{k});
        end
        fclose(fid);

        % HERE WRITE .SH FILE THAT LAUNCHES CUSTOMIZED XCP-D PARAMETERS :)
        % DERIVED VARIABLES
        % [folder, name, ~] = fileparts(custom_confounds_matfile);
        % custom_confounds_tsvfile = fullfile(folder, [name '.tsv']);
        custom_dataset_root = custom_confounds_root_reg;
        custom_desc_file = fullfile(custom_confounds_dir, 'dataset_description.json');
        fid = fopen(custom_desc_file, 'w');
        fprintf(fid, '{\n');
        fprintf(fid, '  "Name": "Custom confounds dataset",\n');
        fprintf(fid, '  "BIDSVersion": "1.9.0",\n');
        fprintf(fid, '  "DatasetType": "derivatives",\n');
        fprintf(fid, '  "GeneratedBy": [\n');
        fprintf(fid, '    {\n');
        fprintf(fid, '      "Name": "MATLAB custom confounds workflow"\n');
        fprintf(fid, '    }\n');
        fprintf(fid, '  ]\n');
        fprintf(fid, '}\n');
        fclose(fid);

        work_dir = fullfile(gngreg_output_dir, sprintf('work_%s_%s_run-%s', participant_label, task_id, run_id));
        if ~exist(gngreg_output_dir, 'dir'); mkdir(gngreg_output_dir); end
        if ~exist(work_dir, 'dir'); mkdir(work_dir); end

        % CREATE BIDS FILTER JSON (NECESSARY TO MAKE XCP-D RUN FOR JUST 1 RUN !!!! ...)
        filter_struct = struct();
        filter_struct.bold = struct( ...
            'session', {{session_id}}, ...
            'task', {{task_id}}, ...
            'run', {{run_id}});
        filter_json_file = fullfile(gngreg_output_dir, sprintf('%s_%s_run-%s_filter.json', participant_label, task_id, run_id));
        fid = fopen(filter_json_file, 'w');
        fprintf(fid, '%s', jsonencode(filter_struct, 'PrettyPrint', true));
        fclose(fid);

        % WRITE LOCAL SHELL SCRIPT THAT RUNS XCP-D
        sh_file = fullfile(gngreg_output_dir, sprintf('run_xcpd_%s_%s_run-%s.sh', participant_label, task_id, run_id));
        script_lines = {
            '#!/usr/bin/env bash'
            'set -euo pipefail'
            ''
            'module purge'
            'module load apptainer'
            ''
            sprintf('XCP_DIR="%s"', XCP_DIR)
            'SIMG="${XCP_DIR}/XCP-D-0.10.7.simg"'
            sprintf('FMRIPREP_DIR="%s"', gngreg_data_dir)
            sprintf('OUTPUT_DIR="%s"', gngreg_output_dir)
            sprintf('ATLAS_DIR="%s"', ATLAS_DIR)
            sprintf('CUSTOM_DATASET="%s"', custom_dataset_root)
            sprintf('NUISANCE_YAML="%s"', yaml_file)
            sprintf('PART="%s"', participant_label)
            sprintf('TASK_ID="%s"', task_id)
            sprintf('SESSION_ID="%s"', session_id)
            sprintf('WORK_DIR="%s"', work_dir)
            sprintf('BIDS_FILTER="%s"', filter_json_file)
            ''
            'mkdir -p "${OUTPUT_DIR}" "${WORK_DIR}"'
            ''
            'NCPUS=2'
            'export OMP_NUM_THREADS="${NCPUS}"'
            'BINDPATHS="/users,/proj,/work"'
            ''
            'apptainer run --cleanenv -B "${BINDPATHS}" \'
            '  "${SIMG}" \'
            '  "${FMRIPREP_DIR}" \'
            '  "${OUTPUT_DIR}" \'
            '  participant \'
            '  --fs-license-file "/work/users/d/g/dga/tools/freesurfer/7.3.2/license.txt" \'
            '  --mode nichart \'
            '  --participant-label "${PART}" \'
            '  --session-id "${SESSION_ID}" \'
            '  --bids-filter-file "${BIDS_FILTER}" \'
            '  -d atlas="${ATLAS_DIR}" custom="${CUSTOM_DATASET}" \'
            '  --atlases Seitzman 4S456Parcels \'
            '  -p "${NUISANCE_YAML}" \'
            '  -t "${TASK_ID}" \'
            '  --file-format nifti \'
            '  --smoothing 0 \'
            '  --motion-filter-type lp \'
            '  --band-stop-min 6 \'
            '  --motion-filter-order 2 \'
            '  -f 0.1 \'
            '  --output-type interpolated \'
            '  --min-time 120 \'
            '  --lower-bpf 0.008 \'
            '  --upper-bpf 0.09 \'
            '  --min-coverage 0.5 \'
            '  --nprocs "${NCPUS}" \'
            '  --clean-workdir \'
            '  -w "${WORK_DIR}"'
            ''            };
        fid = fopen(sh_file, 'w');
        for k = 1:numel(script_lines)
            fprintf(fid, '%s\n', script_lines{k});
        end
        fclose(fid);
        system(sprintf('chmod +x "%s"', sh_file)); % make it executable in terminal :)

        % RUN THE SCRIPT
        command = sprintf('bash "%s"', sh_file);
        [status, cmdout] = system(command);
        disp(status)
        disp(cmdout)
        warning('-------------------------------------------------------------');
        warning('Postproc done with: %s, %s run-%s', sub_id, task_id, run_id);
    end

   %% XCP-D for GNG REWARD runs!!
    countRew_runs = CheckConcurrence.countRew_runsSPd(i,1); task = 'reward'; task_id = 'gngreward';
    for j = 1:countRew_runs
        clear task_idx preproc_data_dir possible_names row1 row2 row3 sub_id sub_task_run new_sub_dir confounds_table;
        indx = j;
        custom_confounds_dir = custom_confounds_root_rew;
        sub_id = ['sub-',num2str(current_id),'W01'];
        run_id = ['0' num2str(indx)];
        session_id = '01';
        sub_custom_confounds_dir = fullfile(custom_confounds_dir,sub_id,'ses-01','func');
        list = dir(sub_custom_confounds_dir);
        names = {list.name}';

        % load only the mat file with custom confounds        
        if indx == 1 && contains(task,'reward') % gngreward1
            task_match = contains(names,'gngreward01')|...
                contains(names,'gngreward_run-01')|...
                contains(names,'gngreward01_run-01')|...
                contains(names,'gngreward01_run-02');
            ext_match = endsWith(names,'.mat');
            task_idx = task_match & ext_match;
        elseif indx == 2 && contains(task,'reward') % gngreward2
            task_match = contains(names,'gngreward02')|...
                contains(names,'gngreward_run-02')|...
                contains(names,'gngreward02_run-01')|...
                contains(names,'gngreward02_run-02');
            ext_match = endsWith(names,'.mat');
            task_idx = task_match & ext_match;
        elseif indx == 3 && contains(task,'reward') % gngreward3
            task_match = contains(names,'gngreward03')|...
                contains(names,'gngreward_run-03')|...
                contains(names,'gngreward02_run-01')|...
                contains(names,'gngreward02_run-02');
            ext_match = endsWith(names,'.mat');
            task_idx = task_match & ext_match;
        elseif indx == 4 && contains(task,'reward') % gngreward4
            task_match = contains(names,'gngreward04')|...
                contains(names,'gngreward_run-04')|...
                contains(names,'gngreward04_run-01')|...
                contains(names,'gngreward04_run-02');
            ext_match = endsWith(names,'.mat');
            task_idx = task_match & ext_match;
        end
        matfile = names{task_idx};

        if ~isempty(matfile) & (size(matfile,1) == 1)
            custom_confounds_matfile = fullfile(sub_custom_confounds_dir,matfile);
            load(custom_confounds_matfile);
        end

        % get confounds name here
        confounds_names = confounds_table.Properties.VariableNames';

        % EDIT AND SAVE THE YAML FILE HERE (32P!! )
        participant_label = erase(sub_id, 'sub-');
        task_names = setdiff(confounds_names(:), fmriprep_names(:), 'stable');
        fmriprep_names = {
            'rot_x'
            'rot_x_derivative1'
            'rot_x_derivative1_power2'
            'rot_x_power2'
            'rot_y'
            'rot_y_derivative1'
            'rot_y_derivative1_power2'
            'rot_y_power2'
            'rot_z'
            'rot_z_derivative1'
            'rot_z_derivative1_power2'
            'rot_z_power2'
            'trans_x'
            'trans_x_derivative1'
            'trans_x_derivative1_power2'
            'trans_x_power2'
            'trans_y'
            'trans_y_derivative1'
            'trans_y_derivative1_power2'
            'trans_y_power2'
            'trans_z'
            'trans_z_derivative1'
            'trans_z_derivative1_power2'
            'trans_z_power2'
            'csf'
            'csf_derivative1'
            'csf_derivative1_power2'
            'csf_power2'
            'white_matter'
            'white_matter_derivative1'
            'white_matter_derivative1_power2'
            'white_matter_power2'};

        yaml_file = fullfile(xcpd_work_dir, sprintf('custom_config_%s_%s_run-%s.yaml', participant_label, task_id, run_id));
        dataset_fmriprep = 'custom';
        dataset_task = 'custom';

        fid = fopen(yaml_file, 'w');
        fprintf(fid, 'name: confounds_custom\n');
        fprintf(fid, 'description: |\n');
        fprintf(fid, '  Nuisance regression including fMRIPrep confounds and HRF-convolved task regressors.\n');
        fprintf(fid, 'confounds:\n');

        fprintf(fid, '  fmriprep_core:\n');
        fprintf(fid, '    dataset: %s\n', dataset_fmriprep);
        fprintf(fid, '    query:\n');
        fprintf(fid, '      space: null\n');
        fprintf(fid, '      cohort: null\n');
        fprintf(fid, '      res: null\n');
        fprintf(fid, '      den: null\n');
        fprintf(fid, '      desc: confounds\n');
        fprintf(fid, '      suffix: timeseries\n');
        fprintf(fid, '      extension: .tsv\n');
        fprintf(fid, '    columns:\n');
        for k = 1:numel(fmriprep_names)
            fprintf(fid, '      - %s\n', fmriprep_names{k});
        end

        fprintf(fid, '\n');
        fprintf(fid, '  task_hrf:\n');
        fprintf(fid, '    dataset: %s\n', dataset_task);
        fprintf(fid, '    query:\n');
        fprintf(fid, '      space: null\n');
        fprintf(fid, '      cohort: null\n');
        fprintf(fid, '      res: null\n');
        fprintf(fid, '      den: null\n');
        fprintf(fid, '      desc: confounds\n');
        fprintf(fid, '      suffix: timeseries\n');
        fprintf(fid, '      extension: .tsv\n');
        fprintf(fid, '    columns:\n');
        for k = 1:numel(task_names)
            fprintf(fid, '      - %s\n', task_names{k});
        end
        fclose(fid);

        % HERE WRITE .SH FILE THAT LAUNCHES CUSTOMIZED XCP-D PARAMETERS :)
        % DERIVED VARIABLES
        % [folder, name, ~] = fileparts(custom_confounds_matfile);
        % custom_confounds_tsvfile = fullfile(folder, [name '.tsv']);
        custom_dataset_root = custom_confounds_root_reg;
        custom_desc_file = fullfile(custom_confounds_dir, 'dataset_description.json');
        fid = fopen(custom_desc_file, 'w');
        fprintf(fid, '{\n');
        fprintf(fid, '  "Name": "Custom confounds dataset",\n');
        fprintf(fid, '  "BIDSVersion": "1.9.0",\n');
        fprintf(fid, '  "DatasetType": "derivatives",\n');
        fprintf(fid, '  "GeneratedBy": [\n');
        fprintf(fid, '    {\n');
        fprintf(fid, '      "Name": "MATLAB custom confounds workflow"\n');
        fprintf(fid, '    }\n');
        fprintf(fid, '  ]\n');
        fprintf(fid, '}\n');
        fclose(fid);

        work_dir = fullfile(gngrew_output_dir, sprintf('work_%s_%s_run-%s', participant_label, task_id, run_id));
        if ~exist(gngrew_output_dir, 'dir'); mkdir(gngrew_output_dir); end
        if ~exist(work_dir, 'dir'); mkdir(work_dir); end

        % CREATE BIDS FILTER JSON (NECESSARY TO MAKE XCP-D RUN FOR JUST 1 RUN !!!! ...)
        filter_struct = struct();
        filter_struct.bold = struct( ...
            'session', {{session_id}}, ...
            'task', {{task_id}}, ...
            'run', {{run_id}});
        filter_json_file = fullfile(gngrew_output_dir, sprintf('%s_%s_run-%s_filter.json', participant_label, task_id, run_id));
        fid = fopen(filter_json_file, 'w');
        fprintf(fid, '%s', jsonencode(filter_struct, 'PrettyPrint', true));
        fclose(fid);

        % WRITE LOCAL SHELL SCRIPT THAT RUNS XCP-D
        sh_file = fullfile(gngrew_output_dir, sprintf('run_xcpd_%s_%s_run-%s.sh', participant_label, task_id, run_id));
        script_lines = {
            '#!/usr/bin/env bash'
            'set -euo pipefail'
            ''
            'module purge'
            'module load apptainer'
            ''
            sprintf('XCP_DIR="%s"', XCP_DIR)
            'SIMG="${XCP_DIR}/XCP-D-0.10.7.simg"'
            sprintf('FMRIPREP_DIR="%s"', gngrew_data_dir)
            sprintf('OUTPUT_DIR="%s"', gngrew_output_dir)
            sprintf('ATLAS_DIR="%s"', ATLAS_DIR)
            sprintf('CUSTOM_DATASET="%s"', custom_dataset_root)
            sprintf('NUISANCE_YAML="%s"', yaml_file)
            sprintf('PART="%s"', participant_label)
            sprintf('TASK_ID="%s"', task_id)
            sprintf('SESSION_ID="%s"', session_id)
            sprintf('WORK_DIR="%s"', work_dir)
            sprintf('BIDS_FILTER="%s"', filter_json_file)
            ''
            'mkdir -p "${OUTPUT_DIR}" "${WORK_DIR}"'
            ''
            'NCPUS=2'
            'export OMP_NUM_THREADS="${NCPUS}"'
            'BINDPATHS="/users,/proj,/work"'
            ''
            'apptainer run --cleanenv -B "${BINDPATHS}" \'
            '  "${SIMG}" \'
            '  "${FMRIPREP_DIR}" \'
            '  "${OUTPUT_DIR}" \'
            '  participant \'
            '  --fs-license-file "/work/users/d/g/dga/tools/freesurfer/7.3.2/license.txt" \'
            '  --mode nichart \'
            '  --participant-label "${PART}" \'
            '  --session-id "${SESSION_ID}" \'
            '  --bids-filter-file "${BIDS_FILTER}" \'
            '  -d atlas="${ATLAS_DIR}" custom="${CUSTOM_DATASET}" \'
            '  --atlases Seitzman 4S456Parcels \'
            '  -p "${NUISANCE_YAML}" \'
            '  -t "${TASK_ID}" \'
            '  --file-format nifti \'
            '  --smoothing 0 \'
            '  --motion-filter-type lp \'
            '  --band-stop-min 6 \'
            '  --motion-filter-order 2 \'
            '  -f 0.1 \'
            '  --output-type interpolated \'
            '  --min-time 120 \'
            '  --lower-bpf 0.008 \'
            '  --upper-bpf 0.09 \'
            '  --min-coverage 0.5 \'
            '  --nprocs "${NCPUS}" \'
            '  --clean-workdir \'
            '  -w "${WORK_DIR}"'
            ''            };
        fid = fopen(sh_file, 'w');
        for k = 1:numel(script_lines)
            fprintf(fid, '%s\n', script_lines{k});
        end
        fclose(fid);
        system(sprintf('chmod +x "%s"', sh_file)); % make it executable in terminal :)

        % RUN THE SCRIPT
        command = sprintf('bash "%s"', sh_file);
        [status, cmdout] = system(command);
        disp(status)
        disp(cmdout)
        warning('-------------------------------------------------------------');
        warning('Postproc done with: %s, %s run-%s', sub_id, task_id, run_id);
    end

end
