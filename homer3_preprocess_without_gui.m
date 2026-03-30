%% homer3_preproc_batch_without_gui_with_mlAct.m
% Batch Homer3-style preprocessing for all .snirf files in a folder
% and save continuous HbO/HbR SNIRF for MNE into:
%    input_dir/preproc_for_mne
%
% Also saves mlAct outputs into:
%    input_dir/preproc_for_mne/mlAct_masks

clear; clc;

%% ---------------- 0) USER SETTINGS -------------------------
input_dir = uigetdir;

if input_dir == 0
    disp('No folder selected');
else
    disp(['Selected folder: ', input_dir]);
end
%input_dir = '/Users/divijnalge/Library/CloudStorage/OneDrive-NanyangTechnologicalUniversity/ntu/Brain Norm/Brain_Norm_Files/preprocessed/tp1arith/sub-bnkg/swapped_files';

out_dir = fullfile(input_dir, 'preproc_for_mne');
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

mask_out_dir = fullfile(out_dir, 'mlAct_masks');
if ~exist(mask_out_dir, 'dir')
    mkdir(mask_out_dir);
end

% PruneChannels
dRange      = [5e-3 9e-1];
SNRthresh   = 3;
SDRange     = [0 45];

windowSec = 3;
hpf_prune = 0.5;
lpf_prune = 2.5;
ScanQualityThresh = 0.5;
SCIThresh = 0.7;
PSPThresh = 0.07;

% MotionArtifactByChannel
tMotion     = 0.5;
tMask       = 0.5;
STDEVthresh = 20;
AMPthresh   = 0.40;

% Spline & Wavelet
p_spline    = 0.99;
iqr_wavelet = 0.72;

% Bandpass
hpf         = 0.01;
lpf         = 0.20;

% OD2Conc
ppf         = [1 1];

% Optional: choose whether to process subfolders too
process_recursively = false;

%% ---------------- 1) FIND INPUT FILES ----------------------
if process_recursively
    snirf_files = dir(fullfile(input_dir, '**', '*.snirf'));
else
    snirf_files = dir(fullfile(input_dir, '*.snirf'));
end

if isempty(snirf_files)
    error('No .snirf files found in: %s', input_dir);
end

fprintf('Found %d SNIRF files.\n', numel(snirf_files));

%% ---------------- 2) INITIALIZE ----------------------------
nRuns = numel(snirf_files);

% mlAct outputs
mlAct_ALL = [];
mlAct_by_run = [];
run_names = strings(nRuns, 1);
valid_run = false(nRuns, 1);

kept_channels_by_run = cell(nRuns, 1);
dropped_channels_by_run = cell(nRuns, 1);
channel_names_ref = {};

%% ---------------- 3) LOOP OVER FILES -----------------------
for k = 1:nRuns

    in_snirf = fullfile(snirf_files(k).folder, snirf_files(k).name);
    run_names(k) = string(snirf_files(k).name);

    % Output filename: <original>_preproc_hb.snirf
    [~, base, ~] = fileparts(snirf_files(k).name);
    out_hb = fullfile(out_dir, sprintf('%s_preproc_hb.snirf', base));

    fprintf('\n=====================================================\n');
    fprintf('(%d/%d) Loading SNIRF: %s\n', k, nRuns, in_snirf);

    try
        %% ---------------- 3) LOAD SNIRF ---------------------------
        snirf_in = SnirfLoad(in_snirf);

        data  = snirf_in.data(1);   % intensity
        probe = snirf_in.probe(1);
        stim  = snirf_in.stim; %#ok<NASGU>  % kept via snirf_out copy

        % Ensure wavelengths exist
        if isempty(probe.wavelengths)
            probe.wavelengths = [760 850];
        end

        nBlk = 1;

        %% ---- initialise manual masks as CELL arrays (important) ---
        mlActMan = cell(nBlk,1);
        tIncMan  = cell(nBlk,1);
        for iBlk = 1:nBlk
            mlActMan{iBlk} = [];
            tIncMan{iBlk}  = [];
        end

        %% ---------------- 4) PRUNE CHANNELS -----------------------
        fprintf('Running hmrR_PruneChannels...\n');
        mlActAuto_basic = hmrR_PruneChannels( ...
             data, probe, mlActMan, tIncMan, ...
             dRange, SNRthresh, SDRange);

        fprintf('Running hmrR_PruneChannelsPlus_equalized...\n');
        mlActAuto = hmrR_PruneChannelsPlus_equalized( ...
            data, probe, mlActMan, tIncMan, windowSec, ...
            hpf_prune, lpf_prune, ScanQualityThresh, SCIThresh, PSPThresh, ...
            dRange, SNRthresh, SDRange);

        pruning_compare = isequal(mlActAuto{1}, mlActAuto_basic{1}); %#ok<NASGU>

        % ----------------------------------------------------------
        % Use equalized pruning output for mlAct exports
        % If you want the old/basic mask instead, change this to:
        %   mlAct_source = mlActAuto_basic{1};
        % ----------------------------------------------------------
        mlAct_source = mlActAuto{1};

        mask_r = logical(mlAct_source(:,3) == 1);   % Nx1 logical

        % Build channel names once from the first successful run
        if isempty(channel_names_ref)
            nCh = size(mlAct_source, 1);
            channel_names_ref = cell(nCh, 1);

            for ii = 1:nCh
                src_idx = mlAct_source(ii, 1);
                det_idx = mlAct_source(ii, 2);
                wl_idx  = mlAct_source(ii, 4);

                if ~isempty(probe.wavelengths) && wl_idx <= numel(probe.wavelengths)
                    wl_val = probe.wavelengths(wl_idx);
                    channel_names_ref{ii} = sprintf('S%d_D%d_%dnm', src_idx, det_idx, round(wl_val));
                else
                    channel_names_ref{ii} = sprintf('S%d_D%d_wl%d', src_idx, det_idx, wl_idx);
                end
            end

            mlAct_ALL = false(nCh,1);
            mlAct_by_run = false(nCh,nRuns);
        else
            if numel(mask_r) ~= size(mlAct_by_run,1)
                error('Run %d has %d channels, expected %d.', ...
                    k, numel(mask_r), size(mlAct_by_run,1));
            end
        end

        % Store kept / dropped names for this run
        kept_channels_by_run{k} = channel_names_ref(mask_r);
        dropped_channels_by_run{k} = channel_names_ref(~mask_r);

        % OR across runs
        mlAct_ALL = mlAct_ALL | mask_r;

        % Store this run
        mlAct_by_run(:,k) = mask_r;
        valid_run(k) = true;

        fprintf('Kept in this run: %d / %d\n', sum(mask_r), numel(mask_r));

        %% ---------------- 5) INTENSITY → OD -----------------------
        fprintf('Running hmrR_Intensity2OD...\n');
        dod = hmrR_Intensity2OD(data);

        %% -------- 6) MOTION ARTIFACT (BY CHANNEL) -----------------
        fprintf('Running hmrR_MotionArtifactByChannel...\n');
        [tIncAuto, tIncChAuto] = hmrR_MotionArtifactByChannel( ...
            dod, probe, mlActMan, mlActAuto, tIncMan, ...
            tMotion, tMask, STDEVthresh, AMPthresh); %#ok<NASGU>

        %% -------- 7) SPLINE MOTION CORRECTION ---------------------
        fprintf('Running hmrR_MotionCorrectSpline...\n');
        turnon_spline = 1;
        dod_spline = hmrR_MotionCorrectSpline( ...
            dod, mlActAuto, tIncChAuto, p_spline, turnon_spline);

        %% -------- 8) WAVELET MOTION CORRECTION --------------------
        fprintf('Running hmrR_MotionCorrectWavelet...\n');
        turnon_wavelet = 1;
        dod_wavelet = hmrR_MotionCorrectWavelet( ...
            dod_spline, mlActMan, mlActAuto, iqr_wavelet, turnon_wavelet);

        %% ---------------- 9) BANDPASS FILTER ----------------------
        fprintf('Running hmrR_BandpassFilt (hpf=%.3f, lpf=%.3f)...\n', hpf, lpf);
        dod_filt = hmrR_BandpassFilt(dod_wavelet, hpf, lpf);

        %% ---------------- 10) OD → HbO/HbR ------------------------
        fprintf('Running hmrR_OD2Conc...\n');
        dc = hmrR_OD2Conc(dod_filt, probe, ppf);

        %% ---------------- 11) SAVE OUTPUT SNIRF -------------------
        fprintf('Saving: %s\n', out_hb);

        snirf_out          = snirf_in;  % copy keeps meta + stim
        snirf_out.data(1)  = dc;        % replace with HbO/HbR
        snirf_out.probe(1) = probe;

        SnirfSave(out_hb, snirf_out);

        fprintf('Done: %s\n', snirf_files(k).name);

    catch ME
        fprintf(2, 'FAILED on: %s\nReason: %s\n', snirf_files(k).name, ME.message);
        continue;
    end
end

%% ---------------- 12) SAVE mlAct OUTPUTS --------------------
if isempty(mlAct_ALL)
    error('No valid runs were processed. No mlAct outputs were generated.');
end

mlAct_ALL = logical(mlAct_ALL);
mlAct_by_run = logical(mlAct_by_run);

kept_channels_all = channel_names_ref(mlAct_ALL);
dropped_channels_all = channel_names_ref(~mlAct_ALL);

% MATLAB file
mat_file = fullfile(mask_out_dir, 'mlAct_masks.mat');
save(mat_file, ...
    'mlAct_ALL', 'mlAct_by_run', 'run_names', 'valid_run', ...
    'channel_names_ref', ...
    'kept_channels_by_run', 'dropped_channels_by_run', ...
    'kept_channels_all', 'dropped_channels_all');

fprintf('\nSaved MATLAB file:\n  %s\n', mat_file);

% Python-ready text file
txt_file = fullfile(mask_out_dir, 'mlAct_masks_python.txt');
fid = fopen(txt_file, 'w');

% -------------------------------------------------
% Save OR mask
% -------------------------------------------------
fprintf(fid, 'mlAct_ALL = [');
for i = 1:length(mlAct_ALL)
    fprintf(fid, '%d', mlAct_ALL(i));
    if i < length(mlAct_ALL)
        fprintf(fid, ', ');
    end
end
fprintf(fid, ']\n\n');

% -------------------------------------------------
% Save per-run masks
% -------------------------------------------------
fprintf(fid, 'mlAct_by_run = [\n');
for k = 1:nRuns
    fprintf(fid, '  [');
    for i = 1:size(mlAct_by_run,1)
        fprintf(fid, '%d', mlAct_by_run(i,k));
        if i < size(mlAct_by_run,1)
            fprintf(fid, ', ');
        end
    end
    fprintf(fid, ']');
    if k < nRuns
        fprintf(fid, ',');
    end

    run_label = char(run_names(k));
    if valid_run(k)
        fprintf(fid, '   # %s\n', run_label);
    else
        fprintf(fid, '   # %s (FAILED)\n', run_label);
    end
end
fprintf(fid, ']\n\n');

% -------------------------------------------------
% Save reference channel names
% -------------------------------------------------
fprintf(fid, 'channel_names_ref = [\n');
for i = 1:numel(channel_names_ref)
    fprintf(fid, '    "%s"', channel_names_ref{i});
    if i < numel(channel_names_ref)
        fprintf(fid, ',');
    end
    fprintf(fid, '\n');
end
fprintf(fid, ']\n\n');

% -------------------------------------------------
% Save kept channels by run
% -------------------------------------------------
fprintf(fid, 'kept_channels_by_run = [\n');
for k = 1:nRuns
    fprintf(fid, '    [');
    kc = kept_channels_by_run{k};
    for i = 1:numel(kc)
        fprintf(fid, '"%s"', kc{i});
        if i < numel(kc)
            fprintf(fid, ', ');
        end
    end
    fprintf(fid, ']');
    if k < nRuns
        fprintf(fid, ',');
    end

    run_label = char(run_names(k));
    if valid_run(k)
        fprintf(fid, '   # %s\n', run_label);
    else
        fprintf(fid, '   # %s (FAILED)\n', run_label);
    end
end
fprintf(fid, ']\n\n');

% -------------------------------------------------
% Save dropped channels by run
% -------------------------------------------------
fprintf(fid, 'dropped_channels_by_run = [\n');
for k = 1:nRuns
    fprintf(fid, '    [');
    dc_run = dropped_channels_by_run{k};
    for i = 1:numel(dc_run)
        fprintf(fid, '"%s"', dc_run{i});
        if i < numel(dc_run)
            fprintf(fid, ', ');
        end
    end
    fprintf(fid, ']');
    if k < nRuns
        fprintf(fid, ',');
    end

    run_label = char(run_names(k));
    if valid_run(k)
        fprintf(fid, '   # %s\n', run_label);
    else
        fprintf(fid, '   # %s (FAILED)\n', run_label);
    end
end
fprintf(fid, ']\n\n');

% -------------------------------------------------
% Save OR-kept and OR-dropped channels across all runs
% -------------------------------------------------
fprintf(fid, 'kept_channels_all = [\n');
for i = 1:numel(kept_channels_all)
    fprintf(fid, '    "%s"', kept_channels_all{i});
    if i < numel(kept_channels_all)
        fprintf(fid, ',');
    end
    fprintf(fid, '\n');
end
fprintf(fid, ']\n\n');

fprintf(fid, 'dropped_channels_all = [\n');
for i = 1:numel(dropped_channels_all)
    fprintf(fid, '    "%s"', dropped_channels_all{i});
    if i < numel(dropped_channels_all)
        fprintf(fid, ',');
    end
    fprintf(fid, '\n');
end
fprintf(fid, ']\n');

fclose(fid);

fprintf('Saved Python-ready text file:\n  %s\n', txt_file);

% CSV files
csv_all = fullfile(mask_out_dir, 'mlAct_ALL.csv');
writematrix(double(mlAct_ALL), csv_all);

csv_runs = fullfile(mask_out_dir, 'mlAct_by_run.csv');
writematrix(double(mlAct_by_run), csv_runs);

fprintf('Saved CSV files:\n  %s\n  %s\n', csv_all, csv_runs);

%% ---------------- 13) SUMMARY ------------------------------
fprintf('\nAll done.\n');
fprintf('Preprocessed SNIRF outputs are in:\n  %s\n', out_dir);
fprintf('mlAct outputs are in:\n  %s\n', mask_out_dir);

fprintf('\nSummary:\n');
fprintf('mlAct_ALL keeps %d / %d channels\n', sum(mlAct_ALL), numel(mlAct_ALL));

for k = 1:nRuns
    if valid_run(k)
        fprintf('Run %d keeps %d / %d channels   (%s)\n', ...
            k, sum(mlAct_by_run(:,k)), size(mlAct_by_run,1), run_names(k));
    else
        fprintf('Run %d FAILED   (%s)\n', k, run_names(k));
    end
end