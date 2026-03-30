%% find_mlAct_masks_only.m
% Find:
%   1) mlAct_ALL   = OR across all runs
%   2) mlAct_by_run = per-run masks (92 x nRuns)
% and save them to files

clear; clc;

%% ---------------- 0) USER SETTINGS -------------------------
input_dir = uigetdir;

if input_dir == 0
    disp('No folder selected');
else
    disp(['Selected folder: ', input_dir]);
end
%input_dir = '/Users/divijnalge/Library/CloudStorage/OneDrive-NanyangTechnologicalUniversity/ntu/Brain Norm/Brain_Norm_Files/preprocessed/tp1arith/sub-bndy/swapped_files';

out_dir = fullfile(input_dir, 'mlAct_masks');
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

% PruneChannels settings
dRange      = [5e-3 9e-1];
SNRthresh   = 3;
SDRange     = [0 45];

windowSec = 3;
hpf_prune = 0.5;
lpf_prune = 2.5;
ScanQualityThresh = 0.5;
SCIThresh = 0.7;
PSPThresh = 0.07;

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

% Will initialize after first run once we know channel count
mlAct_ALL = [];
mlAct_by_run = [];
run_names = strings(nRuns, 1);

kept_channels_by_run = cell(nRuns, 1);
dropped_channels_by_run = cell(nRuns, 1);
channel_names_ref = {};

%% ---------------- 3) LOOP OVER FILES -----------------------
for k = 1:nRuns

    in_snirf = fullfile(snirf_files(k).folder, snirf_files(k).name);
    run_names(k) = string(snirf_files(k).name);

    fprintf('\n=====================================================\n');
    fprintf('(%d/%d) Loading SNIRF: %s\n', k, nRuns, in_snirf);

    try
        %% ---------------- LOAD SNIRF ---------------------------
        snirf_in = SnirfLoad(in_snirf);

        data  = snirf_in.data(1);
        probe = snirf_in.probe(1);

        % Ensure wavelengths exist
        if isempty(probe.wavelengths)
            probe.wavelengths = [760 850];
        end

        %% ---- initialize manual masks as CELL arrays ----------
        nBlk = 1;
        mlActMan = cell(nBlk,1);
        tIncMan  = cell(nBlk,1);
        for iBlk = 1:nBlk
            mlActMan{iBlk} = [];
            tIncMan{iBlk}  = [];
        end

        %% ---------------- PRUNE CHANNELS ----------------------
        fprintf('Running hmrR_PruneChannels_equalized...\n');

        % Choose ONE of the two below.
        % Option A: basic pruning
        %mlActAuto_2 = hmrR_PruneChannels( ...
        %    data, probe, mlActMan, tIncMan, ...
        %    dRange, SNRthresh, SDRange);

        % Option B: equalized pruning (uncomment if you want this instead)
        mlActAuto_2 = hmrR_PruneChannelsPlus_equalized( ...
             data, probe, mlActMan, tIncMan, windowSec, ...
             hpf_prune, lpf_prune, ScanQualityThresh, SCIThresh, PSPThresh, ...
             dRange, SNRthresh, SDRange);

        mask_r = logical(mlActAuto_2{1}(:,3) == 1);   % Nx1 logical

                % Build channel names once from the first run
        ml = mlActAuto_2{1};   % Nx4 matrix: [src det active wl]
        if isempty(channel_names_ref)
            nCh = size(ml, 1);
            channel_names_ref = cell(nCh, 1);

            for ii = 1:nCh
                src_idx = ml(ii, 1);
                det_idx = ml(ii, 2);
                wl_idx  = ml(ii, 4);

                % Use actual wavelength if available, else wl index
                if ~isempty(probe.wavelengths) && wl_idx <= numel(probe.wavelengths)
                    wl_val = probe.wavelengths(wl_idx);
                    channel_names_ref{ii} = sprintf('S%d_D%d_%dnm', src_idx, det_idx, round(wl_val));
                else
                    channel_names_ref{ii} = sprintf('S%d_D%d_wl%d', src_idx, det_idx, wl_idx);
                end
            end
        end

        % Store kept / dropped channel names for this run
        kept_channels_by_run{k} = channel_names_ref(mask_r);
        dropped_channels_by_run{k} = channel_names_ref(~mask_r);

        % Initialize sizes after first run
        if k == 1
            nCh = numel(mask_r);
            mlAct_ALL = false(nCh,1);
            mlAct_by_run = false(nCh,nRuns);
        else
            if numel(mask_r) ~= size(mlAct_by_run,1)
                error('Run %d has %d channels, expected %d.', ...
                    k, numel(mask_r), size(mlAct_by_run,1));
            end
        end

        % OR across runs
        mlAct_ALL = mlAct_ALL | mask_r;

        % Store this run
        mlAct_by_run(:,k) = mask_r;

        fprintf('Kept in this run: %d / %d\n', sum(mask_r), numel(mask_r));

    catch ME
        fprintf(2, 'FAILED on: %s\nReason: %s\n', snirf_files(k).name, ME.message);
        continue;
    end
end

%% ---------------- 4) SAVE FILES ----------------------------
mlAct_ALL = logical(mlAct_ALL);
mlAct_by_run = logical(mlAct_by_run);

mat_file = fullfile(out_dir, 'mlAct_masks.mat');
save(mat_file, 'mlAct_ALL', 'mlAct_by_run', 'run_names');

fprintf('\nSaved MATLAB file:\n  %s\n', mat_file);

%% ---------------- 5) SAVE PYTHON-READY TEXT (0/1 + channel names) ----------------
txt_file = fullfile(out_dir, 'mlAct_masks_python.txt');
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
    fprintf(fid, '   # %s\n', run_names(k));
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
    for i = 1:numel(kept_channels_by_run{k})
        fprintf(fid, '"%s"', kept_channels_by_run{k}{i});
        if i < numel(kept_channels_by_run{k})
            fprintf(fid, ', ');
        end
    end
    fprintf(fid, ']');
    if k < nRuns
        fprintf(fid, ',');
    end
    fprintf(fid, '   # %s\n', run_names(k));
end
fprintf(fid, ']\n\n');

% -------------------------------------------------
% Save dropped channels by run
% -------------------------------------------------
fprintf(fid, 'dropped_channels_by_run = [\n');
for k = 1:nRuns
    fprintf(fid, '    [');
    for i = 1:numel(dropped_channels_by_run{k})
        fprintf(fid, '"%s"', dropped_channels_by_run{k}{i});
        if i < numel(dropped_channels_by_run{k})
            fprintf(fid, ', ');
        end
    end
    fprintf(fid, ']');
    if k < nRuns
        fprintf(fid, ',');
    end
    fprintf(fid, '   # %s\n', run_names(k));
end
fprintf(fid, ']\n\n');

% -------------------------------------------------
% Save OR-kept and OR-dropped channels across all runs
% -------------------------------------------------
kept_channels_all = channel_names_ref(mlAct_ALL);
dropped_channels_all = channel_names_ref(~mlAct_ALL);

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

%% ---------------- 6) OPTIONAL CSV --------------------------
csv_all = fullfile(out_dir, 'mlAct_ALL.csv');
writematrix(double(mlAct_ALL), csv_all);

csv_runs = fullfile(out_dir, 'mlAct_by_run.csv');
writematrix(double(mlAct_by_run), csv_runs);

fprintf('Saved CSV files:\n  %s\n  %s\n', csv_all, csv_runs);

%% ---------------- 7) PRINT SUMMARY -------------------------
fprintf('\nSummary:\n');
fprintf('mlAct_ALL keeps %d / %d channels\n', sum(mlAct_ALL), numel(mlAct_ALL));

for k = 1:nRuns
    fprintf('Run %d keeps %d / %d channels   (%s)\n', ...
        k, sum(mlAct_by_run(:,k)), size(mlAct_by_run,1), run_names(k));
end