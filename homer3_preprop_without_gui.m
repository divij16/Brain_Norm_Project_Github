%% homer3_preproc_batch_without_gui.m
% Batch Homer3-style preprocessing for all .snirf files in a folder
% and save continuous HbO/HbR SNIRF for MNE into:
%    input_dir/preproc_for_mne

clear; clc;

%% ---------------- 0) USER SETTINGS -------------------------
input_dir = '/Users/divijnalge/Library/CloudStorage/OneDrive-NanyangTechnologicalUniversity/ntu/Brain Norm/Brain_Norm_Files/snirf_files_with_correct_stim_events_BNDY';

out_dir = fullfile(input_dir, 'preproc_for_mne');
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
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

%% ---------------- 2) LOOP OVER FILES -----------------------
for k = 1:numel(snirf_files)

    in_snirf = fullfile(snirf_files(k).folder, snirf_files(k).name);

    % Output filename: <original>_preproc_hb.snirf
    [~, base, ~] = fileparts(snirf_files(k).name);
    out_hb = fullfile(out_dir, sprintf('%s_preproc_hb.snirf', base));

    fprintf('\n=====================================================\n');
    fprintf('(%d/%d) Loading SNIRF: %s\n', k, numel(snirf_files), in_snirf);

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

        mlActAuto_2 = hmrR_PruneChannels( ...
             data, probe, mlActMan, tIncMan, ...
             dRange, SNRthresh, SDRange);

        mlActAuto = hmrR_PruneChannelsPlus_equalized( ...
            data, probe, mlActMan, tIncMan, windowSec, ...
            hpf_prune, lpf_prune, ScanQualityThresh, SCIThresh, PSPThresh, ...
            dRange, SNRthresh, SDRange);

        pruning_compare = mlActAuto{1} == mlActAuto_2{1}; %#ok<NASGU>

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
        % Continue to next file instead of stopping the whole batch
        continue;
    end
end

fprintf('\nAll done. Outputs are in:\n  %s\n', out_dir);
