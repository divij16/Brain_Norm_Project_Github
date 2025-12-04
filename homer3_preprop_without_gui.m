%% homer3_preproc_without_gui.m
%
% Homer3-style preprocessing in MATLAB (no GUI) and export
% continuous HbO/HbR SNIRF compatible with MNE.
%
% Order:
%   1) hmrR_PruneChannels        (on intensity)
%   2) hmrR_Intensity2OD
%   3) hmrR_MotionArtifactByChannel
%   4) hmrR_MotionCorrectSpline
%   5) hmrR_MotionCorrectWavelet
%   6) hmrR_BandpassFilt
%   7) hmrR_OD2Conc
%   8) Save .snirf file

clear; clc;

%% ---------------- 0) USER SETTINGS -------------------------

in_snirf = 'sub-bndy_ses-tp1_task-arith_run-1_nirs_with_events.snirf';
out_hb   = 'output_preproc_hb.snirf';

% PruneChannels
dRange      = [5e-3 9e-1];
SNRthresh   = 3;
SDRange     = [0 45];

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

%% ---------------- 1) LOAD SNIRF ---------------------------

fprintf('Loading SNIRF: %s\n', in_snirf);
snirf_in = SnirfLoad(in_snirf);

data  = snirf_in.data(1);   % DataClass (intensity)
probe = snirf_in.probe(1);
stim  = snirf_in.stim;

% Ensure wavelengths exist (edit if your device is different)
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

%% ---------------- 2) PRUNE CHANNELS -----------------------

fprintf('Running hmrR_PruneChannels...\n');

% Your version of Homer3 (from previous output) accepts mlActMan, tIncMan.
% If your hmrR_PruneChannels has an extra "reset" argument, set reset=0
% and add it at the end.
mlActAuto = hmrR_PruneChannels( ...
    data, probe, mlActMan, tIncMan, ...
    dRange, SNRthresh, SDRange);

% mlActAuto is a cell array in newer Homer3 versions.

%% ---------------- 3) INTENSITY → OD -----------------------

fprintf('Running hmrR_Intensity2OD...\n');
dod = hmrR_Intensity2OD(data);   % DataClass, optical density

%% -------- 4) MOTION ARTIFACT (BY CHANNEL) -----------------

fprintf('Running hmrR_MotionArtifactByChannel...\n');

% Signature from your tooltip:
% hmrR_MotionArtifactByChannel(data,probe,mlActMan,mlActAuto,...
%                              tIncMan,tMotion,tMask,std_thresh,amp_thresh)

[tIncAuto, tIncChAuto] = hmrR_MotionArtifactByChannel( ...
    dod, probe, mlActMan, mlActAuto, tIncMan, ...
    tMotion, tMask, STDEVthresh, AMPthresh);

%% -------- 5) SPLINE MOTION CORRECTION ---------------------

fprintf('Running hmrR_MotionCorrectSpline...\n');

% From tooltip: hmrR_MotionCorrectSpline(data_dod,mlAct,tIncCh,p,turnon)
turnon_spline = 1;

dod_spline = hmrR_MotionCorrectSpline( ...
    dod, mlActAuto, tIncChAuto, p_spline, turnon_spline);

%% -------- 6) WAVELET MOTION CORRECTION --------------------

fprintf('Running hmrR_MotionCorrectWavelet...\n');

turnon_wavelet = 1;  % 1 = apply wavelet, 0 = skip

% 5-argument signature:
% hmrR_MotionCorrectWavelet(data_dod, mlActMan, mlActAuto, iqr, turnon)
dod_wavelet = hmrR_MotionCorrectWavelet( ...
    dod_spline, mlActMan, mlActAuto, iqr_wavelet, turnon_wavelet);


%% ---------------- 7) BANDPASS FILTER ----------------------

fprintf('Running hmrR_BandpassFilt (hpf = %.3f, lpf = %.3f)...\n', hpf, lpf);
dod_filt = hmrR_BandpassFilt(dod_wavelet, hpf, lpf);

%% ---------------- 8) OD → HbO/HbR -------------------------

fprintf('Running hmrR_OD2Conc...\n');
dc = hmrR_OD2Conc(dod_filt, probe, ppf);   % DataClass (HbO/HbR)

%% 9) SAVE CONTINUOUS HbO/HbR SNIRF FOR MNE  ---------------------------

fprintf('Saving continuous HbO/HbR SNIRF:\n  %s\n', out_hb);

% Reuse the loaded SnirfClass object and just replace data & probe
snirf_out          = snirf_in;   % copy original (keeps meta/info)
snirf_out.data(1)  = dc;        % replace data block with preprocessed Hb
snirf_out.probe(1) = probe;     % ensure any probe edits (wavelengths) kept

SnirfSave(out_hb, snirf_out);

fprintf('Done.\n');
