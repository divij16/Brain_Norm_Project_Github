%% homer3_preproc_without_gui_glm.m
%
% Homer3-style preprocessing in MATLAB (no GUI) and export:
%   1) continuous HbO/HbR SNIRF (for MNE)
%   2) GLM results from hmrR_GLM saved as CSV

clear; clc;

%% ---------------- 0) USER SETTINGS -------------------------

in_snirf   = 'sub-bndy_ses-tp1_task-arith_run-1_nirs_with_events.snirf';
out_hb     = 'output_preproc_hb.snirf';         % preprocessed Hb SNIRF
out_glmcsv = 'glm_homer3_stats.csv';            % GLM stats CSV

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
stim  = snirf_in.stim;      % StimClass array

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

mlActAuto = hmrR_PruneChannels( ...
    data, probe, mlActMan, tIncMan, ...
    dRange, SNRthresh, SDRange);   % cell array

%% ---------------- 3) INTENSITY → OD -----------------------

fprintf('Running hmrR_Intensity2OD...\n');
dod = hmrR_Intensity2OD(data);   % DataClass, optical density

%% -------- 4) MOTION ARTIFACT (BY CHANNEL) -----------------

fprintf('Running hmrR_MotionArtifactByChannel...\n');

[tIncAuto, tIncChAuto] = hmrR_MotionArtifactByChannel( ...
    dod, probe, mlActMan, mlActAuto, tIncMan, ...
    tMotion, tMask, STDEVthresh, AMPthresh);
% tIncAuto is a cell array (one per block)

%% -------- 5) SPLINE MOTION CORRECTION ---------------------

fprintf('Running hmrR_MotionCorrectSpline...\n');

turnon_spline = 1;
dod_spline = hmrR_MotionCorrectSpline( ...
    dod, mlActAuto, tIncChAuto, p_spline, turnon_spline);

%% -------- 6) WAVELET MOTION CORRECTION --------------------

fprintf('Running hmrR_MotionCorrectWavelet...\n');

turnon_wavelet = 1;
dod_wavelet = hmrR_MotionCorrectWavelet( ...
    dod_spline, mlActMan, mlActAuto, iqr_wavelet, turnon_wavelet);

%% ---------------- 7) BANDPASS FILTER ----------------------

fprintf('Running hmrR_BandpassFilt (hpf = %.3f, lpf = %.3f)...\n', hpf, lpf);
dod_filt = hmrR_BandpassFilt(dod_wavelet, hpf, lpf);

%% ---------------- 8) OD → HbO/HbR -------------------------

fprintf('Running hmrR_OD2Conc...\n');
dc = hmrR_OD2Conc(dod_filt, probe, ppf);   % DataClass (HbO/HbR)

%% ---------------- 9) RUN GLM (hmrR_GLM) -------------------

fprintf('Running hmrR_GLM (GLM_HRF_Drift_SS_Concentration)...\n');

% GLM parameters matching GLM_HRF_Drift_SS_Concentration:
trange            = [-5.0 40.0];     % [tPre tPost] in seconds
glmSolveMethod    = 1;              % 1 = OLS
idxBasis          = 1;              % 1 = Gaussian basis
paramsBasis       = [1.0 1.0];      % [stdev step] (Gaussian width & spacing)
rhoSD_ssThresh    = 15.0;           % short-separation distance (mm)
flagNuisanceRMethod = 1;            % use SS with highest correlation
driftOrder        = 0;              % polynomial drift order
c_vector          = 0;              % no contrast (just condition betas)

% No auxiliary regressors or tCCA, so:
Aaux  = [];                         % no extra regressors
rcMap = [];                         % not used unless flagNuisanceRMethod==3

% data_y is the concentration DataClass (dc) with HbO/HbR
[data_yavg, data_yavgstd, nTrials, data_ynew, data_yresid, ...
 data_ysum2, beta_blks, yR_blks, hmrstats] = ...
    hmrR_GLM(dc, stim, probe, mlActAuto, Aaux, tIncAuto, rcMap, ...
             trange, glmSolveMethod, idxBasis, paramsBasis, ...
             rhoSD_ssThresh, flagNuisanceRMethod, driftOrder, c_vector);

% hmrstats contains:
%   hmrstats.beta_label : cell array of beta labels (Cond1, xDrift, SS, etc.)
%   hmrstats.tval       : [#Betas x #Channels x #Hb] t-stats
%   hmrstats.pval       : [#Betas x #Channels x #Hb] p-values
%   hmrstats.ml         : measurement list (source/det info)

%% 10) EXPORT GLM STATS TO CSV ------------------------------

fprintf('Exporting GLM stats to CSV: %s\n', out_glmcsv);

beta_label = hmrstats.beta_label;    % 1 x nBetas (cell)
tval       = hmrstats.tval;          % nBetas x nCh x nHb
pval       = hmrstats.pval;          % nBetas x nCh x nHb
ml         = hmrstats.ml;            % nCh x (at least 2: src, det)

[nBetas, nCh, nHb] = size(tval);

% Preallocate columns for table
nRows      = nBetas * nCh * nHb;
BetaIndex  = zeros(nRows,1);
BetaLabel  = cell(nRows,1);
HbLabel    = cell(nRows,1);
Src        = zeros(nRows,1);
Det        = zeros(nRows,1);
Tval       = zeros(nRows,1);
Pval       = zeros(nRows,1);

row = 0;
for iB = 1:nBetas
    thisLabel = beta_label{iB};  
    for iCh = 1:nCh
        src = ml(iCh,1);
        det = ml(iCh,2);
        for iHb = 1:nHb
            row = row + 1;
            BetaIndex(row) = iB;
            BetaLabel{row} = thisLabel;

            if iHb == 1
                HbLabel{row} = 'HbO';
            elseif iHb == 2
                HbLabel{row} = 'HbR';
            else
                HbLabel{row} = 'HbT';
            end

            Src(row)  = src;
            Det(row)  = det;
            Tval(row) = tval(iB,iCh,iHb);
            Pval(row) = pval(iB,iCh,iHb);
        end
    end
end

T = table(BetaIndex, BetaLabel, HbLabel, Src, Det, Tval, Pval);

writetable(T, out_glmcsv);

fprintf('GLM CSV written.\n');

%% 11) SAVE CONTINUOUS HbO/HbR SNIRF FOR MNE  ----------------

fprintf('Saving continuous HbO/HbR SNIRF:\n  %s\n', out_hb);

snirf_out          = snirf_in;   % copy original (keeps meta/info)
snirf_out.data(1)  = dc;        % replace data block with preprocessed Hb
snirf_out.probe(1) = probe;     % keep updated probe (wavelengths etc.)

SnirfSave(out_hb, snirf_out);

fprintf('Done.\n');
