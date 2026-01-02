%% homer3_preproc_without_gui_glm_subjectlevel.m
%
% Runs:
%   A) Run-level GLM (hmrR_GLM) for EACH run -> CSV per run
%   B) Subject/session-level GLM (hmrS_GLM) per (ses, task) group -> CSV per group + HRF plots
%
% Needs Homer3 on MATLAB path. :contentReference[oaicite:1]{index=1}

clear; clc;

%% ---------------- 0) USER SETTINGS -------------------------

% Folder that contains your SNIRFs
run_dir = '/Users/divijnalge/Library/CloudStorage/OneDrive-NanyangTechnologicalUniversity/ntu/Brain Norm/Brain_Norm_Files/snirf_files_with_correct_stim_events_BNDY';

% Pattern to match your files
file_pattern = '*_nirs_with_events.snirf';

% Output folder (auto-created)
out_dir = fullfile(run_dir, 'homer3_glm_out');
if ~exist(out_dir, 'dir'); mkdir(out_dir); end

% ---------------- Preproc params (same as your script) ----------------
dRange      = [5e-3 9e-1];
SNRthresh   = 3;
SDRange     = [0 45];

tMotion     = 0.5;
tMask       = 0.5;
STDEVthresh = 20;
AMPthresh   = 0.40;

p_spline    = 0.99;
iqr_wavelet = 0.72;

hpf         = 0.01;
lpf         = 0.20;

ppf         = [1 1];

% ---------------- GLM params (you can auto-switch by task) ----------------
glmSolveMethod       = 1;         % OLS
idxBasis             = 1;         % Gaussian basis
paramsBasis          = [1.0 1.0]; % [step stdev] used by Homer3 Gaussian basis
rhoSD_ssThresh       = 15.0;      % short-sep threshold (mm)
flagNuisanceRMethod  = 1;         % short-sep with highest correlation
c_vector             = 0;         % no contrast

% Drift order:
% - If you KEEP hmrR_BandpassFilt: driftOrder = 0 (common choice)
% - If you REMOVE bandpass filtering: driftOrder = 3 (common choice)
driftOrder = 0;

% No aux regressors / no tCCA in your pipeline
Aaux  = [];
rcMap = [];

%% ---------------- 1) FIND FILES ----------------------------

F = dir(fullfile(run_dir, file_pattern));
if isempty(F)
    error('No SNIRF files found. Check run_dir and file_pattern.');
end

% Group by session+task using filename tokens "ses-..." and "task-..."
groups = containers.Map('KeyType','char','ValueType','any');

for k = 1:numel(F)
    fname = F(k).name;

    sesTok  = regexp(fname, 'ses-([^_]+)', 'tokens', 'once');
    taskTok = regexp(fname, 'task-([^_]+)', 'tokens', 'once');

    if isempty(sesTok);  ses = 'ses-NA';  else; ses  = ['ses-' sesTok{1}];  end
    if isempty(taskTok); task = 'task-NA'; else; task = ['task-' taskTok{1}]; end

    key = [ses '__' task];

    if ~isKey(groups, key)
        groups(key) = {};
    end

    groups(key) = [groups(key), {fullfile(F(k).folder, fname)}];
end

groupKeys = groups.keys;

%% ---------------- 2) LOOP GROUPS ---------------------------

for g = 1:numel(groupKeys)

    key = groupKeys{g};
    run_paths = groups(key);

    fprintf('\n====================================================\n');
    fprintf('GROUP: %s | runs = %d\n', key, numel(run_paths));
    fprintf('====================================================\n');

    % Choose trange based on task in filename (edit if your task names differ)
    if contains(key, 'task-visuowm')
        trange = [-2.0 12.5];
    else
        % default to arith style
        trange = [-5.0 40.0];
    end

    % Collect run-level outputs needed by hmrS_GLM
    nRuns = numel(run_paths);
    dcRuns       = cell(1, nRuns);
    stimRuns     = cell(1, nRuns);
    mlActRuns    = cell(1, nRuns);
    tIncAutoRuns = cell(1, nRuns);
    AauxRuns     = cell(1, nRuns);
    rcMapRuns    = cell(1, nRuns);

    % We'll reuse the probe from the first run (assumes same montage)
    probe_first = [];

    for r = 1:nRuns

        in_snirf = run_paths{r};
        [~, base, ~] = fileparts(in_snirf);

        fprintf('\n-- RUN %d/%d: %s\n', r, nRuns, base);

        %% LOAD
        snirf_in = SnirfLoad(in_snirf);

        data  = snirf_in.data(1);   % intensity
        probe = snirf_in.probe(1);
        stim  = snirf_in.stim;

        if isempty(probe.wavelengths)
            probe.wavelengths = [760 850];
        end

        if isempty(probe_first)
            probe_first = probe;
        end

        nBlk = 1;

        %% init manual masks (cells)
        mlActMan = cell(nBlk,1);
        tIncMan  = cell(nBlk,1);
        for iBlk = 1:nBlk
            mlActMan{iBlk} = [];
            tIncMan{iBlk}  = [];
        end

        %% PRUNE
        mlActAuto = hmrR_PruneChannels( ...
            data, probe, mlActMan, tIncMan, ...
            dRange, SNRthresh, SDRange);

        %% INTENSITY -> OD
        dod = hmrR_Intensity2OD(data);

        %% MOTION ARTIFACT DETECTION
        [tIncAuto, tIncChAuto] = hmrR_MotionArtifactByChannel( ...
            dod, probe, mlActMan, mlActAuto, tIncMan, ...
            tMotion, tMask, STDEVthresh, AMPthresh);

        %% SPLINE
        turnon_spline = 1;
        dod_spline = hmrR_MotionCorrectSpline( ...
            dod, mlActAuto, tIncChAuto, p_spline, turnon_spline);

        %% WAVELET
        turnon_wavelet = 1;
        dod_wavelet = hmrR_MotionCorrectWavelet( ...
            dod_spline, mlActMan, mlActAuto, iqr_wavelet, turnon_wavelet);

        %% BANDPASS
        dod_filt = hmrR_BandpassFilt(dod_wavelet, hpf, lpf);

        %% OD -> CONC (HbO/HbR)
        dc = hmrR_OD2Conc(dod_filt, probe, ppf);

        %% ---------------- RUN-LEVEL GLM (hmrR_GLM) ----------------
        fprintf('   Running run-level GLM (hmrR_GLM)...\n');

        [~, ~, ~, ~, ~, ~, ~, ~, hmrstats_run] = ...
            hmrR_GLM(dc, stim, probe, mlActAuto, Aaux, tIncAuto, rcMap, ...
                     trange, glmSolveMethod, idxBasis, paramsBasis, ...
                     rhoSD_ssThresh, flagNuisanceRMethod, driftOrder, c_vector);

        out_glm_run_csv = fullfile(out_dir, [base '_glm_run.csv']);
        export_hmrstats_csv(hmrstats_run, out_glm_run_csv);

        %% Save preprocessed Hb SNIRF (optional but useful for MNE)
        out_hb = fullfile(out_dir, [base '_hb_preproc.snirf']);
        snirf_out          = snirf_in;
        snirf_out.data(1)  = dc;
        snirf_out.probe(1) = probe;

        % NOTE: If your Homer3 expects reversed args, swap these two:
        SnirfSave(out_hb, snirf_out);

        %% Store for SUBJECT/SESSION GLM
        dcRuns{r}       = dc;
        stimRuns{r}     = stim;
        mlActRuns{r}    = mlActAuto;
        tIncAutoRuns{r} = tIncAuto;

        % hmrS_GLM concatenation code uses cell2mat(AauxRuns{r}), so make it a cell
        AauxRuns{r} = {Aaux};
        rcMapRuns{r} = rcMap;
    end

    %% ---------------- SUBJECT/SESSION GLM (hmrS_GLM) ----------------

    % Example for run r (adjust indexing to your variables)
    ml = dcRuns{r}.GetMeasListSrcDetPairs('reshape');   % rows might be 92
    mlAct_pair = cell2mat(mlActRuns{r});               % currently 46

    [pairList, ~, pairIdx] = unique(ml(:,1:2), 'rows', 'stable');
    mlAct_full = mlAct_pair(pairIdx);                 % now 92

    mlActRuns{r} = {mlAct_full};  % keep Homer3's "cell-of-cell" style

    fprintf('\n*** Running subject/session-level GLM (hmrS_GLM) for group: %s ***\n', key);

    [dcAvg, dcAvgStd, nTrialsS, dcNew, dcResid, dcSum2, betaS, RS, hmrstats_sess] = ...
        hmrS_GLM(dcRuns, stimRuns, probe_first, mlActRuns, AauxRuns, tIncAutoRuns, rcMapRuns, ...
                 trange, glmSolveMethod, idxBasis, paramsBasis, ...
                 rhoSD_ssThresh, flagNuisanceRMethod, driftOrder, c_vector);

    out_glm_sess_csv = fullfile(out_dir, [key '_glm_session.csv']);
    export_hmrstats_csv(hmrstats_sess, out_glm_sess_csv);

    %% ---------------- PLOTS (quick sanity check) ----------------
    % Plot mean HRF across channels for HbO/HbR, each condition
    plot_session_mean_hrf(dcAvg, key, out_dir);

    fprintf('DONE group: %s\n', key);
end

fprintf('\nAll groups complete. Outputs in:\n  %s\n', out_dir);

%% ===================== HELPERS ============================

function export_hmrstats_csv(hmrstats, out_csv)
    beta_label = hmrstats.beta_label;
    tval       = hmrstats.tval;
    pval       = hmrstats.pval;
    ml         = hmrstats.ml;

    [nBetas, nCh, nHb] = size(tval);

    % --- FIX: make beta_label length match nBetas (pad/trim) ---
    if isempty(beta_label)
        beta_label = cell(nBetas,1);
    end

    if numel(beta_label) < nBetas
        warning('beta_label shorter than tval rows (%d < %d). Padding labels.', numel(beta_label), nBetas);
        for k = (numel(beta_label)+1):nBetas
            beta_label{k} = sprintf('Regressor_%03d', k);
        end
    elseif numel(beta_label) > nBetas
        warning('beta_label longer than tval rows (%d > %d). Trimming labels.', numel(beta_label), nBetas);
        beta_label = beta_label(1:nBetas);
    end
    % --- end FIX ---


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
    writetable(T, out_csv);
    fprintf('   Wrote: %s\n', out_csv);

    fprintf('DEBUG export: nBetas=%d, length(beta_label)=%d\n', nBetas, numel(beta_label));

end

function plot_session_mean_hrf(data_yavg, tag, out_dir)
    if isempty(data_yavg)
        warning('No data_yavg to plot.');
        return;
    end

    t = data_yavg(1).GetTime();
    Y = data_yavg(1).GetDataTimeSeries('reshape');

    % Expecting something like: time x Hb x Ch x Cond
    if ndims(Y) < 4
        warning('Unexpected HRF data shape; skipping plot.');
        return;
    end

    nHb   = size(Y,2);
    nCond = size(Y,4);

    HbNames = {'HbO','HbR','HbT'};

    for hb = 1:min(nHb,2) % plot HbO/HbR
        figure('Name',[tag ' ' HbNames{hb}]); hold on;

        for c = 1:nCond
            tmp = squeeze(Y(:,hb,:,c)); % time x ch

            % NaN-safe mean across channels (no toolbox needed)
            denom = sum(~isnan(tmp), 2);
            numer = sum(tmp .* (~isnan(tmp)), 2);
            m = numer ./ max(denom, 1);

            plot(t, m);
        end

        xline(0);
        xlabel('Time (s)');
        ylabel('Concentration (a.u.)');
        legend(compose('Cond %d',1:nCond),'Location','best');
        title([tag ' | mean HRF across channels | ' HbNames{hb}]);

        out_png = fullfile(out_dir, [tag '_' HbNames{hb} '_meanHRF.png']);
        saveas(gcf, out_png);
        close(gcf);

        fprintf('   Saved plot: %s\n', out_png);
    end
end
