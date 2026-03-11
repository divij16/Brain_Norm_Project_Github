clear; clc;

%% ---------------- 0) USER SETTINGS -------------------------
input_dir = '/Users/divijnalge/Library/CloudStorage/OneDrive-NanyangTechnologicalUniversity/ntu/Brain Norm/Brain_Norm_Files/preprocessed/tp1arith/sub-bndy';

out_dir = fullfile(input_dir, 'swapped_files');
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

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

    [~, base, ~] = fileparts(snirf_files(k).name);
    out_snirf = fullfile(out_dir, sprintf('%s_swapped.snirf', base));

    fprintf('\n=====================================================\n');
    fprintf('(%d/%d) Loading SNIRF: %s\n', k, numel(snirf_files), in_snirf);

    try
        %% ---------------- 3) LOAD SNIRF ---------------------------
        snirf_in = SnirfLoad(in_snirf);
        probe = snirf_in.probe(1);

        %% ---------------- 4) SWAP PROBE LABELS + POSITIONS --------------------
        probe = snirf_in.probe(1);
        
        src = probe.sourceLabels(:);
        det = probe.detectorLabels(:);
        
        if numel(det) ~= 16
            error("Expected 16 detector labels, got %d", numel(det));
        end
        if numel(src) ~= 20
            error("Expected 20 source labels, got %d", numel(src));
        end
        
        % --- Read positions if present ---
        srcPos2D = [];
        srcPos3D = [];
        detPos2D = [];
        detPos3D = [];
        
        if ~isempty(probe.sourcePos2D)
            srcPos2D = probe.sourcePos2D;
        end
        if ~isempty(probe.sourcePos3D)
            srcPos3D = probe.sourcePos3D;
        end
        if ~isempty(probe.detectorPos2D)
            detPos2D = probe.detectorPos2D;
        end
        if ~isempty(probe.detectorPos3D)
            detPos3D = probe.detectorPos3D;
        end
        
        % =========================
        % Swap detectors: D1<->D9, ..., D8<->D16
        % =========================
        det_new = det;
        
        if ~isempty(detPos2D)
            detPos2D_new = detPos2D;
        end
        if ~isempty(detPos3D)
            detPos3D_new = detPos3D;
        end
        
        for i = 1:8
            a = i;
            b = i + 8;
        
            % swap labels
            det_new([a b]) = det_new([b a]);
        
            % swap 2D positions
            if ~isempty(detPos2D)
                detPos2D_new([a b], :) = detPos2D_new([b a], :);
            end
        
            % swap 3D positions
            if ~isempty(detPos3D)
                detPos3D_new([a b], :) = detPos3D_new([b a], :);
            end
        end
        
        % =========================
        % Swap sources: S1<->S11, ..., S10<->S20
        % =========================
        src_new = src;
        
        if ~isempty(srcPos2D)
            srcPos2D_new = srcPos2D;
        end
        if ~isempty(srcPos3D)
            srcPos3D_new = srcPos3D;
        end
        
        for i = 1:10
            a = i;
            b = i + 10;
        
            % swap labels
            src_new([a b]) = src_new([b a]);
        
            % swap 2D positions
            if ~isempty(srcPos2D)
                srcPos2D_new([a b], :) = srcPos2D_new([b a], :);
            end
        
            % swap 3D positions
            if ~isempty(srcPos3D)
                srcPos3D_new([a b], :) = srcPos3D_new([b a], :);
            end
        end
        
        % --- Write back ---
        probe.detectorLabels = det_new;
        probe.sourceLabels   = src_new;
        
        if ~isempty(detPos2D)
            probe.detectorPos2D = detPos2D_new;
        end
        if ~isempty(detPos3D)
            probe.detectorPos3D = detPos3D_new;
        end
        if ~isempty(srcPos2D)
            probe.sourcePos2D = srcPos2D_new;
        end
        if ~isempty(srcPos3D)
            probe.sourcePos3D = srcPos3D_new;
        end
        
        % write probe back
        snirf_in.probe(1) = probe;
       
        %% ---------------- 5) TRUE CHANNEL SWAP ---------------------
        % Edit snirf_in.data(iBlk).measurementList entries directly
        nBlocks = length(snirf_in.data);

        for iBlk = 1:nBlocks
            dcBlk = snirf_in.data(iBlk);

            % Get the editable measurementList
            mlObj = dcBlk.measurementList;

            % print mlAct before swapping
            ml_BEFORE = dcBlk.GetMeasList();
            disp("First 5 rows of GetMeasList BEFORE swap:");
            disp(ml_BEFORE(1:5,:));
            disp("Last 5 rows of GetMeasList BEFORE swap:");
            disp(ml_BEFORE(length(ml_BEFORE) - 4:length(ml_BEFORE),:));

            % Case A: measurementList is a cell array
            if iscell(mlObj)
                %disp("MeasurementList is a cell array")
                for i = 1:numel(mlObj)
                    mlObj{i} = swap_one_meas_entry(mlObj{i});
                end

            % Case B: measurementList is a struct or object array
            elseif isstruct(mlObj) || isobject(mlObj)
                %disp("MeasurementList is a struct or object array")
                for i = 1:numel(mlObj)
                    mlObj(i) = swap_one_meas_entry(mlObj(i));
                end

            else
                % The toolbox stores measurementList in an unusual type
                error("measurementList is type %s; cannot safely edit. (iBlk=%d)", class(mlObj), iBlk);
            end

            % write back into the block, then back into snirf
            dcBlk.measurementList = mlObj;
            snirf_in.data(iBlk) = dcBlk;
        end

        dcTest = snirf_in.data(1);
        ml_AFTER = dcTest.GetMeasList();
        disp("First 5 rows of GetMeasList AFTER swap:");
        disp(ml_AFTER(1:5,:));
        disp("Last 5 rows of GetMeasList AFTER swap:");
        disp(ml_AFTER(length(ml_AFTER) - 4:length(ml_AFTER),:));

        %% ---------------- 6) SAVE OUTPUT SNIRF ---------------------
        fprintf('Saving: %s\n', out_snirf);
        SnirfSave(out_snirf, snirf_in);
        fprintf('Done: %s\n', snirf_files(k).name);

    catch ME
        fprintf(2, '\nERROR processing %s\n%s\n', snirf_files(k).name, ME.message);
        % Uncomment for full stack:
        % fprintf(2, '%s\n', getReport(ME, 'extended'));
    end
end

fprintf('\nAll runs finished.\n');

%% ---------------- LOCAL HELPER FUNCTIONS --------------------
function entry = swap_one_meas_entry(entry)
    % Works if entry is a struct OR an object with properties:
    % sourceIndex, detectorIndex (common in SNIRF/Homer3)

    % read
    s = entry.sourceIndex;
    d = entry.detectorIndex;

    % swap source index: 1–10 -> +10, 11–20 -> -10
    if s >= 1 && s <= 10
        s = s + 10;
    elseif s >= 11 && s <= 20
        s = s - 10;
    end

    % swap detector index: 1–8 -> +8, 9–16 -> -8
    if d >= 1 && d <= 8
        d = d + 8;
    elseif d >= 9 && d <= 16
        d = d - 8;
    end

    % write
    entry.sourceIndex   = s;
    entry.detectorIndex = d;
end