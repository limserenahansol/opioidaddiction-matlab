%% ================================================================
%  Longitudinal aggregator for addiction_concate  (run-numbered)
%  Scans: BASE\day1..dayN\<cage>\<mouse>\concat_out_*\combined_pupil_digital.(csv|xlsx)
%  Also reads *.jsonl to attach:
%     - Session_Paradigm (e.g., "Progressive Ratio", "Fixed Ratio", "MicroLogReplay Passive", "Unknown")
%     - isPassive (boolean; 1 if passive by label OR if no req/requirement found anywhere)
%     - RequirementLast (last requirement anywhere in the JSONL)
%     - Trial (per-frame trial index; prefer Injector_TTL; fallback: reward_cmd; fallback2: trial_start)
%     - TrialRequirement (per-frame requirement from jsonl)
%  Adds Immersion latency (Excel/CSV) and TST/HOT summaries (frames & %).
%  Vertcat-safe: harmonizes columns/types/order across sessions.
%  Outputs in: BASE\longitudinal_outputs\run_###
%  NEW: coverage_per_mouse.csv, issues.txt, expanded summary.txt
% ================================================================
clear; clc;

% -------------------- EDIT THIS --------------------
BASE = 'K:\addiction_concate_Dec_2025';
% --------------------------------------------------

% ----- make output root with auto-increment number run_### -----
OUT_BASE = fullfile(BASE, 'longitudinal_outputs');
if ~exist(OUT_BASE,'dir'), mkdir(OUT_BASE); end
OUT_ROOT = nextRunRoot(OUT_BASE);           % e.g. ...\run_003
OUT_PM   = fullfile(OUT_ROOT, 'per_mouse');
if ~exist(OUT_PM,'dir'), mkdir(OUT_PM); end

ISSUES = {};  % collect human-readable issues/warnings

% ---------- collect strictly named day folders day1, day2, ... ----------
dDays = dir(fullfile(BASE, 'day*'));
isGoodDay = arrayfun(@(d) d.isdir && ~ismember(d.name,{'.','..'}) && ...
    ~isempty(regexp(lower(d.name), '^day\d+$','once')), dDays);
if any(~isGoodDay)
    bads = {dDays(~isGoodDay).name};
    for i=1:numel(bads)
        ISSUES{end+1} = sprintf('[name] Nonconforming day folder ignored: %s', bads{i});
    end
end
dDays = dDays(isGoodDay);
getDayIdx = @(nm) sscanf(nm, 'day%d');
[~,ord]   = sort(arrayfun(@(d)getDayIdx(d.name), dDays));
dDays     = dDays(ord);
if isempty(dDays), error('No "dayN" folders found under %s', BASE); end

ALL      = {};
metaRows = {};        % per-session inventory rows
dayRows  = {};        % per-day coverage rows (for per-mouse coverage summary)

for iday = 1:numel(dDays)
    dayName  = dDays(iday).name;
    dayPath  = fullfile(dDays(iday).folder, dayName);
    dayIndex = getDayIdx(dayName);

    % ---- EVERY subfolder is a cage folder ----
    dCages = dir(dayPath);
    isCage = arrayfun(@(d) d.isdir && ~ismember(d.name,{'.','..'}), dCages);
    dCages = dCages(isCage);

    for ic = 1:numel(dCages)
        cageFolderName = dCages(ic).name;
        cagePath       = fullfile(dCages(ic).folder, cageFolderName);
        cageID4        = canonicalCageID(cageFolderName);  % first 4 digits, else first 4 chars

        % ---- any subfolder is a mouse ----
        dMice = dir(cagePath);
        isMouse = arrayfun(@(d) d.isdir && ~ismember(d.name,{'.','..'}), dMice);
        dMice = dMice(isMouse);

        for im = 1:numel(dMice)
            mouseID   = dMice(im).name;
            mouseDir  = fullfile(dMice(im).folder, mouseID);
            mouse_key = string([cageID4 '_' mouseID]);

            % coverage flags for this (mouse, day)
            hasCSVDay  = false;
            hasJSONDay = false;
            hasTSTDay  = false;
            hasHOTDay  = false;
            hasIMMDay  = false;
            hasSTRAUBDay = false;

            % sessions inside mouse: any folder starting with concat_out
            dSess = dir(fullfile(mouseDir, 'concat_out*'));
            isSess = arrayfun(@(d) d.isdir && ~ismember(d.name,{'.','..'}), dSess);
            dSess = dSess(isSess);
            if isempty(dSess)
                ISSUES{end+1} = sprintf('[session] No concat_out* under %s (day=%s, cage=%s, mouse=%s)', ...
                                        mouseDir, dayName, cageFolderName, mouseID);
            end

            % sort sessions (timestamp inside name preferred)
            sessKey = zeros(numel(dSess),1);
            for s=1:numel(dSess)
                tok = regexp(dSess(s).name, '\d{8}_\d{6}','match','once');
                if ~isempty(tok)
                    try, sessKey(s) = datenum(tok,'yyyymmdd_HHMMSS'); catch, sessKey(s)=dSess(s).datenum; end
                else
                    sessKey(s) = dSess(s).datenum;
                end
            end
            [~,so]=sort(sessKey); dSess=dSess(so);

            % day-level Imm/TST/HOT presence (same per session)
            [immLatency_day, immFile_day] = findImmersionLatency(mouseDir);
            hasIMMDay = hasIMMDay || ~isempty(immFile_day);

            [tstCols_day, tstFile_day] = summarizeManualScoringMat(mouseDir, 'ManualScoringResults*tst*.mat', 'TST');
            hasTSTDay = hasTSTDay || ~isempty(tstFile_day);

            [hotCols_day, hotFile_day] = summarizeManualScoringMat(mouseDir, 'ManualScoringResults*hot*.mat', 'HOT');
            hasHOTDay = hasHOTDay || ~isempty(hotFile_day);
            [straubCols_day, straubFile_day] = summarizeManualScoringMat(mouseDir, 'ManualScoringResults*straub*.mat', 'STRAUB');
hasSTRAUBDay = hasSTRAUBDay || ~isempty(straubFile_day);


            nSessKept = 0;

            for isess = 1:numel(dSess)
                sessPath = fullfile(dSess(isess).folder, dSess(isess).name);
                [csvPath, usedName] = findCombinedFile(sessPath);
                if isempty(csvPath)
                    ISSUES{end+1} = sprintf('[csv] Missing combined pupil file under %s', sessPath);
                    continue;
                end
                hasCSVDay = true;

                % Read combined table + normalize columns
                T = safeReadCombined(csvPath);
                if ~hasCombinedColumns(T)
                    ISSUES{end+1} = sprintf('[csv] Combined file at %s is missing expected columns; filled placeholders', csvPath);
                end

                % === JSONL: parse trials, last requirement, reward & trial_start times, paradigm, passive ===
                [reqLast, reqByTrial, rewardElapsed, reqByReward, trialStartElapsed, paradigm, isPassive, jsonPath, parseNotes] = ...
                    findAndParseJSONL(dayPath, cagePath, mouseDir, sessPath, cageID4, mouseID);

                hasJSON = ~isempty(jsonPath);
                hasJSONDay = hasJSONDay || hasJSON;

                if hasJSON && ~isempty(parseNotes)
                    ISSUES{end+1} = sprintf('[jsonl] %s (file=%s)', parseNotes, jsonPath);
                elseif ~hasJSON
                    ISSUES{end+1} = sprintf('[jsonl] No JSONL found near session: %s', sessPath);
                end

                % === Trials: prefer Injector TTL; else reward_cmd; else trial_start ===
                [trialIdx, trialReq, alignMethod] = assignTrialsSmart(T, reqByTrial, rewardElapsed, reqByReward, trialStartElapsed);

                if strcmpi(alignMethod,'none')
                    ISSUES{end+1} = sprintf('[trials] No TTL/reward/trial_start to align trials (session=%s)', sessPath);
                end

                % per-frame adds
                T.Trial            = double(trialIdx);
                T.TrialRequirement = double(trialReq);
                T.TrialRequirement(T.Trial==0) = 0;
                % --- NOTICE + mark TTL-only trials as "extra" instead of NaN ---
% Compare number of trials defined by injector TTL vs JSONL
nTTLtrials = max(T.Trial);
if isempty(nTTLtrials) || ~isfinite(nTTLtrials), nTTLtrials = 0; end

% Prefer reward_cmd count; fallback to reqByTrial; then trial_start
if ~isempty(rewardElapsed)
    nJsonlTrials = numel(rewardElapsed); baseLabel = 'reward_cmd';
elseif ~isempty(reqByTrial)
    nJsonlTrials = numel(reqByTrial);    baseLabel = 'reqByTrial';
elseif ~isempty(trialStartElapsed)
    nJsonlTrials = numel(trialStartElapsed); baseLabel = 'trial_start';
else
    nJsonlTrials = 0; baseLabel = 'jsonl_none';
end

% If injector TTL yields more trials than JSONL, label those as "extra"
if nTTLtrials > nJsonlTrials
    idxExtra = (T.Trial > nJsonlTrials) & (T.Trial > 0) & isnan(T.TrialRequirement);
    if any(idxExtra)
        % convert this column to string for this session and fill "extra"
        T.TrialRequirement = string(T.TrialRequirement);
        T.TrialRequirement(idxExtra) = "extra";
        % keep '0' rows as "0"
        T.TrialRequirement(T.Trial==0) = "0";
        ISSUES{end+1} = sprintf(['[notice] %s | %s | sess=%d: TTL/injector trials=%d ' ...
                                 'vs JSONL(%s)=%d → marked %d rows TrialRequirement="extra".'], ...
                                 dayName, mouseID, isess, nTTLtrials, baseLabel, nJsonlTrials, nnz(idxExtra));
    end
elseif nTTLtrials < nJsonlTrials
    % Just a notice if JSONL has more trials than TTL
    ISSUES{end+1} = sprintf('[notice] %s | %s | sess=%d: TTL/injector trials=%d < JSONL(%s)=%d.', ...
                            dayName, mouseID, isess, nTTLtrials, baseLabel, nJsonlTrials);
end
% --- end NOTICE block ---

                T.RequirementLast  = repmat(double(reqLast), height(T), 1);
                T.Session_Paradigm = repmat(string(paradigm), height(T), 1);
                T.isPassive        = repmat(double(isPassive), height(T), 1);
                nTrialsJsonl = max([ numel(reqByTrial), numel(rewardElapsed), numel(trialStartElapsed) ]);  % record only; no dropping

          
                % ===== Immersion latency (Excel/CSV) — optional =====
                T.Immersion_Latency_s = repmat(double(immLatency_day), height(T), 1);   % NaN if none
                T.Immersion_File      = repmat(string(immFile_day), height(T), 1);      % "" if none

                % ===== TST summaries (frames + percent per behavior) — optional =====
                [T, nTST] = addScalarCols(T, tstCols_day);
                T.TST_File = repmat(string(tstFile_day), height(T), 1);

                % ===== HOT summaries (frames + percent per behavior) — optional =====
                [T, nHOT] = addScalarCols(T, hotCols_day);
                T.HOT_File = repmat(string(hotFile_day), height(T), 1);
                % ===== STRAUB summaries (frames + percent per behavior) — optional =====
[T, nSTRAUB] = addScalarCols(T, straubCols_day);
T.STRAUB_File = repmat(string(straubFile_day), height(T), 1);


                % ---- metadata ----
                T.day_index     = repmat(double(dayIndex), height(T),1);
                T.day_name      = repmat(string(dayName), height(T),1);
                T.cage_folder   = repmat(string(cageFolderName), height(T),1);
                T.cage          = repmat(string(cageID4), height(T),1);
                T.mouse_id      = repmat(string(mouseID), height(T),1);
                T.session_idx   = repmat(double(isess), height(T),1);
                T.session_dir   = repmat(string(sessPath), height(T),1);
                T.source_file   = repmat(string(usedName), height(T),1);
                T.mouse_key     = repmat(mouse_key, height(T),1);

metaCols = {'mouse_key','cage','cage_folder','mouse_id','day_index','day_name', ...
            'session_idx','session_dir','source_file','Session_Paradigm','isPassive', ...
            'RequirementLast','Immersion_Latency_s','Immersion_File','TST_File','HOT_File','STRAUB_File'};

                T = putFirst(T, metaCols);

                ALL{end+1} = T; %#ok<AGROW>
                nSessKept = nSessKept + 1;

                % Per-session inventory (compact)
                rewardStr     = join(string(round(rewardElapsed,6)), ',');
                trialStartStr = join(string(round(trialStartElapsed,6)), ',');
metaRows(end+1,:) = { ...
    mouse_key, dayIndex,dayName,cageID4,cageFolderName,mouseID,isess,height(T), ...
    sessPath, usedName, jsonPath, paradigm, double(isPassive), ...
    reqLast, numel(reqByTrial), char(rewardStr), char(trialStartStr), ...
    double(immLatency_day),immFile_day, nTST, nHOT, nSTRAUB, tstFile_day, hotFile_day, straubFile_day, ...
    hasJSON, true, ~isempty(tstFile_day), ~isempty(hotFile_day), ~isempty(straubFile_day), ~isempty(immFile_day), alignMethod};
 %#ok<AGROW>
            end

            % If no session parsed any JSONL but a JSONL exists up the tree, mark presence (for day coverage)
            if ~hasJSONDay
                jsonPathAny = findJsonlTopOnly(dayPath, cagePath, mouseDir, cageID4, mouseID);
                hasJSONDay = ~isempty(jsonPathAny);
                if hasJSONDay
                    ISSUES{end+1} = sprintf('[jsonl] JSONL present for day but not tied to any session (day=%s, cage=%s, mouse=%s): %s', ...
                                             dayName, cageFolderName, mouseID, jsonPathAny);
                end
            end

            % day coverage row
dayRows(end+1,:) = { ...
    mouse_key, cageID4, cageFolderName, mouseID, dayIndex, dayName, ...
    hasCSVDay, hasJSONDay, hasTSTDay, hasHOTDay, hasSTRAUBDay, hasIMMDay, nSessKept};

        end
    end
end

if isempty(ALL), error('No data collected under %s', BASE); end

% ---------- Harmonize (same vars, order, types) then concatenate ----------
ALL = makeTablesVertcatSafe(ALL);
Grand = vertcat(ALL{:});

% Sort nicely
keyCols = {'mouse_key','day_index','session_idx'};
if ismember('Frame', Grand.Properties.VariableNames)
    Grand = sortrows(Grand, [keyCols, {'Trial','Frame'}]);
else
    Grand = sortrows(Grand, [keyCols, {'Trial'}]);
end

% ---- save outputs ----
if ~exist(OUT_ROOT,'dir'), mkdir(OUT_ROOT); end
if ~exist(OUT_PM,'dir'),   mkdir(OUT_PM);   end

grandCSV = fullfile(OUT_ROOT,'ALL_mice_longitudinal.csv');
writetable(Grand, grandCSV);
fprintf('Wrote grand longitudinal CSV:\n  %s\n', grandCSV);

% per-mouse
mKeys = unique(Grand.mouse_key,'stable');
for i=1:numel(mKeys)
    mk  = string(mKeys(i));
    sub = Grand(Grand.mouse_key==mk,:);
    parts = split(mk,"_",2);
    cg = char(parts(1)); ms = char(parts(2));
    outDir = fullfile(OUT_PM, cg, ms); if ~exist(outDir,'dir'), mkdir(outDir); end
    outFile = fullfile(outDir, sprintf('mouse_%s_ALLDAYS.csv', regexprep(char(mk), '[:\\/*?"<>|]', '-')));
    writetable(sub, outFile);
end
fprintf('Wrote %d per-mouse CSV(s) under:\n  %s\n', numel(mKeys), OUT_PM);

% ---------- build per-session inventory table ----------
S = cell2table(metaRows,'VariableNames',{ ...
    'mouse_key','day_index','day_name','cage','cage_folder','mouse_id','session_idx','rows','session_dir', ...
    'source_file','jsonl_file','Session_Paradigm','isPassive','RequirementLast','nTrialsJsonl', ...
    'RewardTimes_s','TrialStartTimes_s','Immersion_Latency_s','Immersion_File', ...
    'nTSTcols','nHOTcols','nSTRAUBcols','TST_File','HOT_File','STRAUB_File', ...
    'hasJSONL','hasCSV','hasTST','hasHOT','hasSTRAUB','hasIMM','alignMethod'});

invCSV = fullfile(OUT_ROOT,'session_inventory.csv'); writetable(S, invCSV);

% ---------- per-mouse per-day coverage summary ----------
D = cell2table(dayRows,'VariableNames', { ...
    'mouse_key','cage','cage_folder','mouse_id','day_index','day_name', ...
    'hasCSVDay','hasJSONDay','hasTSTDay','hasHOTDay','hasSTRAUBDay','hasIMMDay','nSessions'});

D = sortrows(D, {'day_index','cage','mouse_id'});

% Aggregate per-mouse counts
C = groupsummary(D, 'mouse_key', 'sum', {'hasCSVDay','hasJSONDay','hasTSTDay','hasHOTDay','hasSTRAUBDay','hasIMMDay'});
C.Properties.VariableNames = {'mouse_key','GroupCount','Days_with_CSV','Days_with_JSONL','Days_with_TST','Days_with_HOT','Days_with_STRAUB','Days_with_Immersion'};
U = groupsummary(D, 'mouse_key', @(x) numel(unique(x)), 'day_index');
C.Total_Days = U.fun1_day_index;
C.GroupCount = [];  % drop
covCSV = fullfile(OUT_ROOT,'coverage_per_mouse.csv'); writetable(C, covCSV);

% ---------- write summary ----------
summaryPath = fullfile(OUT_ROOT,'summary.txt');
fid = fopen(summaryPath,'w');
fprintf(fid,'Longitudinal aggregation summary\n');
fprintf(fid,'Base: %s\n', BASE);
fprintf(fid,'Output root: %s\n', OUT_ROOT);
fprintf(fid,'Grand rows: %d\n', height(Grand));
fprintf(fid,'Mice: %d\n\n', numel(mKeys));
fprintf(fid,'Session inventory CSV: %s\n', invCSV);
fprintf(fid,'Per-mouse coverage CSV: %s\n\n', covCSV);

% Sessions block
fprintf(fid,'--- Session inventory (compact) ---\n');
for r=1:height(S)
    fprintf(fid,['day=%d (%s) | cage=%s (%s) | mouse=%s | session=%d | rows=%d | ' ...
                 'paradigm=%s | passive=%d | lastReq=%s | jsonlTrials=%d | align=%s | Immersion=%.3f s | TSTcols=%d | HOTcols=%d\n'], ...
        S.day_index(r), S.day_name{r}, S.cage{r}, S.cage_folder{r}, S.mouse_id{r}, S.session_idx(r), S.rows(r), ...
        S.Session_Paradigm{r}, S.isPassive(r), num2str(S.RequirementLast(r)), S.nTrialsJsonl(r), S.alignMethod{r}, ...
        S.Immersion_Latency_s(r), S.nTSTcols(r), S.nHOTcols(r));
end
fprintf(fid,'\n');

% Coverage block
fprintf(fid,'--- Per-mouse day coverage ---\n');
for i=1:height(C)
    fprintf(fid,'mouse=%s | days=%d | CSV=%d | JSONL=%d | TST=%d | HOT=%d | Immersion=%d\n', ...
        C.mouse_key{i}, C.Total_Days(i), C.Days_with_CSV(i), C.Days_with_JSONL(i), ...
        C.Days_with_TST(i), C.Days_with_HOT(i), C.Days_with_Immersion(i));
end
fprintf(fid,'\n');

% Issues block
issuesPath = fullfile(OUT_ROOT,'issues.txt');
if ~isempty(ISSUES)
    fprintf(fid,'--- Issues / Warnings (%d) ---\n', numel(ISSUES));
    for i=1:numel(ISSUES), fprintf(fid,'%s\n', ISSUES{i}); end
    fprintf(fid,'\nFull issues list saved to: %s\n', issuesPath);
else
    fprintf(fid,'--- Issues / Warnings ---\n(none)\n');
end
fclose(fid);

% Save full issues file
fid2 = fopen(issuesPath,'w');
for i=1:numel(ISSUES), fprintf(fid2,'%s\n', ISSUES{i}); end
fclose(fid2);

fprintf('Wrote summary:\n  %s\n', summaryPath);
fprintf('Wrote issues:\n  %s\n', issuesPath);

%% ================== Helpers ==================
function OUT_ROOT = nextRunRoot(OUT_BASE)
    d = dir(fullfile(OUT_BASE,'run_*')); nums = 0;
    for i=1:numel(d)
        tok = regexp(d(i).name,'^run_(\d+)$','tokens','once');
        if ~isempty(tok), nums = max(nums, str2double(tok{1})); end
    end
    next = nums + 1;
    OUT_ROOT = fullfile(OUT_BASE, sprintf('run_%03d', next));
    if ~exist(OUT_ROOT,'dir'), mkdir(OUT_ROOT); end
end

function id4 = canonicalCageID(folderName)
    m = regexp(folderName,'\d{4}','match','once'); % first 4 consecutive digits
    if ~isempty(m), id4 = m; else, id4 = folderName(1:min(4,numel(folderName))); end
end

function [csvPath, usedName] = findCombinedFile(sessPath)
    cands = {fullfile(sessPath,'combined_pupil_digital.csv'), fullfile(sessPath,'combined_pupil_digital.xlsx')};
    for i=1:numel(cands)
        if exist(cands{i},'file')
            csvPath = cands{i}; [~,nm,ex]=fileparts(csvPath); usedName=[nm ex]; return;
        end
    end
    % fallback: any CSV that looks like the combined file
    csvList = dir(fullfile(sessPath,'*.csv'));
    for k=1:numel(csvList)
        fp = fullfile(csvList(k).folder, csvList(k).name);
        try
            opts = detectImportOptions(fp,'VariableNamingRule','preserve'); opts.DataLines=[2 201];
            Tprev = readtable(fp, opts);
            if hasCombinedColumns(Tprev), csvPath = fp; usedName = csvList(k).name; return;
            end
        catch
        end
    end
    csvPath=''; usedName='';
end

function tf = hasCombinedColumns(T)
    V = lower(string(T.Properties.VariableNames));
    tf = ismember("frame",V) && any(contains(V,"diameter")) && any(contains(V,"camtime"));
end

function T = safeReadCombined(fp)
    try, T = readtable(fp,'VariableNamingRule','preserve'); catch, T = readtable(fp); end
    T = normalizeCombinedColumns(T);
end

function T = normalizeCombinedColumns(T)
    VN  = T.Properties.VariableNames;
    VNl = lower(string(VN));
    canon = struct( ...
        'Frame',             guessName(VN,VNl,["frame"]), ...
        'PupilTimestamp_s',  guessName(VN,VNl,["pupiltimestamp_s","pupiltimestamp","timestamp","time_s","time"]), ...
        'Diameter_px',       guessName(VN,VNl,["diameter_px","smootheddiameter_px","diameter_pixels","diameter"]), ...
        'CenterX_px',        guessName(VN,VNl,["centerx_px","centerx"]), ...
        'CenterY_px',        guessName(VN,VNl,["centery_px","centery"]), ...
        'UsedScoreFallback', guessName(VN,VNl,["usedscorefallback","usedscore","score"]), ...
        'CamTime_s',         guessName(VN,VNl,["camtime_s","camtime"]), ...
        'CamTime_rel_s',     guessName(VN,VNl,["camtime_rel_s","relative"]), ...
        'Camera_TTL',        guessName(VN,VNl,["camera_ttl","camera"]), ...
        'Lick_TTL',          guessName(VN,VNl,["lick_ttl","lick"]), ...
        'Injector_TTL',      guessName(VN,VNl,["injector_ttl","injector","injector_t"]), ...
        'PlotTime_s_30fps',  guessName(VN,VNl,["plottime_s_30fps","plottime_s","plottime"]) ...
    );

    flds = fieldnames(canon);
    for i=1:numel(flds)
        want = flds{i}; have = canon.(want);
        if ~isempty(have) && ~strcmpi(have, want)
            if exist('renamevars','file')==2
                T = renamevars(T, have, want);
            else
                idx = find(strcmpi(T.Properties.VariableNames, have),1);
                if ~isempty(idx), T.Properties.VariableNames{idx} = want; end
            end
        end
    end

    order = {'Frame','PupilTimestamp_s','Diameter_px','CenterX_px','CenterY_px','UsedScoreFallback', ...
             'CamTime_s','CamTime_rel_s','Camera_TTL','Lick_TTL','Injector_TTL','PlotTime_s_30fps'};
    for i=1:numel(order)
        if ~ismember(order{i}, T.Properties.VariableNames)
            T.(order{i}) = nan(height(T),1);
        end
    end
    rest = setdiff(T.Properties.VariableNames, order, 'stable');
    T = T(:, [order, rest]);
end

function nm = guessName(VN, VNl, keys)
    nm = '';
    for k=1:numel(keys)
        targ = lower(string(keys(k)));
        idx = find(VNl == targ, 1); if ~isempty(idx), nm = char(VN(idx)); return; end
        idx = find(contains(VNl, targ), 1); if ~isempty(idx), nm = char(VN(idx)); return; end
    end
end

function T = putFirst(T, names)
    if isstring(names), names = cellstr(names); end
    present = intersect(names, T.Properties.VariableNames, 'stable');
    rest    = setdiff(T.Properties.VariableNames, present, 'stable');
    T = T(:, [present, rest]);
end

%% -------- JSONL & trials (revised) ----------
function [reqLast, reqByTrial, rewardElapsed, reqByReward, trialStartElapsed, paradigm, isPassive, jsonPath, parseNotes] = ...
    findAndParseJSONL(dayPath, cagePath, mouseDir, sessPath, cageID4, mouseID)

    lists = { dir(fullfile(sessPath,'*.jsonl')), ...
              dir(fullfile(mouseDir,'*.jsonl')), ...
              dir(fullfile(cagePath,'*.jsonl')), ...
              dir(fullfile(dayPath,'*.jsonl')) };
    filt = @(L) L(~cellfun('isempty', regexp(lower({L.name}), ...
                   [regexptranslate('escape', lower(mouseID)) '|' regexptranslate('escape', lower(cageID4))], 'once')));
    cand = [];
    for i=1:numel(lists)
        L = lists{i}; if isempty(L), continue; end
        F = filt(L); if isempty(F), F = L; end
        cand = [cand; F(:)]; %#ok<AGROW>
        if ~isempty(F), break; end
    end
    if isempty(cand)
        reqLast = NaN; reqByTrial = []; rewardElapsed = []; reqByReward = [];
        trialStartElapsed = []; paradigm = "Unknown"; isPassive = 0; jsonPath = ""; parseNotes = "No JSONL found";
        return;
    end
    [~,ix] = max([cand.datenum]);
    jsonPath = fullfile(cand(ix).folder, cand(ix).name);
    [reqLast, reqByTrial, rewardElapsed, reqByReward, trialStartElapsed, paradigm, isPassive, parseNotes] = parseJsonlForRequirements(jsonPath);
end

function [lastReq, reqByTrial, rewardElapsed, reqByReward, trialStartElapsed, paradigm, isPassive, parseNotes] = parseJsonlForRequirements(fp)
    lastReq           = NaN;
    reqByTrial        = [];
    rewardElapsed     = [];
    reqByReward       = [];
    trialStartElapsed = [];
    paradigmFound     = "";
    parseNotes        = "";

    reqStart = containers.Map('KeyType','double','ValueType','double');
    reqEnd   = containers.Map('KeyType','double','ValueType','double');
    curTrial = NaN;
    sawReq   = false;

    fid = fopen(fp,'r');
    if fid < 0
        parseNotes = "Could not open JSONL";
        paradigm   = "Unknown";
        isPassive  = 1;
        return;
    end
    cleaner = onCleanup(@() fclose(fid));

    while true
        ln = fgetl(fid);
        if ~ischar(ln), break; end
        ln = strtrim(ln);
        if isempty(ln), continue; end

        try
            obj = jsondecode(ln);
        catch
            % bad line, skip
            continue;
        end

        if isstruct(obj)
            % paradigm label (first one wins)
            if isfield(obj,'paradigm') && paradigmFound == ""
                paradigmFound = string(obj.paradigm);
            end

            % requirement fields
            if isfield(obj,'requirement')
                rv = toNum(obj.requirement);
                if isfinite(rv), lastReq = double(rv); sawReq = true; end
            end
            if isfield(obj,'req')
                rv = toNum(obj.req);
                if isfinite(rv), lastReq = double(rv); sawReq = true; end
            end

            % trial field
            if isfield(obj,'trial')
                tr = toNum(obj.trial);
                if isfinite(tr), curTrial = double(tr); end
            end

            % event tags
            if isfield(obj,'tag')
                tg = string(obj.tag);

                if tg=="trial_start"
                    % time
                    el = NaN;
                    if isfield(obj,'elapsed'), el = toNum(obj.elapsed); end
                    if isfinite(el)
                        trialStartElapsed(end+1,1) = double(el); %#ok<AGROW>
                    end
                    % requirement snapshot at trial start
                    k = curTrial; if isfield(obj,'trial'), k = double(toNum(obj.trial)); end
                    rv = NaN; if isfield(obj,'req'), rv = toNum(obj.req); end
                    if isfinite(k) && isfinite(rv), reqStart(k) = double(rv); end

                elseif tg=="trial_end"
                    k = curTrial; if isfield(obj,'trial'), k = double(toNum(obj.trial)); end
                    rv = NaN;
                    if isfield(obj,'requirement'), rv = toNum(obj.requirement);
                    elseif isfield(obj,'req'),     rv = toNum(obj.req);
                    end
                    if isfinite(k) && isfinite(rv), reqEnd(k) = double(rv); end

                elseif tg=="reward_cmd"
                    el = NaN;
                    if isfield(obj,'elapsed'), el = toNum(obj.elapsed); end
                    if isfinite(el)
                        rewardElapsed(end+1,1) = double(el); %#ok<AGROW>
                        k = curTrial; if isfield(obj,'trial'), k = double(toNum(obj.trial)); end
                        rv = NaN;
                        if isfinite(k)
                            if isKey(reqEnd,k)
                                rv = reqEnd(k);
                            elseif isKey(reqStart,k)
                                rv = reqStart(k);
                            end
                        end
                        if ~isfinite(rv), rv = NaN; end
                        reqByReward(end+1,1) = double(rv); %#ok<AGROW>
                    end
                end
            end
        end
    end

    % Build per-trial requirement vector (k >= 1), prefer trial_end
    ks = unique([cell2mat(keys(reqStart)), cell2mat(keys(reqEnd))]);
    ks = ks(isfinite(ks));
    ksPos = ks(ks >= 1);
    if ~isempty(ksPos)
        maxK = max(ksPos);
        reqByTrial = nan(1, maxK);
        for k = 1:maxK
            if isKey(reqEnd,k)
                reqByTrial(k) = reqEnd(k);
            elseif isKey(reqStart,k)
                reqByTrial(k) = reqStart(k);
            else
                reqByTrial(k) = NaN;
            end
        end
    end

    % Normalize paradigm + passive flag
    paradigm  = normalizeParadigm(paradigmFound);
    isPassive = 0;
    s = lower(paradigm);
    if contains(s,"passive") || contains(s,"replay")
        isPassive = 1;
    end
    if ~sawReq
        isPassive = 1;
        if paradigm == "Unknown"
            parseNotes = 'No requirement/req fields found — treating as passive';
        end
    end
end


function p = normalizeParadigm(praw)
    p = "Unknown";
    if strlength(praw)==0, return; end
    s = lower(string(praw));
    if contains(s,'progressive')
        p = "Progressive Ratio";
    elseif contains(s,'fixed')
        p = "Fixed Ratio";
    elseif contains(s,'replay') || contains(s,'micro')
        if contains(s,'passive'), p = "MicroLogReplay Passive"; else, p = "MicroLogReplay"; end
    else
        p = string(praw);
    end
end

function jp = findJsonlTopOnly(dayPath, cagePath, mouseDir, cageID4, mouseID)
    lists = { dir(fullfile(mouseDir,'*.jsonl')), dir(fullfile(cagePath,'*.jsonl')), dir(fullfile(dayPath,'*.jsonl')) };
    filt = @(L) L(~cellfun('isempty', regexp(lower({L.name}), ...
                   [regexptranslate('escape', lower(mouseID)) '|' regexptranslate('escape', lower(cageID4))], 'once')));
    cand = [];
    for i=1:numel(lists)
        L = lists{i}; if isempty(L), continue; end
        F = filt(L); if isempty(F), F = L; end
        cand = [cand; F(:)]; %#ok<AGROW>
        if ~isempty(F), break; end
    end
    if isempty(cand), jp = ""; else
        [~,ix] = max([cand.datenum]); jp = fullfile(cand(ix).folder, cand(ix).name);
    end
end

function [trialIdx, trialReq, alignMethod] = assignTrialsSmart(T, reqByTrial, rewardElapsed, reqByReward, trialStartElapsed)
    n = height(T); trialIdx = zeros(n,1); trialReq = nan(n,1); alignMethod = 'none';

    % Branch 1: Injector TTL
    if ismember('Injector_TTL', T.Properties.VariableNames)
        inj = T.Injector_TTL; inj(isnan(inj)) = 0; inj = inj > 0.5;
        onset = find(inj & [false; ~inj(1:end-1)]);
        if ~isempty(onset)
            if onset(1) > 1, trialIdx(1:onset(1)-1) = 0; end
            K = numel(onset);
            for k = 1:K
                i1 = onset(k); if k < K, i2 = onset(k+1)-1; else, i2 = n; end
                trialIdx(i1:i2) = k;
                rv = NaN; if ~isempty(reqByTrial) && k <= numel(reqByTrial), rv = reqByTrial(k); end
                trialReq(i1:i2) = rv;
            end
            alignMethod = 'injector_ttl';
            return;
        end
    end

    tb = pickTimebaseSeconds(T);

    % Branch 2: reward_cmd fallback
    if ~isempty(tb) && ~isempty(rewardElapsed)
        el = double(rewardElapsed(:)); rq = double(reqByReward(:));
        good = isfinite(el); el = el(good); rq = rq(good);
        if ~isempty(el)
            [el, so] = sort(el); rq = rq(so);
            K = numel(el);
            for k = 1:K
                i1 = find(tb >= el(k), 1, 'first'); if isempty(i1), continue; end
                if k < K
                    i2 = find(tb >= el(k+1), 1, 'first') - 1; if isempty(i2), i2 = n; end
                else
                    i2 = n;
                end
                if i2 >= i1
                    trialIdx(i1:i2) = k;
                    rv = NaN; if k <= numel(rq), rv = rq(k); end
                    trialReq(i1:i2) = rv;
                end
            end
            alignMethod = 'reward_cmd';
            return;
        end
    end

    % Branch 3: trial_start fallback
    if ~isempty(tb) && ~isempty(trialStartElapsed)
        el = double(trialStartElapsed(:));
        good = isfinite(el); el = el(good);
        if ~isempty(el)
            el = sort(el); K = numel(el);
            for k = 1:K
                i1 = find(tb >= el(k), 1, 'first'); if isempty(i1), continue; end
                if k < K
                    i2 = find(tb >= el(k+1), 1, 'first') - 1; if isempty(i2), i2 = n; end
                else
                    i2 = n;
                end
                if i2 >= i1
                    trialIdx(i1:i2) = k;
                    rv = NaN; if ~isempty(reqByTrial) && k <= numel(reqByTrial), rv = reqByTrial(k); end
                    trialReq(i1:i2) = rv;
                end
            end
            alignMethod = 'trial_start';
            return;
        end
    end
end

function tb = pickTimebaseSeconds(T)
    cands = {'CamTime_rel_s','PupilTimestamp_s','CamTime_s','PlotTime_s_30fps'};
    tb = [];
    for i=1:numel(cands)
        nm = cands{i};
        if ismember(nm, T.Properties.VariableNames)
            v = double(T.(nm));
            if any(isfinite(v)), tb = v; return; end
        end
    end
end

function x = toNum(v)
    if isnumeric(v) && isscalar(v), x = double(v); return; end
    if isstring(v) || ischar(v)
        s = char(string(v)); s = strrep(s, ',', '.');
        s = regexprep(s, '[^\d\.\-+eE]', '');
        if isempty(s), x = NaN; else, x = str2double(s); end
        return;
    end
    x = NaN;
end

%% -------- Immersion/TST/HOT helpers ----------
%% -------- Immersion/TST/HOT helpers (REVISED, STRICT) ----------
%% -------- Immersion/TST/HOT helpers (REVISED, STRICT) ----------
function [lat, usedFile] = findImmersionLatency(mouseDir)
% Strict: read cell B1 (row 1, col 2) only. No guessing/averaging.
    lat = NaN; usedFile = '';
    L = [dir(fullfile(mouseDir,'*immer*.xlsx')); ...
         dir(fullfile(mouseDir,'*immer*.xls'));  ...
         dir(fullfile(mouseDir,'*immer*.xlsm')); ...
         dir(fullfile(mouseDir,'*immer*.csv'))];
    if isempty(L), return; end
    [~,ix] = max([L.datenum]); fp = fullfile(L(ix).folder, L(ix).name); usedFile = fp;

    try
        if endsWith(lower(fp), '.csv')
            C = readcell(fp, 'DatetimeType','text');   % treat like spreadsheet
        else
            C = readcell(fp, 'DatetimeType','text', 'Sheet', 1);
        end
        lat = parseB1(C);
    catch ME
        warning('Immersion read failed (%s): %s', fp, ME.message);
        lat = NaN;
    end

    % Final guardrails; log out-of-range to console (will flow into issues.txt via caller)
    if ~(isfinite(lat) && lat>=0 && lat<=60)
        fprintf('[immersion] ALERT: value %.3f from %s looks suspicious.\n', lat, fp);
    end
end

function val = parseB1(C)
% Parse B1 (row1,col2) from a readcell() matrix C.
    if isempty(C) || size(C,2) < 2
        val = NaN; return;
    end
    v = C{1,2};
    val = toNumStrict(v);
end

function x = toNumStrict(v)
% Robust scalar parser: numeric → double; string like '1.72 s' → 1.72; comma decimals allowed.
    if isnumeric(v) && isscalar(v)
        x = double(v);
        return;
    end
    if isstring(v) || ischar(v)
        s = char(string(v));
        s = strrep(s, ',', '.');                 % comma decimals → dot
        s = regexprep(s, '[^\d\.\-+eE]', '');    % keep only number chars
        if isempty(s), x = NaN; else, x = str2double(s); end
        return;
    end
    x = NaN;
end

function [cols, usedFile] = summarizeManualScoringMat(mouseDir, pattern, prefix)
% Loads latest ManualScoringResults* file for the day (mouseDir), builds frames/% per behavior.
% Expects variables: score (code per frame), behaviorNames (cellstr), numFrames (scalar).
    cols = struct(); usedFile = '';
    L = dir(fullfile(mouseDir, pattern));
    if isempty(L), return; end
    [~,ix] = max([L.datenum]); fp = fullfile(L(ix).folder, L(ix).name); usedFile = fp;

    try
        S = load(fp);
    catch ME
        warning('[%s] Could not load %s: %s', prefix, fp, ME.message); return;
    end

    % Find fields either at top-level or nested one layer down (robust)
    [score, behaviorNames, numFrames] = peelScoringStruct(S);

    if isempty(score)
        warning('[%s] No "score" vector found in %s', prefix, fp);
        return;
    end
    score = score(:);

    if ~isscalar(numFrames) || ~isfinite(numFrames)
        numFrames = numel(score);
    end
    if numFrames ~= numel(score)
        fprintf('[%s] NOTE: numFrames(%d) != length(score)(%d) → using score length.\n', ...
            prefix, numFrames, numel(score));
        numFrames = numel(score);
    end

    % Behavior names
    if isempty(behaviorNames)
        behaviorNames = {'Nonmoving','Licking','Rearing','Flinching','HindlimbLicking','Jump'};
    elseif isstring(behaviorNames) || ischar(behaviorNames)
        behaviorNames = cellstr(behaviorNames);
    end
    behaviorNames = behaviorNames(:).';
    % Sanitize names (valid MATLAB field names)
    for i=1:numel(behaviorNames)
        behaviorNames{i} = matlab.lang.makeValidName(strtrim(behaviorNames{i}));
        if behaviorNames{i} == "", behaviorNames{i} = sprintf('Behavior%d', i-1); end
    end

    % Determine label range and ensure we have names for all codes 0..B
    uniq = unique(score(~isnan(score)));
    uniq = uniq(isfinite(uniq));
    if isempty(uniq), uniq = 0; end
    B = max([uniq(:); numel(behaviorNames)-1]);
    while numel(behaviorNames) < B+1
        behaviorNames{end+1} = sprintf('Behavior%d', numel(behaviorNames));
    end

    % Tally frames + percentage
    for b = 0:B
        nmPretty = behaviorNames{b+1};
        frames = sum(score == b);
        pct = 100 * frames / max(1, numFrames);
        cols.(sprintf('%s_Frames_%s', prefix, nmPretty)) = double(frames);
        cols.(sprintf('%s_Pct_%s',    prefix, nmPretty)) = double(pct);
    end
    cols.(sprintf('%s_TotalFrames', prefix)) = double(numFrames);
end

function [score, behaviorNames, numFrames] = peelScoringStruct(S)
% Extracts score / behaviorNames / numFrames from top level or a one-level nested struct/table.
    score = []; behaviorNames = []; numFrames = [];
    % try direct
    if isfield(S,'score'),          score = S.score; end
    if isfield(S,'behaviorNames'),  behaviorNames = S.behaviorNames; end
    if isfield(S,'numFrames'),      numFrames = S.numFrames; end
    if ~isempty(score), return; end

    % search one level down
    fns = fieldnames(S);
    for i=1:numel(fns)
        v = S.(fns{i});
        try
            if isstruct(v)
                if isfield(v,'score') && isempty(score), score = v.score; end
                if isfield(v,'behaviorNames') && isempty(behaviorNames), behaviorNames = v.behaviorNames; end
                if isfield(v,'numFrames') && isempty(numFrames), numFrames = v.numFrames; end
            elseif istable(v)
                if any(strcmpi(v.Properties.VariableNames,'score')) && isempty(score), score = v.score; end
                if any(strcmpi(v.Properties.VariableNames,'behaviorNames')) && isempty(behaviorNames), behaviorNames = v.behaviorNames; end
                if any(strcmpi(v.Properties.VariableNames,'numFrames')) && isempty(numFrames), numFrames = v.numFrames; end
            end
        catch
            % ignore bad branches
        end
    end
end

function [T, nAdded] = addScalarCols(T, colsStruct)
% Same signature, but safer with types.
    nAdded = 0;
    if isempty(colsStruct), return; end
    f = fieldnames(colsStruct);
    for i=1:numel(f)
        nm = matlab.lang.makeValidName(f{i});
        if exist('matlab.lang.makeUniqueStrings','file')
            nm = matlab.lang.makeUniqueStrings(nm, T.Properties.VariableNames);
        end
        val = colsStruct.(f{i});
        if ~isscalar(val), val = double(val(1)); else, val = double(val); end
        T.(nm) = repmat(val, height(T), 1);
        nAdded = nAdded + 1;
    end
end


%% -------- Vertcat safety harmonizer ----------
function ALL2 = makeTablesVertcatSafe(ALL)
    canon = ALL{1}.Properties.VariableNames;
    for i=2:numel(ALL)
        extra = setdiff(ALL{i}.Properties.VariableNames, canon, 'stable');
        canon = [canon, extra];
    end
    targetType = containers.Map('KeyType','char','ValueType','char');
    for k=1:numel(canon)
        nm = canon{k}; want = 'double';
        for i=1:numel(ALL)
            if ismember(nm, ALL{i}.Properties.VariableNames)
                v = ALL{i}.(nm);
                if isstring(v) || ischar(v) || (iscell(v) && (isempty(v) || ischar(v{1}) || isstring(v)))
                    want = 'string'; break;
                end
            end
        end
        targetType(nm) = want;
    end
    ALL2 = ALL;
    for i=1:numel(ALL2)
        T = ALL2{i};
        miss = setdiff(canon, T.Properties.VariableNames, 'stable');
        for m=1:numel(miss)
            nm = miss{m};
            if strcmp(targetType(nm),'string'), T.(nm) = repmat("", height(T), 1);
            else, T.(nm) = nan(height(T), 1);
            end
        end
        for k=1:numel(canon)
            nm = canon{k};
            if strcmp(targetType(nm),'string')
                if ~isstring(T.(nm)), T.(nm) = string(T.(nm)); end
            else
                if ~isnumeric(T.(nm))
                    try, T.(nm) = double(T.(nm)); catch, T.(nm) = nan(height(T),1); end
                else
                    if ~isa(T.(nm),'double'), T.(nm) = double(T.(nm)); end
                end
            end
        end
        T = T(:, canon);
        ALL2{i} = T;
    end
end
