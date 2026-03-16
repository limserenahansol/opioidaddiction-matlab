% -------------------------------------------------------------------------
% MATLAB Code: Regression of Tail Immersion vs PK (Daily Aggregation)
% 1. Aggregates behavioral data by DAY first (handling duplicate rows).
% 2. Metrics (8 Total): 
%    - Last 3 Days (Avg, Peak, Low)
%    - Last 5 Days (Avg)
%    - Specific Days (10, 9, 8, 7)
% 3. Matches with PK data using Smart ID Matching (Two Sheets).
% 4. Generates 32 Regression Plots.
% -------------------------------------------------------------------------

function regress_Immersion_vs_PK_DailyMetrics
    clear; close all; clc;

    %% 1. Load Behavioral Data
    fprintf('--- Step 1: Loading Behavioral Data ---\n');
    
    % --- UPDATE PATH IF NEEDED ---
    rootTry = 'K:\addiction_concate_Dec_2025\longitudinal_outputs\';
    % -----------------------------
    
    if ~exist(rootTry, 'dir')
        error('Behavioral data folder not found: %s', rootTry);
    end
    
    d = dir(fullfile(rootTry, 'run_*'));
    if isempty(d)
        error('No run_* folders found in %s', rootTry);
    end
    [~, ix] = max([d.datenum]);
    runDir = fullfile(d(ix).folder, d(ix).name);
    csvPath = fullfile(runDir, 'ALL_mice_longitudinal.csv');
    
    if ~exist(csvPath, 'file')
        error('CSV not found: %s', csvPath);
    end
    
    fprintf('Reading Behavior CSV: %s\n', csvPath);
    T_beh = readtable(csvPath, 'VariableNamingRule', 'preserve');
    
    % Find Immersion Column
    cols = T_beh.Properties.VariableNames;
    immColIdx = find(contains(cols, 'Immersion', 'IgnoreCase', true) & contains(cols, 'Latency', 'IgnoreCase', true), 1);
    
    if isempty(immColIdx)
        error('Could not find Immersion Latency column in CSV.');
    end
    immCol = cols{immColIdx};
    
    % --- CALCULATE METRICS PER MOUSE ---
    uniqueMice = unique(T_beh.mouse_key);
    
    % Initialize struct for metrics
    BehData = struct('key', {}, 'number', {}, ...
        'Avg_Last3', {}, 'Peak_Last3', {}, 'Low_Last3', {}, ...
        'Avg_Last5', {}, ...
        'Day10', {}, 'Day9', {}, 'Day8', {}, 'Day7', {});
    
    count = 0;
    for i = 1:numel(uniqueMice)
        mKey = string(uniqueMice(i));
        rows = T_beh(string(T_beh.mouse_key) == mKey, :);
        
        rawVals = rows.(immCol);
        rawDays = rows.day_index;
        
        % Ensure numeric days
        if iscell(rawDays) || isstring(rawDays), rawDays = str2double(string(rawDays)); end
        
        % Remove NaNs
        validMask = ~isnan(rawVals) & ~isnan(rawDays);
        rawVals = rawVals(validMask);
        rawDays = rawDays(validMask);
        
        if isempty(rawVals), continue; end
        
        % --- AGGREGATE BY DAY ---
        % Handle multiple rows for the same day (take mean of that day)
        uniqueDays = unique(rawDays);
        dailyVals = nan(size(uniqueDays));
        for d = 1:length(uniqueDays)
            dayID = uniqueDays(d);
            dailyVals(d) = mean(rawVals(rawDays == dayID));
        end
        
        % Sort by Day (Ascending)
        [uniqueDays, sortIdx] = sort(uniqueDays);
        dailyVals = dailyVals(sortIdx);
        
        % Extract Numeric ID for matching
        digits = regexp(char(mKey), '\d+', 'match');
        if ~isempty(digits)
            count = count + 1;
            BehData(count).key    = lower(char(mKey));
            BehData(count).number = digits{1};
            
            % --- CALCULATE METRICS ---
            
            % 1, 2, 3: Last 3 Days
            if length(dailyVals) >= 3
                last3 = dailyVals(end-2:end);
            else
                last3 = dailyVals;
            end
            BehData(count).Avg_Last3  = mean(last3);
            BehData(count).Peak_Last3 = max(last3);
            BehData(count).Low_Last3  = min(last3);
            
            % 4: Last 5 Days
            if length(dailyVals) >= 5
                last5 = dailyVals(end-4:end);
            else
                last5 = dailyVals;
            end
            BehData(count).Avg_Last5 = mean(last5);
            
            % 5, 6, 7, 8: Specific Days
            BehData(count).Day10 = getDayVal(uniqueDays, dailyVals, 10);
            BehData(count).Day9  = getDayVal(uniqueDays, dailyVals, 9);
            BehData(count).Day8  = getDayVal(uniqueDays, dailyVals, 8);
            BehData(count).Day7  = getDayVal(uniqueDays, dailyVals, 7);
        end
    end
    fprintf('Computed metrics for %d mice.\n', count);

    %% 2. Load PK Data (Two Sheets)
    fprintf('--- Step 2: Loading PK Data ---\n');
    excelPath = 'C:\Users\hsollim\Documents\data_PKassay.xlsx';
    
    if ~exist(excelPath, 'file')
        error('PK Excel file not found: %s', excelPath);
    end
    
    % --- Load Blood (Sheet 1) ---
    fprintf('Reading Sheet 1 (Blood)...\n');
    try
        rawBlood = readcell(excelPath, 'Sheet', 1);
    catch
        error('Could not read Sheet 1. Ensure the Excel file has at least 2 sheets.');
    end
    
    rowIdxBlood = find(contains(string(rawBlood(:,2)), 'Sample Name'), 1);
    if isempty(rowIdxBlood)
        error('Could not find "Sample Name" header in Sheet 1 (Blood).');
    end
    [Names_Blood, Blood_Morphine, Blood_M3G] = extractTableData(rawBlood, rowIdxBlood);
    
    % --- Load Brain (Sheet 2) ---
    fprintf('Reading Sheet 2 (Brain)...\n');
    try
        rawBrain = readcell(excelPath, 'Sheet', 2);
    catch
        error('Could not read Sheet 2. Ensure the Excel file has at least 2 sheets.');
    end
    
    rowIdxBrain = find(contains(string(rawBrain(:,2)), 'Sample Name'), 1);
    if isempty(rowIdxBrain)
        error('Could not find "Sample Name" header in Sheet 2 (Brain).');
    end
    [Names_Brain, Brain_Morphine, Brain_M3G] = extractTableData(rawBrain, rowIdxBrain);
    
    % Validation
    if length(Names_Blood) ~= length(Names_Brain)
        warning('Number of Blood samples (%d) does not match Brain samples (%d). Matching might be misaligned.', length(Names_Blood), length(Names_Brain));
    end

    %% 3. Smart Matching
    fprintf('--- Step 3: Matching Data ---\n');
    
    N = length(Names_Blood);
    MatchedData = table();
    MatchedData.SampleName = Names_Blood;
    MatchedData.Blood_Morphine = Blood_Morphine;
    MatchedData.Blood_M3G = Blood_M3G;
    
    % Handle Brain data mapping (assuming row order might differ or just using Blood index if they match)
    lenBr = length(Names_Brain);
    if lenBr >= N
        MatchedData.Brain_Morphine = Brain_Morphine(1:N);
        MatchedData.Brain_M3G = Brain_M3G(1:N);
    else
        % Pad with NaN if brain has fewer rows
        MatchedData.Brain_Morphine = [Brain_Morphine; nan(N-lenBr,1)];
        MatchedData.Brain_M3G = [Brain_M3G; nan(N-lenBr,1)];
    end
    
    % Initialize metric columns
    MatchedData.Lat_Avg3  = nan(N, 1);
    MatchedData.Lat_Peak3 = nan(N, 1);
    MatchedData.Lat_Low3  = nan(N, 1);
    MatchedData.Lat_Avg5  = nan(N, 1);
    MatchedData.Lat_Day10 = nan(N, 1);
    MatchedData.Lat_Day9  = nan(N, 1);
    MatchedData.Lat_Day8  = nan(N, 1);
    MatchedData.Lat_Day7  = nan(N, 1);
    
    MatchedData.Group = strings(N, 1);
    
    tags = {'black', 'blue', 'norm', 'red', 'green', 'orange', 'white'};
    
    for i = 1:N
        pkName = lower(char(Names_Blood(i)));
        
        % Group
        if contains(pkName, 'active'), MatchedData.Group(i) = "Active";
        else, MatchedData.Group(i) = "Passive"; end
        
        % Extract ID & Tag
        digits = regexp(pkName, '\d+', 'match');
        if isempty(digits), continue; end
        pkID = digits{1};
        
        pkTag = '';
        for t = tags
            if contains(pkName, t{1}), pkTag = t{1}; break; end
        end
        
        % Find Match
        candidatesIdx = find(strcmp({BehData.number}, pkID));
        bestMatchIdx = [];
        
        if length(candidatesIdx) == 1
            bestMatchIdx = candidatesIdx(1);
        elseif length(candidatesIdx) > 1
            for c = candidatesIdx
                behKey = BehData(c).key;
                if ~isempty(pkTag) && contains(behKey, pkTag)
                    bestMatchIdx = c; break;
                elseif strcmp(pkTag, 'norm') && ~contains(behKey, {'black','blue','red','green','orange','white'})
                    bestMatchIdx = c; break;
                end
            end
        end
        
        if ~isempty(bestMatchIdx)
            MatchedData.Lat_Avg3(i)  = BehData(bestMatchIdx).Avg_Last3;
            MatchedData.Lat_Peak3(i) = BehData(bestMatchIdx).Peak_Last3;
            MatchedData.Lat_Low3(i)  = BehData(bestMatchIdx).Low_Last3;
            MatchedData.Lat_Avg5(i)  = BehData(bestMatchIdx).Avg_Last5;
            MatchedData.Lat_Day10(i) = BehData(bestMatchIdx).Day10;
            MatchedData.Lat_Day9(i)  = BehData(bestMatchIdx).Day9;
            MatchedData.Lat_Day8(i)  = BehData(bestMatchIdx).Day8;
            MatchedData.Lat_Day7(i)  = BehData(bestMatchIdx).Day7;
        end
    end
    
    % Keep matched rows (at least one metric valid)
    MatchedData = MatchedData(~isnan(MatchedData.Lat_Avg3), :);
    fprintf('Matched %d samples.\n', height(MatchedData));

    %% 4. Generate 32 Regression Plots
    % Create New Folder
    outDir = fullfile(runDir, 'figs', 'PK_Regression_DailyMetrics');
    if ~exist(outDir, 'dir'), mkdir(outDir); end
    
    % Define the 8 Metrics
    metricsList = {
        'Avg_Last3',    MatchedData.Lat_Avg3,  'Avg Immersion Latency (Last 3 Days)';
        'Peak_Last3',   MatchedData.Lat_Peak3, 'Peak Immersion Latency (Last 3 Days)';
        'Low_Last3',    MatchedData.Lat_Low3,  'Lowest Immersion Latency (Last 3 Days)';
        'Avg_Last5',    MatchedData.Lat_Avg5,  'Avg Immersion Latency (Last 5 Days)';
        'Day10',        MatchedData.Lat_Day10, 'Immersion Latency (Day 10)';
        'Day9',         MatchedData.Lat_Day9,  'Immersion Latency (Day 9)';
        'Day8',         MatchedData.Lat_Day8,  'Immersion Latency (Day 8)';
        'Day7',         MatchedData.Lat_Day7,  'Immersion Latency (Day 7)';
    };
    
    % Define Analytes
    plotsInfo = {
        'Blood', 'Morphine', MatchedData.Blood_Morphine;
        'Blood', 'M3G',      MatchedData.Blood_M3G;
        'Brain', 'Morphine', MatchedData.Brain_Morphine;
        'Brain', 'M3G',      MatchedData.Brain_M3G
    };
    
    % Loop: 8 Metrics x 4 Analytes = 32 Plots
    for m = 1:size(metricsList, 1)
        mName = metricsList{m, 1};
        xData = metricsList{m, 2};
        xLbl  = metricsList{m, 3};
        
        for k = 1:size(plotsInfo, 1)
            tissue = plotsInfo{k, 1};
            analyte = plotsInfo{k, 2};
            yData = plotsInfo{k, 3};
            
            if strcmpi(tissue, 'Blood'), unit = 'ng/mL'; else, unit = 'ng/g'; end
            
            % Generate Plot
            makeRegPlot(xData, yData, MatchedData.Group, MatchedData.SampleName, ...
                        tissue, analyte, unit, mName, xLbl, outDir);
        end
    end
    
    fprintf('Done! 32 Figures saved to: %s\n', outDir);
end

%% Helper Functions

function val = getDayVal(days, vals, targetDay)
    idx = find(days == targetDay, 1);
    if ~isempty(idx)
        val = vals(idx);
    else
        val = NaN;
    end
end

function makeRegPlot(x, y, groups, names, tissue, analyte, unit, metricName, xLbl, outDir)
    % Filter NaNs (Critical for specific days like Day 10 which might be missing)
    mask = ~isnan(x) & ~isnan(y);
    x = x(mask); y = y(mask); groups = groups(mask); names = names(mask);
    
    if isempty(x), return; end

    f = figure('Color','w','Position',[100 100 800 600], 'Visible', 'off'); 
    set(f, 'ToolBar', 'none'); % Suppress toolbar warning
    hold on;
    
    cBlue = [0 0 1]; cRed = [1 0 0];
    gAct = strcmp(groups, "Active");
    gPas = strcmp(groups, "Passive");
    
    % Scatter
    scatter(x(gAct), y(gAct), 100, 'o', 'filled', 'MarkerFaceColor', cBlue, 'MarkerEdgeColor', 'k', 'DisplayName', 'Active');
    scatter(x(gPas), y(gPas), 100, 'o', 'filled', 'MarkerFaceColor', cRed, 'MarkerEdgeColor', 'k', 'DisplayName', 'Passive');
    
    % Labels (Name)
    text(x, y, names, 'VerticalAlignment','bottom', 'HorizontalAlignment','left', ...
        'FontSize', 8, 'Interpreter', 'none', 'Color', [0.3 0.3 0.3]);
    
    % Regression
    if length(x) > 1
        mdl = fitlm(x, y);
        xFit = linspace(min(x)*0.9, max(x)*1.1, 100)';
        yFit = predict(mdl, xFit);
        plot(xFit, yFit, 'k-', 'LineWidth', 2, 'DisplayName', 'Linear Fit');
        
        % Stats
        R2 = mdl.Rsquared.Ordinary;
        pval = mdl.Coefficients.pValue(2);
        
        % Significance
        if pval < 0.05
            sigStr = '(SIGNIFICANT)'; titleCol = 'r';
        else
            sigStr = '(Not Significant)'; titleCol = 'k';
        end
        titleText = sprintf('[%s] %s %s vs Immersion\nR^2 = %.3f, p = %.3f %s', ...
            metricName, tissue, analyte, R2, pval, sigStr);
    else
        titleText = sprintf('[%s] %s %s (Not enough points)', metricName, tissue, analyte);
        titleCol = 'k';
    end
    
    % Labels
    xlabel(xLbl, 'FontSize', 12);
    ylabel(sprintf('%s %s [%s]', tissue, analyte, unit), 'FontSize', 12);
    title(titleText, 'FontSize', 14, 'Interpreter', 'none', 'Color', titleCol);
    
    legend([findobj(gca,'DisplayName','Active'), findobj(gca,'DisplayName','Passive')], 'Location', 'best');
    grid on; box on;
    
    % Save
    fileName = fullfile(outDir, sprintf('Regress_%s_%s_%s.png', tissue, analyte, metricName));
    exportgraphics(f, fileName, 'Resolution', 300);
    close(f);
end

function [sampleNames, morphine, m3g] = extractTableData(rawMatrix, headerRowIdx)
    headers = string(rawMatrix(headerRowIdx, :));
    col_Name = find(contains(headers, 'Sample Name', 'IgnoreCase', true), 1);
    col_M = find(contains(headers, 'Morphine') & ~contains(headers, 'glucuronide'), 1);
    col_3G = find(contains(headers, 'Morphine-3'), 1);
    
    sampleNames = []; morphine = []; m3g = [];
    if isempty(col_Name) || isempty(col_M) || isempty(col_3G), return; end
    
    r = headerRowIdx + 1;
    while r <= size(rawMatrix, 1)
        valName = rawMatrix{r, col_Name};
        
        % STOP CONDITION:
        % Stop if empty, missing, NaN (end of file)
        % OR if we hit the "Sample Name" header of the NEXT table (e.g. Brain table)
        if (isnumeric(valName) && isnan(valName)) || isempty(valName) || any(ismissing(valName))
            break; 
        end
        if contains(string(valName), 'Sample Name', 'IgnoreCase', true)
            break; 
        end
        
        sampleNames = [sampleNames; string(valName)]; %#ok<AGROW>
        morphine = [morphine; cleanNumber(rawMatrix{r, col_M})]; %#ok<AGROW>
        m3g = [m3g; cleanNumber(rawMatrix{r, col_3G})]; %#ok<AGROW>
        r = r + 1;
    end
end

function num = cleanNumber(x)
    if isnumeric(x), num = x;
    elseif ischar(x) || isstring(x)
        if contains(x, 'LLOQ', 'IgnoreCase', true), num = 0;
        else, num = str2double(x); end
    else, num = NaN; end
end