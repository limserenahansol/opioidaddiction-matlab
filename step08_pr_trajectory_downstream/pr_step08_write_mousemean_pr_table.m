function pr_step08_write_mousemean_pr_table(U, outDir)
% pr_step08_write_mousemean_pr_table — per mouse x PR trial means for Step 9 k-means
% Same schema as run_pr_trial_licking_trajectory: MouseMean_LickRate_by_PRLevel.csv
W = U;
if ismember('phase', W.Properties.VariableNames)
    W = W(~strcmp(W.phase, '(drop)'), :);
end
if height(W) < 1
    W = U;
end
if height(W) < 1
    return
end
G = findgroups(W.mouse_key, W.Trial);
summaryMouse = table;
summaryMouse.mouse_key = cellstr(splitapply(@(x) x(1), W.mouse_key, G));
summaryMouse.Trial = splitapply(@(x) x(1), double(W.Trial), G);
summaryMouse.n = splitapply(@numel, W.WithinTrialLickRate, G);
summaryMouse.MeanWithinTrialLickRate = splitapply(@(x) mean(double(x), 'omitnan'), W.WithinTrialLickRate, G);
if ismember('BetweenTrialLickRate', W.Properties.VariableNames)
    summaryMouse.MeanBetweenTrialLickRate = splitapply(@(x) mean(double(x), 'omitnan'), W.BetweenTrialLickRate, G);
else
    summaryMouse.MeanBetweenTrialLickRate = nan(height(summaryMouse), 1);
end
if ismember('Group', W.Properties.VariableNames)
    summaryMouse.Group = cellstr(splitapply(@(x) x(1), W.Group, G));
else
    summaryMouse.Group = repmat({'All'}, height(summaryMouse), 1);
end
fmm = fullfile(outDir, 'MouseMean_LickRate_by_PRLevel.csv');
writetable(summaryMouse, fmm);
fprintf('Wrote: %s\n', fmm);
end
