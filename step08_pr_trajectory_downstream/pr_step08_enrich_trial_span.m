function TT = pr_step08_enrich_trial_span(T, tb, TTin)
% For each per-trial row, add trial_span_s = max(CamTime)-min(CamTime) for that
% (mouse, day, session, Trial) in the longitudinal table (full sampled trial width).
% Vectorized: one O(n) pass on the long table, then ismember to TT (not a row-by-row scan of T).

if isempty(TTin) || height(TTin) < 1
    TT = TTin;
    return
end
T = T;
T.mouse_key = cellstr(string(T.mouse_key));
if ~isnumeric(T.day_index),   T.day_index   = double(T.day_index);   end
if ~isnumeric(T.session_idx), T.session_idx = double(T.session_idx); end
if ~isnumeric(T.Trial),       T.Trial       = double(T.Trial); end
tb = double(tb(:));
nT = height(T);
if numel(tb) ~= nT
    error('Timebase must match longitudinal table height (got %d vs %d).', numel(tb), nT);
end
kLong = strcat(string(T.mouse_key), char(0), string(T.day_index), char(0), string(T.session_idx), char(0), string(T.Trial));
G = findgroups(kLong);
spanG = splitapply(@pr_step08_local_trial_span_range, tb, G);
% One representative key per group, same order 1:max(G) as splitapply
uK = splitapply(@(c) c(1), kLong, G);
TT = TTin;
if ~isnumeric(TT.Trial), TT.Trial = double(TT.Trial); end
kRow = strcat(string(TT.mouse_key), char(0), string(TT.day_index), char(0), string(TT.session_idx), char(0), string(TT.Trial));
[tf, loc] = ismember(kRow, uK);
spanV = nan(height(TT), 1);
ok = tf & (loc > 0);
spanV(ok) = spanG(loc(ok));
TT.trial_span_s = spanV;
end

function s = pr_step08_local_trial_span_range(t)
t = double(t);
t = t(isfinite(t));
if isempty(t)
    s = NaN;
else
    s = max(t) - min(t);
end
end
