function TT = pr_step08_enrich_between_lick(T, tb, TTin)
% For each per-trial row, add BetweenTrialLickRate = licks in (reward_{j-1}, start_j) / duration
% First PR trial in a session has NaN. Requires TTin.time_to_reward, Trial, and T with Trial, Lick_TTL.

if isempty(TTin) || height(TTin) < 1
    TT = TTin;
    return
end
T = T;
T.mouse_key = cellstr(string(T.mouse_key));
if ~isnumeric(T.day_index),   T.day_index   = double(T.day_index);   end
if ~isnumeric(T.session_idx), T.session_idx = double(T.session_idx); end
if ~isnumeric(T.Trial),       T.Trial       = double(T.Trial); end
L = T.Lick_TTL > 0.5;
tb = double(tb(:));
TT = TTin;
TT.BetweenTrialLickRate = nan(height(TT), 1);
TT.BetweenTrialNlicks = nan(height(TT), 1);
TT.mouse_key = cellstr(string(TT.mouse_key));
if ~isnumeric(TT.day_index),   TT.day_index   = double(TT.day_index);   end
if ~isnumeric(TT.session_idx), TT.session_idx = double(TT.session_idx); end
if ~isnumeric(TT.Trial),       TT.Trial       = double(TT.Trial); end
if ~ismember('time_to_reward', TT.Properties.VariableNames)
    warning('per_trial has no time_to_reward; between-trial rate not computed.');
    return
end
g = unique(TT(:, {'mouse_key', 'day_index', 'session_idx'}), 'rows', 'stable');
for i = 1:height(g)
    mk = g.mouse_key(i);
    di = g.day_index(i);
    se = g.session_idx(i);
    r0 = strcmp(string(T.mouse_key), string(mk)) & (T.day_index == di) & (T.session_idx == se);
    if ~any(r0)
        continue
    end
    trList = unique(T.Trial(r0 & isfinite(T.Trial)));
    trList = sort(trList);
    for j = 2:numel(trList)
        tK = trList(j - 1);
        tK1 = trList(j);
        rK = r0 & (T.Trial == tK);
        rK1 = r0 & (T.Trial == tK1);
        if ~any(rK) || ~any(rK1)
            continue
        end
        ttk = tb(rK);
        ttk1 = tb(rK1);
        triS = min(ttk(isfinite(ttk)));
        triS1 = min(ttk1(isfinite(ttk1)));
        rt = TT.time_to_reward(strcmp(TT.mouse_key, mk) & (TT.day_index == di) & (TT.session_idx == se) & (TT.Trial == tK));
        if isempty(rt)
            continue
        end
        rt = double(rt(1));
        rewardK = triS + rt;
        durRaw = triS1 - rewardK;
        if ~(durRaw > 1e-6)
            continue
        end
        % Floor duration so a tiny inter-trial interval does not explode licks/s.
        minDurS = 0.02;
        dur = max(durRaw, minDurS);
        licks = L & (tb > rewardK) & (tb < triS1) & r0;
        nLk = double(sum(licks));
        rate = nLk / dur;
        u = find(strcmp(TT.mouse_key, mk) & (TT.day_index == di) & (TT.session_idx == se) & (TT.Trial == tK1), 1);
        if ~isempty(u)
            TT.BetweenTrialLickRate(u) = rate;
            TT.BetweenTrialNlicks(u) = nLk;
        end
    end
end
end
