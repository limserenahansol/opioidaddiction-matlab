function mk = step09_normalize_mouse_key(x)
% Align with compute_addiction_index_EFA_mousefit_projectDays_v3 normalize_mouse_key
% + spaces -> underscore so "6100 black" matches "6100_black"
if isempty(x)
    mk = x;
    return
end
if iscategorical(x)
    s = string(x);
elseif iscell(x)
    s = string(x);
else
    s = string(x);
end
s = strtrim(s);
s = replace(s, "\", "_");
s = replace(s, "/", "_");
s = replace(s, "-", "_");
s = replace(s, " ", "_");
% collapse duplicate underscores (simple)
for i = 1:numel(s)
    t = s(i);
    while contains(t, "__")
        t = replace(t, "__", "_");
    end
    s(i) = t;
end
mk = s;
for i = 1:numel(mk)
    parts = split(mk(i), "_");
    parts = parts(strlength(parts) > 0);
    if numel(parts) >= 2
        p1 = parts(1);
        p2 = lower(parts(2));
        mk(i) = p1 + "_" + p2;
    end
end
end
