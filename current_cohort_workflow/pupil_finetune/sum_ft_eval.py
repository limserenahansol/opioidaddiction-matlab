"""Sum the leave-mice-out evaluation of the fine-tuned vs original U-Net over the folds (eval_fold*.csv)."""
import sys
from pathlib import Path
import pandas as pd
d = Path(sys.argv[1])
E = pd.concat([pd.read_csv(f) for f in sorted(d.glob("eval_fold*.csv"))], ignore_index=True)
print(f"held-out frames: {len(E)} ({E.session.nunique()} sessions)")
for m in ("ft", "orig"):
    c = E[m].value_counts().to_dict(); lab = E[E.label >= 1]; npk = int(E[m].isin(["correct", "wrong", "false", "blob_unknown"]).sum())
    bad = c.get("wrong", 0) + c.get("false", 0)
    print(f"{m:5s} pupil on contour ({len(lab)}): correct {c.get('correct', 0)}, wrong {c.get('wrong', 0)}, miss {c.get('miss', 0)} "
          f"(recall {(lab[m] == 'correct').mean():.1%}) | not visible ({(E.label == -1).sum()}): none-ok {c.get('none_ok', 0)}, "
          f"false {c.get('false', 0)} | visible, no contour ({(E.label == 0).sum()}): blob {c.get('blob_unknown', 0)}, none {c.get('none_unknown', 0)} "
          f"| bad among known picks {bad}/{npk - c.get('blob_unknown', 0)}")
    r = E.loc[(E.label >= 1) & (E[m] == "correct"), f"{m}_diam_ratio"]
    print(f"      diameter ratio (pick / labelled contour) for correct picks: median {r.median():.2f}, IQR {r.quantile(.25):.2f}-{r.quantile(.75):.2f}")
print(E.groupby("session")[["ft", "orig"]].agg(lambda s: f"{(s == 'correct').sum()}c/{(s == 'wrong').sum()}w/{(s == 'miss').sum()}m/{(s == 'false').sum()}f").to_string())
