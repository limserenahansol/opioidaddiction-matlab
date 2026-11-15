"""24 frames per session with the CURRENT roi.mat (window = ROI expanded 2x), at times offset half a step from the
earlier sample set (other frames). Writes <key>.mat (samples HxWxN uint8, roi, roi_in_crop) and meta.json."""
import json, subprocess, sys
from pathlib import Path
import imageio_ffmpeg, numpy as np, scipy.io as sio
PROJ, OUT = Path(sys.argv[1]), Path(sys.argv[2]); N = int(sys.argv[3]) if len(sys.argv) > 3 else 24
OUT.mkdir(parents=True, exist_ok=True)
FF = imageio_ffmpeg.get_ffmpeg_exe(); W, H = 1920, 1080
def grab(video, t):
    raw = subprocess.run([FF, "-hide_banner", "-loglevel", "error", "-ss", f"{t:.3f}", "-i", str(video), "-frames:v", "1",
                          "-f", "rawvideo", "-pix_fmt", "gray", "-"], capture_output=True, check=True).stdout
    return np.frombuffer(raw[: W * H], np.uint8).reshape(H, W)
def duration(video):
    r = subprocess.run([FF, "-hide_banner", "-i", str(video)], capture_output=True, text=True).stderr
    hh, mm, ss = r.split("Duration: ")[1].split(",")[0].split(":"); return int(hh) * 3600 + int(mm) * 60 + float(ss)
meta = []
for v in sorted(PROJ.glob("video/day*/*/*/*.mp4")):
    day, cage, mouse = v.parts[-4], v.parts[-3], v.parts[-2]; key = f"{day}_{cage}_{mouse}"
    roi = sio.loadmat(v.parent / "roi.mat")["roi"].ravel().astype(float); x, y, w, h = roi
    cx, cy = x - 1 + w / 2, y - 1 + h / 2
    ex0, ey0 = int(max(0, cx - w)), int(max(0, cy - h)); ex1, ey1 = int(min(W, cx + w)), int(min(H, cy + h))
    dur = duration(v); step = (dur - 15) / (N - 1)
    ts = np.linspace(7.5, dur - 7.5, N) + step / 2; ts = np.minimum(ts, dur - 1.0)
    crops = np.stack([grab(v, t)[ey0:ey1, ex0:ex1] for t in ts], axis=-1)
    rx, ry = x - 1 - ex0, y - 1 - ey0
    sio.savemat(OUT / f"{key}.mat", {"samples": crops, "roi": roi[None, :], "roi_in_crop": np.array([[rx, ry, w, h]])},
                do_compression=True)
    meta.append(dict(key=key, video=str(v), roi=roi.tolist(), crop_origin=[ex0, ey0], roi_in_crop=[rx, ry, w, h],
                     duration_container_s=dur, sample_t=ts.tolist()))
    print("done", key, crops.shape, flush=True)
(OUT / "meta.json").write_text(json.dumps(meta, indent=1))
