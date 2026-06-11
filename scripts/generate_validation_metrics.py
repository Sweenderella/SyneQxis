#!/usr/bin/env python3
"""Generate SyneQxis/Q1K validation metrics for the tool-release paper.

Run on Narval from the project root:
    python3 paper_validation_metrics.py --project-path /lustre07/scratch/rsweety/white_paper/wd

Outputs CSVs and a LaTeX snippet under paper_validation_metrics/.
The script computes only metrics supported by available files; unavailable
metrics are reported as NA with a reason.
"""

from __future__ import annotations

import argparse
import math
import re
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd


TASKS = ("GO", "PLR", "VEP")
SITES = ("HSJ", "MHC")


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--project-path", default="/lustre07/scratch/rsweety/white_paper/wd")
    p.add_argument("--outdir", default="paper_validation_metrics")
    p.add_argument("--skip-heavy", action="store_true", help="Skip MNE-based epoch/raw scans")
    return p.parse_args()


def norm_site(x):
    if pd.isna(x):
        return np.nan
    x = str(x).strip().upper()
    return {"MNI": "MHC"}.get(x, x)


def mean_sd(x):
    x = pd.to_numeric(pd.Series(x), errors="coerce").dropna()
    if len(x) == 0:
        return np.nan, np.nan
    return float(x.mean()), float(x.std(ddof=1)) if len(x) > 1 else np.nan


def pct(n, d):
    return np.nan if d == 0 else 100.0 * n / d


def cohen_d_paired(a, b):
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    mask = np.isfinite(a) & np.isfinite(b)
    diff = b[mask] - a[mask]
    if len(diff) < 2 or diff.std(ddof=1) == 0:
        return np.nan
    return float(diff.mean() / diff.std(ddof=1))


def load_demographics(project):
    candidates = [
        project / "participants.tsv",
        project / "q1k_new" / "participants.tsv",
        project / "derivatives" / "init" / "GO" / "participants.tsv",
    ]
    for path in candidates:
        if path.exists():
            df = pd.read_csv(path, sep="\t")
            if len(df) and {"participant_id", "sex"}.issubset(df.columns):
                return df, path
    return pd.DataFrame(), None


def load_inventory(project):
    path = project / "q1k_manager_inventory" / "q1k_master_inventory.csv"
    if path.exists():
        return pd.read_csv(path), path
    return pd.DataFrame(), None


def demographics_table(participants, inventory):
    if participants.empty:
        return pd.DataFrame()
    df = participants.copy()
    df["site_clean"] = df.get("site", pd.Series(index=df.index, dtype=object)).map(norm_site)
    age_col = "eeg_age" if "eeg_age" in df.columns else "age"
    df["age_num"] = pd.to_numeric(df.get(age_col), errors="coerce")
    df["female"] = df.get("sex", "").astype(str).str.lower().str.startswith("female")

    rows = []
    for site in list(SITES) + ["TOTAL"]:
        d = df[df["site_clean"].isin(SITES)] if site == "TOTAL" else df[df["site_clean"] == site]
        age_m, age_s = mean_sd(d["age_num"])
        row = {
            "site": site,
            "N_demographics": len(d),
            "age_mean": age_m,
            "age_sd": age_s,
            "female_n": int(d["female"].sum()) if len(d) else 0,
            "female_pct": pct(int(d["female"].sum()), len(d)),
        }
        if "group" in d.columns:
            for g in ("proband", "mother", "father", "sibling"):
                row[f"group_{g}_n"] = int((d["group"].astype(str).str.lower() == g).sum())
        if "ndd" in d.columns:
            row["ndd_n"] = int((d["ndd"].astype(str).str.lower() == "ndd").sum())
            row["no_ndd_n"] = int((d["ndd"].astype(str).str.lower() == "no_ndd").sum())
        if "asd" in d.columns:
            asd = d["asd"].astype(str).str.lower()
            row["asd_confirmed_or_1_n"] = int(((asd == "1") | asd.str.contains("confirmed", na=False)).sum())
            row["asd_suspected_n"] = int(asd.str.contains("suspected", na=False).sum())
        rows.append(row)

    demo = pd.DataFrame(rows)

    if not inventory.empty:
        inv = inventory[inventory["site"].isin(SITES)].copy()
        for site in list(SITES) + ["TOTAL"]:
            d = inv if site == "TOTAL" else inv[inv["site"] == site]
            eeg_subs = set(d.loc[d["source_mff_count"] > 0, "subject"])
            et_subs = set(d.loc[d["source_et_count"] > 0, "subject"])
            idx = demo["site"] == site
            demo.loc[idx, "unique_subjects_with_source_EEG"] = len(eeg_subs)
            demo.loc[idx, "unique_subjects_with_source_ET"] = len(et_subs)
            demo.loc[idx, "unique_subjects_with_source_EEG_ET_pair"] = len(eeg_subs & et_subs)
    return demo


def source_inventory_table(inventory):
    if inventory.empty:
        return pd.DataFrame()
    inv = inventory[inventory["site"].isin(SITES)].copy()
    rows = []
    for site in list(SITES) + ["TOTAL"]:
        d = inv if site == "TOTAL" else inv[inv["site"] == site]
        rows.append({
            "site": site,
            "expected_site_task_subject_rows": len(d),
            "source_EEG_rows": int((d["source_mff_count"] > 0).sum()),
            "source_ET_rows": int((d["source_et_count"] > 0).sum()),
            "complete_source_EEG_ET_rows": int(((d["source_mff_count"] > 0) & (d["source_et_count"] > 0)).sum()),
            "source_EEG_pct": pct(int((d["source_mff_count"] > 0).sum()), len(d)),
            "source_ET_pct": pct(int((d["source_et_count"] > 0).sum()), len(d)),
            "complete_source_EEG_ET_pct": pct(int(((d["source_mff_count"] > 0) & (d["source_et_count"] > 0)).sum()), len(d)),
        })
    return pd.DataFrame(rows)


def stage_summary(project, inventory):
    p = project / "q1k_manager_inventory" / "q1k_summary_by_site_task.csv"
    if p.exists():
        return pd.read_csv(p)
    if inventory.empty:
        return pd.DataFrame()
    rows = []
    for (site, task), d in inventory.groupby(["site", "task"]):
        rows.append({
            "site": site, "task": task, "subjects": d["subject"].nunique(),
            "init_completed": int((d["init_status"] == "completed").sum()),
            "pylossless_completed": int((d["pylossless_status"] == "completed").sum()),
            "sync_loss_completed": int((d["sync_loss_status"] == "completed").sum()),
            "segment_completed": int((d["segment_status"] == "completed").sum()),
            "autoreject_completed": int((d["autoreject_status"] == "completed").sum()),
        })
    return pd.DataFrame(rows)


def impedance_metrics(project, inventory):
    rows = []
    files = list((project / "derivatives" / "init").glob("*/*/ses-01/eeg/*channels.tsv"))
    files += list((project / "derivatives" / "init").glob("*/*/ses-01/eeg/*electrodes.tsv"))
    per_recording = []
    for f in files:
        try:
            df = pd.read_csv(f, sep="\t")
        except Exception:
            continue
        imp_cols = [c for c in df.columns if "imped" in c.lower()]
        if not imp_cols:
            continue
        vals = pd.to_numeric(df[imp_cols[0]], errors="coerce").dropna()
        if len(vals):
            task = f.parts[f.parts.index("init") + 1] if "init" in f.parts else "UNKNOWN"
            sub = next((x[4:] for x in f.parts if x.startswith("sub-")), "UNKNOWN")
            site = re.match(r"([A-Z]+)", sub).group(1) if re.match(r"([A-Z]+)", sub) else "UNKNOWN"
            per_recording.append({"site": site, "task": task, "subject": sub, "mean_impedance_kohm": float(vals.mean())})
    df = pd.DataFrame(per_recording)
    if df.empty:
        return pd.DataFrame([{"metric": "EEG mean impedance <50 kOhm", "site": s, "value": np.nan, "reason": "No impedance column found in downloaded TSV metadata"} for s in SITES])
    for site in SITES:
        d = df[df.site == site]
        rows.append({"metric": "EEG mean impedance <50 kOhm", "site": site, "n": len(d), "value_pct": pct((d.mean_impedance_kohm < 50).sum(), len(d)), "reason": ""})
    return pd.DataFrame(rows)


def clean_duration_metrics(project):
    try:
        import mne
    except Exception as e:
        return pd.DataFrame([{"metric": "EEG >10 min clean data", "site": s, "value": np.nan, "reason": f"mne import failed: {e}"} for s in SITES])
    rows = []
    recs = []
    for f in (project / "derivatives" / "pylossless").glob("*/*/ses-01/eeg/*_eeg.edf"):
        try:
            raw = mne.io.read_raw_edf(f, preload=False, verbose=False)
            total = raw.n_times / raw.info["sfreq"]
            bad = sum(d for d, desc in zip(raw.annotations.duration, raw.annotations.description) if str(desc).startswith("BAD"))
            clean = max(0.0, total - bad)
            sub = next(x[4:] for x in f.parts if x.startswith("sub-"))
            site = re.match(r"([A-Z]+)", sub).group(1)
            task = f.parts[f.parts.index("pylossless") + 1]
            recs.append({"site": site, "task": task, "subject": sub, "clean_sec": clean, "total_sec": total})
        except Exception:
            continue
    df = pd.DataFrame(recs)
    if df.empty:
        return pd.DataFrame([{"metric": "EEG >10 min clean data", "site": s, "value": np.nan, "reason": "No readable pylossless EDF files found"} for s in SITES])
    for site in SITES:
        d = df[df.site == site]
        rows.append({"metric": "EEG >10 min clean data", "site": site, "n": len(d), "value_pct": pct((d.clean_sec > 600).sum(), len(d)), "reason": "estimated from pylossless EDF duration minus BAD annotations"})
    return pd.DataFrame(rows)


def et_valid_gaze_metrics(project):
    try:
        import mne
    except Exception as e:
        return pd.DataFrame([{"metric": "ET >50% valid gaze samples", "site": s, "value": np.nan, "reason": f"mne import failed: {e}"} for s in SITES])
    recs = []
    for f in (project / "derivatives" / "init").glob("*/*/ses-01/et/*_et.fif"):
        try:
            raw = mne.io.read_raw_fif(f, preload=True, verbose=False)
            gaze = [ch for ch in raw.ch_names if ch.lower() in {"xpos_left", "ypos_left", "xpos_right", "ypos_right"}]
            if not gaze:
                continue
            data = raw.get_data(picks=gaze)
            valid = np.isfinite(data).all(axis=0)
            valid &= np.any(data != 0, axis=0)
            sub = next(x[4:] for x in f.parts if x.startswith("sub-"))
            site = re.match(r"([A-Z]+)", sub).group(1)
            task = f.parts[f.parts.index("init") + 1]
            recs.append({"site": site, "task": task, "subject": sub, "valid_gaze_fraction": float(valid.mean())})
        except Exception:
            continue
    df = pd.DataFrame(recs)
    if df.empty:
        return pd.DataFrame([{"metric": "ET >50% valid gaze samples", "site": s, "value": np.nan, "reason": "No readable init ET fif gaze channels found"} for s in SITES])
    rows = []
    for site in SITES:
        d = df[df.site == site]
        rows.append({"metric": "ET >50% valid gaze samples", "site": site, "n": len(d), "value_pct": pct((d.valid_gaze_fraction > 0.50).sum(), len(d)), "reason": "computed from finite nonzero x/y gaze channels"})
    return pd.DataFrame(rows)


def et_calibration_metrics(project):
    recs = []
    pats = [
        re.compile(r"average\\s+error\\s*[:=]?\\s*([0-9.]+)", re.I),
        re.compile(r"avg\\.?\\s+error\\s*[:=]?\\s*([0-9.]+)", re.I),
        re.compile(r"validation.*?([0-9.]+)\\s*deg", re.I),
    ]
    for f in (project / "source_prime").glob("*/*/et/*/*.asc"):
        try:
            txt = f.read_text(errors="ignore")
        except Exception:
            continue
        vals = []
        for pat in pats:
            vals += [float(x) for x in pat.findall(txt)]
        if vals:
            site = f.parts[f.parts.index("source_prime") + 1]
            recs.append({"site": site, "subject": f.parts[f.parts.index(site) + 1], "calibration_error_deg": min(vals)})
    df = pd.DataFrame(recs)
    if df.empty:
        return pd.DataFrame([{"metric": "ET calibration error <0.5 deg", "site": s, "value": np.nan, "reason": "No parseable calibration error found in ASC files"} for s in SITES])
    rows = []
    for site in SITES:
        d = df[df.site == site]
        rows.append({"metric": "ET calibration error <0.5 deg", "site": site, "n": len(d), "value_pct": pct((d.calibration_error_deg < 0.5).sum(), len(d)), "reason": "parsed from ASC text when available"})
    return pd.DataFrame(rows)


def go_srt_metrics(project):
    rows = []
    subj_rows = []
    labels = {
        "gap": (("dtgc_d", "dtgc"), ("gcgc", "gcgc_d")),
        "overlap": (("dtoc_d", "dtoc"), ("gcoc", "gcoc_d")),
    }
    for f in (project / "derivatives" / "sync_loss" / "GO").glob("*/ses-01/eeg/*_events.tsv"):
        try:
            df = pd.read_csv(f, sep="\t")
        except Exception:
            continue
        if "trial_type" not in df or "onset" not in df:
            continue
        sub = next(x[4:] for x in f.parts if x.startswith("sub-"))
        site = re.match(r"([A-Z]+)", sub).group(1)
        ons = pd.to_numeric(df["onset"], errors="coerce").values
        typ = df["trial_type"].astype(str).values
        out = {"site": site, "subject": sub}
        for cond, (target_labs, gaze_labs) in labels.items():
            srts = []
            target_idx = np.where(np.isin(typ, target_labs))[0]
            gaze_idx = np.where(np.isin(typ, gaze_labs))[0]
            gaze_onsets = ons[gaze_idx]
            for i in target_idx:
                dt = gaze_onsets - ons[i]
                dt = dt[(dt > 0.05) & (dt < 1.5)]
                if len(dt):
                    srts.append(float(dt.min() * 1000))
            out[f"{cond}_srt_ms"] = np.nanmean(srts) if srts else np.nan
            out[f"{cond}_n_trials"] = len(srts)
        subj_rows.append(out)
    sdf = pd.DataFrame(subj_rows)
    if sdf.empty:
        return pd.DataFrame(), pd.DataFrame([{"metric": "GO SRT", "value": np.nan, "reason": "No sync_loss GO events.tsv with dt/gc labels found"}])
    gap = sdf["gap_srt_ms"].dropna()
    overlap = sdf["overlap_srt_ms"].dropna()
    paired = sdf.dropna(subset=["gap_srt_ms", "overlap_srt_ms"])
    rows.append({
        "metric": "GO gap SRT ms", "mean": gap.mean(), "sd": gap.std(ddof=1), "n_subjects": len(gap),
        "reason": "heuristic: next gaze-condition event after target-condition event within 50-1500 ms",
    })
    rows.append({"metric": "GO overlap SRT ms", "mean": overlap.mean(), "sd": overlap.std(ddof=1), "n_subjects": len(overlap), "reason": ""})
    rows.append({
        "metric": "GO overlap-minus-gap effect ms",
        "mean": (paired["overlap_srt_ms"] - paired["gap_srt_ms"]).mean(),
        "sd": (paired["overlap_srt_ms"] - paired["gap_srt_ms"]).std(ddof=1),
        "n_subjects": len(paired),
        "cohen_d_paired": cohen_d_paired(paired["gap_srt_ms"], paired["overlap_srt_ms"]),
        "reason": "positive means overlap slower than gap",
    })
    return sdf, pd.DataFrame(rows)


def vep_metrics(project):
    try:
        import mne
    except Exception as e:
        return pd.DataFrame(), pd.DataFrame([{"metric": "VEP", "value": np.nan, "reason": f"mne import failed: {e}"}])
    subj = []
    for f in (project / "derivatives" / "autoreject" / "epoch_fif_files" / "VEP").glob("*_epo.fif"):
        try:
            ep = mne.read_epochs(f, preload=True, verbose=False)
            ch = "Oz" if "Oz" in ep.ch_names else ("E75" if "E75" in ep.ch_names else None)
            if ch is None:
                continue
            conds = [c for c in ("sv06_d", "sv15_d", "sv06", "sv15") if c in ep.event_id]
            if not conds:
                continue
            data = np.concatenate([ep[c].get_data(picks=[ch]) for c in conds], axis=0).squeeze()
            ev = data.mean(axis=0)
            times = ep.times
            win = (times >= 0.08) & (times <= 0.13)
            if not win.any():
                continue
            idx = np.where(win)[0][np.argmax(ev[win])]
            amp_uv = ev[idx] * 1e6
            lat_ms = times[idx] * 1000
            snr_db = np.nan
            if "sv15_d" in ep.event_id or "sv15" in ep.event_id:
                c = "sv15_d" if "sv15_d" in ep.event_id else "sv15"
                x = ep[c].copy().crop(tmin=0, tmax=min(1.2, ep.times[-1])).get_data(picks=[ch]).squeeze()
                if x.ndim == 1:
                    x = x[None, :]
                sf = ep.info["sfreq"]
                freqs = np.fft.rfftfreq(x.shape[1], 1 / sf)
                spec = np.abs(np.fft.rfft(x, axis=1)).mean(axis=0)
                target = np.argmin(np.abs(freqs - 15))
                noise_mask = ((freqs >= 12) & (freqs <= 18) & (np.abs(freqs - 15) > 1))
                if noise_mask.any() and spec[noise_mask].mean() > 0:
                    snr_db = 20 * np.log10(spec[target] / spec[noise_mask].mean())
            sub = f.name.split("_")[0][4:]
            site = re.match(r"([A-Z]+)", sub).group(1)
            subj.append({"site": site, "subject": sub, "channel": ch, "p100_latency_ms": lat_ms, "p100_amplitude_uv": amp_uv, "snr_15hz_db": snr_db})
        except Exception:
            continue
    sdf = pd.DataFrame(subj)
    if sdf.empty:
        return sdf, pd.DataFrame([{"metric": "VEP", "value": np.nan, "reason": "No readable VEP autoreject epochs with Oz/E75"}])
    rows = []
    for col, name in [("p100_latency_ms", "VEP P100 latency ms"), ("p100_amplitude_uv", "VEP P100 amplitude uV"), ("snr_15hz_db", "VEP 15Hz SNR dB")]:
        vals = sdf[col].dropna()
        rows.append({"metric": name, "mean": vals.mean(), "sd": vals.std(ddof=1), "n_subjects": len(vals), "reason": "computed from autoreject VEP epochs; Oz preferred, E75 fallback"})
    return sdf, pd.DataFrame(rows)


def plr_metrics(project):
    try:
        import mne
        from scipy.optimize import curve_fit
    except Exception as e:
        return pd.DataFrame(), pd.DataFrame([{"metric": "PLR", "value": np.nan, "reason": f"mne/scipy import failed: {e}"}])
    subj = []
    for f in (project / "derivatives" / "autoreject" / "epoch_fif_files" / "PLR").glob("*_epo.fif"):
        try:
            ep = mne.read_epochs(f, preload=True, verbose=False)
            if "plro_d" not in ep.event_id:
                continue
            ep = ep["plro_d"]
            ch = "pupil_left" if "pupil_left" in ep.ch_names else ("pupil_right" if "pupil_right" in ep.ch_names else None)
            if ch is None:
                continue
            y = ep.get_data(picks=[ch]).squeeze()
            if y.ndim == 1:
                y = y[None, :]
            sig = np.nanmean(y, axis=0) * 1e3
            times = ep.times
            base_mask = (times >= -0.5) & (times <= 0)
            constr_mask = (times >= 0.1) & (times <= 1.5)
            if not base_mask.any() or not constr_mask.any():
                continue
            baseline = float(np.nanmean(sig[base_mask]))
            min_i = np.where(constr_mask)[0][np.nanargmin(sig[constr_mask])]
            min_val = float(sig[min_i])
            amp_pct = (baseline - min_val) / baseline * 100 if baseline else np.nan
            threshold = baseline - 0.1 * (baseline - min_val)
            post = np.where((times >= 0) & (sig <= threshold))[0]
            latency_ms = float(times[post[0]] * 1000) if len(post) else np.nan
            tau = np.nan
            fit_mask = (times >= times[min_i]) & (times <= min(2.0, times[-1]))
            if fit_mask.sum() > 10:
                tx = times[fit_mask] - times[min_i]
                yy = sig[fit_mask]
                def model(t, asym, amp, tau_):
                    return asym - amp * np.exp(-t / tau_)
                try:
                    popt, _ = curve_fit(model, tx, yy, p0=[baseline, baseline - min_val, 0.8], bounds=([-np.inf, 0, 0.05], [np.inf, np.inf, 10]), maxfev=5000)
                    tau = float(popt[2])
                except Exception:
                    pass
            sub = f.name.split("_")[0][4:]
            site = re.match(r"([A-Z]+)", sub).group(1)
            subj.append({"site": site, "subject": sub, "pupil_channel": ch, "baseline": baseline, "latency_ms": latency_ms, "constriction_amplitude_pct": amp_pct, "redilation_tau_s": tau})
        except Exception:
            continue
    sdf = pd.DataFrame(subj)
    if sdf.empty:
        return sdf, pd.DataFrame([{"metric": "PLR", "value": np.nan, "reason": "No readable PLR autoreject epochs with plro_d pupil channel"}])
    rows = []
    for col, name in [("latency_ms", "PLR latency ms"), ("constriction_amplitude_pct", "PLR constriction amplitude pct baseline"), ("redilation_tau_s", "PLR redilation tau s")]:
        vals = sdf[col].dropna()
        rows.append({"metric": name, "mean": vals.mean(), "sd": vals.std(ddof=1), "n_subjects": len(vals), "reason": "computed from plro_d autoreject epochs"})
    return sdf, pd.DataFrame(rows)


def write_latex(out, demo, invsrc, stage, quality, task_metrics):
    def fmt(x, digits=1):
        if pd.isna(x):
            return "NA"
        return f"{float(x):.{digits}f}"

    lines = []
    lines.append("% Auto-generated by paper_validation_metrics.py\n")
    if not demo.empty:
        hs = demo[demo.site == "HSJ"].iloc[0]
        mh = demo[demo.site == "MHC"].iloc[0]
        tt = demo[demo.site == "TOTAL"].iloc[0]
        lines.append("\\paragraph{Demographics.} ")
        lines.append(
            f"The demographic file included {int(tt.N_demographics)} HSJ/MHC participants "
            f"(HSJ N={int(hs.N_demographics)}, MHC N={int(mh.N_demographics)}). "
            f"Mean age was {fmt(hs.age_mean)} $\\pm$ {fmt(hs.age_sd)} years at HSJ and "
            f"{fmt(mh.age_mean)} $\\pm$ {fmt(mh.age_sd)} years at MHC. "
            f"Female representation was {fmt(hs.female_pct)}\\% at HSJ and {fmt(mh.female_pct)}\\% at MHC. "
            f"Unique source EEG+ET pairs were available for {int(hs.unique_subjects_with_source_EEG_ET_pair)} HSJ "
            f"and {int(mh.unique_subjects_with_source_EEG_ET_pair)} MHC participants.\n\n"
        )
    if not invsrc.empty:
        tt = invsrc[invsrc.site == "TOTAL"].iloc[0]
        lines.append("\\paragraph{Source inventory.} ")
        lines.append(
            f"Across HSJ and MHC, the inventory contained {int(tt.expected_site_task_subject_rows)} expected "
            f"site-task-subject recordings, including {int(tt.source_EEG_rows)} source EEG recordings, "
            f"{int(tt.source_ET_rows)} source ET recordings, and {int(tt.complete_source_EEG_ET_rows)} complete "
            f"source EEG+ET recording pairs.\n\n"
        )
    lines.append("\\paragraph{Unavailable metrics.} Metrics reported as NA in the CSV outputs were not present in the downloaded metadata and require raw QC exports or task-specific behavioral result files. No unavailable value was imputed.\n")
    (out / "validation_latex_snippet.tex").write_text("".join(lines))


def main():
    args = parse_args()
    project = Path(args.project_path)
    out = project / args.outdir
    out.mkdir(parents=True, exist_ok=True)

    participants, participants_path = load_demographics(project)
    inventory, inventory_path = load_inventory(project)

    demo = demographics_table(participants, inventory)
    invsrc = source_inventory_table(inventory)
    stage = stage_summary(project, inventory)
    quality_parts = [impedance_metrics(project), et_calibration_metrics(project)]
    if not args.skip_heavy:
        quality_parts += [clean_duration_metrics(project), et_valid_gaze_metrics(project)]
    quality = pd.concat(quality_parts, ignore_index=True)

    go_subject, go_metrics = go_srt_metrics(project)
    if args.skip_heavy:
        vep_subject = pd.DataFrame()
        vep_metrics_df = pd.DataFrame([{"metric": "VEP", "value": np.nan, "reason": "skipped by --skip-heavy"}])
        plr_subject = pd.DataFrame()
        plr_metrics_df = pd.DataFrame([{"metric": "PLR", "value": np.nan, "reason": "skipped by --skip-heavy"}])
    else:
        vep_subject, vep_metrics_df = vep_metrics(project)
        plr_subject, plr_metrics_df = plr_metrics(project)
    task_metrics = pd.concat([go_metrics, vep_metrics_df, plr_metrics_df], ignore_index=True)

    outputs = {
        "validation_demographics.csv": demo,
        "validation_source_inventory.csv": invsrc,
        "validation_stage_summary.csv": stage,
        "validation_quality_metrics.csv": quality,
        "validation_task_metrics.csv": task_metrics,
        "validation_go_subject_srt.csv": go_subject,
        "validation_vep_subject_metrics.csv": vep_subject,
        "validation_plr_subject_metrics.csv": plr_subject,
    }
    for name, df in outputs.items():
        df.to_csv(out / name, index=False)

    write_latex(out, demo, invsrc, stage, quality, task_metrics)

    print(f"Wrote validation outputs to: {out}")
    for name, df in outputs.items():
        print(f"{name}: {df.shape[0]} rows")
    print("\nKey files:")
    print(out / "validation_demographics.csv")
    print(out / "validation_quality_metrics.csv")
    print(out / "validation_task_metrics.csv")
    print(out / "validation_latex_snippet.tex")


if __name__ == "__main__":
    main()
