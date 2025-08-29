#!/usr/bin/env python3
import os
import re
import json
import math
import time
import warnings
import csv
import threading  # ### FIX: thread-safe log için
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np
import pandas as pd

# ### FIX: Başsız ortamlar için backend'i pyplot'tan önce ayarla
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from astropy import units as u
try:
    from astroquery.ipac.nexsci.nasa_exoplanet_archive import NasaExoplanetArchive
except Exception:
    from astroquery.nasa_exoplanet_archive import NasaExoplanetArchive
import lightkurve as lk

### ayarlar
OUTPUT_DIR = "/arf/scratch/egitim112/exoplanet_output2"
CACHE_DIR = os.path.join(OUTPUT_DIR, "lk_cache")
PNG_DIR = os.path.join(OUTPUT_DIR, "png")
CSV_DIR = os.path.join(OUTPUT_DIR, "csv")
MANIFEST = os.path.join(OUTPUT_DIR, "manifest.jsonl")
LOGFILE = os.path.join(OUTPUT_DIR, "run_log.jsonl")

# Yerel CSV yolu (senin yüklediğin dosya)
INPUT_FILE = "/arf/scratch/egitim112/transit_data.csv"

MAX_WORKERS = max(4, os.cpu_count() or 4)
START_INDEX = 40   # kaçıncı satırdan başlayacağını burada ayarlarsın
MAX_TARGETS = 20
MISSION_PRIORITY = ["TESS", "Kepler", "K2"]
AUTHOR_PRIORITY = ["SPOC", "QLP", "Kepler", "K2"]
FLATTEN_WINDOW = 301
TIME_BIN = 0.001
RETRY = 3
RETRY_BASE_SLEEP = 2.0
###

warnings.filterwarnings("ignore")
try:
    lk.log.setLevel("ERROR")
except Exception:
    pass

for d in [OUTPUT_DIR, CACHE_DIR, PNG_DIR, CSV_DIR]:
    os.makedirs(d, exist_ok=True)

# ### FIX: log yazımı için kilit
_LOG_LOCK = threading.Lock()

def sanitize(name: str) -> str:
    name = re.sub(r"[^\w\-\.]+", "_", name, flags=re.UNICODE)
    return re.sub(r"_+", "_", name).strip("_")

def save_line(path, rec: dict):
    # ### FIX: thread-safe append
    with _LOG_LOCK:
        with open(path, "a", encoding="utf-8") as f:
            f.write(json.dumps(rec, ensure_ascii=False) + "\n")

def safe_value(x, unit=None):
    if x is None:
        return None
    try:
        if pd.isna(x):
            return None
    except Exception:
        pass
    if hasattr(x, "to") and hasattr(x, "unit"):
        try:
            return x.to(unit).value if unit else x.value
        except Exception:
            try:
                return x.value
            except Exception:
                try:
                    return float(x)
                except Exception:
                    return None
    try:
        return float(x)
    except Exception:
        return None

def get_time_offset(lc) -> float:
    # ### FIX: daha güvenli biçimde format yakala
    fmt = None
    try:
        fmt = getattr(lc, "time_format", None)
    except Exception:
        pass
    if fmt is None:
        try:
            fmt = getattr(getattr(lc, "time", None), "format", None)
        except Exception:
            pass
    if fmt:
        f = str(fmt).lower()
        if f == "btjd":
            return 2457000.0
        if f == "bkjd":
            return 2454833.0
    return 0.0

def detect_delimiter(path, comment_char='#'):
    with open(path, 'r', encoding='utf-8', errors='ignore') as f:
        sample_lines = []
        for _ in range(200):
            line = f.readline()
            if not line:
                break
            if line.lstrip().startswith(comment_char) or line.strip() == "":
                continue
            sample_lines.append(line)
            if len(sample_lines) >= 20:
                break
        sample = ''.join(sample_lines)
    try:
        dialect = csv.Sniffer().sniff(sample, delimiters=[',', '\t', ';', '|', ':'])
        return dialect.delimiter
    except Exception:
        return ','

def _read_csv_robust(path, sep, header=True):
    # ### FIX: bazı pandas sürümlerinde engine='c' + comment sorun çıkarabilir → python'a düş
    try:
        return pd.read_csv(path, sep=sep, header=0 if header else None, engine='c', comment='#', low_memory=False)
    except Exception:
        return pd.read_csv(path, sep=sep, header=0 if header else None, engine='python', comment='#', low_memory=False)

def rows_from_local_file(path, name_col=0, header=True, max_targets=None):
    sep = detect_delimiter(path)
    df = _read_csv_robust(path, sep, header=header)

    needed_cols = {"pl_name", "hostname", "pl_orbper", "pl_tranmid", "pl_trandur"}
    dfcols = set([str(c) for c in df.columns])
    rows = []
    if needed_cols.issubset(dfcols):
        if max_targets:
            df = df.iloc[START_INDEX:START_INDEX + max_targets]
        for _, r in df.iterrows():
            rows.append({k: r[k] for k in ["pl_name", "hostname", "pl_orbper", "pl_tranmid", "pl_trandur"]})
        return rows
    # Eğer tam parametre yoksa isim sütunundan al ve eksik parametreleri NASA'dan sorgula
    names = df.iloc[:, name_col].astype(str).tolist()
    if max_targets:
        names = names[:max_targets]
    for name in names:
        name_clean = str(name).strip()
        if not name_clean:
            continue
        esc = name_clean.replace("'", "''")
        try:
            tbl = NasaExoplanetArchive.query_criteria(
                table="pscomppars",
                select="pl_name,hostname,pl_orbper,pl_tranmid,pl_trandur",
                where=f"pl_name = '{esc}'"
            )
            if len(tbl) > 0:
                r = tbl[0]
                rows.append({k: r[k] for k in tbl.colnames})
            else:
                rows.append({"pl_name": name_clean, "hostname": None, "pl_orbper": None, "pl_tranmid": None, "pl_trandur": None})
        except Exception as e:
            save_line(LOGFILE, {"planet": name_clean, "status": "fetch_params_error", "error": repr(e)})
            rows.append({"pl_name": name_clean, "hostname": None, "pl_orbper": None, "pl_tranmid": None, "pl_trandur": None})
    return rows

def _flatten_or_normalize(lc):
    # ### FIX: kısa seri/NaN durumları için daha yumuşak yaklaşım
    lc2 = lc.remove_nans()
    try:
        wl = int(FLATTEN_WINDOW)
        if wl % 2 == 0:
            wl += 1
        return lc2.flatten(window_length=wl)
    except Exception:
        try:
            return lc2.normalize()
        except Exception:
            return lc2

def search_download_lightcurve(hostname: str):
    # ### FIX: daha çok log ve retry üst katmanda
    for mission in MISSION_PRIORITY:
        try:
            search = lk.search_lightcurve(hostname, mission=mission)
        except Exception as e:
            save_line(LOGFILE, {"host": hostname, "mission": mission, "status": "search_error", "error": repr(e)})
            continue
        if len(search) == 0:
            continue
        for author in AUTHOR_PRIORITY:
            sub = search[search.author == author]
            if len(sub) == 0:
                continue
            try:
                lcc = sub.download_all(download_dir=CACHE_DIR)
                if not lcc or len(lcc) == 0:
                    continue
                lc = lcc.stitch()
                lc = _flatten_or_normalize(lc)
                return lc, mission, author
            except Exception as e:
                save_line(LOGFILE, {"host": hostname, "mission": mission, "author": author, "status": "download_error", "error": repr(e)})
                continue
        # yazar filtrelemeden dene
        try:
            lcc = search.download_all(download_dir=CACHE_DIR)
            if lcc and len(lcc) > 0:
                lc = lcc.stitch()
                lc = _flatten_or_normalize(lc)
                return lc, mission, "auto"
        except Exception as e:
            save_line(LOGFILE, {"host": hostname, "mission": mission, "status": "download_error_auto", "error": repr(e)})
            pass
    return None, None, None

def fold_plot_save(planet, host, P_day, t0_bjd, dur_hr, lc, mission, author):
    base = sanitize(planet)
    png_path = os.path.join(PNG_DIR, f"{base}.png")
    csv_path = os.path.join(CSV_DIR, f"{base}.csv")

    offset = get_time_offset(lc)
    epoch_time = t0_bjd - offset

    folded = lc.fold(period=P_day, epoch_time=epoch_time)
    try:
        folded_binned = folded.bin(time_bin_size=TIME_BIN)
    except Exception:
        folded_binned = folded

    if dur_hr is not None and not (isinstance(dur_hr, float) and math.isnan(dur_hr)):
        dur_days = dur_hr / 24.0
        half_win = min(0.5, 3.0 * (dur_days / P_day))
    else:
        half_win = 0.15

    plt.figure(figsize=(9, 4))
    ax = plt.gca()
    try:
        folded.plot(ax=ax, marker=".", linestyle="none", alpha=0.3, label="Folded")
    except Exception:
        pass
    try:
        folded_binned.plot(ax=ax, marker=".", linestyle="none", label="Binned")
    except Exception:
        pass
    plt.xlim(-half_win, half_win)
    ttl = f"{planet} — Transit (mission={mission}, author={author})"
    plt.title(ttl)
    plt.xlabel("Faz (gün)")
    plt.ylabel("Normalize Akı")
    plt.legend()
    plt.tight_layout()
    plt.savefig(png_path, dpi=150)
    plt.close()

    # ### FIX: flux_err güvenli çıkarım
    try:
        flux_err = getattr(folded, "flux_err", None)
        if flux_err is not None:
            try:
                flux_err = flux_err.value
            except Exception:
                flux_err = np.asarray(flux_err)
    except Exception:
        flux_err = None

    df = pd.DataFrame({
        "phase_day": np.asarray(getattr(folded.time, "value", folded.time)),
        "flux": np.asarray(getattr(folded.flux, "value", folded.flux)),
        "flux_err": flux_err
    })
    df.to_csv(csv_path, index=False)

    save_line(MANIFEST, {
        "planet": planet,
        "host": host,
        "period_day": P_day,
        "t0_bjd": t0_bjd,
        "time_offset_applied": offset,
        "mission": mission,
        "author": author,
        "png": os.path.basename(png_path),
        "csv": os.path.basename(csv_path)
    })

def already_done(planet: str) -> bool:
    base = sanitize(planet)
    return (os.path.exists(os.path.join(PNG_DIR, f"{base}.png")) and
            os.path.exists(os.path.join(CSV_DIR, f"{base}.csv")))

def process_one(row_dict):
    planet = str(row_dict.get("pl_name", "")).strip()
    host = row_dict.get("hostname") if row_dict.get("hostname") is not None else planet
    if host is None:
        host = planet
    host = str(host).strip()

    if already_done(planet):
        save_line(LOGFILE, {"planet": planet, "status": "skip_exists"})
        return planet, "skip_exists"

    P_day = safe_value(row_dict.get("pl_orbper"), u.day)
    t0_bjd = safe_value(row_dict.get("pl_tranmid"), u.day)
    dur_hr = safe_value(row_dict.get("pl_trandur"), u.hour)

    # 1) Eğer eksik parametre varsa önce NASA'dan tekrar dene (satır bazlı sorgu)
    if P_day is None or t0_bjd is None:
        try:
            esc = planet.replace("'", "''")
            tbl = NasaExoplanetArchive.query_criteria(
                table="pscomppars",
                select="pl_name,hostname,pl_orbper,pl_tranmid,pl_trandur",
                where=f"pl_name = '{esc}'"
            )
            if len(tbl) > 0:
                r = tbl[0]
                if P_day is None:
                    P_day = safe_value(r.get("pl_orbper"), u.day)
                if t0_bjd is None:
                    t0_bjd = safe_value(r.get("pl_tranmid"), u.day)
                if dur_hr is None:
                    dur_hr = safe_value(r.get("pl_trandur"), u.hour)
        except Exception as e:
            save_line(LOGFILE, {"planet": planet, "status": "fetch_params_error", "error": repr(e)})

    # 2) Eğer hâlâ pl_tranmid yok ama periyot varsa: lightkurve'den LC indirip kaba bir t0 tahmini dene
    if P_day is not None and t0_bjd is None:
        try:
            lc_try, mission_try, author_try = search_download_lightcurve(host)
            if lc_try is not None:
                try:
                    lc_proc = _flatten_or_normalize(lc_try)
                    t_arr = np.asarray(getattr(lc_proc.time, "value", lc_proc.time))
                    f_arr = np.asarray(getattr(lc_proc.flux, "value", lc_proc.flux))
                    import pandas as _pd
                    s = _pd.Series(f_arr)
                    win = max(3, int(len(s) * 0.01))
                    sm = s.rolling(window=win, center=True, min_periods=1).median().values
                    idx = int(np.argmin(sm))
                    offset = get_time_offset(lc_proc)
                    t0_bjd = float(t_arr[idx]) + float(offset)
                    save_line(LOGFILE, {"planet": planet, "status": "estimated_t0_from_lc", "t0_bjd": t0_bjd, "mission": mission_try, "author": author_try})
                except Exception as e:
                    save_line(LOGFILE, {"planet": planet, "status": "estimate_t0_failed", "error": repr(e)})
        except Exception as e:
            save_line(LOGFILE, {"planet": planet, "status": "lc_fetch_failed_for_t0", "error": repr(e)})

    if P_day is None or t0_bjd is None:
        save_line(LOGFILE, {"planet": planet, "status": "skip_missing_params"})
        return planet, "skip_missing_params"

    last_err = None
    for attempt in range(RETRY):
        try:
            lc, mission, author = search_download_lightcurve(host)
            if lc is None:
                save_line(LOGFILE, {"planet": planet, "status": "no_data"})
                return planet, "no_data"

            fold_plot_save(planet, host, P_day, t0_bjd, dur_hr, lc, mission, author)
            save_line(LOGFILE, {"planet": planet, "status": "ok", "mission": mission, "author": author})
            return planet, "ok"
        except Exception as e:
            last_err = repr(e)
            time.sleep(RETRY_BASE_SLEEP * (2 ** attempt))

    save_line(LOGFILE, {"planet": planet, "status": "error", "error": last_err})
    return planet, "error"


def fetch_table():
    # ### NOT: TRUBA compute node'unda internet yoksa bu çağrı time-out verir
    tbl = NasaExoplanetArchive.query_criteria(
        table="pscomppars",
        select="pl_name,hostname,pl_orbper,pl_tranmid,pl_trandur",
        where="pl_tranmid IS NOT NULL AND pl_orbper IS NOT NULL"
    )
    if MAX_TARGETS:
        tbl = tbl[:MAX_TARGETS]
    rows = []
    for r in tbl:
        rows.append({k: r[k] for k in tbl.colnames})
    return rows

def main():
    print(" Exoplanet Archive sorgulanıyor...")
    if INPUT_FILE and os.path.exists(INPUT_FILE):
        print(f" Local input file kullanılıyor: {INPUT_FILE}")
        rows = rows_from_local_file(INPUT_FILE, name_col=0, header=True, max_targets=MAX_TARGETS)
    else:
        rows = fetch_table()
    total = len(rows)
    print(f" Hedef sayısı: {total}")
    ok = skip = nodata = err = 0
    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as ex:
        futures = [ex.submit(process_one, r) for r in rows]
        for f in as_completed(futures):
            planet, status = f.result()
            if status == "ok":
                ok += 1
            elif status.startswith("skip"):
                skip += 1
            elif status == "no_data":
                nodata += 1
            else:
                err += 1
    summary_msg = f"\n Bitti | OK: {ok} | Skip: {skip} | No-data: {nodata} | Error: {err}"
    path_msg = (f" Çıktı klasörü: {OUTPUT_DIR}\n"
                f"- PNG: {PNG_DIR}\n- CSV: {CSV_DIR}\n- Manifest: {MANIFEST}\n- Log: {LOGFILE}")

# hem log dosyasına yaz
    save_line(LOGFILE, {"status": "summary", "ok": ok, "skip": skip,
                    "no_data": nodata, "error": err})

# hem de ekrana yazmayı dene (ama kapanmışsa sessiz geç)
    try:
        print(summary_msg, flush=True)
        print(path_msg, flush=True)
    except Exception:
        pass

if __name__ == "__main__":
    main()

