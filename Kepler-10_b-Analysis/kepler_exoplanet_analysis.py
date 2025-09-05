import lightkurve as lk
from astroquery.ipac.nexsci.nasa_exoplanet_archive import NasaExoplanetArchive
from lightkurve.periodogram import BoxLeastSquaresPeriodogram
from lightkurve import LightCurveCollection
import matplotlib.pyplot as plt
from astropy import units as u
from lightkurve import search_lightcurve
import numpy as np
from pathlib import Path
import shutil
import glob
import pandas as pd
import os
import gc
from astropy.io import fits
import math
import argparse

fits.Conf.use_memmap = False

# Argument parser
parser = argparse.ArgumentParser(description="Exoplanet processing script")
parser.add_argument("--planetname", type=str, required=True, help="Name of the planet (e.g. 'Kepler-10 b')")
args = parser.parse_args()

# Use the argument instead of input()
# Quotes " " are needed if the planet name contains spaces.
# python kepler-exoplanet-analysis_EOA_v1.py --planetname "Kepler-10 b"
planet = args.planetname

def sanitize_name(name: str) -> str:
    import re
    s = re.sub(r"[^\w\-_\. ]", "_", name).strip()
    return s.replace(" ", "_")

planet_sanitized = sanitize_name(planet)
outdir = Path(f"./{planet_sanitized}")
outdir.mkdir(parents=True, exist_ok=True)
print(f"Çıktılar {outdir} içine kaydedilecek.")

# NASA Exoplanet Archive'den gezegen bilgisi çekme
tbl = NasaExoplanetArchive.query_criteria(
    table="pscomppars",
    select="pl_name,hostname,pl_orbper,st_teff,st_rad,st_mass,sy_dist,sy_vmag,sy_gaiamag",
    where=f"pl_name='{planet}'"
)
if len(tbl) == 0:
    raise RuntimeError(f"{planet} bulunamadı. İsim formatını kontrol et!")

P_catalog = tbl["pl_orbper"][0].to_value("d")  # birimi 'd' (gün) olarak çıkar
host = str(tbl["hostname"][0])
print(f"{planet} için katalog periyodu: {P_catalog} gün, host: {host}")

tpf = lk.search_targetpixelfile(planet, author="Kepler", cadence="long").download_all(download_dir=str(outdir))
print("TPF kaydedildi:")
print("İndirilen TPF sayısı:", len(tpf))

N = len(tpf)        # number of elements
ncols = 5           # fixed number of columns
nrows = math.ceil(N / ncols)  # number of rows needed

fig, axes = plt.subplots(nrows, ncols, figsize=(3*ncols, 3*nrows))
axes = axes.flatten()

for i, ax in enumerate(axes):
    if i < N:
        cadence = tpf[i]
        cadence.plot(aperture_mask=cadence.pipeline_mask, ax=ax, show_colorbar=False)
        ax.set_title(f"Cadence {i}", fontsize=8)
    else:
        ax.axis("off")  # hide unused plots

plt.tight_layout()
plot_path = os.path.join(outdir, "tpf_grid.pdf")
plt.savefig(plot_path, dpi=300, bbox_inches="tight")

lc_collection = []

for i, t in enumerate(tpf):
    lc = t.to_lightcurve(aperture_mask=t.pipeline_mask).flatten(window_length=101).remove_nans().remove_outliers()
    globals()[f"lc_{i}"] = lc
    lc_collection.append(lc)
    print(f"lc_{i} oluşturuldu ve lc_collection'a eklendi.")
    
lc_collection = LightCurveCollection(lc_collection)
lc_stitched   = lc_collection.stitch()

lc_stitched.plot()
plot_path = os.path.join(outdir, "stitched_lightcurve.png")
plt.savefig(plot_path, dpi=300)

fig, ax = plt.subplots(figsize=(20,5))
for lc in lc_collection:
  lc.plot(ax=ax, label=f'Quarter {lc.quarter}');
  
plot_path = os.path.join(outdir, "collection_plot.png")
plt.savefig(plot_path, dpi=300)
  
  
min_period, max_period = 0.5, ((lc_stitched.time[-1].value - lc_stitched.time[0].value) / 3)
print(min_period, max_period)

periods = [] 
for i, lc in enumerate(lc_collection):
    try:
        lc_clean = lc.remove_nans().remove_outliers()
        bls = lc_clean.to_periodogram(method="bls",minimum_period=min_period, maximum_period=max_period)
        bls_period = bls.period_at_max_power.value
        periods.append(bls_period)
        print(f"lc_{i} için bulunan periyot: {bls_period:.5f} d")

    except Exception as e:
        print(f"lc_{i} için hata oluştu: {e}")
        
# ortalama periyot
if periods:
    expected = P_catalog   
    tol = 0.1      

    # filtreleme
    filtered_periods = [p for p in periods if abs(p - expected) < tol]

    if filtered_periods:
        avg_period = np.mean(filtered_periods)
        print("\nBulunan periyotlar:", [f"{p:.5f}" for p in periods])
        print("Filtrelenmiş periyotlar:", [f"{p:.5f}" for p in filtered_periods])
        print(f"Ortalama periyot (filtreli): {avg_period:.5f} d")
    else:
        print("Filtreye uyan periyot bulunamadı.")
else:
    print("Hiç periyot bulunamadı.")
    
bls = lc_stitched.to_periodogram(method="bls", minimum_period=min_period, maximum_period=max_period, frequency_factor=10000)
bls.plot()
plt.semilogx()
bls_period = bls.period_at_max_power.value
print(f"BLS ile bulunan periyot: {bls_period:.5f} d")
plt.axvline(bls_period, color='r', linestyle='dotted', label=f"Period = {bls_period:.4f} d", alpha=0.6)
plt.legend()
plot_path = os.path.join(outdir, "bls_period.png")
plt.savefig(plot_path, dpi=300)

folded_lc = lc_stitched.fold(period=bls_period).bin(time_bin_size=0.001)
folded_lc.plot()
plot_path = os.path.join(outdir, "folded_lightcurve.png")
plt.savefig(plot_path, dpi=300)
write_path = os.path.join(outdir, "binned_lightcurve.csv")
folded_lc.to_table().write(write_path, format='csv', overwrite=True)

print(f"Görseller ve CSV '{outdir}' klasörüne kaydedildi.")

# Özet CSV (periyotlar ve yıldız bilgileri)
summary_data = {
    "pl_name": [planet],
    "hostname": [host],
    "P_catalog_days": [P_catalog],
    "P_bls_days": [bls_period],
    "st_teff_K": [tbl["st_teff"][0] if "st_teff" in tbl.colnames else None],
    "st_rad_Rsun": [tbl["st_rad"][0] if "st_rad" in tbl.colnames else None],
    "st_mass_Msun": [tbl["st_mass"][0] if "st_mass" in tbl.colnames else None],
    "sy_dist_pc": [tbl["sy_dist"][0] if "sy_dist" in tbl.colnames else None],
    "sy_vmag": [tbl["sy_vmag"][0] if "sy_vmag" in tbl.colnames else None],
    "sy_gaiamag": [tbl["sy_gaiamag"][0] if "sy_gaiamag" in tbl.colnames else None],
}

summary_df = pd.DataFrame(summary_data)
summary_path = os.path.join(outdir, "planet_summary.csv")
summary_df.to_csv(summary_path, index=False)

print(f"Özet CSV kaydedildi: {summary_path}")


#mastdowload silmek için
def cleanup_mast(outdir: Path):
    mast_path = outdir / "mastDownload"
    if mast_path.exists():
        try:
            # RAM'deki objeleri serbest bırak
            gc.collect()
            shutil.rmtree(mast_path)
            print(f"{mast_path} klasörü silindi (ham MAST indirmeleri temizlendi).")
        except Exception as e:
            print(f"{mast_path} silinirken hata oluştu: {e}")
    else:
        print("mastDownload klasörü bulunamadı")

cleanup_mast(outdir)

