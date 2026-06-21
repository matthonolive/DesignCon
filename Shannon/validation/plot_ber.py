#!/usr/bin/env python3
"""Plot SER/BER vs SNR from ber_sweep.csv. Usage: plot_ber.py in.csv out.png"""
import sys, csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

inp = sys.argv[1] if len(sys.argv) > 1 else "ber_sweep.csv"
out = sys.argv[2] if len(sys.argv) > 2 else "ber_sweep.png"
snr=[]; mf=[]; mx=[]; dfe=[]; slc=[]
with open(inp) as f:
    for r in csv.DictReader(f):
        snr.append(float(r["snr_db"]))
        mf.append(float(r["ser_mlse_float"])); mx.append(float(r["ser_mlse_fixed"]))
        dfe.append(float(r["ser_dfe"]));        slc.append(float(r["ser_slicer"]))
def clip(v): return [max(x,1e-6) for x in v]
plt.figure(figsize=(7,5))
plt.semilogy(snr, clip(slc), 's--', label="memoryless slicer")
plt.semilogy(snr, clip(dfe), '^-',  label="FFE+DFE (target cancel)")
plt.semilogy(snr, clip(mf),  'o-',  label="FFE+MLSE (float)")
plt.semilogy(snr, clip(mx),  'x:',  label="FFE+MLSE (fixed, firmware)")
plt.gca().invert_xaxis()
plt.xlabel("RX SNR (dB)"); plt.ylabel("Symbol error rate")
plt.title("PAM-4 FFE→PR-target→MLSE vs DFE / slicer")
plt.grid(True, which="both", alpha=0.3); plt.legend()
plt.tight_layout(); plt.savefig(out, dpi=120)
print("saved", out)