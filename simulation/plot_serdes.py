#!/usr/bin/env python3
"""
plot_serdes.py — visualize serdes_sim.c output.

Usage:
    python3 plot_serdes.py [dir] [--osf N] [--show]

Looks in <dir> (default '.') for any of:
    eye_ffe.csv, eye_dfe.csv          eye-diagram sample traces (OSF rows)
    rx_ffe_dfe_train_trace.csv        RX FFE/DFE adaptation trace (has header)
    ctle_sweep_trace.csv              CTLE sweep trace (has header)
    tx_ffe_train_trace.csv            TX FFE adaptation trace (has header)
    channel_taps.txt                  one FIR tap per line

Saves a PNG next to each input it finds; pass --show to also open windows.
"""
import os
import sys
import argparse
import numpy as np

import matplotlib
# default to a headless backend; --show switches to interactive
matplotlib.use("Agg")
import matplotlib.pyplot as plt   # noqa: E402

OSF_DEFAULT = 16


def _save(fig, path, show):
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    print("  wrote", path)
    if not show:
        plt.close(fig)


# ── eye diagram ──────────────────────────────────────────────────────────
def plot_eye(csv, osf, show):
    """CSV is an OSF-row matrix; column-major flatten rebuilds the time series."""
    m = np.loadtxt(csv, delimiter=",")
    if m.ndim == 1:
        m = m.reshape(osf, -1)
    trace = m.ravel(order="F")                  # rows fastest = sample order
    trace = trace[np.isfinite(trace)]

    win = 2 * osf                               # show 2 unit intervals
    n = (len(trace) - win) // osf
    if n <= 0:
        print(f"  [skip] {csv}: not enough samples")
        return
    x = np.arange(win) / osf                    # x axis in UI

    fig, ax = plt.subplots(figsize=(7, 4.5))
    for k in range(n):
        seg = trace[k * osf : k * osf + win]
        ax.plot(x, seg, color="#1f77b4", alpha=0.05, linewidth=0.7)
    ax.axvline(0.5, color="0.6", ls="--", lw=0.8)
    ax.axvline(1.5, color="0.6", ls="--", lw=0.8)
    ax.axhline(0.0, color="0.6", ls="--", lw=0.8)
    ax.set_xlabel("time  [UI]")
    ax.set_ylabel("amplitude")
    ax.set_title(f"Eye diagram — {os.path.basename(csv)}  ({n} symbols)")
    ax.set_xlim(0, 2)
    _save(fig, os.path.splitext(csv)[0] + ".png", show)


# ── RX FFE / DFE adaptation ──────────────────────────────────────────────
def plot_rx_train(csv, show):
    d = np.genfromtxt(csv, delimiter=",", names=True)
    names = d.dtype.names
    ffe = [c for c in names if c.startswith("RX_FFE")]
    dfe = [c for c in names if c.startswith("DFE")]
    pt = d["pt"]

    fig, axs = plt.subplots(2, 1, figsize=(7.5, 7), sharex=True)
    for c in ffe:
        axs[0].plot(pt, d[c], lw=0.9, label=c)
    axs[0].set_title("RX FFE tap adaptation")
    axs[0].set_ylabel("tap value")
    if len(ffe) <= 16:
        axs[0].legend(fontsize=6, ncol=2, loc="best")

    # bit_error magnitude (sampled instants only -> drop NaNs)
    be = d["bit_error"]
    m = np.isfinite(be)
    axs[1].plot(pt[m], np.abs(be[m]), lw=0.7, color="#d62728")
    axs[1].set_title("|error| at decision instants")
    axs[1].set_xlabel("sample index (pt)")
    axs[1].set_ylabel("|bit_error|")
    if dfe:
        ax2 = axs[1].twinx()
        for c in dfe:
            ax2.plot(pt, d[c], lw=0.8, ls=":")
        ax2.set_ylabel("DFE taps (dotted)")
    _save(fig, os.path.splitext(csv)[0] + ".png", show)


# ── CTLE sweep ───────────────────────────────────────────────────────────
def plot_ctle(csv, show):
    d = np.genfromtxt(csv, delimiter=",", names=True)
    pt = d["pt"]
    fig, axs = plt.subplots(2, 1, figsize=(7.5, 7), sharex=True)
    axs[0].plot(pt, d["ctle_A"], lw=1.0, label="ctle_A")
    ax0b = axs[0].twinx()
    ax0b.plot(pt, d["ctle_z"] / 1e9, lw=1.0, color="#ff7f0e", label="ctle_z [GHz]")
    axs[0].set_ylabel("ctle_A")
    ax0b.set_ylabel("ctle_z [GHz]")
    axs[0].set_title("CTLE sweep — gain / zero selection")

    be = d["bit_error"]
    m = np.isfinite(be)
    axs[1].plot(pt[m], be[m] ** 2, lw=0.6, color="#2ca02c")
    axs[1].set_title("squared error during sweep")
    axs[1].set_xlabel("sample index (pt)")
    axs[1].set_ylabel("error^2")
    _save(fig, os.path.splitext(csv)[0] + ".png", show)


# ── channel taps ─────────────────────────────────────────────────────────
def plot_taps(txt, osf, show):
    h = np.loadtxt(txt)
    fig, axs = plt.subplots(2, 1, figsize=(7.5, 7))
    axs[0].stem(np.arange(len(h)), h, basefmt=" ", markerfmt=" ")
    axs[0].set_title(f"Channel impulse response — {os.path.basename(txt)} "
                     f"({len(h)} taps @ Fs)")
    axs[0].set_xlabel("tap index (Fs samples)")
    axs[0].set_ylabel("h[n]")

    # symbol-spaced pulse response: single symbol held OSF samples
    pr = np.convolve(h, np.ones(osf))
    # best sampling phase = phase whose peak symbol-spaced cursor is largest
    phases = [np.max(np.abs(pr[ph::osf])) for ph in range(osf)]
    ph = int(np.argmax(phases))
    s = pr[ph::osf]
    k = int(np.argmax(np.abs(s)))
    idx = np.arange(len(s)) - k
    axs[1].stem(idx, s, basefmt=" ")
    axs[1].set_xlim(-5, 20)
    axs[1].axvline(0, color="0.6", ls="--", lw=0.8)
    axs[1].set_title(f"Symbol-spaced pulse response  (phase {ph}/{osf}, "
                     f"cursor={s[k]:.3f})")
    axs[1].set_xlabel("UI relative to cursor")
    axs[1].set_ylabel("amplitude")
    _save(fig, os.path.splitext(txt)[0] + "_pulse.png", show)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dir", nargs="?", default=".")
    ap.add_argument("--osf", type=int, default=OSF_DEFAULT)
    ap.add_argument("--show", action="store_true")
    args = ap.parse_args()
    if args.show:
        matplotlib.use("TkAgg", force=True)

    D = args.dir
    jobs = [
        ("eye_ffe.csv",                lambda p: plot_eye(p, args.osf, args.show)),
        ("eye_dfe.csv",                lambda p: plot_eye(p, args.osf, args.show)),
        ("rx_ffe_dfe_train_trace.csv", lambda p: plot_rx_train(p, args.show)),
        ("ctle_sweep_trace.csv",       lambda p: plot_ctle(p, args.show)),
        ("tx_ffe_train_trace.csv",     lambda p: plot_rx_train(p, args.show)),
        ("channel_taps.txt",           lambda p: plot_taps(p, args.osf, args.show)),
    ]
    found = False
    for fname, fn in jobs:
        path = os.path.join(D, fname)
        if os.path.exists(path):
            found = True
            print("plotting", fname)
            try:
                fn(path)
            except Exception as e:                       # noqa: BLE001
                print(f"  [error] {fname}: {e}")
    if not found:
        print(f"No serdes_sim output files found in '{D}'.")
        sys.exit(1)
    if args.show:
        plt.show()


if __name__ == "__main__":
    main()
