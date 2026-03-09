#!/usr/bin/env python3
"""
Find and print ONLY the histograms booked by ScoutingHadronicTopMassAnalyzer
inside a DQMIO ROOT file.

Builds the exact expected histogram name list from the C++ bookHistograms()
logic, then matches them against what is actually in the file.

Usage:
    python DQM_Hist_finder.py ScoutingDQM_DQM.root
    python DQM_Hist_finder.py ScoutingDQM_DQM.root --plot
    python DQM_Hist_finder.py ScoutingDQM_DQM.root --filter resonanceZ

Requires: uproot, numpy, awkward  (pip install uproot numpy awkward)
Optional: matplotlib              (pip install matplotlib)  -- for --plot
"""

import sys
import argparse
import numpy as np

try:
    import uproot
except ImportError:
    print("ERROR: 'uproot' is required.  Install with:  pip install uproot awkward numpy")
    sys.exit(1)

# ---------------------------------------------------------------------------
#  Build the exact histogram names that ScoutingHadronicTopMassAnalyzer books
# ---------------------------------------------------------------------------

RESONANCES = ["resonanceZ", "resonanceJ", "resonanceY", "resonanceAll"]

# Per-resonance kinematic histogram suffixes (from bookHistograms_resonance)
KINEMATIC_SUFFIXES = [
    "Probe_sctElectron_Pt_Barrel",
    "Probe_sctElectron_Pt_Endcap",
    "Probe_sctElectron_Eta",
    "Probe_sctElectron_Phi",
    "Probe_sctElectron_HoverE_Barrel",
    "Probe_sctElectron_HoverE_Endcap",
    "Probe_sctElectron_OoEMOoP_Barrel",
    "Probe_sctElectron_OoEMOoP_Endcap",
    "Probe_sctElectron_dPhiIn_Barrel",
    "Probe_sctElectron_dPhiIn_Endcap",
    "Probe_sctElectron_dEtaIn_Barrel",
    "Probe_sctElectron_dEtaIn_Endcap",
    "Probe_sctElectron_SigmaIetaIeta_Barrel",
    "Probe_sctElectron_SigmaIetaIeta_Endcap",
    "Probe_sctElectron_MissingHits_Barrel",
    "Probe_sctElectron_MissingHits_Endcap",
    "Probe_sctElectron_Trackfbrem_Barrel",
    "Probe_sctElectron_Trackfbrem_Endcap",
    "Probe_sctElectron_Track_pt_Barrel",
    "Probe_sctElectron_Track_pt_Endcap",
    "Probe_sctElectron_Track_pMode_Barrel",
    "Probe_sctElectron_Track_pMode_Endcap",
    "Probe_sctElectron_Track_etaMode_Barrel",
    "Probe_sctElectron_Track_etaMode_Endcap",
    "Probe_sctElectron_Track_phiMode_Barrel",
    "Probe_sctElectron_Track_phiMode_Endcap",
    "Probe_sctElectron_Track_qoverpModeError_Barrel",
    "Probe_sctElectron_Track_qoverpModeError_Endcap",
    "Probe_sctElectron_RelEcalIsolation_Barrel",
    "Probe_sctElectron_RelEcalIsolation_Endcap",
    "Probe_sctElectron_RelHcalIsolation_Barrel",
    "Probe_sctElectron_RelHcalIsolation_Endcap",
    "Probe_sctElectron_RelTrackIsolation_Barrel",
    "Probe_sctElectron_RelTrackIsolation_Endcap",
    "sctElectron_Invariant_Mass",
    "Probe_sctElectron_Pt_Barrel_passID",
    "Probe_sctElectron_Pt_Endcap_passID",
    # Leading / sub-leading trigger efficiency (base DST)
    "leading_Pt_Barrel_passBaseDST",
    "leading_Pt_Endcap_passBaseDST",
    "leading_Eta_passBaseDST",
    "subleading_Pt_Barrel_passBaseDST",
    "subleading_Pt_Endcap_passBaseDST",
    "subleading_Eta_passBaseDST",
]
# NOTE: The trigger-path-specific histograms (fireTrigObj) are dynamic and
#       depend on the triggerSelection config.  We can't predict those names
#       statically, so we match them with a prefix pattern below.


def build_expected_names():
    """Return set of histogram short-names (no folder prefix) that the module books."""
    names = {"EventCount"}
    for res in RESONANCES:
        for suffix in KINEMATIC_SUFFIXES:
            names.add(f"{res}_{suffix}")
    return names


def is_module_histogram(full_path, expected_names):
    """
    Check whether a DQMIO FullName belongs to ScoutingHadronicTopMassAnalyzer.
    Matches by the histogram short name (last component of the path).
    Also catches dynamic trigger-efficiency histograms via prefix matching.
    """
    short = full_path.split("/")[-1] if "/" in full_path else full_path
    # Exact match for the known static names
    if short in expected_names:
        return True
    # Dynamic trigger-efficiency histograms: {res}_{leading|subleading}_*_fireTrigObj
    for res in RESONANCES:
        if short.startswith(f"{res}_leading_") and short.endswith("_fireTrigObj"):
            return True
        if short.startswith(f"{res}_subleading_") and short.endswith("_fireTrigObj"):
            return True
    return False


# ---------------------------------------------------------------------------
#  DQMIO reading helpers
# ---------------------------------------------------------------------------

DQMIO_TREES = {
    "TH1Fs": "TH1F", "TH1Ds": "TH1D", "TH1Ss": "TH1S", "TH1Is": "TH1I",
    "TH2Fs": "TH2F", "TH2Ds": "TH2D", "TH2Ss": "TH2S", "TH2Is": "TH2I",
    "TH3Fs": "TH3F", "TProfiles": "TProfile", "TProfile2Ds": "TProfile2D",
}


def scan_dqmio(root_file, expected_names):
    """
    Scan every DQMIO TTree and return:
      found   – list of (tree_name, index, full_path) for module histograms
      others  – total count of histograms NOT belonging to this module
    """
    found, others = [], 0
    for tree_name in DQMIO_TREES:
        if tree_name not in root_file:
            continue
        tree = root_file[tree_name]
        if "FullName" not in tree.keys():
            continue
        names = tree["FullName"].array(library="np")
        for idx, fp in enumerate(names):
            if is_module_histogram(fp, expected_names):
                found.append((tree_name, idx, fp))
            else:
                others += 1
    return found, others


def read_histogram(root_file, tree_name, idx):
    """Deserialise one histogram object from a DQMIO Value branch."""
    return root_file[tree_name]["Value"].array(
        entry_start=idx, entry_stop=idx + 1, library="np"
    )[0]


# ---------------------------------------------------------------------------
#  Display helpers
# ---------------------------------------------------------------------------

def print_hist_summary(full_path, htype, hist_obj):
    """One-line compact summary for a histogram."""
    try:
        values = hist_obj.values()
        total  = np.sum(values)
        edges  = hist_obj.axis().edges()
        nbins  = len(values)
        status = f"entries={total:<10.6g}  bins={nbins}  range=[{edges[0]:.4g},{edges[-1]:.4g}]"
    except Exception:
        status = "(could not read bin info)"
    short = full_path.split("/")[-1] if "/" in full_path else full_path
    tag   = "  FILLED" if "entries" in status and "entries=0" not in status else "   EMPTY"
    print(f"  {tag}  [{htype:>8s}]  {short:60s}  {status}")


def print_hist_detail(full_path, htype, hist_obj):
    """Verbose multi-line info for a single histogram."""
    print(f"\n{'='*90}")
    print(f"  Path : {full_path}")
    print(f"  Type : {htype}")
    try:
        values  = hist_obj.values()
        edges   = hist_obj.axis().edges()
        total   = np.sum(values)
        print(f"  Bins : {len(values)}")
        print(f"  Range: [{edges[0]:.4g}, {edges[-1]:.4g}]")
        print(f"  Sum of bin contents: {total:.6g}")
        if total > 0:
            centers = 0.5 * (edges[:-1] + edges[1:])
            print(f"  Weighted mean: {np.average(centers, weights=values):.4g}")
        nonzero = np.nonzero(values)[0]
        if len(nonzero) > 0:
            print(f"  Non-zero bins: {len(nonzero)}")
            print(f"  {'Bin':>6s}  {'Low':>12s}  {'High':>12s}  {'Content':>12s}")
            for bi in nonzero[:20]:
                print(f"  {bi:>6d}  {edges[bi]:>12.4g}  {edges[bi+1]:>12.4g}  {values[bi]:>12.4g}")
            if len(nonzero) > 20:
                print(f"  ... and {len(nonzero)-20} more non-zero bins")
        else:
            print("  ** All bins are empty **")
    except Exception as e:
        print(f"  (Error: {e})")
    print(f"{'='*90}")


def plot_histogram(full_path, hist_obj, output_dir):
    """Save a 1-D histogram as PNG."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return
    try:
        values = hist_obj.values()
        edges  = hist_obj.axis().edges()
    except Exception:
        return
    if np.sum(values) == 0:
        return
    fig, ax = plt.subplots(figsize=(8, 5))
    centers = 0.5 * (edges[:-1] + edges[1:])
    widths  = edges[1:] - edges[:-1]
    ax.bar(centers, values, width=widths, edgecolor="black", linewidth=0.3, alpha=0.7)
    short = full_path.split("/")[-1]
    ax.set_title(short, fontsize=10)
    ax.set_xlabel("Value"); ax.set_ylabel("Entries")
    import os
    os.makedirs(output_dir, exist_ok=True)
    out = os.path.join(output_dir, short.replace("/", "_") + ".png")
    fig.tight_layout(); fig.savefig(out, dpi=120); plt.close(fig)
    print(f"    -> saved {out}")


# ---------------------------------------------------------------------------
#  Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Find ONLY ScoutingHadronicTopMassAnalyzer histograms in a DQMIO file"
    )
    parser.add_argument("rootfile", nargs="?", default="ScoutingDQM_DQM.root",
                        help="DQMIO ROOT file  (default: ScoutingDQM_DQM.root)")
    parser.add_argument("--filter", "-f", default=None,
                        help="Extra substring filter on histogram name")
    parser.add_argument("--detail", "-d", action="store_true",
                        help="Print detailed per-bin info for each histogram")
    parser.add_argument("--plot", "-p", action="store_true",
                        help="Save PNGs of non-empty 1D histograms")
    parser.add_argument("--plot-dir", default="dqm_plots",
                        help="Output directory for plots  (default: dqm_plots/)")
    args = parser.parse_args()

    # ---- open file ----
    print(f"\nOpening: {args.rootfile}\n")
    try:
        f = uproot.open(args.rootfile)
    except Exception as e:
        print(f"ERROR: {e}"); sys.exit(1)

    # ---- build expected names & scan ----
    expected = build_expected_names()
    found, n_other = scan_dqmio(f, expected)

    # optional extra filter
    if args.filter:
        found = [(t, i, p) for t, i, p in found if args.filter.lower() in p.lower()]

    # ---- report ----
    print(f"  Total histograms in file       : {len(found) + n_other}")
    print(f"  From ScoutingHadronicTopMass   : {len(found)}")
    print(f"  From other modules (ignored)   : {n_other}")

    if not found:
        print("\n  ** No ScoutingHadronicTopMassAnalyzer histograms found! **")
        print("     Possible reasons:")
        print("       - Module did not run or booked 0 histograms")
        print("       - OutputInternalPath does not match expected names")
        print("\n     All DQM paths in the file:")
        for tree_name in DQMIO_TREES:
            if tree_name not in f:
                continue
            tree = f[tree_name]
            if "FullName" not in tree.keys():
                continue
            for name in tree["FullName"].array(library="np"):
                print(f"       {name}")
        sys.exit(0)

    # Determine DQM folder from the first match
    first_path = found[0][2]
    folder = "/".join(first_path.split("/")[:-1]) if "/" in first_path else "(root)"
    print(f"  DQM folder (OutputInternalPath): {folder}")

    # Separate filled vs empty
    filled, empty = [], []
    print(f"\n{'─'*90}")
    print(f"  ScoutingHadronicTopMassAnalyzer histograms")
    print(f"{'─'*90}")
    for tree_name, idx, full_path in found:
        htype = DQMIO_TREES[tree_name]
        try:
            h = read_histogram(f, tree_name, idx)
            total = np.sum(h.values())
        except Exception:
            h, total = None, 0
        print_hist_summary(full_path, htype, h)
        if total > 0:
            filled.append((tree_name, idx, full_path, htype, h))
        else:
            empty.append((tree_name, idx, full_path, htype, h))

    # Quick tally
    print(f"\n  Summary:  {len(filled)} filled,  {len(empty)} empty,  {len(found)} total")

    # Check which expected names are missing from the file entirely
    found_shorts = {p.split("/")[-1] for _, _, p in found}
    missing = sorted(expected - found_shorts)
    if missing:
        print(f"\n  Expected but MISSING from file ({len(missing)}):")
        for m in missing:
            print(f"    - {m}")

    # ---- detail / plot ----
    if args.detail:
        print(f"\n\n{'#'*90}")
        print(f"  DETAILED HISTOGRAM INFO  (only filled histograms)")
        print(f"{'#'*90}")
        for tree_name, idx, full_path, htype, h in filled:
            print_hist_detail(full_path, htype, h)

    if args.plot:
        print(f"\n  Saving plots to {args.plot_dir}/ ...")
        for tree_name, idx, full_path, htype, h in filled:
            if "TH1" in htype:
                plot_histogram(full_path, h, args.plot_dir)
        print("  Done plotting.")

    print()


if __name__ == "__main__":
    main()
