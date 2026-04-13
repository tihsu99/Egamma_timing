import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse
from sklearn.metrics import roc_curve, auc
import os, re

def load_data(file_path):
    return pd.read_parquet(file_path)

def compute_weights(sig_df, bkg_df, pt_bins, eta_bins):
    """Compute 2D pt-eta weights to reweight signal to background"""
    sig_hist, _, _ = np.histogram2d(sig_df['eg-et'], sig_df['eg-eta'], bins=[pt_bins, eta_bins])
    bkg_hist, _, _ = np.histogram2d(bkg_df['eg-et'], bkg_df['eg-eta'], bins=[pt_bins, eta_bins])
    
    weights_2d = np.divide(bkg_hist, sig_hist, out=np.zeros_like(bkg_hist), where=sig_hist>0)
   

    sig_pt_idx = np.clip(np.digitize(sig_df['eg-et'], pt_bins) - 1, 0, len(pt_bins)-2)
    sig_eta_idx = np.clip(np.digitize(sig_df['eg-eta'], eta_bins) - 1, 0, len(eta_bins)-2)

    weights = weights_2d[sig_pt_idx, sig_eta_idx]
    
    return weights

def compute_roc(sig_df, bkg_df, var, eta_cut=None):
    """Compute ROC curve and AUC."""
    sig_region = sig_df.query(eta_cut) if eta_cut else sig_df
    bkg_region = bkg_df.query(eta_cut) if eta_cut else bkg_df

    sig_vals = sig_region[var].to_numpy()
    bkg_vals = bkg_region[var].to_numpy()

    sig_w = sig_region["weight"].to_numpy()
    bkg_w = np.ones(len(bkg_vals))  # unweighted background

    y_true = np.concatenate([np.ones_like(sig_vals), np.zeros_like(bkg_vals)])
    y_score = np.concatenate([sig_vals, bkg_vals])
    weights = np.concatenate([sig_w, bkg_w])

    fpr, tpr, _ = roc_curve(y_true, -y_score, sample_weight=weights)
    roc_auc = auc(fpr, tpr)
    return fpr, tpr, roc_auc


def plot_roc_comparison(sig_df, bkg_df, iso_vars, eta_cut, default_iso, default_cut, title, outname):
    plt.figure(figsize=(8, 6))
    record = dict()

    sig_region = sig_df.query(eta_cut) if eta_cut else sig_df
    bkg_region = bkg_df.query(eta_cut) if eta_cut else bkg_df


    colors = plt.get_cmap("tab10").colors

    for color_idx, iso in enumerate(iso_vars):
        # --- Default (dt_9999.0) ---
        iso_name = iso.replace("-dt_9999.0", "").replace("eg-", "").replace("_iso-", " ").replace("_veto", "-ex")
        fpr, tpr, roc_auc = compute_roc(sig_region, bkg_region, iso, eta_cut=None)
        plt.plot(tpr, 1 - fpr, lw=2, label=f"{iso_name}", color= colors[color_idx])

        # Default cut efficiency point
        if iso == default_iso:
          sig_sel = sig_region[iso] < default_cut
          bkg_sel = bkg_region[iso] < default_cut

          sig_num = sig_region.loc[sig_sel, "weight"].sum()
          sig_den = sig_region["weight"].sum()
          sig_eff_def = sig_num / sig_den if sig_den > 0 else 0

          bkg_num = bkg_sel.sum()
          bkg_den = len(bkg_region)
          bkg_eff_def = bkg_num / bkg_den if bkg_den > 0 else 0

          plt.scatter(sig_eff_def, 1 - bkg_eff_def, color="black", zorder=5,
                    label=f"ref regular cut: {default_cut}")
          plt.axvline(sig_eff_def, color="red", lw=1.5, ls=":")

        # --- Best alternative dt variant ---
        base_prefix = iso.rsplit("-dt_", 1)[0]  # e.g. "eg-layerCluster_iso-cone_veto"
        print(base_prefix)
        dt_candidates = [c for c in sig_df.columns if c.startswith(base_prefix) and c != iso]


        idx = np.argmin(np.abs(tpr - sig_eff_def))
        bkg_rej = 1 - fpr[idx]
        record[iso] = {
          "auc": roc_auc,
          "rej@wp": bkg_rej
        }

        best_iso, best_auc, best_curve = None, -1, None
        print(dt_candidates)
        for dt_var in dt_candidates:
            fpr_dt, tpr_dt, auc_dt = compute_roc(sig_region, bkg_region, dt_var, eta_cut=None)
            idx = np.argmin(np.abs(tpr_dt - sig_eff_def))
            bkg_rej = 1 - fpr_dt[idx]

            record[dt_var] = {
              "auc": auc_dt,
              "rej@wp": bkg_rej
            }

            if best_auc < auc_dt:
                best_auc = auc_dt
                best_iso = dt_var
                best_curve = (tpr_dt, fpr_dt, auc_dt)

        if best_iso is not None:
            tpr_dt, fpr_dt, auc_dt = best_curve
            best_iso_name = best_iso.replace("eg-", "").replace("_iso-", " ").replace("_veto", "-ex")
            best_iso_name = re.sub(r"-dt_([0-9.]+)", lambda m: f" ($\\Delta t = {float(m.group(1)):.2f}$)", best_iso_name)

            plt.plot(tpr_dt, 1 - fpr_dt, lw=2, ls="--", color = colors[color_idx],
                     label=f"{best_iso_name}")

    # Cosmetics
    plt.xlabel("Signal efficiency")
    plt.ylabel("Background rejection")
    plt.title(title)
    plt.legend(fontsize=8, loc="best", ncol=2)
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(outname)
    plt.close()

    return record

def plot_s_bkg(signal, background, bins=50, density=True, labels=("Signal", "Background"),
               xlabel="Variable", ylabel="Events", title="Signal vs Background", logy=False, outname=None):
    """
    Plot signal and background distributions on the same histogram.

    Parameters
    ----------
    signal : array-like
        Array of signal values.
    background : array-like
        Array of background values.
    bins : int or sequence, optional
        Number of bins or bin edges.
    density : bool, optional
        If True, normalize histograms to probability density.
    labels : tuple of str, optional
        Labels for signal and background histograms.
    xlabel : str, optional
        Label for the x-axis.
    ylabel : str, optional
        Label for the y-axis.
    title : str, optional
        Title of the plot.
    logy : bool, optional
        If True, use logarithmic y-axis.
    """

    # Compute common bin edges
    data_all = np.concatenate([signal, background])
    counts, bin_edges = np.histogram(data_all, bins=bins)

    plt.figure(figsize=(8,6))
    plt.hist(signal, bins=bin_edges, density=density, alpha=0.6, label=labels[0], histtype="stepfilled")
    plt.hist(background, bins=bin_edges, density=density, alpha=0.6, label=labels[1], histtype="stepfilled")


    plt.xlabel(xlabel, fontsize=14)
    plt.ylabel(ylabel, fontsize=14)
    plt.title(title, fontsize=16)
    plt.legend(fontsize=12)
    if logy:
        plt.yscale("log")
    plt.grid(alpha=0.3)
    plt.tight_layout()
    if outname is not None:
        plt.savefig(outname)
    plt.show()
    plt.close()

def plot_dt_metric(records, base_prefixes, metric="auc", outname = None):
    """
    Plot a single metric vs. dt in one subplot.
    
    Args:
        records (dict): {iso: {metric: value, ...}}, keys like "prefix_dt_50.0" or "prefix_dt_9999.0".
        base_prefixes (list[str]): prefixes to plot.
        metric (str): which metric to plot, e.g. "auc" or "rej@wp".
    """
    colors = plt.cm.tab10.colors  # distinct colors per prefix
    fig, ax = plt.subplots(figsize=(8, 6))

    for p_idx, prefix in enumerate(base_prefixes):
        dts, vals = [], []
        default_val = None

        for iso, vals_dict in records.items():
            if iso.startswith(prefix):
                try:
                    dt_val = float(iso.split("-dt_")[-1])
                    print(dt_val)
                except ValueError:
                    continue

                if abs(dt_val - 9999.0) < 1:  # default no-cut
                    default_val = vals_dict.get(metric, None)
                else:
                    print(dt_val)
                    dts.append(dt_val)
                    vals.append(vals_dict.get(metric, None))

        if not dts and default_val is None:
            continue

        # sort by dt
        if dts:
            order = np.argsort(dts)
            dts = np.array(dts)[order]
            vals = np.array(vals)[order]

            ax.plot(dts, vals,
                    linestyle="-",
                    color=colors[p_idx % len(colors)],
                    marker="o",
                    label=f"{prefix} (with cut)")

        # plot default as horizontal line
        if default_val is not None:
            ax.axhline(default_val,
                       linestyle="--",
                       color=colors[p_idx % len(colors)],
                       label=f"{prefix} (no cut)")

    ax.set_xlabel("dt")
    ax.set_ylabel(metric)
    ax.set_title(f"{metric} vs dt")
    ax.grid(True, linestyle=":")
    ax.legend()

    plt.tight_layout()
    if outname is not None:
      plt.savefig(outname)
    plt.show()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--inputdir")
    parser.add_argument("--sigproc", required=True, help="Signal parquet file path")
    parser.add_argument("--bkgproc", required=True, help="Background parquet file path")
    parser.add_argument("--plotdir")
    args = parser.parse_args()
   
    os.makedirs(args.plotdir, exist_ok = True)

    sig_df = load_data(f"{args.inputdir}/{args.sigproc}_iso.parquet")
    bkg_df = load_data(f"{args.inputdir}/{args.bkgproc}_iso.parquet")
    
    pt_bins = np.linspace(0, 200, 21)
    eta_bins = np.linspace(-3.0, 3.0, 31)
    
    weights = compute_weights(sig_df, bkg_df, pt_bins, eta_bins)
    sig_df["weight"] = weights 
    eta_regions = {
        "eta<2.0": "(abs(`eg-eta`) < 2.0) & (`eg-et` > 30)",
        "eta>2.0": "(abs(`eg-eta`) > 2.0) & (`eg-et` > 30)"
    }
    
    iso_types = [
        # layerCluster
        "eg-layerCluster_iso-cone_veto-dt_9999.0",
        "eg-layerCluster_iso-seed_veto-dt_9999.0",
        "eg-layerCluster_iso-cluster_veto-dt_9999.0",
        # trackster
        "eg-trackster_iso-cone_veto-dt_9999.0",
        "eg-trackster_iso-seed_veto-dt_9999.0",
        "eg-trackster_iso-cluster_veto-dt_9999.0",
        # rechit
        "eg-HGCRecHit_iso-cone_veto-dt_9999.0",
        "eg-HGCRecHit_iso-seed_veto-dt_9999.0",
        "eg-HGCRecHit_iso-cluster_veto-dt_9999.0"
    ]
    
    default_var = "eg-layerCluster_iso-cone_veto-dt_9999.0"
    default_cuts = {"eta<2.0": 130, "eta>2.0": 340}
    

#    for iso_name in sig_df.columns:
#        plot_s_bkg(sig_df[iso_name], bkg_df[iso_name], outname = os.path.join(args.plotdir, f"{iso_name}.png"))

    for region_name, eta_cut in eta_regions.items():
        plot_name = f"roc_comparison_{region_name.replace('>','gt').replace('<','lt')}.png"

        records = plot_roc_comparison(
            sig_df,
            bkg_df,
            iso_vars=iso_types,             # list of iso variable names
            eta_cut=eta_cut,
            default_iso=default_var,        # e.g. "eg-layerCluster_iso-cone_veto-dt_9999.0"
            default_cut=default_cuts[region_name],  # { "eta<2.0": 140, "eta>2.0": 340 }
            title=f"ROC comparison in {region_name}",
            outname=os.path.join(args.plotdir, plot_name),
        )

        print(f"Plot saved: {plot_name}")


        iso_types_simple = [
          # layerCluster
          "eg-layerCluster_iso-cone_veto",
          "eg-layerCluster_iso-seed_veto",
          "eg-layerCluster_iso-cluster_veto",
          # trackster
          "eg-trackster_iso-cone_veto",
          "eg-trackster_iso-seed_veto",
          "eg-trackster_iso-cluster_veto",
          # rechit
          "eg-HGCRecHit_iso-cone_veto",
          "eg-HGCRecHit_iso-seed_veto",
          "eg-HGCRecHit_iso-cluster_veto"
        ]

        for metric in ["auc", "rej@wp"]:
             for iso in ["layerCluster", "trackster", "HGCRecHit"]:
                 base_prefixes = [i for i in iso_types_simple if iso in i]
                 plot_dt_metric(records, base_prefixes, metric=metric, outname = os.path.join(args.plotdir, plot_name.replace("roc_comparison", f"{iso}-{metric}")))


if __name__ == "__main__":
    main()

