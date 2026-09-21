"""Giotto file-in/file-out adapter for the public ReconST Python API.

The model, training objective, and selection rule live in the independent
ReconST package. This adapter uses the public ReconST 0.1.0 core API.
"""

from contextlib import nullcontext, redirect_stdout
import json
import os
import random
from typing import Optional

import numpy as np
import pandas as pd
import scanpy as sc
import torch
from torch.utils.data import DataLoader

from reconst.data import GeneExpressionDataset, create_data_loader, prepare_common_genes
from reconst.model import FeatureScreeningAutoencoder
from reconst.trainer import evaluate_model, select_genes, train_model


def _resolve_device(device):
    if device is not None:
        return torch.device(device)
    if torch.cuda.is_available():
        return torch.device("cuda")
    if getattr(torch.backends, "mps", None) is not None and torch.backends.mps.is_available():
        return torch.device("mps")
    return torch.device("cpu")


def _set_seed(seed):
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)


def _gene_symbols(adata):
    if "gene_symbol" in adata.var.columns:
        return adata.var["gene_symbol"].astype(str).tolist()
    return adata.var_names.astype(str).tolist()


def run_reconst_pipeline(
    sc_h5ad: str,
    spatial_h5ad: str,
    output_dir: str,
    embedding_size: int = 200,
    num_epochs: int = 20,
    batch_size: int = 256,
    lr: float = 1e-3,
    weight_decay: float = 1e-5,
    l_lambda: float = 1e-4,
    dropout: float = 0.05,
    leaky_slope: float = 0.3,
    train_split: float = 0.8,
    threshold: float = 0.001,
    min_genes_per_cell: Optional[int] = 200,
    min_cells_per_gene: Optional[int] = 100,
    normalize_target_sum: Optional[float] = 1e4,
    seed: int = 42,
    device: Optional[str] = None,
    verbose: bool = True,
    num_threads: int = 1,
) -> dict:
    """Run the full ReconST gene-panel selection pipeline.

    Parameters
    ----------
    sc_h5ad : str
        Path to the scRNA-seq AnnData .h5ad. ``.X`` should hold raw counts
        and ``var_names`` should be unique feature IDs matching ``spatial_h5ad``
        (gene symbols or Ensembl IDs, using the same convention in both).
    spatial_h5ad : str
        Path to the paired spatial AnnData .h5ad (e.g. MERFISH).
    output_dir : str
        Directory to write results into. Created if missing.
    embedding_size, num_epochs, batch_size, lr, weight_decay, l_lambda,
    dropout, leaky_slope, train_split, threshold :
        Model and training hyperparameters.
    min_genes_per_cell, min_cells_per_gene, normalize_target_sum :
        scRNA-seq preprocessing. Pass ``None`` on any to skip that step.
    seed : int
        Seed for Python, NumPy, and PyTorch RNGs so the selected panel is
        reproducible across runs.
    device : str, optional
        Torch device string. Auto-detected (cuda → mps → cpu) when None.
    verbose : bool
        Print per-epoch losses and summary lines.
    num_threads : int
        Torch CPU threads during this call. Defaults to one for compatibility
        with R native libraries. The previous setting is restored afterwards.

    Returns
    -------
    dict
        Paths of every artifact written and a few summary statistics.
        Stable keys: ``results_csv``, ``model_pt``, ``history_csv``,
        ``metrics_json``, ``n_common_genes``, ``n_selected_genes``,
        ``loss_all_genes``, ``loss_selected_genes``.
    """

    if not isinstance(num_threads, int) or num_threads < 1:
        raise ValueError("num_threads must be a positive integer.")
    previous_threads = torch.get_num_threads()
    torch.set_num_threads(num_threads)
    try:
        os.makedirs(output_dir, exist_ok=True)
        _set_seed(seed)
        torch_device = _resolve_device(device)

        if verbose:
            print(f"[ReconST] device: {torch_device}")
            print(f"[ReconST] loading {sc_h5ad}")
        adata_sc = sc.read_h5ad(sc_h5ad)
        if verbose:
            print(f"[ReconST] loading {spatial_h5ad}")
        adata_sp = sc.read_h5ad(spatial_h5ad)

        for label, adata in (("scRNA-seq", adata_sc), ("spatial", adata_sp)):
            if adata.n_obs == 0 or adata.n_vars == 0:
                raise ValueError(f"{label} input has no cells or genes.")
            if not adata.var_names.is_unique or any(not str(x).strip() for x in adata.var_names):
                raise ValueError(f"{label} feature IDs must be nonempty and unique.")
        if not 0 < train_split < 1:
            raise ValueError("train_split must be between 0 and 1.")
        if num_epochs < 1 or batch_size < 1 or embedding_size < 1:
            raise ValueError("num_epochs, batch_size, and embedding_size must be positive.")
        if not np.isfinite(threshold):
            raise ValueError("threshold must be finite.")

        if min_genes_per_cell is not None:
            sc.pp.filter_cells(adata_sc, min_genes=min_genes_per_cell)
        if min_cells_per_gene is not None:
            sc.pp.filter_genes(adata_sc, min_cells=min_cells_per_gene)
        train_size = int(train_split * adata_sc.n_obs)
        if train_size < 1 or train_size >= adata_sc.n_obs or adata_sc.n_vars == 0:
            raise ValueError(
                "Reference filtering/splitting leaves no usable training or test data. "
                "Lower min_genes_per_cell/min_cells_per_gene (or set them to None), "
                "and check train_split."
            )
        if normalize_target_sum is not None:
            sc.pp.normalize_total(adata_sc, target_sum=normalize_target_sum)

        adata_sc_common, adata_sp_common, common_genes = prepare_common_genes(adata_sc, adata_sp)
        if not common_genes:
            raise ValueError("No common feature IDs; both inputs must use the same ID convention.")
        if verbose:
            print(f"[ReconST] common genes: {len(common_genes)}")

        train_loader, test_loader_sc = create_data_loader(
            adata_sc_common, batch_size=batch_size, train_split=train_split
        )

        sp_matrix = (
            adata_sp_common.X.toarray()
            if hasattr(adata_sp_common.X, "toarray")
            else adata_sp_common.X
        )
        sp_dataset = GeneExpressionDataset(torch.tensor(sp_matrix, dtype=torch.float32))
        test_loader_sp = DataLoader(sp_dataset, batch_size=batch_size, shuffle=False)

        model = FeatureScreeningAutoencoder(
            input_size=adata_sc_common.n_vars,
            embedding_size=embedding_size,
            dp=dropout,
            lk=leaky_slope,
        ).to(torch_device)

        # ReconST 0.1.0 does not have a verbose argument. Suppress its stdout only
        # when requested, while preserving the training implementation unchanged.
        with open(os.devnull, "w") as quiet_output:
            with nullcontext() if verbose else redirect_stdout(quiet_output):
                train_losses, test_losses = train_model(
                    model,
                    train_loader,
                    test_loader_sc,
                    num_epochs=num_epochs,
                    lr=lr,
                    weight_decay=weight_decay,
                    l_lambda=l_lambda,
                    device=torch_device,
                )

        genes_mask, feature_importances = select_genes(model, threshold=threshold)

        loss_all = evaluate_model(model, test_loader_sp, gene_mask=None, device=torch_device)
        loss_selected = evaluate_model(
            model, test_loader_sp, gene_mask=genes_mask, device=torch_device
        )

        gene_symbols = _gene_symbols(adata_sc_common)
        results_df = pd.DataFrame(
            {
                "gene_id": adata_sc_common.var_names.astype(str).tolist(),
                "gene_symbol": gene_symbols,
                "importance": feature_importances,
                "selected": genes_mask.astype(bool),
            }
        )

        results_csv = os.path.join(output_dir, "reconst_results.csv")
        model_pt = os.path.join(output_dir, "autoencoder_model.pth")
        history_csv = os.path.join(output_dir, "training_history.csv")
        metrics_json = os.path.join(output_dir, "evaluation_metrics.json")

        results_df.to_csv(results_csv, index=False)
        torch.save(model.state_dict(), model_pt)
        pd.DataFrame(
            {
                "epoch": list(range(1, len(train_losses) + 1)),
                "train_loss": train_losses,
                "test_loss": test_losses,
            }
        ).to_csv(history_csv, index=False)

        summary = {
            "n_common_genes": int(len(common_genes)),
            "n_selected_genes": int(genes_mask.sum()),
            "threshold": float(threshold),
            "loss_all_genes": float(loss_all),
            "loss_selected_genes": float(loss_selected),
            "seed": int(seed),
            "device": str(torch_device),
            "num_threads": num_threads,
        }
        with open(metrics_json, "w") as f:
            json.dump(summary, f, indent=2)

        if verbose:
            print(
                f"[ReconST] selected {summary['n_selected_genes']}/"
                f"{summary['n_common_genes']} genes at threshold {threshold}"
            )
            print(f"[ReconST] spatial loss all={loss_all:.4f} selected={loss_selected:.4f}")
            print(f"[ReconST] wrote results to {output_dir}")

        return {
            "results_csv": results_csv,
            "model_pt": model_pt,
            "history_csv": history_csv,
            "metrics_json": metrics_json,
            **summary,
        }
    finally:
        torch.set_num_threads(previous_threads)
