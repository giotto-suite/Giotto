"""Small real-data-format checks; never touch the saved research results."""
import importlib.util
import tempfile
import unittest
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import torch
from scipy import sparse

adapter_path = Path(__file__).resolve().parents[2] / "inst/python/python_reconst.py"
spec = importlib.util.spec_from_file_location("giotto_reconst", adapter_path)
adapter = importlib.util.module_from_spec(spec)
spec.loader.exec_module(adapter)
run_reconst_pipeline = adapter.run_reconst_pipeline


class PipelineTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="reconst_test_")
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        rng = np.random.default_rng(9)
        self.reference = ad.AnnData(
            sparse.csr_matrix(rng.poisson(3, (12, 4)).astype(np.float32)),
            var=pd.DataFrame(
                {"gene_symbol": ["SharedSymbol", "SharedSymbol", "C", "D"]},
                index=["001", "002", "003", "004"],
            ),
        )
        # Different order and one spatial-only gene exercise ID-based alignment.
        self.spatial = ad.AnnData(
            sparse.csr_matrix(rng.poisson(3, (5, 4)).astype(np.float32)),
            var=pd.DataFrame(index=["004", "001", "002", "spatial_only"]),
        )

    def run_pipeline(self, folder="results", **kwargs):
        sc_path, sp_path = self.root / "sc.h5ad", self.root / "sp.h5ad"
        self.reference.write_h5ad(sc_path)
        self.spatial.write_h5ad(sp_path)
        options = dict(
            sc_h5ad=str(sc_path), spatial_h5ad=str(sp_path),
            output_dir=str(self.root / folder), num_epochs=1,
            batch_size=4, embedding_size=2, min_genes_per_cell=None,
            min_cells_per_gene=None, normalize_target_sum=None,
            seed=42, device="cpu", verbose=False,
        )
        options.update(kwargs)
        return run_reconst_pipeline(**options)

    def test_artifacts_identifiers_and_repeatability(self):
        previous_threads = torch.get_num_threads()
        first = self.run_pipeline("first")
        self.assertEqual(torch.get_num_threads(), previous_threads)
        second = self.run_pipeline("second")
        for name in ("results_csv", "model_pt", "history_csv", "metrics_json"):
            self.assertTrue(Path(first[name]).is_file())
        frame = pd.read_csv(first["results_csv"], dtype={"gene_id": str})
        self.assertEqual(frame.gene_id.tolist(), ["001", "002", "004"])
        self.assertEqual(frame.gene_symbol.tolist(), ["SharedSymbol", "SharedSymbol", "D"])
        self.assertEqual(first["n_common_genes"], 3)
        self.assertTrue(np.isfinite(first["loss_all_genes"]))
        self.assertTrue(np.isfinite(first["loss_selected_genes"]))
        np.testing.assert_allclose(
            frame.importance,
            pd.read_csv(second["results_csv"]).importance,
            rtol=0, atol=0,
        )
        self.assertEqual(int(frame.selected.sum()), first["n_selected_genes"])

    def test_no_common_ids(self):
        self.spatial.var_names = ["x", "y", "z", "w"]
        previous_threads = torch.get_num_threads()
        with self.assertRaisesRegex(ValueError, "No common feature IDs"):
            self.run_pipeline()
        self.assertEqual(torch.get_num_threads(), previous_threads)

    def test_duplicate_ids(self):
        self.reference.var_names = ["001", "001", "003", "004"]
        with self.assertRaisesRegex(ValueError, "nonempty and unique"):
            self.run_pipeline()

    def test_filters_leave_no_training_data(self):
        with self.assertRaisesRegex(ValueError, "filtering/splitting"):
            self.run_pipeline(min_genes_per_cell=200)


if __name__ == "__main__":
    unittest.main()
