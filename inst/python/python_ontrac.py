from optparse import Values

import pandas as pd
from ONTraC.integrate import run_ontrac
from ONTraC.utils import write_version_info


def prepare_ontrac_options(preprocessing_dir: str,
                           GNN_dir: str,
                           NTScore_dir: str,
                           n_cpu: int = 4,
                           n_neighbors: int = 50,
                           n_local: int = 20,
                           device: str = 'cpu',
                           epochs: int = 1000,
                           patience: int = 100,
                           min_delta: float = 0.001,
                           min_epochs: int = 50,
                           batch_size: int = 0,
                           seed: int = 42,
                           lr: float = 0.03,
                           hidden_feats: int = 4,
                           k: int = 6,
                           modularity_loss_weight: float = 0.3,
                           purity_loss_weight: float = 300.0,
                           regularization_loss_weight: float = 0.1,
                           beta: float = 0.03) -> Values:
    """Prepare the options for the ONTraC program.
    :param dataset: str, the dataset to be used.
    :param preprocessing_dir: str, the directory of the preprocessed data.
    :param GNN_dir: str, the directory of the GNN model.
    :param NTScore_dir: str, the directory of the NTScore model.
    :param n_cpu: int, the number of CPUs to be used in niche network construction.
    :param n_neighbors: int, the number of neighbors in niche network construction.
    :param n_local: int, the index of local neighbors used for normalization in niche features calculation.
    :param device: str, the device to be used in GNN training.
    :param epochs: int, the maximum number of epochs in GNN training.
    :param patience: int, the number of epochs to wait before early stopping. 
    :param min_delta: float, the minimum delta to be considered as improvement in early stopping.
    :param min_epochs: int, the minimum epochs to be trained.
    :param batch_size: int, the batch size for each iteration in GNN training.
    :param seed: int, the seed.
    :param lr: float, the learning rate in GNN training.
    :param hidden_feats: int, the hidden features in GNN model.
    :param k: int, the number of niche clusters.
    :param modularity_loss_weight: float, the modularity loss weight.
    :param purity_loss_weight: float, the purity loss weight.
    :param regularization_loss_weight: float, the regularization loss weight.
    :param beta: float, the beta value in softmax function.
    :return: Values, the options.
    """

    options = Values()
    options.preprocessing_dir = preprocessing_dir
    options.GNN_dir = GNN_dir
    options.NTScore_dir = NTScore_dir
    options.n_cpu = n_cpu
    options.n_neighbors = n_neighbors
    options.n_local = n_local
    options.device = device
    options.epochs = epochs
    options.patience = patience
    options.min_delta = min_delta
    options.min_epochs = min_epochs
    options.batch_size = batch_size
    options.seed = seed
    options.lr = lr
    options.hidden_feats = hidden_feats
    options.k = k
    options.modularity_loss_weight = modularity_loss_weight
    options.purity_loss_weight = purity_loss_weight
    options.regularization_loss_weight = regularization_loss_weight
    options.beta = beta

    return options


def ONTraC(ONTraC_input: pd.DataFrame,
           preprocessing_dir: str,
           GNN_dir: str,
           NTScore_dir: str,
           n_cpu: int = 4,
           n_neighbors: int = 50,
           n_local: int = 20,
           device: str = 'cpu',
           epochs: int = 1000,
           patience: int = 100,
           min_delta: float = 0.001,
           min_epochs: int = 50,
           batch_size: int = 0,
           seed: int = 42,
           lr: float = 0.03,
           hidden_feats: int = 4,
           k: int = 6,
           modularity_loss_weight: float = 0.3,
           purity_loss_weight: float = 300.0,
           regularization_loss_weight: float = 0.1,
           beta: float = 0.03) -> None:
    """Run the ONTraC program.
    :param ONTraC_input: pd.DataFrame, the input data.
    :param dataset: str, the dataset to be used.
    :param preprocessing_dir: str, the directory of the preprocessed data.
    :param GNN_dir: str, the directory of the GNN model.
    :param NTScore_dir: str, the directory of the NTScore model.
    :param n_cpu: int, the number of CPUs to be used in niche network construction.
    :param n_neighbors: int, the number of neighbors in niche network construction.
    :param n_local: int, the index of local neighbors used for normalization in niche features calculation.
    :param device: str, the device to be used in GNN training.
    :param epochs: int, the maximum number of epochs in GNN training.
    :param patience: int, the number of epochs to wait before early stopping. 
    :param min_delta: float, the minimum delta to be considered as improvement in early stopping.
    :param min_epochs: int, the minimum epochs to be trained.
    :param batch_size: int, the batch size for each iteration in GNN training.
    :param seed: int, the seed.
    :param lr: float, the learning rate in GNN training.
    :param hidden_feats: int, the hidden features in GNN model.
    :param k: int, the number of niche clusters.
    :param modularity_loss_weight: float, the modularity loss weight.
    :param purity_loss_weight: float, the purity loss weight.
    :param regularization_loss_weight: float, the regularization loss weight.
    :param beta: float, the beta value in softmax function."""

    write_version_info()

    options = prepare_ontrac_options(preprocessing_dir=preprocessing_dir,
                                     GNN_dir=GNN_dir,
                                     NTScore_dir=NTScore_dir,
                                     n_cpu=n_cpu,
                                     n_neighbors=n_neighbors,
                                     n_local=n_local,
                                     device=device,
                                     epochs=epochs,
                                     patience=patience,
                                     min_delta=min_delta,
                                     min_epochs=min_epochs,
                                     batch_size=batch_size,
                                     seed=seed,
                                     lr=lr,
                                     hidden_feats=hidden_feats,
                                     k=k,
                                     modularity_loss_weight=modularity_loss_weight,
                                     purity_loss_weight=purity_loss_weight,
                                     regularization_loss_weight=regularization_loss_weight,
                                     beta=beta)

    run_ontrac(options=options, ori_data_df=ONTraC_input)
