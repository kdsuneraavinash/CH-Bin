import pandas as pd
from sklearn.metrics import adjusted_rand_score


def get_confusion_matrix(df_gt_bin: pd.DataFrame) -> pd.DataFrame:
    return df_gt_bin.groupby(["SPECIES", "BIN"]).size().unstack(fill_value=0)


def get_precision(df_cm: pd.DataFrame, n_binned: int):
    max_s = df_cm.max(axis=0)
    sigma_k_max_s = max_s.sum()
    return sigma_k_max_s / n_binned


def get_recall(df_cm: pd.DataFrame, n_all: int):
    max_k = df_cm.max(axis=1)
    sigma_s_max_k = max_k.sum()
    return sigma_s_max_k / n_all


def get_f1(df_cm: pd.DataFrame, n_binned: int, n_all: int):
    precision = get_precision(df_cm, n_binned=n_binned)
    recall = get_recall(df_cm, n_all=n_all)
    return (2 * precision * recall) / (precision + recall)


def get_ari(df_gt_bin: pd.DataFrame):
    return adjusted_rand_score(labels_true=df_gt_bin.SPECIES, labels_pred=df_gt_bin.BIN)
