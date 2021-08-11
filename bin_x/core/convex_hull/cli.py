# import numpy as np
# import pandas as pd
#
# from bin_x.core.convex_hull.algorithm import fit_cluster
# from bin_x.core.convex_hull.distance import find_distance_matrix
#
#
# def step_01_apply_algorithm(config: _Config, n_clusters: int, df: pd.DataFrame):
#     samples = df.drop(["CONTIG_NAME", "PARENT_NAME", "CLUSTER"], axis=1).values
#     init_bins = df.CLUSTER.values.copy()
#     n_samples = len(samples)
#
#     print(f"Creating a distance matrix of {n_samples}x{n_samples} shape")
#     mat_dist = find_distance_matrix(samples)
#
#     print(f"Performing binning using {config.qp_solver} solver")
#     return fit_cluster(
#         config.n_iter,
#         n_samples,
#         samples,
#         n_clusters,
#         config.n_neighbors,
#         mat_dist,
#         init_bins,
#         metric=config.metric,
#         qp_solver=config.qp_solver,
#     )
#
#
# def step_02_assign_bin(dirs: _DirManager, n_clusters: int, df: pd.DataFrame, convex_labels: np.ndarray):
#     if np.any(convex_labels < 0):
#         raise ValueError("There were some unclustered points left... Aborting.")
#     df_cluster_column = pd.DataFrame({"BIN": convex_labels})
#     df_bin = df.drop("CLUSTER", axis=1)
#     df_bin = pd.concat([df_bin, df_cluster_column], axis=1)
#     parent_groups = df_bin[["PARENT_NAME", "BIN"]].groupby("PARENT_NAME")
#
#     # The dataframe with bin assigned
#     df_dist_bin = parent_groups.BIN.apply(lambda x: np.bincount(x).argmax()).reset_index()
#     df_dist_bin = df_dist_bin.rename(columns={"PARENT_NAME": "CONTIG_NAME"})
#     df_dist_bin.to_csv(dirs.dist_bin_csv, index=False)
#
#     # The dataframe with probabilities
#     dist_prob = {}
#     for i, (contig_name, group_df) in enumerate(parent_groups):
#         bin_counts = np.bincount(group_df.BIN, minlength=n_clusters)
#         dist_prob[contig_name] = bin_counts / bin_counts.sum()
#     df_dist_prob = pd.DataFrame.from_dict(dist_prob, orient="index", columns=range(n_clusters))
#     df_dist_prob.index.name = "CONTIG_NAME"
#     df_dist_prob = df_dist_prob.reset_index()
#     df_dist_prob.to_csv(dirs.dist_prob_csv, index=False)
#
#     return df_dist_bin, df_dist_prob
#
#
# def step_03_visualize(dirs: _DirManager, n_clusters: int, df: pd.DataFrame, df_dist_bin: pd.DataFrame):
#     df_kmer_bins = df_dist_bin.merge(df, left_on="CONTIG_NAME", right_on="PARENT_NAME")
#     df_kmer_counts = df_kmer_bins.drop(["CONTIG_NAME_x", "CONTIG_NAME_y", "PARENT_NAME", "CLUSTER", "BIN"], axis=1)
#     df_2d = reduce_to_2d(df_kmer_counts)
#
#     plt.figure(figsize=(15, 10))
#     plt.title("Convex hull Result (Split)")
#     for i in range(-1, n_clusters):
#         sns.scatterplot(x="X", y="Y", data=df_2d[df_kmer_bins.BIN == i], label=f"Cluster {i}", alpha=0.1)
#     plt.savefig(dirs.split_bin_png)
#
#     plt.figure(figsize=(15, 10))
#     plt.title("Bin Counts")
#     sns.histplot(df_kmer_bins, x="BIN", discrete=True)
#     plt.savefig(dirs.bin_counts_png)
#
#
# def convex_bin(config_file, df_file, output_dir=None):
#     click.secho("Initializing Convex Binning...", fg="blue", bold=True)
#     config = _Config(config_file)
#
#     if not output_dir:
#         timestamp = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
#         output_dir = Path("out") / config.dataset_name / "bin-out" / timestamp
#     dirs = _DirManager(output_dir, config_file)
#     dirs.create_dirs()
#
#     df = pd.read_csv(df_file)
#     n_clusters = df.CLUSTER.max() + 1
#
#     click.secho(f"01. Running algorithm with {config.metric} metric...", bold=True)
#     convex_labels = step_01_apply_algorithm(config, n_clusters, df)
#
#     click.secho(f"02. Assigning the bin via voting...", bold=True)
#     df_dist_bin, df_dist_prob = step_02_assign_bin(dirs, n_clusters, df, convex_labels)
#
#     step_03_visualize(dirs, n_clusters, df, df_dist_bin)
#     click.secho(f"Dumped outputs on {dirs.binning_dir}", fg="green", bold=True)
#
#     return dirs.dist_bin_csv
#
#
# @click.command()
# @click.option(
#     "--config_file", prompt="Configuration file", help="The configuration JSON to use for dataset generation."
# )
# @click.option("--df_file", prompt="Dataframe file", help="The dataframe file to perform the binning operation.")
# def main(config_file, df_file):
#     try:
#         convex_bin(config_file, df_file)
#     except Exception as e:
#         handle_error(e)
#
#
# if __name__ == "__main__":
#     main()
