import sys
import traceback

import click
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


def handle_error(e: Exception):
    click.secho(f"Error: {e}", fg="red", bold=True)
    with open("errors.log", "a") as fa:
        fa.write(f"{e}\n{traceback.format_exc()}\n")
    sys.exit(1)


def reduce_dimensions_to_2d(df):
    df_values = df.values
    scaled_df_values = StandardScaler().fit_transform(df_values)

    pca = PCA(n_components=2)
    principal_components = pca.fit_transform(scaled_df_values)
    df_principal = pd.DataFrame(data=principal_components, columns=["X", "Y"])
    return df_principal
