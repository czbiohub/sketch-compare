import pandas as pd

import dash
# import dash_bio as dashbio
import dash_html_components as html
import dash_core_components as dcc
import dash_table
from scipy.spatial.distance import pdist, squareform
import scipy.cluster.hierarchy as sch
import scipy as scipy


import plotly.graph_objects as go
import plotly.figure_factory as ff

from utils import join_on_ska_labels
from clustermap import create_clustermap
from data_table import create_datatable


app = dash.Dash(__name__)

### Data Manipulation ###
species_id_skeeters = pd.read_csv(
    "~/code/skeeters/data/sample_genus_and_species.csv",
    header=0,
    index_col=0
)

species_id_skeeters.columns = species_id_skeeters.columns.str.replace(" ", "_")
species_id_skeeters.rename(
    columns={

        'corrected_species': 'sourmash_species',
        'corrected_genus': 'sourmash_genus'
    },
    inplace=True
)
final_pivot_w_labels, cluster_map, linkage = join_on_ska_labels(
    "s3://phoenixlogan-sketches/k15_10000000.distances.tsv",
    species_id_skeeters,
    "_001_10000000_15",
)

# columns and rows for app display
columns = list(final_pivot_w_labels.columns.values)
rows = list(final_pivot_w_labels.index.get_level_values(6).unique())

app.layout = html.Div([
    "Species to display",
    dcc.Dropdown(
        id='species-input',
        options=[
            {
                'label': row,
                'value': row
            }
            for row in rows
        ],
        value=rows,
        multi=True
    ),
    dcc.Graph(id='cluster-map'),
    dash_table.DataTable(
        id='outlier-table',
        columns=[
            {"name": i, "id": i}
            for i in [
                "species",
                "sourmash_species",
                "ska_species"
            ]
        ],
    )
])


@app.callback(
    dash.dependencies.Output('cluster-map', 'figure'),
    [dash.dependencies.Input('species-input', 'value')]
)
def update_clustergram(rows):
    print("ROWS: ", rows)
    return create_clustermap(
        final_pivot_w_labels.loc[
            final_pivot_w_labels.index.get_level_values(6).isin(rows)
        ]
    )


@app.callback(
    dash.dependencies.Output('outlier-table', 'data'),
    [dash.dependencies.Input('species-input', 'value')]
)
def update_outliertable(rows):
    data_w_species = final_pivot_w_labels.loc[
        final_pivot_w_labels.index.get_level_values(6).isin(rows)
    ]

    return data_w_species.reset_index()[[
        "species",
        "sourmash_species",
        "ska_species"
    ]].query(
            'species != sourmash_species | ' + 
            'species != ska_species | ' + 
            'sourmash_species != ska_species'
        ).to_dict('records')


if __name__ == '__main__':
    app.run_server(debug=True)
