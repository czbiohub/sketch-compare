#exploring
from collections import namedtuple

import seaborn as sns
import scipy.cluster.hierarchy as sch
from scipy.spatial.distance import pdist, squareform
import plotly.graph_objects as go
import plotly.figure_factory as ff
from numpy import nan


def create_clustermap(pivot_w_labels, complete_metadata):

    ### Create/Add Dendograms ###
    # Initialize figure by creating upper dendrogram
    fig = ff.create_dendrogram(
        pivot_w_labels,
        orientation='bottom',
        distfun=lambda x: pdist(x, metric="correlation"),
        linkagefun=lambda x: sch.linkage(x, "ward", optimal_ordering=True),
        color_threshold=0.5,
        labels=list(pivot_w_labels.index.get_level_values(0)),
    )

    for i in range(len(fig['data'])):
        fig['data'][i]['yaxis'] = 'y2'

    # Create Side Dendrogram
    dendro_side = ff.create_dendrogram(
        pivot_w_labels,
        orientation='right',
        distfun=lambda x: pdist(x, metric="correlation"),
        color_threshold=0.4,
        linkagefun=lambda x: sch.linkage(x, "ward", optimal_ordering=True),    
    )

    for i in range(len(dendro_side['data'])):
        dendro_side['data'][i]['xaxis'] = 'x2'

    # Add Side Dendrogram Data to Figure
    for data in dendro_side['data']:
        fig.add_trace(data)

    ### Add in color bar/cluster labels ###
    color_bar_meta = namedtuple('colorbar_data', [
        'labeling_type',
        'index',
        'color_palette',
        'xaxis_label',
        'yaxis_label'
        ]
    )
    
    color_bar_data = [
        color_bar_meta('species', 2, 'bright', 'x5', 'y5'),
        color_bar_meta('sourmash_species', 4, 'deep', 'x4', 'y4'),
        color_bar_meta('ska_species', 6, 'dark', 'x3', 'y3')
    ]
    
    # sample_id to species lookup dictionary
    species_lookup = pivot_w_labels.reset_index()[[
        "Sample_1",
        "species",
        "sourmash_species",
        "ska_species"
    ]].set_index("Sample_1").to_dict('index')
    
    for species_meta in color_bar_data:
        
        species_coloring = dict(
            zip(
                pivot_w_labels.index.get_level_values(species_meta.index).unique(), 
                sns.color_palette(species_meta.color_palette, 13)
            )
        )
        # incase nan value color grey
        species_coloring[nan] = (220, 220, 220)
        
        species_labels = [
            species_lookup[i][species_meta.labeling_type] 
            for i in list(fig['layout']['xaxis']['ticktext'])
        ]
        
        species_marker_colors = [
            "rgb"+ str(species_coloring[species_lookup[i][species_meta.labeling_type]])
            for i in list(fig['layout']['xaxis']['ticktext'])   
        ]
        
        bar_color_labels = go.Bar(
            x=list(pivot_w_labels.index.get_level_values(0)),
            y=[10 for i in pivot_w_labels.index.get_level_values(0)],
            marker_color=species_marker_colors,
            hovertext=species_labels,
            xaxis=species_meta.xaxis_label,
            yaxis=species_meta.yaxis_label,
        )
        
        fig.add_trace(bar_color_labels)

    ### Create Heatmap ###
    dendro_leaves = dendro_side['layout']['yaxis']['ticktext']
    dendro_leaves = list(map(int, dendro_leaves))

    data_dist = pdist(pivot_w_labels, metric="correlation")
    heat_data = squareform(data_dist)
    heat_data = heat_data[dendro_leaves, :]
    heat_data = heat_data[:, dendro_leaves]

    heatmap = [
        go.Heatmap(
            x=dendro_leaves,
            y=dendro_leaves,
            z=heat_data,
            colorscale='magma',
            #reversescale = True,
        )
    ]

    heatmap[0]['x'] = fig['layout']['xaxis']['tickvals']
    heatmap[0]['y'] = dendro_side['layout']['yaxis']['tickvals']

    # Add Heatmap Data to Figure
    for data in heatmap:
        fig.add_trace(data)
        
    ### Edit Layout ###
    fig.update_layout({
        'width': 1200, 
        'height': 1200,
        'showlegend': False, 
        'hovermode': 'closest',
        'font': dict(size=8),
    })

    # Edit xaxis
    fig.update_layout(
        xaxis={
            'domain': [.15, 1],
            'mirror': False,
            'showgrid': False,
            'showline': False,
            'zeroline': False,
            'showticklabels': False,                     
            #'tickangle': 45,
            'ticks': ""
        }
    )
    # Edit xaxis2
    fig.update_layout(
        xaxis2={
            'domain': [0, .15],
            'mirror': False,
            'showgrid': False,
            'showline': False,
            'zeroline': False,
            'showticklabels': False,
            #'tickangle': 45,
            'ticks': ""
        }
    )

    # Edit xaxis3
    fig.update_layout(
        xaxis3={
            'domain': [0.15, 1],
            'mirror': False,
            'showgrid': False,
            'showline': False,
            'zeroline': False,
            'showticklabels': False,
            'tickangle': 45,
            'ticks': ""
        }
    )

    # Edit xaxis4
    fig.update_layout(
        xaxis4={
            'domain': [0.15, 1],
            'mirror': False,
            'showgrid': False,
            'showline': False,
            'zeroline': False,
            'showticklabels': False,
            'tickangle': 45,
            'ticks': ""
        }
    )
    # Edit xaxis5
    fig.update_layout(
        xaxis5={
            'domain': [0.15, 1],
            'mirror': False,
            'showgrid': False,
            'showline': False,
            'zeroline': False,
            'showticklabels': False,
            'tickangle': 45,
            'ticks': ""
        }
    )

    # Edit yaxis
    fig.update_layout(
        yaxis={
            'domain': [0, .85],
            'mirror': False,
            'showgrid': False,
            'showline': False,
            'zeroline': False,
            'showticklabels': False,
            'ticks': ""
        }
    )
    # Edit yaxis2
    fig.update_layout(
        yaxis2={
            'domain': [.87, .975],
            'mirror': False,
            'showgrid': False,
            'showline': False,
            'zeroline': False,
            'showticklabels': False,
            'ticks': ""
        }
    )

    # Edit yaxis3
    fig.update_layout(
        yaxis3={
            'domain': [0.83, 0.84],
            'mirror': False,
            'showgrid': False,
            'showline': False,
            'zeroline': False,
            'showticklabels': False,
            'ticks': ""
        }
    )
    
    # Edit yaxis4
    fig.update_layout(
        yaxis4={
            'domain': [0.84, 0.85],
            'mirror': False,
            'showgrid': False,
            'showline': False,
            'zeroline': False,
            'showticklabels': False,
            'ticks': ""
        }
    )
    
    # Edit yaxis5
    fig.update_layout(
        yaxis5={
            'domain': [0.85, 0.86],
            'mirror': False,
            'showgrid': False,
            'showline': False,
            'zeroline': False,
            'showticklabels': False,
            'ticks': ""
        }
    )

    return fig
