#!/usr/bin/env python
import networkx as nx
import matplotlib.pyplot as plt
from fa2 import ForceAtlas2
import pandas as pd
import random
import numpy
import glob
import os
from plotly.offline import download_plotlyjs, init_notebook_mode, iplot
import plotly.graph_objs as go
import plotly
import sys

random.seed(246)
numpy.random.seed(4812)

file_path = sys.argv[1]
dir_path = sys.argv[2]
saveas = sys.argv[3]


def create_nets(file_path, dir_path, saveas):
    nams = os.path.splitext(os.path.basename(file_path))[0]
    df = pd.read_csv(file_path, sep="\t", lineterminator="\r")
    graphtype = nx.Graph()
    G = nx.from_pandas_edgelist(
        df,
        source="source",
        target="target",
        edge_attr=["weight", "length"],
        create_using=graphtype,
    )
    # This gets the sample indices
    samples = [
        key for key, val in enumerate(list(G.nodes)) if val in set(list(df.iloc[:, 0]))
    ]
    color_map = ["grey"] * len(list(G.nodes))
    for index in samples:
        color_map[index] = "red"
    forceatlas2 = ForceAtlas2(
        # Behavior alternatives
        outboundAttractionDistribution=True,  # Dissuade hubs
        linLogMode=False,  # NOT IMPLEMENTED
        adjustSizes=False,  # Prevent overlap (NOT IMPLEMENTED)
        edgeWeightInfluence=0.5,
        # Performance
        jitterTolerance=0.1,  # Tolerance
        barnesHutOptimize=False,
        barnesHutTheta=0,
        multiThreaded=False,  # NOT IMPLEMENTED
        # Tuning
        scalingRatio=5,
        strongGravityMode=False,
        gravity=1,
        # Log
        verbose=True,
    )
    positions = forceatlas2.forceatlas2_networkx_layout(G, pos=None, iterations=5000)
    if saveas == "json":
        size_map = [10] * len(list(G.nodes))
        for index in samples:
            size_map[index] = 25
        plotlyize(
            G=G,
            n_position=positions,
            n_color=color_map,
            n_size=size_map,
            plt_title=nams,
            dir_path=dir_path,
        )
    if saveas == "svg":
        size_map = [20] * len(list(G.nodes))
        for index in samples:
            size_map[index] = 60
        svgize(
            G=G,
            n_position=positions,
            n_color=color_map,
            n_size=size_map,
            plt_title=nams,
            dir_path=dir_path,
            df=df,
        )
    if saveas == "all":
        size_map = [10] * len(list(G.nodes))
        for index in samples:
            size_map[index] = 25
        plotlyize(
            G=G,
            n_position=positions,
            n_color=color_map,
            n_size=size_map,
            plt_title=nams,
            dir_path=dir_path,
        )
        size_map = [20] * len(list(G.nodes))
        for index in samples:
            size_map[index] = 60
        svgize(
            G=G,
            n_position=positions,
            n_color=color_map,
            n_size=size_map,
            plt_title=nams,
            dir_path=dir_path,
            df=df,
        )


def run_nets(file_path, dir_path, saveas):
    fils = glob.glob(file_path + "/*.tsv")
    for i in fils:
        create_nets(file_path=i, dir_path=dir_path, saveas=saveas)


def plotlyize(G, n_position, n_color, n_size, plt_title, dir_path):
    edge_trace = go.Scatter(
        x=[], y=[], line=dict(width=0.5, color="#888"), hoverinfo="none", mode="lines"
    )
    for n, p in n_position.items():
        G.nodes[n]["n_position"] = p
    for edge in G.edges():
        x0, y0 = G.nodes[edge[0]]["n_position"]
        x1, y1 = G.nodes[edge[1]]["n_position"]
        edge_trace["x"] += tuple([x0, x1, None])
        edge_trace["y"] += tuple([y0, y1, None])
    node_trace = go.Scatter(
        x=[],
        y=[],
        text=[],
        mode="markers+text",
        textposition="bottom center",
        hoverinfo="none",
        marker=dict(color=[], size=[], line=dict(width=0)),
    )
    for node in G.nodes():
        x, y = G.nodes[node]["n_position"]
        node_trace["x"] += tuple([x])
        node_trace["y"] += tuple([y])
        node_trace["text"] += tuple([node])
        node_trace["marker"]["color"] = tuple(n_color)
        node_trace["marker"]["size"] = tuple(n_size)
    fig = go.Figure(
        data=[edge_trace, node_trace],
        layout=go.Layout(
            title="",
            titlefont=dict(size=16),
            showlegend=False,
            margin=dict(b=20, l=5, r=5, t=40),
            xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        ),
    )
    fig.write_json(dir_path + "/" + plt_title + ".json")


def svgize(G, n_position, n_color, n_size, plt_title, dir_path, df):
    samps = [
        key for key, val in enumerate(list(G.nodes)) if val in set(list(df.iloc[:, 0]))
    ]
    samps = [list(G.nodes)[i] for i in samps]
    samps = {samps[i]: samps[i] for i in range(len(samps))}
    nx.draw_networkx_nodes(
        G, n_position, with_labels=True, node_color=n_color, alpha=0.8, node_size=n_size
    )
    nx.draw_networkx_edges(G, n_position, edge_color="green", alpha=0.05)
    offset = 15
    pos_labels = {}
    keys = n_position.keys()
    for key in keys:
        x, y = n_position[key]
        pos_labels[key] = (x, y + offset)

    nx.draw_networkx_labels(G, pos_labels, labels=samps, font_size=12)
    figure = plt.gcf()
    figure.set_size_inches(20, 20)
    plt.show(block=False)
    plt.savefig(dir_path + "/svg/" + plt_title + ".svg", format="svg")


run_nets(file_path=file_path, dir_path=dir_path, saveas=saveas)
