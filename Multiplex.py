import sys
import pandas as pd
import numpy as np
import random as rn
sys.path.insert(0, '/Users/tkay/pymnet')
from pymnet import *
import matplotlib.style
import matplotlib as mpl

## Create multiplex object
mnet = MultiplexNetwork(couplings = "categorical")

## Add layers
layers = ["Social","Gene Expression", "Behavior", "Microbiota", "Physical Environment"]
for x in layers:
    mnet.add_layer(x)

## Add nodes
nodes  = [2,4,5,6,9,14,15,19,27,28,32,35,42,43,52,53,54,59,63,65,80,99,113,115,118,126,135,139,147,149,153,156,159,161,169,170,173,192,207,212,215,218,221,229,237,238,245,250,255,257,260,262,263,264,265,268,278,326,342,353,356,380,390,391,397,437,462,476,482,518,534,540,556,564,575,593,600,610,614,620,632,635,638,646,663,665]
for x in nodes:
    mnet.add_node(x)

## Add edges
BrainEdges =  pd.read_csv('edges_brain.txt', sep=" ", header=None)
for edge in range(0,len(BrainEdges.index)):
    F = BrainEdges[0][edge]
    T = BrainEdges[1][edge]
    mnet[F,'Gene Expression'][T,'Gene Expression'] = 1

SocialEdges = pd.read_csv('edges_social.txt', sep=" ", header=None)
for edge in range(0,len(SocialEdges.index)):
    F = SocialEdges[0][edge]
    T = SocialEdges[1][edge]
    mnet[F,'Social'][T,'Social'] = 1

SpaceEdges = pd.read_csv('edges_space.txt', sep=" ", header=None)
for edge in range(0,len(SpaceEdges.index)):
    F = SpaceEdges[0][edge]
    T = SpaceEdges[1][edge]
    mnet[F,'Physical Environment'][T,'Physical Environment'] = 1

BehaviorEdges = pd.read_csv('edges_task.txt', sep=" ", header=None)
for edge in range(0,len(BehaviorEdges.index)):
    F = BehaviorEdges[0][edge]
    T = BehaviorEdges[1][edge]
    mnet[F,'Behavior'][T,'Behavior'] = 1

GutEdges = pd.read_csv('edges_gut.txt', sep=" ", header=None)
for edge in range(0,len(GutEdges.index)):
    F = GutEdges[0][edge]
    T = GutEdges[1][edge]
    mnet[F,'Microbiota'][T,'Microbiota'] = 1

## Import separately per layer
BrainCoordinates = pd.read_csv('layout_brain.txt', sep=" ", header=None)
SocialCoordinates = pd.read_csv('layout_social.txt', sep=" ", header=None)
BehaviorCoordinates = pd.read_csv('layout_task.txt', sep=" ", header=None)
GutCoordinates = pd.read_csv('layout_gut.txt', sep=" ", header=None)
SpaceCoordinates = pd.read_csv('layout_space.txt', sep=" ", header=None)

Coord = {}
for position in range(0,len(BrainCoordinates.index)):
    node = BrainCoordinates[0][position]
    Coord[node, 'Gene Expression'] = (BrainCoordinates[1][position],BrainCoordinates[2][position])

for position in range(0,len(SocialCoordinates.index)):
    node = SocialCoordinates[0][position]
    Coord[node, 'Social'] = (SocialCoordinates[1][position],SocialCoordinates[2][position])

for position in range(0,len(SpaceCoordinates.index)):
    node = SpaceCoordinates[0][position]
    Coord[node, 'Physical Environment'] = (SpaceCoordinates[1][position],SpaceCoordinates[2][position])

for position in range(0,len(GutCoordinates.index)):
    node = GutCoordinates[0][position]
    Coord[node, 'Microbiota'] = (GutCoordinates[1][position],GutCoordinates[2][position])

for position in range(0,len(BehaviorCoordinates.index)):
    node = BehaviorCoordinates[0][position]
    Coord[node, 'Behavior'] = (BehaviorCoordinates[1][position],BehaviorCoordinates[2][position])

## Read node colours
NodeColors = pd.read_csv('colours.txt', sep=" ", header=None)
NodeCols = {}
for i in range(0,len(list(mnet))):
    node = list(mnet)[i]
    index = np.where(NodeColors[0] == node)
    NodeCols[node, 'Social'] = NodeColors[1][index[0][0]]
    NodeCols[node, 'Gene Expression'] = NodeColors[1][index[0][0]]
    NodeCols[node, 'Microbiota'] = NodeColors[1][index[0][0]]
    NodeCols[node, 'Physical Environment'] = NodeColors[1][index[0][0]]
    NodeCols[node, 'Behavior'] = NodeColors[1][index[0][0]]

## Layer order
LayerOrder = {0: "Social", 1: "Physical Environment", 2: "Gene Expression", 3: "Behavior", 4: "Microbiota"}

## Plot network
rn.seed(3)
mpl.rcParams['mathtext.rm'] = 'serif'
fig=draw(mnet,
         layergap=0.9, #0.8 looks better but gets cut off
         camera_dist=17,
         layerPadding=0,
         nodeColorDict=NodeCols,
         defaultLayerColor="gainsboro",
         defaultEdgeColor="darkgray",
         edgeWidthRule={"rule":"edgetype","intra":0.05,"inter":0.02},
         edgeStyleRule={"rule":"edgetype","intra":"-","inter":":"},
         nodeLabelRule={},
         nodelayerCoords=Coord,
         nodeLabelColorDict=NodeCols,
         defaultLayerLabelSize=7,
         defaultNodeSize=4,
         defaultLayerLabelStyle="italic",
         defaultNodeLabelSize=2,
         layerOrderDict=LayerOrder  # Not working
         )
fig.savefig("mnet.pdf")
