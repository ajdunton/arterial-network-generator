# Arterial Network Generator 
Author: Aaron Dunton (ajdunton@gmail.com)

This program generates the topology of community-scale civil infrastructure networks that we call *arterial* networks. The total area is broken into subareas and the arterial network connects the subareas. We generate many realizations of the arterial topology. Dunton and Gardoni (2023) uses this stochastic arterial network generator to evaluate the probability of disconnection of subareas given localized damage to an unknown network (i.e., the location of damage is known, but the components of the infrastructure at that location and the criticality of those components are unknown). Inside each subarea, there is a *capillary* network, and we have another generator procedure for those networks (https://github.com/ajdunton/capillary-network-generator). The case study in Dunton and Gardoni (2023) and included in this repository is for wastewater networks, but modifications can be made to generate other types of civil infrastructure networks.
<p align="center">
  <img src="https://github.com/ajdunton/arterial-network-generator/assets/147078788/ad5881aa-b283-4163-88cd-96fb77bf9736" width="600">
</p>

# Inputs
The inputs, in the input directory, are listed below. Included in this repository are the inputs for a case study location in Irving, Texas. These inputs require significant manual preparation; see Dunton and Gardoni (2023) for details. Note that all geospatial inputs need to be in the same projected coordinate reference system, with units in feet.
- box_bounds: shapefile with a single polygon feature, the total area that the infrastructure network covers, to be divided into subareas
- dem_ft: digital elevation model in feet, used in evaluating edge weights for generating wastewater networks
- streets: shapefile of streets, including the indicator attributes *main* and *known*. *main* indicates whether each street segment is a major road that will be a boundary between subareas and a potential arterial infrastructure line (*main* = 1) or a smaller street segment that is inside a subarea and a potential capillary infrastructure line (*main* = 0). *known* indicates whether a main street segment is known to have an arterial infrastructure line along it (*known* = 1) or not (*known* = 0) (i.e., allowing us to integrate some known lines into the procedure, see the Additional Functionality section below).
- wwtp_loc: shapefile with a single point attribute at the location from which the infrastructure distribute resources or to which the infrastructure collect waste (e.g., the wastewater treatment plant for wastewater networks)

# Outputs
The outputs define the generated arterial networks and the subareas as follows:
- edges, edge_1, edges_2, ... : shapefiles of linestring geometries representing edges for each realization of the arterial network, with attributes defining the source node (*node_0*) and destination node (*node_1*) of each edge
- outlets: shapefile with point geometries representing the outlets for each of the subareas, also nodes in the arterial network, with an attribute defining the node index (*node_ind*)
- subareas: shapefile with polygon geometries representing the subareas, each with a node index (*node_ind*) attribute that identifies the corresponding outlet/node
- termnode_ind: text document with the node index of the terminal node (i.e., node at wwtp_loc)

# Procedure
The following is a brief overview of the 3-step procedure that is implemented in script.py. Note that params.py and funcs.py are modules that are imported into script.py. These modules contain parameters/constants and functions, respectively.

## Step 1. Identify Subareas
Preliminary subareas are identified based on the mainline street network. However, the non-mainline street network in each preliminary subarea may be disconnected, either from the raw data or from manual pre-processing as described in Dunton and Gardoni (2023). Each connected component of the non-mainline street network corresponds to a different capillary infrastructure network. We determine one final subarea for each connected component of the non-mainline street network in each preliminary subarea.

## Step 2. Identify Arterial Nodes (i.e., Outlets for Wastewater Networks)
We identify the location where the capillary infrastructure in each subarea meet the arterial infrastructure by selecting from the points where the non-mainline streets in the subarea intersect the mainline streets. We include two options for how to determine this point, specified in params.py as *outlet_method*.
- outlet_method = 'elev': choose the point with the lowest elevation
- outlet_method = 'wwtp': choose the point closest to the wastewater treatment plant (based on distance or weight, as is described in Step 3)

## Step 3. Identify Arterial Network Topology
The arterial network is a subtree of the mainline tree network, spanning the outlet nodes. In particular, we connect each outlet node to the wastewater treatment plant via the shortest path. We include two options for how to define the edge weights that are used in the shortest path tree, specified in params.py as *ls_weight_method*.
- ls_weight_method = 'distance': the weight of each edge is the distance of the corresponding street segment
- ls_weight_method = 'hyrbid-elev': the weight of each directed edge is the distance multiplied by a factor that accounts for the elevation of the destination node (i.e., downward flow is preferred). See Dunton and Gardoni (2023) for the formulation of this weight.

# Generating Multiple Networks for Risk Analysis of Unknown Infrastructure to Localized Damage
We generate multiple realizations of the network topology by applying multiplicative random error to the edge weights in Step 3. The number of realizations to be generated is specified in params.py as *number_realizations*. We sample the multiplicative random error independently for each edge from a lognormal distribution with mean = 1 and various values of standard deviation (i.e., sensitivity analysis). See Dunton and Gardoni (2023) for a description of why this method is appropriate for estimating the probability of disconnection of each of the subareas given localized damage, e.g., producing the following map:
<p align="center">
  <img src="https://github.com/ajdunton/arterial-network-generator/assets/147078788/ed6ba23a-8bfd-423e-8cc2-4653e6fbeb73" width="600">
</p>

# Additional Functionality
The arterial network generator can be used for additional analysis that was not included in Dunton and Gardoni (2023). We include some description of additional analysis here.

## Known Infrastructure Lines
We consider the effect of known infrastructure lines, specifically if a portion of the network extending from the treatment plant is known. We reduce the weights of the edges that are known arterial infrastructure lines to almost 0 by multiplying by a factor specified in params.py as *known_weight_ratio*. For example, when a mainline (shown below in green) is known to extend to the west of the study area, the extent of disconnection is likely smaller than the extent shown above.
<p align="center">
  <img src="https://github.com/ajdunton/arterial-network-generator/assets/147078788/e1de4f9a-8740-4d79-a88d-47b3ee541637" width="250">
</p>

## Integrating with the Capillary Network Generator
We integrate this arterial network generator procedure with our previously-published capillary network generator (https://github.com/ajdunton/capillary-network-generator). For each subarea, we clip the input data needed for the capillary network generator to the subarea extent. We then generate the capillary topology. We assess how the results for the risk analysis change when we consider the capillary topology. In general, for the case study location and a single network generated, only 62 additional buildings of the 42,605 in the study area were disconnected when capillary networks were generated. The figure below demonstrates a specific location where there was some difference in the building-level results. Overall, this result indicates that if the subareas are defined at the resolution that we used for the case study, the effect of the capillary networks on the building-level results is negligible, and the above analysis is valid.
<p align="center">
  <img src="https://github.com/ajdunton/arterial-network-generator/assets/147078788/1663ba4d-1396-497d-ba6f-c069cbfa2532" width="700">
</p>

# Reference
See the following paper for further details about this procedure and case study:

Dunton, A., and Gardoni, P. (2023). Developing digital twins of infrastructure for risk analysis. 14th International Conference on Applications of Statistics and Probability in Civil Engineering (ICASP14), Dublin, Ireland, 2023.
