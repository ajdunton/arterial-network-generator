# Arterial Network Generator 
Author: Aaron Dunton (ajdunton@gmail.com)

This program generates the topology of community-scale civil infrastructure networks that we call *arterial* networks. The total area is broken into subareas and the arterial network connects the subareas. We generate many realizations of the arterial topology. Dunton and Gardoni (2023) uses this stochastic arterial network generator to evaluate the probability of disconnection of subareas given localized damage to an unknown network (i.e., the location of damage is known, but the components of the infrastructure at that location and the criticality of those components are unknown). Inside each subarea, there is a *capillary* network, and we have another generator procedure for those networks (https://github.com/ajdunton/capillary-network-generator). The case study in Dunton and Gardoni (2023) and included in this repository is for wastewater networks, but modifications can be made to generator other types of civil infrastructure networks.

# Inputs
The inputs, in the input directory, are listed below. Included in this repository are the inputs for a case study location in Irving, Texas. Note that all geospatial inputs need to be in the same projected coordinate reference system, with units in feet.
- box_bounds: shapefile with a single polygon feature, the total area that the infrastructure network covers
- dem_ft: digital elevation model in feet
- streets: shapefile of streets, including the indicator attributes *main* and *known*. *main* indicates whether each street segment is a major road that will be a boundary between subareas and a potential arterial infrastructure line (*main*=1) or a smaller street segment that is inside a subarea and a potential capillary infrastructure line (*main*=0). *known* indicates whether a main street segment is known to have an arterial infrastructure line along it (*known*=1) or not (*known*=0) (i.e., allowing us to integrate some known lines into the procedure). See Dunton and Gardoni (2023) for more details about preparation of this input.
- wwtp_loc: shapefile with a single point attribute at the location from which the infrastructure distribute resources or to which the infrastructure collect waste (e.g., the wastewater treatment plant for wastewater networks)

# Outputs

# Procedure
The following is a brief overview of the 3-step procedure that is implemented in script.py. Note that params.py and funcs.py are modules that are imported into script.py. These modules contain parameters/constants and functions, respectively.

## 1. 

## 2.

## 3.

# Reference
See the following paper for further details about this procedure and results for the case study:

Dunton, A., and Gardoni, P. (2023). Developing digital twins of infrastructure for risk analysis. 14th International Conference on Applications of Statistics and Probability in Civil Engineering (ICASP14), Dublin, Ireland, 2023.

# Additional Analysis
In addition to the analysis in Dunton and Gardoni (2023), the arterial network generator can be used for additional analysis. We include some description of additional analysis here.

## Known Infrastructure Lines

## Integrating with the Capillary Network Generator
Note that this program interfaces well with the capillary network generator available at __
