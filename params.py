# Author: Aaron Dunton (ajdunton@gmail.com)

# How to determine the outlet for each subarea: 'elev' or 'wwtp'
# elev: candidate point with lowest elevation
# wwtp: closest candidate point to the wwtp
outlet_method = 'wwtp' 

# How to evaluate the weights for each edge: 'distance' or 'hybrid-elev'
ls_weight_method = 'hybrid-elev'

# Number of edge realizatinos to generate
number_realizations = 100

# How to modify the weight of edges that are known mainlines?
# Value between 0 and 1, not exactly 0 to prevent possibility of 0-weight cycle
known_weight_ratio = 0.01 