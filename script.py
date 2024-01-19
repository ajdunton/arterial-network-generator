# Author: Aaron Dunton (ajdunton@gmail.com)

import os
import geopandas as gpd
import rasterio
from shapely.geometry import Point
from shapely.ops import substring
import numpy as np
import networkx as nx
import json

import funcs as fn
import params as pr

# %% Determine sub-areas
# If streets_clean exists, open directly
if os.path.exists('./inter/streets_clean/streets_clean.shp'):
    streets = gpd.read_file('./inter/streets_clean/streets_clean.shp')

# If streets_clean does not exist, call function to clean the street shapefile
# Also saves the cleaned streets and returns the GeoDataFrame
else: 
    streets = fn.street_cleaning()

coref = streets.crs

# Main streets and total area geodataframes
main_streets = streets.loc[streets.main==1]
tot_area = gpd.read_file('./inputs/box_bounds/box_bounds.shp')

# Get first sub-areas based on the mainline geometries
subareas = fn.first_cut(main_streets.geometry, tot_area)

# Further refine sub-areas if the street network in each area is not connected
# Output list of new area geometries
all_new_areas = []
# Output list of indices to drop from subareas
drop_inds = []

# Loop over each of the subareas already identified based on the mainlines
for ind, geom_orig in zip(subareas.index, subareas.geometry):  
    # Clip streets to area and get graph
    s_in_area = streets.clip(geom_orig, keep_geom_type=True).loc[streets.main
                                                                 ==0]
    
    if len(s_in_area) > 0:
        s_in_area.to_file('./inter/cur_sts')
        cur_g, cur_pos = fn.graph_from_shp('./inter/cur_sts/cur_sts.shp')

        # If the graph is not connected, refine the area and add data to output 
        # lists
        if not nx.is_connected(cur_g):
            new_area_list = fn.refine_area(s_in_area, geom_orig, cur_g, 
                                           cur_pos)       
            for i in new_area_list: all_new_areas.append(i)
            drop_inds.append(ind)
            
    else: drop_inds.append(ind)

# Drop areas that were split
subareas = subareas.drop(index=drop_inds)

# Add new areas
all_new_areas_gdf = gpd.GeoDataFrame(geometry=all_new_areas).set_crs(coref)
subareas = subareas.append(all_new_areas_gdf).reset_index(drop=True)
subareas['com_id'] = subareas.index
        
# %% Determine outlets
with rasterio.open('inputs/dem_ft/dem_ft.tif') as dem:
    ls_graph, ls_pos, tn_id = fn.main_graph(main_streets, dem)

subareas['outlets'] = subareas.apply(fn.identify_outlet, axis=1, sts=streets, 
                                     graph=ls_graph, pos=ls_pos)

outlet_node_inds = [fn.identify_ind(pt,ls_pos) for pt in subareas.outlets]
subareas['node_inds'] = outlet_node_inds

# %% Determine path between outlets
# Shortest path tree
ls_spst = fn.main_spst(ls_graph, tn_id, outlet_node_inds)

# %% Save large-scale topology outputs 
# Make outputs directory if it does not exist
if not os.path.isdir('outputs'):
    os.mkdir('outputs')

# Save subareas
out = gpd.GeoDataFrame({'node_ind': subareas['node_inds'],
                        'geometry': subareas['geometry']}).set_crs(coref)
out.to_file('outputs/subareas')

# Save outlets
out = gpd.GeoDataFrame({'node_ind': subareas['node_inds'],
                        'geometry': subareas['outlets']}).set_crs(coref)
out.to_file('outputs/outlets')

# Save edges
spst_edge_geoms = []
nodes_0 = []
nodes_1 = []
for i,j in ls_spst.edges():
    pt_i = Point(ls_pos[i])
    pt_j = Point(ls_pos[j])
    ln = main_streets.loc[(main_streets.intersects(pt_i)) &
                             (main_streets.intersects(pt_j))].geometry
    if pt_i.buffer(1).intersects(ln[ln.index[0]].boundary[0]):
        spst_edge_geoms.append(ln[ln.index[0]])
    else:     
        spst_edge_geoms.append(substring(ln[ln.index[0]], 
                                         ln[ln.index[0]].length,
                                         0))
    nodes_0.append(i)
    nodes_1.append(j)
    
out = gpd.GeoDataFrame({'node_0': nodes_0, 'node_1': nodes_1,
                        'geometry':spst_edge_geoms}).set_crs(coref)

out.to_file('outputs/edges')

# Save termninal node index
with open('./outputs/termnode_ind','w') as handle:
    json.dump(tn_id,handle)

# %% Multiple Network Realizations and Save Outputs
if pr.number_realizations > 1:
    for realiz_number in range(1,pr.number_realizations): 
        # Re-evaluate weight of edges with multiplicative noise
        graph_with_eps = ls_graph.copy()
        n_edges = graph_with_eps.number_of_edges()
        eps_array = np.random.lognormal(mean=0, 
                                        sigma=0.05,
                                        size=n_edges) 
        for ind,(i,j) in enumerate(graph_with_eps.edges()):
            graph_with_eps[i][j]['weight'] = \
                graph_with_eps[i][j]['weight']*eps_array[ind]

        # Shortest path tree
        noise_spst = fn.main_spst(graph_with_eps, tn_id, outlet_node_inds)
        
        # Save edges
        spst_edge_geoms = []
        nodes_0 = []
        nodes_1 = []
        for i,j in noise_spst.edges():
            pt_i = Point(ls_pos[i])
            pt_j = Point(ls_pos[j])
            ln = main_streets.loc[(main_streets.intersects(pt_i)) &
                                  (main_streets.intersects(pt_j))].geometry 
            if pt_i.buffer(1).intersects(ln[ln.index[0]].boundary[0]):
                spst_edge_geoms.append(ln[ln.index[0]])
            else:
                spst_edge_geoms.append(substring(ln[ln.index[0]], 
                                                 ln[ln.index[0]].length,
                                                 0))
            nodes_0.append(i)
            nodes_1.append(j)
            
        out = gpd.GeoDataFrame({'node_0': nodes_0, 'node_1': nodes_1,
                                'geometry':spst_edge_geoms}).set_crs(coref)
        out.to_file('outputs/edges_'+str(realiz_number))