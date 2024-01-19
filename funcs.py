# Author: Aaron Dunton (ajdunton@gmail.com)

import geopandas as gpd
import networkx as nx
import numpy as np
import math
import shapely
import os
from shapely.ops import split, unary_union, polygonize, nearest_points
from shapely.geometry import Point, LineString
from geovoronoi import voronoi_regions_from_coords
from statistics import mode
import igraph as ig
import leidenalg as la

import params as pr

def coord_dist(p0,p1):
    """
    Calculates the distance between two points.

    Parameters
    ----------
    p0 : Tuple
        Coordinate tuple of point 0.
    p1 : TYPE
        Coordinate tuple of point 1.

    Returns
    -------
    float
        Distance between two points.
    """
    return math.sqrt((p0[0] - p1[0])**2 + (p0[1] - p1[1])**2)

def split_in_half(l):
    """
    Function to split a shapely LineString in half.

    Parameters
    ----------
    l : LineString
        Line to split in half.

    Returns
    -------
    list
        List of two shapely LineStrings, original line split in half.

    """
    pt = l.interpolate(0.5,normalized=True)

    delx = l.coords[0][0]-l.coords[-1][0]
    dely = l.coords[0][1]-l.coords[-1][1]

    perpdelx = 1/(math.sqrt(delx**2+dely**2))*dely
    perpdely = -1/(math.sqrt(delx**2+dely**2))*delx

    pt1 = Point(pt.coords[0][0]+perpdelx*0.01*l.length, 
                pt.coords[0][1]+perpdely*0.01*l.length)
    pt2 = Point(pt.coords[0][0]-perpdelx*0.01*l.length, 
                pt.coords[0][1]-perpdely*0.01*l.length)
    breakline = LineString((pt1,pt2))

    return [x for x in split(l,breakline)]

def graph_from_shp(path):
    """
    Creates an undirected graph and position dictionary from a shapefile, with
    integer node numbers.
    
    Parameters
    ----------
    path : string
        Path to the shapefile.

    Returns
    -------
    X : networkx graph 
        Undirected graph corresponding to the input shapefile.
    pos : dictionary
        Position dictionary for the graph {node #: (coord_x,coord_y)}.
    """
    G = nx.read_shp(path).to_undirected()
    pos = {k: v for k,v in enumerate(G.nodes())}
    X = nx.Graph()
    X.add_nodes_from(pos.keys())
    co_edges = [set(x) for x in G.edges()]
    int_edges = [tuple(x for x,y in pos.items() if y in edg) 
                 for edg in co_edges]
    X.add_edges_from(int_edges)
    return X, pos

def digraph_from_shp(path):
    """
    Creates a directed graph and position dictioary from a shapefile, with 
    integer node numbers.
    
    There are two directed edges for each line in the shapefile.

    Parameters
    ----------
    path : string
        Path to the shapefile.

    Returns
    -------
    X : networkx digraph
        Directed graph corresponding to the input shapefile.
    pos : dictionary
        Position dictionary for the graph {node #: (coord_x,coord_y)}.

    """
    G = nx.read_shp(path)
    pos = {k: v for k,v in enumerate(G.nodes())}
    X = nx.DiGraph() 
    X.add_nodes_from(pos.keys())
    l = list(G.edges())
    edg = [tuple(k for k,v in pos.items() if v in sl) for sl in l]
    edg_copy = edg.copy()
    for x,y in edg_copy: edg.append((y,x))
    X.add_edges_from(edg)
    
    return X, pos

def street_cleaning():
    """
    Cleans the street network shapefile as inputs/streets/streets.shp and saves
    the resulting GeoDataFrame as streets_clean/streets_clean.shp, 4 steps:
       (1) Removes streets that are not LineStrings
       (2) Removes streets that begin and end at the same point
       (3) If two streets have the same endpoints, splits the longer
       (4) If the corresponding graph is not connected, remove streets 
           corresponding to the smaller component

    Returns
    -------
    coref : Coordinate reference system
        Coordinate reference system of the street shapefile.

    """
    # Open streets as a geodataframe and preliminary cleaning
    streets = gpd.read_file("inputs/streets/streets.shp")
    streets = streets.loc[streets['geometry'].is_valid, :]
    # Need to sort before dropping duplicates so that main=1 row stays when 
    # there is a duplicate with main=0
    streets = streets.sort_values('main',ascending=False)
    streets = streets.drop_duplicates(subset=['geometry'])
    streets = streets.reset_index(drop=True)[["main","known","geometry"]]
    
    # Remove streets that are not LineStrings
    notlinestrings = [i for i,geom in zip(streets.index,streets.geometry)
                      if type(geom)!=LineString]
    streets = streets.drop(index=notlinestrings).reset_index(drop=True)
    
    # Remove streets that begin and end at the same point (e.g., roundabouts)
    samepoints = [i for i,geom in zip(streets.index,streets.geometry)
                  if geom.coords[0]==geom.coords[-1]]
    streets = streets.drop(index=samepoints).reset_index(drop=True)
    
    # If two streets the same endpoints, split the longer
    multiples = []
    visited_rows = set()
    for ind_i, geom_i in zip(streets.index, streets.geometry):
        # Spatial query: get the streets that are intersect the current street,
        # with 1-foot buffer
        inter_streets = streets.loc[streets.intersects(geom_i.buffer(1))]
        visited_rows.add(ind_i)
        for ind_j, geom_j in zip(inter_streets.index,inter_streets.geometry):
            if ind_j not in visited_rows:
                # Calculate the separation between start points
                sep_0 = coord_dist(geom_i.coords[0], geom_j.coords[0])
                if sep_0 < 1:
                    # If start points are same, are end points the same?
                    sep_1 = coord_dist(geom_i.coords[-1], geom_j.coords[-1])
                    # Add to multiples if both start and end points are same
                    if sep_1 < 1:
                        multiples.append((ind_i,ind_j))         
                # Same computation, but case that start/end is reversed
                else:
                    sep_0 = coord_dist(geom_i.coords[0], geom_j.coords[-1])
                    if sep_0 < 1:
                        sep_1 = coord_dist(geom_i.coords[-1], geom_j.coords[0])
                        if sep_1 < 1:
                            multiples.append((ind_i,ind_j))    
                            
    # Remove streets that will be split and add split lines
    lines_and_main = [] #Save tuples of lines to split in half, main, known
    inds_to_drop = []
    for i in multiples:
        if streets.iloc[i[0]].geometry.length > \
        streets.iloc[i[1]].geometry.length:
            lines_and_main.append((streets.iloc[i[0]].geometry, 
                                   streets.iloc[i[0]].main,
                                   streets.iloc[i[0]].known))
            inds_to_drop.append(i[0])
        else:
            lines_and_main.append((streets.iloc[i[1]].geometry, 
                                   streets.iloc[i[1]].main,
                                   streets.iloc[i[1]].known))
            inds_to_drop.append(i[1])
    # Drop streets that will be split
    streets = streets.drop(index=inds_to_drop)
    # Add split streets back
    for geom,main_ind,known_ind in lines_and_main:
        try:
            new_rows = gpd.GeoDataFrame(data={'main':[main_ind,main_ind],
                                              'known':[known_ind,known_ind],
                                              'geometry':split_in_half(geom)})  
            streets = streets.append(new_rows)
        except ValueError:
            pass
    streets = streets.reset_index(drop=True)
    
    # Save cleaned streets to re-open with networkx
    if not os.path.isdir('inter'):
        os.mkdir('inter')
    streets.to_file('./inter/streets_clean')
    
    # Remove all but largest component if corresponding graph is not connected
    graph_check, pos_check = \
        graph_from_shp('./inter/streets_clean/streets_clean.shp')
    if not nx.is_connected(graph_check):
        nodes_to_delete = set(graph_check.nodes()) - \
            max(nx.connected_components(graph_check),key=len)
        s_inds = set()
        for i in nodes_to_delete:
            node_point = Point(pos_check[i]).buffer(1)
            for j in streets[streets.intersects(node_point)].index:
                s_inds.add(j)
        streets = streets.drop(index=s_inds).reset_index(drop=True)
        streets.to_file('./inter/streets_clean')  

    return streets

def identify_ind(pt,pos):
    """
    Identifies the graph index corresponding to a point geometry.

    Parameters
    ----------
    pt : point
        Point geometry of the target node.
    pos : dicitonary
        Position dictionary of the graph.

    Returns
    -------
    ind : int
        Index of the target node.
    """
    dists = {i:coord_dist(pt.coords[0],xy) for i,xy in pos.items()}
    ind = [i for i,dist in dists.items() if dist==min(dists.values())][0]
    return ind

def first_cut(lines,area):
    """
    Function that returns sub-areas given lines that might split the input area
    (i.e., they split the input area if they form a closed area, including the 
     boundary of the input area).

    Parameters
    ----------
    lines : GeoSeries
        Lines that may split the input area.
    area : GeoDataFrame
        Total area to be split into sub-areas according to any closed 
        boundaries indicated by input lines.

    Returns
    -------
    sub_areas : GeoDataFrame
        Resulting sub-areas, covering the entire input area.

    """
    # Get boundary from area and split at line intersections
    boundary = area.geometry[area.index[0]].buffer(-1).boundary
    split_list = [i for i in split(boundary,lines.unary_union)]
    split_gs = gpd.GeoSeries(split_list)
    lines = lines.clip(area.buffer(-1))
    lines_and_boundary = lines.append(split_gs).reset_index(drop=True)
    
    # Polygonize and create subareas geodataframe
    geom_list = [i for i in polygonize(lines_and_boundary)]
    sub_areas = gpd.GeoDataFrame(geometry=geom_list).set_crs(lines.crs)
    
    return sub_areas

def closest_line_ind(poly,lines):
    """
    Function that takes in a polygon geometry and a line gdf and returns 
    the index of the closest line to the polygon.

    Parameters
    ----------
    poly : shapely polygon
        Polygon from which we are determining the nearest line.
    lines : geopandas geodataframe
        Linestrings to which we are determining the nearest line.

    Returns
    -------
    line_ind : ing
        The index of the closest line
    """
    links = [None]*len(lines)
    for j, geom in enumerate(lines['geometry']):
        pt = nearest_points(poly,geom)
        links[j] = ((pt[0].x,pt[0].y),(pt[1].x,pt[1].y))
    zipped_list = zip([coord_dist(x,y) for x,y in links], list(lines.index))
    sorted_list = sorted(zipped_list)
    return sorted_list[0][1]

def points_on_line_and_endpoints(l,max_space):
    """
    Returns a list of points along a line that are equally spaced at no greater
    spacing that specified as input, including the endpoints.

    Parameters
    ----------
    l : LineString
        Line along which points are determined.
    max_space : float/int
        Maximum spacing of points along the line.

    Returns
    -------
    list
        List of Points.

    """
    # Choose number of points based on line length and max spacing, at least 1
    npts = math.ceil(l.length/max_space)+2
    pts = [l.interpolate((j)/(npts-1), normalized=True) for j in range(npts)]
    return pts


def refine_area(sts, total_area, graph, pos):
    """
    Function that breaks an area into a sub-areas based on connected components
    of a graph (i.e., the sub-areas contain connected graph components), where 
    the area is defined based on Voronoi polygons associated with points with a
    maximum spacing of 50 feet along the streets in each component.

    Parameters
    ----------
    sts : GeoDataFrame
        Streets that are in the area; the output sub-areas have connected 
        street networks.
    toal_area : Polygon
        Total area to be divided into sub-areas with connected street networks.
    graph : NetworkX Graph
        Graph to define the components of the street network in the area.
    pos : dictionary
        Position dictionary for the input graph.

    Returns
    -------
    new_areas : List of Polygons
        List of new areas, each with a connected street network inside.

    """
    try: # First try with voronoi points at 50 foot separation    
        # List of tuples: (Point, Component ID)
        points_and_comp_id = []
        
        # Loop over components
        for comp_i, g_inds in enumerate(nx.connected_components(graph)):
            # List and geodataframe of node points
            node_pts_list = [Point(pos[i]) for i in g_inds]
            node_pts = \
            gpd.GeoDataFrame(geometry=node_pts_list).set_crs(sts.crs)
            
            # Get streets for this component based on the node points gdf
            comp_streets = \
            sts.loc[sts.intersects(node_pts.buffer(1).unary_union)]
            
            # Loop over the streets in this component                       
            for street in comp_streets.geometry:
                # Get points on each street with maximum spacing of 50 feet                    
                for pt in points_on_line_and_endpoints(street,50):
                    # Add each point and the current component id to the list
                    points_and_comp_id.append((pt,comp_i))
        
        # Unpack the data
        points_data = {'comp_ind':[ind for _,ind in points_and_comp_id],
                       'geometry':[pt for pt,_ in points_and_comp_id]}
        
        # Create a gdf of points from which the voronoi polygons will be 
        # determined
        vor_pts_gdf = gpd.GeoDataFrame(points_data).set_crs(sts.crs)
        
        # Voronoi polygons and unpack results, save as column voroni_poly 
        polys, points = voronoi_regions_from_coords(list(vor_pts_gdf.geometry),
                                                    total_area)
    
    except RuntimeError: #try again with 100 if there is some numerical issue
        # List of tuples: (Point, Component ID)
        points_and_comp_id = []
        
        # Loop over components
        for comp_i, g_inds in enumerate(nx.connected_components(graph)):
            # List and geodataframe of node points
            node_pts_list = [Point(pos[i]) for i in g_inds]
            node_pts = \
            gpd.GeoDataFrame(geometry=node_pts_list).set_crs(sts.crs)
            
            # Get streets for this component based on the node points gdf
            comp_streets = \
            sts.loc[sts.intersects(node_pts.buffer(1).unary_union)]
            
            # Loop over the streets in this component                       
            for street in comp_streets.geometry:
                # Get points on each street with maximum spacing of 50 feet                    
                for pt in points_on_line_and_endpoints(street,100):
                    # Add each point and the current component id to the list
                    points_and_comp_id.append((pt,comp_i))
        
        # Unpack the data
        points_data = {'comp_ind':[ind for _,ind in points_and_comp_id],
                       'geometry':[pt for pt,_ in points_and_comp_id]}
        
        # Create a gdf of points from which the voronoi polygons will be 
        # determined
        vor_pts_gdf = gpd.GeoDataFrame(points_data).set_crs(sts.crs)
        
        # Voronoi polygons and unpack results, save as column voroni_poly
        polys, points = voronoi_regions_from_coords(list(vor_pts_gdf.geometry),
                                                    total_area)        
    
    geom_list = [None]*len(vor_pts_gdf)
    for x,y in polys.items():
        geom_list[points[x][0]] = y
    vor_pts_gdf['voroni_poly'] = geom_list
    vor_pts_gdf.dropna(subset=['voroni_poly'],inplace=True)
    
    # Create list of new areas by looping over components and unionizing voroni
    # polygons
    new_areas = []
    for i in set(vor_pts_gdf.comp_ind):
        pts_in_comp = vor_pts_gdf.loc[vor_pts_gdf.comp_ind==i]
        new_areas.append(unary_union(pts_in_comp.voroni_poly))
    
    return new_areas

def refine_area_bld(blds, sts, tot_a_cur, graph, pos):
    """
    Function that breaks an area into a sub-areas based on connected components
    of a graph (i.e., the sub-areas contain connected graph components), where 
    the area is defined based on Voronoi polygons on buildings.

    Parameters
    ----------
    blds : GeoDataFrame
        Buildings that are in the area, used to define the voronoi polygons 
        that define the output geometry.
    sts : GeoDataFrame
        Streets that are in the area; the output sub-areas have connected 
        street networks.
    tot_a_cur : Polygon
        Total area to be divided into sub-areas with connected street networks.
    graph : NetworkX Graph
        Graph to define the components of the street network in the area.
    pos : dictionary
        Position dictionary for the input graph.

    Returns
    -------
    new_areas : List of Polygons
        List of new areas, each with a connected street network inside.

    """
    # Index of the closest street to each building
    street_index_for_bld = [None]*len(blds)   
    for i,bld_geom in enumerate(blds['geometry']):
        street_index_for_bld[i] = closest_line_ind(bld_geom,sts)
    blds['street_index'] = street_index_for_bld
    
    # Create Voronoi polygon associated with each building
    polys, points = voronoi_regions_from_coords(blds.centroid, tot_a_cur)
    geom_list = [None]*len(blds)
    for x,y in polys.items():
        geom_list[points[x][0]] = y
    blds['voroni_poly'] = geom_list
    
    # Output list of new areas
    new_areas = []
    
    # Go through each connected component of the graph
    for inds in nx.connected_components(graph):
        # Identify gdf of streets for each component
        pts_list = [Point(pos[i]) for i in inds]
        pts = gpd.GeoDataFrame(geometry=pts_list).set_crs(blds.crs)
        s_new = sts.loc[sts.intersects(pts.buffer(1).unary_union)]
        
        # Identify gdf of buildings for each component
        b_new = blds.loc[blds.street_index.isin(s_new.index)]

        # Combine Voronoi polygons associated with buildings in the new area
        new_areas.append(unary_union(b_new.voroni_poly))
    
    return new_areas

def edges_for_gdf(lines,pos):
    """
    Identifies corresponding edge tuples for a GeoSeries of LineStrings.

    Parameters
    ----------
    lines : GeoSeries
        GeoSeries of lines, each of which will be identified as an edge.
    pos : Dictionary
        Position dictionary of the corresponding graph.

    Returns
    -------
    list_of_edges : List
        List of edge tuples.

    """
    i_of_edges = np.empty(len(lines),dtype=int)
    j_of_edges = np.empty(len(lines),dtype=int)
    for ind,l in enumerate(lines):
        i_coords = l.boundary[0].coords[0]
        j_coords = l.boundary[1].coords[0]
        
        for node_ind,node_coords in pos.items():
            dist = coord_dist(i_coords,node_coords)
            if dist<1:
                break
        i_of_edges[ind] = node_ind
        
        for node_ind,node_coords in pos.items():
            dist = coord_dist(j_coords,node_coords)
            if dist<1:
                break
        j_of_edges[ind] = node_ind
        
    return [(i,j) for i,j in zip(i_of_edges,j_of_edges)]

def sample_line(l,ras):
    """
    Funciton to get the profile of a raster along a line, with the number of 
    sample points based on cell size and line length.

    Parameters
    ----------
    l : LineString geometry (or GeoSeries)
        Line along which raster should be sampled.
    ras : rasterio.io.DatasetReader
        Raster to sample.

    Returns
    -------
    profile : List
        Values along the raster.

    """
    if isinstance(l,gpd.geoseries.GeoSeries):
        min_size = min((abs(ras.transform[0]),abs(ras.transform[4]))) 
        npts = max(math.ceil(2*l.length[l.index[0]]/min_size), 2) 
        samplepts = [l.interpolate(j/(npts-1), normalized=True) 
                     for j in range(npts)]
        profile = [[i for i in ras.sample([(j[l.index[0]].x, 
                                            j[l.index[0]].y)])][0][0] 
                   for j in samplepts]
   
    if isinstance(l,shapely.geometry.linestring.LineString):
        min_size = min((abs(ras.transform[0]),abs(ras.transform[4])))
        npts = max(math.ceil(2*l.length/min_size), 2)
        samplepts = [l.interpolate(j/(npts-1), normalized=True) 
                     for j in range(npts)]
        profile = [[i for i in ras.sample([(j.x, j.y)])][0][0] 
                   for j in samplepts]

    return profile

def eval_weight(row,dem):
    """
    Function that evaluates the weights for the street edges.

    Parameters
    ----------
    row : Row
        Row of the streets GeoDataFrame, from .apply(axis=1).
    dem : rasterio.io.DatasetReader
        Digital elevation model.

    Returns
    -------
    wij : Float/Int
        Weight of the edge in the forward direction (as indicated by the order 
        of the points in the linestring).
    wji : Float/Int
        Weight of the edge in the reverse direction (as indicated by the order 
        of the points in the linestring).

    """
    prof = sample_line(row.geometry,dem)
    h_max = max(prof)
    hi = prof[0]
    hj = prof[-1]
    
    length = row.geometry.length
    
    if (hi-hj)/length > -pr.scr:
        wij = min(0.5**(((hj-hi)+(h_max-max(hi,hj)))/pr.hpt5), 5)
    else: 
        wij = 0
        
    if (hj-hi)/length > -pr.scr:
        wji = min(0.5**(((hi-hj)+(h_max-max(hi,hj)))/pr.hpt5), 5)
    else:
        wji = 0
        
    return wij, wji

def leiden(graph):
    """
    Runs the Leiden algorithm on the input graph, considering edge weights, and 
    does some post-processing of the results.

    Parameters
    ----------
    graph : networkx graph
        Graph with weights as edge attributes.

    Returns
    -------
    node_assignments : dictionary
        Community for each node {Node # : Community #}.

    """
    ig_graph = ig.Graph.from_networkx(graph)

    partition = la.ModularityVertexPartition(ig_graph, weights='weight')
    opt = la.Optimiser()
    opt.consider_comms = la.RAND_NEIGH_COMM
    opt.refine_consider_comms = la.RAND_NEIGH_COMM
    _ = opt.optimise_partition(partition)
    
    # Assigned communities for each node {node number: community}
    # Assumes that partition.membership in same order as streetgraph.nodes()
    node_assignments = dict(zip(graph.nodes(), partition.membership))
    
    # Post-processing: assign to mode among neighbor's communities if none of 
    # neighbors are assigned to same community
    for i in graph.nodes():    
        adj_coms = []
        for j in graph.neighbors(i):
            adj_coms.append(node_assignments[j])

        new_com = None
        if all(com!=node_assignments[i] for com in adj_coms):
            new_com = mode(adj_coms)
            print(f'Node {i} moved from {node_assignments[i]} to {new_com}.')
            node_assignments[i] = new_com
            
    return node_assignments

def points_on_line(l,max_space):
    """
    Returns a list of points along a line that are equally spaced at no greater
    spacing that specified as input.

    Parameters
    ----------
    l : LineString
        Line along which points are determined.
    max_space : float/int
        Maximum spacing of points along the line.

    Returns
    -------
    list
        List of Points.

    """
    # Choose number of points based on line length and max spacing, at least 1
    npts = max(math.ceil(l.length/max_space), 1)
    pts = [l.interpolate((j+1)/(npts+1), normalized=True) for j in range(npts)]
    return pts

def get_subareas(nodes, streets):
    """
    Returns a GeoDataFrame of sub-areas based on community assignments in 
    nodes, also adding intermediate nodes when getting the voroni polygons to
    smooth out the borders.

    Parameters
    ----------
    nodes : GeoDataFrame
        Nodes GDF with 'com_id', 'geometry', and 'node_id' for each row/node.
    streets : GeoDataFrame
        Streets GDF with 'graph_edge' and 'geometry' for each row/street.

    Returns
    -------
    subareas_gdf : GeoDataFrame
        GDF of the sub-areas as determined by voronoi polygons, also including 
        the community id and associated node ids for each row/subarea.

    """
    # Create copy of nodes to which extra nodes will be added
    vor_pts_gdf = nodes[['com_id','geometry']].copy()
    
    # Add new points if community of two adjacent nodes is the same
    new_rows = []
    for inds,geom in zip(streets['graph_edge'],streets['geometry']):
        com_i = int(nodes.loc[nodes['node_id']==inds[0]].com_id)
        com_j = int(nodes.loc[nodes['node_id']==inds[1]].com_id)
        
        if  com_i == com_j:
            for k in points_on_line(geom,50):
                new_rows.append({'com_id':com_i, 'geometry':k})
                
    vor_pts_gdf = vor_pts_gdf.append(gpd.GeoDataFrame(new_rows))
                
    # Open the bounding box, total area that all Voronoi polygons will cover
    bound_box = gpd.read_file("inputs/box_bounds/box_bounds.shp")


    # Voronoi polygons for every node
    polys, points = voronoi_regions_from_coords(list(vor_pts_gdf["geometry"]),
                                                bound_box["geometry"][0])
    points_reversed = {y[0]:x for x,y in points.items()}
    for vor_id, input_id in points.items():
        if len(input_id) > 1:
            for x in input_id[1:]: points_reversed[x] = vor_id
    vor_pts_gdf['vor'] = [polys[points_reversed[i]] 
                          for i in range(len(vor_pts_gdf))]

    # Combine polygons that belong to the same sub-area and create subareas gdf
    subs = {i:unary_union(vor_pts_gdf.loc[vor_pts_gdf['com_id']==i]['vor']) 
            for i in set(nodes['com_id'])}
    
    node_ids = [None]*len(subs)
    for i,j in enumerate(subs.keys()):
        node_ids[i] = [ind for ind,com in zip(nodes['node_id'],
                                              nodes['com_id'])
                       if com==j]
    
    subareas_gdf = gpd.GeoDataFrame({'com_id':subs.keys(),
                                     'node_ids':node_ids,
                                     'geometry':subs.values()})
    return subareas_gdf

def identify_outlet(row, sts, graph, pos):
    """
    Identifies outlet for each subarea using the specified method.

    Parameters
    ----------
    row : Row
        Row of subareas GeoDataFrame, from .apply(axis=1).
    sts : GeoDataFrame
        GeoDataFrame of the streets.
    graph : NetworkX graph
        Mainline street network graph, with attributes.
    pos : dictionary
        Mainline street network position dictionary.

    Returns
    -------
    outlet : Point
        Outlet location for each row/subarea.
    """
    streets_in_area = sts.clip(row.geometry, keep_geom_type=True)
    non_main = streets_in_area.loc[streets_in_area.main==0]
    all_mainlines = sts.loc[sts.main==1]
    
    candidates = non_main.unary_union.intersection(all_mainlines.unary_union)
    if isinstance(candidates,Point):
        outlet = candidates
    else:
        values = []
        for pt in candidates:
            pt_ind = identify_ind(pt,pos)
            if pr.outlet_method=='elev':
                values.append(graph.nodes[pt_ind]['elev'])
            elif pr.outlet_method=='wwtp':
                values.append(graph.nodes[pt_ind]['path_weight_to_wwtp'])
        outlet = [pt for pt,val in zip(candidates,values) 
                  if val==min(values)][0]
    return outlet

def main_graph(mains, dem):
    """
    Function to get full graph of mainline streets; samples DEM at each node 
    and stores as node attribute 'elev'; evaluates the weight of each edge 
    using the specified method and stores as edge attribute 'weight'; and 
    determines the shortest path distance from each node to the treatment plant
    node and stores as node attribute 'path_weight_to_wwtp'.

    Parameters
    ----------
    mains : GeoDataFrame
        GDF of main streets.
    dem : raster
        Digital elevation model.

    Returns
    -------
    graph : NetworkX DiGraph
        Full street mainline graph; node attributes 'elev' and 
        'path_weight_to_wwtp'; edge attribute 'weight'.
    pos : dictionary
        Position dictionary of the nodes in graph.
    term_id : int
        Node index for the treatment plant (i.e., terminal) node.
    """
    mains.to_file('./inter/main_streets')
    graph, pos = digraph_from_shp('./inter/main_streets/main_streets.shp')
    
    # Sample elevation at nodes
    elevations = dict.fromkeys(graph.nodes())
    for i in graph.nodes():
        graph.nodes[i]['elev'] = next(dem.sample([pos[i]]))[0]
        elevations[i] = graph.nodes[i]['elev']
    
    # Evaluate weight of edges
    if pr.ls_weight_method == 'distance':
        for i,j in graph.edges(): 
            lij = coord_dist(pos[i],pos[j])
            graph[i][j]['weight'] = lij
            
    elif pr.ls_weight_method == 'hybrid-elev':
        max_elev = max(elevations.values())
        min_elev = min(elevations.values())
        for i,j in graph.edges():
            lij = coord_dist(pos[i],pos[j])
            hj = elevations[j]
            graph[i][j]['weight'] = (hj-min_elev)/(max_elev-min_elev)*lij
            
    known_mains = mains[mains.known==1]
    for line in known_mains.geometry:
        pt_i,pt_j = line.boundary
        i = identify_ind(pt_i,pos)
        j = identify_ind(pt_j,pos)
        graph[i][j]['weight'] = pr.known_weight_ratio*graph[i][j]['weight']
        graph[j][i]['weight'] = pr.known_weight_ratio*graph[j][i]['weight']
            
    # For each node, evaluate minimum distance to the treatment plant node
    wwtp_loc = gpd.read_file('inputs/wwtp_loc/wwtp_loc.shp')
    term_id = identify_ind(wwtp_loc.geometry[wwtp_loc.index[0]],pos)
    dist, path = nx.single_source_dijkstra(graph,term_id)
    for i in graph.nodes():
        graph.nodes[i]['path_weight_to_wwtp'] = dist[i]
    
    return graph, pos, term_id

def max_min(val, max_val, min_val):
    """
    Returns value scaled using min-max scaling.

    Parameters
    ----------
    val : float
        Value to be scaled.
    max_val : float
        Maximum value in data (i.e., max_min(max_val) = 1).
    min_val : float
        Minimum value in data (i.e., max_min(min_val) = 0).

    Returns
    -------
    float
        Scaled value, between 0 and 1.

    """
    return (val-min_val)/(max_val-min_val)

def main_spst(graph,target,needed_nodes,keepweights=False):
    """
    Determines the shortest path spanning tree of an input graph considering a 
    target, only requiring to connect to needed_nodes, using Dijkstra's 
    algorithm.

    Parameters
    ----------
    graph : networkx graph
        Input graph.
    target : int
        Node number corresponding to the target node.
    needed_nodes : list
        List of node indices that should be included in the tree (i.e., other 
        nodes may or may not be included in the tree, depending on if they are 
        in one of the shortest paths to connect the needed_nodes). 
    keepweights : boolean
        If true, store weights from the input graph in the output tree 
        (default False).

    Returns
    -------
    t : networkx directed graph
        Shortest path spanning tree, where edge direction indicates the path 
        from each node to the target node.
    """
    dist, path = nx.single_source_dijkstra(graph,target)
    node_set = set()
    edge_set = set()
    for n_id in needed_nodes:
        path_from_n_to_tn = path[n_id].copy()
        path_from_n_to_tn.reverse()
        if len(path_from_n_to_tn) > 1:
            for i,j in zip(path_from_n_to_tn,path_from_n_to_tn[1:]):
                edge_set.add((i,j))
                node_set.add(i)
                node_set.add(j)
    t = nx.DiGraph()
    t.add_nodes_from(node_set)
    t.add_edges_from(edge_set)
    if keepweights:
        for i,j in t.edges():
            t[i][j]['weight'] = graph[i][j]['weight']
    return t

def assign_to_cell(row,lspop):
    """
    Assigns each a building to a cell of the vectorized landscan population 
    data.

    Parameters
    ----------
    row : row
        Row of the buildings GeoDataFrame, from .apply(axis=1).
    lspop : GeoDataFrame
        Vectorized landscan population data.

    Returns
    -------
    cell : int
        Index of the corresponding lspop cell/row for each building/row.

    """
    cell_list = [ind 
            for ind,geom 
            in zip(lspop.index,lspop.geometry) 
            if row.geometry.centroid.within(geom)]
    if len(cell_list)>0: cell = int(cell_list[0])
    else: cell = None
    return cell

def tot_bld_area(row, buildings):
    """
    Calculates the total building area associated with each cell/row of the 
    vectorized population data.

    Parameters
    ----------
    row : row
        Row of the population geodataframe, from .apply(axis=1).
    buildings : GeoDataFrame
        Buildings GeoDataFrame, with corresponding cell/row index of the 
        population GDF as lspop_id.

    Returns
    -------
    float
        Sum of the area of all of the building fooprints assigned to each 
        cell/row.

    """
    return sum([geom.area 
                for lspop_id, geom 
                in zip(buildings.lspop_id, buildings.geometry) 
                if lspop_id==row.name])

def population_column(row, lspop):
    """
    Calculates the poplation assigned to each building.

    Parameters
    ----------
    row : row
        Row of the buildings GeoDataFrame, from .apply(axis=1).
    lspop : GeoDataFrame
        Vectorized landscan population data, with total building area for each 
        cell/row = tot_bld_area, and population for each cell/row = pop.

    Returns
    -------
    TYPE
        DESCRIPTION.

    """
    a_ind = row.geometry.area
    tot_area = lspop.iloc[int(row.lspop_id)].tot_bld_area
    tot_pop = lspop.iloc[int(row.lspop_id)]['pop']
    return tot_pop*a_ind/tot_area

def detect_boundary(seg,area,outlet):
    """
    Detects whether each line in the clipped streets GDF intersects the 
    boundary, and which of these intersects the outlet.

    Parameters
    ----------
    seg : row 
        Row of the streets GeoDataFrame, from .apply(axis=1).
    area : polygon
        Area for corresponding sub-area.
    outlet : point
        Outlet location for corresponding sub-area.

    Returns
    -------
    result : string
        String indicating if each row/street segment is 'on boundary, outlet'; 
        'on boundary, no outlet'; or 'not on boundary'.

    """
    if seg.geometry.intersects(area.boundary.buffer(1)):
        if seg.geometry.intersects(outlet.buffer(1)):
            result = "on boundary, outlet"
        else: result = "on boundary, no outlet"
    else: result = "not on boundary"
    return result