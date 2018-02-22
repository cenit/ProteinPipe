# -*- coding: utf-8 -*-
"""
Created on Thu Feb  1 15:31:12 2018

@author: nico
"""

import networkx as nx # graph and layout
import pandas as pd # edges dataframe 
import sys # exit
import os # exists file
import argparse # parse command line
import numpy as np # append numpy array
import scipy # squareform, pdist, csgraph, eigsh
import itertools # combination
import scipy.stats
# blender python
import bpy 
from math import acos, degrees, pi
from mathutils import Vector
from matplotlib import colors as mcolors

ALL_C = dict(mcolors.BASE_COLORS, **mcolors.CSS4_COLORS)
ALL_C = {mcolors.to_rgb(ash) : name for name, ash in ALL_C.items()}

#sns.color_palette("hls", 20)
all_colors = {
    "ALA" : (0.85999999999999999, 0.37119999999999997, 0.33999999999999997),
    "ARG" : (0.85999999999999999, 0.5272, 0.33999999999999997),
    "ASN" : (0.85999999999999999, 0.68320000000000003, 0.33999999999999997),
    "ASP" : (0.85999999999999999, 0.83920000000000017, 0.33999999999999997),
    "CYS" : (0.72479999999999989, 0.85999999999999999, 0.33999999999999997),
    "GLN" : (0.56880000000000008, 0.85999999999999999, 0.33999999999999997),
    "GLU" : (0.41279999999999994, 0.85999999999999999, 0.33999999999999997),
    "GLY" : (0.33999999999999997, 0.85999999999999999, 0.42320000000000013),
    "HIS" : (0.33999999999999997, 0.85999999999999999, 0.57920000000000016),
    "ILE" : (0.33999999999999997, 0.85999999999999999, 0.73520000000000008),
    "LEU" : (0.33999999999999997, 0.82879999999999987, 0.85999999999999999),
    "LYS" : (0.33999999999999997, 0.67279999999999973, 0.85999999999999999),
    "MET" : (0.33999999999999997, 0.51679999999999948, 0.85999999999999999),
    "PHE" : (0.33999999999999997, 0.36079999999999973, 0.85999999999999999),
    "PRO" : (0.47520000000000029, 0.33999999999999997, 0.85999999999999999),
    "SER" : (0.63119999999999976, 0.33999999999999997, 0.85999999999999999),
    "THR" : (0.7871999999999999, 0.33999999999999997, 0.85999999999999999),
    "TRP" : (0.85999999999999999, 0.33999999999999997, 0.77679999999999927),
    "TYR" : (0.85999999999999999, 0.33999999999999997, 0.62079999999999991),
    "VAL" : (0.85999999999999999, 0.33999999999999997, 0.46479999999999977)
}


def blend_rec(graph, tcoords, gcoords, colors, label, node_size=3, edge_thickness=0.25):
    # edge obj
    bpy.data.materials.new(name="light_gray")
    bpy.data.materials["light_gray"].diffuse_color = (0.4627450980392157, 0.4627450980392157, 0.4627450980392157)
    bpy.data.materials["light_gray"].specular_intensity = 0.5
    # sphere obj
    for color, aa in zip(colors, label):
        bpy.data.materials.new(name=aa)
        bpy.data.materials[aa].diffuse_color = color
        bpy.data.materials[aa].specular_intensity = 0.5

        # Transparency parameters
        bpy.data.materials[aa].use_transparency = True
        bpy.data.materials[aa].transparency_method = "RAYTRACE"
        bpy.data.materials[aa].alpha = 0.95
        bpy.data.materials[aa].raytrace_transparency.fresnel = 0.1
        bpy.data.materials[aa].raytrace_transparency.ior = 1.15

    # sphere rec obj
    bpy.data.materials.new(name="red")
    bpy.data.materials["red"].diffuse_color = (1, 0, 0)
    bpy.data.materials["red"].specular_intensity = 0.5

    # text obj
    text = bpy.data.objects.new("label", bpy.data.curves.new(type="FONT", name="curve"))
    bpy.data.materials.new(name="black")
    bpy.data.materials["black"].diffuse_color = (0, 0, 0)
    bpy.data.materials["black"].specular_intensity = 0.5
    
    # Set scene, light and alpha_mode
    scene = bpy.context.scene
    scene.render.engine = 'BLENDER_RENDER' # 'CYCLE'
    scene.render.alpha_mode = 'TRANSPARENT' # remove background

    area = next(area for area in bpy.context.screen.areas if area.type == 'VIEW_3D')
    area.spaces[0].region_3d.view_perspective = 'CAMERA'

    bpy.data.worlds["World"].light_settings.use_ambient_occlusion = True
    bpy.data.worlds["World"].light_settings.samples = 10

    # camera position
    if(len(bpy.data.cameras) == 1):
        obj = bpy.data.objects['Camera'] # bpy.types.Camera

        mean_x = np.mean([el[0] for el in tcoords.values()])
        mean_y = np.mean([el[1] for el in tcoords.values()])
        mean_z = np.mean([el[2] for el in tcoords.values()])

        obj.location.x = mean_x
        obj.location.y = mean_y
        min_x = min([el[0] for el in tcoords.values()])
        max_x = max([el[0] for el in tcoords.values()])
        min_y = min([el[1] for el in tcoords.values()])
        max_y = max([el[1] for el in tcoords.values()])
        obj.location.z = mean_z + max(max_x - min_x, max_y - min_y)/(2 * np.tan(np.radians(20)/2))
        obj.rotation_euler = (0.0, 0.0, np.radians(90))
        obj.keyframe_insert(data_path="location", frame=10.0)

    # Add some mesh primitives
    bpy.ops.object.select_all(action='DESELECT')
    bpy.ops.mesh.primitive_uv_sphere_add()
    sphere = bpy.context.object
    bpy.ops.mesh.primitive_cylinder_add()
    cylinder = bpy.context.object
    cylinder.active_material = bpy.data.materials["light_gray"]

    bpy.ops.mesh.primitive_uv_sphere_add()
    gsphere = bpy.context.object
    gsphere.active_material = bpy.data.materials["red"]
    
    # Keep references to all nodes and edges
    shapes = []
    # Keep separate references to shapes to be smoothed
    shapes_to_smooth = []
    
    # Draw nodes
    nx.set_node_attributes(graph, [], "color")
    for node, aa in zip(graph.nodes(), label):
        # Coloring rule for nodes. Edit this to suit your needs!
        col = colors[node]
        
        # Copy mesh primitive and edit to make node
        # (You can change the shape of drawn nodes here)
        node_sphere = sphere.copy()
        node_sphere.data = sphere.data.copy()
        node_sphere.location = tcoords[node]
        node_sphere.dimensions = [node_size] * 3
        node_sphere.active_material = bpy.data.materials[aa]
        bpy.context.scene.objects.link(node_sphere)
        shapes.append(node_sphere)
        shapes_to_smooth.append(node_sphere)
        
        #lbl = text.copy()
        #lbl.data = text.data.copy()
        #lbl.data.body = label[node]
        #lbl.rotation_mode = "AXIS_ANGLE"
        #lbl.rotation_euler = (0.0, 0.0, np.radians(90))
        #lbl.active_material = bpy.data.materials["black"]
        #lbl.dimensions = [.00000001] * 3
        #lbl.location = np.asarray(tcoords[node]) + [0., 0., node_size/2]
        #bpy.context.scene.objects.link(lbl)

        # guess atom
        gues_sphere = gsphere.copy()
        gues_sphere.data = gsphere.data.copy()
        gues_sphere.dimensions = [node_size] * 3
        gues_sphere.location = gcoords[node]
        gues_sphere.active_material = bpy.data.materials["red"]
        bpy.context.scene.objects.link(gues_sphere)
        shapes.append(gues_sphere)
        shapes_to_smooth.append(gues_sphere)


    # Draw edges
    for source, target in graph.edges():
        # Get source and target locations by drilling down into data structure
        source_loc = tcoords[source]
        target_loc = tcoords[target]

        diff = [c2 - c1 for c2, c1 in zip(source_loc, target_loc)]
        cent = [(c2 + c1) / 2 for c2, c1 in zip(source_loc, target_loc)]
        mag = sum([(c2 - c1) ** 2 for c1, c2 in zip(source_loc, target_loc)]) ** 0.5
        
        # Euler rotation calculation
        v_axis = Vector(diff).normalized()
        v_obj = Vector((0, 0, 1))
        v_rot = v_obj.cross(v_axis)
        angle = acos(v_obj.dot(v_axis))

        # Copy mesh primitive to create edge
        edge_cylinder = cylinder.copy()
        edge_cylinder.data = cylinder.data.copy()
        edge_cylinder.dimensions = [edge_thickness] * 2 + [mag - node_size]
        edge_cylinder.location = cent
        edge_cylinder.rotation_mode = "AXIS_ANGLE"
        edge_cylinder.rotation_axis_angle = [angle] + list(v_rot)
        bpy.context.scene.objects.link(edge_cylinder)
        shapes.append(edge_cylinder)
        shapes_to_smooth.append(edge_cylinder)
        
    # Remove primitive meshes
    bpy.ops.object.select_all(action='DESELECT')
    sphere.select = True
    gsphere.select = True
    cylinder.select = True
    
    # If the starting cube is there, remove it
    if "Cube" in bpy.data.objects.keys():
        bpy.data.objects.get("Cube").select = True
    bpy.ops.object.delete()

    # Smooth specified shapes
    for shape in shapes_to_smooth:
        shape.select = True
    bpy.context.scene.objects.active = shapes_to_smooth[0]
    bpy.ops.object.shade_smooth()

    # Join shapes
    for shape in shapes:
        shape.select = True
    bpy.context.scene.objects.active = shapes[0]
    bpy.ops.object.join()

    # Center object origin to geometry
    bpy.ops.object.origin_set(type="ORIGIN_GEOMETRY", center="MEDIAN")

    # Refresh scene
    bpy.context.scene.update()

if __name__ == "__main__":
    argv = sys.argv
    argv = argv[argv.index("--") + 1:]  # get all args after "--"
    
    description = "BlenderPDBRec"
    parser = argparse.ArgumentParser(description = description)
    parser.add_argument("-T", required=True, dest="tcoord", action="store", help="True Protein coordinates filename (.xyz)", default="")
    parser.add_argument("-G", required=True, dest="gcoord", action="store", help="Guessed Protein coordinates filename (.guess)", default="")
    parser.add_argument("-s", required=False, dest="nsize", action="store", help="node size", default=3)
    parser.add_argument("-l", required=False, dest="esize", action="store", help="edge thickness", default=.25)
    parser.add_argument("-t", required=False, dest="thr", action="store", help="cmap threashold", default=8)
    
    if len(sys.argv) <= 1:
        parser.print_help()
        sys.exit(1)
    else:
        args = parser.parse_args(argv)
    
    tfile = args.tcoord
    gfile = args.gcoord
    node_size = float(args.nsize)
    edge_thickness = float(args.esize)
    thr = float(args.thr)

    true_protein = pd.read_csv(tfile, sep="\t", header=None, names=["atoms", "x", "y", "z"])
    gues_protein = pd.read_csv(gfile, sep="\t", header=None, names=["atoms", "x", "y", "z"])

    cmap = nx.from_numpy_matrix( np.asarray(scipy.spatial.distance.squareform(scipy.spatial.distance.pdist(true_protein.iloc[:,1:], metric="euclidean") < thr), dtype=np.float) )

    # Translation
    true_protein.iloc[:, 1:] -= true_protein.iloc[:, 1:].mean(axis = 0)
    gues_protein.iloc[:, 1:] -= gues_protein.iloc[:, 1:].mean(axis = 0)

    # Scaling
    true_protein.iloc[:, 1:] /= np.linalg.norm(true_protein.iloc[:, 1:], axis=0)
    gues_protein.iloc[:, 1:] /= np.linalg.norm(gues_protein.iloc[:, 1:], axis=0)

    # find right permutation of axis 
    combo = list(itertools.permutations(["x", "y", "z"]))
    permutation = list(combo[ np.argmax( [ sum([ scipy.stats.pearsonr(x=true_protein.iloc[:, 1:][true], y=gues_protein.iloc[:, 1:][guess])[0] 
                                           for true, guess in zip(["x","y","z"], list(comb))
                                          ]) 
                                        for comb in combo ]
                                        )])

    tot_dist =  np.sum( np.sqrt( np.sum( (true_protein.iloc[:, 1:] - gues_protein.iloc[:, 1:][permutation])**2, axis=1) ) )
    print("Total distance : ", tot_dist )
    print("RMSD : ", tot_dist / len(true_protein) )

    colors = [all_colors[aa] for aa in list(true_protein.atoms)]
    names = list(true_protein.atoms)
    true_position = {i : (x, y, z) for i, (x,y,z) in true_protein[["x", "y", "z"]].iterrows()}
    gues_position = {i : (x, y, z) for i, (x,y,z) in gues_protein[permutation].iterrows()}

    blend_rec(graph = cmap, tcoords = true_position, gcoords = gues_position, colors = colors, label = names, node_size = node_size, edge_thickness = edge_thickness )

