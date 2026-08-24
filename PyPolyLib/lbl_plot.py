"""
lbl_plot.py — Visualizing Z-polyhedra with matplotlib.
Typical usage:
    from pypolylib import LBL
    a = LBL("{(i, j) | 1 <= i <= 10, 1 <= j <= 10}")
    a.plot()
"""
import sys
import pypolylib_core as pl
import pypolylib

import math
import numpy as np

# plotting libs:
from scipy.spatial import ConvexHull
import colorsys
import pyvista as pv

# # 2D -> not used anymore
# import matplotlib.pyplot as plt
# from shapely.geometry import MultiPoint

# this variable is used to artificially bound an unbounded lbl being plotted,
# it is the ray/line multiplier to get an upper/lower bound
ZOOM = 4

# allow pyvista to plot polyhedra that have no integer points:
pv.global_theme.allow_empty_mesh = True


# this class is used to store a list of LBLs to be plotted.
# Contains everything that's needed to plot the right object
# (bounded or not, bounding box, etc.)
class LBL_view:
    lbl = None          # the LBL it comes from
    color = None
    show_points:bool
    lat:pl.Matrix
    poly_coordinate:pl.Polyhedron
    poly_convex_hull:pl.Polyhedron
    is_bounded:bool
    # build python lists of FP vertices/rays/lines of hull once:
    vrl:tuple


# my class to plot single LBLs with pyvista's Plotter and set object colors
# inherits from pyvista's Plotter class:
class MyWindow(pv.Plotter):
    hue:float
    _3D:bool
    LBLs:list           # list of LBL_view's
    all_bounded:bool    # they are all bounded
    my_origin: list

    def __init__(self):
        self.hue = 0.0
        self._3D = False
        self.LBLs = []
        self.all_bounded = True
        self.my_origin = []
        super().__init__()

    # round robbin hue, this is a pretty nice generator of various colors
    def color(self):
        self.hue = (self.hue + 0.618033988749895) % 1.0
        return colorsys.hsv_to_rgb(self.hue, .55, .9)


# this is not very nice, but it's the best way to do it.
# _lbl_plot_window[-1] is set to a MyWindow instance, and a new window will
# be added to the list, unless subplot = True.
_lbl_plot_window = [MyWindow()]


# main plotting function
def lbl_plot(
        lbl, *args,
        show_points=True, subplot=False, show_grid="custom", color=None,
        **kwargs
    ):
    """
    Plots the LBL using the pyvista plotting library

    Optional arguments:
    - show_points (boolean): draw points that are part of the LBL
      (default: True)
    - subplot (boolean): set this to True if you want to plot another LBL
        in the same window -it will be rendered later in that case
        (default: False)
    - show_grid (string): 'custom', 'full box', 'pyvista' or None
        can be extended in future versions (default: 'custom')
    - color (string): color used to plot this object (default: random for each
        sub-LBL of the union)
    - supplementary arguments are transmitted to the plotter function
      (pyvista's plotter.show(**kwargs))
      Example usage: you can set option "interactive_update=True" to
      continue running python and get another plot window. Rendering will be
      ensured by the latest call to LBL.plot(). All currently opened plot
      windows will be closed when one of them is closed.
    """
    global _lbl_plot_window

    node = lbl._lbl
    if node is None:
        return

    # get the last opened window
    plotter = _lbl_plot_window[-1]  # inherits from pv.plotter


    # build a python list of single LBL views and add it to plotter.LBLs:
    while node is not None:
        p = node.P
        if node.Lat.nbrows - 1 >= 3:
            plotter._3D = True
        lat = _lat3D(node.Lat) # always project to 3D, one way or another!
        # loop over polyhedra inside node
        while p is not None:
            # add each single LBL view to the list
            view = LBL_view()
            view.lbl = lbl
            view.color = color
            view.lat = lat
            view.poly_coordinate = p
            view.poly_convex_hull = p.image(lat)
            view.show_points = show_points
            # keep a local list of vertices, rays, lines in the view for later use
            view.vrl = _get_vertices(view.poly_convex_hull)
            if view.vrl[1] or view.vrl[2]:
                view.is_bounded = False
                plotter.all_bounded = False
            else:
                view.is_bounded = True
            plotter.LBLs.append(view)

            p = p.next
        node = node.next


    # if subplot is True do not render anything yet: the next call will add
    # some LBLs in the same window and do the rendering
    if not subplot:
        # this is the main part that will render all single LBLs
        _render_LBLs(plotter)

        if not plotter._3D:
            # this is a 2D plot in a 3D scene, disable z-axis:
            plotter.enable_2d_style()
            plotter.view_xy()

        if show_grid == "custom":
            _show_grid(plotter)
        elif show_grid == "full box":
            _show_grid(plotter, full_bbox=True)
        elif show_grid == "pyvista":
            plotter.show_grid()
        else:
            print(f"unknown option: show_grid = '{show_grid}'")
        plotter.reset_camera()
        plotter.enable_parallel_projection()

        # set the same scale in all opened windows (if there are more than one)
        if len(_lbl_plot_window) > 1:
            _set_global_scale()

        # get in main graphic loop
        # (blocking unless you set "interactive_update=True")
        plotter.show(**kwargs)

        # you cannot keep older windows open once you leave the main loop
        if "interactive_update" not in kwargs:
            # ensure everything is closed properly
            pv.close_all()
            _lbl_plot_window = []

        # will open a new window in a subsequent call
        # (not done if subplot==True)
        _lbl_plot_window.append(MyWindow())



def _render_LBLs(plotter):
    # If some lbls are unbounded, compute a bounding box globally
    # (modify the plotter.LBLs)
    bbox = None
    if not plotter.all_bounded:
        # get the matrix of constraints of a bounding box
        bbox = _bounding_box_lbl(plotter.LBLs)

    # loop over lattices and polyhedral domains to plot them
    for view in plotter.LBLs:
        lat, poly_coordinate, poly_convex_hull = view.lat, view.poly_coordinate, view.poly_convex_hull
        if bbox:
            poly_convex_hull = poly_convex_hull.add_constraints(bbox)
            bounded_coord = poly_convex_hull.preimage(lat)
            poly_coordinate = poly_coordinate.intersect(bounded_coord)
        if view.color is None:
            view.color = plotter.color()
        mesh = _poly2pyvista(poly_convex_hull, view.poly_convex_hull, view.vrl)
        if mesh:
            plotter.add_mesh(mesh,
                            show_edges=True,
                            opacity=.25, color=view.color)
            # draw also lines and rays originating from vertices -> does not render nice, suppressed for now

        if view.show_points:
            points = list(pypolylib._iter_single_lbl(lat, poly_coordinate, set()))
            plotter.add_mesh(
                pv.PolyData(np.asarray(points, dtype=float)),
                opacity=.5, color=view.color,
                point_size=16, render_points_as_spheres=True
            )


def _set_global_scale():

    maxscale = max(p.camera.parallel_scale for p in _lbl_plot_window)

    for p in _lbl_plot_window:
        p.camera.parallel_scale = maxscale
    for p in _lbl_plot_window[:-1]:
        p.update()


def _darker(color, factor=0.7):
    """
    Convert a color to darker (or lighter) version.
    """
    c = pv.Color(color)
    h, l, s = colorsys.rgb_to_hls(*c.float_rgb)
    l *= factor
    return colorsys.hls_to_rgb(h, l, s)


def _show_grid(window, color="red", full_bbox=False):
    """
    Draw the coordinates of the origin, axes along i, j(, k), and
    grids in the basis planes.
    """
    def _line(p1, p2):
        plotter.add_mesh(pv.Line(p1, p2), line_width=1,
                        show_scalar_bar=False, color=color)

    plotter = window
    _3D = window._3D
    text_color = _darker(color, .4)

    # set bounds
    xmin, xmax, ymin, ymax, zmin, zmax = plotter.bounds
    xmin = np.floor(xmin); ymin = np.floor(ymin); zmin = np.floor(zmin)
    xmax = np.ceil(xmax); ymax = np.ceil(ymax); zmax = np.ceil(zmax)

    if len(_lbl_plot_window) <= 1:
        origin = np.array([xmin, ymin, zmin])
        if xmin < 0 and xmax > 0: origin[0] = 0.
        if ymin < 0 and ymax > 0: origin[1] = 0.
        if zmin < 0 and zmax > 0: origin[2] = 0.
        # save the origin of this plot, will be reused
        # by all plots if there are multiple ones
        plotter.my_origin = origin
    else:
        origin = _lbl_plot_window[0].my_origin

    # draw axes:
    plotter.add_mesh(
        pv.Line((xmin, origin[1], origin[2]), (xmax, origin[1], origin[2])),
        line_width=3, show_scalar_bar=False, color=color,)
    plotter.add_mesh(
        pv.Line((origin[0], ymin, origin[2]), (origin[0], ymax, origin[2])),
        line_width=3, show_scalar_bar=False, color=color,)
    if _3D:
        plotter.add_mesh(
            pv.Line((origin[0], origin[1], zmin), (origin[0], origin[1], zmax)),
            line_width=3, show_scalar_bar=False, color=color,)

        # and origin coordinates
        plotter.add_point_labels([origin],
            [f"({origin[0]:g}, {origin[1]:g}, {origin[2]:g}) "],
            always_visible=True,
            text_color=text_color, shape_opacity=0.0,
            justification_horizontal = "right", justification_vertical = "top",
        )
    else:
        # 2D, only (i, j):
        plotter.add_point_labels([origin],
            [f"({origin[0]:g}, {origin[1]:g}) "],
            always_visible=True,
            text_color=text_color, shape_opacity=0.0,
            justification_horizontal = "right", justification_vertical = "top",
        )

    # i/j/k directions
    if xmax != origin[0]:
        plotter.add_point_labels(
            [(xmax, origin[1], origin[2])],
            [f"i={int(xmax)}"],
            always_visible=True,
            text_color=text_color, shape_opacity=0.0,
            justification_horizontal = "right", justification_vertical = "top",
        )
    if ymax != origin[1]:
        plotter.add_point_labels(
            [(origin[0], ymax, origin[2])],
            [f"j={int(ymax)}"],
            always_visible=True,
            text_color=text_color, shape_opacity=0.0,
            justification_horizontal = "right", justification_vertical = "top",
        )
    if _3D:
        if zmax != origin[2]:
            plotter.add_point_labels(
                [(origin[0], origin[1], zmax)],
                [f"k={int(zmax)}"],
                always_visible=True,
                text_color=text_color, shape_opacity=0.0,
                justification_horizontal = "right",
                justification_vertical = "top",
            )

    # scale *10 if too many points
    LARGE = 15
    xscale = yscale = zscale = 1
    while (xmax-xmin)/xscale > LARGE: xscale *= 10
    while (ymax-ymin)/yscale > LARGE: yscale *= 10
    if _3D:
        while (zmax-zmin)/zscale > LARGE: zscale *= 10

    if not full_bbox:
        # grid at origin in x/y/z planes:
        for x in range(int(xmin), int(xmax) + 1, xscale):
            # in plane z = origin[2]
            _line((x, ymin, origin[2]), (x, ymax, origin[2]))
            if _3D:
                # in plane y = origin[1]
                _line((x, origin[1], zmin), (x, origin[1], zmax))

        for y in range(int(ymin), int(ymax) + 1, yscale):
            # in plane z = origin[2]
            _line((xmin, y, origin[2]), (xmax, y, origin[2]))
            if _3D:
                # in plane x = origin[0]
                _line((origin[0], y, zmin),(origin[0], y, zmax))

        if _3D:
            for z in range(int(zmin), int(zmax) + 1, zscale):
                # in plane x = origin[0]
                _line((origin[0], ymin, z), (origin[0], ymax, z))
                # in plane y = origin[1]
                _line((xmin, origin[1], z), (xmax, origin[1], z))

    else: # if full_bbox:
        # grid at min/max in x/y/z planes:
        for x in range(int(xmin), int(xmax) + 1, xscale):
            # in plane z = origin[2]->zmin and zmax
            _line((x, ymin, zmin), (x, ymax, zmin))
            _line((x, ymin, zmax), (x, ymax, zmax))
            if _3D:
                # in plane y = origin[1]->ymin and ymax
                _line((x, ymin, zmin), (x, ymin, zmax))
                _line((x, ymax, zmin), (x, ymax, zmax))

        for y in range(int(ymin), int(ymax) + 1, yscale):
            # in plane z = origin[2]->zmin and zmax
            _line((xmin, y, zmin), (xmax, y, zmin))
            _line((xmin, y, zmax), (xmax, y, zmax))
            if _3D:
                # in plane x = origin[0]->xmin and xmax
                _line((xmin, y, zmin),(xmin, y, zmax))
                _line((xmax, y, zmin),(xmax, y, zmax))

        if _3D:
            for z in range(int(zmin), int(zmax) + 1, zscale):
                # in plane x = origin[0]->
                _line((xmin, ymin, z), (xmin, ymax, z))
                _line((xmax, ymin, z), (xmax, ymax, z))
                # in plane y = origin[1]->
                _line((xmin, ymin, z), (xmax, ymin, z))
                _line((xmin, ymax, z), (xmax, ymax, z))

def _bounding_box_lbl(lbl_list):
    """Compute a bounding box of an unbounded lbl (only necessary constraints).

    Get all vertices and rays/lines, then compute a box of
      vertices + ZOOM*rays +- ZOOM*lines.

    lbl_list is a list of: (lat, poly_coordinate, poly_convex_hull, is_bounded)
    """
    if not lbl_list:
        # empty lbl?
        return None
    # all poly's (convex hulls) have the same dimension
    dim = lbl_list[0].poly_convex_hull.dimension

    # all vertices/rays/lines of all hulls
    vertices = []
    rays = []
    lines = []
    for view in lbl_list:
        vertices.extend(view.vrl[0])
        rays.extend(view.vrl[1])
        lines.extend(view.vrl[2])

    # vertices is a non-empty list of tuples of floats
    # rays+lines is a non-empty list of tuples of integers
    mini = [0] * dim
    maxi = [0] * dim
    for d in range(dim):
        mini[d] = min(vertices[v][d] for v in range(len(vertices)))
        maxi[d] = max(vertices[v][d] for v in range(len(vertices)))
    # we now have a min and max value of vertices coordinates

    # scan the rays to add them to the mini/maxi values if needed
    for r in rays:
        for d in range(dim):
            mini[d] = min(mini[d], mini[d] + ZOOM * r[d])
            maxi[d] = max(maxi[d], maxi[d] + ZOOM * r[d])
    # scan the lines to add/substract them to the mini/maxi values if needed
    for l in lines:
        for d in range(dim):
            mini[d] = min(mini[d], mini[d] + ZOOM * l[d])
            mini[d] = min(mini[d], mini[d] - ZOOM * l[d])
            maxi[d] = max(maxi[d], maxi[d] + ZOOM * l[d])
            maxi[d] = max(maxi[d], maxi[d] - ZOOM * l[d])

    # compute the constraints to be added to the hull polyhedra:
    constraints_str = f"{dim * 2} {dim + 2}\n"
    for i in range(dim):
        cons = ["1"] + ["0"] * (dim + 1)
        cons[dim + 1] = str(math.floor(-mini[i]))
        cons[i + 1] = "1"    # x >= mini[d]
        constraints_str += " ".join(cons) + "\n"
        cons[dim + 1] = str(math.ceil(maxi[i]))
        cons[i + 1] = "-1"    # x <= maxi[d]
        constraints_str += " ".join(cons) + "\n"
    bbox_constraints = pl.matrix_read_from_string(constraints_str)

    return bbox_constraints


def _poly2pyvista(poly, poly_unbounded, vrl):
    """Transform a PolyLib 3D Polyhedron into a PyVista polyhedron."""

    # vertices (of the bounded poly): as the list of tuples of FP coordinates
    vertices, _, _ = _get_vertices(poly)
    if not vertices:
        return []
    vertices = np.asarray(vertices, dtype=float)

    # build faces: as the (flat) list of faces
    # [num_vertices, vertex0, vertex1..., num_vertices2, vertex2, vertex3...]
    faces = []
    dim = poly.dimension
    # constraints from the unbounded polyhedron, rays from the bounded one.
    constraint = poly_unbounded.constraint

    ray = poly.ray
    # scan constraints to build faces:
    for c in range(poly_unbounded.nbconstraints):
        # ignore the constraints from the bounding box to draw an open
        # polyhedron if unbounded.

        normal = tuple(constraint[c, d+1] for d in range(dim))
        face = []
        for r in range(poly.nbrays):
            sat = sum(constraint[c, d+1]*ray[r, d+1] for d in range(dim+1))
            # does ray r saturate constraint c?
            if sat == 0:
                face.append(r)
        if len(face) >= 3:
            face = _order_face_indices(vertices, face, normal)
            faces.extend([len(face)] + face)
    if faces:
        # # thought that was a good idea, but it's not so nice (not visible):
        # # add lines around vertices: in the directions of + ray and +- line.
        # infinite_directions = []
        # for v in vertices:
        #     for r in vrl[1]:
        #         # rays
        #         infinite_directions.append(
        #             pv.Line((v[0], v[1], v[2]),
        #                     (v[0]+r[0], v[1]+r[1], v[2]+r[2]))
        #         )
        #     for l in vrl[2]:
        #         # lines
        #         infinite_directions.append(
        #             pv.Line((v[0], v[1], v[2]),
        #                     (v[0]+l[0], v[1]+l[1], v[2]+l[2]))
        #         )
        #         infinite_directions.append(
        #             pv.Line((v[0], v[1], v[2]),
        #                     (v[0]-l[0], v[1]-l[1], v[2]-l[2]))
        #         )

        return pv.PolyData(vertices, faces) # + infinite_directions
    return []


def _order_face_indices(vertices, face_indices, normal):
    face_indices = np.asarray(face_indices, dtype=int)
    ordered_local = _order_face(vertices[face_indices], normal)
    return [face_indices[i] for i in ordered_local]


def _order_face(vertices, normal):
    """Order coplanar 3D vertices cyclically."""
    normal = np.asarray(normal, dtype=float)
    # Normalize normal
    normal /= np.linalg.norm(normal)
    if len(vertices) == 3:
        order = np.array([0,1,2])
    else:
        # Build an orthonormal basis (u,v,n): choose a vector not parallel to n
        if abs(normal[0]) < 0.9:
            tmp = np.array([1., 0., 0.])
        else:
            tmp = np.array([0., 1., 0.])
        u = np.cross(normal, tmp)
        u /= np.linalg.norm(u)
        v = np.cross(normal, u)

        # Project points onto the face plane
        center = vertices.mean(axis=0)
        centered = vertices - center
        pts2d = np.column_stack((centered @ u, centered @ v))

        # Convex hull in 2D
        hull = ConvexHull(pts2d)
        order = hull.vertices

    # Fix orientation
    p0, p1, p2 = vertices[order[:3]]
    face_normal = np.cross(p1 - p0, p2 - p0)
    if np.dot(face_normal, normal) < 0:
        order = order[::-1]
    return order


def _get_vertices(poly):
    """Return the vertices/rays/lines of poly as lists of tuples of floats.
    """
    vertices = []
    rays = []
    lines = []

    dim = poly.dimension
    poly_ray = poly.ray
    for v in range(poly.nbrays):
        if poly_ray[v, 0] == 0:
            # line
            lines.append(tuple(poly_ray[v, x+1] for x in range(dim)))
        else:
            v_div = poly_ray[v, dim + 1]
            if v_div == 0:
                # ray
                rays.append(tuple(poly_ray[v, x+1] for x in range(dim)))
            else:
                # vertex
                vertices.append(tuple(poly_ray[v, x+1]/v_div for x in range(dim)))

    return vertices, rays, lines


def _lat3D(node_lat):
    """Change a Lattice to be a 3D representation/projection.

    - if lower than 3D -> set the plane z=0 (and y=0 if 1D)
    - if greater than 3D -> project out extra dimensions (TODO: could be improved)
    """
    if node_lat.nbrows == 4:
        # 3D is ok
        return node_lat
    else:
        # extend node_lat to be a 3D projection!
        # - if larger than 3D: project out extra dimensions
        # - if smaller than 3D: set z=0 (and y=0 if 1D)
        lat = pl.matrix_alloc(4, node_lat.nbcolumns)
        for i in range(min(3, node_lat.nbrows - 1)):
            for j in range(node_lat.nbcolumns):
                lat[i, j] = node_lat[i, j]
        lat[3, lat.nbcolumns - 1] = 1
        return lat
