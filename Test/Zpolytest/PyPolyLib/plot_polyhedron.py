import numpy as np
import matplotlib.pyplot as plt
from shapely.geometry import Point, MultiPoint
# from shapely.geometry.polygon import Polygon

# BIG = 1000

# ----------------------
#   FONCTIONS INTERNES
# ----------------------

def solve_intersection(line1, line2, tol=1e-12):
    a1, b1, c1 = line1
    a2, b2, c2 = line2
    A = np.array([[a1, b1], [a2, b2]])
    det = np.linalg.det(A)
    if abs(det) < tol:
        return None
    return tuple(np.linalg.solve(A, np.array([-c1, -c2])).tolist())


def compute_polyhedron(inequalities, tol=1e-8):
    """
    Retourne le polygone Shapely correspondant,
    ou None si vide ou non-borné.
    """
    vertices = []
    n = len(inequalities)
    for i in range(n):
        for j in range(i+1, n):
            p = solve_intersection(inequalities[i], inequalities[j])
            if p is None:
                continue
            x, y = p
            ok = True
            for (a, b, c) in inequalities:
                if a*x + b*y + c < -tol:
                    ok = False
                    break
            if ok:
                vertices.append((x, y))

    if len(vertices) == 0:
        return None  # vide ou non-borné

    mp = MultiPoint(vertices)
    region = mp.convex_hull

    if region.is_empty or region.geom_type not in ("Polygon", "MultiPolygon"):
        return None

    if region.geom_type == "MultiPolygon":
        region = max(region.geoms, key=lambda p: p.area)

    return region


# ----------------------
#   FONCTION PRINCIPALE
# ----------------------

def plot_polyhedra(list_of_polyhedra, show_integer_points=True, grid_margin=1):
    """
    list_of_polyhedra : liste de polyèdres
       un polyèdre = liste de tuples (a,b,c) pour ax+by+c>=0

    Exemple :
    [
        [(a1,b1,c1),(a2,b2,c2),...],
        [(d1,e1,f1),(d2,e2,f2),...],
        ...
    ]
    """

    fig, ax = plt.subplots()
    colors = plt.cm.tab10(np.linspace(0, 1, len(list_of_polyhedra)))

    for k, inequalities in enumerate(list_of_polyhedra):

        region = compute_polyhedron(inequalities)
        if region is None:
            print(f"Polyèdre {k+1} : vide ou non borné.")
            continue

        px, py = region.exterior.xy

        # Polygone
        ax.plot(px, py, color=colors[k], linewidth=2)
        ax.fill(px, py, color=colors[k], alpha=0.15)

        # Points entiers
        if show_integer_points:
            minx, miny, maxx, maxy = region.bounds
            minx = int(np.floor(minx))
            miny = int(np.floor(miny))
            maxx = int(np.ceil(maxx))
            maxy = int(np.ceil(maxy))

            integer_points = []
            for x in range(minx, maxx+1):
                for y in range(miny, maxy+1):
                    p = Point(x, y)
                    if region.contains(p) or region.touches(p):
                        integer_points.append((x, y))

            if integer_points:
                xs, ys = zip(*integer_points)
                ax.scatter(xs, ys, s=20, color=colors[k])

    ax.set_aspect('equal')
    ax.grid(True, linestyle='--', alpha=0.4)
    plt.show()

