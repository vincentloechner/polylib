"""
lbl_plot_nD_axes.py — Visualizing LBLs:
Handle an axes interactive frame in the corner of the window.
"""
import numpy as np
import pyvista as pv

from vtkmodules.vtkFiltersSources import vtkArrowSource
from vtkmodules.vtkRenderingCore import (
    vtkActor,
    vtkPolyDataMapper,
    vtkPropPicker,
    vtkRenderer,
)


class AxesnD:
    """
    Small interactive 4D (test) coordinate frame displayed in a corner
    of a PyVista render window.

    The four 4D basis vectors are projected onto a 2D plane:

                  i
                  ↑
             j ←──┼──→ k
                  │
                  ↓
                  l
    """

    def __init__(
        self,
        plotter,
        viewport=(0.02, 0.02, 0.22, 0.22),
    ):
        print("Init Axes4D")
        self.plotter = plotter
        self.viewport = viewport

        self.axis_names = ["i", "j", "k", "l"]

        # Projection of the 4D basis vectors onto the widget plane.
        # Each column is the 2D projection of one 4D basis vector.
        self.projection = np.array([
            [ 0.0, -1.0,  1.0,  0.0],
            [ 1.0,  0.0,  0.0, -1.0],
        ])

        # Dedicated renderer.
        self.renderer = vtkRenderer()
        self.renderer.SetViewport(*viewport)
        self.renderer.SetLayer(1)

        # Transparent background.
        self.renderer.SetBackgroundAlpha(0.0)
        # self.renderer.SetInteractive(False)
        self.renderer.SetInteractive(True)

        # Main renderer is layer 0, our renderer is layer 1.
        self.plotter.render_window.SetNumberOfLayers(2)
        self.plotter.render_window.AddRenderer(self.renderer)

        self.actors = {}
        self.picker = vtkPropPicker()

        self._create_axes()
        self._setup_camera()

        # Handle clicks:
        self.dragged_axis = None
        self.plotter.iren.add_observer(
            "LeftButtonPressEvent",
            self._on_left_button_press,
        )
        self.plotter.iren.add_observer(
            "MouseMoveEvent",
            self._on_mouse_move,
        )
        self.plotter.iren.add_observer(
            "LeftButtonReleaseEvent",
            self._on_left_button_release,
        )

    # ------------------------------------------------------------------
    # Axes
    # ------------------------------------------------------------------

    def _create_axes(self):
        print(f"creating axes")
        for i, name in enumerate(self.axis_names):
            # direction = self.projection[:, i]
            print(f"  {i}: {name}")

            actor = self._create_axis()

            self.renderer.AddActor(actor)
            self.actors[name] = actor

            self._update_axis(i)

    def _create_axis(self):
        arrow = vtkArrowSource()

        arrow.SetTipLength(0.20)
        arrow.SetTipRadius(0.07)
        arrow.SetShaftRadius(0.025)

        arrow.SetShaftResolution(12)
        arrow.SetTipResolution(12)

        mapper = vtkPolyDataMapper()
        mapper.SetInputConnection(arrow.GetOutputPort())

        actor = vtkActor()
        actor.SetMapper(mapper)

        actor.SetScale(0.75)

        return actor

    # ------------------------------------------------------------------
    # Camera
    # ------------------------------------------------------------------
    def _setup_camera(self):
        camera = self.renderer.GetActiveCamera()

        camera.ParallelProjectionOn()

        camera.SetPosition(0.0, 0.0, 10.0)
        camera.SetFocalPoint(0.0, 0.0, 0.0)
        camera.SetViewUp(0.0, 1.0, 0.0)

        camera.SetParallelScale(1.8)

    # ------------------------------------------------------------------
    # Picking
    # ------------------------------------------------------------------
    def _on_left_button_press(self, caller, event):
        x, y = caller.GetEventPosition()

        width, height = self.plotter.render_window.GetSize()

        xmin, ymin, xmax, ymax = self.viewport

        if not (
            xmin * width <= x <= xmax * width
            and ymin * height <= y <= ymax * height
        ):
            return

        self.picker.Pick(x, y, 0, self.renderer)

        picked_actor = self.picker.GetActor()

        for name, actor in self.actors.items():
            if picked_actor == actor:
                self.dragged_axis = name

                print(f"Start dragging {name}")
                break

    def _on_mouse_move(self, caller, event):
        if self.dragged_axis is None:
            return

        x, y = caller.GetEventPosition()

        self.axis_dragged(
            self.dragged_axis,
            x,
            y,
        )

        self.plotter.render()

    def _on_left_button_release(self, caller, event):
        if self.dragged_axis is not None:
            print(f"Stop dragging {self.dragged_axis}")

        self.dragged_axis = None

    # ------------------------------------------------------------------
    # User action
    # ------------------------------------------------------------------
    def axis_dragged(self, name, x, y):
        i = self.axis_names.index(name)

        width, height = self.plotter.render_window.GetSize()

        xmin, ymin, xmax, ymax = self.viewport

        # Position du centre du widget en pixels.
        cx = 0.5 * (xmin + xmax) * width
        cy = 0.5 * (ymin + ymax) * height

        # Taille du widget.
        widget_width = (xmax - xmin) * width
        widget_height = (ymax - ymin) * height

        # Position de la souris par rapport au centre,
        # normalisée par la taille du widget.
        vx = (x - cx) / widget_width * 2.0
        vy = (cy - y) / widget_height * 2.0

        # Limite éventuellement la longueur de l'axe.
        v = np.array([vx, vy])

        length = np.linalg.norm(v)

        if length > 0.9:
            v *= 0.9 / length

        self.projection[:, i] = v

        self._update_axis(i)

    def _update_axis(self, i):
        name = self.axis_names[i]
        actor = self.actors[name]

        v = self.projection[:, i]

        length = np.linalg.norm(v)

        if length < 1e-6:
            actor.SetVisibility(False)
            return

        actor.SetVisibility(True)

        angle = np.degrees(np.arctan2(v[1], v[0]))

        actor.SetPosition(0.0, 0.0, 0.0)
        actor.SetOrientation(0.0, 0.0, angle)
        actor.SetScale(length)

        self.renderer.ResetCameraClippingRange()