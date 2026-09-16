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


class Axes4D:
    """
    Small interactive 4D coordinate frame displayed in a corner
    of a PyVista render window.

    The four 4D basis vectors are projected onto a 2D plane:

                  X
                  ↑
             Y ←──┼──→ Z
                  │
                  ↓
                  W
    """

    def __init__(
        self,
        plotter,
        viewport=(0.02, 0.02, 0.22, 0.22),
    ):
        self.plotter = plotter
        self.viewport = viewport

        # Projection of the 4D basis vectors onto the widget plane.
        self.directions = {
            "X": np.array([0.0, 1.0]),
            "Y": np.array([-1.0, 0.0]),
            "Z": np.array([1.0, 0.0]),
            "W": np.array([0.0, -1.0]),
        }

        # Dedicated renderer.
        self.renderer = vtkRenderer()
        self.renderer.SetViewport(*viewport)
        self.renderer.SetLayer(1)

        # Transparent background.
        self.renderer.SetBackgroundAlpha(0.0)
        self.renderer.SetInteractive(1)

        # Main renderer is layer 0, our renderer is layer 1.
        self.plotter.render_window.SetNumberOfLayers(2)
        self.plotter.render_window.AddRenderer(self.renderer)

        self.actors = {}

        self._create_axes()
        self._setup_camera()
        self._setup_picker()

        # Handle clicks:
        self.dragged_axis = None
        self.last_mouse_position = None
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
        for name, direction in self.directions.items():
            actor = self._create_axis(direction)
            self.renderer.AddActor(actor)
            self.actors[name] = actor

    def _create_axis(self, direction):
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

        # vtkArrowSource points along +X.
        #
        # Rotate it in the widget plane so that +X becomes
        # the desired 2D direction.
        angle = np.degrees(
            np.arctan2(direction[1], direction[0])
        )

        actor.RotateZ(angle)

        # ArrowSource has length 1.
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

    def _setup_picker(self):
        self.picker = vtkPropPicker()

        self.plotter.iren.add_observer(
            "LeftButtonPressEvent",
            self._on_left_button_press,
        )

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
                self.last_mouse_position = (x, y)

                print(f"Start dragging {name}")
                break

    def _on_mouse_move(self, caller, event):
        if self.dragged_axis is None:
            return

        x, y = caller.GetEventPosition()

        old_x, old_y = self.last_mouse_position

        dx = x - old_x
        dy = y - old_y

        self.last_mouse_position = (x, y)

        self.axis_dragged(
            self.dragged_axis,
            dx,
            dy,
            x,
            y,
        )

    def _on_left_button_release(self, caller, event):
        if self.dragged_axis is not None:
            print(f"Stop dragging {self.dragged_axis}")

        self.dragged_axis = None
        self.last_mouse_position = None

    # ------------------------------------------------------------------
    # User action
    # ------------------------------------------------------------------
    def axis_dragged(self, name, dx, dy, x, y):
        print(
            f"{name}: dx={dx}, dy={dy}, "
            f"position=({x}, {y})"
        )

    def axis_clicked(self, name):
        print(f"Clicked on 4D axis: {name}")
