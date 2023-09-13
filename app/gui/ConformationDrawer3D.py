import plotly.graph_objects as go

from ..src.Models.Conformation3D import Conformation3D
from ..src.Models.Polarity import Polarity
from .ConformationDrawer import ConformationDrawer


class ConformationDrawer3D(ConformationDrawer):
    """Class used to represent protein conformations in 3D."""

    def __init__(self, conformation: Conformation3D, colors: dict[str, str]) -> None:
        """Constructor for the ConformationDrawer3D class.

        Parameters
        ----------
        conformation : Conformation2D
            Conformation to be drawn.
        """
        super().__init__(conformation, colors)

    def draw(self) -> None:
        """Draws the conformation."""
        # Points to be drawn
        i = 0
        points = []
        links = []
        colors = []

        for coords, residue in self._conformation.amino_acid_coordinates.items():
            point = {}
            point["x"] = coords[0]
            point["y"] = coords[1]
            point["z"] = coords[2]

            if residue.polarity == Polarity.POLAR:
                point["color"] = self._colors["P"]
                colors.append(self._colors["P"])
            else:
                point["color"] = self._colors["H"]
                colors.append(self._colors["H"])

            # index of the residue in the protein
            found_index = None

            for index, amino_acid in enumerate(self.conformation.protein.sequence):
                if amino_acid.id == residue.id:
                    found_index = index
                    break  # Stop searching after the first match

            if found_index is None:
                point["label"] = "?"
            else:
                point["label"] = str(found_index)

            points.append(point)

            if i > 0:
                link = (str(i - 1), str(i))
                links.append(link)

            i += 1

        fig = go.Figure()

        for point in points:
            fig.add_trace(
                go.Scatter3d(
                    x=[point["x"]],
                    y=[point["y"]],
                    z=[point["z"]],
                    text=[point["label"]],
                    mode="markers+text",
                    marker=dict(color=point["color"]),
                )
            )

        for line in links:
            start_point = [point for point in points if point["label"] == line[0]][0]
            end_point = [point for point in points if point["label"] == line[1]][0]
            fig.add_trace(
                go.Scatter3d(
                    x=[start_point["x"], end_point["x"]],
                    y=[start_point["y"], end_point["y"]],
                    z=[start_point["z"], end_point["z"]],
                    mode="lines",
                    line=dict(color="blue"),
                )
            )

        fig.update_layout(
            scene=dict(
                xaxis_title="",
                yaxis_title="",
                zaxis_title="",
            ),
            title="3D Representation of the protein",
            showlegend=False,
        )

        fig.update_traces(
            textfont_size=14, textposition="top left", marker=dict(size=20)
        )

        fig.update_layout(
            scene=dict(
                xaxis_showspikes=False,  # Hide x-axis
                yaxis_showspikes=False,  # Hide y-axis
                zaxis_showspikes=False,  # Hide z-axis
            )
        )

        return fig
