import plotly.express as px

from ..src.Models.Conformation2D import Conformation2D
from ..src.Models.Polarity import Polarity
from .ConformationDrawer import ConformationDrawer


class ConformationDrawer2D(ConformationDrawer):
    """Class used to represent protein conformations in 2D."""

    def __init__(self, conformation: Conformation2D, colors: dict[str, int]) -> None:
        """Constructor for the ConformationDrawer2D class.

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
                link = {}
                link["start"] = str(i - 1)
                link["end"] = str(i)
                links.append(link)

            i += 1

        fig = px.scatter()

        for point in points:
            fig.add_trace(
                px.scatter(
                    x=[point["x"]],
                    y=[point["y"]],
                    text=[point["label"]],
                ).data[0]
            )

        for line in links:
            start_point = [
                point for point in points if point["label"] == line["start"]
            ][0]
            end_point = [point for point in points if point["label"] == line["end"]][0]
            fig.add_trace(
                px.line(
                    x=[start_point["x"], end_point["x"]],
                    y=[start_point["y"], end_point["y"]],
                ).data[0]
            )

        fig.update_layout(
            xaxis_title="",
            yaxis_title="",
            title="2D Structure of the conformation",
            showlegend=False,
        )

        colorscale = [(0, "black"), (128 / 255, "green"), (1, "white")]
        fig.update_traces(
            textfont_size=14,
            textposition="top left",
            marker=dict(size=30, color=colors, colorscale=colorscale, opacity=1),
        )

        fig.update_xaxes(showticklabels=False)
        fig.update_yaxes(showticklabels=False)

        return fig
