from dataclasses import dataclass
from typing import Tuple

from .TopoCoordinates import TopoCoordinates


@dataclass(slots=True)
class Coordinates2D(TopoCoordinates):
    """Coordinates2D is a class that represents coordinates of an element in NxN."""

    _coordinates: Tuple[int, int]  # Coordinates of the element in NxN

    @property
    def coordinates(self) -> Tuple[int, int]:
        """Getter for the attribute coordinates of the element in NxN.

        Returns
        -------
        Tuple[int, int]
            Coordinates of the element in NxN.
        """
        return self._coordinates

    @coordinates.setter
    def coordinates(self, coordinates: Tuple[int, int]) -> None:
        """Setter for the attribute coordinates of the element in NxN.

        Parameters
        ----------
        coordinates : Tuple[int, int]
            Coordinates of the element in NxN to be assigned.
        """
        self._coordinates = coordinates

    def topo_coordinates_model(self) -> str:
        return "NxN"
