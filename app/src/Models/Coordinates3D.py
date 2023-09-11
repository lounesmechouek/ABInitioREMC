from dataclasses import dataclass
from typing import Tuple

from .TopoCoordinates import TopoCoordinates


@dataclass(slots=True)
class Coordinates3D(TopoCoordinates):
    """Coordinates3D is a class that represents coordinates of an element in NxNxN."""

    _coordinates: Tuple[int, int, int]  # Coordinates of the element in NxNxN

    @property
    def coordinates(self) -> Tuple[int, int, int]:
        """Getter for the attribute coordinates of the element in NxNxN.

        Returns
        -------
        Tuple[int, int, int]
            Coordinates of the element in NxNxN.
        """
        return self._coordinates

    @coordinates.setter
    def coordinates(self, coordinates: Tuple[int, int, int]) -> None:
        """Setter for the attribute coordinates of the element in NxN.

        Parameters
        ----------
        coordinates : Tuple[int, int, int]
            Coordinates of the element in NxNxN to be assigned.
        """
        self._coordinates = coordinates

    def topo_coordinates_model(self) -> str:
        return "NxNxN"
