from abc import ABC, abstractmethod
from dataclasses import dataclass


@dataclass(slots=True)
class TopoCoordinates(ABC):
    """TopoCoordinates is an abstract class that represents topological coordinates in a given space."""

    # * We use inheritance to avoid code duplication if we want to add more topological coordinates models (polar, etc.).
    # TODO: Define more abstract high-level methods for the topological coordinates if necessary.

    @abstractmethod
    def topo_coordinates_model(self) -> str:
        """Returns the type of the topological coordinates.

        Returns
        -------
        str
            Type of the topological coordinates.
        """
        pass
