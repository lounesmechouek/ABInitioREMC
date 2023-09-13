from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Tuple

from .TopoCoordinates import TopoCoordinates


@dataclass(slots=True)
class Lattice(ABC):
    """Abstract class that represents a lattice."""

    _dimensions: Tuple[int, ...]  # Dimensions of the lattice (x, y)

    _cell_values: dict[
        Tuple[int, ...], bool
    ]  # Values of the cells of the lattice (False : empty, True : occupied

    @property
    def dimensions(self) -> Tuple[int, ...]:
        """Getter for the attribute dimensions of the lattice.

        Returns
        -------
        tuple[int, ...]
            Dimensions of the lattice.
        """
        return self._dimensions

    @dimensions.setter
    def dimensions(self, dimensions: Tuple[int, ...]) -> None:
        """Setter for the attribute dimensions of the lattice.

        Parameters
        ----------
        dimensions : tuple[int, ...]
            Dimensions of the lattice to be assigned.
        """
        self._dimensions = dimensions

    @property
    def cell_values(self) -> dict[Tuple[int, ...], bool]:
        """Getter for the attribute cell_values of the lattice.

        Returns
        -------
        dict[Tuple[int, ...], bool]
            Values of the cells of the lattice.
        """
        return self._cell_values

    @cell_values.setter
    def cell_values(self, cell_values: dict[Tuple[int, ...], bool]) -> None:
        """Setter for the attribute cell_values of the lattice.

        Parameters
        ----------
        cell_values : dict[Tuple[int, ...], bool]
            Values of the cells of the lattice to be assigned.
        """
        self._cell_values = cell_values

    @abstractmethod
    def are_adjacent(self, cell1: TopoCoordinates, cell2: TopoCoordinates) -> bool:
        """Checks if two cells are adjacent.

        Parameters
        ----------
        cell1 : TopoCoordinates
            First cell.
        cell2 : TopoCoordinates
            Second cell.

        Returns
        -------
        bool
            True if the cells are adjacent. False otherwise.
        """
        pass

    @abstractmethod
    def compute_end_moves(self, cell: TopoCoordinates) -> list[TopoCoordinates]:
        """Computes the end moves of a cell.

        Parameters
        ----------
        cell : TopoCoordinates
            Cell.

        Returns
        -------
        list[TopoCoordinates]
            New possible positions of the cell after end moves.
        """
        pass

    @abstractmethod
    def compute_corner_moves(self, cell: TopoCoordinates) -> list[TopoCoordinates]:
        """Computes the corner moves of a cell.

        Parameters
        ----------
        cell : TopoCoordinates
            Cell.

        Returns
        -------
        list[TopoCoordinates]
            New possible positions of the cell after corner moves.
        """
        pass

    @abstractmethod
    def compute_crankshaft_moves(self, cell: TopoCoordinates) -> list[TopoCoordinates]:
        """Computes the crankshaft moves of a cell.

        Parameters
        ----------
        cell : TopoCoordinates
            Cell.

        Returns
        -------
        list[TopoCoordinates]
            New possible positions of the cell after crankshaft moves.
        """
        pass

    @abstractmethod
    def compute_pull_moves(self, cell: TopoCoordinates) -> list[TopoCoordinates]:
        """Computes the pull moves of a cell.

        Parameters
        ----------
        cell : TopoCoordinates
            Cell.

        Returns
        -------
        list[TopoCoordinates]
            New possible positions of the cell after pull moves.
        """
        pass
