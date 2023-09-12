from abc import ABC, abstractmethod
from dataclasses import dataclass

from .TopoCoordinates import TopoCoordinates


class Lattice(ABC):
    """Abstract class that represents a lattice."""

    dimensions: list[int]  # Dimensions of the lattice
    cell_values: dict[
        TopoCoordinates, bool
    ]  # Values of the cells of the lattice (0 : empty, 1 : occupied)

    @property
    def dimensions(self) -> list[int]:
        """Getter for the attribute dimensions of the lattice.

        Returns
        -------
        list[int]
            Dimensions of the lattice.
        """
        return self._dimensions

    @dimensions.setter
    def dimensions(self, dimensions: list[int]) -> None:
        """Setter for the attribute dimensions of the lattice.

        Parameters
        ----------
        dimensions : list[int]
            Dimensions of the lattice to be assigned.
        """
        self._dimensions = dimensions

    @property
    def cell_values(self) -> dict[TopoCoordinates, bool]:
        """Getter for the attribute cell_values of the lattice.

        Returns
        -------
        dict[TopoCoordinates, bool]
            Values of the cells of the lattice.
        """
        return self._cell_values

    @cell_values.setter
    def cell_values(self, cell_values: dict[TopoCoordinates, bool]) -> None:
        """Setter for the attribute cell_values of the lattice.

        Parameters
        ----------
        cell_values : dict[TopoCoordinates, bool]
            Values of the cells of the lattice to be assigned.
        """
        self._cell_values = cell_values

    @abstractmethod
    def are_topological_neighbours(
        self, cell1: TopoCoordinates, cell2: TopoCoordinates
    ) -> bool:
        """Checks if two cells are topological neighbours.

        Parameters
        ----------
        cell1 : TopoCoordinates
            First cell.
        cell2 : TopoCoordinates
            Second cell.

        Returns
        -------
        bool
            True if the cells are topological neighbours, False otherwise.
        """
        pass

    @abstractmethod
    def get_topological_neighbours(
        self, cell: TopoCoordinates
    ) -> list[TopoCoordinates]:
        """Returns the topological neighbours of a cell.

        Parameters
        ----------
        cell : TopoCoordinates
            Cell.

        Returns
        -------
        list[TopoCoordinates]
            Topological neighbours of the cell.
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
