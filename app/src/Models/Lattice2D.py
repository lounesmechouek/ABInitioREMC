import random
from dataclasses import dataclass
from typing import Tuple

from .AminoAcidHP import AminoAcidHP
from .Coordinates2D import Coordinates2D
from .Lattice import Lattice
from .ProteinHP import ProteinHP


@dataclass(slots=True)
class Lattice2D(Lattice):
    """Class that represents a 2D lattice."""

    def __init__(self, dimensions: Tuple[int, int]) -> None:
        """Constructor for the Lattice2D class.

        Parameters
        ----------
        dimensions : tuple[int, int]
            Dimensions of the lattice.
        """
        if len(dimensions) != 2:
            raise ValueError("Dimensions must be a tuple of two integers.")

        self._dimensions = dimensions
        values = {}
        for i in range(dimensions[0]):
            for j in range(dimensions[1]):
                values[(i, j)] = False
        self._cell_values = values

    def reset_lattice(self) -> None:
        """Resets the lattice to its initial state."""
        for i in range(self._dimensions[0]):
            for j in range(self._dimensions[1]):
                self._cell_values[(i, j)] = False

    def set_cell_value(self, cell: Coordinates2D, value: bool) -> None:
        """Sets the value of a single cell.

        Parameters
        ----------
        cell : Coordinates2D
            Cell.
        value : bool
            Value to be assigned to the cell.
        """
        if (
            cell.coordinates[0] >= 0
            and cell.coordinates[1] >= 0
            and cell.coordinates[0] < self.dimensions[0]
            and cell.coordinates[1] < self.dimensions[1]
        ):
            self._cell_values[cell.coordinates] = value

        else:
            raise ValueError("Cell coordinates out of bounds.")

    def are_adjacent(self, cell1: Coordinates2D, cell2: Coordinates2D) -> bool:
        """Checks if two cells are topological neighbours.

        Parameters
        ----------
        cell1 : Coordinates2D
            First cell.
        cell2 : Coordinates2D
            Second cell.

        Returns
        -------
        bool
            True if the cells are topological neighbours, False otherwise.
        """
        # Points are not adjacent if they are the same point
        if cell1.coordinates == cell2.coordinates:
            return False

        # Check if the cells are adjacent (same row or column)
        adjacent = (
            abs(cell1.coordinates[0] - cell2.coordinates[0]) == 0
            and abs(cell1.coordinates[1] - cell2.coordinates[1]) == 1
        ) or (
            abs(cell1.coordinates[1] - cell2.coordinates[1]) == 0
            and abs(cell1.coordinates[0] - cell2.coordinates[0]) == 1
        )
        return adjacent

    def get_random_adjacent_cell(
        self, cell: Coordinates2D, exclude: list[Coordinates2D]
    ) -> Tuple[int, int]:
        """Returns a random adjacent cell.

        Parameters
        ----------
        cell : Coordinates2D
            Cell.
        exclude : list[Coordinates2D]
            List of cells to exclude from the sampling.

        Returns
        -------
        Tuple[int, int]
            Random adjacent cell.
        """
        candidates = []
        # Each cell has at most 4 adjacent cells
        if (cell.coordinates[0] - 1 >= 0) and (
            (cell.coordinates[0] - 1, cell.coordinates[1]) not in exclude
        ):
            candidates.append((cell.coordinates[0] - 1, cell.coordinates[1]))
        if (cell.coordinates[0] + 1 < self._dimensions[0]) and (
            (cell.coordinates[0] + 1, cell.coordinates[1]) not in exclude
        ):
            candidates.append((cell.coordinates[0] + 1, cell.coordinates[1]))
        if (cell.coordinates[1] - 1 >= 0) and (
            (cell.coordinates[0], cell.coordinates[1] - 1) not in exclude
        ):
            candidates.append((cell.coordinates[0], cell.coordinates[1] - 1))
        if (cell.coordinates[1] + 1 < self._dimensions[1]) and (
            (cell.coordinates[0], cell.coordinates[1] + 1) not in exclude
        ):
            candidates.append((cell.coordinates[0], cell.coordinates[1] + 1))

        if len(candidates) == 0:
            raise ValueError("Cell has no adjacent cells.")
        elif len(candidates) == 1:
            return candidates[0]
        else:
            random_index = random.randint(0, len(candidates) - 1)
            return candidates[random_index]

    def get_all_adjacent_cells(self, cell: Coordinates2D) -> list[Tuple[int, int]]:
        """Returns all adjacent cells.

        Parameters
        ----------
        cell : Coordinates2D
            Cell.

        Returns
        -------
        list[Coordinates2D]
            List of adjacent cells.
        """
        candidates = []
        # Each cell has at most 4 adjacent cells
        if cell.coordinates[0] - 1 >= 0:
            candidates.append((cell.coordinates[0] - 1, cell.coordinates[1]))
        if cell.coordinates[0] + 1 < self._dimensions[0]:
            candidates.append((cell.coordinates[0] + 1, cell.coordinates[1]))
        if cell.coordinates[1] - 1 >= 0:
            candidates.append((cell.coordinates[0], cell.coordinates[1] - 1))
        if cell.coordinates[1] + 1 < self._dimensions[1]:
            candidates.append((cell.coordinates[0], cell.coordinates[1] + 1))

        return candidates

    def compute_end_moves(
        self,
        cell: Coordinates2D,
        protein: ProteinHP,
        amino_acid_coords: dict[Tuple[int, int], AminoAcidHP],
    ) -> list[Tuple[int, int]]:
        """Computes the end moves for a given cell.

        Parameters
        ----------
        cell : Coordinates2D
            The cell to compute the end moves for.
        protein : ProteinHP
            The studied protein.
        amino_acid_coords: dict[Tuple[int, int], AminoAcidHP]
            Dictionary containing the coordinates of the amino acids in the conformation.

        Returns
        -------
        list[Tuple[int, int]]
            List of possible new positions.
        """
        # We get all the adjacent cells:
        adjacent_cells = self.get_all_adjacent_cells(cell)

        if len(adjacent_cells) == 0:
            raise ValueError("Cell has no adjacent cells.")

        # We make sure there is exactly one connected neighbour to the cell that is occupied:
        neighbour = None
        neighbour_counter = 0

        for adjacent_cell in adjacent_cells:
            if self._cell_values[adjacent_cell]:
                if protein.are_neighbours(
                    amino_acid_coords[adjacent_cell],
                    amino_acid_coords[cell.coordinates],
                ):
                    if neighbour_counter == 0:
                        neighbour = adjacent_cell
                        neighbour_counter += 1
                    else:
                        raise ValueError(
                            "Cell has more than one occupied adjacent cell. This is not an end cell."
                        )
        if neighbour is None:
            raise ValueError("Cell has no occupied adjacent cell. No move is possible.")

        # An end move pivots the residue to a free position adjacent to its connected neighbour
        # We get all the adjacent cells of the neighbour:
        neighbour_adjacent_cells = self.get_all_adjacent_cells(Coordinates2D(neighbour))

        if len(neighbour_adjacent_cells) == 0:
            raise ValueError(
                "Neighbour has no adjacent cells. No end move is possible."
            )

        new_positions = []
        for neighbour_adjacent_cell in neighbour_adjacent_cells:
            if not self._cell_values[neighbour_adjacent_cell]:
                new_positions.append(neighbour_adjacent_cell)

        if len(new_positions) == 0:
            raise ValueError(
                "Neighbour has no free adjacent cells. No end move is possible."
            )

        return new_positions

    def compute_corner_moves(
        self,
        cell: Coordinates2D,
        protein: ProteinHP,
        amino_acid_coords: dict[Tuple[int, int], AminoAcidHP],
    ) -> list[Tuple[int, int]]:
        """Computes the corner moves for a given cell.

        Parameters
        ----------
        cell : Coordinates2D
            The cell to compute the corner moves for.
        protein : ProteinHP
            The studied protein.
        amino_acid_coords: dict[Tuple[int, int], AminoAcidHP]
            Dictionary containing the coordinates of the amino acids in the conformation.

        Returns
        -------
        list[Tuple[int, int]]
            List of possible new positions.
        """
        # We get all the adjacent cells:
        adjacent_cells = self.get_all_adjacent_cells(cell)

        if len(adjacent_cells) == 0:
            raise ValueError("Cell has no adjacent cells.")

        # We make sure there are exactly two connected neighbours:
        neighbours = []
        for adjacent_cell in adjacent_cells:
            if self._cell_values[adjacent_cell]:
                if protein.are_neighbours(
                    amino_acid_coords[adjacent_cell],
                    amino_acid_coords[cell.coordinates],
                ):
                    neighbours.append(adjacent_cell)

        if len(neighbours) != 2:
            raise ValueError(
                "Cell has more or less than two occupied adjacent cells. This is not a corner cell."
            )

        # A corner move pivots the residue to a free position adjacent to its connected neighbours
        # We get all the adjacent cells of the neighbours:
        neighbour_adjacent_cells = []
        for neighbour in neighbours:
            neighbour_adjacent_cells.append(
                self.get_all_adjacent_cells(Coordinates2D(neighbour))
            )

        # If one of the neighbours has no adjacent cells, no corner move is possible
        if len(neighbour_adjacent_cells) == 0:
            raise ValueError(
                "Neighbours have no adjacent cells. No corner move is possible."
            )

        for neighbour_adjacent_cell in neighbour_adjacent_cells:
            if len(neighbour_adjacent_cell) == 0:
                raise ValueError(
                    "Neighbour has no adjacent cells. No corner move is possible."
                )

        # We get the common adjacent cells of the neighbours:
        new_positions = []
        for coordinates in neighbour_adjacent_cells[0]:
            if coordinates in neighbour_adjacent_cells[1]:
                if not self._cell_values[coordinates]:
                    new_positions.append(coordinates)

        if len(new_positions) == 0:
            raise ValueError(
                "Neighbours have no free adjacent cells. No corner move is possible."
            )

        return new_positions

    def search_u_structures(self) -> dict[int, list[Coordinates2D]]:
        """Searches the lattice for U structures.

        Returns
        -------
        dict[int, list[Coordinates2D]]
            Dictionary containing the U structures found in the lattice. The keys are ids and the values are lists of
        """
        pass

    def compute_crankshaft_moves(self, cell: Coordinates2D) -> list[Coordinates2D]:
        pass

    def compute_pull_moves(self, cell: Coordinates2D) -> list[Coordinates2D]:
        pass
