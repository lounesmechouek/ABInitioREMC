from abc import ABC, abstractmethod
from dataclasses import dataclass

from .AminoAcid import AminoAcid


@dataclass(slots=True)
class Protein(ABC):
    """Protein is an abstract class that represents a protein regardless of the model used to describe it."""

    _name: str  # Name of the protein.
    _sequence: list[AminoAcid]  # Sequence of amino acids that compose the protein.

    @property
    def name(self) -> str:
        """Getter for the attribute name of the protein.

        Returns
        -------
        str
            Name of the protein.
        """
        return self._name

    @name.setter
    def name(self, name: str) -> None:
        """Setter for the attribute name of the protein.

        Parameters
        ----------
        name : str
            Name of the protein to be assigned.
        """
        self._name = name

    @property
    def sequence(self) -> list[AminoAcid]:
        """Getter for the attribute sequence of the protein.

        Returns
        -------
        list[AminoAcid]
            Sequence of amino acids that compose the protein.
        """
        return self._sequence

    @sequence.setter
    def sequence(self, sequence: list[AminoAcid]) -> None:
        """Setter for the attribute sequence of the protein.

        Parameters
        ----------
        sequence : list[AminoAcid]
            Sequence of amino acids that compose the protein to be assigned.
        """
        self._sequence = sequence

    @abstractmethod
    def protein_model(self) -> str:
        """Returns the type of the protein.

        Returns
        -------
        str
            Name of the protein model.
        """
        pass
