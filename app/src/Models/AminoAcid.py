from abc import ABC, abstractmethod
from dataclasses import dataclass


@dataclass(slots=True)
class AminoAcid(ABC):
    """AminoAcid is an abstract class that represents an amino acid regardless of the model used to describe it."""

    _id: int  # ID of the amino acid
    _name: str  # Name of the amino acid. Ex: Alanine
    _abbreviation: str  # Abbreviation of the amino acid (3 letters). Ex: Ala

    @property
    def id(self) -> int:
        """Getter for the attribute id of the amino acid.

        Returns
        -------
        int
            ID of the amino acid.
        """
        return self._id

    @id.setter
    def id(self, id: int) -> None:
        """Setter for the attribute id of the amino acid.

        Parameters
        ----------
        id : int
            ID of the amino acid to be assigned.
        """
        self._id = id

    @property
    def name(self) -> str:
        """Getter for the attribute name of the amino acid.

        Returns
        -------
        str
            Name of the amino acid.
        """
        return self._name

    @name.setter
    def name(self, name: str) -> None:
        """Setter for the attribute name of the amino acid.

        Parameters
        ----------
        name : str
            Name of the amino acid to be assigned.
        """
        self._name = name

    @property
    def abbreviation(self) -> str:
        """Getter for the attribute abbreviation of the amino acid.

        Returns
        -------
        str
            Abbreviation of the amino acid.
        """
        return self._abbreviation

    @abbreviation.setter
    def abbreviation(self, abbreviation: str) -> None:
        """Setter for the attribute abbreviation of the amino acid.

        Parameters
        ----------
        abbreviation : str
            Abbreviation of the amino acid to be assigned.
        """
        self._abbreviation = abbreviation

    @abstractmethod
    def acide_model(self) -> str:
        """Returns the type of the amino acid.

        Returns
        -------
        str
            Type of the amino acid.
        """
        pass
