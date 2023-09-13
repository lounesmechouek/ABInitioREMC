from abc import ABC, abstractmethod

from ..src.Models.Conformation import Conformation


class ConformationDrawer(ABC):
    """Class used to represent protein conformations in HP-Model."""

    _conformation: Conformation  # Conformation to be drawn
    _colors: dict  # Colors of the amino acids (key: polarity "H"-"P", value: color)

    def __init__(self, conformation: Conformation, colors: dict) -> None:
        """Constructor for the ConformationDrawer class.

        Parameters
        ----------
        conformation : Conformation
            Conformation to be drawn.
        """
        self._conformation = conformation
        self._colors = colors

    @property
    def conformation(self) -> Conformation:
        """Getter for the attribute conformation of the conformation drawer.

        Returns
        -------
        Conformation
            Conformation to be drawn.
        """
        return self._conformation

    @conformation.setter
    def conformation(self, conformation: Conformation) -> None:
        """Setter for the attribute conformation of the conformation drawer.

        Parameters
        ----------
        conformation : Conformation
            Conformation to be drawn.
        """
        self._conformation = conformation

    @property
    def colors(self) -> dict:
        """Getter for the attribute colors of the conformation drawer.

        Returns
        -------
        dict
            Colors of the amino acids (key: polarity "H"-"P", value: color)
        """
        return self._colors

    @colors.setter
    def colors(self, colors: dict) -> None:
        """Setter for the attribute colors of the conformation drawer.

        Parameters
        ----------
        colors : dict
            Colors of the amino acids (key: polarity "H"-"P", value: color)
        """
        self._colors = colors

    @abstractmethod
    def draw(self) -> None:
        """Draws the conformation."""
        pass
