import json

from app.src.Models.Protein import Protein

from ..Exceptions.EmptyFile import EmptyFile
from ..Models.AminoAcidHP import AminoAcidHP
from ..Models.Polarity import Polarity
from ..Models.Protein import Protein
from ..Models.ProteinHP import ProteinHP
from ..Models.ProteinModel import ProteinModel
from .ProteinIO import ProteinIO


class JSONProteinIO(ProteinIO):
    """This is an interface used for reading and writing protein data from JSON files."""

    _filename: str

    def __init__(self, filename: str) -> None:
        """Constructor for the JSONProteinIO class.

        Parameters
        ----------
        filename : str
            Name of the JSON file to be read from and written to.
        """
        self._filename = filename

    @property
    def filename(self) -> str:
        """Getter for the attribute filename of the JSONProteinIO.

        Returns
        -------
        str
            Name of the JSON file to be read from and written to.
        """
        return self._filename

    @filename.setter
    def filename(self, filename: str) -> None:
        """Setter for the attribute filename of the JSONProteinIO.

        Parameters
        ----------
        filename : str
            Name of the JSON file to be read from and written to.
        """
        self._filename = filename

    def _get_json_data(self) -> dict:
        """Gets the data from the JSON file.

        Returns
        -------
        dict
            Data from the JSON file.
        """
        try:
            with open(self._filename, "r") as file:
                data = json.load(file)
                if data is None:
                    raise EmptyFile(f"No content was found in {self._filename}.")
                return data
        except FileNotFoundError:
            raise FileNotFoundError(f"File {self._filename} was not found.")
        except json.JSONDecodeError as e:
            raise Exception(f"File {self._filename} is not a valid JSON file.") from e
        except Exception as e:
            raise Exception(
                f"An unknown error occurred while reading {self._filename}."
            ) from e
        finally:
            if "file" in locals():
                file.close()

    def _save_json_data(self, json_data: dict) -> None:
        """Saves data to the JSON file.

        Parameters
        ----------
        json_data : dict
            Data to be saved to the JSON file.
        """
        try:
            with open(self._filename, "w") as file:
                json.dump(json_data, file, indent=4)
        except Exception as e:
            raise e

    def _json_to_proteins_hp(self, json_data: dict) -> list[ProteinHP]:
        """Converts a JSON object to a list of proteins.

        Parameters
        ----------
        json_data : dict
            JSON object to be converted.

        Returns
        -------
        list[ProteinHP]
            List of proteins converted from a JSON object.
        """
        proteins = []
        for protein in json_data:
            seq_amino_acids = []
            for amino_acid in protein["sequence"]:
                pol = None
                if amino_acid["polarity"] == "H":
                    pol = Polarity.HYDROPHOBIC
                elif amino_acid["polarity"] == "P":
                    pol = Polarity.POLAR
                else:
                    raise Exception(
                        f"Invalid polarity {amino_acid['polarity']} at protein {protein['name']}. Please make sure that the polarity is either H or P."
                    )

                try:
                    seq_amino_acids.append(
                        AminoAcidHP(
                            amino_acid["name"],
                            amino_acid["abbreviation"],
                            pol,
                        )
                    )
                except Exception as e:
                    # TODO : Check for the specific exceptions that can be raised during the creation of an AminoAcidHP object.
                    raise e

            try:
                proteins.append(
                    ProteinHP(
                        protein["name"],
                        seq_amino_acids,
                        int(protein["e_star"]),
                        int(protein["recommended_dimension"]),
                    )
                )
            except Exception as e:
                # TODO : Check for the specific exceptions that can be raised during the creation of a ProteinHP object.
                raise e

        return proteins

    def load_proteins(self, protein_model: ProteinModel) -> list[Protein]:
        """Loads proteins from a JSON file.

        Parameters
        ----------
        protein_model : ProteinModel
            Model of the protein to be loaded.

        Returns
        -------
        list[Protein]
            List of proteins loaded from a JSON file.
        """
        try:
            protein_data = self._get_json_data()
        except Exception as e:
            raise e

        if protein_model == ProteinModel.HYDROPHOBIC_POLAR:
            try:
                proteins = self._json_to_proteins_hp(protein_data)
                return proteins
            except Exception as e:
                raise e
        else:
            raise Exception(f"Protein model {protein_model} is not supported.")

    def _proteins_hp_to_json(self, proteins: list[ProteinHP]) -> dict:
        """Converts a list of proteins to a JSON object.

        Parameters
        ----------
        proteins : list[ProteinHP]
            List of proteins to be converted.

        Returns
        -------
        dict
            JSON object converted from a list of proteins.
        """
        json_data = []
        for protein in proteins:
            seq = []
            for amino_acid in protein.sequence:
                pol = None
                if amino_acid.polarity == Polarity.HYDROPHOBIC:
                    pol = "H"
                else:
                    pol = "P"
                seq.append(
                    {
                        "name": amino_acid.name,
                        "abbreviation": amino_acid.abbreviation,
                        "polarity": pol,
                    }
                )
            json_data.append(
                {
                    "name": protein.name,
                    "sequence": seq,
                    "e_star": protein.e_star,
                    "recommended_dimension": protein.recommended_dimension,
                }
            )
        return json_data

    def save_proteins(
        self, proteins: list[Protein], protein_model: ProteinModel
    ) -> None:
        """Saves proteins to a JSON file.

        Parameters
        ----------
        proteins : list[Protein]
            List of proteins to be saved to a JSON file.
        protein_model : ProteinModel
            Model of the protein to be saved.
        """
        if protein_model == ProteinModel.HYDROPHOBIC_POLAR:
            try:
                json_data = self._proteins_hp_to_json(proteins)
                self._save_json_data(json_data)
            except Exception as e:
                raise e
        else:
            raise Exception(f"Protein model {protein_model} is not supported.")
