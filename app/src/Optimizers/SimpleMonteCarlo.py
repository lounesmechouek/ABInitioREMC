import copy
import math
import random

from ..Controllers.ConformationManager import ConformationManager
from ..Models.Conformation import Conformation


class SimpleMonteCarlo:
    """Class for Monte Carlo optimization algorithm in the AB-Initio context."""

    _phi: int  # Number of search steps.

    @property
    def phi(self) -> int:
        """Getter for the attribute phi of the MonteCarlo class.

        Returns
        -------
        int
            Number of search steps.
        """
        return self._phi

    @phi.setter
    def phi(self, phi: int) -> None:
        """Setter for the attribute phi of the MonteCarlo class.

        Parameters
        ----------
        phi : int
            Number of search steps to be assigned.
        """
        self._phi = phi

    def __init__(self, phi: int) -> None:
        """Constructor for the MonteCarlo class.

        Parameters
        ----------
        phi : int
            Number of search steps.
        rho : float, optional
            Probability to use pull moves, by default 0.0
        """
        self._phi = phi

    def optimize(
        self,
        conformation: Conformation,
        temperature: float,
        conf_manager: ConformationManager,
    ) -> Conformation:
        """Optimizes a conformation using the Monte Carlo algorithm.

        Parameters
        ----------
        conformation : Conformation
            Conformation to be optimized.
        temperature : float
            Temperature of the replica.
        conf_manager : ConformationManager
            Conformation manager that is used to compute the neigbourhood.

        Returns
        -------
        Conformation
            Optimized conformation.
        """
        optimal_conformation = copy.deepcopy(conformation)

        for i in range(self._phi):
            neighbourhood = []

            # We compute the neighbourhood of the conformation
            try:
                neighbourhood = conf_manager.compute_vhsd_neighbourhood(
                    optimal_conformation
                )
                # TODO : Compute pull moves.
            except Exception as e:
                raise e

            if len(neighbourhood) == 0:
                return optimal_conformation

            # We select a random conformation from the neighbourhood.
            random_conformation = copy.deepcopy(random.choice(neighbourhood))

            # We compute the energy of the conformations.
            try:
                random_conformation.compute_energy()
                optimal_conformation.compute_energy()
            except Exception as e:
                raise e

            if (
                random_conformation.computed_energy
                < optimal_conformation.computed_energy
            ):
                optimal_conformation = copy.deepcopy(random_conformation)
            else:
                q = random.uniform(0, 1)
                threshold = math.exp(
                    -(
                        random_conformation.computed_energy
                        - optimal_conformation.computed_energy
                    )
                    / temperature
                )
                if q > threshold:
                    optimal_conformation = copy.deepcopy(random_conformation)

        return optimal_conformation
