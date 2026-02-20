class FurutaParams:
    def __init__(self, mP=0.0):
        # Konstanten
        self.g = 9.81

        # Geometrie
        self._L1 = 0.155
        self._L2 = 0.240
        self._lP = 0.1917

        # Massen
        self._mS = 0.0274
        self._mP = mP

        # Trägheitsmomente (Basis)
        self._J1_hat = 2.89e-2

        # Reibungsparameter
        self.mu_V1 = 50.91e-2
        self.mu_H1 = 93.58e-3
        self.mu_V2 = 1.332e-5
        self.mu_H2 = 3.018e-4

        self.epsilon = 0.01
        self.tau_max = 4.0

        self._recompute()

    # ------------------
    # Properties (inputs)
    # ------------------
    @property
    def L1(self):
        return self._L1

    @L1.setter
    def L1(self, value):
        self._L1 = value
        self._recompute()

    @property
    def L2(self):
        return self._L2

    @L2.setter
    def L2(self, value):
        self._L2 = value
        self._recompute()

    @property
    def lP(self):
        return self._lP

    @lP.setter
    def lP(self, value):
        self._lP = value
        self._recompute()

    @property
    def mS(self):
        return self._mS

    @mS.setter
    def mS(self, value):
        self._mS = value
        self._recompute()

    @property
    def mP(self):
        return self._mP

    @mP.setter
    def mP(self, value):
        self._mP = value
        self._recompute()

    @property
    def J1_hat(self):
        return self._J1_hat

    @J1_hat.setter
    def J1_hat(self, value):
        self._J1_hat = value
        self._recompute()

    # ------------------
    # Derived parameters
    # ------------------
    def _recompute(self):
        self.m2 = self._mS + self._mP
        self.J0_hat = self._J1_hat + self.m2 * self._L1**2

        self.l2 = self._L2 / 2.0
        self.J2_hat = (
            (1 / 12) * self._mS * self._L2**2
            + self._mS * self.l2**2
            + self._mP * self._lP**2
        )