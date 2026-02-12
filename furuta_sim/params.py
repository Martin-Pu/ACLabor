# -*- coding: utf-8 -*-

# ==========================================================
# Parameter aus Tabelle 1 (Labor – AC – Furuta-Pendel)
# Alle Größen in SI-Einheiten
# ==========================================================

# Physikalische Konstanten
g = 9.81                    # m/s^2

# Geometrie
L1 = 0.155                  # m   (155 mm)
L2 = 0.240                  # m   (240 mm)
lP = 0.1917                 # m   (191.7 mm)

# Massen
mS = 0.0274                 # kg  (Stabmasse, 27.4 g)
mP = 0 # 0.0510                 # kg  (Pendelmasse, 51.0 g)

# Trägheitsmomente
J1_hat = 2.89e-2            # kg*m^2

# Reibungsparameter Rotor
mu_V1 = 50.91e-2            # N*m*s
mu_H1 = 93.58e-3            # N*m

# Reibungsparameter Pendel
mu_V2 = 1.332e-5            # N*m*s
mu_H2 = 3.018e-4            # N*m

# Epsilon für Reibungsmodel
epsilon = 0.01

# Aktuatorgrenzen
tau_max = 4.0               # N*m

#Berechnete Werte aus Aufgabe 3.2
m2 = mP + mS
J0_hat = J1_hat + m2*lP
l2 = L2/2 # für pendelmasse mP = 0
J2_hat = 1/12*mS*L2**2 + ms*l2**2 + mP*lP**2
