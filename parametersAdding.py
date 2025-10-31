#just starting to write in the parameters, can copy and paste later into something 

# --- Import libraries ---
from pyomo.environ import *
import math

# --- Create model ---
model = ConcreteModel()
#need to make a set 

model.c = Set(initialize=['CH3', 'CH2', 'CH', 'C', '=CH', '=C', 'NH2', 'NH', 'N', '=N', '-OH', 'O', '=O'])

#set for 313.15K cp group contributions in [J/K*mol], Rayer paper Table 4 GAA method
cp_313_data={'CH3': 43.56, 'CH2': 31.40, 'CH': 48.19, 'C': 10.04, '=CH': 21.24, '=C':-0.81, 'NH2':58.40, 
              'NH':42.27, 'N':2.00, '=N':34.81, 'OH':57.51, 'O':22.50, '=O':29.03
            }

#set for 393.15K cp group contributions in [J/K*mol] Rayer paper Table 4 GAA method
cp_393_data={'CH3': 58.02, 'CH2': 30.10, 'CH': 22.93, 'C': -6.66, '=CH': 38.17, '=C':-0.17, 'NH2':68.48, 
              'NH':57.87, 'N':9.68, '=N':13.04, 'OH':70.33, 'O':30.57, '=O':36.65
            }

print("done")