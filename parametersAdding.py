#just starting to write in the parameters, can copy and paste later into something 

# --- Import libraries ---
from pyomo.environ import *
import math

# --- Create model ---
model = ConcreteModel()
#need to make a set 

#model.c is the set of chemical species 
model.c = Set(initialize=['CH3', 'CH2', 'CH', 'C', '=CH', '=C', 'NH2', 'NH', 'N', '=N', '-OH', 'O', '=O'])

#parameters 

#set for 313.15K cp group contributions in [J/K*mol], Rayer paper Table 4 GAA method
cp_313_data={'CH3': 43.56, 'CH2': 31.40, 'CH': 48.19, 'C': 10.04, '=CH': 21.24, '=C':-0.81, 'NH2':58.40, 
              'NH':42.27, 'N':2.00, '=N':34.81, 'OH':57.51, 'O':22.50, '=O':29.03
            }
model.cp313 = Param(model.c, initialize=cp_313_data)
#set for 393.15K cp group contributions in [J/K*mol] Rayer paper Table 4 GAA method
cp_393_data={'CH3': 58.02, 'CH2': 30.10, 'CH': 22.93, 'C': -6.66, '=CH': 38.17, '=C':-0.17, 'NH2':68.48, 
              'NH':57.87, 'N':9.68, '=N':13.04, 'OH':70.33, 'O':30.57, '=O':36.65
            }

model.cp393 = Param(model.c, initialize=cp_393_data)

#the calculation for cps
 


print("done")


'''
# Gibbs standard free energy
gibbs_data = {
    'H': -10.021, 'H2': -21.096, 'H2O': -37.986, 'N': -9.846, 'N2': -28.653,
    'NH': -18.918, 'NO': -28.032, 'O': -14.64, 'O2': -30.594, 'OH': -26.11
}
model.g = Param(model.c, initialize=gibbs_data)
#So now, inside Pyomo you can reference model.g['H2O'] and get -37.986

# Gibbs standard free energy + pressure term
gibbs_p_data = {c : gibbs_data[c] + math.log(750 * 0.07031) for c in gibbs_data.keys()}
model.g_p = Param(model.c, initialize=gibbs_p_data)
'''