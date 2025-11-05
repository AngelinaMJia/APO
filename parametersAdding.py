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

#here maybe a group count, like 
#group_count = {'CH3': 1, 'CH2': 1, 'CH': 1, etc}

#the calculation for cps
#do the summation as a loop 
model.n_groups = Param(model.c, initialize={g: group_count.get(g, 0) for g in model.c})

#define both as constraints since calculated and is not known before optimization begins 
model.cp313_constraint_rule = Var(bounds = (0, None))
model.cp393_constraint_rule = Var(bounds = (0, None))

#from the "pop-quiz equations" slide 
def cp313_constraint_rule(model):
    return model.cp_313_total == sum(model.n_groups[g] * model.cp313[g] for g in model.c)
model.cp313_constraint = Constraint(rule=cp313_constraint_rule)

def cp393_constraint_rule(model):
    return model.cp_393_total == sum(model.n_groups[g] * model.cp393[g] for g in model.c)
model.cp393_constraint = Constraint(rule=cp393_constraint_rule)


print(f"Cp_313 (J/mol·K): {value(model.cp_313_total):.2f}")
print(f"Cp_393 (J/mol·K): {value(model.cp_393_total):.2f}")
print("done")
