import pyomo.environ as pe

model = pe.ConcreteModel(name="Solvent_CAMD_Step1")

# define 1st-order groups
model.G1 = pe.Set(initialize=[
    'CH3', 'CH2', 'CH', 'C', 'NH2', 'NH', 'N', 'OH'
])

# define 2nd-order groups
model.G2 = pe.Set(initialize=[
    'CH2-OH', 'CH-OH', 'C-OH',
    'CH2-NH2', 'CH-NH2',
    'CH2-NH', 'CH-NH',
])

# universal constants (boiling point, melting point, molar volume constant)
model.Tb0 = pe.Param(initialize=244.7889)