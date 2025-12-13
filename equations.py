import pyomo.environ as pe
import math

# Model Creation ---
model = pe.ConcreteModel(name="Solvent_CAMD")

# Defining Sets ---

# --- (Loris) Set of 1st-Order Groups (G1_Hukkerikar) ---
model.G1_Hukkerikar = pe.Set(initialize=[
    # Alkyl Groups
    'CH3',     # Group #1
    'CH2',     # Group #2
    'CH',      # Group #3
    'C',       # Group #4
    
    # Primary Amines
    'CH2NH2',  # Group #54
    'CHNH2',   # Group #55
    'CNH2',    # Group #56 (C-NH2)
    
    # Secondary Amines
    'CH3NH',   # Group #57
    'CH2NH',   # Group #58
    'CHNH',    # Group #59
    
    # Tertiary Amines
    'CH3N',    # Group #60
    'CH2N',    # Group #61
    
    # Hydroxyl Group
    'OH',      # Group #29 (Terminal -OH group)
    
    # Groups to Cover Heat Capacity Method Groups
    'CH2O',    # Group #50 ('O')
    'CH3CO',   # Group #33 ('=O')
    'CH2CO',   # Group #34 ('=O')
    'CH2=CH',  # Group #5  ('=CH')
    'CH=CH',   # Group #6  ('=CH', '=C')
    'CH=N',    # Group #66 ('=N')
])

#limit on number of unsaturated bonds
unsaturated= ['CH2=CH', 'CH=CH', 'CH=N']


# --- Heat Capacity Group Set ---
# ONLY for heat capacity equations, groups differ for other method
model.Cp_Groups = pe.Set(initialize=[
    'CH3', 'CH2', 'CH', 'C', '=CH', '=C', 'NH2', 'NH', 'N', '=N', 'OH', 'O', '=O'
])

print("--- All Sets Defined ---")
# print statements are safe here because sets are defined
# model.G1_Hukkerikar.pprint() 


# Defining Parameters ---
UB_num_groups= 15 #upper bound on the number of groups

# --- Heat Capacity Cp Parameters ---
cp_313_data={'CH3': 43.56, 'CH2': 31.40, 'CH': 48.19, 'C': 10.04, '=CH': 21.24, '=C':-0.81, 'NH2':58.40,
             'NH':42.27, 'N':2.00, '=N':34.81, 'OH':57.51, 'O':22.50, '=O':29.03}
model.cp313 = pe.Param(model.Cp_Groups, initialize=cp_313_data)

cp_393_data={'CH3': 58.02, 'CH2': 30.10, 'CH': 22.93, 'C': -6.66, '=CH': 38.17, '=C':-0.17, 'NH2':68.48,
             'NH':57.87, 'N':9.68, '=N':13.04, 'OH':70.33, 'O':30.57, '=O':36.65}
model.cp393 = pe.Param(model.Cp_Groups, initialize=cp_393_data)


# --- (Loris) Universal Constants (from Hukkerikar, Table 5) ---
model.Tb0 = pe.Param(initialize=244.7889, doc="Boiling point constant (K)")
model.Tm0 = pe.Param(initialize=144.0977, doc="Melting point constant (K)")
model.Vm0 = pe.Param(initialize=0.0123, doc="Molar volume constant (m^3/kmol)")

# --- (Loris) External Data (CO2 Hansen Parameters http://kinampark.com/PL/files/Books/Hansen%20Solubility%20Parameters%202000.pdf)
# Need to double check these as I don't even think they're consistent within the one report I got them from
model.delta_D_CO2 = pe.Param(initialize=14.7, doc="Hansen D for CO2 (MPa^0.5)")
model.delta_P_CO2 = pe.Param(initialize=3.9, doc="Hansen P for CO2 (MPa^0.5)")
model.delta_H_CO2 = pe.Param(initialize=6.7, doc="Hansen H for CO2 (MPa^0.5)")

# --- (Loris) Structural Parameters (for G1_Hukkerikar) ---
valency_data = {
    'CH3': 1, 'CH2': 2, 'CH': 3, 'C': 4,
    'CH2NH2': 1, 'CHNH2': 2, 'CNH2': 3,
    'CH3NH': 1, 'CH2NH': 2, 'CHNH': 3,
    'CH3N': 2, 'CH2N': 3,
    'OH': 1,
    'CH2O': 2, 'CH3CO': 1, 'CH2CO': 2, 'CH2=CH': 1, 'CH=CH': 2, 'CH=N': 1
}
mw_data_kg_kmol = {
    'CH3': 15.035, 'CH2': 14.026, 'CH': 13.018, 'C': 12.011,
    'CH2NH2': 30.059, 'CHNH2': 29.051, 'CNH2': 28.043,
    'CH3NH': 29.051, 'CH2NH': 28.043, 'CHNH': 27.035,
    'CH3N': 43.078, 'CH2N': 42.070,
    'OH': 17.008,
    'CH2O': 30.026, 'CH3CO': 43.042, 'CH2CO': 42.036, 'CH2=CH': 27.043, 'CH=CH': 26.037, 'CH=N': 27.025
}
for g1 in model.G1_Hukkerikar:
    valency_data.setdefault(g1, 2)
    mw_data_kg_kmol.setdefault(g1, 30.0)

model.valency = pe.Param(model.G1_Hukkerikar, initialize=valency_data, doc="Valency of G1 groups")
model.MW = pe.Param(model.G1_Hukkerikar, initialize=mw_data_kg_kmol, doc="Molecular weight of G1 groups (kg/kmol)")


# --- (Loris) G1_Hukkerikar Contributions (C_i) ---
C_Tb_data = {'CH3': 0.9218, 'CH2': 0.5780, 'CH': -0.1189, 'C': -0.6495, 'CH2NH2': 2.2640, 'CHNH2': 1.4372, 'CNH2': 0.8863, 'CH3NH': 1.9860, 'CH2NH': 1.2690, 'CHNH': 0.5940, 'CH3N': 0.9990, 'CH2N': 0.3324, 'OH': 2.2476, 'CH2O': 0.9750, 'CH3CO': 2.6907, 'CH2CO': 1.9665, 'CH2=CH': 1.4953, 'CH=CH': 1.2001, 'CH=N': 1.1099}
C_Tm_data = {'CH3': 0.7555, 'CH2': 0.2966, 'CH': -0.5960, 'C': -0.3679, 'CH2NH2': 3.3490, 'CHNH2': 36.2974, 'CNH2': 11.5011, 'CH3NH': 2.7394, 'CH2NH': 2.0378, 'CHNH': 1.3226, 'CH3N': 0.8482, 'CH2N': -0.4084, 'OH': 3.2424, 'CH2O': 0.6741, 'CH3CO': 3.2535, 'CH2CO': 2.8589, 'CH2=CH': 1.0430, 'CH=CH': 0.6600, 'CH=N': 9.2631}
C_Vm_data = {'CH3': 0.0238, 'CH2': 0.0166, 'CH': 0.0084, 'C': -0.0015, 'CH2NH2': 0.0262, 'CHNH2': 0.0214, 'CNH2': 0.0144, 'CH3NH': 0.0279, 'CH2NH': 0.0246, 'CHNH': 0.0182, 'CH3N': 0.0265, 'CH2N': 0.0190, 'OH': 0.0042, 'CH2O': 0.0228, 'CH3CO': 0.0347, 'CH2CO': 0.0283, 'CH2=CH': 0.0333, 'CH=CH': 0.0244, 'CH=N': 0.0}
C_dD_data = {'CH3': 7.5983, 'CH2': -0.0023, 'CH': -7.5390, 'C': -15.6455, 'CH2NH2': 8.1995, 'CHNH2': -0.3812, 'CNH2': 0, 'CH3NH': 7.7307, 'CH2NH': 0.0223, 'CHNH': -7.5377, 'CH3N': 0.3494, 'CH2N': -6.7232, 'OH': 8.0503, 'CH2O': 0.1706, 'CH3CO': 8.1107, 'CH2CO': 0.5371, 'CH2=CH': 7.7504, 'CH=CH': 0.4284, 'CH=N': 0.0}
C_dP_data = {'CH3': 2.3037, 'CH2': -0.1664, 'CH': -3.3851, 'C': -5.1979, 'CH2NH2': 5.2101, 'CHNH2': 0.5616, 'CNH2': 0, 'CH3NH': 2.5065, 'CH2NH': -0.7159, 'CHNH': -4.6694, 'CH3N': -0.4250, 'CH2N': -0.7354, 'OH': 5.2379, 'CH2O': 0.5137, 'CH3CO': 6.3823, 'CH2CO': 1.2706, 'CH2=CH': 3.6752, 'CH=CH': 3.0492, 'CH=N': 0.0}
C_dH_data = {'CH3': 2.2105, 'CH2': -0.2150, 'CH': -2.6826, 'C': -6.4821, 'CH2NH2': 6.7984, 'CHNH2': 2.8953, 'CNH2': 0, 'CH3NH': 7.2551, 'CH2NH': 1.4183, 'CHNH': -2.2824, 'CH3N': 2.4585, 'CH2N': -7.3014, 'OH': 11.8005, 'CH2O': 0.8246, 'CH3CO': 3.4394, 'CH2CO': -0.0788, 'CH2=CH': 2.7673, 'CH=CH': 0.8631, 'CH=N': 0.0}
# NO VALUES IN SOURCE FOR DISPERSION PARAMETERS OF CNH2 OR CH=N, SET AS 0 FOR NOW but need to either remove or estimate (from discussion board)

for g1 in model.G1_Hukkerikar:
    C_Tb_data.setdefault(g1, 1.0); C_Tm_data.setdefault(g1, 1.0); C_Vm_data.setdefault(g1, 0.015)
    C_dD_data.setdefault(g1, 0.0); C_dP_data.setdefault(g1, 0.0); C_dH_data.setdefault(g1, 0.0)

model.C_Tb = pe.Param(model.G1_Hukkerikar, initialize=C_Tb_data)
model.C_Tm = pe.Param(model.G1_Hukkerikar, initialize=C_Tm_data)
model.C_Vm = pe.Param(model.G1_Hukkerikar, initialize=C_Vm_data)
model.C_dD = pe.Param(model.G1_Hukkerikar, initialize=C_dD_data)
model.C_dP = pe.Param(model.G1_Hukkerikar, initialize=C_dP_data)
model.C_dH = pe.Param(model.G1_Hukkerikar, initialize=C_dH_data)

# Variables

# Decision Variables
model.N = pe.Var(model.G1_Hukkerikar, domain=pe.NonNegativeIntegers, bounds=(0, 10), doc="Count of 1st-order groups")

# Calculated Property Variables
model.Tb_calc = pe.Var(domain=pe.Reals, bounds=(200, 1000))
model.Tm_calc = pe.Var(domain=pe.Reals, bounds=(100, 600))
model.rho_calc = pe.Var(domain=pe.Reals, bounds=(500, 2000))
model.RED_calc = pe.Var(domain=pe.Reals, bounds=(0, 10000))
model.Cp_313_calc = pe.Var(bounds=(0, None))
model.Cp_393_calc = pe.Var(bounds=(0, None))

# Linking Groups from Both Reports (Need to double check this is complete/correct)
def N_Cp_mapping_rule(model, cp_group):
    if cp_group == 'CH3':
        return model.N['CH3'] + model.N['CH3NH'] + model.N['CH3N'] + model.N['CH3CO']
    if cp_group == 'CH2':
        return model.N['CH2'] + model.N['CH2NH2'] + model.N['CH2NH'] + model.N['CH2N'] + model.N['CH2O'] + model.N['CH2CO']
    if cp_group == 'CH':
        return model.N['CH'] + model.N['CHNH2'] + model.N['CHNH']
    if cp_group == 'C':
        return model.N['C'] + model.N['CNH2']
    if cp_group == 'NH2':
        return model.N['CH2NH2'] + model.N['CHNH2'] + model.N['CNH2']
    if cp_group == 'NH':
        return model.N['CH3NH'] + model.N['CH2NH'] + model.N['CHNH']
    if cp_group == 'N':
        return model.N['CH3N'] + model.N['CH2N']
    if cp_group == 'OH':
        return model.N['OH']
    if cp_group == 'O': 
        return model.N['CH2O']
    if cp_group == '=O': 
        return model.N['CH3CO'] + model.N['CH2CO']
    if cp_group == '=CH':
        return model.N['CH2=CH'] + 2 * model.N['CH=CH'] + model.N['CH=N']
    if cp_group == '=C':
        return model.N['CH2=CH']
    if cp_group == '=N':
        return model.N['CH=N']
    return 0
model.N_Cp = pe.Expression(model.Cp_Groups, rule=N_Cp_mapping_rule)


# Constraints

# Heat Capacity
def cp313_calc_rule(model):
    return model.Cp_313_calc == sum(model.N_Cp[i] * model.cp313[i] for i in model.Cp_Groups)
model.cp313_sum = pe.Constraint(rule=cp313_calc_rule)

def cp393_calc_rule(model):
    return model.Cp_393_calc == sum(model.N_Cp[i] * model.cp393[i] for i in model.Cp_Groups)
model.cp393_sum = pe.Constraint(rule=cp393_calc_rule)

# Structural Feasibility
# Octet Rule for Acyclic Chains: sum(N*(v-2)) = -2
model.c_octet = pe.Constraint(
    expr = sum((model.valency[i] - 2) * model.N[i] for i in model.G1_Hukkerikar) == -2
)

# Bonding Rule
def bonding_rule(model, j):
    total_groups = sum(model.N[i] for i in model.G1_Hukkerikar)
    # N_j * (Valency_j - 1) + 2 <= Total_Number_of_Groups
    return model.N[j] * (model.valency[j] - 1) + 2 <= total_groups
model.c_bonding = pe.Constraint(model.G1_Hukkerikar, rule=bonding_rule)

# Upper Bound on size
model.c_num_of_groups = pe.Constraint(
    expr = sum(model.N[i] for i in model.G1_Hukkerikar) <= UB_num_groups
)

#constaint on saturated groups
model.max_unsat = pe.Constraint(
    expr = sum(model.N[g] for g in unsaturated) <= 1
)



# Property Calculations
small_epsilon = 1e-6

model.f_Tb_expr = pe.Expression(expr=sum(model.N[i] * model.C_Tb[i] for i in model.G1_Hukkerikar))
model.f_Tm_expr = pe.Expression(expr=sum(model.N[i] * model.C_Tm[i] for i in model.G1_Hukkerikar))

model.c_Tb_calc = pe.Constraint(expr= model.Tb_calc == model.Tb0 * pe.log(model.f_Tb_expr + small_epsilon))
model.c_Tm_calc = pe.Constraint(expr= model.Tm_calc == model.Tm0 * pe.log(model.f_Tm_expr + small_epsilon))

model.MW_calc_expr = pe.Expression(expr=sum(model.N[i] * model.MW[i] for i in model.G1_Hukkerikar))
model.Vm_calc_expr = pe.Expression(expr=model.Vm0 + sum(model.N[i] * model.C_Vm[i] for i in model.G1_Hukkerikar))

model.c_rho_calc = pe.Constraint(expr= model.rho_calc * model.Vm_calc_expr == model.MW_calc_expr)

model.f_dD_expr = pe.Expression(expr=sum(model.N[i] * model.C_dD[i] for i in model.G1_Hukkerikar))
model.f_dP_expr = pe.Expression(expr=sum(model.N[i] * model.C_dP[i] for i in model.G1_Hukkerikar))
model.f_dH_expr = pe.Expression(expr=sum(model.N[i] * model.C_dH[i] for i in model.G1_Hukkerikar))

dD_term = model.f_dD_expr - model.delta_D_CO2
dP_term = model.f_dP_expr - model.delta_P_CO2
dH_term = model.f_dH_expr - model.delta_H_CO2
Ra_squared = 4*(dD_term**2) + (dP_term**2) + (dH_term**2)
model.c_RED_calc = pe.Constraint(expr= (model.RED_calc)**2 == Ra_squared)
#Need to divide this by R0_CO2 but I am not sure where to get value from? I'll find at some point.
#It would multiply with RED_calc, equation written backwards essentially to avoid square rooting. 

# Temperature Bounds
model.c_Tb_min = pe.Constraint(expr= model.Tb_calc >= 393.0)
model.c_Tm_max = pe.Constraint(expr= model.Tm_calc <= 313.0)


# Objective Function
def objective_rule(m):
    return 1.0 * model.RED_calc + 1.0 * model.Cp_313_calc + 1.0 * model.Cp_393_calc - 0.001 * model.rho_calc
model.obj = pe.Objective(rule=objective_rule, sense=pe.minimize)

print("\nModel built successfully!")

solver = pe.SolverFactory('gams')
solver.options['solver'] = 'baron'
print("\nStarting solver...")
try:
    results = solver.solve(model, tee=True)
    print("Termination Condition:", results.solver.termination_condition)
    
    # Print Results
    print("\n--- Optimal Molecule ---")
    for i in model.G1_Hukkerikar:
        if pe.value(model.N[i]) > 0.5:
            print(f"{i}: {pe.value(model.N[i])}")
            
    print("\n--- Properties ---")
    print(f"Tb: {pe.value(model.Tb_calc):.2f} K")
    print(f"Tm: {pe.value(model.Tm_calc):.2f} K")
    print(f"Density: {pe.value(model.rho_calc):.2f} kg/m3")
    print(f"RED: {pe.value(model.RED_calc):.2f}")

except Exception as e:
    print(f"Solver failed: {e}")