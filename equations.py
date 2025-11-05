import pyomo.environ as pe
import math

# Model Creation
model = pe.ConcreteModel(name="Solvent_CAMD")

# Defining Sets
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

# --- (Loris) Set of 2nd-Order Groups (G2_Hukkerikar) ---
model.G2_Hukkerikar = pe.Set(initialize=[
    '(CH3)2CH',  # Correction for 2x CH3 + 1x CH
    '(CH3)3C',   # Correction for 3x CH3 + 1x C
    'CHOH',      # Correction for 1x CH + 1x OH
    'COH',       # Correction for 1x C + 1x OH
])

# --- Heat Capacity Group Set ---
# ONLY for heat capacity equations, groups differ for other method
model.Cp_Groups = pe.Set(initialize=[
    'CH3', 'CH2', 'CH', 'C', '=CH', '=C', 'NH2', 'NH', 'N', '=N', 'OH', 'O', '=O'
])

print("All Sets Defined")
print("G1_Hukkerikar (Decision Variables):")
model.G1_Hukkerikar.pprint()
print("G2_Hukkerikar (Correction Groups):")
model.G2_Hukkerikar.pprint()
print("Cp_Groups (for Heat Capacity):")
model.Cp_Groups.pprint()


# Defining Parameters

# --- Heat Capacity Cp Parameters ---
cp_313_data={'CH3': 43.56, 'CH2': 31.40, 'CH': 48.19, 'C': 10.04, '=CH': 21.24, '=C':-0.81, 'NH2':58.40,
             'NH':42.27, 'N':2.00, '=N':34.81, 'OH':57.51, 'O':22.50, '=O':29.03
            }
model.cp313 = pe.Param(model.Cp_Groups, initialize=cp_313_data)
cp_393_data={'CH3': 58.02, 'CH2': 30.10, 'CH': 22.93, 'C': -6.66, '=CH': 38.17, '=C':-0.17, 'NH2':68.48,
             'NH':57.87, 'N':9.68, '=N':13.04, 'OH':70.33, 'O':30.57, 'O=':36.65
            }
model.cp393 = pe.Param(model.Cp_Groups, initialize=cp_393_data)


# --- (Loris) Universal Constants (from Hukkerikar, Table 5) ---
model.Tb0 = pe.Param(initialize=244.7889, doc="Boiling point constant (K)")
model.Tm0 = pe.Param(initialize=144.0977, doc="Melting point constant (K)")
model.Vm0 = pe.Param(initialize=0.0123, doc="Molar volume constant (m^3/kmol)")

# --- (Loris) External Data (CO2 Hansen Parameters - NEED TO FIND ACTUAL VALUES) ---
model.delta_D_CO2 = pe.Param(initialize=15.0, doc="Hansen D for CO2 (MPa^0.5)")
model.delta_P_CO2 = pe.Param(initialize=5.0, doc="Hansen P for CO2 (MPa^0.5)")
model.delta_H_CO2 = pe.Param(initialize=6.0, doc="Hansen H for CO2 (MPa^0.5)")

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
    'CH3N': 43.078, 'CH2N': 42.070, # Note: CH3N is (CH3)-N<, CH2N is >CH2-N<
    'OH': 17.008,
    'CH2O': 30.026, 'CH3CO': 43.042, 'CH2CO': 42.036, 'CH2=CH': 27.043, 'CH=CH': 26.037, 'CH=N': 27.025
}
model.valency = pe.Param(model.G1_Hukkerikar, initialize=valency_data, doc="Valency of G1 groups")
model.MW = pe.Param(model.G1_Hukkerikar, initialize=mw_data_kg_kmol, doc="Molecular weight of G1 groups (kg/kmol)")


# --- (Loris) G1_Hukkerikar Contributions (C_i) for Tb (from Table S5) ---
C_Tb_data = {
    'CH3': 0.9218, 'CH2': 0.5780, 'CH': -0.1189, 'C': -0.6495,
    'CH2NH2': 2.2640, 'CHNH2': 1.4372, 'CNH2': 0.8863,
    'CH3NH': 1.9860, 'CH2NH': 1.2690, 'CHNH': 0.5940,
    'CH3N': 0.9990, 'CH2N': 0.3324,
    'OH': 2.2476,
    'CH2O': 0.9750, 'CH3CO': 2.6907, 'CH2CO': 1.9665,
    'CH2=CH': 1.4953, 'CH=CH': 1.2001, 'CH=N': 1.1099
}
model.C_Tb = pe.Param(model.G1_Hukkerikar, initialize=C_Tb_data, doc="Boiling T 1st-order contributions")

# --- (Loris) G1_Hukkerikar Contributions (C_i) for Tm (from Table S5) ---
C_Tm_data = {
    'CH3': 0.7555, 'CH2': 0.2966, 'CH': -0.5960, 'C': -0.3679,
    'CH2NH2': 3.3490, 'CHNH2': 36.2974, 'CNH2': 11.5011,
    'CH3NH': 2.7394, 'CH2NH': 2.0378, 'CHNH': 1.3226,
    'CH3N': 0.8482, 'CH2N': -0.4084,
    'OH': 3.2424,
    'CH2O': 0.6741, 'CH3CO': 3.2535, 'CH2CO': 2.8589,
    'CH2=CH': 1.0430, 'CH=CH': 0.6600, 'CH=N': 9.2631
}
model.C_Tm = pe.Param(model.G1_Hukkerikar, initialize=C_Tm_data, doc="Melting T 1st-order contributions")

# --- (Loris) G1_Hukkerikar Contributions (C_i) for Vm (from Table S5) ---
C_Vm_data = {
    'CH3': 0.0238, 'CH2': 0.0166, 'CH': 0.0084, 'C': -0.0015,
    'CH2NH2': 0.0262, 'CHNH2': 0.0214, 'CNH2': 0.0144,
    'CH3NH': 0.0279, 'CH2NH': 0.0246, 'CHNH': 0.0182,
    'CH3N': 0.0265, 'CH2N': 0.0190,
    'OH': 0.0042,
    'CH2O': 0.0228, 'CH3CO': 0.0347, 'CH2CO': 0.0283,
    'CH2=CH': 0.0333, 'CH=CH': 0.0244, 'CH=N': 0.0  # CH=N not in Table S5 Vm column
}
model.C_Vm = pe.Param(model.G1_Hukkerikar, initialize=C_Vm_data, doc="Volume 1st-order contributions")

# --- (Loris) G1_Hukkerikar Contributions (C_i) for dD (from Table S5) ---
C_dD_data = {
    'CH3': 7.5983, 'CH2': -0.0023, 'CH': -7.5390, 'C': -15.6455,
    'CH2NH2': 8.1995, 'CHNH2': -0.3812, 'CNH2': -6.8570,
    'CH3NH': 7.7307, 'CH2NH': 0.0223, 'CHNH': -7.5377,
    'CH3N': 0.3494, 'CH2N': -6.7232,
    'OH': 8.0503,
    'CH2O': 0.1706, 'CH3CO': 8.1107, 'CH2CO': 0.5371,
    'CH2=CH': 7.7504, 'CH=CH': 0.4284, 'CH=N': 0.0 # Not in table S5
}
model.C_dD = pe.Param(model.G1_Hukkerikar, initialize=C_dD_data, doc="Hansen-D 1st-order contributions")

# --- (Loris) G1_Hukkerikar Contributions (C_i) for dP (from Table S5) ---
C_dP_data = {
    'CH3': 2.3037, 'CH2': -0.1664, 'CH': -3.3851, 'C': -5.1979,
    'CH2NH2': 5.2101, 'CHNH2': 0.5616, 'CNH2': 0.3445,
    'CH3NH': 2.5065, 'CH2NH': -0.7159, 'CHNH': -4.6694,
    'CH3N': -0.4250, 'CH2N': -0.7354,
    'OH': 5.2379,
    'CH2O': 0.5137, 'CH3CO': 6.3823, 'CH2CO': 1.2706,
    'CH2=CH': 3.6752, 'CH=CH': 3.0492, 'CH=N': 0.0 # Not in table S5
}
model.C_dP = pe.Param(model.G1_Hukkerikar, initialize=C_dP_data, doc="Hansen-P 1st-order contributions")

# --- (Loris) G1_Hukkerikar Contributions (C_i) for dH (from Table S5) ---
C_dH_data = {
    'CH3': 2.2105, 'CH2': -0.2150, 'CH': -2.6826, 'C': -6.4821,
    'CH2NH2': 6.7984, 'CHNH2': 2.8953, 'CNH2': -0.5293,
    'CH3NH': 7.2551, 'CH2NH': 1.4183, 'CHNH': -2.2824,
    'CH3N': 2.4585, 'CH2N': -7.3014,
    'OH': 11.8005,
    'CH2O': 0.8246, 'CH3CO': 3.4394, 'CH2CO': -0.0788,
    'CH2=CH': 2.7673, 'CH=CH': 0.8631, 'CH=N': 0.0 # Not in table S5
}
model.C_dH = pe.Param(model.G1_Hukkerikar, initialize=C_dH_data, doc="Hansen-H 1st-order contributions")


# --- (Loris) G2_Hukkerikar Contributions (D_j) for Tb (from Table S6) ---
D_Tb_data = {
    '(CH3)2CH': 0.0563,  # Row #1
    '(CH3)3C': 0.0460,   # Row #2
    'CHOH': -0.1849,   # Row #16
    'COH': -0.3371,    # Row #17
}
model.D_Tb = pe.Param(model.G2_Hukkerikar, initialize=D_Tb_data, doc="Boiling T 2nd-order contributions")

# --- (Loris) G2_Hukkerikar Contributions (D_j) for Tm (from Table S6) ---
D_Tm_data = {
    '(CH3)2CH': 0.1542,  # Row #1
    '(CH3)3C': -0.1090,  # Row #2
    'CHOH': 0.1715,   # Row #16
    'COH': 0.6164,    # Row #17
}
model.D_Tm = pe.Param(model.G2_Hukkerikar, initialize=D_Tm_data, doc="Melting T 2nd-order contributions")

# --- (Loris) G2_Hukkerikar Contributions (D_j) for Vm (from Table S6) ---
D_Vm_data = {
    '(CH3)2CH': 0.0021,  # Row #1
    '(CH3)3C': 0.0045,   # Row #2
    'CHOH': 0.0001,   # Row #16
    'COH': 0.0038,    # Row #17
}
model.D_Vm = pe.Param(model.G2_Hukkerikar, initialize=D_Vm_data, doc="Volume 2nd-order contributions")

# --- (Loris) G2_Hukkerikar Contributions (D_j) for dD (from Table S6) ---
D_dD_data = {
    '(CH3)2CH': -0.2581, # Row #1
    '(CH3)3C': 0.1777,   # Row #2
    'CHOH': 0.0394,   # Row #16
    'COH': -0.0212,    # Row #17
}
model.D_dD = pe.Param(model.G2_Hukkerikar, initialize=D_dD_data, doc="Hansen-D 2nd-order contributions")

# --- (Loris) G2_Hukkerikar Contributions (D_j) for dP (from Table S6) ---
D_dP_data = {
    '(CH3)2CH': 0.2698, # Row #1
    '(CH3)3C': -0.0817,   # Row #2
    'CHOH': 0.1176,   # Row #16
    'COH': -1.0511,    # Row #17
}
model.D_dP = pe.Param(model.G2_Hukkerikar, initialize=D_dP_data, doc="Hansen-P 2nd-order contributions")

# --- (Loris) G2_Hukkerikar Contributions (D_j) for dH (from Table S6) ---
D_dH_data = {
    '(CH3)2CH': -0.1884, # Row #1
    '(CH3)3C': 2.3617,   # Row #2
    'CHOH': 0.6216,   # Row #16
    'COH': -0.8848,    # Row #17
}
model.D_dH = pe.Param(model.G2_Hukkerikar, initialize=D_dH_data, doc="Hansen-H 2nd-order contributions")


print("\n--- All Parameters Defined ---")
print("Heat Capacity Cp (313K) Parameters:")
model.cp313.pprint()
print("\n(Loris) Tm (G1) Parameters:")
model.C_Tm.pprint()
print("\n(Loris) Tm (G2) Parameters:")
model.D_Tm.pprint()


# Defining Variables
model.N = pe.Var(model.G1_Hukkerikar, domain=pe.NonNegativeIntegers, bounds=(0, 10), doc="Count of 1st-order groups")
model.M = pe.Var(model.G2_Hukkerikar, domain=pe.NonNegativeIntegers, bounds=(0, 10), doc="Count of 2nd-order groups")


#--Angelina heat capcity vairables---
model.cp313_rule = pe.Var(bounds = (0, None)) #assumes real domain 
model.cp393_rule = pe.Var(bounds = (0, None))


print("\n--- Decision Variables Defined ---")
model.N.pprint()

# Defining Constraints 

#---Angelina Cp calculations---
#will neeed to rename/intialize something for this, but basically have something with the number of each group
#^needs to be the same key (group name like "CH3" as what is defined in the cp313_data)
model.n_groups = pe.Param(model.c, initialize={g: group_count.get(g, 0) for g in model.c})


def cp313_rule(model):
    return model.cp_313_total == sum(model.n_groups[g] * model.cp313[g] for g in model.c)
model.cp313_constraint = pe.Constraint(rule=cp313_rule)

def cp393_rule(model):
    return model.cp_393_total == sum(model.n_groups[g] * model.cp393[g] for g in model.c)
model.cp393_constraint = pe.Constraint(rule=cp393_rule)
