"""
@author: Kenneth Y. Wertheim

Version 13.9 is the PRIMAGE-facing version of the neuroblastoma model.
It supports precise action of chemotherapeutics on intracellular protein species.

This file defines the initial attributes of each Schwann cell agent.
"""

"""
Load the necessary libraries.
"""
import numpy as np
import math as math

"""
This class combines the parameters and initial conditions defined in Environment.py and 
converts the combined into various initial agent attributes.
"""
class Schwann:
    """
    Initial agent attributes.
    1. Spatial coordinates.
    2. Mechanically related attributes.
    3. Stressors affecting the agent.
    4. Attributes about the agent's progress towards cycling, apoptosis, and necrosis.
    """    
    def __init__(self, R_tumour=0, cycle=-1, apop=-1, apop_signal=-1, necro=-1, necro_signal=-1, telo_count=-1):                
        """
        The agent's spatial coordinates, which are continuous variables.
        By default, it's randomly assigned as the imaging biomarkers are not resolved at this scale.
        If the initial domain is a sphere rather than a cube, this algorithm ensures that the agents are uniformly distributed throughout 
        the tumour. https://karthikkaranth.me/blog/generating-random-points-in-a-sphere/#using-normally-distributed-random-numbers
        """
        self.x = np.random.uniform(-R_tumour, R_tumour)
        self.y = np.random.uniform(-R_tumour, R_tumour)
        self.z = np.random.uniform(-R_tumour, R_tumour)

        """
        Mechanically related attributes.
        1. Fx, Fy, and Fz: Forces in independent directions (kg s-2 micron).
        2. overlap: The agent's overlap with its neighbouring agents.
        3. neighbours: The number of agents within the agent's search distance.
        4. mobile: This Boolean variable indicates whether the agent is mobile.
        """  
        self.Fx = 0
        self.Fy = 0
        self.Fz = 0
        self.overlap = 0
        self.neighbours = 0        
        self.mobile = 1

        """
        Stressors affecting the agent.
        1. hypoxia: This Boolean variable indicates whether the agent is hypoxic.
        2. nutrient: This Boolean variable indicates whether the agent has sufficient nutrients (other than O2).
        3. ATP: This Boolean variable indicates whether the agent has sufficient energy.
        4. telo_count: The total number of telomere units in the agent, a discrete variable ranging from zero to 60.
        5. DNA_damage: This Boolean variable indicates whether the agent has damaged DNA.
        6. DNA_unreplicated: This Boolean variable indicates whether the agent has unreplicated DNA.
        """
        self.hypoxia = 0
        self.nutrient = 1
        self.ATP = 1
        if telo_count == -1:
            self.telo_count = np.random.randint(35, 46)
        else:
            self.telo_count = telo_count        
        self.DNA_damage = 0
        self.DNA_unreplicated = 0        

        """
        Attributes about the agent's progress towards cycling, apoptosis, and necrosis.
        1. cycle: This flag (continuous, 0 to 4) indicates the agent's position in the cell cycle.        
        1. apop: This Boolean variable indicates if the agent is apoptotic.
        2. apop_signal: The total number of apoptotic signals in the agent.
        3. necro: This Boolean variable indicates if the agent is necrotic.
        4. necro_signal: The total number of necrotic signals in the agent.
        5. necro_critical: The number of necrotic signals required to trigger necrosis in the agent.
        """
        if cycle == -1:
            self.cycle = np.random.uniform(0, 4)
        else:
            self.cycle = cycle                
        if apop == -1:
            self.apop = 0
        else:
            self.apop = apop
        if apop_signal == -1:
            self.apop_signal = 0
        else:
            self.apop_signal = apop_signal
        if necro == -1:
            self.necro = 0
        else:
            self.necro = necro
        if necro_signal == -1:
            self.necro_signal = 0
        else:
            self.necro_signal = necro_signal
        self.necro_critical = np.random.randint(3, 169) #It is between 3 and 168, inclusive (Warren et al., 2016).