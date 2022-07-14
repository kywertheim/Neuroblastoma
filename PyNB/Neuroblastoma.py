"""
@author: Kenneth Y. Wertheim

Version 13.9 is the PRIMAGE-facing version of the neuroblastoma model.
It supports precise action of chemotherapeutics on intracellular protein species.

This file defines the initial/fixed attributes of each neuroblastic agent.
"""

"""
Load the necessary libraries.
"""
import numpy as np
import math as math

"""
This class combines the parameters and initial conditions defined in Environment.py and 
converts the combined into various initial/fixed agent attributes.
"""
class Neuroblastoma:
    """
    Fixed agent attributes.
    1. Mutation profile.
    2. Functional activities of various intracellular species.

    Initial agent attributes.
    1. Spatial coordinates.
    2. Mechanically related attributes.
    3. Stressors affecting the agent.
    4. Attributes about the agent's cycling progress.
    5. Attributes about the agent's progress towards apoptosis and necrosis.
    6. Attributes indicating the status of each intracellular species.
    """
    def __init__(self, R_tumour=0, histology=0, gradiff=0, x=-1, y=-1, z=-1, MYCN_amp=-1, TERT_rarngm=-1, ATRX_inact=-1, ALT=-1, ALK=-1, MYCN_fn00=-1, MYCN_fn10=-1, MYCN_fn01=-1, MYCN_fn11=-1, MAPK_RAS_fn00=-1, MAPK_RAS_fn10=-1, MAPK_RAS_fn01=-1, MAPK_RAS_fn11=-1, p53_fn=-1, CHK1_fn=-1, p21_fn=-1, p27_fn=-1, p73_fn=-1, CDC25C_fn=-1, CDS1_fn=-1, ID2_fn=-1, IAP2_fn=-1, HIF_fn=-1, BNIP3_fn=-1, JAB1_fn=-1, Bcl2_Bclxl_fn=-1, BAK_BAX_fn=-1, CAS_fn=-1, VEGF_fn=-1, cycle=-1, apop=-1, apop_signal=-1, necro=-1, necro_signal=-1, telo_count=-1):       
        """
        The agent's mutation profile.
        Each mutation is represented by a Boolean variable except ALK, which can take three discrete values.
        """        
        if MYCN_amp == -1:
            self.MYCN_amp = (np.random.uniform(0, 1) < 0.5)
        else:
            self.MYCN_amp = MYCN_amp
        if TERT_rarngm == -1:
            self.TERT_rarngm = (np.random.uniform(0, 1) < 0.5)
        else:
            self.TERT_rarngm = TERT_rarngm
        if ATRX_inact == -1:
            self.ATRX_inact = (np.random.uniform(0, 1) < 0.5)
        else:
            self.ATRX_inact = ATRX_inact
        if ALT == -1:
            self.ALT = (np.random.uniform(0, 1) < 0.5)
        else:
            self.ALT = ALT
        if ALK == -1:
            self.ALK = np.random.randint(0,3) #For ALK, 0 means wild type, 1 means ALK amplification or activation, and 2 means other RAS mutations.
        else:
            self.ALK = ALK
        
        """
        Functional activities of various intracellular species.
        All variables are continuous (0 to 1).
        Note that MYCN_fn and MAPK_RAS_fn depend on the mutation profile.
        """
        if MYCN_fn00 == -1:
            self.MYCN_fn00 = np.random.uniform(0, 1)
        else:
            self.MYCN_fn00 = MYCN_fn00        
        if MYCN_fn10 == -1:
            self.MYCN_fn10 = np.random.uniform(0, 1)
        else:
            self.MYCN_fn10 = MYCN_fn10
        if MYCN_fn01 == -1:
            self.MYCN_fn01 = np.random.uniform(0, 1)
        else:
            self.MYCN_fn01 = MYCN_fn01
        if MYCN_fn11 == -1:
            self.MYCN_fn11 = np.random.uniform(0, 1)
        else:
            self.MYCN_fn11 = MYCN_fn11
        if MAPK_RAS_fn00 == -1:
            self.MAPK_RAS_fn00 = np.random.uniform(0, 1)
        else:
            self.MAPK_RAS_fn00 = MAPK_RAS_fn00        
        if MAPK_RAS_fn10 == -1:
            self.MAPK_RAS_fn10 = np.random.uniform(0, 1)
        else:
            self.MAPK_RAS_fn10 = MAPK_RAS_fn10
        if MAPK_RAS_fn01 == -1:
            self.MAPK_RAS_fn01 = np.random.uniform(0, 1)
        else:
            self.MAPK_RAS_fn01 = MAPK_RAS_fn01
        if MAPK_RAS_fn11 == -1:
            self.MAPK_RAS_fn11 = np.random.uniform(0, 1)
        else:
            self.MAPK_RAS_fn11 = MAPK_RAS_fn11
        if self.MYCN_amp == 0 and self.ALK ==0:
            self.MYCN_fn = self.MYCN_fn00
            self.MAPK_RAS_fn = self.MAPK_RAS_fn00
        elif self.MYCN_amp == 1 and self.ALK == 0:
            self.MYCN_fn = self.MYCN_fn10
            self.MAPK_RAS_fn = self.MAPK_RAS_fn10
        elif self.MYCN_amp == 0 and self.ALK == 1:
            self.MYCN_fn = self.MYCN_fn01
            self.MAPK_RAS_fn = self.MAPK_RAS_fn01
        elif self.MYCN_amp == 1 and self.ALK == 1:
            self.MYCN_fn = self.MYCN_fn11
            self.MAPK_RAS_fn = self.MAPK_RAS_fn11
        elif self.MYCN_amp == 0 and self.ALK == 2:
            self.MYCN_fn = self.MYCN_fn00
            self.MAPK_RAS_fn = self.MAPK_RAS_fn01        
        else:
            self.MYCN_fn = self.MYCN_fn10
            self.MAPK_RAS_fn = self.MAPK_RAS_fn11
        if p53_fn == -1:
            self.p53_fn = np.random.uniform(0, 1)
        else:
            self.p53_fn = p53_fn
        if p73_fn == -1:
            self.p73_fn = np.random.uniform(0, 1)
        else:
            self.p73_fn = p73_fn
        if HIF_fn == -1:
            self.HIF_fn = np.random.uniform(0, 1)
        else:
            self.HIF_fn = HIF_fn
        if CHK1_fn == -1:
            self.CHK1_fn = np.random.uniform(0, 1)
        else:
            self.CHK1_fn = CHK1_fn
        if p21_fn == -1:
            self.p21_fn = np.random.uniform(0, 1)
        else:
            self.p21_fn = p21_fn
        if p27_fn == -1:
            self.p27_fn = np.random.uniform(0, 1)
        else:
            self.p27_fn = p27_fn
        if CDC25C_fn == -1:
            self.CDC25C_fn = np.random.uniform(0, 1)
        else:
            self.CDC25C_fn = CDC25C_fn
        if CDS1_fn == -1:
            self.CDS1_fn = np.random.uniform(0, 1)
        else:
            self.CDS1_fn = CDS1_fn
        if ID2_fn == -1:
            self.ID2_fn = np.random.uniform(0, 1)
        else:
            self.ID2_fn = ID2_fn
        if IAP2_fn == -1:
            self.IAP2_fn = np.random.uniform(0, 1)
        else:
            self.IAP2_fn = IAP2_fn
        if BNIP3_fn == -1:
            self.BNIP3_fn = np.random.uniform(0, 1)
        else:
            self.BNIP3_fn = BNIP3_fn
        if JAB1_fn == -1:
            self.JAB1_fn = np.random.uniform(0, 1)
        else:
            self.JAB1_fn = JAB1_fn
        if Bcl2_Bclxl_fn == -1:
            self.Bcl2_Bclxl_fn = np.random.uniform(0, 1)
        else:
            self.Bcl2_Bclxl_fn = Bcl2_Bclxl_fn
        if BAK_BAX_fn == -1:
            self.BAK_BAX_fn = np.random.uniform(0, 1)
        else:
            self.BAK_BAX_fn = BAK_BAX_fn
        if CAS_fn == -1:
            self.CAS_fn = np.random.uniform(0, 1)
        else:
            self.CAS_fn = CAS_fn
        if VEGF_fn == -1:
            self.VEGF_fn = np.random.uniform(0, 1)
        else:
            self.VEGF_fn = VEGF_fn

        """
        The agent's spatial coordinates, which are continuous variables.
        By default, it's randomly assigned as the imaging biomarkers are not resolved at this scale.
        If the initial domain is a sphere rather than a cube, this algorithm ensures that the agents are uniformly distributed throughout 
        the tumour. https://karthikkaranth.me/blog/generating-random-points-in-a-sphere/#using-normally-distributed-random-numbers
        """
        if x == -1:
            self.x = np.random.uniform(-R_tumour, R_tumour)
        else:
            self.x = x
        if y == -1:
            self.y = np.random.uniform(-R_tumour, R_tumour)
        else:
            self.y = y
        if z == -1:
            self.z = np.random.uniform(-R_tumour, R_tumour)
        else:
            self.z = z

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
        elif telo_count == 1:
            self.telo_count = np.random.randint(47, 61)        
        elif telo_count == 2:
            self.telo_count = np.random.randint(21, 34)
        elif telo_count == 3:
            self.telo_count = np.random.randint(34, 47)
        else:
            self.telo_count = telo_count        
        self.DNA_damage = 0
        self.DNA_unreplicated = 0

        """
        Attributes about the agent's cycling progress.
        1. cycle: This flag (continuous, 0 to 4) indicates the agent's position in the cell cycle.
        2. degdiff: Degree of differentiation (continuous between 0 and 1).
        3. cycdiff: Probability that the agent enters the cell cycle (G0 to G1).
        """
        if cycle == -1:
            self.cycle = np.random.uniform(0, 4)
        else:
            self.cycle = cycle
        if histology == 0:
            if gradiff == 0:
                self.degdiff = 0
            elif gradiff == 1:
                self.degdiff = np.random.uniform(0, 0.2)
            elif gradiff == 2:
                self.degdiff = np.random.uniform(0.2, 0.4)
        elif histology == 2:
            dummy = np.random.uniform(0, 1)
            if dummy < 0.33:
                if gradiff == 0:
                    self.degdiff = 0
                elif gradiff == 1:
                    self.degdiff = np.random.uniform(0, 0.2)
                elif gradiff == 2:
                    self.degdiff = np.random.uniform(0.2, 0.4)
            elif dummy < 0.66:
                self.degdiff = np.random.uniform(0.4, 0.6)
            else:
                self.degdiff = np.random.uniform(0.6, 0.8)
        elif histology == 3:
                self.degdiff = np.random.uniform(0.4, 0.6) 
        elif histology == 5:
            self.degdiff = np.random.uniform(0.6, 0.8)
        elif histology == 6:
            self.degdiff = np.random.uniform(0.8, 1)
        self.cycdiff = 1 - self.degdiff #Assumed to be time-independent, i.e. at an equilibrium.

        """
        Attributes about the agent's progress towards apoptosis and necrosis.
        1. apop: This Boolean variable indicates if the agent is apoptotic.
        2. apop_signal: The total number of apoptotic signals in the agent.
        3. necro: This Boolean variable indicates if the agent is necrotic.
        4. necro_signal: The total number of necrotic signals in the agent.
        5. necro_critical: The number of necrotic signals required to trigger necrosis in the agent.
        """
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

        """
        Attributes indicating the status of each intracellular species.
        These Boolean variables are informed by the mutation profile and functional activities of the species.
        """            
        if self.MYCN_amp == 1 or self.TERT_rarngm == 1:
            self.telo = 1
        else:
            self.telo = 0
        if self.ALT == 0 and self.ATRX_inact == 1:
            self.ALT = 1            
        self.MYCN = 0
        self.MAPK_RAS = 0
        self.JAB1 = 0
        self.CHK1 = 0
        self.CDS1 = 0
        self.CDC25C = 0
        self.ID2 = 0
        self.IAP2 = 0
        self.HIF = 0
        self.BNIP3 = 0
        self.VEGF = 0
        self.p53 = 0
        self.p73 = 0
        self.p21 = 0
        self.p27 = 0
        self.Bcl2_Bclxl = 0
        self.BAK_BAX = 0
        self.CAS = 0        
        if (np.random.uniform(0, 1) < self.MYCN_fn):
            self.MYCN = 1
        else:
            self.MYCN = 0
        if (np.random.uniform(0, 1) < self.MAPK_RAS_fn):
            self.MAPK_RAS = 1
        else:
            self.MAPK_RAS = 0
        if (np.random.uniform(0, 1) < self.JAB1_fn):
            self.JAB1 = 1
        else:
            self.JAB1 = 0
        if (np.random.uniform(0, 1) < self.CHK1_fn and self.DNA_damage == 1):
            self.CHK1 = 1
        elif (np.random.uniform(0, 1) < self.CHK1_fn and self.DNA_damage == 1 and self.MYCN == 1):
            self.CHK1 = 1
        else:
            self.CHK1 = 0
        if (np.random.uniform(0, 1) < self.CDS1_fn and self.DNA_unreplicated == 1):
            self.CDS1 = 1
        else:
            self.CDS1 = 0
        if (np.random.uniform(0, 1) < self.CDC25C_fn):
            self.CDC25C = 1
        else:
            self.CDC25C = 0
        if self.CDC25C == 1:
            if self.CDS1 == 1:
                self.CDC25C = 0
            elif self.CHK1 == 1:
                self.CDC25C = 0
        if (np.random.uniform(0, 1) < self.ID2_fn and self.MYCN == 1):
            self.ID2 = 1
        else:
            self.ID2 = 0
        if (np.random.uniform(0, 1) < self.IAP2_fn and self.hypoxia == 1):
            self.IAP2 = 1
        else:
            self.IAP2 = 0
        if (np.random.uniform(0, 1) < self.HIF_fn and self.hypoxia == 1):
            self.HIF = 1
        elif (np.random.uniform(0, 1) < self.HIF_fn and self.hypoxia == 1 and self.JAB1 == 1):
            self.HIF = 1
        else:
            self.HIF = 0
        if self.HIF == 1: #Note that the effects of p53 and p73 on HIF are considered before p53 and p73 are updated.
            if self.p53 == 1:
                self.HIF = 0
            elif self.p73 == 1:
                self.HIF = 0
        if (np.random.uniform(0, 1) < self.BNIP3_fn and self.HIF == 1):
            self.BNIP3 = 1
        else:
            self.BNIP3 = 0
        if (np.random.uniform(0, 1) < self.VEGF_fn and self.HIF == 1):
            self.VEGF = 1
        else:
            self.VEGF = 0
        if (np.random.uniform(0, 1) < self.p53_fn and self.DNA_damage == 1):
            self.p53 = 1
        elif (np.random.uniform(0, 1) < self.p53_fn and self.DNA_damage == 1 and self.MYCN == 1):
            self.p53 = 1                        
        elif (np.random.uniform(0, 1) < self.p53_fn and self.HIF == 1):
            self.p53 = 1
        elif (np.random.uniform(0, 1) < self.p53_fn and self.HIF == 1 and self.MYCN == 1):
            self.p53 = 1
        else:
            self.p53 = 0
        if self.p53 == 1:
            if (self.MYCN == 1 and self.MYCN_amp == 1): #MYCN's effects on p53 depend on the former's amplification status (Tang et al., 2006).
                self.p53 = 0
        if (np.random.uniform(0, 1) < self.p73_fn and self.CHK1 == 1):
            self.p73 = 1
        elif (np.random.uniform(0, 1) < self.p73_fn and self.HIF == 1):
            self.p73 = 1
        else:
            self.p73 = 0
        if (np.random.uniform(0, 1) < self.p21_fn and self.HIF == 1):
            self.p21 = 1
        elif (np.random.uniform(0, 1) < self.p21_fn and self.p53 == 1):
            self.p21 = 1
        else:
            self.p21 = 0
        if self.p21 == 1:
            if self.MAPK_RAS == 1:
                self.p21 = 0
            elif self.MYCN == 1:
                self.p21 = 0
        if (np.random.uniform(0, 1) < self.p27_fn and self.HIF == 1):
            self.p27 = 1
        elif (np.random.uniform(0, 1) < self.p27_fn and self.p53 == 1):
            self.p27 = 1
        else:
            self.p27 = 0
        if self.p27 == 1:
            if self.MAPK_RAS == 1:
                self.p27 = 0
            elif self.MYCN == 1:
                self.p27 = 0
        if (np.random.uniform(0, 1) < self.Bcl2_Bclxl_fn):
            self.Bcl2_Bclxl = 1
        else:
            self.Bcl2_Bclxl = 0
        if self.Bcl2_Bclxl == 1:
            if self.BNIP3 == 1:
                self.Bcl2_Bclxl = 0
            elif self.p53 == 1:
                self.Bcl2_Bclxl = 0
            elif self.p73 == 1:
                self.Bcl2_Bclxl = 0
        if (np.random.uniform(0, 1) < self.BAK_BAX_fn and self.hypoxia == 1):
            self.BAK_BAX = 1
        elif (np.random.uniform(0, 1) < self.BAK_BAX_fn and self.p53 == 1):
            self.BAK_BAX = 1
        elif (np.random.uniform(0, 1) < self.BAK_BAX_fn and self.p73 == 1):
            self.BAK_BAX = 1
        else:
            self.BAK_BAX = 0
        if self.BAK_BAX == 1:
            if self.Bcl2_Bclxl == 1:
                self.BAK_BAX = 0
            elif self.IAP2 == 1:
                self.BAK_BAX = 0
        if (np.random.uniform(0, 1) < self.CAS_fn and self.BAK_BAX == 1 and self.ATP == 1):
            self.CAS = 1
        elif (np.random.uniform(0, 1) < self.CAS_fn and self.hypoxia == 1 and self.ATP == 1):
            self.CAS = 1
        else:
            self.CAS = 0