"""
@author: Kenneth Y. Wertheim

Version 13.9 is the PRIMAGE-facing version of the neuroblastoma model.
It supports precise action of chemotherapeutics on intracellular protein species.

This file contains all the parameters and initial conditions of the model, 
either defined directly or derived.
"""

"""
Load the necessary libraries.
"""
import math as math
import numpy as np

"""
Parameters and initial conditions.
"""
class Environment:
    """
    If the JSON load leads to a version warning, this value will be set.
    This is a static variable, so it will affect all instances.
    """
    VERSION_WARNING = False
    
    """
    This function returns the model's semantic version [official release, major, minor].
    This manually updated constant allows JSON exports to show the model version.
    """
    def VERSION():
        return [0, 13, 9, Environment.VERSION_WARNING]

    """
    The following constants are defined directly.
    1. Integration with data from consortium partners.
    2. Physical parameters (neuroblasts and Schwann cells).
    3. Cell cycle parameters (neuroblasts and Schwann cells).
    4. Stress-Related parameters (neuroblasts and Schwann cells).
    5. Cell death parameters (neuroblasts and Schwann cells).
    6. NB-SC crosstalk parameters.
    7. Mechanical model parameters.
    8. Microenvironment parameters.
    9. Initial conditions (neuroblasts).
    10. Initial conditions (Schwann cells).
    """    
    def __init__(self):
        """   
        Data layer 0 (integration with imaging biomarkers), part 1.
        1. V_tumour: Initial volume of tumour (cubic microns).
        2. boundary_max: Factor by which the tumour can expand in all three directions.                
        3. x_displace: Displacement magnitude in the boundary (x-direction).
        4. y_displace: Displacement magnitude in the boundary (y-direction).
        5. z_displace: Displacement magnitude in the boundary (z-direction).
        """
        self.V_tumour = 2e6 #The volume of a finite element is around 2000**3 cubic microns, but the Python model cannot handle more than a few hundred agents.
        self.boundary_max = 1.26 #According to Aherne and Buck (1971), the volume doubling time is around 60 days.
        self.x_displace = 2
        self.y_displace = 2
        self.z_displace = 2

        """   
        Data layer 0 (integration with imaging biomarkers), part 2.
        1. cellularity: Initial cellularity in the tumour (continuous, 0.05 to 0.95).
        2. histology_init and histology: Histological type.
        3. gradiff: Grade of differentiation for neuroblastoma.
        4. theta_sc: Fraction of Schwann cells in the cell population (continuous, 0.05 to 0.95), 
        derived from the histological type and grade of differentiation.
        """
        self.cellularity = np.random.uniform(0.05, 0.95)        
        self.histology_init = 0 #0 is neuroblastoma, 1 is ganglioneuroblastoma, 2 is nodular ganglioneuroblastoma, 3 is intermixed ganglioneuroblastoma, 4 is ganglioneuroma, 5 is maturing ganglioneuroma, and 6 is mature ganglioneuroma.
        if self.histology_init == 1: #If it is 1 or 4, assign the subtype stochastically.
            if np.random.uniform(0, 1) < 0.5:
                self.histology = 2
            else:
                self.histology = 3
        elif self.histology_init == 4: #If it is 1 or 4, assign the subtype stochastically.
            if np.random.uniform(0, 1) < 0.5:
                self.histology = 5
            else:
                self.histology = 6
        else:
            self.histology = self.histology_init
        self.gradiff = np.random.randint(0, 3) #0 is undifferentiated, 1 is pooly differentiated, and 2 is differentiating.
        if self.histology == 0:
            if self.gradiff == 0:
                self.theta_sc = np.random.uniform(0.05, 0.17)
            elif self.gradiff == 1:
                self.theta_sc = np.random.uniform(0.17, 0.33)
            elif self.gradiff == 2:
                self.theta_sc = np.random.uniform(0.33, 0.5)              
        elif self.histology == 2:
            dummy = np.random.uniform(0, 1)
            if dummy < 0.33:
                if self.gradiff == 0:
                    self.theta_sc = np.random.uniform(0.05, 0.17)
                elif self.gradiff == 1:
                    self.theta_sc = np.random.uniform(0.17, 0.33)
                elif self.gradiff == 2:             
                    self.theta_sc = np.random.uniform(0.33, 0.5)
            elif dummy < 0.66:
                self.theta_sc = np.random.uniform(0.5, 0.67)
            else:
                self.theta_sc = np.random.uniform(0.67, 0.83)                
        elif self.histology == 3:
                self.theta_sc = np.random.uniform(0.5, 0.67)
        elif self.histology == 5:
            self.theta_sc = np.random.uniform(0.67, 0.83)
        elif self.histology == 6:
            self.theta_sc = np.random.uniform(0.83, 0.95)

        """   
        Data layer 0 (integration with imaging biomarkers), part 3.
        1. O2: Initial oxygen level (continuous, 0 to 1), scaled by the level in the kidney, 72 mmHg (Carreau et al., 2011).                
        2. chemo_start: Time points at which chemo cycles begin.
        3. chemo_end: Time points at which chemo cycles end.
        4. chemo_effects: Probabilities that CHK1, JAB1, HIF, MYCN, TEP1, and p53 are inhibited by chemotherapy.
        """
        self.O2 = np.random.uniform(2/72, 32/72) #Oxygen level in hypoxic tumours = 2 to 32 mmHg (McKeown, 2014).
        self.chemo_start = [0, 240] #First two-week cycle of Rapid COJEC.
        self.chemo_end = [96, 336] #First two-week cycle of Rapid COJEC.
        self.chemo_effects = [0.6026663265, 0.6026663265, 0, 0.6026663265, 0.6026663265, 0] #Rapid COJEC inhibits these species.

        """
        Data Layer 1 (integration with genetic/molecular biomarkers of neuroblasts).
        1. MYCN_amp: MYCN amplification status (categorical, 0 or 1), default (-1) means unknown.
        2. TERT_rarngm: TERT rearrangement status (categorical, 0 or 1), default (-1) means unknown.
        3. ATRX_inact: ATRX inactivation status (categorical, 0 or 1), default (-1) means unknown.
        While this mutation switches ALT on, it may also have unrelated functions, but this model does not consider them.
        4. ALT: Alternative lengthening of telomeres status (categorical, 0 or 1), default (-1) means unknown.                
        5. ALK: ALK amplification or activating mutation status (discrete, 0 or 1 or 2), default (-1) means unknown.
        """
        self.MYCN_amp=np.random.randint(0, 2)
        self.TERT_rarngm=np.random.randint(0, 2)
        self.ATRX_inact=np.random.randint(0, 2)
        self.ALT=np.random.randint(0, 2)
        self.ALK=np.random.randint(0, 3) #0 means wild type, 1 means ALK amplification or activation, and 2 means other RAS mutations.
        
        """
        Data Layer 2 (integration with genetic/molecular biomarkers of neuroblasts). 
        1. MYCN_fnxx: Functional activity of MYCN (continuous, 0 to 1), default (-1) means unknown.
        This value depends on MYCN_amp and ALK.
        2. MAPK_RAS_fnxx: Functional activity of MAPK/RAS signalling (continuous, 0 to 1), default (-1) means unknown.
        This value depends on MYCN_amp and ALK.
        3. p53_fn: Functional activity of p53 signalling (continuous, 0 to 1), default (-1) means unknown.
        4. p73_fn: Functional activity of p73 signalling (continuous, 0 to 1), default (-1) means unknown.
        5. HIF_fn: Functional activity of HIF signalling (continuous, 0 to 1), default (-1) means unknown.
        """
        self.MYCN_fn11=0.942648156 #LHC, index 564.
        self.MYCN_fn00=0.8*0.71*self.MYCN_fn11  
        self.MYCN_fn10=0.71*self.MYCN_fn11
        self.MYCN_fn01=0.8*self.MYCN_fn11
        self.MAPK_RAS_fn11=0.377627766 #LHC, index 564.
        self.MAPK_RAS_fn10=0.77*self.MAPK_RAS_fn11
        self.MAPK_RAS_fn01=0.003161348 #LHC, index 564.
        self.MAPK_RAS_fn00=0.77*self.MAPK_RAS_fn01
        self.p53_fn=0.198089528 #LHC, index 564.
        self.p73_fn=0.141041534 #LHC, index 564.
        self.HIF_fn=0.591769646 #LHC, index 564.
        
        """
        Data Layer 3 (integration with genetic/molecular biomarkers of neuroblasts).
        Activity levels of various species/pathways (continuous, 0 to 1), default (-1) means unknown.
        Assumed to be one and not selected for calibration.
        """        
        self.CHK1_fn=1
        self.p21_fn=1
        self.p27_fn=1
        self.CDC25C_fn=1
        self.CDS1_fn=1
        self.ID2_fn=1
        self.IAP2_fn=1
        self.BNIP3_fn=1
        self.JAB1_fn=1
        self.Bcl2_Bclxl_fn=1
        self.BAK_BAX_fn=1
        self.CAS_fn=1
        self.VEGF_fn=1

        """   
        Physical parameters (neuroblasts and Schwann cells).
        1. rho_tumour: Initial cell density of the cellular region in the tumour (cells per cubic micron).
        2. R_cell: Cell radius in microns at the beginning of the cell cycle.
        3. R_voxel: Half of a voxel's side length in microns.
        """
        self.rho_tumour = 9.39e-05
        self.R_cell = 11/2 #The radii of most animal cells range from 5 to 15 microns (Del Monte, 2009).
        self.R_voxel = 15 #This must be at least as big as R_cell.
        
        """
        Cell cycle parameters (neuroblasts and Schwann cells).
        1. cycle_stages: Durations of G1, S, G2, and M in hours (Harper and Brooks, 2005).
        2. glycoEff: Efficiency of glycolysis compared to oxidative phosphorylation.
        3. P_cycle_nb: Basal probability of cycling for neuroblasts.
        4. P_cycle_sc: Basal probability of cycling for Schwann cells.
        """
        self.cycle_stages = [12, 6, 4, 2]
        self.glycoEff = 1/15 #du Plessis et al. (2015).
        self.P_cycle_nb = 0.0457 #LHC_Cal6, index 9.
        self.P_cycle_sc = 0.0325 #LHC_Cal6, index 9.

        """
        Stress-Related parameters (neuroblasts and Schwann cells).
        1. C50_necro: Concentration (M) or partial pressure of oxygen (mmHg) at which 50 % of the tumour cell population die through necrosis.
        However, this is used to calculate the probability that a living cell is hypoxic. This is time-independent, i.e at equilibrium.
        2. telo_maximum: Maximum number of telomere units in a cell.
        3. telo_critical: Maximum number of telomere units in a senescent cell.
        4. P_telorp: Probability of gaining one unit of telomere in an hour, when telomerase or ALT is active.
        5. P_apopChemo: Probability of gaining DNA damage in an hour due to chemotherapy.
        6. P_DNA_damageHypo: Probability of gaining DNA damage in an hour due to hypoxia.
        7. P_DNA_damagerp: Probability of repairing DNA damage in an hour.
        8. P_unrepDNA: Probability of gaining unreplicated DNA in an hour.
        9. P_unrepDNAHypo: Probability of gaining unreplicated DNA in an hour due to hypoxia.
        10. P_unrepDNArp: Probability of repairing unreplicated DNA in an hour.
        """
        self.C50_necro = 1.2/2.2779e-4/32 #Given in mmHg, converted to g per dm3, then to moles per dm3 or M (Warren and Partridge, 2016).
        self.telo_maximum = 60 #A normal human foetal cell population will divide between 40 and 60 times before entering a senescence phase due to shortening telomeres (Hayflick et al., 1961).
        self.telo_critical = 20 #A normal human foetal cell population will divide between 40 and 60 times before entering a senescence phase due to shortening telomeres (Hayflick et al., 1961).
        self.P_telorp = 0.08895382 #LHC, index 564.
        self.P_apopChemo = 0.644 #LHC_Cal4, index 754.
        self.P_DNA_damageHypo = 0.772947675 #LHC, index 564.
        self.P_DNA_damagerp = 0.771497002 #LHC, index 564.
        self.P_unrepDNA = 0 #Assumed to be insignificant compared to hypoxia-induced events.
        self.P_unrepDNAHypo = 0.434578817 #LHC, index 564.
        self.P_unrepDNArp = 0.890953082 #LHC, index 564.

        """
        Cell death parameters (neuroblasts and Schwann cells).
        1. P_DNA_damage_pathways: Probability of DNA damages triggering CAS-independent pathways to induce an apoptotic signal in an hour.
        2. apop_critical: Number of apoptotic signals needed to kill the cell.
        3. P_apoprp: Probability of losing an apoptotic signal in an unstressed cell in an hour.
        4. P_2ndnecro: Probability of secondary necrosis in an hour.
        5. P_necroIS: Probability of the immune system triggering a necrotic signal in a living cell per necrotic cell present per hour.
        6. P_necrorp: Probability of losing a necrotic signal in an unstressed cell in an hour.
        """
        self.P_DNA_damage_pathways = 0.256 #LHC_Cal4, index 754.
        self.apop_critical = 3 #Elmore (2007).
        self.P_apoprp = 0.957831979 #LHC, index 564.
        self.P_2ndnecro = 0.2 #Dunster, Byrne, and King (2014).
        self.P_necroIS = 0.57675841 #LHC, index 564.
        self.P_necrorp = 0.98970852 #LHC, index 564.

        """
        NB-SC crosstalk parameters.
        1. scpro_jux: Scaling factor for the influence of neuroblasts on Schwann cell proliferation, juxtacrine.
        2. nbdiff_jux: Scaling factor for the influence of Schwann cells on neuroblast differentiation, juxtacrine.
        3. nbdiff_amount: Amount of neuroblast differentiation achieved in an hour, triggered by Schwann cells.
        4. nbapop_jux: Scaling factor for the influence of Schwann cells on neuroblast apoptosis, juxtacrine.
        5. scpro_para: Scaling factor for the influence of neuroblasts on Schwann cell proliferation, paracrine.             
        6. nbdiff_para: Scaling factor for the influence of Schwann cells on neuroblast differentiation, paracrine.
        7. nbapop_para: Scaling factor for the influence of Schwann cells on neuroblast apoptosis, paracrine.            
        """
        self.scpro_jux = 0.0374*0.1 #LHC_Cal3b, index 163, modified.
        self.nbdiff_jux = 0.000521 #LHC_Cal6, index 9.
        self.nbdiff_amount = 0.01 #Assumed.
        self.nbapop_jux = 0.0304 #LHC_Cal3b, index 163.
        self.scpro_para = self.scpro_jux/10 #Assumed.
        self.nbdiff_para = self.nbdiff_jux/10 #Assumed.
        self.nbapop_para = self.nbapop_jux/10 #Assumed. 
        
        """
        Mechanical model parameters.
        1. min_overlap: Minimum overlap below which two cells cannot interact (microns).
        2. k1: Linear force law parameter in N m-1.
        3. R_neighbours: A cell's search distance in the search for neighbours (microns).
        4. N_neighbours: Number of other cells allowed within a cell's search distance before contact inhibition activates.
        5. k_locom: Factor by which a cell magnifies the force acting on it upon contact inhibition.
        6. mu: Viscosity in N s m-1.
        7. dt: Time step for the mechanical model in seconds.
        """
        self.min_overlap = -4e-6*1e6*0 #Set to zero so that when two cells are just touching, they stop interacting immediately, i.e. no bouncing off.
        self.k1 = 2.2e-3 #Pathmanathan et al. (2009).
        self.R_neighbours = 1.05*(2*1.5*self.R_cell) #Calibrated by trial and error.
        self.N_neighbours = 2 #Assumed.
        self.k_locom = 2 #Assumed.
        self.mu = 0.4 #Pathmanathan et al. (2009).
        self.dt = 36 #Pathmanathan et al. (2009).

        """
        Microenvironment parameters.
        1. P_O20: Production rate of oxygen given in moles per cell per hour.
        2. Cs_O2: Concentration scale or maximum partial pressure of oxygen in moles per dm3 (M).
        3. staticO2: If this is on, the initial oxygen level is taken as the equilibrium and the vasculature will return the O2 level to this equilibrium at the end of every time step.
        4. ang_critical: Number of angiogenic signals needed to update the vasculature.
        5. P_matrix: Production rate of matrix by one living Schwann cell (cubic microns per hour).
        6. P_lysis: Probability of an apoptotic or necrotic cell being engulfed by an immune cell in an hour.
        """
        self.P_O20 = -250e-12*60/80000 #Negative because it is actually consumption (Grimes et al., 2014).
        self.Cs_O2 = 72/2.2779e-4/32 #Oxygen level in the kidney (Carreau et al., 2011).
        self.staticO2 = 0
        self.ang_critical = 100 #In an experiment (Utzinger et al., 2015), it took 100 hours for two microvessel fragments to inosculate to the vascular network.
        self.P_matrix = 3.125e-12*0.12*1.89e12 #Protein production rate (Conlon et al., 2003), proportion of collagen in the output (DeClerck et al., 1987), and the volume of hydrated collagen I (Levick, 1987).
        self.P_lysis = 0.35 #Jagiella et al. (2016).
        
        """   
        Initial conditions (neuroblasts).
        1. cycle: A flag (continuous, 0 to 4) indicating the cell's position in the cell cycle, default (-1) means random initialisation.
        2. apop: A flag (Boolean variable) indicating if the cell is apoptotic, default (-1) means it is not apoptotic (zero).
        3. apop_signal: Number of apoptotic signals (categorical, 0, 1, 2, and 3), default (-1) means setting it to zero.
        4. necro: A flag (Boolean variable) indicating if the cell is necrotic, default (-1) means it is not necrotic (zero).
        5. necro_signal: Number of necrotic signals (categorical, integers from 0 to 168), default (-1) means setting it to zero.
        6. telo_count: Number of telomere units (categorical, integers from 0 to 60).      
        """
        self.cycle=-1
        self.apop=-1
        self.apop_signal=-1
        self.necro=-1
        self.necro_signal=-1 #Maximum is 168 (Warren et al., 2016).
        self.telo_count=-1 #-1 means random initialisation between 35 and 45, inclusive.
        if (self.ALT == 1 or self.ATRX_inact == 1) and (self.MYCN_amp == 1 or self.TERT_rarngm == 1):
            self.telo_count = np.random.randint(1, 3)
        elif self.ALT == 1 or self.ATRX_inact == 1:
            self.telo_count = 1 #1 means random initialisation between 47 and 60, inclusive.
        elif self.MYCN_amp == 1 or self.TERT_rarngm == 1:    
            self.telo_count = 2 #2 means random initialisation between 21 and 33, inclusive.
        else:
            self.telo_count = 3 #3 means random initialisation between 34 and 46, inclusive.
        
        """
        Initial conditions (Schwann cells).
        1. cycle_sc: A flag (continuous, 0 to 4) indicating the cell's position in the cell cycle, default (-1) means random initialisation.
        2. apop_sc: A flag (Boolean variable) indicating if the cell is apoptotic, default (-1) means it is not apoptotic (zero).
        3. apop_signal_sc: Number of apoptotic signals (categorical, 0, 1, 2, and 3), default (-1) means setting it to zero.
        4. necro_sc: A flag (Boolean variable) indicating if the cell is necrotic, default (-1) means it is not necrotic (zero).
        5. necro_signal_sc: Number of necrotic signals (categorical, integers from 0 to 168), default (-1) means setting it to zero.
        6. telo_count_sc: Number of telomere units (categorical, integers from 0 to 60).              
        """
        self.cycle_sc=-1
        self.apop_sc=-1
        self.apop_signal_sc=-1
        self.necro_sc=-1
        self.necro_signal_sc=-1 #Maximum is 168 (Warren et al., 2016).
        self.telo_count_sc=-1 #-1 means random initialisation between 35 and 45, inclusive.

    """
    Half of the initial tumour length in one dimension in microns (derived).
    """       
    def R_tumour(self):
        return 0.5*(self.V_tumour)**(1/3)

    """
    Tumour boundary in microns (derived).
    """       
    def x_bc_minus(self):
        return -self.boundary_max*self.R_tumour()
    
    def x_bc_plus(self):
        return self.boundary_max*self.R_tumour()
    
    def y_bc_minus(self):
        return -self.boundary_max*self.R_tumour()
    
    def y_bc_plus(self): 
        return self.boundary_max*self.R_tumour() 
    
    def z_bc_minus(self):
        return -self.boundary_max*self.R_tumour()
    
    def z_bc_plus(self):
        return self.boundary_max*self.R_tumour()

    """
    Initial number of neuroblasts (derived).
    """               
    def N_cell (self):
        return np.ceil(self.rho_tumour*self.V_tumour*self.cellularity*(1-self.theta_sc))

    """
    Initial number of Schwann cells (derived).
    """    
    def N_scell (self):
        return np.ceil(self.rho_tumour*self.V_tumour*self.cellularity*self.theta_sc)

    """   
    Voxel volume in cubic microns (derived).
    """
    def V_grid(self):
        return (2*self.R_voxel)**3

    """   
    Voxel side area in microns squared (derived).
    """    
    def A_grid(self):
        return (2*self.R_voxel)**2

    """   
    Initial conditions (extracellular environment, matrix).
    Matrix volume fraction in each voxel (continuous, 0 to 1).
    """
    def matrix_grid(self):
        N_units = math.ceil(self.R_tumour()/self.R_voxel) #Number of voxels in each dimension.
        if (N_units % 2) == 0:
            N_units += 1 #Number of voxels in each dimension.
        matrix_grid = (1-self.cellularity)*np.ones((N_units, N_units, N_units))
        return matrix_grid