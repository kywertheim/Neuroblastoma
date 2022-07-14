# -*- coding: utf-8 -*-
"""
Created on Mon Oct  7 12:44:20 2019

@author: Kenneth Y. Wertheim
"""

"""
PARAMETERS.
1. Extracellular Environment.
    a. EC_EVM_MMP: 
        Degrades the extracellular matrix and promotes metastasis
    b. EC_EVM_VEGF:
        Promotes angiogenesis and metastasi
2. NeuroBlastoma
    a. NB_count (int):
            Number of NeuroBlastoma cells
    b. NB_position_mean (float[3]):
            Mean of NeuroBlastoma positions, for Gaussian distribution
    c. NB_position_variance (float[3]):
            Variance of NeuroBlastoma positions, for Gaussian distribution
    d. NB_dummy_force (float):
            Mean of all NeuroBlastoma force magnitudes
    e. NB_neighbour_histogram (int[]):
            Histogram of NeuroBlastoma neighbourhood volumes
            Note: The length of this should equal max neighbours+1
"""
class ValidationData: 
    def __init__(self
        , EC_EVM_MMP=0
        , EC_EVM_VEGF=0
        , NB_count=0
        , NB_position_mean=[0] * 3
        , NB_position_variance=[0] * 3
        , NB_dummy_force=0
        , NB_neighbour_histogram=[0] * 7
    ):
        self.EC_EVM_MMP = EC_EVM_MMP
        self.EC_EVM_VEGF = EC_EVM_VEGF
        self.NB_count = NB_count
        self.NB_position_mean = NB_position_mean
        self.NB_position_variance = NB_position_variance
        self.NB_dummy_force = NB_dummy_force
        self.NB_neighbour_histogram = NB_neighbour_histogram