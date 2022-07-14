"""
@author: Kenneth Y. Wertheim

Version 13.9 is the PRIMAGE-facing version of the neuroblastoma model.
It supports precise action of chemotherapeutics on intracellular protein species.

This file contains the functions used to update the model's state vector
and orders the updates to simulate cancer progression.
"""

"""
Load the necessary libraries.
"""
import itertools as itt
import numpy as np
import math as math
import statistics as stat
import random as rn
import copy as cp
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation
import time
from PyNB.Neuroblastoma import Neuroblastoma
from PyNB.Schwann import Schwann
from PyNB.Environment import Environment
from PyNB.JSON import load
from PyNB.JSON import save
from PyNB.ValidationData import ValidationData
import argparse

"""
Define a parser for configuration.
"""
def getParser():
    """
    The parser allows the user to configure a simulation directly.
    """
    parser = argparse.ArgumentParser(description='Configure initial state export.')
    parser.add_argument('-rng', dest='seed', type=int, nargs=1, default=[0],
                        help='Random seed, r', required=False);
    """
    Override the default initialisation routine for neuroblast agents.
    """
    parser.add_argument('-x', dest='x', type=float, nargs=1, default=[-1],
                       help='Coordinate of the cell, x', required=False);
    parser.add_argument('-y', dest='y', type=float, nargs=1, default=[-1],
                       help='Coordinate of the cell, y', required=False);
    parser.add_argument('-z', dest='z', type=float, nargs=1, default=[-1],
                       help='Coordinate of the cell, z', required=False);
    parser.add_argument('-MYCN_amp', dest='MYCN_amp', type=int, nargs=1, default=[envn.MYCN_amp],
                       help='MYCN amplification status, MYCN_amp', required=False);                        
    parser.add_argument('-TERT_rarngm', dest='TERT_rarngm', type=int, nargs=1, default=[envn.TERT_rarngm],
                       help='TERT rearrangement status, TERT_rarngm', required=False);
    parser.add_argument('-ATRX_inact', dest='ATRX_inact', type=int, nargs=1, default=[envn.ATRX_inact],
                       help='ATRX inactivation status, ATRX_inact', required=False);
    parser.add_argument('-ALT', dest='ALT', type=int, nargs=1, default=[envn.ALT],
                       help='Alternative lengthening of telomeres status, ALT', required=False);
    parser.add_argument('-ALK', dest='ALK', type=int, nargs=1, default=[envn.ALK],
                       help='ALK amplification or activating mutation status, ALK', required=False);
    parser.add_argument('-MYCN_fn00', dest='MYCN_fn00', type=float, nargs=1, default=[envn.MYCN_fn00],
                       help='Function of MYCN (MYCN_amp off and ALK off), MYCN_fn00', required=False);
    parser.add_argument('-MYCN_fn10', dest='MYCN_fn10', type=float, nargs=1, default=[envn.MYCN_fn10],
                   help='Function of MYCN (MYCN_amp on and ALK off), MYCN_fn10', required=False);
    parser.add_argument('-MYCN_fn01', dest='MYCN_fn01', type=float, nargs=1, default=[envn.MYCN_fn01],
                   help='Function of MYCN (MYCN_amp off and ALK on), MYCN_fn01', required=False);
    parser.add_argument('-MYCN_fn11', dest='MYCN_fn11', type=float, nargs=1, default=[envn.MYCN_fn11],
                   help='Function of MYCN (MYCN_amp on and ALK on), MYCN_fn11', required=False);
    parser.add_argument('-MAPK_RAS_fn00', dest='MAPK_RAS_fn00', type=float, nargs=1, default=[envn.MAPK_RAS_fn00],
                       help='Function of MAPK/RAS signalling (MYCN_amp off and ALK off), MAPK_RAS_fn00', required=False);
    parser.add_argument('-MAPK_RAS_fn10', dest='MAPK_RAS_fn10', type=float, nargs=1, default=[envn.MAPK_RAS_fn10],
                   help='Function of MAPK/RAS signalling (MYCN_amp on and ALK off), MAPK_RAS_fn10', required=False);
    parser.add_argument('-MAPK_RAS_fn01', dest='MAPK_RAS_fn01', type=float, nargs=1, default=[envn.MAPK_RAS_fn01],
                   help='Function of MAPK/RAS signalling (MYCN_amp off and ALK on), MAPK_RAS_fn01', required=False);
    parser.add_argument('-MAPK_RAS_fn11', dest='MAPK_RAS_fn11', type=float, nargs=1, default=[envn.MAPK_RAS_fn11],
                   help='Function of MAPK/RAS signalling (MYCN_amp on and ALK on), MAPK_RAS_fn11', required=False);
    parser.add_argument('-p53_fn', dest='p53_fn', type=float, nargs=1, default=[envn.p53_fn],
                       help='Function of p53 signalling, p53_fn', required=False);
    parser.add_argument('-CHK1_fn', dest='CHK1_fn', type=float, nargs=1, default=[envn.CHK1_fn],
                       help='Function of CHK1 signalling, CHK1_fn', required=False);
    parser.add_argument('-p21_fn', dest='p21_fn', type=float, nargs=1, default=[envn.p21_fn],
                       help='Function of p21 signalling, p21_fn', required=False);
    parser.add_argument('-p27_fn', dest='p27_fn', type=float, nargs=1, default=[envn.p27_fn],
                       help='Function of p27 signalling, p27_fn', required=False);
    parser.add_argument('-p73_fn', dest='p73_fn', type=float, nargs=1, default=[envn.p73_fn],
                       help='Function of p73 signalling, p73_fn', required=False);
    parser.add_argument('-CDC25C_fn', dest='CDC25C_fn', type=float, nargs=1, default=[envn.CDC25C_fn],
                       help='Function of CDC25C signalling, CDC25C_fn', required=False);
    parser.add_argument('-CDS1_fn', dest='CDS1_fn', type=float, nargs=1, default=[envn.CDS1_fn],
                       help='Function of CDS1 signalling, CDS1_fn', required=False);    
    parser.add_argument('-ID2_fn', dest='ID2_fn', type=float, nargs=1, default=[envn.ID2_fn],
                       help='Function of ID2 signalling, ID2_fn', required=False);
    parser.add_argument('-IAP2_fn', dest='IAP2_fn', type=float, nargs=1, default=[envn.IAP2_fn],
                       help='Function of IAP2 signalling, IAP2_fn', required=False);
    parser.add_argument('-HIF_fn', dest='HIF_fn', type=float, nargs=1, default=[envn.HIF_fn],
                       help='Function of HIF signalling, HIF_fn', required=False);
    parser.add_argument('-BNIP3_fn', dest='BNIP3_fn', type=float, nargs=1, default=[envn.BNIP3_fn],
                       help='Function of BNIP3 signalling, BNIP3_fn', required=False);
    parser.add_argument('-JAB1_fn', dest='JAB1_fn', type=float, nargs=1, default=[envn.JAB1_fn],
                       help='Function of JAB1 signalling, JAB1_fn', required=False);
    parser.add_argument('-Bcl2_Bclxl_fn', dest='Bcl2_Bclxl_fn', type=float, nargs=1, default=[envn.Bcl2_Bclxl_fn],
                       help='Function of Bcl2_Bclxl signalling, Bcl2_Bclxl_fn', required=False);
    parser.add_argument('-BAK_BAX_fn', dest='BAK_BAX_fn', type=float, nargs=1, default=[envn.BAK_BAX_fn],
                       help='Function of BAK_BAX signalling, BAK_BAX_fn', required=False);
    parser.add_argument('-CAS_fn', dest='CAS_fn', type=float, nargs=1, default=[envn.CAS_fn],
                       help='Function of CAS signalling, CAS_fn', required=False);
    parser.add_argument('-VEGF_fn', dest='VEGF_fn', type=float, nargs=1, default=[envn.VEGF_fn],
                       help='Function of VEGF signalling, VEGF_fn', required=False);
    parser.add_argument('-cycle', dest='cycle', type=float, nargs=1, default=[envn.cycle],
                       help='Cell cycle progress, cycle', required=False);
    parser.add_argument('-apop', dest='apop', type=int, nargs=1, default=[envn.apop],
                       help='Apoptotic status, apop', required=False);    
    parser.add_argument('-apop_signal', dest='apop_signal', type=int, nargs=1, default=[envn.apop_signal],
                       help='Number of apoptotic signals, apop_signal', required=False);
    parser.add_argument('-necro', dest='necro', type=int, nargs=1, default=[envn.necro],
                       help='Necrotic status, necro', required=False); 
    parser.add_argument('-necro_signal', dest='necro_signal', type=int, nargs=1, default=[envn.necro_signal],
                       help='Number of necrotic signals, necro_signal', required=False);                 
    parser.add_argument('-telo_count', dest='telo_count', type=int, nargs=1, default=[envn.telo_count],
                       help='Number of telomere units, telo_count', required=False);
    """
    Override the default initialisation routine for Schwann cell agents.
    """
    parser.add_argument('-cycle_sc', dest='cycle_sc', type=float, nargs=1, default=[envn.cycle_sc],
                       help='Cell cycle progress in a Schwann cell, cycle_sc', required=False);
    parser.add_argument('-apop_sc', dest='apop_sc', type=int, nargs=1, default=[envn.apop_sc],
                       help='Apoptotic status in a Schwann cell, apop_sc', required=False);
    parser.add_argument('-apop_signal_sc', dest='apop_signal_sc', type=int, nargs=1, default=[envn.apop_signal_sc],
                       help='Number of apoptotic signals in a Schwann cell, apop_signal_sc', required=False);
    parser.add_argument('-necro_sc', dest='necro_sc', type=int, nargs=1, default=[envn.necro_sc],
                       help='Necrotic status in a Schwann cell, necro_sc', required=False);
    parser.add_argument('-necro_signal_sc', dest='necro_signal_sc', type=int, nargs=1, default=[envn.necro_signal_sc],
                       help='Number of necrotic signals in a Schwann cell, necro_signal_sc', required=False);              
    parser.add_argument('-telo_count_sc', dest='telo_count_sc', type=int, nargs=1, default=[envn.telo_count_sc],
                       help='Number of telomere units in a Schwann cell, telo_count_sc', required=False);
    """
    Override environmental properties.
    """
    env = Environment();
    for key, val in env.__dict__.items():
        if isinstance(val, list):
            parser.add_argument('-env_%s'%(key), metavar='_', dest=key, type=float, nargs=len(val), default=val,
                       help='Environment.descriptions[key]', 
                       required=False);
        else:
            parser.add_argument('-env_%s'%(key), metavar='_', dest=key, type=float, nargs=1, default=[val],
                       help='Environment.descriptions[key]', 
                       required=False);             
    #Output file
    #Filetype required
    #parser.add_argument('out', nargs=1,
    #                   help='The output file.', required=False);      
    return parser;

"""
Define the cell type–independent functions.
"""
def vasculature(P_O2v, ang_signal, matrix_grid, Nnbl_grid, Nscl_grid, N_vf, t):
    """
    Case 1:
    Set up the initial vasculature in terms of the amount of oxygen it supplies in one time step.
    The assumption is that it can supply the amount consumed by the initial population of living neuroblasts and Schwann cells in one time step.
    Case 2:
    When there are more living VEGF-producing neuroblasts than living Schwann cells, an angiogenic signal is produced.
    When there are enough angiogenic signals, calculate the amount of oxygen consumed by the current population of living cells.
    Take this as the new oxygen supply rate if it exceeds the old rate.
    """
    if t == 0:
        [span_xr, span_yr, span_zr] = np.shape(matrix_grid) 
        N_voxel = span_xr*span_yr*span_zr
        P_O2v = -1*(envn.N_cell()+envn.N_scell())*envn.P_O20*step_size/(N_voxel*envn.V_grid()*1e-15)/envn.Cs_O2
        ang_signal = 0
    elif N_vf > np.sum(Nscl_grid):
        ang_signal += 1*step_size
    if ang_signal == envn.ang_critical:
        [span_xr, span_yr, span_zr] = np.shape(matrix_grid) 
        N_voxel = span_xr*span_yr*span_zr
        dummy_P_O2v = -1*(np.sum(Nnbl_grid)+np.sum(Nscl_grid))*envn.P_O20*step_size/(N_voxel*envn.V_grid()*1e-15)/envn.Cs_O2
        if dummy_P_O2v > P_O2v:
            P_O2v = dummy_P_O2v
        ang_signal = 0
    return P_O2v, ang_signal

def fresolve(nblist, sclist, matrix_grid):
    """
    Resolve the repulsive forces within the entire cell population by cell migration.
    """
    iteration_resforce = 0
    total_overlap = 0
    list_overlap = []
    list_neighbours = []
    total_force = 0
    #max_displace = 0
    for i in nblist:
        [i.Fx, i.Fy, i.Fz] = [0, 0, 0]
        i.overlap = 0
        i.neighbours = 0
    for i in sclist:
        [i.Fx, i.Fy, i.Fz] = [0, 0, 0]
        i.overlap = 0
        i.neighbours = 0
    for i in nblist:
        if i.cycle < 1:
            Ri = i.cycle*envn.cycle_stages[0]/(envn.cycle_stages[0] + envn.cycle_stages[2])
        elif i.cycle < 2:
            Ri = envn.cycle_stages[0]/(envn.cycle_stages[0] + envn.cycle_stages[2])
        elif i.cycle < 3:
            Ri = (envn.cycle_stages[0] + (i.cycle-2)*envn.cycle_stages[2])/(envn.cycle_stages[0] + envn.cycle_stages[2])
        else:
            Ri = (envn.cycle_stages[0] + envn.cycle_stages[2])/(envn.cycle_stages[0] + envn.cycle_stages[2])
        for j in nblist:
            displacement_ij = np.array([i.x-j.x, i.y-j.y, i.z-j.z])
            if np.linalg.norm(displacement_ij) < envn.R_neighbours:
                i.neighbours += 1
            if j.cycle < 1:
                Rj = j.cycle*envn.cycle_stages[0]/(envn.cycle_stages[0] + envn.cycle_stages[2])
            elif j.cycle < 2:
                Rj = envn.cycle_stages[0]/(envn.cycle_stages[0] + envn.cycle_stages[2])
            elif j.cycle < 3:
                Rj = (envn.cycle_stages[0] + (j.cycle-2)*envn.cycle_stages[2])/(envn.cycle_stages[0] + envn.cycle_stages[2])
            else:
                Rj = (envn.cycle_stages[0] + envn.cycle_stages[2])/(envn.cycle_stages[0] + envn.cycle_stages[2])
            overlap_ij = 2*envn.R_cell + (Ri+Rj)*envn.R_cell - np.linalg.norm(displacement_ij)
            if overlap_ij == 2*envn.R_cell + (Ri+Rj)*envn.R_cell:
                overlap_ij = 0
                force_vector = [0, 0, 0]
            elif overlap_ij < envn.min_overlap:
                Fij = 0
                force_vector = Fij*displacement_ij/np.linalg.norm(displacement_ij)
            else:
                Fij = envn.k1*overlap_ij
                force_vector = Fij*displacement_ij/np.linalg.norm(displacement_ij)
                i.overlap += overlap_ij
            i.Fx += force_vector[0]
            i.Fy += force_vector[1]
            i.Fz += force_vector[2]
        for j in sclist:
            displacement_ij = np.array([i.x-j.x, i.y-j.y, i.z-j.z])
            if np.linalg.norm(displacement_ij) < envn.R_neighbours:
                i.neighbours += 1            
            if j.cycle < 1:
                Rj = j.cycle*envn.cycle_stages[0]/(envn.cycle_stages[0] + envn.cycle_stages[2])
            elif j.cycle < 2:
                Rj = envn.cycle_stages[0]/(envn.cycle_stages[0] + envn.cycle_stages[2])
            elif j.cycle < 3:
                Rj = (envn.cycle_stages[0] + (j.cycle-2)*envn.cycle_stages[2])/(envn.cycle_stages[0] + envn.cycle_stages[2])
            else:
                Rj = (envn.cycle_stages[0] + envn.cycle_stages[2])/(envn.cycle_stages[0] + envn.cycle_stages[2])
            overlap_ij = 2*envn.R_cell + (Ri+Rj)*envn.R_cell - np.linalg.norm(displacement_ij)
            if overlap_ij < envn.min_overlap:
                Fij = 0
                force_vector = Fij*displacement_ij/np.linalg.norm(displacement_ij)
            else:
                Fij = envn.k1*overlap_ij
                force_vector = Fij*displacement_ij/np.linalg.norm(displacement_ij)
                i.overlap += overlap_ij
            i.Fx += force_vector[0]
            i.Fy += force_vector[1]
            i.Fz += force_vector[2]         
        if i.neighbours > envn.N_neighbours and i.mobile == 1:
            i.Fx = envn.k_locom*i.Fx
            i.Fy = envn.k_locom*i.Fy
            i.Fz = envn.k_locom*i.Fz
        total_overlap += i.overlap
        total_force += np.linalg.norm(np.array([i.Fx, i.Fy, i.Fz]))
    for i in sclist:
        if i.cycle < 1:
            Ri = i.cycle*envn.cycle_stages[0]/(envn.cycle_stages[0] + envn.cycle_stages[2])
        elif i.cycle < 2:
            Ri = envn.cycle_stages[0]/(envn.cycle_stages[0] + envn.cycle_stages[2])
        elif i.cycle < 3:
            Ri = (envn.cycle_stages[0] + (i.cycle-2)*envn.cycle_stages[2])/(envn.cycle_stages[0] + envn.cycle_stages[2])
        else:
            Ri = (envn.cycle_stages[0] + envn.cycle_stages[2])/(envn.cycle_stages[0] + envn.cycle_stages[2])
        for j in sclist:
            displacement_ij = np.array([i.x-j.x, i.y-j.y, i.z-j.z])
            if np.linalg.norm(displacement_ij) < envn.R_neighbours:
                i.neighbours += 1            
            if j.cycle < 1:
                Rj = j.cycle*envn.cycle_stages[0]/(envn.cycle_stages[0] + envn.cycle_stages[2])
            elif j.cycle < 2:
                Rj = envn.cycle_stages[0]/(envn.cycle_stages[0] + envn.cycle_stages[2])
            elif j.cycle < 3:
                Rj = (envn.cycle_stages[0] + (j.cycle-2)*envn.cycle_stages[2])/(envn.cycle_stages[0] + envn.cycle_stages[2])
            else:
                Rj = (envn.cycle_stages[0] + envn.cycle_stages[2])/(envn.cycle_stages[0] + envn.cycle_stages[2])
            overlap_ij = 2*envn.R_cell + (Ri+Rj)*envn.R_cell - np.linalg.norm(displacement_ij)
            if overlap_ij == 2*envn.R_cell + (Ri+Rj)*envn.R_cell:
                overlap_ij = 0
                force_vector = [0, 0, 0]
            elif overlap_ij < envn.min_overlap:
                Fij = 0
                force_vector = Fij*displacement_ij/np.linalg.norm(displacement_ij)
            else:
                Fij = envn.k1*overlap_ij
                force_vector = Fij*displacement_ij/np.linalg.norm(displacement_ij)
                i.overlap += overlap_ij
            i.Fx += force_vector[0]
            i.Fy += force_vector[1]
            i.Fz += force_vector[2]
        for j in nblist:
            displacement_ij = np.array([i.x-j.x, i.y-j.y, i.z-j.z])
            if np.linalg.norm(displacement_ij) < envn.R_neighbours:
                i.neighbours += 1            
            if j.cycle < 1:
                Rj = j.cycle*envn.cycle_stages[0]/(envn.cycle_stages[0] + envn.cycle_stages[2])
            elif j.cycle < 2:
                Rj = envn.cycle_stages[0]/(envn.cycle_stages[0] + envn.cycle_stages[2])
            elif j.cycle < 3:
                Rj = (envn.cycle_stages[0] + (j.cycle-2)*envn.cycle_stages[2])/(envn.cycle_stages[0] + envn.cycle_stages[2])
            else:
                Rj = (envn.cycle_stages[0] + envn.cycle_stages[2])/(envn.cycle_stages[0] + envn.cycle_stages[2])
            overlap_ij = 2*envn.R_cell + (Ri+Rj)*envn.R_cell - np.linalg.norm(displacement_ij)
            if overlap_ij < envn.min_overlap:
                Fij = 0
                force_vector = Fij*displacement_ij/np.linalg.norm(displacement_ij)
            else:
                Fij = envn.k1*overlap_ij
                force_vector = Fij*displacement_ij/np.linalg.norm(displacement_ij)
                i.overlap += overlap_ij
            i.Fx += force_vector[0]
            i.Fy += force_vector[1]
            i.Fz += force_vector[2]
        if i.neighbours > envn.N_neighbours and i.mobile == 1:
            i.Fx = envn.k_locom*i.Fx
            i.Fy = envn.k_locom*i.Fy
            i.Fz = envn.k_locom*i.Fz          
        total_overlap += i.overlap      
        total_force += np.linalg.norm(np.array([i.Fx, i.Fy, i.Fz]))
    for i in nblist:
        list_overlap.append(i.overlap)
        list_neighbours.append(i.neighbours)
    for i in sclist:
        list_overlap.append(i.overlap)
        list_neighbours.append(i.neighbours)
    max_overlap = max(list_overlap)
    min_overlap_log = min(list_overlap)
    max_neighbours = max(list_neighbours)
    median_neighbours = stat.median(list_neighbours)
    mean_neighbours = stat.mean(list_neighbours)
    min_neighbours = min(list_neighbours)
    print('Maximum overlap per cell: {}.'.format(max_overlap))
    print('Minimum overlap per cell: {}.'.format(min_overlap_log))
    print('Maximum number of neighbours per cell, including itself: {}.'.format(max_neighbours))
    print('Median number of neighbours per cell, including itself: {}.'.format(median_neighbours))
    print('Mean number of neighbours per cell, including itself: {}.'.format(mean_neighbours))
    print('Minimum number of neighbours per cell, including itself: {}.'.format(min_neighbours))
    print('========================================================================')
    while True:
        dummy_overlap = 0
        list_overlap = []
        list_neighbours = []
        dummy_force = 0
        [span_xr, span_yr, span_zr] = np.shape(matrix_grid)
        span_x = span_xr*envn.R_voxel*2
        span_y = span_yr*envn.R_voxel*2
        span_z = span_zr*envn.R_voxel*2
        for i in nblist:  
            x_unit = math.floor((i.x + span_x/2)/envn.R_voxel/2)
            y_unit = math.floor((i.y + span_y/2)/envn.R_voxel/2)
            z_unit = math.floor((i.z + span_z/2)/envn.R_voxel/2)  
            dummy_matrix = matrix_grid[x_unit, y_unit, z_unit]
            mu_eff = (1+dummy_matrix)*envn.mu
            i.x += i.Fx*envn.dt/mu_eff
            i.y += i.Fy*envn.dt/mu_eff         
            i.z += i.Fz*envn.dt/mu_eff
            #dummy_displace = math.sqrt((i.Fx*envn.dt/mu_eff)**2+(i.Fy*envn.dt/mu_eff)**2+(i.Fz*envn.dt/mu_eff)**2)
            #if dummy_displace > max_displace:
            #    max_displace = dummy_displace
            boundary_conditions(i)
            [i.Fx, i.Fy, i.Fz] = [0, 0, 0]   
            i.overlap = 0
            i.neighbours = 0
        for i in sclist:
            x_unit = math.floor((i.x + span_x/2)/envn.R_voxel/2)
            y_unit = math.floor((i.y + span_y/2)/envn.R_voxel/2)
            z_unit = math.floor((i.z + span_z/2)/envn.R_voxel/2)  
            dummy_matrix = matrix_grid[x_unit, y_unit, z_unit]
            mu_eff = (1+dummy_matrix)*envn.mu       
            i.x += i.Fx*envn.dt/mu_eff
            i.y += i.Fy*envn.dt/mu_eff         
            i.z += i.Fz*envn.dt/mu_eff
            #dummy_displace = math.sqrt((i.Fx*envn.dt/mu_eff)**2+(i.Fy*envn.dt/mu_eff)**2+(i.Fz*envn.dt/mu_eff)**2)
            #if dummy_displace > max_displace:
            #    max_displace = dummy_displace
            boundary_conditions(i)
            [i.Fx, i.Fy, i.Fz] = [0, 0, 0]
            i.overlap = 0
            i.neighbours = 0
        iteration_resforce += 1            
        [x, y, z, xs, ys, zs, xa, ya, za, xsa, ysa, zsa, xn, yn, zn, xsn, ysn, zsn, xl, yl, zl, xsl, ysl, zsl, xsl_col, ysl_col, zsl_col, N_vf, V, Vs] = locate(nblist, sclist)
        [matrix_grid, N_grid, Nnb_grid, Nnba_grid, Nnbn_grid, Nnbl_grid, Nsc_grid, Nsca_grid, Nscn_grid, Nscl_grid, Nscl_col_grid] = CAexpand(matrix_grid, x, y, z, xs, ys, zs, xa, ya, za, xsa, ysa, zsa, xn, yn, zn, xsn, ysn, zsn, xl, yl, zl, xsl, ysl, zsl, xsl_col, ysl_col, zsl_col, V, Vs)
        for i in nblist:
            if i.cycle < 1:
                Ri = i.cycle*envn.cycle_stages[0]/(envn.cycle_stages[0] + envn.cycle_stages[2])
            elif i.cycle < 2:
                Ri = envn.cycle_stages[0]/(envn.cycle_stages[0] + envn.cycle_stages[2])
            elif i.cycle < 3:
                Ri = (envn.cycle_stages[0] + (i.cycle-2)*envn.cycle_stages[2])/(envn.cycle_stages[0] + envn.cycle_stages[2])
            else:
                Ri = (envn.cycle_stages[0] + envn.cycle_stages[2])/(envn.cycle_stages[0] + envn.cycle_stages[2])            
            for j in nblist:
                displacement_ij = np.array([i.x-j.x, i.y-j.y, i.z-j.z])
                if np.linalg.norm(displacement_ij) < envn.R_neighbours:
                    i.neighbours += 1                
                if j.cycle < 1:
                    Rj = j.cycle*envn.cycle_stages[0]/(envn.cycle_stages[0] + envn.cycle_stages[2])
                elif j.cycle < 2:
                    Rj = envn.cycle_stages[0]/(envn.cycle_stages[0] + envn.cycle_stages[2])
                elif j.cycle < 3:
                    Rj = (envn.cycle_stages[0] + (j.cycle-2)*envn.cycle_stages[2])/(envn.cycle_stages[0] + envn.cycle_stages[2])
                else:
                    Rj = (envn.cycle_stages[0] + envn.cycle_stages[2])/(envn.cycle_stages[0] + envn.cycle_stages[2])
                overlap_ij = 2*envn.R_cell + (Ri+Rj)*envn.R_cell - np.linalg.norm(displacement_ij)
                if overlap_ij == 2*envn.R_cell + (Ri+Rj)*envn.R_cell:
                    overlap_ij = 0
                    force_vector = [0, 0, 0]
                elif overlap_ij < envn.min_overlap:
                    Fij = 0
                    force_vector = Fij*displacement_ij/np.linalg.norm(displacement_ij)
                else:
                    Fij = envn.k1*overlap_ij
                    force_vector = Fij*displacement_ij/np.linalg.norm(displacement_ij)
                    i.overlap += overlap_ij
                i.Fx += force_vector[0]
                i.Fy += force_vector[1]
                i.Fz += force_vector[2]
            for j in sclist:
                displacement_ij = np.array([i.x-j.x, i.y-j.y, i.z-j.z])
                if np.linalg.norm(displacement_ij) < envn.R_neighbours:
                    i.neighbours += 1  
                if j.cycle < 1:
                    Rj = j.cycle*envn.cycle_stages[0]/(envn.cycle_stages[0] + envn.cycle_stages[2])
                elif j.cycle < 2:
                    Rj = envn.cycle_stages[0]/(envn.cycle_stages[0] + envn.cycle_stages[2])
                elif j.cycle < 3:
                    Rj = (envn.cycle_stages[0] + (j.cycle-2)*envn.cycle_stages[2])/(envn.cycle_stages[0] + envn.cycle_stages[2])
                else:
                    Rj = (envn.cycle_stages[0] + envn.cycle_stages[2])/(envn.cycle_stages[0] + envn.cycle_stages[2])
                overlap_ij = 2*envn.R_cell + (Ri+Rj)*envn.R_cell - np.linalg.norm(displacement_ij)
                if overlap_ij < envn.min_overlap:
                    Fij = 0
                    force_vector = Fij*displacement_ij/np.linalg.norm(displacement_ij)
                else:
                    Fij = envn.k1*overlap_ij
                    force_vector = Fij*displacement_ij/np.linalg.norm(displacement_ij)
                    i.overlap += overlap_ij
                i.Fx += force_vector[0]
                i.Fy += force_vector[1]
                i.Fz += force_vector[2]
            if i.neighbours > envn.N_neighbours and i.mobile == 1:
                i.Fx = envn.k_locom*i.Fx
                i.Fy = envn.k_locom*i.Fy
                i.Fz = envn.k_locom*i.Fz
            dummy_overlap += i.overlap           
            dummy_force += np.linalg.norm(np.array([i.Fx, i.Fy, i.Fz]))
        for i in sclist:
            if i.cycle < 1:
                Ri = i.cycle*envn.cycle_stages[0]/(envn.cycle_stages[0] + envn.cycle_stages[2])
            elif i.cycle < 2:
                Ri = envn.cycle_stages[0]/(envn.cycle_stages[0] + envn.cycle_stages[2])
            elif i.cycle < 3:
                Ri = (envn.cycle_stages[0] + (i.cycle-2)*envn.cycle_stages[2])/(envn.cycle_stages[0] + envn.cycle_stages[2])
            else:
                Ri = (envn.cycle_stages[0] + envn.cycle_stages[2])/(envn.cycle_stages[0] + envn.cycle_stages[2])
            for j in sclist:
                displacement_ij = np.array([i.x-j.x, i.y-j.y, i.z-j.z])
                if np.linalg.norm(displacement_ij) < envn.R_neighbours:
                    i.neighbours += 1                  
                if j.cycle < 1:
                    Rj = j.cycle*envn.cycle_stages[0]/(envn.cycle_stages[0] + envn.cycle_stages[2])
                elif j.cycle < 2:
                    Rj = envn.cycle_stages[0]/(envn.cycle_stages[0] + envn.cycle_stages[2])
                elif j.cycle < 3:
                    Rj = (envn.cycle_stages[0] + (j.cycle-2)*envn.cycle_stages[2])/(envn.cycle_stages[0] + envn.cycle_stages[2])
                else:
                    Rj = (envn.cycle_stages[0] + envn.cycle_stages[2])/(envn.cycle_stages[0] + envn.cycle_stages[2])
                overlap_ij = 2*envn.R_cell + (Ri+Rj)*envn.R_cell - np.linalg.norm(displacement_ij)
                if overlap_ij == 2*envn.R_cell + (Ri+Rj)*envn.R_cell:
                    overlap_ij = 0
                    force_vector = [0, 0, 0]
                elif overlap_ij < envn.min_overlap:
                    Fij = 0
                    force_vector = Fij*displacement_ij/np.linalg.norm(displacement_ij)
                else:
                    Fij = envn.k1*overlap_ij
                    force_vector = Fij*displacement_ij/np.linalg.norm(displacement_ij)
                    i.overlap += overlap_ij
                i.Fx += force_vector[0]
                i.Fy += force_vector[1]
                i.Fz += force_vector[2]
            for j in nblist:
                displacement_ij = np.array([i.x-j.x, i.y-j.y, i.z-j.z])
                if np.linalg.norm(displacement_ij) < envn.R_neighbours:
                    i.neighbours += 1                  
                if j.cycle < 1:
                    Rj = j.cycle*envn.cycle_stages[0]/(envn.cycle_stages[0] + envn.cycle_stages[2])
                elif j.cycle < 2:
                    Rj = envn.cycle_stages[0]/(envn.cycle_stages[0] + envn.cycle_stages[2])
                elif j.cycle < 3:
                    Rj = (envn.cycle_stages[0] + (j.cycle-2)*envn.cycle_stages[2])/(envn.cycle_stages[0] + envn.cycle_stages[2])
                else:
                    Rj = (envn.cycle_stages[0] + envn.cycle_stages[2])/(envn.cycle_stages[0] + envn.cycle_stages[2])
                overlap_ij = 2*envn.R_cell + (Ri+Rj)*envn.R_cell - np.linalg.norm(displacement_ij)
                if overlap_ij < envn.min_overlap:
                    Fij = 0
                    force_vector = Fij*displacement_ij/np.linalg.norm(displacement_ij)
                else:
                    Fij = envn.k1*overlap_ij
                    force_vector = Fij*displacement_ij/np.linalg.norm(displacement_ij)
                    i.overlap += overlap_ij
                i.Fx += force_vector[0]
                i.Fy += force_vector[1]
                i.Fz += force_vector[2]               
            if i.neighbours > envn.N_neighbours and i.mobile == 1:
                i.Fx = envn.k_locom*i.Fx
                i.Fy = envn.k_locom*i.Fy
                i.Fz = envn.k_locom*i.Fz           
            dummy_overlap += i.overlap
            dummy_force += np.linalg.norm(np.array([i.Fx, i.Fy, i.Fz]))            
        for i in nblist:
            list_overlap.append(i.overlap)
            list_neighbours.append(i.neighbours)
        for i in sclist:
            list_overlap.append(i.overlap)
            list_neighbours.append(i.neighbours)
        max_overlap = max(list_overlap)
        min_overlap_log = min(list_overlap)
        max_neighbours = max(list_neighbours)
        median_neighbours = stat.median(list_neighbours)
        mean_neighbours = stat.mean(list_neighbours)
        min_neighbours = min(list_neighbours)
        print('Maximum overlap per cell: {}.'.format(max_overlap))
        print('Minimum overlap per cell: {}.'.format(min_overlap_log))
        print('Maximum number of neighbours per cell, including itself: {}.'.format(max_neighbours))
        print('Median number of neighbours per cell, including itself: {}.'.format(median_neighbours))
        print('Mean number of neighbours per cell, including itself: {}.'.format(mean_neighbours))
        print('Minimum number of neighbours per cell, including itself: {}.'.format(min_neighbours))
        print('========================================================================')            
        if (len(nblist) + len(sclist)) < 2:
            #print(max_displace)
            break
        elif max_neighbours <= envn.N_neighbours:
            #print(max_displace)
            break
        elif max_overlap < 0.15*envn.R_cell:
            #print(max_displace)
            break
        elif iteration_resforce > step_size*3600/envn.dt:
            #print(max_displace)
            break
        else:
            total_overlap = dummy_overlap

def locate(nblist, sclist):
    """
    Extract key metrics from the cell populations.
    1. Locate the three types of neuroblasts and Schwann cells (living, apoptotic, and necrotic), 
    as well as the two populations as a whole and the collagen-producing Schwann cells.
    2. Record the volume of each cell from each population.
    3. Count the VEGF-producing neuroblasts.
    """
    x = []
    y = []
    z = []   
    xa = []
    ya = []
    za = []
    xn = []
    yn = []
    zn = []
    xl = []
    yl = []
    zl = []    
    xs = []
    ys = []
    zs = []    
    xsa = []
    ysa = []
    zsa = []
    xsn = []
    ysn = []
    zsn = [] 
    xsl = []
    ysl = []
    zsl = []
    xsl_col = []
    ysl_col = []
    zsl_col = []
    N_vf = 0
    V = []
    Vs = []
    for i in nblist:
        x.append(i.x)
        y.append(i.y)
        z.append(i.z)
        if i.cycle < 1:
            Ri = i.cycle*envn.cycle_stages[0]/(envn.cycle_stages[0] + envn.cycle_stages[2])
        elif i.cycle < 2:
            Ri = envn.cycle_stages[0]/(envn.cycle_stages[0] + envn.cycle_stages[2])
        elif i.cycle < 3:
            Ri = (envn.cycle_stages[0] + (i.cycle-2)*envn.cycle_stages[2])/(envn.cycle_stages[0] + envn.cycle_stages[2])
        else:
            Ri = (envn.cycle_stages[0] + envn.cycle_stages[2])/(envn.cycle_stages[0] + envn.cycle_stages[2])
        dummy_R = (1+Ri)*envn.R_cell
        V.append((4*math.pi*dummy_R**3)/3)
        if i.apop == 1:
            xa.append(i.x)
            ya.append(i.y)
            za.append(i.z)
        elif i.necro == 1:
            xn.append(i.x)
            yn.append(i.y)
            zn.append(i.z)
        else:
            xl.append(i.x)
            yl.append(i.y)
            zl.append(i.z)
            if i.VEGF == 1:
                N_vf += 1
    for i in sclist:
        xs.append(i.x)
        ys.append(i.y)
        zs.append(i.z)
        if i.cycle < 1:
            Ri = i.cycle*envn.cycle_stages[0]/(envn.cycle_stages[0] + envn.cycle_stages[2])
        elif i.cycle < 2:
            Ri = envn.cycle_stages[0]/(envn.cycle_stages[0] + envn.cycle_stages[2])
        elif i.cycle < 3:
            Ri = (envn.cycle_stages[0] + (i.cycle-2)*envn.cycle_stages[2])/(envn.cycle_stages[0] + envn.cycle_stages[2])
        else:
            Ri = (envn.cycle_stages[0] + envn.cycle_stages[2])/(envn.cycle_stages[0] + envn.cycle_stages[2])
        dummy_R = (1+Ri)*envn.R_cell
        Vs.append((4*math.pi*dummy_R**3)/3)
        if i.apop == 1:
            xsa.append(i.x)
            ysa.append(i.y)
            zsa.append(i.z)
        elif i.necro == 1:
            xsn.append(i.x)
            ysn.append(i.y)
            zsn.append(i.z)
        else:
            xsl.append(i.x)
            ysl.append(i.y)
            zsl.append(i.z)
            if i.neighbours <= envn.N_neighbours:
                xsl_col.append(i.x)
                ysl_col.append(i.y)
                zsl_col.append(i.z)                
    return x, y, z, xs, ys, zs, xa, ya, za, xsa, ysa, zsa, xn, yn, zn, xsn, ysn, zsn, xl, yl, zl, xsl, ysl, zsl, xsl_col, ysl_col, zsl_col, N_vf, V, Vs

def CAexpand(matrix_grid, x, y, z, xs, ys, zs, xa, ya, za, xsa, ysa, zsa, xn, yn, zn, xsn, ysn, zsn, xl, yl, zl, xsl, ysl, zsl, xsl_col, ysl_col, zsl_col, V, Vs):
    """
    Expand the grids to reflect tumour growth.
    1. Check if the tumour has expanded in each direction.
    2. If so, expand the matrix grid in that direction, using the average quantity in the existing grid to initialise the new voxels.
    3. Expand the grids for various cell types.
    """
    [oldspan_xr, oldspan_yr, oldspan_zr] = np.shape(matrix_grid) 
    oldspan_x = math.ceil(oldspan_xr/2)    
    oldspan_y = math.ceil(oldspan_yr/2) 
    oldspan_z = math.ceil(oldspan_zr/2)
    newspan_x = math.ceil((max(np.abs(x+xs))/envn.R_voxel/2) + 0.5)
    newspan_y = math.ceil((max(np.abs(y+ys))/envn.R_voxel/2) + 0.5)
    newspan_z = math.ceil((max(np.abs(z+zs))/envn.R_voxel/2) + 0.5)
    matrix_dummy = np.mean(matrix_grid)
    if newspan_x > oldspan_x:
        length_pad = newspan_x-oldspan_x
        matrix_grid = np.pad(matrix_grid, ((length_pad, length_pad), (0,0), (0,0)), 'constant', constant_values=matrix_dummy)
    if newspan_y > oldspan_y:
        length_pad = newspan_y-oldspan_y
        matrix_grid = np.pad(matrix_grid, ((0,0), (length_pad, length_pad), (0,0)), 'constant', constant_values=matrix_dummy)
    if newspan_z > oldspan_z:
        length_pad = newspan_z-oldspan_z
        matrix_grid = np.pad(matrix_grid, ((0,0), (0,0), (length_pad, length_pad)), 'constant', constant_values=matrix_dummy)
    Nnb_grid = CAcounting(x, y, z, matrix_grid)
    Nnba_grid = CAcounting(xa, ya, za, matrix_grid)
    Nnbn_grid = CAcounting(xn, yn, zn, matrix_grid)
    Nnbl_grid = CAcounting(xl, yl, zl, matrix_grid)    
    Nsc_grid = CAcounting(xs, ys, zs, matrix_grid)
    Nsca_grid = CAcounting(xsa, ysa, zsa, matrix_grid)
    Nscn_grid = CAcounting(xsn, ysn, zsn, matrix_grid)
    Nscl_grid = CAcounting(xsl, ysl, zsl, matrix_grid)
    Nscl_col_grid = CAcounting(xsl_col, ysl_col, zsl_col, matrix_grid)
    N_grid = Nnb_grid + Nsc_grid  
    return matrix_grid, N_grid, Nnb_grid, Nnba_grid, Nnbn_grid, Nnbl_grid, Nsc_grid, Nsca_grid, Nscn_grid, Nscl_grid, Nscl_col_grid

def CAcounting(x, y, z, grid):
    """
    Build the grid for a particular cell type by locating and counting cells of this type.
    """
    dummy_grid = np.zeros(np.shape(grid))
    [span_xr, span_yr, span_zr] = np.shape(grid)
    span_x = span_xr*envn.R_voxel*2
    span_y = span_yr*envn.R_voxel*2
    span_z = span_zr*envn.R_voxel*2
    for i in range(len(x)):    
        x_unit = math.floor((x[i] + span_x/2)/envn.R_voxel/2)
        y_unit = math.floor((y[i] + span_y/2)/envn.R_voxel/2)
        z_unit = math.floor((z[i] + span_z/2)/envn.R_voxel/2)  
        dummy_grid[x_unit, y_unit, z_unit] += 1
    return dummy_grid

def record(t, x, y, z, xs, ys, zs, xa, ya, za, xsa, ysa, zsa, xn, yn, zn, xsn, ysn, zsn, xl, yl, zl, xsl, ysl, zsl, O2, matrix_grid, N_grid, Nnb_grid, Nnba_grid, Nnbn_grid, Nnbl_grid, Nsc_grid, Nsca_grid, Nscn_grid, Nscl_grid, ang_signal):
    """
    Record the outputs after t time steps.
    """
    print('========================================================================')
    print('Time: {} hours.'.format(t))
    print('1. Total neuroblast count: {}.'.format(len(x)))
    print('-Apoptotic neuroblast count: {}.'.format(len(xa)))
    print('-Necrotic neuroblast count: {}.'.format(len(xn)))
    print('-Living neuroblast count: {}.'.format(len(xl)))
    print('2. Total Schwann cell count: {}.'.format(len(xs)))
    print('-Apoptotic Schwann cell count: {}.'.format(len(xsa)))
    print('-Necrotic Schwann cell count: {}.'.format(len(xsn)))
    print('-Living Schwann cell count: {}.'.format(len(xsl)))    
    print('3. Oxygen level, scaled: {}.'.format(O2))
    print('4. Average matrix volume fraction: {}.'.format(np.sum(matrix_grid)/np.size(matrix_grid)))
    print('5. Minimum cell count per voxel: {}.'.format(np.min(N_grid)))
    print('6. Maximum cell count per voxel: {}.'.format(np.max(N_grid)))
    print('7. Grid dimensions: {}.'.format(np.shape(matrix_grid)))
    [dummyx, dummyy, dummyz] = np.shape(matrix_grid) 
    V = dummyx*dummyy*dummyz*envn.V_grid()
    print('8. Tumour volume, naive estimate: {} cubic microns.'.format(V))
    print('9. Tumour volume, naive estimate: {} cubic mm.'.format(V*1e-9))
    print('10. Number of angiogenic signals: {}.'.format(ang_signal))
    print('========================================================================')
    return t, x, y, z, xs, ys, zs, xa, ya, za, xsa, ysa, zsa, xn, yn, zn, xsn, ysn, zsn, xl, yl, zl, xsl, ysl, zsl, O2, matrix_grid, N_grid, Nnb_grid, Nnba_grid, Nnbn_grid, Nnbl_grid, Nsc_grid, Nsca_grid, Nscn_grid, Nscl_grid, V

def remove(self):
    """
    Remove an apoptotic or necrotic cell if it is engulfed by an immune cell.
    """
    if self.apop == 1 or self.necro == 1:
        if np.random.uniform(0, 1) < envn.P_lysis*step_size:
            return True
        else:
            return False
    else:
        return False

def drift(self, R_cell):
    """
    Randomly and slightly move a generic cell.
    """
    dummy_dir = [np.random.uniform(-1, 1), np.random.uniform(-1, 1) , np.random.uniform(-1, 1)]
    norm_dir = dummy_dir/np.linalg.norm(dummy_dir)
    self.x += 2*R_cell * norm_dir[0]
    self.y += 2*R_cell * norm_dir[1]
    self.z += 2*R_cell * norm_dir[2]

def boundary_conditions(self):
    """
    If a cell migrates into the tumour boundary, it will be moved in the opposite direction.
    Assume fixed boundary conditions throughout the simulation.
    """
    if self.x < envn.x_bc_minus():
        #print('Invasion into the -x boundary.')
        self.x += envn.x_displace*(envn.x_bc_minus()-self.x) #The magnitude depends on the extent of invasion and the field's sensitivity to the extent.
    elif self.x > envn.x_bc_plus():
        #print('Invasion into the +x boundary.')        
        self.x -= envn.x_displace*(self.x-envn.x_bc_plus()) #The magnitude depends on the extent of invasion and the field's sensitivity to the extent.
    if self.y < envn.y_bc_minus():
        #print('Invasion into the -y boundary.')        
        self.y += envn.y_displace*(envn.y_bc_minus()-self.y) #The magnitude depends on the extent of invasion and the field's sensitivity to the extent.
    elif self.y > envn.y_bc_plus():
        #print('Invasion into the +y boundary.')            
        self.y -= envn.y_displace*(self.y-envn.y_bc_plus()) #The magnitude depends on the extent of invasion and the field's sensitivity to the extent.
    if self.z < envn.z_bc_minus():
        #print('Invasion into the -z boundary.')        
        self.z += envn.z_displace*(envn.z_bc_minus()-self.z) #The magnitude depends on the extent of invasion and the field's sensitivity to the extent.
    elif self.z > envn.z_bc_plus():
        #print('Invasion into the +z boundary.')        
        self.z -= envn.z_displace*(self.z-envn.z_bc_plus()) #The magnitude depends on the extent of invasion and the field's sensitivity to the extent.

def alter(O2, matrix_grid, Nnbl_grid, Nscl_grid, Nscl_col_grid):
    """
    Update the tumour microenvironment by considering the agents collectively.
    """
    if envn.staticO2 == 0:
        """
        If the oxygen level is not static, calculate the amount consumed by all living cells in a time step, add the amount supplied by the vasculature, and the result will be the new O2 level.
        Assuming a diffusivity of 1.75e-5 cm2 s-1 (Grote et al., 1977), the diffusion length in an hour is 0.5 cm or 5 mm, so diffusion is not limiting; it is instantaneous.
        """
        [span_xr, span_yr, span_zr] = np.shape(Nnbl_grid) 
        N_voxel = span_xr*span_yr*span_zr
        O2 += (np.sum(Nnbl_grid)+np.sum(Nscl_grid))*envn.P_O20*step_size/(N_voxel*envn.V_grid()*1e-15)/envn.Cs_O2
        O2 += P_O2v
        if O2 < 0:
            O2 = 0
        elif O2 > 1:
            O2 = 1         
    matrix_grid += Nscl_col_grid*envn.P_matrix*step_size/envn.V_grid() #Calculate the volume of matrix produced by living Schwann cells in a time step.
    return O2, matrix_grid

"""
Define the neuroblastic agent functions.
"""
def sense(self, O2, Nnbl_grid, Nscl_grid, Nnbn_grid, Nscn_grid):
    """
    Let a neuroblastoma cell respond and potentially adapt to the extracellular environment.
    """
    if self.apop == 0 and self.necro == 0:
        """
        Detect the presence of stressors in the microenvironment.
        """
        if np.random.uniform(0, 1) < (1-O2*envn.Cs_O2/(O2*envn.Cs_O2+envn.C50_necro)): #The probability is derived from the hypoxic equilibrium, so it is time-independent.
            self.hypoxia = 1
        else:
            self.hypoxia = 0       
    
        if np.random.uniform(0, 1) < (1-O2*envn.Cs_O2/(O2*envn.Cs_O2+envn.C50_necro)):
            self.nutrient = 0
        else:
            self.nutrient = 1
               
        chemo = 0
        for i in range(len(chemo_start)):
            if t>=chemo_start[i] and t<=chemo_end[i] and np.random.uniform(0, 1)<sum(chemo_effects)/len(chemo_effects):
                chemo = 1 #Drug delivery is assumed to be instantaneous.
                break

        """
        Repair any damaged DNA and induce new damages to DNA.
        1. p53 and p73 repair damaged DNA.
        2. Shortened telomeres make the chromosomes unstable.
        3. Hypoxia damages DNA.
        4. Chemotherapy damages DNA.
        """
        if self.DNA_damage == 1:
            if self.p53 == 1 or self.p73 == 1: #Effects of p53 and p73 on DNA_damage are considered before the proteins are updated.
                self.DNA_damage = 0
        else:
            if np.random.uniform(0, 1) < (1-self.telo_count/envn.telo_critical)*step_size:
                self.DNA_damage = 1
            elif (np.random.uniform(0, 1) < envn.P_DNA_damageHypo*step_size and self.hypoxia == 1):
                self.DNA_damage = 1
            elif chemo == 1 and 1 < self.cycle < 2 and np.random.uniform(0, 1) < envn.P_apopChemo*step_size:
                self.DNA_damage = 1        

        """
        Update necrotic signals.
        1. Signals from necrotic cells.
        2. Glycolysis produces lactic acid even if it is successful.
        3. Lack of nutrients other than oxygen.
        """
        [span_xr, span_yr, span_zr] = np.shape(Nnbl_grid)
        span_x = span_xr*envn.R_voxel*2
        span_y = span_yr*envn.R_voxel*2
        span_z = span_zr*envn.R_voxel*2
        x_unit = math.floor((self.x + span_x/2)/envn.R_voxel/2)
        y_unit = math.floor((self.y + span_y/2)/envn.R_voxel/2)
        z_unit = math.floor((self.z + span_z/2)/envn.R_voxel/2)
        
        dummy_Nn = Nnbn_grid[x_unit, y_unit, z_unit] + Nscn_grid[x_unit, y_unit, z_unit]
        if (x_unit-1) >= 0:
            dummy_Nn += Nnbn_grid[x_unit-1, y_unit, z_unit] + Nscn_grid[x_unit-1, y_unit, z_unit]
        if (x_unit+1) < span_xr:
            dummy_Nn += Nnbn_grid[x_unit+1, y_unit, z_unit] + Nscn_grid[x_unit+1, y_unit, z_unit]
        if (y_unit-1) >= 0:
            dummy_Nn += Nnbn_grid[x_unit, y_unit-1, z_unit] + Nscn_grid[x_unit, y_unit-1, z_unit]
        if (y_unit+1) < span_yr:
            dummy_Nn += Nnbn_grid[x_unit, y_unit+1, z_unit] + Nscn_grid[x_unit, y_unit+1, z_unit]
        if (z_unit-1) >= 0:
            dummy_Nn += Nnbn_grid[x_unit, y_unit, z_unit-1] + Nscn_grid[x_unit, y_unit, z_unit-1]
        if (z_unit+1) < span_zr:
            dummy_Nn += Nnbn_grid[x_unit, y_unit, z_unit+1] + Nscn_grid[x_unit, y_unit, z_unit+1]
        
        stress = 0
        
        for j in range(0, int(dummy_Nn)):
            if np.random.uniform(0, 1) < envn.P_necroIS*step_size:
                self.necro_signal += 1*step_size
                stress = 1
  
        if self.nutrient == 1:
            if self.hypoxia == 0:
                self.ATP = 1
                if self.necro_signal > 0 and (np.random.uniform(0, 1) < envn.P_necrorp*step_size) and stress == 0:
                    self.necro_signal -= 1*step_size
            elif (np.random.uniform(0, 1) < envn.glycoEff): #The contribution of glycolysis to necrosis is time-independent.
                self.ATP = 1
                self.necro_signal += 1*step_size
            else:
                self.ATP = 0
                self.necro_signal += 2*step_size
        else:
            self.ATP = 0
            self.necro_signal += 1*step_size
        
        """
        Update most of the intracellular species.
        Note that the effects of p53 and p73 on HIF are considered before p53 and p73 are updated.
        Chemotherapy inhibits CHK1, JAB1, HIF, MYCN, and p53. It is assumed that drug delivery is instantaneous.
        """
        if (np.random.uniform(0, 1) < self.MYCN_fn):
            self.MYCN = 1
        else:
            self.MYCN = 0
        if self.MYCN == 1:
            for i in range(len(chemo_start)):
                if t>=chemo_start[i] and t<=chemo_end[i] and np.random.uniform(0, 1) < chemo_effects[3]:
                    self.MYCN = 0
                    break

        if (np.random.uniform(0, 1) < self.MAPK_RAS_fn):
            self.MAPK_RAS = 1
        else:
            self.MAPK_RAS = 0
            
        if (np.random.uniform(0, 1) < self.JAB1_fn):
            self.JAB1 = 1
        else:
            self.JAB1 = 0
        if self.JAB1 == 1:
            for i in range(len(chemo_start)):
                if t>=chemo_start[i] and t<=chemo_end[i] and np.random.uniform(0, 1) < chemo_effects[1]:
                    self.JAB1 = 0
                    break
                
        if (np.random.uniform(0, 1) < self.CHK1_fn and self.DNA_damage == 1):
            self.CHK1 = 1
        elif (np.random.uniform(0, 1) < self.CHK1_fn and self.DNA_damage == 1 and self.MYCN == 1):
            self.CHK1 = 1
        else:
            self.CHK1 = 0
        if self.CHK1 == 1:
            for i in range(len(chemo_start)):
                if t>=chemo_start[i] and t<=chemo_end[i] and np.random.uniform(0, 1) < chemo_effects[0]:
                    self.CHK1 = 0
                    break     

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
        if self.HIF == 1:
            if self.p53 == 1:
                self.HIF = 0
            elif self.p73 == 1:
                self.HIF = 0   
        if self.HIF == 1:
            for i in range(len(chemo_start)):
                if t>=chemo_start[i] and t<=chemo_end[i] and np.random.uniform(0, 1) < chemo_effects[2]:
                    self.HIF = 0
                    break 

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
            if (self.MYCN == 1 and self.MYCN_amp == 1):
                self.p53 = 0
        if self.p53 == 1:
            for i in range(len(chemo_start)):
                if t>=chemo_start[i] and t<=chemo_end[i] and np.random.uniform(0, 1) < chemo_effects[5]:
                    self.p53 = 0
                    break         
        
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
    
        """
        Try to repair any unreplicated DNA and consider the effects of unreplicated DNA on the intracellular species.
        """
        if self.DNA_unreplicated == 1:
            if self.p53 == 1 or self.p73 == 1:
                self.DNA_unreplicated = 0

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

        """
        Differentiation due to stimulation from Schwann cells.
        """
        dummy_Nnbl = Nnbl_grid[x_unit, y_unit, z_unit]
        dummy_Nscl = Nscl_grid[x_unit, y_unit, z_unit]
        if (x_unit-1) >= 0:
            dummy_Nnbl += Nnbl_grid[x_unit-1, y_unit, z_unit]
            dummy_Nscl += Nscl_grid[x_unit-1, y_unit, z_unit]
        if (x_unit+1) < span_xr:
            dummy_Nnbl += Nnbl_grid[x_unit+1, y_unit, z_unit]
            dummy_Nscl += Nscl_grid[x_unit+1, y_unit, z_unit]
        if (y_unit-1) >= 0:
            dummy_Nnbl += Nnbl_grid[x_unit, y_unit-1, z_unit]
            dummy_Nscl += Nscl_grid[x_unit, y_unit-1, z_unit]
        if (y_unit+1) < span_yr:
            dummy_Nnbl += Nnbl_grid[x_unit, y_unit+1, z_unit]
            dummy_Nscl += Nscl_grid[x_unit, y_unit+1, z_unit]
        if (z_unit-1) >= 0:
            dummy_Nnbl += Nnbl_grid[x_unit, y_unit, z_unit-1]
            dummy_Nscl += Nscl_grid[x_unit, y_unit, z_unit-1]
        if (z_unit+1) < span_zr:
            dummy_Nnbl += Nnbl_grid[x_unit, y_unit, z_unit+1]
            dummy_Nscl += Nscl_grid[x_unit, y_unit, z_unit+1]
        
        if np.random.uniform(0, 1) < step_size*envn.nbdiff_jux*dummy_Nscl/(dummy_Nnbl+dummy_Nscl):
            self.degdiff += envn.nbdiff_amount*step_size
            if self.degdiff > 1:
                self.degdiff = 1
        elif np.random.uniform(0, 1) < step_size*envn.nbdiff_para*np.sum(Nscl_grid)/(np.sum(Nnbl_grid)+np.sum(Nscl_grid)):
            self.degdiff += envn.nbdiff_amount*step_size
            if self.degdiff > 1:
                self.degdiff = 1
        self.cycdiff = 1 -self.degdiff                

        """
        Update apoptotic signals.
        1. CAS triggers apoptosis.
        2. CAS-independent pathways linking damaged DNA to apoptosis.
        3. Schwann cells induce apoptosis.
        """        
        stress = 0
        
        if self.CAS == 1:
            self.apop_signal += 1*step_size
            stress = 1
        elif self.DNA_damage == 1 and self.ATP == 1 and np.random.uniform(0, 1) < step_size*envn.P_DNA_damage_pathways:
            self.apop_signal += 1*step_size
            stress = 1            

        if np.random.uniform(0, 1) < step_size*envn.nbapop_jux*dummy_Nscl/(dummy_Nnbl+dummy_Nscl):
            self.apop_signal += 1*step_size
            stress = 1
        elif np.random.uniform(0, 1) < step_size*envn.nbapop_para*np.sum(Nscl_grid)/(np.sum(Nnbl_grid)+np.sum(Nscl_grid)):
            self.apop_signal += 1*step_size
            stress = 1

        if self.apop_signal > 0 and (np.random.uniform(0, 1) < envn.P_apoprp*step_size) and stress == 0:
            self.apop_signal -= 1*step_size

        """
        Check if a cell is living, apoptotic, or necrotic.
        """
        if self.apop_signal >= envn.apop_critical:
            self.telo = 0
            self.ALT = 0
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
            self.mobile = 0
            self.ATP = 0
            self.apop = 1
        elif self.necro_signal >= self.necro_critical:
            self.telo = 0
            self.ALT = 0
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
            self.mobile = 0
            self.ATP = 0
            self.necro = 1
    elif self.apop == 1 and self.necro == 0:
        if np.random.uniform(0, 1) < envn.P_2ndnecro*step_size: #Secondary necrosis.
            self.apop = 0
            self.necro = 1
    
def cell_cycle(self, cycle_stages):
    """
    Progress through the neuroblastoma cell cycle.    
    In the cell cycle, 0 = G0, 1 = G1/S, 2 = S/G2, 3 = G2/M, 4 = division.
    """
    if np.random.uniform(0, 1) < envn.P_cycle_nb and self.neighbours <= envn.N_neighbours and self.ATP == 1 and self.apop == 0 and self.necro == 0:
        """
        Regardless of the cell cycle stage, the neuroblast is subjected to several regulatory mechanisms.
        1. Basal cycling rate.
        2. Contact inhibition.
        3. ATP availability.
        4. Not being apoptotic/necrotic.
        """
        if self.cycle < 1:
            if self.cycle == 0:
                if (np.random.uniform(0, 1) < self.cycdiff):
                    """
                    Fully differentiated cells cannot exit the resting phase.
                    """                    
                    if ((self.MAPK_RAS == 1 or self.MYCN == 1) and self.p21 == 0 and self.p27 == 0) or (self.ID2 == 1):
                        self.cycle += step_size/cycle_stages[0]
            elif ((self.MAPK_RAS == 1 or self.MYCN == 1) and self.p21 == 0 and self.p27 == 0) or (self.ID2 == 1):
                self.cycle += step_size/cycle_stages[0]
                if self.cycle >= 1 and ((self.MAPK_RAS == 1 and self.p21 == 0 and self.p27 == 0) or (self.ID2 == 1)) == 0:
                    self.cycle -= step_size/cycle_stages[0]
        elif self.cycle < 2:
            if self.p21 == 0 and self.p27 == 0:
                self.cycle += step_size/cycle_stages[1]
                if self.DNA_unreplicated == 0:
                    """
                    The cell's DNA replicates in the S phase.
                    """  
                    if np.random.uniform(0, 1) < envn.P_unrepDNA*step_size:
                        self.DNA_unreplicated = 1
                    elif (np.random.uniform(0, 1) < envn.P_unrepDNAHypo*step_size and self.hypoxia == 1):
                        self.DNA_unreplicated = 1
        elif self.cycle < 3:
            self.cycle += step_size/cycle_stages[2]
            if self.cycle >= 3 and self.CDC25C == 0:  
                self.cycle -= step_size/cycle_stages[2]
        elif self.cycle < 4:
            self.cycle += step_size/cycle_stages[3]

def divide(self):
    """
    Telomere repair and division of a living neuroblast.
    """
    if self.apop == 0 and self.necro == 0:
        if self.telo_count < envn.telo_maximum:
            """
            If the number of telomere units is below the maximum value, make an attempt to repair it.
            """            
            if self.telo == 1 and (np.random.uniform(0, 1) < envn.P_telorp*step_size):
                telo_dummy = 1
                for i in range(len(chemo_start)):
                    """
                    Telomerase: Chemotherapy inhibits TEP1, which is important for telomerase activity.
                    It is assumed that drug delivery is instantaneous.
                    """  
                    if t>=chemo_start[i] and t<=chemo_end[i] and np.random.uniform(0, 1) < chemo_effects[4]:
                        telo_dummy = 0
                        break
                if telo_dummy == 1:
                    self.telo_count += 1
            elif self.ALT == 1 and (np.random.uniform(0, 1) < envn.P_telorp*step_size):
                """
                ALT: This mechanism acts independently of telomerase.
                """                 
                self.telo_count += 1                
        if self.cycle >= 4:
            """
            At the end of the cell cycle, the cell divides.
            """            
            self.cycle = 0 #Back to the beginning of the cell cycle.
            if self.telo_count > 0:
                self.telo_count -= 1 #If it has at least one telomere unit, shorten its telomere.
            return True
        else:
            return False
    else:
        return False

"""
Define the Schwann cells' agent functions.
"""
def sc_sense(self, O2, Nnbn_grid, Nscn_grid):
    """
    Let a Schwann cell respond and potentially adapt to the extracellular environment.
    """
    if self.apop == 0 and self.necro == 0:
        """
        Detect the presence of stressors in the microenvironment.
        The probability for necrosis due to hypoxia is given in a paper (Warren and Partridge, 2016).
        """
        if np.random.uniform(0, 1) < (1-O2*envn.Cs_O2/(O2*envn.Cs_O2+envn.C50_necro)): #The probability is derived from the hypoxic equilibrium, so it is time-independent.
            self.hypoxia = 1
        else:
            self.hypoxia = 0       
    
        if np.random.uniform(0, 1) < (1-O2*envn.Cs_O2/(O2*envn.Cs_O2+envn.C50_necro)):
            self.nutrient = 0
        else:
            self.nutrient = 1

        chemo = 0
        for i in range(len(chemo_start)):
            if t>=chemo_start[i] and t<=chemo_end[i] and np.random.uniform(0, 1)<sum(chemo_effects)/len(chemo_effects):
                chemo = 1
                break
           
        """
        Repair any damaged DNA and induce new damages to DNA.
        1. Shortened telomeres make the chromosomes unstable.
        2. Hypoxia damages DNA.
        3. Chemotherapy damages DNA.
        4. In the absence of chemotherapy, let the cell repair its DNA.
        """
        if self.DNA_damage == 0:
            if np.random.uniform(0, 1) < (1-self.telo_count/envn.telo_critical)*step_size:
                self.DNA_damage = 1
            elif (np.random.uniform(0, 1) < envn.P_DNA_damageHypo*step_size and self.hypoxia == 1):
                self.DNA_damage = 1
            elif chemo == 1 and 1 < self.cycle < 2 and np.random.uniform(0, 1) < envn.P_apopChemo*step_size:
                self.DNA_damage = 1
        elif chemo == 0 and np.random.uniform(0, 1) < envn.P_DNA_damagerp*step_size:
            self.DNA_damage = 0
                
        """
        Update necrotic signals.
        1. Signals from necrotic cells.
        2. Glycolysis produces lactic acid even if it is successful.
        3. Lack of nutrients other than oxygen.
        """
        [span_xr, span_yr, span_zr] = np.shape(Nnbn_grid)
        span_x = span_xr*envn.R_voxel*2
        span_y = span_yr*envn.R_voxel*2
        span_z = span_zr*envn.R_voxel*2
        x_unit = math.floor((self.x + span_x/2)/envn.R_voxel/2)
        y_unit = math.floor((self.y + span_y/2)/envn.R_voxel/2)
        z_unit = math.floor((self.z + span_z/2)/envn.R_voxel/2)        
        
        dummy_Nn = Nnbn_grid[x_unit, y_unit, z_unit] + Nscn_grid[x_unit, y_unit, z_unit]
        if (x_unit-1) >= 0:
            dummy_Nn += Nnbn_grid[x_unit-1, y_unit, z_unit] + Nscn_grid[x_unit-1, y_unit, z_unit]
        if (x_unit+1) < span_xr:
            dummy_Nn += Nnbn_grid[x_unit+1, y_unit, z_unit] + Nscn_grid[x_unit+1, y_unit, z_unit]
        if (y_unit-1) >= 0:
            dummy_Nn += Nnbn_grid[x_unit, y_unit-1, z_unit] + Nscn_grid[x_unit, y_unit-1, z_unit]
        if (y_unit+1) < span_yr:
            dummy_Nn += Nnbn_grid[x_unit, y_unit+1, z_unit] + Nscn_grid[x_unit, y_unit+1, z_unit]
        if (z_unit-1) >= 0:
            dummy_Nn += Nnbn_grid[x_unit, y_unit, z_unit-1] + Nscn_grid[x_unit, y_unit, z_unit-1]
        if (z_unit+1) < span_zr:
            dummy_Nn += Nnbn_grid[x_unit, y_unit, z_unit+1] + Nscn_grid[x_unit, y_unit, z_unit+1]

        stress = 0
        
        for j in range(0, int(dummy_Nn)):
            if np.random.uniform(0, 1) < envn.P_necroIS*step_size:
                self.necro_signal += 1*step_size
                stress = 1
       
        if self.nutrient == 1:
            if self.hypoxia == 0:
                self.ATP = 1
                if self.necro_signal > 0 and (np.random.uniform(0, 1) < envn.P_necrorp*step_size) and stress == 0:
                    self.necro_signal -= 1*step_size
            elif (np.random.uniform(0, 1) < envn.glycoEff): #The contribution of glycolysis to necrosis is time-independent.
                self.ATP = 1
                self.necro_signal += 1*step_size
            else:
                self.ATP = 0
                self.necro_signal += 2*step_size
        else:
            self.ATP = 0
            self.necro_signal += 1*step_size

        """
        If chemotherapy is inactive, let the cell repair its unreplicated DNA.
        """
        if self.DNA_unreplicated == 1:
            if chemo == 0 and np.random.uniform(0, 1) < envn.P_unrepDNArp*step_size:
                self.DNA_unreplicated = 0

        """
        Update apoptotic signals.
        1. CAS is activated by hypoxia. Since CAS is not considered for Schwann cells, the inputs are linked to the apoptotic signals directly.
        2. CAS is activated by DNA damages, but chemo disrupts the pathways linking these damages to CAS. Since CAS is not considered for Schwann cells, the inputs are linked to the apoptotic signals directly.
        3. Missing DNA damage response pathways.
        """
        stress = 0
        
        if self.hypoxia == 1 and self.ATP == 1:
            self.apop_signal += 1*step_size
            stress = 1       
        elif chemo == 0 and self.DNA_damage == 1 and self.ATP == 1:
            self.apop_signal += 1*step_size
            stress = 1
        elif self.DNA_damage == 1 and self.ATP == 1 and np.random.uniform(0, 1) < step_size*envn.P_DNA_damage_pathways:
            self.apop_signal += 1*step_size
            stress = 1

        if self.apop_signal > 0 and (np.random.uniform(0, 1) < envn.P_apoprp*step_size) and stress == 0:
            self.apop_signal -= 1*step_size

        """
        Check if a cell is living, apoptotic, or necrotic.
        """
        if self.apop_signal >= envn.apop_critical:
            self.mobile = 0
            self.ATP = 0
            self.apop = 1
        elif self.necro_signal >= self.necro_critical:
            self.mobile = 0
            self.ATP = 0
            self.necro = 1
    
    elif self.apop == 1 and self.necro == 0:
        if np.random.uniform(0, 1) < envn.P_2ndnecro*step_size: #Secondary necrosis.
            self.apop = 0
            self.necro = 1

def sc_cell_cycle(self, cycle_stages, Nnbl_grid, Nscl_grid):
    """
    Progress through the Schwann cell's cell cycle.    
    In the cell cycle, 0 = G0, 1 = G1/S, 2 = S/G2, 3 = G2/M, 4 = division.
    """

    """
    The Schwann cell's ability to cycle depends on stimulation from the neuroblastic population.
    Therefore, count the neuroblasts and Schwann cells in the 3D von Neumann neighbourhood.
    Then, decide whether the stimulation is sufficient. If not, the Schwann cell has to rely on its basal cycling rate.
    """
    [span_xr, span_yr, span_zr] = np.shape(Nnbl_grid)
    span_x = span_xr*envn.R_voxel*2
    span_y = span_yr*envn.R_voxel*2
    span_z = span_zr*envn.R_voxel*2
    x_unit = math.floor((self.x + span_x/2)/envn.R_voxel/2)
    y_unit = math.floor((self.y + span_y/2)/envn.R_voxel/2)
    z_unit = math.floor((self.z + span_z/2)/envn.R_voxel/2)    
    dummy_Nnbl = Nnbl_grid[x_unit, y_unit, z_unit]
    dummy_Nscl = Nscl_grid[x_unit, y_unit, z_unit]
    if (x_unit-1) >= 0:
        dummy_Nnbl += Nnbl_grid[x_unit-1, y_unit, z_unit]
        dummy_Nscl += Nscl_grid[x_unit-1, y_unit, z_unit]
    if (x_unit+1) < span_xr:
        dummy_Nnbl += Nnbl_grid[x_unit+1, y_unit, z_unit]
        dummy_Nscl += Nscl_grid[x_unit+1, y_unit, z_unit]
    if (y_unit-1) >= 0:
        dummy_Nnbl += Nnbl_grid[x_unit, y_unit-1, z_unit]
        dummy_Nscl += Nscl_grid[x_unit, y_unit-1, z_unit]
    if (y_unit+1) < span_yr:
        dummy_Nnbl += Nnbl_grid[x_unit, y_unit+1, z_unit]
        dummy_Nscl += Nscl_grid[x_unit, y_unit+1, z_unit]
    if (z_unit-1) >= 0:
        dummy_Nnbl += Nnbl_grid[x_unit, y_unit, z_unit-1]
        dummy_Nscl += Nscl_grid[x_unit, y_unit, z_unit-1]
    if (z_unit+1) < span_zr:
        dummy_Nnbl += Nnbl_grid[x_unit, y_unit, z_unit+1]
        dummy_Nscl += Nscl_grid[x_unit, y_unit, z_unit+1]  
    if np.random.uniform(0, 1) < step_size*envn.scpro_jux*dummy_Nnbl/(dummy_Nnbl+dummy_Nscl):
        """
        Juxtacrine signallings.
        """
        dummy_scpro = 1
    elif np.random.uniform(0, 1) < step_size*envn.scpro_para*np.sum(Nnbl_grid)/(np.sum(Nnbl_grid)+np.sum(Nscl_grid)):
        """
        Paracrine signallings.
        """
        dummy_scpro = 1
    else:
        dummy_scpro = 0
    if dummy_scpro == 1 or np.random.uniform(0, 1) < envn.P_cycle_sc:
        """
        Basal cycling rate provides a last chance.
        """
        dummy_scycle = 1
    else:
        dummy_scycle = 0

    if dummy_scycle == 1 and self.neighbours <= envn.N_neighbours and self.ATP == 1 and self.apop == 0 and self.necro == 0:
        """
        Regardless of the cell cycle stage, the Schwann cell is subjected to several regulatory mechanisms.
        1. Basal cycling rate and stimulation from Schwann cells.
        2. Contact inhibition.
        3. ATP availability.
        4. Not being apoptotic/necrotic.
        """        
        if self.cycle < 1:
            if self.cycle == 0:
                if self.DNA_damage == 0 and self.hypoxia == 0:
                    """
                    DNA damage and hypoxia must both be off during G1 and S because each can switch on p21 and p27 to arrest cycling.
                    """
                    self.cycle += step_size/cycle_stages[0]
            elif self.DNA_damage == 0 and self.hypoxia == 0:
                """
                DNA damage and hypoxia must both be off during G1 and S because each can switch on p21 and p27 to arrest cycling.
                """                
                self.cycle += step_size/cycle_stages[0]
        elif self.cycle < 2:
            if self.DNA_damage == 0 and self.hypoxia == 0:
                """
                DNA damage and hypoxia must both be off during G1 and S because each can switch on p21 and p27 to arrest cycling.
                """
                self.cycle += step_size/cycle_stages[1]
                if self.DNA_unreplicated == 0:
                    """
                    The cell's DNA replicates in the S phase.
                    """  
                    if np.random.uniform(0, 1) < envn.P_unrepDNA*step_size:
                        self.DNA_unreplicated = 1
                    elif (np.random.uniform(0, 1) < envn.P_unrepDNAHypo*step_size and self.hypoxia == 1):
                        self.DNA_unreplicated = 1
        elif self.cycle < 3:
            self.cycle += step_size/cycle_stages[2]
            if self.cycle >= 3 and (self.DNA_damage == 1 or self.DNA_unreplicated == 1):  
                """
                Either DNA damage or unreplicated DNA can switch off CDC25C to arrest G2/M transition.
                """
                self.cycle -= step_size/cycle_stages[2]
        elif self.cycle < 4:
            self.cycle += step_size/cycle_stages[3]

def sc_divide(self):
    """
    Division of a living Schwann cell.
    """
    if self.apop == 0 and self.necro == 0:
        if self.cycle >= 4: #End of the cell cycle.
            self.cycle = 0 #Back to the beginning of the cell cycle.
            if self.telo_count > 0:
                self.telo_count -= 1 #If it has at least one telomere unit, shorten its telomere.
            return True
        else:
            return False
    else:
        return False
                        
"""
Parameterisation.
1. Define a class instance which contains the default parameters.
2. Define the attributes to be initialised in each neuroblast and Schwann cell.
3. Specify the number of updates (steps) and the physiological time represented by each update in hours (step_size). 
4. Schedule chemotherapy.
5. Handle Dynamic environment args.
6. Configure random seed
"""
envn = Environment()
args = getParser().parse_args();
#steps = 3024
steps = 336
step_size = 1
chemo_start = envn.chemo_start
chemo_end = envn.chemo_end
chemo_effects = envn.chemo_effects
for key, val in envn.__dict__.items():
    if isinstance(val, list):
        setattr(envn, key, getattr(args, key, getattr(envn, key)));
    else:
        setattr(envn, key, getattr(args, key, [getattr(envn, key)])[0]);
if args.seed[0]!=0:
    np.random.seed(args.seed[0]);
    rn.seed(args.seed[0]);

"""
Initialisation.
1. Start timing the simulation.
2. Create the initial populations of neuroblasts and Schwann cells.
3. Create the initial microenvironment (oxygen level and matrix volume fraction).
4. Optimise the initial spatial distribution of all cells.
5. Locate the cells after optimising their spatial distribution.
6. Expand the grids modelling the microenvironment.
7. Create the initial vasculature.
8. Record the initial conditions.
"""
t0 = time.time()
nblist = [Neuroblastoma(envn.R_tumour(), envn.histology, envn.gradiff, args.x[0], args.y[0], args.z[0], args.MYCN_amp[0], args.TERT_rarngm[0], args.ATRX_inact[0], args.ALT[0], args.ALK[0], args.MYCN_fn00[0], args.MYCN_fn10[0], args.MYCN_fn01[0], args.MYCN_fn11[0], args.MAPK_RAS_fn00[0], args.MAPK_RAS_fn10[0], args.MAPK_RAS_fn01[0], args.MAPK_RAS_fn11[0], args.p53_fn[0], args.CHK1_fn[0], args.p21_fn[0], args.p27_fn[0], args.p73_fn[0], args.CDC25C_fn[0], args.CDS1_fn[0], args.ID2_fn[0], args.IAP2_fn[0], args.HIF_fn[0], args.BNIP3_fn[0], args.JAB1_fn[0], args.Bcl2_Bclxl_fn[0], args.BAK_BAX_fn[0], args.CAS_fn[0], args.VEGF_fn[0], args.cycle[0], args.apop[0], args.apop_signal[0], args.necro[0], args.necro_signal[0], args.telo_count[0]) for i in range(int(envn.N_cell()))]
sclist = [Schwann(envn.R_tumour(), args.cycle_sc[0], args.apop_sc[0], args.apop_signal_sc[0], args.necro_sc[0], args.necro_signal_sc[0], args.telo_count_sc[0]) for i in range(int(envn.N_scell()))]
O2 = envn.O2
matrix_grid = envn.matrix_grid()
fresolve(nblist, sclist, matrix_grid)
[x, y, z, xs, ys, zs, xa, ya, za, xsa, ysa, zsa, xn, yn, zn, xsn, ysn, zsn, xl, yl, zl, xsl, ysl, zsl, xsl_col, ysl_col, zsl_col, N_vf, V, Vs] = locate(nblist, sclist)
[matrix_grid, N_grid, Nnb_grid, Nnba_grid, Nnbn_grid, Nnbl_grid, Nsc_grid, Nsca_grid, Nscn_grid, Nscl_grid, Nscl_col_grid] = CAexpand(matrix_grid, x, y, z, xs, ys, zs, xa, ya, za, xsa, ysa, zsa, xn, yn, zn, xsn, ysn, zsn, xl, yl, zl, xsl, ysl, zsl, xsl_col, ysl_col, zsl_col, V, Vs)
[P_O2v, ang_signal, t] = [0, 0, 0]
[P_O2v, ang_signal] = vasculature(P_O2v, ang_signal, matrix_grid, Nnbl_grid, Nscl_grid, N_vf, t)
data = []
data.append(record(t, x, y, z, xs, ys, zs, xa, ya, za, xsa, ysa, zsa, xn, yn, zn, xsn, ysn, zsn, xl, yl, zl, xsl, ysl, zsl, O2, matrix_grid, N_grid, Nnb_grid, Nnba_grid, Nnbn_grid, Nnbl_grid, Nsc_grid, Nsca_grid, Nscn_grid, Nscl_grid, ang_signal))

"""
Simulate neuroblastoma progression by iterating these steps until the end of the simulation.
1. Evaluate each neuroblast.
    (a) Let it sense the microenvironment.
    (b) Try to move it forward in the cell cycle.
    (c) Try to remove it from the system.
    (d) Try to replicate it and perturb the daughter cell's position.
2. Evaluate each Schwann cell.
    (a) Let it sense the microenvironment.
    (b) Try to move it forward in the cell cycle.
    (c) Try to remove it from the system.
    (d) Try to replicate it and perturb the daughter cell's position.
3. Update the spatial distribution of cells using the mechanical model and then update the grids.
4. Update the vasculature.
5. Update the microenvironment.
6. Record the conditions at the end of the time step.
"""
for t in range(steps):
    nblist_new = []
    for i in nblist:
        sense(i, O2, Nnbl_grid, Nscl_grid, Nnbn_grid, Nscn_grid)
        cell_cycle(i, envn.cycle_stages)
        if remove(i):
            pass
        elif divide(i):
            nblist_new.append(i)
            obj_dummy = cp.deepcopy(i)            
            drift(obj_dummy, envn.R_cell)
            nblist_new.append(obj_dummy)
        else:
            nblist_new.append(i)
    del nblist
    nblist = nblist_new
    sclist_new = []
    for i in sclist:
        sc_sense(i, O2, Nnbn_grid, Nscn_grid)
        sc_cell_cycle(i, envn.cycle_stages, Nnbl_grid, Nscl_grid)
        if remove(i):
            pass
        elif sc_divide(i):
            sclist_new.append(i)
            obj_dummy = cp.deepcopy(i)            
            drift(obj_dummy, envn.R_cell)
            sclist_new.append(obj_dummy)            
        else:
            sclist_new.append(i)
    del sclist
    sclist = sclist_new
    [x, y, z, xs, ys, zs, xa, ya, za, xsa, ysa, zsa, xn, yn, zn, xsn, ysn, zsn, xl, yl, zl, xsl, ysl, zsl, xsl_col, ysl_col, zsl_col, N_vf, V, Vs] = locate(nblist, sclist)
    [matrix_grid, N_grid, Nnb_grid, Nnba_grid, Nnbn_grid, Nnbl_grid, Nsc_grid, Nsca_grid, Nscn_grid, Nscl_grid, Nscl_col_grid] = CAexpand(matrix_grid, x, y, z, xs, ys, zs, xa, ya, za, xsa, ysa, zsa, xn, yn, zn, xsn, ysn, zsn, xl, yl, zl, xsl, ysl, zsl, xsl_col, ysl_col, zsl_col, V, Vs)
    if len(xl) == 0:
        print('Tumour regressed, time: {} hours.'.format(t+1))
        break
    else: 
        fresolve(nblist, sclist, matrix_grid)
        [x, y, z, xs, ys, zs, xa, ya, za, xsa, ysa, zsa, xn, yn, zn, xsn, ysn, zsn, xl, yl, zl, xsl, ysl, zsl, xsl_col, ysl_col, zsl_col, N_vf, V, Vs] = locate(nblist, sclist)
        [matrix_grid, N_grid, Nnb_grid, Nnba_grid, Nnbn_grid, Nnbl_grid, Nsc_grid, Nsca_grid, Nscn_grid, Nscl_grid, Nscl_col_grid] = CAexpand(matrix_grid, x, y, z, xs, ys, zs, xa, ya, za, xsa, ysa, zsa, xn, yn, zn, xsn, ysn, zsn, xl, yl, zl, xsl, ysl, zsl, xsl_col, ysl_col, zsl_col, V, Vs)
        [P_O2v, ang_signal] = vasculature(P_O2v, ang_signal, matrix_grid, Nnbl_grid, Nscl_grid, N_vf, t+1)
        [O2, matrix_grid] = alter(O2, matrix_grid, Nnbl_grid, Nscl_grid, Nscl_col_grid)
        data.append(record(t+1, x, y, z, xs, ys, zs, xa, ya, za, xsa, ysa, zsa, xn, yn, zn, xsn, ysn, zsn, xl, yl, zl, xsl, ysl, zsl, O2, matrix_grid, N_grid, Nnb_grid, Nnba_grid, Nnbn_grid, Nnbl_grid, Nsc_grid, Nsca_grid, Nscn_grid, Nscl_grid, ang_signal))

"""
Exit routine.
1. Stop timing the simulation.
2. Save the data.
"""
t1 = time.time()
simulation_t = t1-t0
print('Simulation time: {}.'.format(simulation_t))
np.save('data.npy', data, allow_pickle=True)