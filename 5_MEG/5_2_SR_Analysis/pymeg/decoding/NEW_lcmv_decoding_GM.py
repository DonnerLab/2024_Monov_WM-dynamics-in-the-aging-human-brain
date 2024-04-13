
# Some imports:
import logging
import mne
import numpy as np
import os
import scipy.io as sio
import time
os.environ["PYMEG_CACHE_DIR"] = "/mnt/homes/home024/pmurphy/tmp"
from joblib import Memory # Provides caching of results

from os import makedirs
from os.path import join
from glob import glob

from pymeg import lcmv as pymeglcmv
from pymeg import source_reconstruction as pymegsr
from pymeg import decoding_GM

from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.decomposition import PCA
from sklearn import linear_model

import datetime
import pandas as pd

from pymeg.lcmv_mci import get_stim_epoch
from pymeg.source_reconstruction import (
    get_leadfield,
    make_trans,   # NB: source_reconstruction also has get_trans_epoch, get_ctf_trans - not sure which we want, NW has get_trans
)

# Setup some paths:
memory = Memory(cachedir='/mnt/homes/home024/pmurphy/tmp/')
subjects_dir = "/home/pmurphy/meg_data/MCI/MRIs/fs_converted" # freesurfer subject dirs
trans_dir = "/home/pmurphy/meg_data/MCI/MRIs/trans_mats" # transofrmation matrices

# We need to restrict the number of threads that we use for the cluster
def set_n_threads(n):
    os.environ["OPENBLAS_NUM_THREADS"] = str(n)
    os.environ["MKL_NUM_THREADS"] = str(n)
    os.environ["OMP_NUM_THREADS"] = str(n)

# make dict of subjects/sessions/recordings
# subjects = {'01': [()],
#            '02': [()],
#            '03': [()],
#            '04': [()],
#            '05': [()],
#            '06': [()],
#            '07': [()],
#            '08': [()],
#            '09': [()],
#            '10': [()],
#            '14': [()],
#            '15': [()],
#            '16': [()],
#            '17': [()],
#            '18': [()],
#            '19': [()],
#            '20': [()],
#            '21': [()],
#            '22': [()],
#            '23': [()],
#            '25': [()],
#            '26': [()],
#            '27': [()],
#            '28': [()],
#            '29': [()],
#            '30': [()],
#            '31': [()],
#            '33': [()],
#            '35': [()],
#            '36': [()],
#            '37': [()],
#            '38': [()],
#            '39': [()],
#            '40': [()],
#            '41': [()],
#            '42': [()],
#            '43': [()],
#            '45': [()],
#            '46': [()],
#            '47': [()],
#            '48': [()],
#            '49': [()],
#            '50': [()],
#            '51': [()],
#            '52': [()],
#            '53': [()],
#            '54': [()]}#,# Subjects for which the oct 6 spacing delivers better results
# subjects = {'11': [()],
#             '24': [()],
#             '32': [()]}#,,# Subjects for which the oct 7 spacing delivers better results

# subjects = {'04': [()]}#,,# Test one example

subjects = {'18': [()],
            '51': [()],
            '25': [()],
            '46': [()],
            '19': [()],
            '43': [()],
            '45': [()]}#,,# Subjects for which the fsaverage head is now used due to failure in the SR


mem_demand = {'01': 3, '02': 2, '03': 2, '04': 2, '05': 2, '06': 2,   # SETTINGS FOR VERTEX-AVERAGED DECODING
              '07': 2, '08': 2, '09': 2, '10': 2, '11': 4, '14': 2,
              '15': 2, '16': 3, '17': 2, '18': 2, '19': 3,
              '20': 2, '21': 2, '22': 2, '23': 2, '24': 2, '25': 2,
              '26': 2, '27': 2, '28': 2, '29': 2, '30': 2, '31': 2,
              '32': 10, '33': 2, '35': 2, '36': 2, '37': 2, '38': 2,
              '39': 2, '40': 2, '41': 2, '42': 2, '43': 2, '45': 2,
              '46': 2, '47': 2, '48': 2, '49': 2, '50': 2, '51': 2,
              '52': 2, '53': 2, '54': 2, 'P03': 2}

# typical processing demands (in # vertices) per area (not used but useful reference fo setting ntasks at submit)
mem_area = {'vfcPrimary': 3650, 'vfcEarly': 9610, 'vfcV3ab': 3280,
            'vfcIPS01': 4200, 'vfcIPS23': 2200, 'JWG_aIPS': 870,
            'JWG_IPS_PCeS': 4600, 'JWG_M1': 2900, 'HCPMMP1_premotor': 9900}



#Specify settings on which to run the anaylsis
mode = 'super_roi' # choose between all_parcels, coarse, or super_roi or complete_glasser or whole_cortex
min_freq = 5 #choose lower bound for included frequency spectrum
max_freq = 35 # choose upper bound for included frequency
delay_duration = 1 # select delay duration (1, 3 or 9 s)
#!!!!!!! PC settings for perforing the decoding have to be changed manually!!!!

def a2lab_all_source_areas():
   labels_new_glasser = [
           "V1", "V2", "V3", "V4", "V6", "V3A", "V7", "IPS1", "V3B", "V6A",
           "MST", "LO1", "LO2", "MT","V4t", "FST", "LO3", "V3CD", "PH", "V8",
           "FFC", "PIT", "VMV1", "VMV3", "VMV2", "VVC", "7Pm", "7AL", "7Am", "7PL",
           "7PC", "LIPv", "VIP", "MIP", "LIPd", "AIP", "6a", "6d",  "FEF",
           "6v", "6r", "EC", "PreS", "H", "PeEc", "PHA1", "PHA3", "PEF",
           "TF", "PHA2", "PCV", "7m", "POS1", "v23ab", "d23ab", "31pv",  "31pd", "31a",
           "RSC", "POS2", "ProS", "DVT", "23d", "23c", "SFL", "8Av", "8Ad", "8BL",
           "9p", "8C", "p9-46v", "46", "a9-46v", "9-46d", "9a", "i6-8", "s6-8",
       ]
   glasser_new_clusters = {
           l:['L_{}_ROI-lh'.format(l), 'R_{}_ROI-rh'.format(l)] for l in labels_new_glasser
       }

   jwg_M1 = {
           "JWG_lr_M1": ['lh.JWDG.lr_M1-lh', 'rh.JWDG.lr_M1-rh',
           ],
       }
   #concatenate rois of interest

   glasser_new_clusters.update(jwg_M1)

   areas_to_labels = glasser_new_clusters
   return areas_to_labels

def a2lab_all_glasser_parcels():
        labels = [
            "V1",	"MST",	"V6",	"V2",	"V3",	"V4",	"V8",	"4"	, "3b",	"FEF",
           	"PEF", "55b",	"V3A",	"RSC",	"POS2",	"V7",	"IPS1",	"FFC",	"V3B",	"LO1",
           	"LO2",	"PIT",	"MT",	"A1",	"PSL",	"SFL",	"PCV",	"STV",	"7Pm",	"7m",
            "POS1",	"23d",	"v23ab",	"d23ab",	"31pv",	"5m",	"5mv",	"23c",	"5L",	"24dd",
           	"24dv",	"7AL",	"SCEF",	"6ma",	"7Am",	"7PL",	"7PC",	"LIPv",	"VIP",	"MIP",
           	"1",	"2",	"3a",	"6d",	"6mp",	"6v",	"p24pr",	"33pr",	"a24pr",	"p32pr",
           	"a24",	"d32",	"8BM",	"p32",	"10r",	"47m",	"8Av",	"8Ad",	"9m",	"8BL",
           	"9p",	"10d",	"8C",	"44",	"45",	"47l",	"a47r",	"6r",	"IFJa",	"IFJp",
           	"IFSp",	"IFSa",	"p9-46v",	"46",	"a9-46v",	"9-46d",	"9a",	"10v",	"a10p",
           	"10pp",	"11l",	"13l",	"OFC",	"47s",	"LIPd",	"6a",	"i6-8",	"s6-8",	"43",
           	"OP4",	"OP1",	"OP2-3",	"52",	"RI",	"PFcm",	"PoI2",	"TA2",	"FOP4",	"MI",	"Pir",
           	"AVI",	"AAIC",	"FOP1",	"FOP3",	"FOP2",	"PFt",	"AIP",	"EC",	"PreS",	"H", "ProS",
           	"PeEc",	"STGa",	"PBelt",	"A5",	"PHA1",	"PHA3",	"STSda",	"STSdp",	"STSvp",	"TGd",
           	"TE1a",	"TE1p",	"TE2a",	"TF",	"TE2p",	"PHT",	"PH",	"TPOJ1",	"TPOJ2",	"TPOJ3",
           	"DVT",	"PGp",	"IP2",	"IP1",	"IP0",	"PFop",	"PF",	"PFm",	"PGi",	"PGs",
           	"V6A",	"VMV1",	"VMV3",	"PHA2",	"V4t",	"FST",	"V3CD",	"LO3",	"VMV2",	"31pd",	"31a",
           	"VVC",	"25",	"s32",	"pOFC",	"PoI1",	"Ig",	"FOP5",	"p10p",	"p47r",	"TGv",
           	"MBelt",	"LBelt",	"A4",	"STSva",	"TE1m",	"PI",	"a32pr",	"p24",
        ]
        # fmt: on
        areas_to_labels = {l:['L_{}_ROI-lh'.format(l), 'R_{}_ROI-rh'.format(l)] for l in labels}
        return areas_to_labels


def a2lab_coarse(): # define areas as vertices of all rois in the area
    areas_to_labels = {
      "V1": [
           "L_V1_ROI-lh", "R_V1_ROI-rh",
      ],
      "V2-V4": [
         "L_V2_ROI-lh", "R_V2_ROI-rh",
         "L_V3_ROI-lh", "R_V3_ROI-rh",
         "L_V4_ROI-lh", "R_V4_ROI-rh",
      ],
      "Dorsal_visual": [
         "L_V6_ROI-lh", "R_V6_ROI-rh",
         "L_V3A_ROI-lh", "R_V3A_ROI-rh",
         "L_V7_ROI-lh", "R_V7_ROI-rh",
         "L_IPS1_ROI-lh", "R_IPS1_ROI-rh",
         "L_V3B_ROI-lh", "R_V3B_ROI-rh",
         "L_V6A_ROI-lh", "R_V6A_ROI-rh",
      ],
      "MT": [
         "L_MST_ROI-lh", "R_MST_ROI-rh",
         "L_LO1_ROI-lh", "R_LO1_ROI-rh",
         "L_LO2_ROI-lh", "R_LO2_ROI-rh",
         "L_MT_ROI-lh", "R_MT_ROI-rh",
         "L_V4t_ROI-lh", "R_V4t_ROI-rh",
         "L_FST_ROI-lh", "R_FST_ROI-rh",
         "L_LO3_ROI-lh", "R_LO3_ROI-rh",
         "L_V3CD_ROI-lh", "R_V3CD_ROI-rh",
         "L_PH_ROI-lh", "R_PH_ROI-rh",
      ],
      "Ventral_visual": [
         "L_V8_ROI-lh", "R_V8_ROI-rh",
         "L_FFC_ROI-lh", "R_FFC_ROI-rh",
         "L_PIT_ROI-lh", "R_PIT_ROI-rh",
         "L_VMV1_ROI-lh", "R_VMV1_ROI-rh",
         "L_VMV2_ROI-lh", "R_VMV2_ROI-rh",
         "L_VMV3_ROI-lh", "R_VMV3_ROI-rh",
         "L_VVC_ROI-lh", "R_VVC_ROI-rh",
      ],
      "Superior_parietal": [
         "L_7Pm_ROI-lh", "R_7Pm_ROI-rh",
         "L_7AL_ROI-lh", "R_7AL_ROI-rh",
         "L_7Am_ROI-lh", "R_7Am_ROI-rh",
         "L_7PL_ROI-lh", "R_7PL_ROI-rh",
         "L_7PC_ROI-lh", "R_7PC_ROI-rh",
         "L_LIPv_ROI-lh", "R_LIPv_ROI-rh",
         "L_VIP_ROI-lh", "R_VIP_ROI-rh",
         "L_MIP_ROI-lh", "R_MIP_ROI-rh",
         "L_LIPd_ROI-lh", "R_LIPd_ROI-rh",
         "L_AIP_ROI-lh", "R_AIP_ROI-rh",
      ],
      "PMd": [
         "L_6a_ROI-lh", "R_6a_ROI-rh",
         "L_6d_ROI-lh", "R_6d_ROI-rh",
      ],
      "Eye_fields": [
         "L_FEF_ROI-lh", "R_FEF_ROI-rh",
         "L_PEF_ROI-lh", "R_PEF_ROI-rh",
      ],
      "PMv": [
         "L_6v_ROI-lh", "R_6v_ROI-rh",
         "L_6r_ROI-lh", "R_6r_ROI-rh",
      ],
      "Medial_temporal": [
         "L_EC_ROI-lh", "R_EC_ROI-rh",
         "L_PreS_ROI-lh", "R_PreS_ROI-rh",
         "L_H_ROI-lh", "R_H_ROI-rh",
         "L_PeEc_ROI-lh", "R_PeEc_ROI-rh",
         "L_PHA1_ROI-lh", "R_PHA1_ROI-rh",
         "L_PHA3_ROI-lh", "R_PHA3_ROI-rh",
         "L_TF_ROI-lh", "R_TF_ROI-rh",
         "L_PHA2_ROI-lh", "R_PHA2_ROI-rh",
      ],
      "Posterior_cingulate": [
         "L_PCV_ROI-lh", "R_PCV_ROI-rh",
         "L_7m_ROI-lh", "R_7m_ROI-rh",
         "L_POS1_ROI-lh", "R_POS1_ROI-rh",
         "L_v23ab_ROI-lh", "R_v23ab_ROI-rh",
         "L_d23ab_ROI-lh", "R_d23ab_ROI-rh",
         "L_31pv_ROI-lh", "R_31pv_ROI-rh",
         "L_31pd_ROI-lh", "R_31pd_ROI-rh",
         "L_31a_ROI-lh", "R_31a_ROI-rh",
         "L_RCS_ROI-lh", "R_RSC_ROI-rh",
         "L_POS2_ROI-lh", "R_POS2_ROI-rh",
         "L_ProS_ROI-lh", "R_ProS_ROI-rh",
         "L_DVT_ROI-lh", "R_DVT_ROI-rh",
         "L_23d_ROI-lh", "R_23d_ROI-rh",
         "L_23c_ROI-lh", "R_23c_ROI-rh",
      ],
      "DLPFC": [
         "L_SFL_ROI-lh", "R_SFL_ROI-rh",
         "L_8Av_ROI-lh", "R_8Av_ROI-rh",
         "L_8Ad_ROI-lh", "R_8Ad_ROI-rh",
         "L_8BL_ROI-lh", "R_8BL_ROI-rh",
         "L_9p_ROI-lh", "R_9p_ROI-rh",
         "L_8C_ROI-lh", "R_8C_ROI-rh",
         "L_p9-46v_ROI-lh", "R_p9-46v_ROI-rh",
         "L_46_ROI-lh", "R_46_ROI-rh",
         "L_a9-46v_ROI-lh", "R_a9-46v_ROI-rh",
         "L_9-46d_ROI-lh", "R_9-46d_ROI-rh",
         "L_9a_ROI-lh", "R_9a_ROI-rh",
         "L_i6-8_ROI-lh", "R_i6-8_ROI-rh",
         "L_s6-8_ROI-lh", "R_s6-8_ROI-rh",
      ],
      "JWG_M1": ["lh.JWDG.lr_M1-lh", "rh.JWDG.lr_M1-rh"],
      }
    return areas_to_labels

def a2lab_super_roi(): # define areas as vertices of all rois in the area
    areas_to_labels = {

      "super_cluster": [
         "L_V2_ROI-lh", "R_V2_ROI-rh",
         "L_V3_ROI-lh", "R_V3_ROI-rh",
         "L_V4_ROI-lh", "R_V4_ROI-rh",
         "L_V6_ROI-lh", "R_V6_ROI-rh",
         "L_V3A_ROI-lh", "R_V3A_ROI-rh",
         "L_V7_ROI-lh", "R_V7_ROI-rh",
         "L_IPS1_ROI-lh", "R_IPS1_ROI-rh",
         "L_V3B_ROI-lh", "R_V3B_ROI-rh",
         "L_V6A_ROI-lh", "R_V6A_ROI-rh",
         "L_FEF_ROI-lh", "R_FEF_ROI-rh",
         "L_PEF_ROI-lh", "R_PEF_ROI-rh",
         "L_SFL_ROI-lh", "R_SFL_ROI-rh",
         "L_8Av_ROI-lh", "R_8Av_ROI-rh",
         "L_8Ad_ROI-lh", "R_8Ad_ROI-rh",
         "L_8BL_ROI-lh", "R_8BL_ROI-rh",
         "L_9p_ROI-lh", "R_9p_ROI-rh",
         "L_8C_ROI-lh", "R_8C_ROI-rh",
         "L_p9-46v_ROI-lh", "R_p9-46v_ROI-rh",
         "L_46_ROI-lh", "R_46_ROI-rh",
         "L_a9-46v_ROI-lh", "R_a9-46v_ROI-rh",
         "L_9-46d_ROI-lh", "R_9-46d_ROI-rh",
         "L_9a_ROI-lh", "R_9a_ROI-rh",
         "L_i6-8_ROI-lh", "R_i6-8_ROI-rh",
         "L_s6-8_ROI-lh", "R_s6-8_ROI-rh",
      ],
      }
    return areas_to_labels



def a2lab_whole_cortex(): # define areas as vertices of all rois in the area
    areas_to_labels = {

      "whole_cortex": [
         "L_V1_ROI-lh",	"L_MST_ROI-lh",	"L_V6_ROI-lh",	"L_V2_ROI-lh",	"L_V3_ROI-lh",	"L_V4_ROI-lh",	"L_V8_ROI-lh",	"L_4_ROI-lh"	, "L_3b_ROI-lh",	"L_FEF_ROI-lh",
               "L_PEF_ROI-lh", "L_55b_ROI-lh",	"L_V3A_ROI-lh",	"L_RSC_ROI-lh",	"L_POS2_ROI-lh",	"L_V7_ROI-lh",	"L_IPS1_ROI-lh",	"L_FFC_ROI-lh",	"L_V3B_ROI-lh",	"L_LO1_ROI-lh",
               "L_LO2_ROI-lh",	"L_PIT_ROI-lh",	"L_MT_ROI-lh",	"L_A1_ROI-lh",	"L_PSL_ROI-lh",	"L_SFL_ROI-lh",	"L_PCV_ROI-lh",	"L_STV_ROI-lh",	"L_7Pm_ROI-lh",	"L_7m_ROI-lh",
               "L_POS1_ROI-lh",	"L_23d_ROI-lh",	"L_v23ab_ROI-lh",	"L_d23ab_ROI-lh",	"L_31pv_ROI-lh",	"L_5m_ROI-lh",	"L_5mv_ROI-lh",	"L_23c_ROI-lh",	"L_5L_ROI-lh",	"L_24dd_ROI-lh",
               "L_24dv_ROI-lh",	"L_7AL_ROI-lh",	"L_SCEF_ROI-lh",	"L_6ma_ROI-lh",	"L_7Am_ROI-lh",	"L_7PL_ROI-lh",	"L_7PC_ROI-lh",	"L_LIPv_ROI-lh",	"L_VIP_ROI-lh",	"L_MIP_ROI-lh",
               "L_1_ROI-lh",	"L_2_ROI-lh",	"L_3a_ROI-lh",	"L_6d_ROI-lh",	"L_6mp_ROI-lh",	"L_6v_ROI-lh",	"L_p24pr_ROI-lh",	"L_33pr_ROI-lh",	"L_a24pr_ROI-lh",	"L_p32pr_ROI-lh",
               "L_a24_ROI-lh",	"L_d32_ROI-lh",	"L_8BM_ROI-lh",	"L_p32_ROI-lh",	"L_10r_ROI-lh",	"L_47m_ROI-lh",	"L_8Av_ROI-lh",	"L_8Ad_ROI-lh",	"L_9m_ROI-lh",	"L_8BL_ROI-lh",
               "L_9p_ROI-lh",	"L_10d_ROI-lh",	"L_8C_ROI-lh",	"L_44_ROI-lh",	"L_45_ROI-lh",	"L_47l_ROI-lh",	"L_a47r_ROI-lh",	"L_6r_ROI-lh",	"L_IFJa_ROI-lh",	"L_IFJp_ROI-lh",
               "L_IFSp_ROI-lh",	"L_IFSa_ROI-lh",	"L_p9-46v_ROI-lh",	"L_46_ROI-lh",	"L_a9-46v_ROI-lh",	"L_9-46d_ROI-lh",	"L_9a_ROI-lh",	"L_10v_ROI-lh",	"L_a10p_ROI-lh",
               "L_10pp_ROI-lh",	"L_11l_ROI-lh",	"L_13l_ROI-lh",	"L_OFC_ROI-lh",	"L_47s_ROI-lh",	"L_LIPd_ROI-lh",	"L_6a_ROI-lh",	"L_i6-8_ROI-lh",	"L_s6-8_ROI-lh",	"L_43_ROI-lh",
               "L_OP4_ROI-lh",	"L_OP1_ROI-lh",	"L_OP2-3_ROI-lh",	"L_52_ROI-lh",	"L_RI_ROI-lh",	"L_PFcm_ROI-lh",	"L_PoI2_ROI-lh",	"L_TA2_ROI-lh",	"L_FOP4_ROI-lh",	"L_MI_ROI-lh",	"L_Pir_ROI-lh",
               "L_AVI_ROI-lh",	"L_AAIC_ROI-lh",	"L_FOP1_ROI-lh",	"L_FOP3_ROI-lh",	"L_FOP2_ROI-lh",	"L_PFt_ROI-lh",	"L_AIP_ROI-lh",	"L_EC_ROI-lh",	"L_PreS_ROI-lh",	"L_H_ROI-lh", "L_ProS_ROI-lh",
               "L_PeEc_ROI-lh",	"L_STGa_ROI-lh",	"L_PBelt_ROI-lh",	"L_A5_ROI-lh",	"L_PHA1_ROI-lh",	"L_PHA3_ROI-lh",	"L_STSda_ROI-lh",	"L_STSdp_ROI-lh",	"L_STSvp_ROI-lh",	"L_TGd_ROI-lh",
               "L_TE1a_ROI-lh",	"L_TE1p_ROI-lh",	"L_TE2a_ROI-lh",	"L_TF_ROI-lh",	"L_TE2p_ROI-lh",	"L_PHT_ROI-lh",	"L_PH_ROI-lh",	"L_TPOJ1_ROI-lh",	"L_TPOJ2_ROI-lh",	"L_TPOJ3_ROI-lh",
               "L_DVT_ROI-lh",	"L_PGp_ROI-lh",	"L_IP2_ROI-lh",	"L_IP1_ROI-lh",	"L_IP0_ROI-lh",	"L_PFop_ROI-lh",	"L_PF_ROI-lh",	"L_PFm_ROI-lh",	"L_PGi_ROI-lh",	"L_PGs_ROI-lh",
               "L_V6A_ROI-lh",	"L_VMV1_ROI-lh",	"L_VMV3_ROI-lh",	"L_PHA2_ROI-lh",	"L_V4t_ROI-lh",	"L_FST_ROI-lh",	"L_V3CD_ROI-lh",	"L_LO3_ROI-lh",	"L_VMV2_ROI-lh",	"L_31pd_ROI-lh",	"L_31a_ROI-lh",
               "L_VVC_ROI-lh",	"L_25_ROI-lh",	"L_s32_ROI-lh",	"L_pOFC_ROI-lh",	"L_PoI1_ROI-lh",	"L_Ig_ROI-lh",	"L_FOP5_ROI-lh",	"L_p10p_ROI-lh",	"L_p47r_ROI-lh",	"L_TGv_ROI-lh",
               "L_MBelt_ROI-lh",	"L_LBelt_ROI-lh",	"L_A4_ROI-lh",	"L_STSva_ROI-lh",	"L_TE1m_ROI-lh",	"L_PI_ROI-lh",	"L_a32pr_ROI-lh",	"L_p24_ROI-lh",
               "R_V1_ROI-rh",	"R_MST_ROI-rh",	"R_V6_ROI-rh",	"R_V2_ROI-rh",	"R_V3_ROI-rh",	"R_V4_ROI-rh",	"R_V8_ROI-rh",	"R_4_ROI-rh"	, "R_3b_ROI-rh",	"R_FEF_ROI-rh",
               "R_PEF_ROI-rh", "R_55b_ROI-rh",	"R_V3A_ROI-rh",	"R_RSC_ROI-rh",	"R_POS2_ROI-rh",	"R_V7_ROI-rh",	"R_IPS1_ROI-rh",	"R_FFC_ROI-rh",	"R_V3B_ROI-rh",	"R_LO1_ROI-rh",
               "R_LO2_ROI-rh",	"R_PIT_ROI-rh",	"R_MT_ROI-rh",	"R_A1_ROI-rh",	"R_PSL_ROI-rh",	"R_SFL_ROI-rh",	"R_PCV_ROI-rh",	"R_STV_ROI-rh",	"R_7Pm_ROI-rh",	"R_7m_ROI-rh",
               "R_POS1_ROI-rh",	"R_23d_ROI-rh",	"R_v23ab_ROI-rh",	"R_d23ab_ROI-rh",	"R_31pv_ROI-rh",	"R_5m_ROI-rh",	"R_5mv_ROI-rh",	"R_23c_ROI-rh",	"R_5L_ROI-rh",	"R_24dd_ROI-rh",
               "R_24dv_ROI-rh",	"R_7AL_ROI-rh",	"R_SCEF_ROI-rh",	"R_6ma_ROI-rh",	"R_7Am_ROI-rh",	"R_7PL_ROI-rh",	"R_7PC_ROI-rh",	"R_LIPv_ROI-rh",	"R_VIP_ROI-rh",	"R_MIP_ROI-rh",
               "R_1_ROI-rh",	"R_2_ROI-rh",	"R_3a_ROI-rh",	"R_6d_ROI-rh",	"R_6mp_ROI-rh",	"R_6v_ROI-rh",	"R_p24pr_ROI-rh",	"R_33pr_ROI-rh",	"R_a24pr_ROI-rh",	"R_p32pr_ROI-rh",
               "R_a24_ROI-rh",	"R_d32_ROI-rh",	"R_8BM_ROI-rh",	"R_p32_ROI-rh",	"R_10r_ROI-rh",	"R_47m_ROI-rh",	"R_8Av_ROI-rh",	"R_8Ad_ROI-rh",	"R_9m_ROI-rh",	"R_8BL_ROI-rh",
               "R_9p_ROI-rh",	"R_10d_ROI-rh",	"R_8C_ROI-rh",	"R_44_ROI-rh",	"R_45_ROI-rh",	"R_47l_ROI-rh",	"R_a47r_ROI-rh",	"R_6r_ROI-rh",	"R_IFJa_ROI-rh",	"R_IFJp_ROI-rh",
               "R_IFSp_ROI-rh",	"R_IFSa_ROI-rh",	"R_p9-46v_ROI-rh",	"R_46_ROI-rh",	"R_a9-46v_ROI-rh",	"R_9-46d_ROI-rh",	"R_9a_ROI-rh",	"R_10v_ROI-rh",	"R_a10p_ROI-rh",
               "R_10pp_ROI-rh",	"R_11l_ROI-rh",	"R_13l_ROI-rh",	"R_OFC_ROI-rh",	"R_47s_ROI-rh",	"R_LIPd_ROI-rh",	"R_6a_ROI-rh",	"R_i6-8_ROI-rh",	"R_s6-8_ROI-rh",	"R_43_ROI-rh",
               "R_OP4_ROI-rh",	"R_OP1_ROI-rh",	"R_OP2-3_ROI-rh",	"R_52_ROI-rh",	"R_RI_ROI-rh",	"R_PFcm_ROI-rh",	"R_PoI2_ROI-rh",	"R_TA2_ROI-rh",	"R_FOP4_ROI-rh",	"R_MI_ROI-rh",	"R_Pir_ROI-rh",
               "R_AVI_ROI-rh",	"R_AAIC_ROI-rh",	"R_FOP1_ROI-rh",	"R_FOP3_ROI-rh",	"R_FOP2_ROI-rh",	"R_PFt_ROI-rh",	"R_AIP_ROI-rh",	"R_EC_ROI-rh",	"R_PreS_ROI-rh",	"R_H_ROI-rh", "R_ProS_ROI-rh",
               "R_PeEc_ROI-rh",	"R_STGa_ROI-rh",	"R_PBelt_ROI-rh",	"R_A5_ROI-rh",	"R_PHA1_ROI-rh",	"R_PHA3_ROI-rh",	"R_STSda_ROI-rh",	"R_STSdp_ROI-rh",	"R_STSvp_ROI-rh",	"R_TGd_ROI-rh",
               "R_TE1a_ROI-rh",	"R_TE1p_ROI-rh",	"R_TE2a_ROI-rh",	"R_TF_ROI-rh",	"R_TE2p_ROI-rh",	"R_PHT_ROI-rh",	"R_PH_ROI-rh",	"R_TPOJ1_ROI-rh",	"R_TPOJ2_ROI-rh",	"R_TPOJ3_ROI-rh",
               "R_DVT_ROI-rh",	"R_PGp_ROI-rh",	"R_IP2_ROI-rh",	"R_IP1_ROI-rh",	"R_IP0_ROI-rh",	"R_PFop_ROI-rh",	"R_PF_ROI-rh",	"R_PFm_ROI-rh",	"R_PGi_ROI-rh",	"R_PGs_ROI-rh",
               "R_V6A_ROI-rh",	"R_VMV1_ROI-rh",	"R_VMV3_ROI-rh",	"R_PHA2_ROI-rh",	"R_V4t_ROI-rh",	"R_FST_ROI-rh",	"R_V3CD_ROI-rh",	"R_LO3_ROI-rh",	"R_VMV2_ROI-rh",	"R_31pd_ROI-rh",	"R_31a_ROI-rh",
               "R_VVC_ROI-rh",	"R_25_ROI-rh",	"R_s32_ROI-rh",	"R_pOFC_ROI-rh",	"R_PoI1_ROI-rh",	"R_Ig_ROI-rh",	"R_FOP5_ROI-rh",	"R_p10p_ROI-rh",	"R_p47r_ROI-rh",	"R_TGv_ROI-rh",
               "R_MBelt_ROI-rh",	"R_LBelt_ROI-rh",	"R_A4_ROI-rh",	"R_STSva_ROI-rh",	"R_TE1m_ROI-rh",	"R_PI_ROI-rh",	"R_a32pr_ROI-rh",	"R_p24_ROI-rh"
      ],
      }
    return areas_to_labels




if mode in ('super_roi'):
   a2lab = a2lab_super_roi
elif mode in ('coarse'):
   a2lab = a2lab_coarse
elif mode in ('all_parcels'):
   a2lab = a2lab_all_source_areas
elif mode in ('complete_glasser'):
   a2lab = a2lab_all_glasser_parcels
elif mode in ('whole_cortex'):
   a2lab = a2lab_whole_cortex

def select_delay_trials2(all_trial_ids,epochs):
    import pandas as pd
    trial_no = pd.DataFrame(all_trial_ids)
    log_ind = trial_no.isin(epochs.events[:,2])
    log_ind = log_ind[log_ind[0] == True]
    selected_trials = log_ind.index

    return selected_trials


def full_length_get_stim_epoch(subject, session, recording):
    from pymeg import preprocessing as pymegprepr
    globstring = '/home/gmonov/meg_analysis/source_reconstruction/Conv2mne/%s_%i_1*fif.gz' % (
        subject, session)  # only loading the 1s file consisting of beginning of each clean trial for estimating the covariance matrix
    filenames = glob(globstring)[0]
    epochs = mne.read_epochs(filenames)
    # epochs.times = epochs.times - 1  # PM: this was somehow necessary for initial pipeline, but *NOT* for induced
    epochs = epochs.pick_channels(
        [x for x in epochs.ch_names if x.startswith('M')])
    id_time = (-0.25 <= epochs.times) & (epochs.times <= 0)
    means = epochs._data[:, :, id_time].mean(-1)
    epochs._data -= means[:, :, np.newaxis]
    min_time, max_time = epochs.times.min() + 0.35, epochs.times.max() - 0.3
    data_cov = pymeglcmv.get_cov(epochs, tmin=min_time, tmax=max_time)
    all_trial_ids = epochs.events[:,2]

    globstring = '/home/gmonov/meg_analysis/source_reconstruction/full_length_Conv2mne/%s_%i_%i*fif.gz' % (
        subject, session, recording)  # only loading the 1s file consisting of beginning of each clean trial for estimating the covariance matrix
    filenames = glob(globstring)[0]
    epochs = mne.read_epochs(filenames)
    # epochs.times = epochs.times - 1  # PM: this was somehow necessary for initial pipeline, but *NOT* for induced
    epochs = epochs.pick_channels(
        [x for x in epochs.ch_names if x.startswith('M')])

    selected_trials = select_delay_trials2(
           all_trial_ids,epochs)


    epochs._data -= means[selected_trials, :, np.newaxis]   # need to select only trials that match those for current delay duration

    return data_cov, epochs, filenames


#Create list of all areas of interest
areas_to_labels = a2lab()

# Submit to cluster. One job per ROI and subject
def submit(only_glasser=False):
    from pymeg import parallel
    from itertools import product
    from NEW_lcmv_decoding_GM import a2lab
    areas_to_labels = a2lab()
    areas_to_labels = {'V1': [()]} # run only one single area
    h=list(areas_to_labels.keys())
    for area in list(areas_to_labels.keys()):
        for subject, tasks in subjects.items():
            for rec in tasks:
                print("Submitting %s -> %s delay %i" % (subject, area, rec))
                parallel.pmap(
                    decode,
                    [(subject, area, rec)],
                    cluster = 'SLURM',
                    walltime="200:00:00",
                    memory= mem_demand[subject]*10 + 10,
                    nodes=1,
                    tasks= mem_demand[subject] + 1,
                    env="mne",
                    name="dcd_" + area + subject,
                )
                time.sleep(1)

# Submit to cluster like hame: wait until jobs are finished to submit more
def submit_like_hame(area,only_glasser=False):
    from pymeg import parallel
    from itertools import product
    jobid=[]
    rec = delay_duration #specify which delay period
    #from NEW_lcmv_decoding_GM import a2lab
    #areas_to_labels = a2lab()
    #areas_to_labels = {'V1': [()]} # run only one single area
    #for area in list(areas_to_labels.keys()):
    for subject in subjects.keys():
         print("Submitting %s -> %s delay %i" % (subject, area, rec))
         out=parallel.pmap(
             decode,
             [(subject, area, rec)],
             cluster = 'SLURM',
             walltime="200:00:00",
             memory= mem_demand[subject]*10 + 10,
             nodes=1,
             tasks= mem_demand[subject] + 1,
             env="mne",
             name="dcd_" + area + subject,
         )
         jobid.append(int(out[0][20:-1]))
         time.sleep(1)
    return jobid




# This function returns labels for one subject
def get_labels(subject, only_glasser):
    if not only_glasser:
        labels = pymegsr.get_labels(
            subject=subject,
            filters=["*wang*.label", "*JWDG*.label"],
            annotations=["HCPMMP1"],
            sdir=subjects_dir,
        )
        labels = pymegsr.labels_exclude(
            labels=labels,
            exclude_filters=[
                "wang2015atlas.IPS4",
                "wang2015atlas.IPS5",
                "wang2015atlas.SPL",
                "JWDG_lat_Unknown",
            ],
        )
       # labels = pymegsr.labels_remove_overlap(
       #     labels=labels, priority_filters=["wang", "JWDG"]
       # )
    else:
        labels = pymegsr.get_labels(
            subject=subject,
            filters=["select_nothing"],
            annotations=["HCPMMP1"],
            sdir=subjects_dir,
        )
    return labels


#@memory.cache
def decode(
    subject,
    area,
    rec,
    epoch_type="stimulus",
    only_glasser=False,
    BEM="three_layer",
    debug=False,
    target="response",
):
    # Only show warning and above:
    mne.set_log_level("WARNING")
    pymeglcmv.logging.getLogger().setLevel(logging.INFO)

    set_n_threads(1)

    # Get all labels for this subject
    if subject in ['02','10','14','31','39','40','P03','18','51','25','46','19','43','45']:  # subjects without MRI, use fsaverage
       labels = pymegsr.get_labels('fsaverage')
    else:
       labels = get_labels(subject, only_glasser)
    # And keep only the ones that we are interested in (area parameter of this function) -> needs to be in areas_to_labels dict defined above
    areas_to_labels = a2lab()
    labels = [x for x in labels if any([cl for cl in areas_to_labels[area] if cl in x.name])]
    print(labels)

    if len(labels) < 1:
        raise RuntimeError('Expecting at least two labels')

    # Turn individual labels into one big label that contains all vertices that
    # belong to one ROI
    label = labels.pop()
    for l in labels:
        label += l

    print('Selecting this label for area %s:'%area, label)

    data=[]; fwds=[]; bems=[]; sources=[];
    #for sess, rec in sessinfo:
    sess = 1 #only considering 1st session
    logging.info("Reading data for %s, delay %i "% (subject, rec))
    if rec == 1: # 1s delay -- all trials:
        data_cov,epoch,epoch_filename = get_stim_epoch(subject, sess, rec)
    else: # 3s and 9s trials are fully concatenated for the decoding analysis
        data_cov,epoch,epoch_filename = full_length_get_stim_epoch(subject, sess, rec)

    data.append((data_cov, epoch)); ##### PM: N.B. data may in fact need to be output currently called 'epochs'

    logging.info("Setting up source space and forward model")
    megdirs = '/home/gmonov/meg_data/'

    if subject in ['17','18','19','20','21','22','23','26','27','28','30','32','33','35','36','38','46','47','48','49','50','51','52','53','54','P03']: # hcs
              raw_filename = glob(megdirs + str(subject) + '_MCI' + '*_01' + '.ds')
    elif subject=='11' and sess==2:  # session has not been added properly in file name
              raw_filename = glob(megdirs + '11_MCI' + '*_01' + '.ds')
    elif subject=='14' and sess==1:  # Patients where second recording is the right one
              raw_filename = glob(megdirs + str(subject) + '-' + str(sess) + '*_02' + '.ds')
    elif subject=='43' and sess==2:  # Patients where second recording is the right one
             raw_filename = glob(megdirs + str(subject) + '-' + str(sess) + '*_02' + '.ds')
    elif subject=='31': # HC where second recording is the right one
              raw_filename = glob(megdirs + str(subject) + '_MCI' + '*_02' + '.ds')
    elif subject=='37': # HC where second recording is the right one
              raw_filename = glob(megdirs + str(subject) + '_MCI' + '*_02' + '.ds')
    elif subject=='24' and sess==1:  # Naming mistake and second recording had to be started
              raw_filename = glob(megdirs + '01C_MCI' + '*_02' + '.ds')
    elif subject=='25' and sess==1:  # Naming mistake
              raw_filename = glob(megdirs + '02C_MCI' + '*_01' + '.ds')
    else:
              raw_filename = glob(megdirs + str(subject) + '-' + str(sess) + '*_01' + '.ds')


    if subject in ['01','20','22','24','26','27','52']:  # these subjects have intersecting BEM layers
        conductivity=(0.3,)  # setting len(conductivity)==1 is a trick to make mne use only the inner skull layer (i.e. a 1-layer BEM)
    else:
        conductivity=(0.3, 0.006, 0.3) # otherwise, use the default settings


    raw_filename = raw_filename[0]
    trans_filename = glob('/home/pmurphy/meg_data/MCI/MRIs/trans_mats/%s_1*fif' % (
            subject))[0]
    if subject in ['02','10','14','31','39','40','P03','18','51','25','46','19','43','45']:  # subjects without MRI, use fsaverage
       forward, bem, source = get_leadfield(
           'fsaverage', raw_filename, epoch_filename, trans_filename, bem_sub_path='bem_ft', conductivity=conductivity)
       fwds.append(forward); # bems.append(bem); sources.append(source);
    else:
       forward, bem, source = get_leadfield(
           subject, raw_filename, epoch_filename, trans_filename, bem_sub_path='bem_ft', conductivity=conductivity)
       fwds.append(forward); # bems.append(bem); sources.append(source);


    # Define TFR parameters
    #fois = np.arange(8, 14, 1)   # PM: np.arange(36, 162, 4)
    lfois = np.arange(min_freq, max_freq+1, 1) #,np.arange(12, 31, 2)])   # lower resolution than typical TF plots, to lower processing demands
    #lfois = np.hstack([np.arange(1, 7, 1),np.arange(16, 31, 2)])   # NB: IF COMMENTED IN, THIS SETTING EXCLUDES 7-15 Hz!!!!
    tfr_params = {
        #"HF": {
        #    "foi": fois,
        #    "cycles": fois * 0.25,
        #    "time_bandwidth": 6 + 1,
        #    "n_jobs": 1,
        #    "est_val": fois,
        #    "est_key": "HF",
        #    "sf": 400,
        #    "decim": 20,
       # },
        "LF": {
            "foi": lfois,
            "cycles": lfois * 0.4,
            "time_bandwidth": 1 + 1,
            "n_jobs": 1,
            "est_val": lfois,
            "est_key": "LF",
            "sf": 400,
            "decim": 20,
        },
    }

    events = [d[1].events[:, 2] for d in data]
    events = np.hstack(events)

    # Compute LCMV filters for each session
    filters = []
    for (data_cov, epochs), forward in zip(data, fwds):
        filters.append(
            pymeglcmv.setup_filters(epochs.info, forward, data_cov, None, [label])
        )
    set_n_threads(1)

    # Specify vertex -> hemisphere mapping array --- COMMENT IN IF WANT TO AVERAGE ACROSS VERTICES
    f = filters[0][label.name]
    avg_vertices = np.zeros((len(f['vertices'][0]) + len(f['vertices'][1]))).astype(bool)
    avg_vertices[:len(f['vertices'][0])] = True

    # specify decoding settings
    clf = Pipeline(
        [
           ("Scaling", StandardScaler()),
           ("PCA", PCA(n_components=0.80, svd_solver='full')),
           ("RidgeReg", linear_model.Ridge(alpha=1)),
        ]
    )

    # specify sample onsets and window for decoding
    smpon = np.array([0])   # vector of sample onsets (s)
    if rec == 1:
        smpwin = [-0.10001, 1.50001]    # window for decoding, relative to sample onset (s) - going marginally outside desired bnds important to catch all times
    elif rec == 3:
        smpwin = [-0.10001, 3.50001]    # window for decoding, relative to sample onset (s) - going marginally outside desired bnds important to catch all times
    elif rec == 9:
        smpwin = [-0.10001, 9.50001]    # window for decoding, relative to sample onset (s) - going marginally outside desired bnds important to catch all times

    # load to-be-decoded variables & check that trials are appropiately aligned
    matname = ('/home/gmonov/meg_analysis/Decoding/behav4decoding/%s_1_%i.mat' % (subject,rec))
    mat = sio.loadmat(matname)

    mat_events = np.int64(np.concatenate(mat["tIDs"]))  # convert matlab events to same type as python events
    assert np.array_equal(events,mat_events[:len(events)])

    # Perform source reconstruction, using for each session the appropriate filter
    # Iterates over sample positions to mitigate memory demands
    all_smp = []  # inialize DataFrame containing all decoding results, across sample positions/variables
    for smp in range(len(smpon)):

        #fname = "/home/gmonov/meg_analysis/Decoding/decodeAv_alpha/%s_%s_%s_avTF.hdf" % (subject, area, str(rec))
        fname = "/home/gmonov/meg_analysis/Decoding/decode_%i_%i_%s/%s_%s_%i_%i_%i.hdf" % (min_freq, max_freq, mode, subject, area, rec, min_freq, max_freq)
        # the try: except: block implements caching, if output is already there, don't do it again.
        try:
            all_s = pd.read_hdf(fname)
        except FileNotFoundError:
            # perform source reconstruction of TF data
            #HF_tfrdata, events, HF_freq, times = decoding_GM.get_lcmv(   # PM: padding by 0.2s (max TF win / 2) for accurate TF estimation
            #    tfr_params["HF"], [d[1].copy().crop(smpon[smp]+smpwin[0]-0.2,smpon[smp]+smpwin[1]+0.2) for d in data], filters, njobs=6    #d[1].copy().crop() pulls out sample-aligned data
           # )
            LF_tfrdata, events, LF_freq, times = decoding_GM.get_lcmv(
                tfr_params["LF"], [d[1].copy().crop(smpon[smp]+smpwin[0]-0.2,smpon[smp]+smpwin[1]+0.2) for d in data], filters, njobs=6
            )

            # Concatenate data
            #tfrdata = np.hstack([HF_tfrdata, LF_tfrdata])
            tfrdata = np.hstack([LF_tfrdata])
            #del LF_tfrdata, HF_tfrdata
            del LF_tfrdata
            #freq = np.concatenate([HF_freq, LF_freq])
            freq = np.concatenate([LF_freq])

            # loop through variables to be decoded
            ctimes = (smpwin[0]+smpon[smp] <= times) & (times <= smpwin[1]+smpon[smp])  # indices of time-points without padding
            all_s = []
            for target in ["mem_loc"]:
                # pull variable to be decoded
                target_vals = mat[target]   #  target_vals will be a numpy ndarray, ntrials*nsamples
                target_vals = target_vals[:len(events),:]
                print('Size:', tfrdata.shape)

                # perform decoding
                dcd = decoding_GM.Decoder(target_vals[:,smp],("RidgeReg",clf))
                k = dcd.classify(
                    tfrdata[:,:,:,ctimes], times[ctimes]-smpon[smp], freq, events, area,   # feeding in times aligned to smp onset
                    average_vertices=False, use_phase=False            ####### PM: NB, set average_vertices to False if want to preserve vertices as separate features
                )
                k.loc[:, "target"] = target  # include target_val label in dataframe
                all_s.append(k)   #  append decoding results for this target_val combo

            all_s = pd.concat(all_s)         # concatenate all target_vals
            all_s.loc[:, 'ROI'] = area       # and include ROI/sample position labels
            all_s.loc[:, "sample"] = str(smp+1)
            all_s.to_hdf(fname, "df")  # save once all target_vals have been iterated over

        all_smp.append(all_s)  # append decoding results for this sample position

    all_smp = pd.concat(all_smp)  # concatenate all sample positions
    all_smp.to_csv("/home/gmonov/meg_analysis/Decoding/decode_%i_%i_%s/%s_%s_%i_%i_%i.csv" % (min_freq, max_freq, mode, subject, area, rec, min_freq, max_freq))
