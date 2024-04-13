# make imports
import numpy as np
from sklearn import svm
from sklearn.model_selection import (
    cross_validate,
    cross_val_predict,
    RandomizedSearchCV,
)
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.decomposition import PCA
from sklearn import linear_model
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
import time
import os
import pandas as pd
import scipy.io as sio
import decode_stim_loc


def decode(
    subject, sampling_rate
):
    clf = Pipeline(
        [
           ("Scaling", StandardScaler()),
           ("RidgeReg", linear_model.Ridge(alpha=1)),
        ]
    )

# load data from subject
# load sample stimulus location
    print(subject)

# load preprocessed eye tracking data

    matname = ('/home/gmonov/eye_tracking_data_analysis/data4decoding/trials4decoding/%s_trial_data_x_y_%i.mat' % (subject,sampling_rate))
    mat = sio.loadmat(matname)
    data = mat["trial_data_x_y"]
    times = mat["times"]
    target_vals = mat["stim_loc"]
    target_vals = target_vals[:,:]
    all_smp = []
# perform decoding
    dcd = decode_stim_loc.Decoder(target_vals[:,0],("RidgeReg",clf))
    k = dcd.classify(
    data, times[:,0]
    )
    all_smp = k
    del data
    del target_vals

    all_smp.to_csv("/home/gmonov/eye_tracking_data_analysis/decoding_results/decode_%s_%i_ndr_x_y.csv" % (subject,sampling_rate))
