import logging
import mne
import numpy as np
import os

from joblib import Memory

from os import makedirs
from os.path import join
os.environ["PYMEG_CACHE_DIR"] = "/mnt/homes/home024/pmurphy/tmp"
from glob import glob

from pymeg import lcmv as pymeglcmv
from pymeg import source_reconstruction as pymegsr

memory = Memory(cachedir=os.environ['PYMEG_CACHE_DIR'])
path = '/home/gmonov/meg_analysis/source_reconstruction/Conv2mne/'


def set_n_threads(n):
    import os
    os.environ['OPENBLAS_NUM_THREADS'] = str(n)
    os.environ['MKL_NUM_THREADS'] = str(n)
    os.environ['OMP_NUM_THREADS'] = str(n)

def select_delay_trials(all_trial_ids,epochs):
    import pandas as pd
    trial_no = pd.DataFrame(all_trial_ids)
    log_ind = trial_no.isin(epochs.events[:,2])
    log_ind = log_ind[log_ind[0] == True]
    selected_trials = log_ind.index

    return selected_trials


#subjects = {'01': [(1, 1), (1, 3), (2, 1), (2, 3), (2, 9)]}#,
#subjects = {'02': [(1, 1)]}#,, (1, 3), (1, 9), (2, 1), (2, 3), (2, 9)]}#,
#subjects = {'03': [(1, 1), (1, 3), (1, 9)]}#,
#subjects = {'04': [(1, 1), (1, 3), (1, 9), (2, 1), (2, 3), (2, 9)],
#subjects = {'05': [(1, 1)]}#,, (1, 3), (1, 9), (2, 1), (2, 3), (2, 9)]}#,
#subjects = {'06': [(1, 1), (1, 3), (1, 9), (2, 1), (2, 3), (2, 9)]}#,
#subjects = {'07': [(1, 1), (1, 3), (1, 9), (2, 1), (2, 3), (2, 9)]}#,,
#subjects = {'08': [(1, 1), (1, 3), (1, 9), (2, 1), (2, 3), (2, 9)]}#,,
#subjects = {'09': [(1, 1), (1, 3), (1, 9), (2, 1), (2, 3), (2, 9)]}#,,
#subjects = {'11': [(1, 1), (1, 3), (1, 9)]}#,,
#subjects = {'18': [(1, 1), (1, 3), (1, 9)]}#,,
#subjects = {'19': [(1, 1), (1, 3), (1, 9)]}#,,
#subjects = {'24': [(1, 1), (1, 3), (1, 9)]}#,,,
#subjects = {'25': [(1, 1), (1, 3), (1, 9)]}#,,,
#subjects = {'32': [(1, 1)]}#,,,
#subjects = {'43': [(1, 1),(1, 3), (1, 9)]}#,,,
#subjects = {'45': [(1, 1), (1, 3), (1, 9)]}#,,,
#subjects = {'46': [(1, 1), (1, 3), (1, 9)]}#,,,
#subjects = {'51': [(1, 1), (1, 3), (1, 9)]}#,,
#subjects = {'11': [(1, 1), (1, 3), (1, 9), (2, 1), (2, 3), (2, 9)]}#,,
#subjects = {'14': [(1, 1), (1, 3), (1, 9)]}#,,
#subjects = {'15': [(1, 1), (1, 3), (1, 9)]}#,,
#subjects = {'16': [(1, 1), (1, 3), (1, 9), (2, 1), (2, 3), (2, 9)]}#,,
#subjects = {'17': [(1, 1), (1, 3), (1, 9)],#,,
subjects = {'43': [(1, 1), (1, 3), (1, 9)],#,,
            '45': [(1, 1), (1, 3), (1, 9)]}#,,
#subjects = {'19': [(1, 1), (1, 3), (1, 9)]}#,,
#subjects = {'20': [(1, 1), (1, 3), (1, 9)]}#,,
#subjects = {'21': [(1, 1), (1, 3), (1, 9)]}#,,
#subjects = {'22': [(1, 1), (1, 3), (1, 9)]}#,,
#subjects = {'23': [(1, 1), (1, 3), (1, 9)]}#,,
#subjects = {'24': [(1, 1), (1, 3), (1, 9)]}#,,
#subjects = {'25': [(1, 1), (1, 3), (1, 9)]}#,,
#subjects = {'26': [(1, 1), (1, 3), (1, 9)]}#,,
#subjects = {'27': [(1, 1), (1, 3), (1, 9)]}#,,
#subjects = {'28': [(1, 1), (1, 3), (1, 9)]}#,,
#subjects = {'29': [(1, 1), (1, 3), (1, 9), (2, 1), (2, 3), (2, 9)]}#,,
#subjects = {'30': [(1, 1), (1, 3), (1, 9)]}#,,
#subjects = {'31': [(1, 1), (1, 3), (1, 9)]}#,,,
#subjects = {'32': [(1, 1), (1, 3), (1, 9)]}#,,,
#subjects = {'33': [(1, 1), (1, 3), (1, 9)]}#,,,
#subjects = {'35': [(1, 1), (1, 3), (1, 9)]}#,,,
#subjects = {'36': [(1, 1), (1, 3), (1, 9)]}#,,,
#subjects = {'37': [(1, 1), (1, 3), (1, 9)]}#,,,
#subjects = {'38': [(1, 1), (1, 3), (1, 9)]}#,,
#subjects = {'39': [(1, 1), (1, 3), (1, 9)]}#,
#subjects = {'40': [(1, 1), (1, 3), (1, 9)]}#,,
#subjects = {'41': [(1, 1), (1, 3), (1, 9), (2, 1), (2, 3), (2, 9)]
#subjects = {'42': [(1, 1), (1, 3), (1, 9)]}#,,,,
#subjects = {'43': [(1, 1), (1, 3), (1, 9), (2, 1), (2, 3), (2, 9)]}#,,,,
#subjects = {'45': [(1, 1), (1, 3), (1, 9)]}#,,,,
#subjects = {'46': [(1, 1), (1, 3), (1, 9)]}#,,,,
#subjects = {'47': [(1, 1), (1, 3), (1, 9)]}#,,,,
#subjects =  {'48': [(1, 1), (1, 3), (1, 9)]}#,,,,
#subjects =  {'49': [(1, 1), (1, 3), (1, 9)]}#,,,,
#subjects =  {'50': [(1, 1), (1, 3), (1, 9)]}
#subjects =  {'51': [(1, 1), (1, 3), (1, 9)]}
#subjects =  {'52': [(1, 1), (1, 3), (1, 9)]}#,,,,
#subjects =  {'53': [(1, 1), (1, 3), (1, 9)]}
#subjects =  {'54': [(1, 1), (1, 3), (1, 9)]}
#subjects =  {'P03': [(1, 1), (1, 3), (1, 9)]}



def submit():
    from pymeg import parallel
    for subject, tasks in subjects.items():
        for session, recording in tasks:
            #for signal in ['LF','HF']:
            for signal in ['BB', 'HF', 'LF']:
                parallel.pmap(
                    extract, [(subject, session, recording, signal)],
                    cluster = 'SLURM', walltime='500:00:00', memory=8000, nodes=1, tasks=8,
                    name='sr' + str(subject) + '_' + str(session) + str(recording) + '_' + signal, env='mne')


def lcmvfilename(subject, session, signal, recording, chunk=None):
    try:
        makedirs(path)
    except:
        pass
    if chunk is None:
        filename = '%s-SESS%i-%i-%s-lcmv.hdf' % (
            subject, session, recording, signal)
    else:
        filename = '%s-SESS%i-%i-%s-chunk%i-lcmv.hdf' % (
            subject, session, recording, signal, chunk)
    return join(path, filename)


def get_stim_epoch(subject, session, recording):
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

    globstring = '/home/gmonov/meg_analysis/source_reconstruction/Conv2mne/%s_%i_%i*fif.gz' % (
        subject, session, recording)  # only loading the 1s file consisting of beginning of each clean trial for estimating the covariance matrix
    filenames = glob(globstring)[0]
    epochs = mne.read_epochs(filenames)
    # epochs.times = epochs.times - 1  # PM: this was somehow necessary for initial pipeline, but *NOT* for induced
    epochs = epochs.pick_channels(
        [x for x in epochs.ch_names if x.startswith('M')])

    selected_trials = select_delay_trials(
           all_trial_ids,epochs)


    epochs._data -= means[selected_trials, :, np.newaxis]   # need to select only trials that match those for current delay duration

    return data_cov, epochs, filenames


def extract(subject, session, recording, signal_type='BB',
            BEM='three_layer', debug=False, chunks=100, njobs=4):
    mne.set_log_level('WARNING')
    pymeglcmv.logging.getLogger().setLevel(logging.INFO)
    set_n_threads(1)

    logging.info('Reading stimulus data')
    data_cov, epochs, epochs_filename = get_stim_epoch(
        subject, session, recording)



    megdirs = '/home/gmonov/meg_data/'

    if subject in ['17','18','19','20','21','22','23','26','27','28','30','32','33','35','36','38','46','47','48','49','50','51','52','53','54','P03']: # hcs
          raw_filename = glob(megdirs + str(subject) + '_MCI' + '*_01' + '.ds')
    elif subject=='11' and session==2:  # session has not been added properly in file name
          raw_filename = glob(megdirs + '11_MCI' + '*_01' + '.ds')
    elif subject=='14' and session==1:  # Patients where second recording is the right one
          raw_filename = glob(megdirs + str(subject) + '-' + str(session) + '*_02' + '.ds')
    elif subject=='43' and session==2:  # Patients where second recording is the right one
          raw_filename = glob(megdirs + str(subject) + '-' + str(session) + '*_02' + '.ds')
    elif subject=='31': # HC where second recording is the right one
          raw_filename = glob(megdirs + str(subject) + '_MCI' + '*_02' + '.ds')
    elif subject=='37': # HC where second recording is the right one
          raw_filename = glob(megdirs + str(subject) + '_MCI' + '*_02' + '.ds')
    elif subject=='24' and session==1:  # Naming mistake and second recording had to be started
          raw_filename = glob(megdirs + '01C_MCI' + '*_02' + '.ds')
    elif subject=='25' and session==1:  # Naming mistake
          raw_filename = glob(megdirs + '02C_MCI' + '*_01' + '.ds')
    else:
          raw_filename = glob(megdirs + str(subject) + '-' + str(session) + '*_01' + '.ds')



    assert len(raw_filename) == 1
    raw_filename = raw_filename[0]
    print(epochs_filename)
    #epochs_filename = glob('/home/gmonov/meg_analysis/source_reconstruction/Conv2mne/' + str (subject) + '_' + str(session) + '_' + str(recording)+ '*fif.gz') # redefine filename for 3s and 9s delay files
    #epochs_filename = epochs_filename[0]
    #print(epochs_filename)
    trans_filename = glob('/home/pmurphy/meg_data/MCI/MRIs/trans_mats/%s_%i*fif' % (
        subject, session))[0]
    #if subject in ['51']:  # use the resting state trans_mat for this subject
    #    trans_filename = glob('/mnt/homes/home032/aarazi/MCI/data/trans_mat/%s_%i*fif' % (
    #        subject, session))[0]
    #else:
    #    trans_filename = glob('/home/pmurphy/meg_data/MCI/MRIs/trans_mats/%s_%i*fif' % (
    #        subject, session))[0]
    logging.info('Setting up source space and forward model')

    if subject in ['01','20','22','24','26','27','52']:  # these subjects have intersecting BEM layers, removed 18 51, 19, 25,45 since they run with fsaverage model
        conductivity=(0.3,)  # setting len(conductivity)==1 is a trick to make mne use only the inner skull layer (i.e. a 1-layer BEM)
    else:
        conductivity=(0.3, 0.006, 0.3) # otherwise, use the default settings

    if subject in ['02','10','14','31','39','40','P03','18','51','25','46','19','43','45']:  # subjects without MRI, use fsaverage # added 18 and 51 due to poor mri data quality
        forward, bem, source = pymegsr.get_leadfield(
            'fsaverage', raw_filename, epochs_filename, trans_filename, bem_sub_path='bem_ft', conductivity=conductivity)
        labels = pymegsr.get_labels('fsaverage')
    else:
        forward, bem, source = pymegsr.get_leadfield(
            subject, raw_filename, epochs_filename, trans_filename, bem_sub_path='bem_ft', conductivity=conductivity)
        labels = pymegsr.get_labels(subject)

    labels = pymegsr.labels_exclude(labels,
                                    exclude_filters=['wang2015atlas.IPS4',
                                                     'wang2015atlas.IPS5',
                                                     'wang2015atlas.SPL',
                                                     'JWDG_lat_Unknown'])
    #labels = pymegsr.labels_remove_overlap(
    #    labels, priority_filters=['wang', 'JWDG'],)

    # Now chunk Reconstruction into blocks of ~100 trials to save Memory
    fois_h = np.arange(36, 162, 4)
    fois_l = np.arange(1, 36, 1)
    tfr_params = {
        'HF': {'foi': fois_h, 'cycles': fois_h * 0.25, 'time_bandwidth': 6 + 1,
               'n_jobs': njobs, 'est_val': fois_h, 'est_key': 'HF', 'sf': 400,
               'decim': 20},
        'LF': {'foi': fois_l, 'cycles': fois_l * 0.4, 'time_bandwidth': 1 + 1,
               'n_jobs': njobs, 'est_val': fois_l, 'est_key': 'LF', 'sf': 400,
               'decim': 20}
    }

    events = epochs.events[:, 2]
    filters = pymeglcmv.setup_filters(epochs.info, forward, data_cov, None, labels)

    set_n_threads(1)

    for i in range(0, len(events), chunks):
        filename = lcmvfilename(
            subject, session, signal_type, recording, chunk=i)
        if os.path.isfile(filename):
            continue
        if signal_type == 'BB':
            logging.info('Starting reconstruction of BB signal')
            M = pymeglcmv.reconstruct_broadband(
                filters, epochs.info, epochs._data[i:i + chunks],
                events[i:i + chunks],
                epochs.times, njobs=1)
        else:
            logging.info('Starting reconstruction of TFR signal')
            M = pymeglcmv.reconstruct_tfr(
                filters, epochs.info, epochs._data[i:i + chunks],
                events[i:i + chunks], epochs.times,
                est_args=tfr_params[signal_type],
                njobs=4)
        M.to_hdf(filename, 'epochs')
    set_n_threads(njobs)
