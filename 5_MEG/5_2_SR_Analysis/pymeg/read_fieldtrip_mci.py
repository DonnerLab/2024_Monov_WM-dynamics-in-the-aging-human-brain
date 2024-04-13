
import h5py
import numpy as np
import os
import glob

from joblib import Memory

import mne

from mne import create_info
from mne.epochs import EpochsArray


memory = Memory(cachedir=os.environ['PYMEG_CACHE_DIR'], verbose=0)


def fix_chs(rawinfo, einfo):
    ch_names = []
    for k in range(len(einfo['chs'])):
        name = einfo['chs'][k]['ch_name']
        newchan = [x for x in rawinfo['chs']
                   if name in x['ch_name']][0]
        einfo['chs'][k] = newchan
        ch_names.append(newchan['ch_name'])
    einfo['ch_names'] = ch_names
    return einfo


@memory.cache
def get_info_for_epochs(rawname):
    raw = mne.io.ctf.read_raw_ctf(rawname)
    return raw.info


def read_ft_epochs(fname, rawinfo, cachedir=os.environ['PYMEG_CACHE_DIR'],
                   trialinfo_col=-1):
    '''
    Read and cache the output of fieldtrip epochs.

    This function reads a matlab file that contains a 'data' struct with the
    following fields:

        trialinfo: matrix
            Dim is ntrials x nr_meta, the columns contain meta information
            about each trial.
        label: list of strings
            Channel names
        sampleinfo: matrix
            Dim is ntrials x 2, the first column contains the start sample
            of each epoch in the raw data.
        time: array
            Contains time points for epochs.
        trial: array
            Dim is time x channels x trials, contains the actial data

    This data is parsed into an MNE Epochs object. To correctly assign
    channel locations, types etc. the info structure from the raw data
    that generated the fieldtrip epochs is used. The channel names in
    the fieldtrip structure should still be relatable to the raw
    channel names, relatable here means that a fieldtrip channel name
    must be contained in the raw channel name.

    Args
        fname: str
            Path to .mat file to load
        rawinfo: mne info structure
            Info structure with correct channel locations etc. This
            should be obtained by reading the raw data corresponding
            to the epochs with MNE.
        cachedir: str
            Path where the epochs are saved on disk. If this is
            None the epochs are returned.
        trialinfo_col: int
            Column in trialinfo which contains trial identifier.

    Output
        Returns path to saved epochs if cachedir is not None, else
        it returns the epochs
    '''
    if cachedir is None:
        return _load_ft_epochs(fname, rawinfo, trialinfo_col=trialinfo_col)
    epochs_path = os.path.join(cachedir, fname + '-epo.fif.gz')
    if not os.path.exists(epochs_path):
        epochs = _load_ft_epochs(fname, rawinfo, trialinfo_col=trialinfo_col)
        epochs.save(epochs_path)
    return epochs_path


def _load_ft_epochs(fname, rawinfo, trialinfo_col=-1):
    # load Matlab/Fieldtrip data
    f = h5py.File(fname)
    list(f.keys())
    ft_data = f['data']
    ft_data.keys()

    trialinfo = ft_data['trialinfo']
    channels = ft_data['label']
    sampleinfo = ft_data['sampleinfo']
    time = ft_data['time']
    sfreq = np.around(1 / np.diff(time[:].ravel()), 2)
    assert(len(np.unique(sfreq)) == 1)
    n_time, n_chans, n_trial = ft_data['trial'].shape

    data = np.zeros((n_trial, n_chans, n_time))
    transposed_data = np.transpose(ft_data['trial'])
    for trial in range(n_trial):
        data[trial, :, :] = transposed_data[trial]

    data = data[:, range(n_chans), :]

    chan_names = []
    for i in range(n_chans):
        st = channels[0][i]
        obj = ft_data[st]
        chan_names.append(''.join(chr(j) for j in obj[:]))
    #ch_names = [x + '-3705' for x in chan_names]

    info = create_info(chan_names, sfreq[0])
    events = np.zeros((n_trial, 3), int)
    events[:, 2] = trialinfo[trialinfo_col]
    events[:, 0] = sampleinfo[0]

    epochs = EpochsArray(data, info, events=events, tmin = min(time)[0], verbose=False)
    epochs.info = fix_chs(rawinfo, epochs.info)
    return epochs


#subjects = {'01': [(1, 3)],
#            '02': [(1, 3), (1, 9)],
#            '03': [(1, 3), (1, 9)],
#            '04': [(1, 3), (1, 9)],
#            '05': [(1, 3), (1, 9)],
#            '06': [(1, 3), (1, 9)],
#            '07': [(1, 3), (1, 9)],
#            '08': [(1, 3), (1, 9)],
#            '09': [(1, 3), (1, 9)],
#            '10': [(1, 3), (1, 9)],
#            '11': [(1, 3), (1, 9)]}
#subjects = {'14': [(1, 3), (1, 9)],
#            '15': [(1, 3), (1, 9)],
#            '16': [(1, 3), (1, 9)],
#            '17': [(1, 3), (1, 9)],
#            '18': [(1, 3), (1, 9)],
#            '19': [(1, 3), (1, 9)],
#            '20': [(1, 3), (1, 9)],
#            '21': [(1, 3), (1, 9)],
#            '22': [(1, 3), (1, 9)],
#            '23': [(1, 3), (1, 9)],
#            '24': [(1, 3), (1, 9)],
#            '25': [(1, 3), (1, 9)],
#            '26': [(1, 3), (1, 9)],
#            '27': [(1, 3), (1, 9)]}
subjects = {'28': [(1, 3), (1, 9)],
            '29': [(1, 3), (1, 9)],
            '30': [(1, 3), (1, 9)],
            '31': [(1, 3), (1, 9)],
            '32': [(1, 3), (1, 9)],
            '33': [(1, 3), (1, 9)],
            '35': [(1, 3), (1, 9)],
            '36': [(1, 3), (1, 9)],
            '37': [(1, 3), (1, 9)],
            '38': [(1, 3), (1, 9)],
            '39': [(1, 3), (1, 9)],
            '40': [(1, 3), (1, 9)],
            '41': [(1, 3), (1, 9)],
            '42': [(1, 3), (1, 9)],
            '43': [(1, 3), (1, 9)],
            '45': [(1, 3), (1, 9)],
            '46': [(1, 3), (1, 9)],
            '47': [(1, 3), (1, 9)],
            '48': [(1, 3), (1, 9)],
            '49': [(1, 3), (1, 9)],
            '50': [(1, 3), (1, 9)],
            '51': [(1, 3), (1, 9)],
            '52': [(1, 3), (1, 9)],
            '53': [(1, 3), (1, 9)],
            '54': [(1, 3), (1, 9)],
            'P03': [(1, 3), (1, 9)]}
full_length = 1 #set to 0 for separated times and 1 for full length long delay trials
if full_length == 0:
   ftdir = '/home/gmonov/meg_analysis/sr_separated_data/'
elif full_length == 1:
   ftdir = '/home/gmonov/meg_analysis/full_length_3_9s/'
megdirs = '/home/gmonov/meg_data/'
if full_length == 0:
    savedir = '/home/gmonov/meg_analysis/source_reconstruction/Conv2mne/'
elif full_length == 1:
    savedir = '/home/gmonov/meg_analysis/source_reconstruction/full_length_Conv2mne/'

for subj, tasks in subjects.items():
    for sess, rec in tasks:
        ftname = ftdir + str(subj) + '_' + str(sess) + '/' + str(subj) + '_' + str(sess) +'_' + str(rec) + '.mat';
        print(ftname)
        if subj in ['17','18','19','20','21','22','23','26','27','28','30','32','33','35','36','38','46','47','48','49','50','51','52','53','54','P03']: # hcs
            megname = glob.glob(megdirs + str(subj) + '_MCI' + '*_01' + '.ds')
        elif subj=='11' and sess==2:  # session has not been added properly in file name
            megname = glob.glob(megdirs + '11_MCI' + '*_01' + '.ds')
        elif subj=='14' and sess==1:  # Patients where second recording is the right one
            megname = glob.glob(megdirs + str(subj) + '-' + str(sess) + '*_02' + '.ds')
        elif subj=='43' and sess==2:  # Patients where second recording is the right one
            megname = glob.glob(megdirs + str(subj) + '-' + str(sess) + '*_02' + '.ds')
        elif subj=='31': # HC where second recording is the right one
            megname = glob.glob(megdirs + str(subj) + '_MCI' + '*_02' + '.ds')
        elif subj=='37': # HC where second recording is the right one
            megname = glob.glob(megdirs + str(subj) + '_MCI' + '*_02' + '.ds')
        elif subj=='24' and sess==1:  # Naming mistake and second recording had to be started
            megname = glob.glob(megdirs + '01C_MCI' + '*_02' + '.ds')
        elif subj=='25' and sess==1:  # Naming mistake
            megname = glob.glob(megdirs + '02C_MCI' + '*_01' + '.ds')
        else:
            megname = glob.glob(megdirs + str(subj) + '-' + str(sess) + '*_01' + '.ds')
        print(megname)
        megname = str(megname)[1:-1]
        megname = megname[1:-1]
        rawinfo = get_info_for_epochs(megname)
        epochs = _load_ft_epochs(ftname, rawinfo)
        epochs.save(savedir + str(subj) + '_' + str(sess) + '_' + str(rec) + '_preproc4mne_induced.mat' + '-epo.fif.gz')
        del epochs
