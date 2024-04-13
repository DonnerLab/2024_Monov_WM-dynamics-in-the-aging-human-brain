
import os
os.environ["PYMEG_CACHE_DIR"] = "/mnt/homes/home024/pmurphy/tmp"

#subjects = {'01': [(1, 1), (1, 3), (2, 1), (2, 3), (2, 9)]}#,,
#subjects = {'02': [(1, 1), (1, 3), (1, 9), (2, 1), (2, 3), (2, 9)],
#            '03': [(1, 1), (1, 3), (1, 9)],
#            '04': [(1, 1), (1, 3), (1, 9), (2, 1), (2, 3), (2, 9)],
#            '05': [(1, 1), (1, 3), (1, 9), (2, 1), (2, 3), (2, 9)],
#            '06': [(1, 1), (1, 3), (1, 9), (2, 1), (2, 3), (2, 9)]}#,,
#subjects = {'07': [(1, 1), (1, 3), (1, 9), (2, 1), (2, 3), (2, 9)],
#            '08': [(1, 1), (1, 3), (1, 9), (2, 1), (2, 3), (2, 9)],
#            '09': [(1, 1), (1, 3), (1, 9), (2, 1), (2, 3), (2, 9)]}#,,
#subjects = {'10': [(1, 1), (1, 3), (1, 9)]}#,,
#subjects = {'11': [(1, 1), (1, 3), (1, 9), (2, 1), (2, 3), (2, 9)],
#           '14': [(1, 1), (1, 3), (1, 9)],
#            '15': [(1, 1), (1, 3), (1, 9)],
#            '16': [(1, 1), (1, 3), (1, 9), (2, 1), (2, 3), (2, 9)],
#            '17': [(1, 1), (1, 3), (1, 9)],
#            '18': [(1, 1), (1, 3), (1, 9)],
#            '19': [(1, 1), (1, 3), (1, 9)]}#,,
#subjects = {'20': [(1, 1), (1, 3), (1, 9)],
#            '21': [(1, 1), (1, 3), (1, 9)],
#            '22': [(1, 1), (1, 3), (1, 9)],
#            '23': [(1, 1), (1, 3), (1, 9)],
#            '24': [(1, 1), (1, 3), (1, 9)],
#            '25': [(1, 1), (1, 3), (1, 9)],
#            '26': [(1, 1), (1, 3), (1, 9)]}#,,
#subjects = {'27': [(1, 1), (1, 3), (1, 9)],
#            '28': [(1, 1), (1, 3), (1, 9)],
#            '29': [(1, 1), (1, 3), (1, 9), (2, 1), (2, 3), (2, 9)],
#            '30': [(1, 1), (1, 3), (1, 9)],
#            '31': [(1, 1), (1, 3), (1, 9)],
#            '32': [(1, 1), (1, 3), (1, 9)],
#            '33': [(1, 1), (1, 3), (1, 9)],
#            '35': [(1, 1), (1, 3), (1, 9)],
#            '36': [(1, 1), (1, 3), (1, 9)]}#,,,
#subjects = {'37': [(1, 1), (1, 3), (1, 9)],
#            '38': [(1, 1), (1, 3), (1, 9)],
#            '39': [(1, 1), (1, 3), (1, 9)],
#            '40': [(1, 1), (1, 3), (1, 9)],
#            '41': [(1, 1), (1, 3), (1, 9), (2, 1), (2, 3), (2, 9)],
#            '42': [(1, 1), (1, 3), (1, 9)],
#            '43': [(1, 1), (1, 3), (1, 9), (2, 1), (2, 3), (2, 9)],
#            '45': [(1, 1), (1, 3), (1, 9)],
#            '46': [(1, 1), (1, 3), (1, 9)],
#            '47': [(1, 1), (1, 3), (1, 9)]}#,,,
#subjects = {'48': [(1, 1), (1, 3), (1, 9)],
#            '49': [(1, 1), (1, 3), (1, 9)],
#            '50': [(1, 1), (1, 3), (1, 9)],
#subjects = {'51': [(1, 1), (1, 3), (1, 9)]}
#            '52': [(1, 1), (1, 3), (1, 9)],
#            '53': [(1, 1), (1, 3), (1, 9)],
#            '54': [(1, 1), (1, 3), (1, 9)],
#            'P03': [(1, 1), (1, 3), (1, 9)]}
#subjects = {'11': [(1, 1), (1, 3), (1, 9)]}#,,
#subjects = {'18': [(1, 1), (1, 3), (1, 9)]}#,,
#subjects = {'19': [(1, 1), (1, 3), (1, 9)]}#,,
#subjects = {'24': [(1, 1), (1, 3), (1, 9)]}#,,
#subjects = {'25': [(1, 1), (1, 3), (1, 9)]}#,,
#subjects = {'32': [(1, 1), (1, 3), (1, 9)]}#,,
#subjects = {'43': [(1, 1), (1, 3), (1, 9)]}#,,
#subjects = {'45': [(1, 1), (1, 3), (1, 9)]}#,,
#subjects = {'46': [(1, 1), (1, 3), (1, 9)]}#,,
subjects = {'45': [(1, 1), (1, 3), (1, 9)]}#,,

def submit_aggregates(cluster='uke'):
    from pymeg import parallel
    from os.path import join
    for subject, tasks in subjects.items():
        for session, delay in tasks:
        #for sessnum in range(1,session+1):
            #for runnum in range(1,3):
                #for datatype in ['F']:
                for datatype in ['F','BB']:
                    parallel.pmap(aggregate, [(subject, session, delay, datatype)],
                                  cluster='SLURM', walltime='24:00:00', memory=8000, nodes = 1, tasks=8,
                                  name='agg' + str(session) + '-' + str(delay)+ '-' + str(subject) + datatype, env='mne')

def aggregate(subject, session, delay, datatype):
    from pymeg import aggregate_sr as asr
    from os.path import join
    bl = 1 # take baseline only from 1s trials
    data = (
        '/home/gmonov/meg_analysis/source_reconstruction/Conv2mne/%s-SESS%i-%i*%s*-lcmv.hdf' % (
            subject, session, delay, datatype))
    data_baseline=(
        '/home/gmonov/meg_analysis/source_reconstruction/Conv2mne/%s-SESS%i-%i*%s*-lcmv.hdf' % (
            subject, session, bl, datatype))

    if datatype == 'F':    # time-frequency
        agg = asr.aggregate_files(data, data_baseline, (-0.4, -0.2), to_decibels=True)
    elif datatype == 'BB':    # broadband
        agg = asr.aggregate_files(data, data_baseline, (-0.2, 0), to_decibels=False)
    elif datatype == 'G':    # gamma (optimized)
        agg = asr.aggregate_files(data, data, (-0.1, -0.05), to_decibels=True)

    filename = join(
        '/home/gmonov/meg_analysis/source_reconstruction/agg/',
        'S%s_SESS%i_%i_%s_agg.hdf' % (subject, session, delay, datatype))
    asr.agg2hdf(agg, filename)


if __name__=="__main__":
    submit_aggregates()
