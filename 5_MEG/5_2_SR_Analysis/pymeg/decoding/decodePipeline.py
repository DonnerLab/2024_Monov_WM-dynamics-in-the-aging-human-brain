"""
Functions implimenting various stages of a decoding pipeline

COMMON INPUTS
Many functions require the same input...

dirs: dict. The values are full filenames (including path) to use for saving
    and loading.

HISTORY
Setup, 2021, Joshua Calder-Travis, j.calder.travis@gmail.com
"""

import time
import numpy as np
import mne
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
import sklearn.model_selection as modelSelec
import sklearn.svm as svm
import sklearn.neural_network as NN
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import pandas as pd
from joblib import Parallel, delayed
import copy
from pymeg import helpers as helpKit


def crossValEvalDecoder(predictors, outcome, timePoints,
                            numFolds, pipeline,
                            trainTimeIdx, timeInvWindow,
                            splitVar=None,
                            evalMode='atAllTimes',
                            njobs=1):
    """ Evaluate the performance of a decoder through cross-validation

    INPUT
    predictors: numpy array. Shape must be
        (num trials or cases, num features, num time points)
    outcome: numpy vector. Contains the values we are aiming to predict. Shape
        should be (num trials or cases,)
    timePoints: numpy vector. Contains time each time point corresponds to.
        Hence, should have shape (numy time points,)
    numFolds: The number of cross validaton folds to use.
    pipeline: str. Sprecifies the scikit-learn decoder pipeline to use. For
        options see findDecodePipeline. Options include both classification
        and regression.
    trainTimeIdx: str or int. 'all' to train seperate decoders on data from
        each and every time point. 'inv' to train a single decoder on data
        from many time points. Alternatively, int specifying the index of
        the time point to use.
    timeInvWindow: 2-length tuple of scalars. If using the time invarient
        decoder (i.e. option 'inv' for trainTimeIdx), which time points
        should be used. This input specifies the window. Code uses the input
        'timePoints' to work out which data falls in the window.
    splitVar: None | array. If None has no effect. Otherwise shape must match
        the shape of outcome. Gives for each case a category (numeric or
        string) to which that case belongs. Data from different categories
        will be weighted during training to ensure each category is taking
        into account equally during traing, and cross validation will be
        stratified to ensure a consistent ratio of cases from the two
        categories. Accuracy data will be evaluated seperately for the two
        categories.
    evalMode: str. If 'atTrainTime' only evaluate the performace of a decoder
        using data from the same time point that was used for training (not
        avaliable if trainTimeIdx='inv'). If 'atAllTimes' evaluate each
        decoder using data from all time points, seperately for each time
        point.
    njobs: positive int. Number of jobs to run in parallel when using
        parallelised computing.

    OUTPUT
    combinedScoreData: Pandas data frame. Contains a single column of decoder
        performance values (accuracy for classifiers, and the coeficient of
        determination for regressions), and a multi-level index. The levels
        are the time point used for training the decoder ('Time_trained'), and
        the time point used for evaluating the decoder ('Time_evaluated'), and
        the cross validation fold from which the data are from ('CV_fold'). If
        splitVar is provided then accuracy is evalated seperately for each
        category and there is an extra index level 'Split_category' to give the
        cateogry.
    """
    # Create cross validation folds and loop over these
    if splitVar is None:
        kFold = modelSelec.KFold(n_splits=numFolds, shuffle=True)
        categories = None
    else:
        kFold = modelSelec.StratifiedKFold(n_splits=numFolds, shuffle=True)
        categories = splitVar

    jobs = []
    tic = time.perf_counter()
    for foldNum, (trainIdxs, testIdxs) in enumerate(
                                        kFold.split(predictors, categories)):
        jobs.append(delayed(runOneFold)(
            foldNum, trainIdxs, testIdxs, predictors, outcome, timePoints,
               pipeline, trainTimeIdx, timeInvWindow, splitVar, evalMode
        ))

    allScoreData = Parallel(n_jobs=njobs, verbose=1)(jobs)
    assert len(allScoreData) == numFolds
    toc = time.perf_counter()
    print('All CV folds ({}) evlauated using {} workers in {} '
          'seconds.'.format(numFolds, njobs, toc - tic))

    combinedScoreData = pd.concat(allScoreData, verify_integrity=True)

    expectedIndex = ['Time_trained', 'Time_evaluated', 'CV_fold']
    if splitVar is not None:
        expectedIndex = ['Split_category'] + expectedIndex
    helpKit.checkDfLevels(combinedScoreData, indexLvs=expectedIndex)
    assert (set(combinedScoreData.columns) == set(['Accuracy'])) or (
        set(combinedScoreData.columns) == set(['correlation_coefficient']))

    return combinedScoreData


def runOneFold(foldNum, trainIdxs, testIdxs, predictors, outcome, timePoints,
               pipeline, trainTimeIdx, timeInvWindow, splitVar, evalMode):
    print('Beginning CV fold: {}'.format(foldNum))
    scoreData, _ = allTimesDecode(predictors, outcome, timePoints,
                                    trainIdxs, testIdxs, pipeline,
                                    trainTimeIdx=trainTimeIdx,
                                    timeInvWindow=timeInvWindow,
                                    splitVar=splitVar,
                                    evalPerformance=evalMode)

    scoreData.loc[:, 'CV_fold'] = foldNum
    scoreData = scoreData.set_index('CV_fold', append=True)
    return scoreData


def allTimesDecode(predictors, outcome, timePoints, trainIdxs, testIdxs,
                    pipeName='default', trainTimeIdx='all',
                    timeInvWindow=(-float("Inf"), float("Inf")),
                    evalPerformance='atAllTimes',
                    splitVar=None):
    """ Function performs a decoding analysis on data with a time dimension.
    Trains seperate decoders for each time point, or a single decoder, and
    evaluates performance at all time points.

    INPUT
    predictors: numpy array. Shape must be
        (num trials or cases, num features, num time points)
    outcome: numpy vector. Contains the values we are aiming to predict. Shape
        should be (num trials or cases,). The values may be strings.
    timePoints: numpy vector. Contains time each time point corresponds to.
        Hence, should have shape (numy time points,)
    trainIdxs: numpy vector. Contains the indicies of the trials/cases to use
        for training all the decoders
    testIdxs: numpy vector. Contains the indicies of the trials/cases to use
        for evaluating the performance of all the decoders
    pipeName: str. Sprecifies the scikit-learn decoder pipeline to use. For
        options see trainAtOneTimePoint. Options include both classification
        and regression.
    trainTimeIdx: str or int. 'all' to train seperate decoders on data from
        each and every time point. 'inv' to train a single decoder on data
        from many time points. Alternatively, int specifying the index of
        the time point to use.
    timeInvWindow: 2-length tuple of scalars. If using the time invarient
        decoder (i.e. option 'inv' for trainTimeIdx), which time points
        should be used. This input specifies the window. Code uses the input
        'timePoints' to work out which data falls in the window.
    evalPerformance: False | 'atTrainTime' | 'atAllTimes'. Whether or not to
        evaluate the performance of the trained decoder(s), and if so, how.
        If false, None will be returned in place of the usual output df. If
        'atTrainTime' only evaluate the performace of a decoder using data
        from the same time point that was used for training (not avaliable
        if trainTimeIdx='inv'). If 'atAllTimes' evaluate each decoder using
        data from all time points, seperately for each time point.
    splitVar: None | array. If None has no effect. Otherwise shape must match
        the shape of outcome. Gives for each case a category (numeric or
        string) to which that case belongs. Data from different categories
        will be weighted during training to ensure each category is taking
        into account equally during traing, and accuracy data will be
        evaluated seperately for the two categories.

    OUTPUT
    df: data frame containing a single column of performance values (accuracy
        for classifiers, and the coeficient of determination for regressions),
        and a multi-level index. The levels are the time point used
        for training the decoder ('Time_trained'), and the time point used
        for evaluating the decoder ('Time_evaluated'). If splitVar is provided
        then there is an additional index level 'Split_category' because
        the accuracy for different categories in splitVar are evaluated
        seperately.
    decoder: If only a single decoder was trained it is also returned, else
        this output is None
    """
    if len(timePoints) != predictors.shape[2]:
        raise AssertionError('Provided number of time points in timePoints, '+
                             'does not match the number of time points in '+
                             'predictors.')

    if evalPerformance == 'atTrainTime':
        if trainTimeIdx == 'inv':
            raise NotImplementedError('Case is not currently covered.')
        generalise = False

    elif evalPerformance == 'atAllTimes':
        generalise = True
    else:
        assert not evalPerformance

    allDecoders, decTrainTimes = allTimeTrain(predictors, outcome, timePoints,
                                                trainIdxs, pipeName,
                                                trainTimeIdx, timeInvWindow,
                                                splitVar=splitVar)

    if evalPerformance:
        df = allTimeEval(allDecoders, decTrainTimes,
                            predictors, outcome,
                            timePoints, testIdxs,
                            splitVar=splitVar,
                            generalise=generalise)
    else:
        df = None

    if len(allDecoders) == 1:
        decoder = allDecoders[0]
    else:
        decoder = None

    return df, decoder


def allTimeTrain(predictors, outcome, timePoints, trainIdxs, pipeName,
                    trainTimeIdx, timeInvWindow, splitVar=None):
    """ Trains seperate decoders for each time point, or a single
    time-invarient decoder.

    INPUT
    predictors: numpy array. Shape must be
        (num trials or cases, num features, num time points)
    outcome: numpy vector. Contains the values we are aiming to predict. Shape
        should be (num trials or cases,). The values may be strings.
    timePoints: numpy vector. Contains time each time point corresponds to.
        Hence, should have shape (numy time points,)
    trainIdxs: numpy vector. Contains the indicies of the trials/cases to use
        for training all the decoders
    pipeName: str. Sprecifies the scikit-learn decoder pipeline to use. For
        options see trainAtOneTimePoint. Options include both classification
        and regression.
    trainTimeIdx: str or int. 'all' to train seperate decoders on data from
        each and every time point. 'inv' to train a single decoder on data
        from many time points. Alternatively, int specifying the index of
        the time point to use.
    timeInvWindow: 2-length tuple of scalars. If using the time invarient
        decoder (i.e. option 'inv' for trainTimeIdx), which time points
        should be used. This input specifies the window. Code uses the input
        'timePoints' to work out which data falls in the window.
    splitVar: None | array. Shape must match the shape of outcomes. Gives for
        each case a category (numeric or string) to which that case belongs.
        Data from different categories will be weighted during training to
        ensure each category is taken into account equally during traing.

    OUTPUT
    allDecoders: list of trained decoders. If trainTimeIdx is 'inv' or a scalar
        then a 1-element list is returned. Otherwise a list as long as the
        number of time points is returned containing decoders trained on
        data from the corresponding time point.
    decTrainTimes: list as long as allDecoders. Gives the time of the
        data that each decoder was trained on.
    """
    if splitVar is not None:
        assert splitVar.shape == outcome.shape
        assert splitVar.shape == (predictors.shape[0],)
        catNames = np.unique(splitVar)
        print('Training decoder while weighting {} different categories '
                'so that they have a balanced effect.'.format(len(catNames)))
        print('The categories are: {}'.format(catNames))

    allDecoders = []
    decTrainTimes = []
    numTPoints = predictors.shape[2]
    assert numTPoints == len(timePoints)

    if trainTimeIdx == 'inv':
        trainTimeIdx = ['inv']
    elif trainTimeIdx == 'all':
        trainTimeIdx = np.arange(numTPoints)
    elif trainTimeIdx.ndim == 0:
        trainTimeIdx = [trainTimeIdx]
    else:
        raise ValueError('Unexpected input')

    trainPreds = predictors[trainIdxs, :, :]
    trainOutcomes = outcome[trainIdxs]
    if splitVar is not None:
        trainSplitVar = splitVar[trainIdxs]
    else:
        theseSplitVar = None

    tic = time.perf_counter()
    for iTIdx in trainTimeIdx:
        if iTIdx == 'inv':
            incTimePoints = np.logical_and(timePoints >= timeInvWindow[0],
                                            timePoints <= timeInvWindow[1])
            thesePredictors = trainPreds[:, :, incTimePoints]

            condensedPreds = [thesePredictors[:, :, ind] for ind in
                                np.arange(thesePredictors.shape[2])]
            numTPoints = len(condensedPreds)
            thesePredictors = np.concatenate(condensedPreds, axis=0)

            assert thesePredictors.shape[1] == predictors.shape[1]
            assert thesePredictors.shape[0] == (len(trainIdxs)*numTPoints)
            assert thesePredictors.ndim == 2

            print(('Training time invarient decoder with data from {} ' +
                    'time points.').format(numTPoints))

            theseOutcomes = np.tile(trainOutcomes, numTPoints)
            if splitVar is not None:
                theseSplitVar = np.tile(trainSplitVar, numTPoints)
            thisTrainTime = 'Time_invarient'
        else:
            thesePredictors = trainPreds[:, :, iTIdx]
            assert(thesePredictors.shape ==
                        (len(trainIdxs), predictors.shape[1]))
            theseOutcomes = trainOutcomes
            if splitVar is not None:
                theseSplitVar = trainSplitVar
            thisTrainTime = timePoints[iTIdx]

            if np.mod(iTIdx, 20) == 0:
                print('Training decoder on time step {}'.format(iTIdx))
                toc = time.perf_counter()
                if iTIdx > 0:
                    print('Averate train time per time step: {}'.format(
                                                        (toc - tic)/iTIdx))

        thisPipeline = trainAtOneTimePoint(thesePredictors,
                                            theseOutcomes,
                                            pipeName,
                                            splitVar=theseSplitVar)
        allDecoders.append(thisPipeline)
        decTrainTimes.append(thisTrainTime)

    assert len(allDecoders) == len(trainTimeIdx)
    return allDecoders, decTrainTimes


def allTimeEval(allDecoders, decTrainTimes,
                predictors, outcome, timePoints,
                testIdxs, splitVar=None,
                generalise=True):
    """ Evaluate the performance of one or several trained decoders on data
    from all time points (with data from each time point consdered
    individually).

    INPUT
    allDecoders: list of trained decoders.
    decTrainTimes: list | array. As long as allDecoders giving the time
        at which each decoder was trained. Can be strings not scalars.
    predictors: numpy array. Shape must be
        (num trials or cases, num features, num time points)
    outcome: numpy vector. Contains the values we are aiming to predict. Shape
        should be (num trials or cases,)
    timePoints: numpy vector. Contains time each time point corresponds to.
        Hence, should have shape (num time points,)
    testIdxs: numpy vector. Contains the indicies of the trials/cases to use
        for evaluating the performance of all the decoders
    splitVar: None | array. If None has no effect. Otherwise shape must match
        the shape of outcome. Gives for each case a category (numeric or
        string) to which that case belongs. Accuracy data will be evaluated
        seperately for the two categories.
    generalise: bool. If true evaluate each decoder at all time points. If
        false evaluate each decoder only using the time point that was used
        for training that decoder. (Cannot be False if a time invarient
        decoder was trained.)

    OUTPUT
    df: data frame containing a single column of performance values (accuracy
        for classifiers, and the coeficient of determination for regressions),
        and a multi-level index. The levels are the time point used
        for training the decoder ('Time_trained'), and the time point used
        for evaluating the decoder ('Time_evaluated'). If splitVar is provided
        then there is an additional index level 'Split_category' because
        the accuracy for different categories in splitVar are evaluated
        seperately.
    """
    assert predictors.shape[0] == outcome.shape[0]
    if splitVar is not None:
        assert predictors.shape[0] == splitVar.shape[0]

    # Trim down to only the test trials (i.e. remove all trials used for
    # training)
    assert len(predictors.shape) == 3
    predictors = predictors[testIdxs, :, :]
    outcome = outcome[testIdxs]
    if splitVar is not None:
        splitVar = splitVar[testIdxs]

    if splitVar is None:
        df = allTimeEvalSingleSplit(allDecoders, decTrainTimes,
                                    predictors, outcome,
                                    timePoints,
                                    generalise=generalise)
    else:
        assert splitVar.shape == outcome.shape
        wasUsed = np.full(outcome.shape, fill_value=False)
        catNames = np.unique(splitVar)

        allAccDfs = []
        for thisCat in catNames:
            relIdxs = np.isin(splitVar, [thisCat])
            assert not np.any(wasUsed[relIdxs])
            wasUsed[relIdxs] = True
            assert len(predictors.shape) == 3
            assert predictors.shape[0] == outcome.shape[0]
            thesePredictors = predictors[relIdxs, :, :]
            theseOutcomes = outcome[relIdxs]
            assert thesePredictors.shape == (theseOutcomes.shape[0],
                                                *predictors.shape[1:])

            df = allTimeEvalSingleSplit(allDecoders, decTrainTimes,
                                        thesePredictors, theseOutcomes,
                                        timePoints,
                                        generalise=generalise)
            allAccDfs.append(df)
        assert np.all(wasUsed)

        df = pd.concat(allAccDfs, keys=catNames, names=['Split_category'],
                        verify_integrity=True)

    if splitVar is None:
        expectedIndex = ['Time_trained', 'Time_evaluated']
    else:
        expectedIndex = ['Split_category', 'Time_trained', 'Time_evaluated']
    helpKit.checkDfLevels(df, indexLvs=expectedIndex)
    assert (list(df.columns) == ['Accuracy']) or (
        list(df.columns) == ['correlation_coefficient'])
    #assert (list(df.columns) == ['correlation_coefficient'])

    return df


def allTimeEvalSingleSplit(allDecoders, decTrainTimes, testPreds, testOut,
                            timePoints, generalise=True):
    """ Evaluate the performance of one or several trained decoders on all data
    from all time points (time points consdered individually). Evaluates the
    performance on all test data without applying any kind of splits (compare
    to allTimeEval).

    INPUT
    allDecoders: list of trained decoders.
    decTrainTimes: list | array. As long as allDecoders giving the time
        at which each decoder was trained. Entries can be scalars or the
        string 'Time_invarient'.
    testPreds: numpy array. Shape must be
        (num trials or cases, num features, num time points). All provided
        data will be used to evaluate accuracy, therefore if cross-validation
        has been perfomred, no data should be passed to this function that was
        used for training. (I.e. only pass the "test" section of the data
        to this function.)
    testOut: numpy vector. Contains the values we are aiming to predict. Shape
        should be (num trials or cases,). All provided data will be used to
        evaluate accuracy, therefore if cross-validation has been perfomred,
        no data should be passed to this function that was used for training.
        (I.e. only pass the "test" section of the data to this dataset.)
    timePoints: numpy vector. Contains time each time point corresponds to.
        Hence, should have shape (num time points,)
    generalise: bool. If true evaluate each decoder at all time points. If
        false evaluate each decoder only using the time point that was used
        for training that decoder. (Cannot be False if a time invarient
        decoder was trained.)

    OUTPUT
    df: data frame containing a single column of performance values (accuracy
        for classifiers, and the coeficient of determination for regressions),
        and a multi-level index. The levels are the time point used
        for training the decoder ('Time_trained'), and the time point used
        for evaluating the decoder ('Time_evaluated').
    """
    timeTrained = []
    timeEvaluated = []
    score = []

    tic = time.perf_counter()
    for thisDecoder, thisDecTime in zip(allDecoders, decTrainTimes):

        if (thisDecTime != 'Time_invarient') and np.mod(thisDecTime, 20) == 0:
            print('Evaluating decoder that was trained on time step '+
                    '{}'.format(thisDecTime))
            if len(timeTrained) > 0:
                toc = time.perf_counter()
                numDone = len(np.unique(timeTrained))
                print('Average time to evaluate a train time: {}'.format(
                    (toc - tic)/numDone))

        evalTimePoints = np.asarray(timePoints)
        assert evalTimePoints.ndim == 1
        evalTimeIdxs = np.arange(evalTimePoints.shape[0])

        if generalise:
            pass
        elif thisDecTime == 'Time_invarient':
            raise ValueError('For Time_invarient decoders generalise '+
                             'must be set to True.')
        else:
            match = np.isin(evalTimePoints, [thisDecTime])
            evalTimePoints = evalTimePoints[match]
            evalTimeIdxs = evalTimeIdxs[match]

        # Loop through all time points evaluating the performance of this
        # decoder at each one
        for iTestTimeIdx, thisTestTime in zip(evalTimeIdxs, evalTimePoints):
            # Find the relevant data
            theseTestPreds = testPreds[:, :, iTestTimeIdx]
            assert theseTestPreds.shape == testPreds.shape[:2]

            predictions = thisDecoder.predict(theseTestPreds)
            thisScore = np.corrcoef(predictions,testOut)[0, 1]
            timeTrained.append(thisDecTime)
            timeEvaluated.append(thisTestTime)
            score.append(thisScore)

    # scoreName = [findDecoderScore(thisDec) for thisDec in allDecoders]
    # scoreName = np.unique(scoreName)
    # assert len(scoreName) == 1
    # scoreName = scoreName[0]

    scoreName = 'correlation_coefficient'

    df = pd.DataFrame({
        'Time_trained': timeTrained,
        'Time_evaluated': timeEvaluated,
        scoreName: score
    })
    df = df.set_index(['Time_trained', 'Time_evaluated'],
                    verify_integrity=True)
    return df


def trainAtOneTimePoint(predictors, outcomes, pipeName='default',
                        splitVar=None):
    """ Setup a new decoder and train on all the data provided, which should
    only be for a single timepoint.

    INPUT
    pipeName: str. Sprecifies the scikit-learn decoder pipeline to use.
        Options include both classification and regression.
    predictors: numpy array. Shape must be (num trials or cases, num features)
    outcome: numpy vector. Contains the values we are aiming to predict. Shape
        should be (num trials or cases,). The values may be strings.
    splitVar: None | array. Shape must match the shape of outcomes. Gives for
        each case a category (numeric or string) to which that case belongs.
        Data from different categories will be weighted during training to
        ensure each category is taken into account equally during traing.

    OUTPUT
    thisPipeline: Decoding pipeline trained on all cases provided
    """
    if splitVar is not None:
        assert pipeName in ['default', 'pca'], 'Options not compatible'

    thisPipeline = findDecodePipeline(pipeName)

    if splitVar is not None:
        catNames = np.unique(splitVar)
        casesInCat = np.full(catNames.shape, fill_value=np.nan)
        assert len(splitVar.shape) == 1
        for iCat, thisCat in enumerate(catNames):
            casesInCat[iCat] = np.sum(np.isin(splitVar, [thisCat]))

        assert not np.any(np.isnan(casesInCat))
        assert not np.any(casesInCat == 0)
        assert np.sum(casesInCat) == splitVar.shape[0]

        # Set the weights so that the sum of the weights over all cases sums
        # to the number of categories
        normaliser = np.sum(casesInCat) / len(catNames)
        weightsPerCat = (1/casesInCat) * normaliser

        weights = np.full(splitVar.shape, fill_value=np.nan)
        assert len(weights.shape) == 1
        for thisCat, thisWeight in zip(catNames, weightsPerCat):
            relIdxs = np.isin(splitVar, [thisCat])
            assert np.all(np.isnan(weights[relIdxs]))
            weights[relIdxs] = thisWeight

        assert not np.any(np.isnan(weights))
        assert np.sum(weights) == len(outcomes)
    else:
        weights = None

    # Train the decoder
    assert(predictors.ndim == 2)
    assert(outcomes.ndim == 1)
    assert(predictors.shape[0] == outcomes.shape[0])
    if weights is None:
        thisPipeline.fit(predictors, outcomes)
    else:
        assert(weights.shape == outcomes.shape)
        thisPipeline.fit(predictors, outcomes, sample_weight=weights)
    return thisPipeline

def findDecodePipeline(pipeName):
    from sklearn import linear_model
    if pipeName == 'default':
        thisPipeline = Pipeline(steps=[
                            ( "standardise", StandardScaler()),
                            ("PCA", PCA(n_components=0.80, svd_solver='full')),
                            ("RidgeReg", linear_model.Ridge(alpha=1))])


    return thisPipeline
