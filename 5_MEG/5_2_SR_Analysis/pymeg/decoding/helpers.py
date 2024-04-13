import numpy as np


def checkDfLevels(df, indexLvs=None, colLvs=None, ignoreOrder=False):
    """ Check that the the levels of the index, or levels of the columns, of a
    pandas dataframe, match those that are expected.

    INPUT
    df: Dataframe to check. Can also pass a series as long as colLvs is None.
    indexLvs: list. Expected levels of the index. If None, is not checked
    colLvs: list. Expected levels of the columns. If None, is not checked
        ignoreOrder: boolean. If true, ignore the order of the index and column
        levels
    """
    levelInfo = {
        'exptIdxLvs': indexLvs,
        'realIdxLvs': df.index.names
    }
    if colLvs is not None:
        colLevelInfo = {
            'exptColLvs': colLvs,
            'realColLvs': df.columns.names
        }
        levelInfo.update(colLevelInfo)

    for key, val in levelInfo.items():
        levelInfo[key] = operateIfNotNone(np.asarray, val)

        if ignoreOrder:
            levelInfo[key] = operateIfNotNone(np.sort, val)

    if indexLvs != None:
        if not np.array_equal(levelInfo['realIdxLvs'],
                                levelInfo['exptIdxLvs']):
            raise mkDfLvsException('index', levelInfo['exptIdxLvs'],
                                levelInfo['realIdxLvs'])
    if colLvs != None:
        if not np.array_equal(levelInfo['realColLvs'],
                                levelInfo['exptColLvs']):
            raise mkDfLvsException('columns', levelInfo['exptColLvs'],
                                levelInfo['realColLvs'])

    if (indexLvs == None) and (colLvs == None):
        raise Exception("No check was requested")


def operateIfNotNone(operation, thisInput):
    if thisInput is None:
        return thisInput
    else:
        return operation(thisInput)
