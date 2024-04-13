import numpy as np
from functools import partial
from itertools import product
from scipy.stats import uniform
from sklearn import svm
from sklearn.model_selection import (
    cross_validate,
    cross_val_predict,
    RandomizedSearchCV,
)
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline

from sklearn.metrics.scorer import make_scorer
from sklearn.metrics import roc_auc_score, mean_squared_error
from sklearn.utils.multiclass import type_of_target
from sklearn.feature_selection import SelectFromModel, SelectKBest

from sklearn.model_selection import StratifiedShuffleSplit
from sklearn.linear_model import LinearRegression

from sklearn.decomposition import PCA
from sklearn.utils import shuffle
import pandas as pd


class Decoder(object):
    def __init__(self, target, classifier=None):
        """

        Args:
            classifier: (name, Sklearn classifier object)
                A tuple that contains a name at first position and a
                sklearn classifier at second position.
            target: target values (to be predcited variable)

        """
        self.classifier = classifier
        self.target = target

        if classifier is None:
            self.clf_name = "SVClin"
            self.clf = None
        else:
            self.clf_name, self.clf = classifier


    def classify(self, data, time):
            """Perform decoding on source reconstructed data.

            Decoding is carried out for each time point separately.

            Args:
                data: ndarray
                    3d: ntrials x channels (X and Y gaze position) x time
                time: ndarray
                    time points that match last dimension of data
            Returns:
                A pandas DataFrame that contains decoding accuracies
            """
            n_trials = data.shape[0]

            target_vals = self.target    # PM: original was target.loc[trial]
            idnan = np.isnan(target_vals)
            scores = []
            print(data.shape)

            for idx_t, tp in enumerate(time):
                X = data[:, :, idx_t]
                print('Time:', tp, 'Size:', X.shape)
                clf = self.clf
                s = categorize(clf, target_vals, X)
                s["latency"] = tp

                scores.append(s)
            return pd.DataFrame(scores)



def categorize(clf, target, data, njobs=1):
    """
    Expects a pandas series and a pandas data frame.
    Both need to be indexed with the same index.
    """
    from sklearn.metrics import make_scorer

    corr_scorer = make_scorer(lambda x, y: np.corrcoef(x, y)[0, 1])
    metrics = {"correlation":corr_scorer}

    score = cross_validate(
        clf, data, target, cv=10, scoring=metrics, return_train_score=False, n_jobs=njobs
    )
    del score["fit_time"]
    del score["score_time"]
    score = {k: np.mean(v) for k, v in list(score.items())}
    print(score)

    return score
