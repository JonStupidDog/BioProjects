print(__doc__)

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import spline
from sklearn.naive_bayes import GaussianNB
from sklearn.svm import SVC
from sklearn.datasets import load_digits
from sklearn.model_selection import learning_curve
from sklearn.model_selection import ShuffleSplit


def plot_learning_curve(ylim=None, train_sizes=np.linspace(10, 500, 50)):
    """
    Generate a simple plot of the test and training learning curve.

    Parameters
    ----------
    estimator : object type that implements the "fit" and "predict" methods
        An object of that type which is cloned for each validation.

    title : string
        Title for the chart.

    X : array-like, shape (n_samples, n_features)
        Training vector, where n_samples is the number of samples and
        n_features is the number of features.

    y : array-like, shape (n_samples) or (n_samples, n_features), optional
        Target relative to X for classification or regression;
        None for unsupervised learning.

    ylim : tuple, shape (ymin, ymax), optional
        Defines minimum and maximum yvalues plotted.

    cv : int, cross-validation generator or an iterable, optional
        Determines the cross-validation splitting strategy.
        Possible inputs for cv are:
          - None, to use the default 3-fold cross-validation,
          - integer, to specify the number of folds.
          - :term:`CV splitter`,
          - An iterable yielding (train, test) splits as arrays of indices.

        For integer/None inputs, if ``y`` is binary or multiclass,
        :class:`StratifiedKFold` used. If the estimator is not a classifier
        or if ``y`` is neither binary nor multiclass, :class:`KFold` is used.

        Refer :ref:`User Guide <cross_validation>` for the various
        cross-validators that can be used here.

    n_jobs : int or None, optional (default=None)
        Number of jobs to run in parallel.
        ``None`` means 1 unless in a :obj:`joblib.parallel_backend` context.
        ``-1`` means using all processors. See :term:`Glossary <n_jobs>`
        for more details.

    train_sizes : array-like, shape (n_ticks,), dtype float or int
        Relative or absolute numbers of training examples that will be used to
        generate the learning curve. If the dtype is float, it is regarded as a
        fraction of the maximum size of the training set (that is determined
        by the selected validation method), i.e. it has to be within (0, 1].
        Otherwise it is interpreted as absolute sizes of the training sets.
        Note that for classification the number of samples usually have to
        be big enough to contain at least one sample from each class.
        (default: np.linspace(0.1, 1.0, 5))
    """

    plt.figure()
    plt.title("Learning Curves")
    if ylim is not None:
        plt.ylim(*ylim)
    plt.xlabel("Training epochs")
    plt.ylabel("Loss/ACC")

    train_loss = np.zeros((50, 10), dtype=float)
    train_acc = np.zeros((50, 10), dtype=float)
    val_loss = np.zeros((50, 10), dtype=float)
    val_acc = np.zeros((50, 10), dtype=float)
    for k in range(1, 3):
        with open('./data/multi_labels/loss'+str(k)+'.txt', 'r') as f:
            lines = f.readlines()
            for i in range(50):
                tmp = lines[i+1].strip().split()
                train_loss[i][k] = float(tmp[1])
                train_acc[i][k] = float(tmp[2])
                val_loss[i][k] = float(tmp[3])
                val_acc[i][k] = float(tmp[4])
        # with open('./data/multi_labels/loss2.txt', 'r') as f:
        #     lines = f.readlines()
        #     for i in range(50):
        #         tmp = lines[i + 1].strip().split()
        #         train_loss[i][0] = float(tmp[1])
        #         train_acc[i][0] = float(tmp[2])
        #         val_loss[i][1] = float(tmp[3])
        #         val_acc[i][1] = float(tmp[4])

    # train_loss_mean = np.mean(train_loss, axis=1)
    # train_loss_std = np.std(train_loss, axis=1)
    # train_acc_mean = np.mean(train_acc, axis=1)
    # train_acc_std = np.std(train_acc, axis=1)

    val_loss_mean = np.mean(val_loss, axis=1)
    val_loss_std = np.std(val_loss, axis=1)
    val_acc_mean = np.mean(val_acc, axis=1)
    val_acc_std = np.std(val_acc, axis=1)

    plt.grid()

    # plt.fill_between(train_sizes, train_loss_mean - train_loss_std,
    #                  train_loss_mean + train_loss_std, alpha=0.1,
    #                  color="r")
    # plt.fill_between(train_sizes, train_acc_mean - train_acc_std,
    #                  train_acc_mean + train_acc_std, alpha=0.1, color="g")
    plt.fill_between(train_sizes, val_loss_mean - val_loss_std,
                     val_loss_mean + val_loss_std, alpha=0.1,
                     color="b")
    plt.fill_between(train_sizes, val_acc_mean - val_acc_std,
                     val_acc_mean + val_acc_std, alpha=0.1, color="y")
    # plt.plot(train_sizes, train_loss_mean, '-', color="r",
    #          label="train acc")
    # plt.plot(train_sizes, train_acc_mean, '-', color="g",
    #          label="train loss")
    plt.plot(train_sizes, val_loss_mean, '-', color="b",
             label="val acc")
    plt.plot(train_sizes, val_acc_mean, '-', color="y",
             label="val loss")

    plt.legend(loc="best")
    return plt

def loss_curve(ylim=None):
    """
    Generate a simple plot of the test and training learning curve.

    Parameters
    ----------
    estimator : object type that implements the "fit" and "predict" methods
        An object of that type which is cloned for each validation.

    title : string
        Title for the chart.

    X : array-like, shape (n_samples, n_features)
        Training vector, where n_samples is the number of samples and
        n_features is the number of features.

    y : array-like, shape (n_samples) or (n_samples, n_features), optional
        Target relative to X for classification or regression;
        None for unsupervised learning.

    ylim : tuple, shape (ymin, ymax), optional
        Defines minimum and maximum yvalues plotted.

    cv : int, cross-validation generator or an iterable, optional
        Determines the cross-validation splitting strategy.
        Possible inputs for cv are:
          - None, to use the default 3-fold cross-validation,
          - integer, to specify the number of folds.
          - :term:`CV splitter`,
          - An iterable yielding (train, test) splits as arrays of indices.

        For integer/None inputs, if ``y`` is binary or multiclass,
        :class:`StratifiedKFold` used. If the estimator is not a classifier
        or if ``y`` is neither binary nor multiclass, :class:`KFold` is used.

        Refer :ref:`User Guide <cross_validation>` for the various
        cross-validators that can be used here.

    n_jobs : int or None, optional (default=None)
        Number of jobs to run in parallel.
        ``None`` means 1 unless in a :obj:`joblib.parallel_backend` context.
        ``-1`` means using all processors. See :term:`Glossary <n_jobs>`
        for more details.

    train_sizes : array-like, shape (n_ticks,), dtype float or int
        Relative or absolute numbers of training examples that will be used to
        generate the learning curve. If the dtype is float, it is regarded as a
        fraction of the maximum size of the training set (that is determined
        by the selected validation method), i.e. it has to be within (0, 1].
        Otherwise it is interpreted as absolute sizes of the training sets.
        Note that for classification the number of samples usually have to
        be big enough to contain at least one sample from each class.
        (default: np.linspace(0.1, 1.0, 5))
    """

    plt.figure()
    plt.title("Learning Curves")
    if ylim is not None:
        plt.ylim(*ylim)
    plt.xlabel("Training epochs")
    plt.ylabel("Loss/ACC")

    train_sizes = []
    train_loss = []
    val_auc1 = []
    val_auc2 = []
    val_auc3 = []
    with open('result/YK16_3.5_early_stop_with_disorder_wind7_hmm_ss_aa_multi_bilstm_fc_loss.txt', 'r') as f:
        lines = f.readlines()
        for i in range(150):
            if i == 0:
                train_loss.append(2.6)
                val_auc1.append(0.5)
                val_auc2.append(0.5)
                val_auc3.append(0.5)
                train_sizes.append(0.0)
            else:
                train_loss.append(float(lines[i].strip().split()[1]))
                val_auc1.append(float(lines[i].strip().split()[6]))
                val_auc2.append(float(lines[i].strip().split()[7]))
                val_auc3.append(float(lines[i].strip().split()[8]))
                train_sizes.append(float(lines[i].strip().split()[0]))


    # x_new = np.linspace(min(train_sizes), max(train_sizes), 50)
    #
    # train_loss1 = spline(train_sizes, train_loss, x_new)
    # y_smooth1 = spline(train_sizes, val_auc1, x_new)
    # y_smooth2 = spline(train_sizes, val_auc2, x_new)
    # y_smooth3 = spline(train_sizes, val_auc3, x_new)


    plt.grid()
    plt.plot(train_sizes, train_loss, '-', color="r", label="training loss", )
    plt.plot(train_sizes, val_auc1, '-', color="b", label="val auc1")
    plt.plot(train_sizes, val_auc2, '-', color="y", label="val auc2")
    plt.plot(train_sizes, val_auc3, '-', color="gray", label="val auc3")
    plt.legend(loc="best")
    return plt


title = "Learning Curves (Naive Bayes)"
loss_curve()

plt.show()