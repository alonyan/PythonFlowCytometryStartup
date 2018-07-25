import numpy as np
from scipy.optimize import fmin

def fit_hill(x, y, variance=None):
    '''
    Input:
    x = numpy array of independent variable
    y = numpy array of dependent variable
    variance = numpy array of the variance, or standard error, of each data point y

    Returns
    numpy array of maximum likelihood estimates of the parameters
    '''
    theta_mle = fmin(obj,
                     get_initial_params(x, y),
                     args=(x, y, variance))
    return theta_mle


def hill_model(theta, x):
    '''
    Input:
    Model Parameters: theta = python list, or numpy array
    theta[0] = amplitude
    theta[1] = hill coefficient
    theta[2] = half effective concentration
    theta[3] = background

    Returns:
    Numpy array, model evaluations
    '''
    return theta[0] * x**theta[1] / (x**theta[1] + theta[2]**theta[1]) + theta[3]


# ========================================================
# ========================================================

def obj(theta, x, y, var):
    '''
    Find the maximum likelihood solution of parameters of the hill model given the data.  Hard constraint for the amplitude to be larger than the background.
    Input:
    theta: Model parameters
    x : numpy array of x data
    y : numpy array of y data
    var : numpy array of y data variance

    Returns:
    Numpy array of model parameters
    '''
    penalty = 0
    if theta[0] < theta[3]:
        penalty = np.inf
    if var is None:
        chi = np.sum(0.5 * (y - hill_model(theta, x))**2) + penalty
    else:
        chi = np.sum(0.5 * (y - hill_model(theta, x))**2 / var) + penalty
    return chi


# ========================================================
# ========================================================


def get_initial_params(x, y, percentile=0.25):
    '''
    Make / get initial estimates of parameters

    Input:
    x: numpy array of x data
    y: numpy array of y data
    percentile : Thresholds for finding the dynamic range, see 'get_dyanmic_range' function for more details.

    Returns:
    Python list of parameters
    '''
    temp = [y.max() - y.min(),
            get_hill_coef(x, y, percentile),
            get_sensitivity(x, y, percentile),
            y.min()]
    return temp

# ========================================================
# ========================================================

def get_hill_coef(x, y, p):
    '''
    Estimate the hill coefficient by computing the d log(f) / d log(x)

    Input:
    x : numpy array of x data
    y : numpy array of y data
    p : scalar, the percentile used to find the dynamic range

    Returns:
    scalar estimate of the hill coefficient
    '''
    idx = get_dynamic_range(x, y, percentile=p)
    data_idxs = range(idx[0], idx[1]+1)
    y_temp = y[data_idxs]
    x_temp = x[data_idxs][y_temp != 0]
    y_temp = y_temp[y_temp != 0]
    cov_matrix = np.cov([np.log(x_temp),
                         np.log(y_temp)])
    return 2 * cov_matrix[0, 1] / cov_matrix[0, 0]


# ========================================================
# ========================================================



def get_sensitivity(x, y, p):
    '''
    Estimate the EC50 or IC50 by computing the intersection of a horizontal line with at y = 0.5 * amplitude + background and the linear expansion of the hill function about the EC50.  The linear expansion is estimated by fitting a line to the 'dynamic' range of the hill function.

    Input:
    x : numpy array of x data
    y : numpy array of y data
    p : scalar of percentile, same as in hill coef estimate

    Return:
    scalar estimate of the EC50 or IC50
    '''
    idx = get_dynamic_range(x, y, percentile=p)
    data_idxs = range(idx[0], idx[1] + 1)
    cov_matrix = np.cov([x[data_idxs],
                         y[data_idxs]])
    m2 = cov_matrix[0, 1] / cov_matrix[0, 0]
    b2 = np.mean(y[data_idxs]) - m2*np.mean(x[data_idxs])
    b1 = 0.5*(y.max() - y.min()) + y.min()
    return (b1 - b2) / m2

# ========================================================
# ========================================================

def get_dynamic_range(x, y, percentile=0.25):
    '''
    Find the points that make up the dynamic range of the hill function.  Find this by finding all points y_i such that,
    percentile*amplitude + b < y_i < (1 - percentile) * amplitude + b,
    where b is background.

    Input:
    x : numpy array of x data
    y : numpy array of y data

    Return :
    python list of indexes that bound the dynamic range.
    '''
    # initialize
    b = y.min()
    amp = y.max() - b
    # if data is in ascending order
    step = 1
    # if data is in descendeing order
    if y.argmax() < y.argmin():
        step = -1
    # find upper bound
    count = y.argmax()
    while y[count]-b > (1 - percentile)*amp:
        count -= step
    max_idx = count + step
    # find lower bound
    count = y.argmin()
    while y[count]-b < percentile*amp:
        count += step
    min_idx = count - step
    return_idx = [min_idx, max_idx]
    if max_idx < min_idx:
        return_idx = [max_idx, min_idx]
    return return_idx
