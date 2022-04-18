# -*- coding: utf-8 -*-
"""
Created on Wed Jul 29 16:24:07 2020

We aim to predict the buzz, using supervised machine learning, 
based the Random Forest and Logistic Regression models based on 27 features, 
feature engineered from accelerometer data and depth data,
then make several plots to compare the quality between the predictions and ground truths.

"""

#%% Functions for Average of Peak Frequency
# https://reader.elsevier.com/reader/sd/pii/S1877050914008643
# https://ieeexplore-ieee-org.ep.fjernadgang.kb.dk/document/5634532

import numpy as np
from numba import jit, njit, float64


#    Find local maxima in a 1D array (from SciPy).
#    "This function finds all local maxima in a 1D array and returns the indices
#    for their edges and midpoints (rounded down for even plateau sizes)"".
# https://github.com/scipy/scipy/blob/master/scipy/signal/_peak_finding_utils.pyx#L19
@njit
def _local_maxima_1d(x):
    """
    Find local maxima in a 1D array.
    This function finds all local maxima in a 1D array and returns the indices
    for their edges and midpoints (rounded down for even plateau sizes).
    Parameters
    ----------
    x : ndarray
        The array to search for local maxima.
    Returns
    -------
    midpoints : ndarray
        Indices of midpoints of local maxima in `x`.
    left_edges : ndarray
        Indices of edges to the left of local maxima in `x`.
    right_edges : ndarray
        Indices of edges to the right of local maxima in `x`.
    Notes
    -----
    - Compared to `argrelmax` this function is significantly faster and can
      detect maxima that are more than one sample wide. However this comes at
      the cost of being only applicable to 1D arrays.
    - A maxima is defined as one or more samples of equal value that are
      surrounded on both sides by at least one smaller sample.
    .. versionadded:: 1.1.0
    """
    # cdef:
    #     np.intp_t[::1] midpoints, left_edges, right_edges
    #     np.intp_t m, i, i_ahead, i_max

    # Preallocate, there can't be more maxima than half the size of `x`
    midpoints = np.empty(x.shape[0] // 2, dtype=np.intp)
    # left_edges = np.empty(x.shape[0] // 2, dtype=np.intp)
    # right_edges = np.empty(x.shape[0] // 2, dtype=np.intp)
    m = 0  # Pointer to the end of valid area in allocated arrays

    # with nogil:
    i = 1  # Pointer to current sample, first one can't be maxima
    i_max = x.shape[0] - 1  # Last sample can't be maxima
    while i < i_max:
        # Test if previous sample is smaller
        if x[i - 1] < x[i]:
            i_ahead = i + 1  # Index to look ahead of current sample

            # Find next sample that is unequal to x[i]
            while i_ahead < i_max and x[i_ahead] == x[i]:
                i_ahead += 1

            # Maxima is found if next unequal sample is smaller than x[i]
            # if x[i_ahead] < x[i] and (i + 2 <= i_ahead) :
            if x[i_ahead] < x[i] :
                # left_edges[m] = i
                # right_edges[m] = i_ahead - 1
                # midpoints[m] = (left_edges[m] + right_edges[m]) // 2
                midpoints[m] = (i + i_ahead - 1) // 2
                m += 1
                # Skip samples that can't be maximum
                i = i_ahead
        i += 1

    # # Keep only valid part of array memory.
    # midpoints.base.resize(m, refcheck=False)
    # left_edges.base.resize(m, refcheck=False)
    # right_edges.base.resize(m, refcheck=False)

    # return midpoints.base, left_edges.base, right_edges.base
    return midpoints[:m]


# Find peaks inside a signal based on peak properties (SciPy).
# "This function takes a 1-D array and finds all local maxima by simple comparison of neighboring values. 
# Optionally, a subset of these peaks can be selected by specifying conditions for a peak’s properties".
# https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.find_peaks.html
@njit
def find_peaks_scipy(x, has_height:int=0, hmin:np.float64=0.0, hmax:np.float64=0.0) -> np.ndarray:
  
    # _argmaxima1d expects array of dtype 'float64'
    x = x.astype(np.float64)
    if x.ndim != 1:
        raise ValueError('`x` must be a 1-D array')

    # peaks, left_edges, right_edges = _local_maxima_1d(x)
    peaks = _local_maxima_1d(x)

    if has_height :
        # Evaluate height condition
        peak_heights = x[peaks]
        # hmin, hmax = _unpack_condition_args(height, x, peaks)

        # keep = _select_by_property(peak_heights, hmin, hmax)
        # https://stackoverflow.com/a/21174962
        keep = np.full(peak_heights.size, True)
        if has_height == 1:
            keep &= (hmin <= peak_heights)
        if has_height == 2:
            keep &= (peak_heights <= hmax)

        peaks = peaks[keep]
        # peaks_l = peaks_l[keep]
        # peaks_r = peaks_r[keep]

    return peaks

@njit
def peak_freq(peak): 
    return len(peak)

#%% Function 'feature_extract' 

# Features engineering for Random Forest and Logistic Regresion, based on the article
# 

@jit(float64[:,:](float64[:, :], float64, float64, float64), nopython=True)
def feature_extract(A:np.ndarray, rate:float64, n_sec:float64, overlap:float64) -> np.ndarray :

  Ax = A[:,0]
  Ay = A[:,1]
  Az = A[:,2]
  depth = A[:,3]
  
  divestate = A[:,4]
  onehot = np.zeros(4)

  B = A[:,-1]
  W_l = rate*n_sec
  n = A.shape[0]
  
  df = np.full( (round(A.shape[0]//(W_l*(1-overlap))), 24), np.inf, dtype=np.float64)
  
  r = i = n
  r = i = 0
  while i+W_l <= n:
    x = Ax[i : i+W_l]
    y = Ay[i : i+W_l]
    z = Az[i : i+W_l]
    d = depth[i : i+W_l]
    
    state = divestate[i : i+W_l]
    onehot[0] = onehot[1] = onehot[2] = onehot[3] = 0
    if np.mean(state) > (state[0] + state[-1])/2 :
      onehot[round(state[-1])] = 1
    else :
      onehot[round(state[0])] = 1
    
    
    buzz = B[i : i+W_l]
    a = np.sqrt(x**2 + y**2 + z**2)
    
    x_m = np.mean(x)
    y_m = np.mean(y)
    z_m = np.mean(z)
    
    # a_m = np.mean(a)
    a_minmax = max(a) - min(a)
    a_std = np.std(a)
    a_rms = np.sqrt( np.mean(a**2) )
    
    x_std = np.std(x)
    y_std = np.std(y)
    z_std = np.std(z)
  
    x_rms = np.sqrt( np.mean(x**2) )
    y_rms = np.sqrt( np.mean(y**2) )
    z_rms = np.sqrt( np.mean(z**2) )
    
    # https://stackoverflow.com/a/39098306
    r_xy = np.cov(x, y, bias=True)[0][1]/(x_std*y_std)
    r_yz = np.cov(z, y, bias=True)[0][1]/(y_std*z_std)
    r_zx = np.cov(x, z, bias=True)[0][1]/(z_std*x_std)
  
    x_minmax = max(x) - min(x)
    y_minmax = max(y) - min(y)
    z_minmax = max(z) - min(z)
    
    depth_mean = np.mean(d)
    
    
    C = 0.0
    if buzz.sum() >= W_l/2 :
      C = 1.0
    
    # https://stackoverflow.com/a/33511352
    L = [ x_m, y_m, z_m, a_minmax, a_std, a_rms, x_std, y_std, z_std, x_rms, y_rms, z_rms,
            r_xy, r_yz, r_zx, x_minmax, y_minmax, z_minmax, depth_mean, 
            onehot[0], onehot[1], onehot[2], onehot[3], C ]
    # L = [ x_m, y_m, z_m, a_minmax, a_std, a_rms, x_std, y_std, z_std, x_rms, y_rms, z_rms,
    #        r_xy, r_yz, r_zx, x_minmax, y_minmax, z_minmax, C ]
    df[r,:] = np.array(L)
    
    # # sma = np.mean(abs(x - x_m)) + np.mean(abs(y - y_m)) + np.mean(abs(z - z_m))
    # # H = sum( abs(a-a_m) * np.log10(abs(a-a_m)) )
    # # sk = sum( (a-a_m)**3 ) / (W_l*np.std(a)**3)
    # # kurtosis = sum( (a-a_m)**4 ) / (W_l*np.std(a)**4)
    # # E_a = abs(np.mean(scipy.fft.fft(a)))
    
    i += round(W_l*(1-overlap))
    r += 1
    
  return df


#%% Peak extract function

@jit(float64[:,:](float64[:, :], float64, float64, float64), nopython=True)
def peak_extract(A:np.ndarray, rate:float64, n_sec:float64, overlap:float64) -> np.ndarray :

  n = A.shape[0]
  Ax = A[:,0]
  Ay = A[:,1]
  Az = A[:,2]
  
  Px = find_peaks_scipy(Ax)
  Py = find_peaks_scipy(Ay)
  Pz = find_peaks_scipy(Az)

  r = i = n
  r = i = 0
  W_l = round(rate*n_sec)
  
  
  df = np.full( (round(A.shape[0]//(W_l*(1-overlap))), 7) , np.inf, dtype=np.float64)
  
  while i+W_l <= n:    
    # find Px s.t. i <= Px <= i + W_l
    p = np.searchsorted( Px, np.array([i,i+W_l-1]) )
    if p[0] == p[1] :
      x_apf = 0
    else :
      peak_x = Px[ p[0] : p[1] ]        
      x_apf = peak_freq(peak_x)
      if len(peak_x) > 1 :
        elapse_time_x = np.diff(peak_x).mean()
      else :
        elapse_time_x = 0
    # find Py s.t. i <= Py <= i + W_l
    p = np.searchsorted( Py, np.array([i,i+W_l-1]) )
    if p[0] == p[1] :
      y_apf = 0
    else :
      peak_y = Py[ p[0] : p[1] ]
      y_apf = peak_freq(peak_y)
      if len(peak_y) > 1 :
        elapse_time_y = np.diff(peak_y).mean()
      else :
        elapse_time_y = 0
    # find Pz s.t. i <= Pz <= i + W_l
    p = np.searchsorted( Pz, np.array([i,i+W_l-1]) )
    if p[0] == p[1] :
      z_apf = 0
    else :
      peak_z = Pz[ p[0] : p[1] ]
      z_apf = peak_freq(peak_z)
      if len(peak_z) > 1 :
        elapse_time_z = np.diff(peak_z).mean()
      else :
        elapse_time_z = 0
    
    apf_mean = (x_apf+y_apf+z_apf)/3.0
    VarAPF = ( (x_apf-apf_mean)**2 + (y_apf-apf_mean)**2 + (z_apf-apf_mean)**2 )/2.0
  
    L = [ x_apf, y_apf, z_apf, VarAPF, elapse_time_x, elapse_time_y, elapse_time_z ]
    df[r,:] = np.array(L)
    
    i += round(W_l*(1-overlap))
    r += 1
    
  print('r = ', r)
  return df


#%% Make datasets to train and test for Random Forest and Logistic Regression

import pandas as pd
pd.set_option('display.max_columns', None)

whale = ['Asgeir','Helge18','Kyrri','Nemo','Siggi']
Acc = [None]*len(whale)

rate = 100
# The ratio for training/test set is 80%:20%
train_ratio = 0.8 # no validation set used since we do not have hyperparmeter selection
test_ratio = 0.8

X_train = [None]*len(whale)
X_test = [None]*len(whale)
y_train = [None]*len(whale)
y_test = [None]*len(whale)

for i in range(len(whale)) :
    Acc[i] = pd.read_parquet('data/accel-' + whale[i] + '.csv.parquet')
    accel = Acc[i]
    
    A = accel.iloc[:,[1,2,3,4,9,-1]].to_numpy()
    n = A.shape[0]
    
    a_train = feature_extract(A[:round(train_ratio*n)], rate, 1, 0.5)
    a_train = a_train[np.isfinite(a_train).all(axis=1)]
    
    a_test = feature_extract(A[round(test_ratio*n):], rate, 1, 0)
    a_test = a_test[np.isfinite(a_test).all(axis=1)]
  
    b_train = peak_extract(A[:round(train_ratio*n),:3], rate, 1, 0.5)
    b_train = b_train[np.isfinite(b_train).all(axis=1)]
    
    b_test = peak_extract(A[round(test_ratio*n):,:3], rate, 1, 0)
    b_test = b_test[np.isfinite(b_test).all(axis=1)]
  
    X_train[i] = np.c_[a_train[:,:-1], b_train]
    X_test[i]  = np.c_[a_test[:,:-1], b_test]
    
    y_train[i] = a_train[:,-1]
    y_test[i]  = a_test[:,-1]


#%% some functions

# Show the range of array x
def _range(x):
  return [min(x), max(x)]

# Divide sequence y to a list of continous sequences
def consecutive(y, thresh=0.0, stepsize=1): # https://stackoverflow.com/a/7353335
    data = np.array(np.where(y > thresh ))[0]
    return np.split(data, np.where(np.diff(data) != stepsize)[0]+1)

from numba import njit

# Calculate the intersections of sequences begin with s0, end with e0
# with sequences begin with s1, end with e1
@njit
def intersect(s0,e0, s1,e1) :
  m = len(s0)
  n = len(s1)
  if n == 0 :
    return None, None
  d = np.zeros( (m, n), dtype=np.int64)
  ovl = np.zeros( (m, n), dtype=np.float32)
  i = j = 0
  
  for j in range(m) :
    for i in range(n) :
      # https://scicomp.stackexchange.com/a/26260
      if s0[j] > e1[i] or s1[i] > e0[j] : # (1) & (2)
        if s0[j] > e1[i] : # (1)
          d[j][i] = e1[i] - s0[j]
        if s1[i] > e0[j] : # (2)
          d[j][i] = s1[i] - e0[j]
      else: # overlap
        s_ovl = max(s0[j], s1[i])
        e_ovl = min(e0[j], e1[i])
        d[j][i] = 0
        ovl[j][i] = (e_ovl-s_ovl)/(e0[j]-s0[j])*100
  return d, ovl

# Compute if y and y_t overlap, if not how far between them
def hit_rate(y, y_t, dist=100, ovl_ptg=50, accel=None) :
    C = consecutive(y, 0.5)
    C_t = consecutive(y_t)
    if len(C[0]) == 0 :
        return 0, 0, None, None
    s1 = np.ravel([ x[0] for x in C ])
    e1 = np.ravel([ x[-1] for x in C ]) + 1
    s0 = np.ravel([ x[0] for x in C_t ])
    e0 = np.ravel([ x[-1] for x in C_t ]) + 1

    a, ovl = intersect(s0,e0, s1,e1)
    hit = np.where( abs(a).min(axis=1) <= dist )[0]
    missed = np.sort(np.array( list(set(range(len(s0)))-set(hit)) ))
    hit_ovl = np.where( ovl.max(axis=1) > ovl_ptg )[0]

    return len(hit)/len(s0), len(hit_ovl)/len(s0), hit, missed



#%% Train with Random Forest
# [https://mljar.com/blog/random-forest-overfitting/]
# [https://stackoverflow.com/a/35012011]
# [https://datascience.stackexchange.com/a/6430]
# [https://stats.stackexchange.com/a/376785]

x_tr = np.vstack(X_train)
y_tr = np.concatenate(y_train)

# Select number of tree for Random Forest
# n_tree = 500
# n_tree = 1000
n_tree = 2000

# Fit the model with sklearn
from sklearn.ensemble import RandomForestClassifier
clf = RandomForestClassifier(n_estimators = n_tree, random_state=0, 
                             class_weight='balanced_subsample', n_jobs=16)
clf.fit( x_tr, y_tr )


#%% Save/Load model [https://stackoverflow.com/a/33500427]

import joblib

# Save the fitted model
with open('data/all_OneHot_' + str(n_tree) + '.rf' , 'wb') as f:
    joblib.dump(clf, f, 5)
# Load the fitted model
with open('data/all_OneHot_' + str(n_tree) + '.rf', 'rb') as f:
    clf = joblib.load(f)


#%% Confusion matrix of RF for each whale

from sklearn.metrics import confusion_matrix

for i in range(len(whale)) :
    y_pred = clf.predict(X_test[i])
    print ( i, ':', confusion_matrix(y_test[i], y_pred).ravel() )

# (tn, fp, fn, tp):
# 0 : [45654    16   234    11]
# 1 : [45357    12  1463    48]
# 2 : [54876     9   756    32]
# 3 : [13653     1    64     0]
# 4 : [13052     0   124     1]


# 0 : [22888     0    69     0]
# 1 : [23242     0   197     1]
# 2 : [27721     1   111     3]
# 3 : [6800    1   58    0]
# 4 : [6511    0   76    1]


# 500 trees
# 0 : [22888     0    69     0]
# 1 : [23242     0   197     1]
# 2 : [27718     4   103    11]
# 3 : [6799    2   58    0]
# 4 : [6509    2   76    1]

# 1000 trees
# 0 : [22888     0    69     0]
# 1 : [23242     0   197     1]
# 2 : [27718     4   102    12]
# 3 : [6799    2   58    0]
# 4 : [6509    2   76    1]

# 2000 trees
# 0 : [22888     0    69     0]
# 1 : [23242     0   197     1]
# 2 : [27718     4   100    14]
# 3 : [6799    2   58    0]
# 4 : [6509    2   76    1]


#%% Count predicted buzz and ground true buzz per dive
# Also calculate the sum of length of predicted buzz vs length of ground true buzz,
# and write to a dataframe

df = [None]*len(whale)
for i in range(len(whale)) :
  y_pred = clf.predict(X_test[i])
  C = consecutive(y_pred)
  C_t = consecutive(y_test[i])
  
  if len(C[0]) :
    s1 = np.ravel([ x[0] for x in C ])
    e1 = np.ravel([ x[-1] for x in C ]) + 1
    s1 *= rate; e1 *= rate
  if len(C_t[0]) :
    s0 = np.ravel([ x[0] for x in C_t ])
    e0 = np.ravel([ x[-1] for x in C_t ]) + 1
    s0 *= rate; e0 *= rate
  
  dt = Acc[i].iloc[round( test_ratio*Acc[i].shape[0] ):,:] # [1371800 rows x 13 columns]
  
  dif_dive = np.diff(dt['Dive_no'].to_numpy())
  dive_st = np.where(dif_dive > 0)[0] + 1
  dive_end = np.where(dif_dive < 0)[0] + 1
  if len(dive_st) < len(dive_end) :
    dive_st = np.r_[0, dive_st]
  if len(dive_st) > len(dive_end) :
    dive_end = np.r_[dive_end, dt.shape[0]]
  
  
  buzz_dive = np.zeros( (len(dive_st), 8) )
  for k in range(len(dive_st)) :
    depth = dt.iloc[dive_st[k]:dive_end[k],4].to_numpy()
    buzz_dive[k, 0] = i
    buzz_dive[k, 1] = k
    buzz_dive[k, 2] = depth.max()
    buzz_dive[k, 3] = depth.shape[0]/100
    
    if len(C_t[0]) :
      buzz_dive_0 = np.where( (s0 >= dive_st[k]) & (s0 < dive_end[k]) )[0]
      buzz_dive[k, 4] = len(buzz_dive_0)
      buzz_dive[k, 6] = (e0[buzz_dive_0] - s0[buzz_dive_0]).sum()
    if len(C[0]) :
      buzz_dive_1 = np.where( (s1 >= dive_st[k]) & (s1 < dive_end[k]) )[0]
      buzz_dive[k, 5] = len(buzz_dive_1)
      buzz_dive[k, 7] = (e1[buzz_dive_1] - s1[buzz_dive_1]).sum()
  
  df[i] = pd.DataFrame(data = buzz_dive , 
                       columns = ['ID','Dive_no','MaxDepth','Duration',
                                  'True','Prediction','Len_true','Len_pred'] )
  df[i]['Name'] = pd.Series(np.repeat(whale[i], len(dive_st)))


Dive_res = pd.concat(df)
from collections import Counter
# https://stackoverflow.com/a/57739795
Dive_res['hash'] = pd.util.hash_pandas_object(Dive_res.loc[:,['True','Prediction']], 
                                              index=False)
Dive_res['Num_point'] = Dive_res['hash'].map(dict(Counter(Dive_res['hash'].to_numpy())))
Dive_res['Diff'] = Dive_res['Prediction'].to_numpy() - Dive_res['True'].to_numpy()
Dive_res['Diff_len'] = Dive_res['Len_pred'].to_numpy() - Dive_res['Len_true'].to_numpy()
Dive_res.to_csv( 'data/Dive_res_RF.csv', index=False )


#%% For each whale

Res = np.zeros( (6,9))
for i in range(len(whale)) :
  y_pred = clf.predict(X_test[i])
  L0 = []
  for p in [99, 50, 25, 0] :
    _, ptg, _, _ = hit_rate(y_pred, y_test[i], ovl_ptg=p)
    if ptg is not None :
      L0.append( round(ptg, 2) )
    else :
      L0.append( None )
  for t in [1, 2, 3, 4, 5] :
    R, _, hit, missed = hit_rate(y_pred, y_test[i], t)
    if R is not None :
      L0.append( round(R, 2) )
    else :
      L0.append( None )
  Res[i,:] = np.ravel(L0)

pd.DataFrame(data=Res.T, columns = whale,
             index=np.ravel(['> 99%', '> 50%', '> 25%', '> 0%', 
                             '1 s', '2 s', '3 s', '4 s', '5 s']) ).to_html()

#%% We calcuate for all whales the overlap/distance between predicted buzzes and ground true buzzes

y_pred_all = clf.predict(np.vstack(X_test))
y_test_all = np.hstack(y_test)

L0 = []
for p in [99, 50, 25, 0] :
  _, ptg, _, _ = hit_rate(y_pred_all, y_test_all, ovl_ptg=p)
  if ptg is not None :
    L0.append( round(ptg, 2) )
  else :
    L0.append( None )
for t in [1, 2, 3, 4, 5] :
  R, _, hit, missed = hit_rate(y_pred_all, y_test_all, t)
  if R is not None :
    L0.append( round(R, 2) )
  else :
    L0.append( None )

Res[-1,:] = np.ravel(L0)
res = pd.DataFrame(data= {'index': np.tile(np.ravel(['> 99%', '> 50%', '> 25%', '> 0%', 
                                           '1 s', '2 s', '3 s', '4 s', '5 s']), Res.shape[0] ) ,
                        'val': Res.ravel() , 
                        'whale': np.repeat(whale + ['All'],Res.shape[1])  })
res['method'] = 'Random Forest'
# Save the result to file 
res.to_csv( 'B:/Codes/Data/Fieldwork 2018/res_RF.csv', index=False )

# Then we can plot Fig.7 with plotnine
from plotnine import *
ggplot(res) + geom_line(aes(x='index',y='val',color='whale', group=1))


#%% Plot Fig. 11 (Random Forest)

# https://matplotlib.org/3.1.0/api/_as_gen/matplotlib.pyplot.figure.html
import matplotlib
import matplotlib.pyplot as plt
matplotlib.rcParams["figure.dpi"] = 300
matplotlib.rcParams["figure.figsize"] = [11, 5]
matplotlib.rcParams['svg.fonttype'] = 'none'
plt.style.use('ggplot')
import matplotlib.patches as patches

# Select 2nd whale (ID 20158) and its 6th dive to plot the prediction and ground true buzz on the same dive
i = 2; k = 6

y_pred = clf.predict(X_test[i])
C = consecutive(y_pred)
C_t = consecutive(y_test[i])

if len(C[0]) :
  s1 = np.ravel([ x[0] for x in C ])
  e1 = np.ravel([ x[-1] for x in C ]) + 1
  s1 *= rate; e1 *= rate
if len(C_t[0]) :
  s0 = np.ravel([ x[0] for x in C_t ])
  e0 = np.ravel([ x[-1] for x in C_t ]) + 1
  s0 *= rate; e0 *= rate

dt = Acc[i].iloc[round( test_ratio*Acc[i].shape[0] ):,:].copy() # [1371800 rows x 13 columns]

dif_dive = np.diff(dt['Dive_no'].to_numpy())
dive_st = np.where(dif_dive > 0)[0] + 1
dive_end = np.where(dif_dive < 0)[0] + 1
if len(dive_st) < len(dive_end) :
  dive_st = np.r_[0, dive_st]
if len(dive_st) > len(dive_end) :
  dive_end = np.r_[dive_end, dt.shape[0]]


depth = dt.iloc[dive_st[k]:dive_end[k],4].to_numpy()

if len(C_t[0]) :
  buzz_dive_0 = np.where( (s0 >= dive_st[k]) & (s0 < dive_end[k]) )[0]
if len(C[0]) :
  buzz_dive_1 = np.where( (s1 >= dive_st[k]) & (s1 < dive_end[k]) )[0]
  
  
fig, ax = plt.subplots()
ax.tick_params(axis='both', which='major', labelsize=15)
plt.gca().invert_yaxis()
X_axis = np.arange(dive_st[k], dive_end[k])[::10]
tick_locs = np.arange(X_axis[0], X_axis[-1], 2*60*rate)
tick_lbls = ((tick_locs-X_axis[0])/(60*rate)).round().astype(int)
plt.xticks(tick_locs.tolist(), [])
ax.tick_params(axis='x', length=0, width=2)


plt.plot(np.arange(dive_st[k], dive_end[k])[::10], 
         dt.iloc[dive_st[k]:dive_end[k],4].to_numpy()[::10], 
        label='_nolegend_', c='k', alpha=0.25, linewidth=6)
for j in buzz_dive_0 :
  plt.plot( np.arange(s0[j], e0[j])[::10], 
           dt.iloc[ s0[j]:e0[j], 4 ].to_numpy()[::10], c='#D55E00', linewidth=6)
for j in buzz_dive_1 :
  plt.plot( np.arange(s1[j], e1[j])[::10], 
           dt.iloc[ s1[j]:e1[j], 4 ].to_numpy()[::10], c='#0072B2', alpha=0.5, linewidth=6)

offset1 = 6000
offset2 = 24000
st = min(s0[min(buzz_dive_0)],s1[min(buzz_dive_1)])
en = max(e0[max(buzz_dive_0)],e1[max(buzz_dive_1)])
st += offset1
en -= offset2

plt.plot( [st, en], [440, 440], linestyle=':', color='black')
plt.plot( [st, en], [530, 530], linestyle=':', color='black')
plt.plot( [st, st], [440, 530], linestyle=':', color='black')
plt.plot( [en, en], [440, 530], linestyle=':', color='black')

h = 0.62
t = 0.17
plt.plot( [st, int( (1-t)*dive_st[k] + t*dive_end[k] )], 
         [440, int(h*max(dt.iloc[dive_st[k]:dive_end[k],4].to_numpy()))], 
         linestyle=':', color='black')
t2 = t + 0.64
plt.plot( [en, int( (1-t2)*dive_st[k] + t2*dive_end[k] )], 
         [440, int(h*max(dt.iloc[dive_st[k]:dive_end[k],4].to_numpy()))], 
         linestyle=':', color='black')

plt.xlabel('   ', fontsize=20, labelpad=20)
plt.ylabel('Depth (m)', fontsize=20, labelpad=20)
title = ax.set_title('Random Forest', fontsize=20, position=(.5, 1.01),
           backgroundcolor='#D9D9D9', color='black',
           verticalalignment="bottom", horizontalalignment="center")
title._bbox_patch._mutation_aspect = 0.04
title.get_bbox_patch().set_boxstyle("square", pad=11.5)

a = plt.axes([.28, .45, .45, .4])

plt.gca().invert_yaxis()
plt.plot(np.arange( st, en )[::10], 
         dt.iloc[st:en,4].to_numpy()[::10], 
         label='_nolegend_ ', c='k', alpha=0.25, linewidth=4)
for j in buzz_dive_0[1:] :
  if (st <= s0[j]) and (e0[j] <= en) :
    plt.plot( np.arange(s0[j], e0[j])[::10],
             dt.iloc[ s0[j]:e0[j], 4 ].to_numpy()[::10], c='#D55E00', linewidth=4)
for j in buzz_dive_1[1:] :
  if (st <= s1[j]) and (e1[j] <= en) :
    plt.plot( np.arange(s1[j], e1[j])[::10],
             dt.iloc[ s1[j]:e1[j], 4 ].to_numpy()[::10], c='#0072B2', alpha=0.5, linewidth=4)


plt.plot( np.arange( st, en )[::10],
          np.repeat( max(dt.iloc[ st:en,4 ].to_numpy()[::10]) + 3, (en-st)//10) , 
          c='k', alpha=0.25, linewidth=6)
for j in buzz_dive_0[1:] :
  if (st <= s0[j]) and (e0[j] <= en) :
    plt.plot( np.arange(s0[j], e0[j])[::10],
              np.repeat( max(dt.iloc[ st:en,4 ].to_numpy()[::10]) + 3, (e0[j]-s0[j])//10) ,  c='#D55E00', linewidth=6)


plt.plot( np.arange( st, en )[::10],
          np.repeat( max(dt.iloc[ st:en,4 ].to_numpy()[::10]) + 6, (en-st)//10) , 
          c='k', alpha=0.25, linewidth=6 )
for j in buzz_dive_1 :
  if (st <= s1[j]) and (e1[j] <= en) :
    plt.plot( np.arange(s1[j], e1[j])[::10],
              np.repeat( max(dt.iloc[ st:en,4 ].to_numpy()[::10]) + 6, (e1[j]-s1[j])//10) ,  c='#0072B2', linewidth=6)

plt.xticks([])
plt.yticks([])
plt.ylim(max(dt.iloc[st:en,4].to_numpy()[::10]) + 8, 
         min(dt.iloc[st:en,4].to_numpy()[::10]) - 5)
plt.savefig('plot/pred_RF.svg', bbox_inches='tight')
plt.show()
  


#%% We do the same procedure of fitting and plotting with logistic regression (LR)

from sklearn.linear_model import LogisticRegression

# Note that 'sag' and 'saga' fast convergence is only guaranteed on features with approximately the same scale. 
# You can preprocess the data with a scaler from sklearn.preprocessing.

# LR = LogisticRegression(solver='saga', class_weight='balanced', max_iter=500, n_jobs=16)
# The problem is not specifically the rarity of events, but rather the possibility of a small number of cases on the rarer of the two outcomes.  
# If you have a sample size of 1000 but only 20 events, you have a problem. 
# If you have a sample size of 10,000 with 200 events, you may be OK. 
# If your sample has 100,000 cases with 2000 events, you’re golden.
#  [https://stats.stackexchange.com/a/227082]
LR = LogisticRegression(solver='saga', max_iter=500)

# https://datascience.stackexchange.com/a/27616
# https://datascience.stackexchange.com/a/12346
from sklearn.preprocessing import StandardScaler
scaler = StandardScaler()
x_tr = scaler.fit_transform( x_tr )
for w in range(len(whale)) :
    X_test[w] = scaler.transform( X_test[w] )
LR.fit(x_tr, y_tr)


#%% Save/Load model [https://stackoverflow.com/a/33500427]

import joblib

with open('data/all_OneHot_no-correct.LR' , 'wb') as f:
    joblib.dump(LR, f, 5)
with open('data/all_OneHot_no-correct.LR', 'rb') as f:
    LR = joblib.load(f)

with open('data/all_OneHot.LR' , 'wb') as f:
    joblib.dump(LR, f, 5)

with open('data/all_OneHot.LR', 'rb') as f:
    LR = joblib.load(f)



#%% Confusion matrix of LR

from sklearn.metrics import confusion_matrix

for i in range(len(whale)) :
  y_pred = LR.predict(X_test[i])
  print ( i, ':', confusion_matrix(y_test[i], y_pred).ravel() )

# (tn, fp, fn, tp): no rare event correction
# 0 : [22874    14    69     0]
# 1 : [23211    31   186    12]
# 2 : [27711    11    82    32]
# 3 : [6800    1   58    0]
# 4 : [6511    0   76    1]

  
# vs. rare event correction
# (tn, fp, fn, tp):
# 0 : [36044  9626    14   231]
# 1 : [35937  9432   377  1134]
# 2 : [50119  4766    62   726]
# 3 : [12302  1352    21    43]
# 4 : [10621  2431    16   109]


#%% Count buzz per dive and plot

df = [None]*len(whale)
for i in range(len(whale)) :
  y_pred = LR.predict(X_test[i])
  C = consecutive(y_pred)
  C_t = consecutive(y_test[i])
  
  if len(C_t[0]) :
    s0 = np.ravel([ x[0] for x in C_t ])
    e0 = np.ravel([ x[-1] for x in C_t ]) + 1
    s0 *= rate; e0 *= rate
  if len(C[0]) :
    s1 = np.ravel([ x[0] for x in C ])
    e1 = np.ravel([ x[-1] for x in C ]) + 1
    s1 *= rate; e1 *= rate
  
  dt = Acc[i].iloc[round( test_ratio*Acc[i].shape[0] ):,:]
  
  dif_dive = np.diff(dt.iloc[:,-3].to_numpy())
  dive_st = np.where(dif_dive > 0)[0] + 1
  dive_end = np.where(dif_dive < 0)[0] + 1
  if len(dive_st) < len(dive_end) :
    dive_st = np.r_[0, dive_st]
  if len(dive_st) > len(dive_end) :
    dive_end = np.r_[dive_end, dt.shape[0]]
  
  
  buzz_dive = np.zeros( (len(dive_st), 8) )
  for k in range(len(dive_st)) :
    depth = dt.iloc[dive_st[k]:dive_end[k],4].to_numpy()
    buzz_dive[k, 0] = i
    buzz_dive[k, 1] = k
    buzz_dive[k, 2] = depth.max()
    buzz_dive[k, 3] = depth.shape[0]/100
    
    if len(C_t[0]) :
      buzz_dive_0 = np.where( (s0 >= dive_st[k]) & (s0 < dive_end[k]) )[0]
      buzz_dive[k, 4] = len(buzz_dive_0)
      buzz_dive[k, 6] = (e0[buzz_dive_0] - s0[buzz_dive_0]).sum()
    if len(C[0]) :
      buzz_dive_1 = np.where( (s1 >= dive_st[k]) & (s1 < dive_end[k]) )[0]
      buzz_dive[k, 5] = len(buzz_dive_1)
      buzz_dive[k, 7] = (e1[buzz_dive_1] - s1[buzz_dive_1]).sum()
  
  df[i] = pd.DataFrame(data = buzz_dive , 
                       columns = ['ID','Dive_no','MaxDepth','Duration',
                                  'True','Prediction','Len_true','Len_pred'] )
  df[i]['Name'] = pd.Series(np.repeat(whale[i], len(dive_st)))


Dive_res = pd.concat(df)

from collections import Counter

# https://stackoverflow.com/a/57739795
Dive_res['hash'] = pd.util.hash_pandas_object(Dive_res.loc[:,['True','Prediction']], 
                                              index=False)
# https://stackoverflow.com/a/41678874
Dive_res['Num_point'] = Dive_res['hash'].map(dict(Counter(Dive_res['hash'].to_numpy())))
Dive_res['Diff'] = Dive_res['Prediction'].to_numpy() - Dive_res['True'].to_numpy()
Dive_res['Diff_len'] = Dive_res['Len_pred'].to_numpy() - Dive_res['Len_true'].to_numpy()


#%% For each whale

Res = np.zeros( (6,9))
for i in range(len(whale)) :
  y_pred = LR.predict(X_test[i])
  L0 = []
  for p in [99, 50, 25, 0] :
    _, ptg, _, _ = hit_rate(y_pred, y_test[i], ovl_ptg=p)
    if ptg is not None :
      L0.append( round(ptg, 2) )
    else :
      L0.append( None )
  for t in [1, 2, 3, 4, 5] :
    R, _, hit, missed = hit_rate(y_pred, y_test[i], t)
    if R is not None :
      L0.append( round(R, 2) )
    else :
      L0.append( None )
  Res[i,:] = np.ravel(L0)


# https://codebeautify.org/htmlviewer/
pd.DataFrame(data=Res.T, columns = whale,
             index=np.ravel(['> 99%', '> 50%', '> 25%', '> 0%', 
                             '1 s', '2 s', '3 s', '4 s', '5 s']) ).to_html()

#%% For all whales

y_pred_all = LR.predict(np.vstack(X_test))
y_test_all = np.hstack(y_test)

L0 = []
for p in [99, 50, 25, 0] :
  _, ptg, _, _ = hit_rate(y_pred_all, y_test_all, ovl_ptg=p)
  if ptg is not None :
    L0.append( round(ptg, 2) )
  else :
    L0.append( None )
for t in [1, 2, 3, 4, 5] :
  R, _, hit, missed = hit_rate(y_pred_all, y_test_all, t)
  if R is not None :
    L0.append( round(R, 2) )
  else :
    L0.append( None )
Res[-1,:] = np.ravel(L0)

res = pd.DataFrame(data= {'index': np.tile(np.ravel(['> 99%', '> 50%', '> 25%', '> 0%', 
                                           '1 s', '2 s', '3 s', '4 s', '5 s']), Res.shape[0] ) ,
                        'val': Res.ravel() , 
                        'whale': np.repeat(whale + ['All'],Res.shape[1])  })
res['method'] = 'Logistic Regression'

# https://matplotlib.org/3.1.0/api/_as_gen/matplotlib.pyplot.figure.html
import matplotlib
import matplotlib.pyplot as plt
matplotlib.rcParams["figure.dpi"] = 300
matplotlib.rcParams["figure.figsize"] = [11, 5]
matplotlib.rcParams['svg.fonttype'] = 'none'
plt.style.use('ggplot')


i = 2; k = 6

y_pred = LR.predict(X_test[i])
C = consecutive(y_pred)
C_t = consecutive(y_test[i])
  
if len(C[0]) :
  s1 = np.ravel([ x[0] for x in C ])
  e1 = np.ravel([ x[-1] for x in C ]) + 1
  s1 *= rate; e1 *= rate
if len(C_t[0]) :
  s0 = np.ravel([ x[0] for x in C_t ])
  e0 = np.ravel([ x[-1] for x in C_t ]) + 1
  s0 *= rate; e0 *= rate

dt = Acc[i].iloc[round( test_ratio*Acc[i].shape[0] ):,:].copy() # [1371800 rows x 13 columns]

from collections import Counter as cnt
dif_dive = np.diff(dt.iloc[:,-3].to_numpy())

pd.DataFrame.from_dict(cnt(dif_dive), orient='index').reset_index()

dive_st = np.where(dif_dive > 0)[0] + 1
dive_end = np.where(dif_dive < 0)[0] + 1
if len(dive_st) < len(dive_end) :
  dive_st = np.r_[0, dive_st]
if len(dive_st) > len(dive_end) :
  dive_end = np.r_[dive_end, dt.shape[0]]

if len(C_t[0]) :
  buzz_dive_0 = np.where( (s0 >= dive_st[k]) & (s0 < dive_end[k]) )[0]
if len(C[0]) :
  buzz_dive_1 = np.where( (s1 >= dive_st[k]) & (s1 < dive_end[k]) )[0]
  

fig, ax = plt.subplots()
ax.tick_params(axis='both', which='major', labelsize=15)

plt.gca().invert_yaxis()
# https://stackoverflow.com/a/1144137
X_axis = np.arange(dive_st[k], dive_end[k])[::10]
tick_locs = np.arange(X_axis[0], X_axis[-1], 2*60*rate)
tick_lbls = ((tick_locs-X_axis[0])/(60*rate)).round().astype(int)
plt.xticks(tick_locs.tolist(), tick_lbls.tolist() )


plt.plot(np.arange(dive_st[k], dive_end[k])[::10], 
         dt.iloc[dive_st[k]:dive_end[k],4].to_numpy()[::10], 
        label='_nolegend_', c='k', alpha=0.25, linewidth=6)
for j in buzz_dive_0 :
  plt.plot( np.arange(s0[j], e0[j])[::10], 
           dt.iloc[ s0[j]:e0[j], 4 ].to_numpy()[::10], c='#D55E00', linewidth=6)
for j in buzz_dive_1 :
  plt.plot( np.arange(s1[j], e1[j])[::10], 
           dt.iloc[ s1[j]:e1[j], 4 ].to_numpy()[::10], c='#0072B2', alpha=0.5, linewidth=6)

plt.xlabel('Time (minutes)', fontsize=20, labelpad=20)
plt.ylabel('   ', fontsize=20, labelpad=20)
# plt.title('Logistic Regression', fontsize=20 )
title = ax.set_title('Logistic Regression', fontsize=20, position=(.5, 1.01),
           backgroundcolor='#D9D9D9', color='black',
           verticalalignment="bottom", horizontalalignment="center")
title._bbox_patch._mutation_aspect = 0.04
title.get_bbox_patch().set_boxstyle("square", pad=10.46)


offset1 = 6000
offset2 = 24000
st = min(s0[min(buzz_dive_0)],s1[min(buzz_dive_1)])
en = max(e0[max(buzz_dive_0)],e1[max(buzz_dive_1)])
st += offset1
en -= offset2

# plt.axline((en-10, 520), (en, 520))
plt.plot( [st, en], [440, 440], linestyle=':', color='black')
plt.plot( [st, en], [530, 530], linestyle=':', color='black')
plt.plot( [st, st], [440, 530], linestyle=':', color='black')
plt.plot( [en, en], [440, 530], linestyle=':', color='black')


h = 0.62
t = 0.17
plt.plot( [st, int( (1-t)*dive_st[k] + t*dive_end[k] )], 
         [440, int(h*max(dt.iloc[dive_st[k]:dive_end[k],4].to_numpy()))], 
         linestyle=':', color='black')
t2 = t + 0.64
plt.plot( [en, int( (1-t2)*dive_st[k] + t2*dive_end[k] )], 
         [440, int(h*max(dt.iloc[dive_st[k]:dive_end[k],4].to_numpy()))], 
         linestyle=':', color='black')

a = plt.axes([.28, .45, .45, .4])
plt.gca().invert_yaxis()
plt.plot(np.arange( st, en )[::10], 
         dt.iloc[st:en,4].to_numpy()[::10], 
         label='_nolegend_ ', c='k', alpha=0.25, linewidth=4)
for j in buzz_dive_0[1:] :
  if (st <= s0[j]) and (e0[j] <= en) :
    plt.plot( np.arange(s0[j], e0[j])[::10],
             dt.iloc[ s0[j]:e0[j], 4 ].to_numpy()[::10], c='#D55E00', linewidth=4)
for j in buzz_dive_1[1:] :
  if (st <= s1[j]) and (e1[j] <= en) :
    plt.plot( np.arange(s1[j], e1[j])[::10],
             dt.iloc[ s1[j]:e1[j], 4 ].to_numpy()[::10], c='#0072B2', alpha=0.5, linewidth=4)


plt.plot( np.arange( st, en )[::10],
          np.repeat( max(dt.iloc[ st:en,4 ].to_numpy()[::10]) + 3, (en-st)//10) , 
          c='k', alpha=0.25, linewidth=6)
for j in buzz_dive_0[1:] :
  if (st <= s0[j]) and (e0[j] <= en) :
    plt.plot( np.arange(s0[j], e0[j])[::10],
              np.repeat( max(dt.iloc[ st:en,4 ].to_numpy()[::10]) + 3, (e0[j]-s0[j])//10) ,  c='#D55E00', linewidth=6)


plt.plot( np.arange( st, en )[::10],
          np.repeat( max(dt.iloc[ st:en,4 ].to_numpy()[::10]) + 6, (en-st)//10) , 
          c='k', alpha=0.25, linewidth=6 )
for j in buzz_dive_1 :
  if (st <= s1[j]) and (e1[j] <= en) :
    plt.plot( np.arange(s1[j], e1[j])[::10],
              np.repeat( max(dt.iloc[ st:en,4 ].to_numpy()[::10]) + 6, (e1[j]-s1[j])//10) ,  c='#0072B2', linewidth=6)

plt.xticks([])
plt.yticks([])
plt.ylim(max(dt.iloc[st:en,4].to_numpy()[::10]) + 8, 
         min(dt.iloc[st:en,4].to_numpy()[::10]) - 5)
plt.savefig('plot/pred_LR.svg', bbox_inches='tight')
plt.show()

#%% Finally, we have 3 figures from 3 models U-Net, Random Forest and Logistic Regression
# We can concatenate them together to have Figure 11 in the article

import svgutils.transform as sg
import sys 

plots = [sg.fromfile(x).getroot() for x in ['plot/pred_UNet.svg','plot/pred_RF.svg', 
                                            'plot/pred_LR.svg']]
plots[1].moveto(0, 340, scale=1)
plots[2].moveto(0, 340*2, scale=1)

#create new SVG figure
fig = sg.SVGFigure("19cm", "28cm")
fig.append(plots)
fig.save("plot/fig_final.svg")
