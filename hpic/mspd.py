"""
Multiscale peak detection method based on continuous wavelet transform 
"""
from __future__ import division
import numpy as np
from scipy.signal import fftconvolve
from scipy.stats import scoreatpercentile, mode
from collections import deque
from math import ceil


def mexican_hat(points, a):
    A = 2 / (np.sqrt(3 * a) * (np.pi ** 0.25))
    wsq = a ** 2
    vec = np.arange(0, points) - (points - 1.0) / 2
    tsq = vec ** 2
    mod = (1 - tsq / wsq)
    gauss = np.exp(-tsq / (2 * wsq))
    total = A * mod * gauss
    return total

def cwt(data, wavelet, widths):
    output = np.zeros([len(widths), len(data)])
    for ind, width in enumerate(widths):
        wavelet_data = wavelet(min(10 * width, len(data)), width)
        output[ind, :] = fftconvolve(data, wavelet_data,
                                     mode='same')
    return output

def local_extreme(data, comparator,
                  axis=0, order=1, mode='clip'):
    if (int(order) != order) or (order < 1):
        raise ValueError('Order must be an int >= 1')
    datalen = data.shape[axis]
    locs = np.arange(0, datalen)
    results = np.ones(data.shape, dtype=bool)
    main = data.take(locs, axis=axis, mode=mode)
    for shift in range(1, order + 1):
        plus = data.take(locs + shift, axis=axis, mode=mode)
        minus = data.take(locs - shift, axis=axis, mode=mode)
        results &= comparator(main, plus)
        results &= comparator(main, minus)
    return results

def ridge_detection(local_max, row_best, col, n_rows, n_cols, minus=True, plus=True):
    cols = deque()
    rows = deque()
    cols.append(col)
    rows.append(row_best)
    col_plus = col
    col_minus = col
    for i in range(1, n_rows):
        row_plus = row_best + i
        row_minus = row_best - i
        segment_plus = 1
        segment_minus = 1
        if minus and row_minus > 0 and segment_minus < col_minus < n_cols - segment_minus - 1:
            if local_max[row_minus, col_minus + 1]:
                col_minus += 1
            elif local_max[row_minus, col_minus - 1]:
                col_minus -= 1
            elif local_max[row_minus, col_minus]:
                col_minus = col_minus
            else:
                col_minus = -1
            if col_minus != -1:
                rows.appendleft(row_minus)
                cols.appendleft(col_minus)
        if plus and row_plus < n_rows and segment_plus < col_plus < n_cols - segment_plus - 1:
            if local_max[row_plus, col_plus + 1]:
                col_plus += 1
            elif local_max[row_plus, col_plus - 1]:
                col_plus -= 1
            elif local_max[row_plus, col_plus]:
                col_plus = col_plus
            else:
                col_plus = -1
            if col_plus != -1:
                rows.append(row_plus)
                cols.append(col_plus)
        if (minus and False == plus and col_minus == -1) or \
                (False == minus and True == plus and col_plus == -1) or \
                (True == minus and True == plus and col_plus == -1 and col_minus == -1):
            break
    return rows, cols

def peaks_position(vec, ridges, cwt2d, wnd=2):
    n_cols = cwt2d.shape[1]
    negs = cwt2d < 0
    local_minus = local_extreme(cwt2d, np.less, axis=1, order=1)
    #zero_crossing = np.abs(np.diff(np.sign(cwt2d))) / 2
    negs |= local_minus
    negs[:, [0, n_cols - 1]] = True
    ridges_select = []
    peaks = []
 
    for ridge in ridges:
        inds = np.where(cwt2d[ridge[0, :], ridge[1, :]] > 0)[0]
        if len(inds) > 0:
            col = int(mode(ridge[1, inds])[0][0])
            rows = ridge[0, :][(ridge[1, :] == col)]
            row = rows[0]
            cols_start = max(col - np.where(negs[row, 0:col][::-1])[0][0], 0)
            cols_end = min(col + np.where(negs[row, col:n_cols])[0][0], n_cols)
            if cols_end > cols_start:
                inds = range(cols_start, cols_end)
                peaks.append(inds[np.argmax(vec[inds])])
                ridges_select.append(ridge)
        elif ridge.shape[1] > 2: # local wavelet coefficients < 0
            cols_accurate = ridge[1, 0:int(ceil(ridge.shape[1] / 2))]
            cols_start = max(np.min(cols_accurate) - 3, 0)
            cols_end = min(np.max(cols_accurate) + 4, n_cols - 1)
            inds = range(cols_start, cols_end)
            if len(inds) > 0:
                peaks.append(inds[np.argmax(vec[inds])])
                ridges_select.append(ridge)
    peaks_1 = np.unique(peaks)
    if len(peaks)> peaks_1.shape[0]:
        ridges_len = np.array([ridge.shape[1] for ridge in ridges_select])
        ridges_refine = []
        for peak in peaks_1:
            inds = np.where(peaks == peak)[0]
            ridge = ridges_select[inds[np.argmax(ridges_len[inds])]]
            ridges_refine.append(ridge)
        return list(peaks_1),ridges_refine
    else:
        return peaks,ridges_select

def ridges_detection(cwt2d, vec):
    n_rows = cwt2d.shape[0]
    n_cols = cwt2d.shape[1]
    local_max = local_extreme(cwt2d, np.greater, axis=1, order=1)
    ridges = []
    rows_init = np.array(range(1, 6))
    cols_small_peaks = np.where(np.sum(local_max[rows_init, :], axis=0) > 0)[0]
    for col in cols_small_peaks:
        best_rows = rows_init[np.where(local_max[rows_init, col])[0]]
        rows, cols = ridge_detection(local_max, best_rows[0], col, n_rows, n_cols, True, True)
        staightness = 1 - float(sum(abs(np.diff(cols)))) / float(len(cols))
        if len(rows) >=4 and staightness >=0 and \
            not(
            len(ridges) > 0 and
            np.array_equal(np.array(cols),ridges[-1][1])
        ):
            ridges.append(np.array([rows, cols], dtype=np.int32))
    return ridges

def signal_noise_ratio(cwt2d, ridges, peaks):
    n_cols = cwt2d.shape[1]
    row_one = cwt2d[0, :]
    row_one_del = np.delete(row_one, np.where(abs(row_one) < 10e-5))
    t = 3 * np.median(np.abs(row_one_del - np.median(row_one_del))) / 0.67
    row_one[row_one > t] = t
    row_one[row_one < -t] = -t
    noises = np.zeros(len(peaks))
    signals = np.zeros(len(peaks))
    for ind, val in enumerate(peaks):
        hf_window = ridges[ind].shape[1]
        window = range(int(max([val - hf_window, 0])), int(min([val + hf_window, n_cols])))
        noises[ind] = scoreatpercentile(np.abs(row_one[window]), per=90)
        signals[ind] = np.max(cwt2d[ridges[ind][0, :], ridges[ind][1, :]])
    return np.abs((signals + np.finfo(float).eps) / (noises + np.finfo(float).eps)), signals

def peaks_detection(data, scales, min_snr=3, intensity=200):
    """
    Peak detection for pure ion chromatogram of LC-MS

    Arguments:
        data: string
            path to the dataset locally
        scales: vector
            1-D array of widths to use for calculating the CWT matrix. 
        min_snr: float
            minimum signal noise ratio.
        intensity:
            minimum intensity of the peak

    Returns:
        Numpy array: shape = (n, 9).
            n is the number of the detected peaks
            The meaning of each column is as following
                # 0 -> m/z value 
                # 1 -> retention time at apex
                # 2 -> intensity at apex
                # 3 -> start of retention time 
                # 4 -> end of retention time  
                # 5 -> intensity in CWT space
                # 6 -> SNR in CWT space
                # 7 -> index of the peak
                # 8 -> length of the peak
    """
    vec = data[:,2]
    cwt2d = cwt(vec, mexican_hat, scales)
    ridges = ridges_detection(cwt2d, vec)
    peaks, ridges = peaks_position(vec, ridges, cwt2d)
    peak_list = np.zeros((len(peaks),9))
    peak_list[:,:3],peak_list[:,3],peak_list[:,4],peak_list[:,7],peak_list[:,8]= data[peaks][:,[1,0,2]],data[0,0],data[-1,0],peaks,data.shape[0]
    peak_list[:,6], peak_list[:,5] = signal_noise_ratio(cwt2d, ridges, peaks)
    peak_list = peak_list[(peak_list[:,6]>=min_snr)&(peak_list[:,2]>=intensity)]
    if peak_list.shape[0]>1:
        peak_list=peak_list[np.argsort(peak_list[:,7])]
        for i in range(peak_list.shape[0]-1):
            rt1_index = int(np.argmin(vec[int(peak_list[i,7]):int(peak_list[i+1,7])])+peak_list[i,7])
            peak_list[i,4],peak_list[i+1,3] = data[rt1_index,0],data[rt1_index,0] 
    return peak_list
