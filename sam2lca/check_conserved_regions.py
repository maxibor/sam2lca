import numpy as np
from scipy.stats import binned_statistic


def compute_coverage(pysam_cov):
    """Computes average coverage of a reference
    Args:
        pysam_cov (np.array): Four dimensional array of coverage for each base
    Returns:
        np.array: 1D int array of coverage at each base
    """
    return np.sum(pysam_cov, axis=0)


def zscore(cov):
    """Compute zscore

    Args:
        cov (np.array): 1D numpy array of coverage (float)
    Returns:
        np.array: zscore(float) of coverage
    """
    mean = np.mean(cov)
    stdev = np.std(cov)
    if mean == stdev == 0:
        return cov
    else:
        return (cov - mean) / stdev


def flag_conserved_regions(cov_array, window_size=500, zscore_thresh=1.65):
    """Flag ultra-conserved regions by checking coverage with zscore

    Args:
        cov_array (np.array): 1D int array of coverage at each base
        window_size(int): size of sliding window
        zscore_thresh(float): zscore threshold
    Returns:
        list: list of start and end positions of windows flagged as
            conserved [[start, end],[start,end]]

    """
    nb_windows = int(cov_array.size / window_size)
    cov_bin_median = binned_statistic(
        np.arange(cov_array.size), cov_array, statistic="median", bins=nb_windows
    )
    cov_bin_zscore = zscore(cov_bin_median[0])
    is_conserved = cov_bin_zscore > zscore_thresh
    conserved_regions = cov_bin_median[1][:-1].astype(int)[is_conserved]
    cons_range = []
    for i in conserved_regions:
        cons_range.append((i, min(i + window_size, cov_array.size - 1)))

    return cons_range


def is_in_conserved(read, cons_ranges):
    """Check if read is in a conserved region

    Args:
        read (pysam read): Read class from PySam
        cons_ranges (list): list of start and end positions of windows flagged as
            conserved [[start, end],[start,end]]
    """
    for r in cons_ranges:
        if read.reference_start > r[1]:
            continue
        if read.reference_start > r[0] and read.reference_end < r[1]:
            return True
    return False
