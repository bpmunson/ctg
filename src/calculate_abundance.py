import numpy  as np
import pandas as pd
from scipy.interpolate import UnivariateSpline
from scipy.stats import gaussian_kde
from scipy.interpolate import CubicSpline

def calculate_abundance(counts_file, min_counts_threshold=10, **kwargs):
    if isinstance(counts_file, str):
        counts_df = read_counts_file(counts_file, **kwargs)

    else:
        #TODO: Validation code needed
        counts_df = counts_file

    counts_df[counts_df == 0] = 1

    thresholds = []
    for col in counts_df.columns:
        thresholds.append(analyzeCountsDist(counts_df[col], min_counts_threshold=min_counts_threshold))

    # return pd.DataFrame(np.vstack([counts_df.columns, thresholds]).T,
    #                    columns=['sampleName', 'log2CountsThresh'],
    #                    index=range(1,len(counts_df.columns) + 1))

    thresholds_df = pd.DataFrame(thresholds,
                       columns=['log2CountsThresh'],
                       index=counts_df.columns)

    thresholds_df.index.name = 'sampleName'

    return thresholds_df

def read_counts_file(fileName, index_col=0, data_col_start=4, sep='\t'):
    counts_df = pd.read_csv(fileName, index_col=index_col, sep=sep)

    return counts_df.iloc[:,data_col_start:]

def analyzeCountsDist(counts_series, binwidth=0.05, min_counts_threshold=10, plot=False):
    log_counts = np.log2(counts_series)


    x = np.arange(-binwidth, log_counts.values.max(), binwidth)
    counts, bin_edges = np.histogram(np.asarray(log_counts), x)

    bin_mids = _averages(bin_edges)

    scaleFactor = np.sum(counts)*binwidth
    density = gaussian_kde(np.asarray(log_counts), bw_method='silverman').evaluate(x)*scaleFactor

    hist = np.vstack([bin_mids, counts])
    nonZeroLog2Counts = hist[:, hist[1,:] > 0]
    f = CubicSpline(nonZeroLog2Counts[0,:], nonZeroLog2Counts[1,:])

    xs = np.arange(np.min(nonZeroLog2Counts[0,:]), np.max(nonZeroLog2Counts[0,:]), 0.5)
    log2CountsSpline = f(xs)

    min_counts = smallestMin(xs, log2CountsSpline, min_counts_threshold)

    if plot:
        plt.bar(bin_mids, counts)
        plt.plot(xs, log2CountsSpline, color='red')

    if check_threshold(min_counts, f, (xs[0], xs[-1])):
        return min_counts

    else:
        return min_counts_threshold


def smallestMin(counts, density, min_counts_threshold):
    '''Smallest-local-minimum-in-valley method'''
    log_min_counts = np.log2(min_counts_threshold)

    min_indices = findMinIndices(density)
    max_indices = findMaxIndices(density)

    #print(counts[max_indices])
    #print(counts[min_indices])

    max_counts = counts[max_indices]
    max_counts_masked = np.ma.array(max_counts, mask = ~(max_counts > log_min_counts))
    max_counts = np.max(max_counts_masked)

    #print(max_counts)

    if max_counts is np.ma.core.masked:
        return None

    min_counts = counts[min_indices]
    min_counts_masked = np.ma.array(min_counts,
                                    mask = ~((min_counts >= log_min_counts) \
                                        & (min_counts < max_counts)))

    #print(min_counts, log_min_counts)

    min_counts = np.min(min_counts_masked)

    if min_counts is np.ma.core.masked:
        return None

    return min_counts

def check_threshold(value, func, bounds, threshold = 0.95):
    auc = func.integrate(bounds[0], bounds[1])

    if func.integrate(value, bounds[1])/auc > 1-threshold:
        return True
    else:
        return False

# def nearest_point(counts, density, min_counts_threshold):
#     '''Near-point-of-spline-and-density method used'''
#     log_min_counts = np.log2(min_counts_threshold)


'''Helper functions'''
def findMinIndices(xs):
    tmp = np.diff(np.sign(np.diff(np.asarray(xs))))
    return np.where(tmp == 2)[0] + 1

def findMaxIndices(xs):
    tmp = np.diff(np.sign(np.diff(np.asarray(xs))))
    return np.where(tmp == -2)[0] + 1

def _averages(series):
    series = np.asarray(series)

    return (series[1:] + series[:-1])/2

if __name__ == "__main__":
    counts_file = '/cellar/users/samsonfong/side_projects/17_12_14_skewed-GI/ctg_output/full/testing_20171214113528/testing_20171214113528_timepoint_counts.txt'

    thresholds = calculate_abundance(counts_file)
    print(thresholds)
