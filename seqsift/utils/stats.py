#! /usr/bin/env python

"""
A collection of classes and functions for stats operations.

Many of the classes and functions are stolen from `PyMsBayes` (Copywright Jamie
Oaks; licensed under the GPL; https://github.com/joaks1/PyMsBayes).
"""

import sys
import os
import math
import operator
import decimal
import fractions
try:
    from cStringIO import StringIO
except ImportError:
    from io import StringIO

from seqsift.utils import GLOBAL_RNG, iteritems, itervalues
from seqsift.utils.messaging import get_logger

_LOG = get_logger(__name__)

def get_counts(elements):
    keys = sorted(set(elements))
    counts = dict(zip(keys, [0 for i in range(len(keys))]))
    for el in elements:
        counts[el] += 1
    return counts

def get_freqs(elements):
    counts = get_counts(elements)
    total = float(sum(itervalues(counts)))
    freqs = {}
    for k, v in iteritems(counts):
        freqs[k] = v / total
    return freqs

def median(samples):
    """
    Return the median value from a list of samples.
    """
    s = sorted(samples)
    n = len(s)
    if n < 1:
        raise ValueError('empty samples')
    mdn = None
    if n % 2 == 0:
        mdn = (s[((n // 2) - 1)] + s[(n // 2)]) / 2
    else:
        mdn = s[((n - 1) // 2)]
    return mdn

def has_floats(samples):
    for s in samples:
        if isinstance(s, float) or isinstance(s, fractions.Fraction) or \
                isinstance(s, decimal.Decimal):
            return True
    return False

def get_bin_width(samples, algorithm = 'freedman-diaconis'):
    """
    Return "best" histogram bin width for samples using specified `algorithm`.
    Options for `algorithm` argument include:

    `algorithm = 'f'|'fd'|'freedman-diaconis'`
        Use the Freedman-Diaconis rule to calculate bin widths as
        2*IQR(samples) / n^(1/3), where n is the number of samples.

    `algorithm = 'c'|'custom'`
        Use custom (hacked) rescaled version of the Freedman-Diaconis rule to
        calculate. It returns the Freedman-Diaconis bin width multiplied by
        (n^(8/7) / (2 * n).
        This is very similar to the vanilla Freedman-Diaconis bin width at
        small sample sizes, but returns a more conservative (wider) bin width
        as sample sizes get large. I have found the F-D bin widths to be too
        narrow at large samples sizes (e.g., n > 10000) and this adjustment can
        allow more consistent estimation of the mode, yet be more precise than
        Sturges` and Doane's formulae.

    `algorithm = 's'|'sturges'`
        Use Sturges' formula to calculate the number of bins as k = LOG2(n) + 1
        where n is the number of samples, then return sample_range / k as the
        bin width.

    `algorithm = 'r'|'rice'`
        Use Rice Rule to estimate the number of bins as k = 2n^(1/3), where n
        is the number of samples, then return sample_range / k as the bin
        width.

    `algorithm = 'd'|'doane'`
        Use Doane's formula to estimate k:
            k = 1 + Log2(n) + Log2(1 + (|skewness|/sigma))
        where n is the number of samples, and
            sigma = sqrt((6(n-2)/((n+1)(n+3))))
        Then return sample_range / k as the bin width.
    """
    if not samples:
        return None
    if len(samples) < 2:
        return math.fabs(samples[0])
    a = algorithm.strip().lower()
    n = len(list(samples))
    if a in ['c', 'custom']:
        scaler = (n ** (float(8)/7)) / (2 * n)
        return scaler * get_bin_width(samples, 'freedman-diaconis')
    if a in ['f', 'fd', 'freedman-diaconis']:
        iqr = quantile(samples, 0.75) - quantile(samples, 0.25)
        return 2 * (float(iqr) / (n ** (float(1)/3)))
    elif a in ['s', 'sturges']:
        k = math.ceil(math.log(n, 2) + 1)
        return (max(samples) - min(samples)) / k
    elif a in ['d', 'doane']:
        if samples < 3:
            return get_bin_width(samples, 'freedman-diaconis')
        sigma = math.sqrt((6 * (n - 2)) / float((n + 1) * (n + 3)))
        ss = SampleSummarizer(samples)
        k = 1 + math.log(n, 2) + math.log((1 + (math.fabs(
                ss.skewness) / sigma)), 2)
        return (max(samples) - min(samples)) / k
    elif a in ['r', 'rice']:
        k = math.ceil(2 * (n ** (float(1)/3)))
        return (max(samples) - min(samples)) / k
    raise ValueError('unsupported `a` argument: {0!r}'.format(a))

def mode_list(samples, bin_width = 'auto', zero_value = 'boundary'):
    """
    Return a list of modes, or mode bins, from a list of values.

    Arguments include:

        `samples` is an iterable set of values, which can be integers, strings
        or floats.

        `bin_width` controls the behavior of the mode estimation, with the
        following options:
        
            `bin_width = 'a'|'auto'` - The default. The function automatically
            determines whether to treat the values as discrete or continuous by
            checking for floating point numbers in the sample. If there are no
            floats the samples are treated as discrete and a list of the most
            common values is returned. If there are floats, a bin width is
            determined by calling `get_bin_width(samples, algorithm='custom')`.
            The values are then binned using this bin width and the
            `zero_value` argument, and a list of tuples is returned, each tuple
            represents the lower and upper bounds of the most common bins.

            `bin_width = None|0` - The samples are treated as
            discrete and a list of the most common values is returned.

            `bin_width = <NUMBER OTHER THAN ZERO>` - The samples are treated as
            floats and are binned into categories of width `bin_width` to
            determine the mode.

            `bin_width =
                'c'|'custom'
                'f'|'fd'|'freedman-diaconis'|
                's'|'sturges'|
                'r'|'rice'|
                'd'|'doane'`
            The 'best' bin width is determined using the specified algorithm
            (see `get_bin_width` function for details regarding the algorithm
            options).
 
        `zero_value` zero always corresponds to a bin, and this option controls
        whether zero is at the center of a bin or at the edge. Options include:

            `zero_value = 'b'|'boundary'` - zero is a boundary between bins.
            In most cases choosing between 'b' and 'c' will be arbitrary, but
            if the samples are bounded by zero (i.e., the parameter is either
            strictly positive or negative, zero should be set as a boundary.

            `zero_value = 'c'|'center'` - zero is at the center of a bin.  If
            the samples are bounded by zero, use 'boundary' instead.  However,
            this option can be useful if the samples span zero and are
            suspected to be centered at zero.

    The function returns:

        If values are treated as discrete (i.e., `bin_width = None`), a list
        containing the most common values is returned. The list will contain
        multiple values if there is a tie.

        If values are treated as floats (i.e. `bin_width != None`), a list of
        tuples containing the lower and upper boundaries of the most common
        bins is returned. The list will contain multiple tuples each
        representing a most common bin, if there is a tie.

    Some examples:
        >>> from pymsbayes.utils.stats import mode_list
        >>> x = range(10) + [2]
        >>> mode_list(x)  # treat values as discrete by default
        [2]
        >>> x += [6]
        >>> mode_list(x)  # a tie!
        [2, 6]
        >>> x = ['a', 'b', 'b', 'c', 'c', 'b']
        >>> # strings work too when treated as discrete
        >>> mode_list(x)
        ['b']
        >>> import random
        >>> x = [random.Random().expovariate(1) for i in range(10000)]
        >>> # specify bin width for continuous values
        >>> mode_list(x, bin_width='auto')
        [(0.0, 0.10405355148832289)]
        >>> x = [random.Random().normalvariate(1, 1) for i in range(10000)]
        >>> mode_list(x, bin_width='auto')
        [(0.8910191831744725, 1.0183076379136828)]
        >>> x = [random.Random().normalvariate(0, 1) for i in range(10000)]
        >>> # zero is a bin boundary by default
        >>> mode_list(x, bin_width='auto') 
        [(-0.1263661814981197, 0.0)]
        >>> # specify zero_value as bin center to get mode that spans zero
        >>> mode_list(x, bin_width='auto', zero_value='center')
        [(-0.06318309074905985, 0.06318309074905985)]

    The beginnings of this function were based on the mode function in DendroPy
    (Copyright Jeet Sukumaran and Mark T. Holder; licensed under BSD License;
    http://pythonhosted.org/DendroPy/):

    Sukumaran, J. and M. T. Holder. 2010. DendroPy: a Python library
    for phylogenetic computing. Bioinformatics 26: 1569-1571.
    """
    if not samples:
        raise ValueError('empty samples')
    if len(list(samples)) == 1:
        return list(samples)
    zero_value = zero_value.strip().lower()
    discrete = False
    if not bin_width:
        discrete = True
    elif hasattr(bin_width, 'lower'):
        bin_width = bin_width.strip().lower()
        if bin_width in ['a', 'auto']:
            discrete = not has_floats(samples)
            if discrete:
                bin_width = None
            else:
                bin_width = 'c'
        else:
            discrete = False
        if not discrete:
            bin_width = get_bin_width(samples, bin_width)
            if bin_width == 0.0:
                bin_width = (max(samples) - min(samples)) / float(10)
                if bin_width == 0.0:
                    bin_width = 0.001
    if not discrete:
        bw = float(bin_width)
    counts = {}
    for s in samples:
        if discrete:
            index = s
        else:
            if zero_value in ['b', 'boundary']:
                index = int(math.floor(s / bw))
                bounds = (0.0, bw)
            elif zero_value in ['c', 'center']:
                index = int(math.floor((s / bw) + 0.5))
                bounds = ((bw / 2), (bw / 2))
            else:
                raise ValueError('unsupported `zero_value` argument: '
                        '{0!r}'.format(zero_value))
        counts[index] = counts.get(index, 0) + 1
    count_tups = sorted(iteritems(counts), key = operator.itemgetter(1),
            reverse = True)
    max_count = count_tups[0][1]
    if discrete:
        return [val for val, cnt in count_tups if cnt >= max_count]
    return [((val * bin_width) - bounds[0], (val * bin_width) + bounds[1]) \
            for val, cnt in count_tups if cnt >= max_count]

def get_hpd_interval(samples, interval_prob = 0.95):
    """
    Return tuple of highest posterior density (HPD).

    The interval is estimated via a sliding window to find the narrowest
    interval that contains the specified proportion of samples.
    """
    if not samples:
        raise ValueError('empty samples')
    if interval_prob <= 0.0:
        raise ValueError('hpd probability interval is zero')
    samples = sorted(samples)
    tail_prob = 1 - interval_prob
    n = len(samples)
    num_exclude = int(round((n * tail_prob) - 0.5))
    if num_exclude < 1:
        num_exclude = 0
    widths = []
    # sliding window to find possible interval widths
    for i in range(num_exclude+1):
        lower = samples[i]
        upper = samples[(n - 1) - num_exclude + i]
        widths.append(upper - lower)
    min_width = min(widths)
    min_index = widths.index(min_width)
    return(samples[min_index], samples[(n - 1) - num_exclude + min_index])

def quantile(samples, p): 
    """
    Return quantile associated with probability `p`.

    Modified from code by Wai Yip Tung, licensed under PSF license and available
    here:
    http://code.activestate.com/recipes/511478-finding-the-percentile-of-the-values/
    """
    if not samples:
        raise ValueError('empty samples')
    s = sorted(samples)
    k = (len(s) - 1) * p
    f = math.floor(k)
    c = math.ceil(k)
    if f == c:
        return s[int(k)]
    d0 = s[int(f)] * (c - k)
    d1 = s[int(c)] * (k - f)
    return d0 + d1

def quantile_95(samples):
    """
    Return tuple of interval of 2.5% and 97.5% quantiles.
    """
    return (quantile(samples, 0.025), quantile(samples, 0.975))

def get_summary(samples, bin_width = 'auto'):
    """
    Return a dict of summaries calculated from the samples.

    The dict has the following items:
        'n': sample_size
        'mean': mean
        'median': median
        'modes': mode (tuple if binning)
        'variance': variance
        'range': range
        'hpdi_95': tuple of 95% highest posterior density interval
        'qi_95': tuple of 2.5% to 97.5% quantile interval
    """
    ss = SampleSummarizer()
    ss.update_samples(samples)
    return {'n': ss.n,
            'mean': ss.mean,
            'median': median(samples),
            'modes': mode_list(samples, bin_width),
            'variance': ss.variance,
            'range': (min(samples), max(samples)),
            'hpdi_95': get_hpd_interval(samples, 0.95),
            'qi_95': quantile_95(samples)}
         
class SampleSummarizer(object):
    count = 0
    def __init__(self, samples = None, tag = None):
        self.__class__.count += 1
        self.name = self.__class__.__name__ + '-' + str(self.count)
        self.tag = tag
        self._min = None
        self._max = None
        self._n = 0
        self._mean = 0.0
        self._sum_devs_2 = 0.0
        self._sum_devs_3 = 0.0
        self._sum_devs_4 = 0.0
        if samples:
            self.update_samples(samples)
    
    def add_sample(self, x):
        n = self._n + 1
        d = x - self._mean
        d_n = d / n
        d_n2 = d_n * d_n
        self._mean = self._mean + d_n
        first_term =  d * d_n * self._n
        self._sum_devs_4 += (first_term * d_n2 * ((n * n) - (3 * n) + 3)) + \
                (6 * d_n2 * self._sum_devs_2) - (4 * d_n * self._sum_devs_3)
        self._sum_devs_3 += (first_term * d_n * (n - 2)) - \
                (3 * d_n * self._sum_devs_2)
        self._sum_devs_2 += first_term
        self._n = n
        if not self._min:
            self._min = x
        elif x < self._min:
            self._min = x
        if not self._max:
            self._max = x
        elif x > self._max:
            self._max = x

    def update_samples(self, x_iter):
        for x in x_iter:
            self.add_sample(x)

    def _get_n(self):
        return self._n
    
    n = property(_get_n)

    def _get_min(self):
        return self._min

    def _get_max(self):
        return self._max

    minimum = property(_get_min)
    maximum = property(_get_max)
    
    def _get_mean(self):
        if self._n < 1:
            return None
        return self._mean
    
    def _get_variance(self):
        if self._n < 1:
            return None
        if self._n == 1:
            return float('inf')
        return (self._sum_devs_2 / (self._n - 1))

    def _get_std_dev(self):
        if self._n < 1:
            return None
        return math.sqrt(self._get_variance())

    def _get_pop_variance(self):
        if self._n < 1:
            return None
        return (self._sum_devs_2 / self._n)

    mean = property(_get_mean)
    variance = property(_get_variance)
    pop_variance = property(_get_pop_variance)
    std_deviation = property(_get_std_dev)

    def _get_skewness(self):
        return ((self._sum_devs_3 * math.sqrt(self._n)) / \
                (self._sum_devs_2 ** (float(3)/2)))
    def _get_kurtosis(self):
        return (((self._n * self._sum_devs_4) / (self._sum_devs_2 ** 2)) - 3)

    skewness = property(_get_skewness)
    kurtosis = property(_get_kurtosis)

    def __str__(self):
        s = StringIO()
        s.write('name = {0}\n'.format(self.name))
        s.write('sample size = {0}\n'.format(self._n))
        s.write('min = {0}\nmax = {1}\n'.format(self._min, self._max))
        s.write('mean = {0}\n'.format(self.mean))
        s.write('variance = {0}\n'.format(self.variance))
        s.write('pop variance = {0}\n'.format(self.pop_variance))
        s.write('skewness = {0}\n'.format(self.skewness))
        s.write('kurtosis = {0}\n'.format(self.kurtosis))
        return s.getvalue()


class SampleSummary(object):
    def __init__(self, sample_size = 0, mean = 0.0, variance = 0.0):
        self.n = sample_size
        self.mean = mean
        self.variance = variance
        if self.n < 1:
            self.n = 0
            self.mean = 0.0
            self.variance = 0.0

    def _get_std_dev(self):
        if self.n < 1:
            return None
        return math.sqrt(self.variance)
    std_deviation = property(_get_std_dev)

    def update(self, sample_summary):
        s1 = self
        s2 = sample_summary
        n = s1.n + s2.n
        mean = ((s1.n / float(n)) * s1.mean) + \
               ((s2.n / float(n)) * s2.mean)
        v = float(((s1.n**2) * s1.variance) + ((s2.n**2) * s2.variance) - \
                (s1.n * s1.variance) - (s1.n * s2.variance) - \
                (s2.n * s1.variance) - (s2.n * s2.variance) + \
                (s1.n * s2.n * s1.variance) + (s1.n * s2.n * s2.variance) + \
                (s1.n * s2.n * ((s1.mean - s2.mean)**2))) / \
                ((s1.n + s2.n - 1) * (s1.n + s2.n))
        self.n = n
        self.mean = mean
        self.variance = v

    def update_iter(self, sample_summaries):
        for ss in sample_summaries:
            self.update(ss)

def mean_squared_error(x, y):
    if not len(x) == len(y):
        raise ValueError('x and y must be the same length')
    sse = 0.0
    for i in range(len(x)):
        sse += (x[i] - y[i]) ** 2
    return sse / len(x)

def root_mean_square_error(x, y):
    return mean_squared_error(x, y) ** 0.5

