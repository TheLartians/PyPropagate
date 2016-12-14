
# TODO: convert to c routines, to reduce N^2 memory requirement

cache={}

def hankel(f, xmax=None, n=0, kmax=None):
    '''
     As in: Theory and operational rules for the discreteHankel transform
     by Natalie Baddour* and Ugo Chouinard
     https://www.researchgate.net/publication/272942976_Theory_and_Operational_Rules_for_the_Discrete_Hankel_Transform#pf3
    '''

    import scipy as sci
    import scipy.special
    import numpy as np

    N = len(f) + 1

    def create_T():
        return T

    jkey = ('j', n)
    if jkey not in cache or len(cache[jkey]) < N:
        cache[jkey] = np.array(sci.special.jn_zeros(n, N + 1))

    jn = cache[jkey]

    Ykey = ('Y', n, N)
    if Ykey not in cache:
        k, m = np.meshgrid(np.arange(N - 1), np.arange(N - 1))
        jN = jn[N - 1]
        # T = 2 * sci.special.jn(n,jn[m]*jn[k]/jn[N-1]) / (sci.special.jn(n+1,jn[m]) * sci.special.jn(n+1,jn[k]) * jn[N-1])
        Y = 2 / (jN * sci.special.jn(n + 1, jn[k]) ** 2) * sci.special.jn(n, jn[m] * jn[k] / jN)
        cache[Ykey] = Y

    if xmax is not None:
        if kmax is not None:
            raise ValueError('provide either xmax an kmax')
    else:
        if kmax is not None:
            xmax = jn[N - 1] / kmax
        else:
            xmax = 1
    return np.dot(cache[Ykey], f) * xmax ** 2 / jn[N - 1]


def hankel_samples(N, xmax=1, n=0):
    import scipy as sci
    import scipy.special
    import numpy as np

    jn = np.array(sci.special.jn_zeros(n, N + 1))
    return jn[:-1] * xmax / jn[N]


def hankel_freq(N, xmax=1, n=0):
    import scipy as sci
    import scipy.special
    import numpy as np

    jn = np.array(sci.special.jn_zeros(n, N))
    return jn / xmax


def hankel_resample_matrix(N1, new, xmax=None, kmax=None, n=0, cache_key=None):
    '''
    As in: Reconstruction of optical fields with the Quasi-discrete Hankel transform
    By: Andrew W. Norfolk and Edward J. Grace
    https://www.osapublishing.org/DirectPDFAccess/31AA535C-9410-DDAC-67361BE3117DB160_199359/oe-18-10-10551.pdf?da=1&id=199359&seq=0&mobile=no
    '''

    if cache_key is not None:
        if cache_key in cache:
            return cache[cache_key]

    import scipy as sci
    import scipy.special
    import numpy as np

    N2 = len(new)

    jn = np.array(sci.special.jn_zeros(n, N1 + 1))
    jN = jn[N1]

    if xmax is not None:
        if kmax is not None:
            raise ValueError('provide either xmax an kmax')
        kmax = jN / xmax
    else:
        if kmax is not None:
            xmax = jN / kmax
        else:
            xmax = 1
            kmax = jN / xmax

    samples = jn[:-1] / kmax

    K = kmax

    k, m = np.meshgrid(np.arange(N1), np.arange(N2))

    samenm = samples[k] == abs(new[m])
    samen = ((samenm).sum(axis=1) > 0)
    same = samen[m]

    S = 2 * samples[k] * sci.special.jn(n, K * new[m]) / (
    K * (samples[k] ** 2 - new[m] ** 2 + samenm) * sci.special.jn(n + 1, K * samples[k]))
    S = S * (same == False) + same * samenm

    if cache_key is not None:
        cache[('S',cache_key)] = S

    return S


def hankel_resample(f, new, xmax=None, kmax=None, n=0):
    import numpy as np
    S = hankel_resample_matrix(len(f), new, xmax, kmax, n)
    return np.dot(S, f)