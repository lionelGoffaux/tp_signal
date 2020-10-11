from __future__ import print_function
from __future__ import division
from __future__ import unicode_literals


import math, cmath

def lporder (freq1, freq2, delta_p, delta_s):
    '''
    FIR lowpass filter length estimator.  freq1 and freq2 are
    normalized to the sampling frequency.  delta_p is the passband
    deviation (ripple), delta_s is the stopband deviation (ripple).
    Note, this works for high pass filters too (freq1 > freq2), but
    doesn't work well if the transition is near f == 0 or f == fs/2
    From Herrmann et al (1973), Practical design rules for optimum
    finite impulse response filters.  Bell System Technical J., 52, 769-99
    '''
    df = abs (freq2 - freq1)
    ddp = math.log10 (delta_p)
    dds = math.log10 (delta_s)

    a1 = 5.309e-3
    a2 = 7.114e-2
    a3 = -4.761e-1
    a4 = -2.66e-3
    a5 = -5.941e-1
    a6 = -4.278e-1

    b1 = 11.01217
    b2 = 0.5124401

    t1 = a1 * ddp * ddp
    t2 = a2 * ddp
    t3 = a4 * ddp * ddp
    t4 = a5 * ddp

    dinf=((t1 + t2 + a3) * dds) + (t3 + t4 + a6)
    ff = b1 + b2 * (ddp - dds)
    n = dinf / df - ff * df + 1
    return n


def remezord (fcuts, mags, devs, fsamp = 2):
    '''
    FIR order estimator (lowpass, highpass, bandpass, mulitiband).
    (n, fo, ao, w) = remezord (f, a, dev)
    (n, fo, ao, w) = remezord (f, a, dev, fs)
    (n, fo, ao, w) = remezord (f, a, dev) finds the approximate order,
    normalized frequency band edges, frequency band amplitudes, and
    weights that meet input specifications f, a, and dev, to use with
    the remez command.
    * f is a sequence of frequency band edges (between 0 and Fs/2, where
      Fs is the sampling frequency), and a is a sequence specifying the
      desired amplitude on the bands defined by f. The length of f is
      twice the length of a, minus 2. The desired function is
      piecewise constant.
    * dev is a sequence the same size as a that specifies the maximum
      allowable deviation or ripples between the frequency response
      and the desired amplitude of the output filter, for each band.
    Use remez with the resulting order n, frequency sequence fo,
    amplitude response sequence ao, and weights w to design the filter b
    which approximately meets the specifications given by remezord
    input parameters f, a, and dev:
    b = remez (n, fo, ao, w)
    (n, fo, ao, w) = remezord (f, a, dev, Fs) specifies a sampling frequency Fs.
    Fs defaults to 2 Hz, implying a Nyquist frequency of 1 Hz. You can
    therefore specify band edges scaled to a particular applications
    sampling frequency.
    In some cases remezord underestimates the order n. If the filter
    does not meet the specifications, try a higher order such as n+1
    or n+2.
    '''
    # get local copies
    fcuts = fcuts[:]
    mags = mags[:]
    devs = devs[:]

    for i in range (len (fcuts)):
        fcuts[i] = float (fcuts[i]) / fsamp

    nf = len (fcuts)
    nm = len (mags)
    nd = len (devs)
    nbands = nm

    if nm != nd:
        raise ValueError("Length of mags and devs must be equal")

    if nf != 2 * (nbands - 1):
        raise ValueError("Length of f must be 2 * len (mags) - 2")

    for i in range (len (mags)):
        if mags[i] != 0:                        # if not stopband, get relative deviation
            devs[i] = devs[i] / mags[i]

    # separate the passband and stopband edges
    f1 = fcuts[0::2]
    f2 = fcuts[1::2]

    n = 0
    min_delta = 2
    for i in range (len (f1)):
        if f2[i] - f1[i] < min_delta:
            n = i
            min_delta = f2[i] - f1[i]

    if nbands == 2:
        # lowpass or highpass case (use formula)
        l = lporder (f1[n], f2[n], devs[0], devs[1])
    else:
        # bandpass or multipass case
        # try different lowpasses and take the worst one that
        #  goes through the BP specs
        l = 0
        for i in range (1, nbands-1):
            l1 = lporder (f1[i-1], f2[i-1], devs[i], devs[i-1])
            l2 = lporder (f1[i], f2[i], devs[i], devs[i+1])
            l = max (l, l1, l2)

    n = int (math.ceil (l)) - 1               # need order, not length for remez

    # cook up remez compatible result
    ff = [0] + fcuts + [1]
    for i in range (1, len (ff) - 1):
        ff[i] *= 2

    aa = []
    for a in mags:
        aa = aa + [a, a]

    max_dev = max (devs)
    wts = [1] * len(devs)
    for i in range (len (wts)):
        wts[i] = max_dev / devs[i]

    return (n, ff, aa, wts)