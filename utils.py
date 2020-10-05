import numpy as np


def pulse(n):
    pulse = np.concatenate((np.ones(1), np.zeros(n-1)))
    return pulse


def get_sinusoid(amp, freq, sample_freq, n_sample):
    n = np.arange(0, n_sample)
    return n/sample_freq, amp * np.sin(freq/sample_freq*2*np.pi*n)


def get_squarewave(amp, freq, sample_freq, n_sample):  # TODO: refactor
    x, sin = get_sinusoid(1, freq, sample_freq, n_sample)
    sin[sin >= 0] = 1
    sin[sin < 0] = -1
    return x, amp * sin

def to_db(h):
    return 20*np.log10(np.maximum(np.abs(h), 1e-5))