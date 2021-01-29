'''a few simple test cases for the find_best_interval function.

run with: python -m pytest test_best_intervals.py
'''
import numpy as np
import pytest

from best_interval import find_best_interval


def _test_find_best_interval(vals, threshold, expected_start, expected_stop):
    '''runs find_best_interval on the given inputs and compares to the
    expected result.
    '''
    start, stop = find_best_interval(vals, threshold)
    assert start == expected_start
    assert stop == expected_stop


def test_all_flat_pos():
    num_samples = 10
    vals = np.ones(num_samples)
    _test_find_best_interval(vals,
                             0.0,
                             expected_start=0,
                             expected_stop=num_samples)


def test_all_flat_neg():
    num_samples = 10
    vals = -np.ones(num_samples)
    _test_find_best_interval(vals, 0.0, expected_start=0, expected_stop=0)


def test_simple_interval():
    vals = np.asarray([-1, -1, 1, 1, 1, -1])
    _test_find_best_interval(vals, 0.0, expected_start=2, expected_stop=5)


def test_simple_rhs():
    vals = np.asarray([-1, -1, -1, 1, 1])
    _test_find_best_interval(vals,
                             0.0,
                             expected_start=3,
                             expected_stop=len(vals))


def test_simple_lhs():
    vals = np.asarray([1, 1, -1, -1, -1])
    _test_find_best_interval(vals, 0.0, expected_start=0, expected_stop=2)


def test_small_gap():
    vals = np.asarray([-1, 1, 1, -1, 1, 1, -1, -1])
    _test_find_best_interval(vals, 0.0, expected_start=1, expected_stop=6)


def test_large_gap():
    vals = np.asarray([-1, 1, -1, -1, 1, 1, -1, -1])
    _test_find_best_interval(vals, 0.0, expected_start=4, expected_stop=6)


@pytest.mark.parametrize('threshold', [-7.0, -0.0001, 13.3, 999])
def test_threshold(threshold):
    vals = np.asarray([-1, 1, -1, -1, 1, 1, -1, -1]) + threshold
    _test_find_best_interval(vals,
                             threshold,
                             expected_start=4,
                             expected_stop=6)
