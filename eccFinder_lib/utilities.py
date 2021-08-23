#!/usr/bin/env python

import subprocess
import operator
import math
import time
import sys

""" A collection of various helper functions"""

complements = str.maketrans("ACGTNURYSWKMBVDHacgtnuryswkmbvdh", "TGCANAYRSWMKVBHDtgcanayrswmkvbhd")


def get_eccFinder_version():
    return 'v1.0.0'


def reverse_complement(seq):
    return seq.translate(complements)[::-1]


def run(cmd):
    if not isinstance(cmd, list):
        raise TypeError("'run' expects a list")

    log("INFO", "Running: %s" % " ".join(cmd))
    if subprocess.call(cmd) != 0:
        raise RuntimeError("Failed : %s" % " ".join(cmd))
    log("INFO", "Finished running : %s" % " ".join(cmd))


def run_oe(cmd, out, err):
    """ Run a command and redirect stdout/stderr. """

    if not isinstance(out, str) or not isinstance(err, str):
        raise TypeError("out/err should be file names (strings)")

    f_out = open(out, "w")
    f_err = open(err, "w")

    log("INFO", "Running: %s > %s 2> %s" % (" ".join(cmd), out, err))
    if subprocess.call(cmd, stdout=f_out, stderr=f_err) != 0:
        raise RuntimeError("Failed : %s > %s 2> %s" % (" ".join(cmd), out, err))

    log("INFO", "Finished running : %s > %s 2> %s" % (" ".join(cmd), out, err))

    f_out.close()
    f_err.close()

def run_e(cmd, err):
    """ Run a command and redirect stderr but not stdout. """

    if not isinstance(err, str):
        raise TypeError("err should be a file name (string)")

    f_err = open(err, "w")

    log("INFO", "Running: %s 2> %s" % (" ".join(cmd), err))
    if subprocess.call(cmd, stderr=f_err) != 0:
        raise RuntimeError('Failed : %s 2> %s' % (" ".join(cmd), err))

    log("INFO", "Finished running : %s > %s" % (" ".join(cmd), err))

    f_err.close()

def run_oae(cmd, out, err):
    """ Run a command and redirect stdout/stderr. Append rather than write to err"""

    if not isinstance(out, str) or not isinstance(err, str):
        raise TypeError("out/err should be file names (strings)")

    f_out = open(out, "w")
    f_err = open(err, "a")

    log("INFO", "Running: %s > %s 2> %s" % (" ".join(cmd), out, err))
    if subprocess.call(cmd, stdout=f_out, stderr=f_err) != 0:
        raise RuntimeError('Failed : %s > %s 2> %s. Check stderr file for details.' % (" ".join(cmd), out, err))

    log("INFO", "Finished running : %s > %s 2> %s" % (" ".join(cmd), out, err))

    f_out.close()
    f_err.close()


def log(level, message):
    """ Log messages to standard error. """
    level = level.upper()
    if level not in {"VERSION", "CMD", "INFO", "WARNING", "DEBUG"}:
        raise ValueError("Invalid logging level: {}".format(level))

    sys.stderr.write(time.ctime() + " --- " + level + ": " + message + "\n")
    sys.stderr.flush()

def p2q(p):
    """ Convert Pr(incorrect) to MAPQ. """
    return round(-10 * math.log(p, 10))


def q2p(q):
    """ Convert MAPQ to Pr(incorrect). """
    return 10**(q/-10)
