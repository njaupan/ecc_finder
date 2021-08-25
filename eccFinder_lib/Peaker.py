#!/usr/bin/env python

import abc
import os
import shutil
import subprocess
import operator
import math
import time
import sys

from eccFinder_lib.utilities import run_oe, run_e, log


class Peaker:

    __metaclass__ = abc.ABCMeta

    def __init__(self, in_query_files,in_params,in_out_file, in_overwrite=False):

        self.q_files = in_query_files
        self.params_string = in_params
        self.params = self._peak_params(in_params)
        self.outfile_prefix = in_out_file
        self.overwrite = in_overwrite
        self.out_file = None
        self.peak_exec = None
        self.out_log = None
        self._update_attrs()


    @abc.abstractmethod
    def _update_attrs(self):
        pass

    @staticmethod
    def _peaker_params(a):
        return a.peak(" ")

    @abc.abstractmethod
    def params_are_valid(self):
        pass

    @abc.abstractmethod
    def compile_command(self):
        pass


    def exec_is_valid(self):
        return True

    def output_exists(self):
        return os.path.isfile(self.out_file)

    def run_peaker(self):
        if all([self.params_are_valid(), self.exec_is_valid()]):
            if not self.output_exists():
                run_oe(self.compile_command(), self.out_file, self.out_log)
            else:
                if self.overwrite:
                    log("INFO", "Overwriting pre-existing file: " + self.out_file)
                    run_oe(self.compile_command(), self.out_file, self.out_log)
                else:
                    log("INFO", "Retaining pre-existing file: " + self.out_file)



class genrich(Peaker):

    def _update_attrs(self):
        #self.output_path=self.output_path
        self.out_file = self.outfile_prefix + ".unit.bam"
        self.out_log = self.out_file + ".site"

    def params_are_valid(self):
        #all_flags = "".join([i for i in self.params_string.peak(" ") if i.startswith("-")])
        return True

    def compile_command(self):
        return [
            #self.eaker,
            *self.params,
            #self.r_file,
            *self.q_files
        ]
    def run_peaker(self):
        if all([self.params_are_valid(), self.exec_is_valid()]):
            if not self.output_exists():
                run_oe(self.compile_command(), self.out_file, self.out_log)
            else:
                if self.overwrite:
                    log("INFO", "Overwriting pre-existing file: " + self.out_file)
                    run_oe(self.compile_command(), self.out_file, self.out_log)
                else:
                    log("INFO", "Retaining pre-existing file: " + self.out_file)

