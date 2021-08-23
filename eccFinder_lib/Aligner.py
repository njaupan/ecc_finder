#!/usr/bin/env python

import abc
import os
import shutil

from eccFinder_lib.utilities import run_oe, run_e, log


class Aligner:

    __metaclass__ = abc.ABCMeta

    def __init__(self, in_ref_file, in_query_files, in_aligner, in_params, in_out_file, in_overwrite=False):
        self.r_file = in_ref_file
        self.q_files = in_query_files
        self.aligner = in_aligner
        self.params_string = in_params
        self.params = self._split_params(in_params)
        self.outfile_prefix = in_out_file
        self.overwrite = in_overwrite
        self.out_file = None
        self.aligner_exec = None
        self.out_log = None
        self._update_attrs()

    @abc.abstractmethod
    def _update_attrs(self):
        pass

    @staticmethod
    def _split_params(a):
        return a.split(" ")

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

    def run_aligner(self):
        if all([self.params_are_valid(), self.exec_is_valid()]):
            if not self.output_exists():
                run_oe(self.compile_command(), self.out_file, self.out_log)
            else:
                if self.overwrite:
                    log("INFO", "Overwriting pre-existing file: " + self.out_file)
                    run_oe(self.compile_command(), self.out_file, self.out_log)
                else:
                    log("INFO", "Retaining pre-existing file: " + self.out_file)


class bwaAligner(Aligner):

    def _update_attrs(self):
        """ Update class attributes for a specific aligner. """
        self.aligner_exec = "bwa mem"
        self.out_file = self.outfile_prefix + ".sam"
        self.out_log = self.out_file + ".log"

    def params_are_valid(self):
        if "-p" in self.params_string:
            raise ValueError("eccFinder names its own alignment files.")

        return True

    def compile_command(self):
        return [
            self.aligner,
            *self.params,
            self.outfile_prefix,
            self.r_file,
            *self.q_files
        ]

    def run_aligner(self):
        """ Run the aligner. """
        if all([self.params_are_valid(), self.exec_is_valid()]):
            if not self.output_exists():
                run_e(self.compile_command(), self.out_log)
            else:
                if self.overwrite:
                    log("INFO", "Overwriting pre-existing file: " + self.out_file)
                    run_e(self.compile_command(), self.out_log)
                else:
                    log("INFO", "Retaining pre-existing file: " + self.out_file)

class Minimap2Aligner(Aligner):

    def _update_attrs(self):
        self.aligner_exec = "minimap2"
        self.out_file = self.outfile_prefix + ".paf"
        self.out_log = self.out_file + ".log"

    def params_are_valid(self):
        all_flags = "".join([i for i in self.params_string.split(" ") if i.startswith("-")])
        if "a" in all_flags:
            raise ValueError("Alignments must not be in SAM format (-a).")

        if "c" in all_flags:
            log("WARNING", "Computing base-alignments (-c) will slow down Minimap2 alignment.")

        return True

    def compile_command(self):
        return [
            self.aligner,
            *self.params,
            self.r_file,
            *self.q_files
        ]

class Minimap2SAMAligner(Aligner):

    def _update_attrs(self):
        self.aligner_exec = "minimap2"
        self.out_file = self.outfile_prefix + ".sam"
        self.out_log = self.out_file + ".log"

    def params_are_valid(self):
        all_flags = "".join([i for i in self.params_string.split(" ") if i.startswith("-")])
        if "a" not in all_flags:
            raise ValueError("Alignments must be in SAM format (-a).")
        return True

    def compile_command(self):
        return [
            self.aligner,
            *self.params,
            self.r_file,
            *self.q_files
        ]