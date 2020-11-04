# -*- coding: utf-8 -*-
"""
Created on Mon Nov  2 01:28:58 2020
"""

import unittest
from read_motif_v6 import find_overlaps


class TestMatchFull(unittest.TestCase):

    def test_single_match(self):
        rv = find_overlaps(motif="GGA", read="TACGGGACGAT",
                           min_overlap=1, full_contain=False)
        self.assertTrue(len(rv) == 1)
        self.assertTrue(rv[0] == (4, 1))

    def test_single_match_start(self):
        rv = find_overlaps(motif="ACGGG", read="ACGGGACGAT",
                           min_overlap=1, full_contain=False)
        self.assertTrue(len(rv) == 1)
        self.assertTrue(rv[0] == (0, 1))

    def test_single_match_end(self):
        rv = find_overlaps(motif="CGA", read="TACGGGACGA",
                           min_overlap=1, full_contain=False)
        self.assertTrue(len(rv) == 1)
        self.assertTrue(rv[0] == (7, 1))

    def test_multi_match_not_overlapping(self):
        rv = find_overlaps(motif="CGA", read="TATTCGATTAGCGAAT",
                           min_overlap=1, full_contain=False)
        self.assertTrue(len(rv) == 2)
        self.assertTrue(rv[0] == (4, 1))
        self.assertTrue(rv[1] == (11, 1))

    def test_multi_match_overlapping(self):
        rv = find_overlaps(motif="CGACGA", read="TATTCGACGACGATTAGCGAAT",
                           min_overlap=1, full_contain=False)
        self.assertTrue(len(rv) == 2)
        self.assertTrue(rv[0] == (4, 1))
        self.assertTrue(rv[1] == (7, 1))


class TestMatchPartial(unittest.TestCase):

    def test_match_start(self):
        rv = find_overlaps(motif="GTA", read="TACGGGACGA",
                           min_overlap=2, full_contain=False)
        self.assertTrue(len(rv) == 1)
        self.assertTrue(rv[0] == (0, 2/3))

    def test_match_end(self):
        rv = find_overlaps(motif="GTAAA", read="TACGGGACGAGT",
                           min_overlap=2, full_contain=False)
        self.assertTrue(len(rv) == 1)
        self.assertTrue(rv[0] == (10, 2/5))


class TestMatchMixed(unittest.TestCase):

    def test_multi_match_start(self):
        rv = find_overlaps(motif="GTA", read="TACGGGTAGA",
                           min_overlap=2, full_contain=False)
        self.assertTrue(len(rv) == 2)
        self.assertTrue(rv[0] == (0, 2/3))
        self.assertTrue(rv[1] == (5, 1))

    def test_multi_match_end(self):
        rv = find_overlaps(motif="GTAAA", read="ACGGTAAAAGT",
                           min_overlap=2, full_contain=False)
        self.assertTrue(len(rv) == 2)
        self.assertTrue(rv[0] == (3, 1))
        self.assertTrue(rv[1] == (9, 2/5))
# further to implement. to give stress and see if it works well.


class NoMatch(unittest.TestCase):
    def test_noMatch(self):
        rv = find_overlaps(motif="GTA", read="AAAAAAAAAAA",
                           min_overlap=2, full_contain=False)
        self.assertTrue(len(rv) == 0)


class TestLongReads(unittest.TestCase):

    def test_file_1(self):
        pass


if __name__ == '__main__':
    unittest.main()
