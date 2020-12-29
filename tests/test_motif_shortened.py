# -*- coding: utf-8 -*-
"""
Created on Mon Nov  2 01:28:58 2020
"""

import unittest
from htsinfer.read_motif_v6 import find_overlaps


class Test_Inputs(unittest.TestCase):

    def test_call(self):
        # 0 arguments, not acceptable
        with self.assertRaises(TypeError):
            find_overlaps()
        # 1 argument, not acceptable
        with self.assertRaises(TypeError):
            find_overlaps("")
        # 2 arguments, not acceptable
        with self.assertRaises(TypeError):
            find_overlaps("", "")
        # 3 or 4 or arguments, acceptable
        find_overlaps("A", "A", 1)
        find_overlaps("A", "A", 1, True)
        # 5 arguments, not acceptable
        with self.assertRaises(TypeError):
            find_overlaps("", "", 0, True, 0)

    def test_positional_arguments(self):
        # Check that argument 1 is of type string
        with self.assertRaises(TypeError):
            find_overlaps(0, "", 1, True)
        # Check that argument 2 is of type string
        with self.assertRaises(TypeError):
            find_overlaps("", 0, 1, True)
        # Check that argument 3 is of type int
        with self.assertRaises(TypeError):
            find_overlaps("", "", "", True)
        # Check that argument 4 is of type bool
        with self.assertRaises(TypeError):
            find_overlaps("", "", 1, "")

    def test_arguments_types(self):
        # Check that motif is of type string
        with self.assertRaises(TypeError):
            find_overlaps(motif=0, read="", min_overlap=1, full_contain=False)
        # Check that read is of type string
        with self.assertRaises(TypeError):
            find_overlaps(motif="", read=0, min_overlap=1, full_contain=False)
        # Check that min_overlap is of type int
        with self.assertRaises(TypeError):
            find_overlaps(motif="", read="",
                          min_overlap="", full_contain=False)
        # Check that use_n is of type bool
        with self.assertRaises(TypeError):
            find_overlaps(motif="", read="", min_overlap=1, full_contain=6)

    def test_argument_range(self):
        # Check that motif is not accepted if it is an empty string
        with self.assertRaises(ValueError):
            find_overlaps(motif="", read="A",
                          min_overlap=1, full_contain=False)
        # Check that read is not accepted if it is an empty string
        with self.assertRaises(ValueError):
            find_overlaps(motif="A", read="",
                          min_overlap=1, full_contain=False)
        # Check that min_overlap is not accepted if smaller than 1
        with self.assertRaises(ValueError):
            find_overlaps(motif="A", read="A",
                          min_overlap=0, full_contain=False)
        # check if read is longer than motif
        with self.assertRaises(ValueError):
            find_overlaps(motif="AAAAAA", read="A",
                          min_overlap=1, full_contain=False)
        # check if motif contains small character
        with self.assertRaises(ValueError):
            find_overlaps(motif="a", read="AAAAAA",
                          min_overlap=1, full_contain=False)
        # check if read contains small character
        with self.assertRaises(ValueError):
            find_overlaps(motif="A", read="Aa",
                          min_overlap=1, full_contain=False)

    def test_return_value(self):
        rv = find_overlaps(motif="G", read="AAAA",
                           min_overlap=1, full_contain=False)
        self.assertTrue(isinstance(rv, list))


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

    def test_return_value_full(self):
        rv = find_overlaps(motif="AAA", read="AAAA",
                           min_overlap=1, full_contain=True)
        self.assertTrue(rv[0] == (0, 1))


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
