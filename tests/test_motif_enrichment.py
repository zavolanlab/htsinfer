"""Unit tests for motif_enrichment.py"""

import pytest
from htsinfer.motif_enrichment import motif_enrichment

FOREGROUND_SAME = {"UGAUUC": 5, "UAAACC": 3, "AAGUUACCU": 2,
                   "AAGCCUU": 1, "AGUUCUA": 1, "UUUCCCG": 5}
BACKGROUND_SAME = {"UGAUUC": 9, "UAAACC": 1, "AAGUUACCU": 3,
                   "AAGCCUU": 1, "AGUUCUA": 4, "UUUCCCG": 3}
FINAL_SAME = {'UGAUUC': [0.72, 0.9302721574455114],
              'UAAACC': [2.4, 0.015461966703500418],
              'AAGUUACCU': [1.0, 0.0],
              'AAGCCUU': [1.1, 0.26809988836485743],
              'AGUUCUA': [0.44, 0.9053908005999413],
              'UUUCCCG': [1.65, 0.032318974952034396]}


FOREGROUND_STR = {"UGAUUC": "test", "UAAACC": 3, "AAGUUACCU": 1,
                  "AAGCCUU": 1, "AGUUCUA": 1, "UUUCCCG": 5}

FOREGROUND_DIFF = {"UGAUUC": 5, "CCUUAA": 3, "AAGUUACCU": 2,
                   "AAGCCUU": 1, "AGUUCUA": 1, "UUUCCCG": 5}
BACKGROUND_DIFF = {"UGAUUC": 9, "UAAACC": 1, "AAGUUACCU": 3,
                   "AAGCCUU": 1, "AGUUCUA": 4, "UUUAAG": 3}
FINAL_DIFF =   {'UGAUUC': [0.7090909090909091, 0.9138690964544254],
                'CCUUAA': [4.7272727272727275, 0.0008378831771956863],
                'AAGUUACCU': [1.0, 0.0],
                'AAGCCUU': [1.1, 0.26809988836485743],
                'AGUUCUA': [0.44, 0.9053908005999413],
                'UUUCCCG': [1.65, 0.032318974952034396],
                'UAAACC': [0.5909090909090909, 0.5224004421250878]} 

class TestMotifEnrichment:
    """Tests for 'motif_enrichment' function."""

    def test_no_dictionaries(self):
        "No arguments passed."
        with pytest.raises(TypeError):
            motif_enrichment()

    def test_one_dictionary(self):
        "Missing one required argument."
        with pytest.raises(TypeError):
            motif_enrichment(FOREGROUND_SAME)

    def test_str_as_value(self):
        "Unsupported operand type 'str' as value in dictionary."
        with pytest.raises(TypeError):
            motif_enrichment(FOREGROUND_STR, BACKGROUND_SAME)

    def test_valid_dictionaries_same_motifs(self):
        "Valid dictionaries passed with same motifs."
        assert motif_enrichment(FOREGROUND_SAME, BACKGROUND_SAME) == FINAL_SAME

    def test_valid_dictionaries_different_motifs(self):
        "Valid dictionaries passed."
        assert motif_enrichment(FOREGROUND_DIFF, BACKGROUND_DIFF) == FINAL_DIFF
