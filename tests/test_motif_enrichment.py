"""Unit tests for motif_enrichment.py"""

import pytest
from htsinfer.motif_enrichment import motif_enrichment

FOREGROUND_SAME = {"UGAUUC": 5, "UAAACC": 3, "AAGUUACCU": 2,
                   "AAGCCUU": 1, "AGUUCUA": 1, "UUUCCCG": 5}
BACKGROUND_SAME = {"UGAUUC": 9, "UAAACC": 1, "AAGUUACCU": 3,
                   "AAGCCUU": 1, "AGUUCUA": 4, "UUUCCCG": 3}
FINAL_SAME = {'UGAUUC': [0.6944444444444444, 0.96190821],
              'UAAACC': [3.75, 0.005024350000000011],
              'AAGUUACCU': [1.0, 0.0],
              'AAGCCUU': [1.1428571428571428, 0.21460819244384766],
              'AGUUCUA': [0.2857142857142857, 0.9375],
              'UUUCCCG': [1.9047619047619047, 0.013209342956542969]}

FOREGROUND_STR = {"UGAUUC": "test", "UAAACC": 3, "AAGUUACCU": 1,
                  "AAGCCUU": 1, "AGUUCUA": 1, "UUUCCCG": 5}

FOREGROUND_DIFF = {"UGAUUC": 5, "CCUUAA": 3, "AAGUUACCU": 2,
                   "AAGCCUU": 1, "AGUUCUA": 1, "UUUCCCG": 5}
BACKGROUND_DIFF = {"UGAUUC": 9, "UAAACC": 1, "AAGUUACCU": 3,
                   "AAGCCUU": 1, "AGUUCUA": 4, "UUUAAG": 3}
FINAL_DIFF = {'UGAUUC': [0.9583333333333334, 0.47129086020755195],
              'CCUUAA': [8.625, 2.3720532859039523e-06],
              'AAGUUACCU': [1.0, 0.0],
              'AAGCCUU': [0.8, 0.4744071960449219],
              'AGUUCUA': [0.32, 0.992154236882925],
              'UUUCCCG': [4.8, 4.0452927350997925e-05],
              'UAAACC': [0.359375, 0.7945355308633791],
              'UUUAAG': [0.23958333333333334, 0.9472567499162279]}


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
