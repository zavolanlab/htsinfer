"""Unit tests for motif_enrichment.py"""

import pytest
from htsinfer.htsinfer import motif_enrichment

FOREGROUND_DICT = {"UGAUUC": 5, "UAAACC": 3, "AAGUUACCU": 1, "AAGCCUU": 1, "AGUUCUA": 1,
                   "UUUCCCG": 5}
BACKGROUND_DICT = {"UGAUUC": 3, "UAAACC": 5, "AAGCCUUAU": 1, "AGUUCUA": 1, "UUUCCCG": 5,
                   "UUGGAA": 7}
FOREGROUND_STR = {"UGAUUC": "test", "UAAACC": 3, "AAGUUACCU": 1, "AAGCCUU": 1, "AGUUCUA": 1,
                  "UUUCCCG": 5}
FINAL_DICT = {"TGATTC": [0.7333333333333334, 0.01646090534979423],
              "TAAACC": [0.7333333333333334, 0.21947873799725662],
              "AAGTTACCT": [1.375, 0.26337448559670795],
              "AAGCCTT": [1.6041666666666667, 0.26337448559670795],
              "AGTTCTA": [1.6041666666666667, 0.26337448559670795],
              "TTTCCCG": [1.6041666666666667, 0.01646090534979423]}


class TestMotifEnrichment:
    """Tests for 'motif_enrichment' function."""

    def test_no_dictionaries(self):
        "No arguments passed."
        with pytest.raises(TypeError):
            motif_enrichment()

    def test_one_dictionary(self):
        "Missing one required argument."
        with pytest.raises(TypeError):
            motif_enrichment(FOREGROUND_DICT)

    def test_str_as_value(self):
        "Unsupported operand type 'str' as value in dictionary."
        with pytest.raises(TypeError):
            motif_enrichment(FOREGROUND_STR, BACKGROUND_DICT)

    def test_valid_dictionaries(self):
        "Valid dictionaries passed."
        assert motif_enrichment(FOREGROUND_DICT, BACKGROUND_DICT) == FINAL_DICT
