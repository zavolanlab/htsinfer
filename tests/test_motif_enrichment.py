"""Unit tests for motif_enrichment.py"""

import pytest
from htsinfer.motif_enrichment import motif_enrichment

FOREGROUND_DICT = {"UGAUUC": 5, "UAAACC": 3, "AAGUUACCU": 2,
                   "AAGCCUU": 1, "AGUUCUA": 1, "UUUCCCG": 5}
BACKGROUND_DICT = {"UGAUUC": 3, "UAAACC": 5, "AAGCCUUAU": 1,
                   "AGUUCUA": 1, "UUUCCCG": 5, "UUGGAA": 7}
FOREGROUND_STR = {"UGAUUC": "test", "UAAACC": 3, "AAGUUACCU": 1,
                    "AAGCCUU": 1, "AGUUCUA": 1, "UUUCCCG": 5}
FINAL_DICT = {'UGAUUC': [2.4545454545454546, 0.0037400679729227004],
                'UAAACC': [1.0909090909090908, 0.28899727345086257],
                'AAGUUACCU': [3.2, 0.0009765625],
                'AAGCCUU': [1.8, 0.09059653269912638],
                'AGUUCUA': [0.9, 0.3899216801618356],
                'UUUCCCG': [0.9, 0.5592643397856017],
                'AAGCCUUAU': [0.26666666666666666, 0.984375],
                'UUGGAA': [0.20454545454545456, 0.9847514510038162]}


class Test_motif_enrichment:
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
