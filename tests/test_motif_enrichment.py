import pytest

# Test parameters
from htsinfer.htsinfer import motif_enrichment

foreground = {"UGAUUC": 5, "UAAACC": 3, "AAGUUACCU": 1, "AAGCCUU": 1, "AGUUCUA": 1, "UUUCCCG": 5}
background = {"UGAUUC": 3, "UAAACC": 5, "AAGCCUUAU": 1, "AGUUCUA": 1, "UUUCCCG": 5, "UUGGAA": 7}
foreground_str = {"UGAUUC": "test", "UAAACC": 3, "AAGUUACCU": 1, "AAGCCUU": 1, "AGUUCUA": 1, "UUUCCCG": 5}
final_dict = {"TGATTC": [0.7333333333333334, 0.01646090534979423],
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
            motif_enrichment(foreground)

    def test_str_as_value(self):
        "Unsupported operand type 'str' as value in dictionary."
        with pytest.raises(TypeError):
            motif_enrichment(foreground_str, background)

    def test_valid_dictionaries(self):
        "Valid dictionaries passed."
        assert motif_enrichment(foreground, background) == final_dict
