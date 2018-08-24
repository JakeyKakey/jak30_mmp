from bin.tax import TaxMatch
import unittest

class TestTaxMatch(unittest.TestCase):

    def test_name_cleaning(self):
        string = ["gi|551368999|emb|BX284682.8|=Zebrafish, DNA sequence from. clone DKEYP-74C5; in linkage group_22, " \
                 "complete sequence scaffold chromosome DNA group linkage assembly strain " \
                  " PREDICTED: hypothetical partial protein from transcript variant clone in genome "]

        expected = ["gi551368999embBX284682 8 Zebrafish   DKEYP-74C5  22       "]

        actual = TaxMatch.format_results(string)

        self.assertEqual(actual,expected)