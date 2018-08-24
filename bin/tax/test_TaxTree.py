


#Sample tree.
#9606	|	Homo sapiens	|		|	Homo	|	Hominidae	|	Primates	|	Mammalia	|	Chordata	|	Metazoa	|	Eukaryota	|
#63221	|	Homo sapiens neanderthalensis	|	Homo sapiens	|	Homo	|	Hominidae	|	Primates	|	Mammalia	|	Chordata	|	Metazoa	|	Eukaryota	|
#741158	|	Homo sapiens subsp. 'Denisova'	|	Homo sapiens	|	Homo	|	Hominidae	|	Primates	|	Mammalia	|	Chordata	|	Metazoa	|	Eukaryota	|
#11676	|	Human immunodeficiency virus 1	|		|	Lentivirus	|	Retroviridae	|		|		|		|		|	Viruses	|
#7157	|	Culicidae	|		|		|		|	Diptera	|	Insecta	|	Arthropoda	|	Metazoa	|	Eukaryota	|
#842271	|	Lepidoptera sp. BOLD:AAE5441	|		|		|		|	Lepidoptera	|	Insecta	|	Arthropoda	|	Metazoa	|	Eukaryota	|
#1350382	|	Scenedesmus sp. LUCC 015	|		|	Scenedesmus	|	Scenedesmaceae	|	Sphaeropleales	|	Chlorophyceae	|	Chlorophyta	|	Viridiplantae	|	Eukaryota	|
#1289305	|	Chlorophyceae sp. WJT24VFNP4	|		|		|		|		|	Chlorophyceae	|	Chlorophyta	|	Viridiplantae	|	Eukaryota	|
#1427098	|	Pseudomonas sp. 161(2013)	|		|	Pseudomonas	|	Pseudomonadaceae	|	Pseudomonadales	|	Gammaproteobacteria	|	Proteobacteria	|		|	Bacteria	|
#1042	|	Erythrobacter sp.	|		|	Erythrobacter	|	Erythrobacteraceae	|	Sphingomonadales	|	Alphaproteobacteria	|	Proteobacteria	|		|	Bacteria	|

#                       / Chordata  Mammalia  Primates  Hominidae  Homo - Homo sapiens
#              / Metazoa
#             |        |                     / Diptera - Culicidae
#             |         \ Arthropoda  Insecta
#   / Eukaryota                              \ Lepidoptera - Lepidoptera sp. BOLD:AAE5441
#  |          |
#  |          |                                           / Sphaeropleales  Scenedesmaceae  Scenedesmus - Scenedesmus sp. LUCC 015
#  |           \ Viridiplantae  Chlorophyta  Chlorophyceae
#--|                                                      \- Chlorophyceae sp. WJT24VFNP4
#  |
#  |- Viruses  Retroviridae  Lentivirus - Human immunodeficiency virus 1
#  |
#  |                          / Gammaproteobacteria  Pseudomonadales  Pseudomonadaceae  Pseudomonas - Pseudomonas sp. 161(2013)
#   \ Bacteria  Proteobacteria
#                             \ Alphaproteobacteria  Sphingomonadales  Erythrobacteraceae  Erythrobacter - Erythrobacter sp.

from bin.tax import TaxTree
import unittest
import pickle
import io
from contextlib import redirect_stdout

species_indexes = ["9606","63221","741158","11676","7157","842271","1350382","1289305","1427098","1042"]

class TestTaxTree(unittest.TestCase):

    def test_extract_data(self):
        test_lineage = "1805293	|	Candidatus Pacearchaeota archaeon CG1_02_30_18	|		|		|		|		|		|		|		|	Archaea	|"
        test_list = ['1805293 ', ' Candidatus Pacearchaeota archaeon CG1_02_30_18 ', ' ', ' ', ' ', ' ', ' ', ' ',
                     ' ', ' Archaea ']

        result_list = TaxTree.extract_data(test_lineage)
        self.assertEqual(result_list, test_list)

    def test_lineage_fix(self):

        original_lineage = [['Name', ' Human immunodeficiency virus 1 '], ['Species', ' '], ['Genus', ' Lentivirus '],
                        ['Family', ' Retroviridae '], ['Order', ' '], ['Class', ' '], ['Phylum', ' '], ['Kingdom', ' '],
                        ['Superkingdom', ' Viruses ']]
        fixed_lineage = [['Name', ' Human immunodeficiency virus 1 '], ['Species', ' Human immunodeficiency virus 1 '],
                         ['Genus', ' Lentivirus '], ['Family', ' Retroviridae '], ['Superkingdom', ' Viruses ']]

        test_lineage = TaxTree.fixLineage(original_lineage)

        self.assertEqual(test_lineage,fixed_lineage)

    def test_most(self):
        with open('test_data/most_example.p', 'rb') as file:
            most_test = pickle.load(file)


        TaxTree.ranked_lineage_filepath = '../../bin/database/taxonomy/rankedlineage.dmp'

        rootNode = TaxTree.TaxNode("ROOT")
        TaxTree.loadRankedLineage(species_indexes,rootNode)
        rootNode.countChildren()

        ranks = ['Species', 'Genus', 'Family', 'Order',
                 'Class', 'Phylum', 'Kingdom', 'Superkingdom', 'All'] 

        f = io.StringIO()
        with redirect_stdout(f):
            for r in ranks:
                rootNode.calculateCommonnessByRank(r,5,rootNode)
        output = f.getvalue()

        self.assertEqual(most_test,output)

    def test_least(self):
        with open('test_data/least_example.p', 'rb') as file:
            least_test = pickle.load(file)

        TaxTree.ranked_lineage_filepath = '../../bin/database/taxonomy/rankedlineage.dmp'

        rootNode = TaxTree.TaxNode("ROOT")
        TaxTree.loadRankedLineage(species_indexes, rootNode)
        rootNode.countChildren()

        ranks = ['Species', 'Genus', 'Family', 'Order',
                 'Class', 'Phylum', 'Kingdom', 'Superkingdom', 'All']

        f = io.StringIO()
        with redirect_stdout(f):
            for r in ranks:
                rootNode.calculateCommonnessByRank(r, 5, rootNode, False)
        output = f.getvalue()

        self.assertEqual(least_test, output)

    def test_ETETree(self):

        with open('test_data/ETETree_example.p', 'rb') as file:
            ETETree_example = pickle.load(file)
        TaxTree.ranked_lineage_filepath = '../../bin/database/taxonomy/rankedlineage.dmp'
        rootNode = TaxTree.TaxNode("ROOT")
        TaxTree.loadRankedLineage(species_indexes, rootNode)
        rootNode.countChildren()
        ETETree = TaxTree.create_ETE_tree(rootNode)

        f = io.StringIO()
        with redirect_stdout(f):
            ETETree.describe()
        ETE_describe = f.getvalue()

        t = io.StringIO()
        with redirect_stdout(t):
            ETETree_example.describe()
        ETE_test_describe = t.getvalue()

        self.assertEqual(ETE_describe,ETE_test_describe)


    def test_tree_integrity_ascii(self):
        with open('test_data/tree_example.p', 'rb') as file:
            test_tree = pickle.load(file)


        TaxTree.ranked_lineage_filepath = '../../bin/database/taxonomy/rankedlineage.dmp'

        rootNode = TaxTree.TaxNode("ROOT")
        TaxTree.loadRankedLineage(species_indexes,rootNode)
        rootNode.countChildren()
        ETETree = TaxTree.create_ETE_tree(rootNode)
        ascii = TaxTree.get_ASCII_tree(ETETree)

        self.assertEqual(ascii, test_tree)