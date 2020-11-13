import unittest
import pandas as pd
from ProteinClassification import GOAnotation
from ProteinClassification.DefaultUse import PreconfiguredProteinClasses


class TestingGOAnnotations(unittest.TestCase):
    def test_correct_download(self):
        kinase_function_go_id = "GO:0016301"
        kinase_go_function = GOAnotation.GoMolecularFunction.from_id(kinase_function_go_id)
        self.assertEqual(kinase_function_go_id, kinase_go_function.go_id)

    def test_id_correction(self):
        replaced_id = "GO:0048365"
        with self.assertWarns(Warning):
            go_function = GOAnotation.GoMolecularFunction.from_id(replaced_id)
        self.assertEqual("GO:0031267", go_function.go_id)


class TestingProteinAnnotation(unittest.TestCase):
    def test_exemplary_proteins(self):
        protein_list = ["Q16512", "P30085", "P25774"]
        protein_mangager = PreconfiguredProteinClasses()
        for uniprot_id in protein_list:
            protein_mangager.add_protein(uniprot_id)
        self.assertEqual({"Transferase (phosphorus-containing groups)", "Transcription regulator"},
                         protein_mangager.get_matching_classes("Q16512"))
        self.assertEqual({"Transferase (phosphorus-containing groups)"},
                         protein_mangager.get_matching_classes("P30085"))
        self.assertEqual({"Peptidase"},
                         protein_mangager.get_matching_classes("P25774"))


if __name__ == '__main__':
    unittest.main()
