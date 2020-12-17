import unittest

from go_protein_annotation import function_extraction
from go_protein_annotation.default_use import DefaultAnnotation

default_annotation = DefaultAnnotation()


class TestingGOAnnotations(unittest.TestCase):
    def test_correct_download(self):
        kinase_function_go_id = "GO:0016301"
        kinase_go_function = function_extraction.GoMolecularFunction.from_id(kinase_function_go_id)
        self.assertEqual(kinase_function_go_id, kinase_go_function.go_id)

    def test_id_correction(self):
        replaced_id = "GO:0048365"
        with self.assertWarns(Warning):
            go_function = function_extraction.GoMolecularFunction.from_id(replaced_id)
        self.assertEqual("GO:0031267", go_function.go_id)


class TestingProteinAnnotation(unittest.TestCase):
    def test_exemplary_proteins(self):
        self.assertEqual({"Kinase", "Transcription regulator"},
                         set(default_annotation.get_protein_functions("Q16512").functions))
        self.assertEqual({"Kinase"},
                         set(default_annotation.get_protein_functions("P30085").functions))
        self.assertEqual({"Peptidase"},
                         set(default_annotation.get_protein_functions("P25774").functions))


if __name__ == '__main__':
    unittest.main()
