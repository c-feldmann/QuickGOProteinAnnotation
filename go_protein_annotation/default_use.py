from .function_extraction import SpecialFunctionAnnotation

# Define default classes:
# First argument is a set of GO functions which a protein must have in order to be classified as this class.
# The second argument defines all GO functions a protein must not have in order to be classified as this class.
# Third argument is the given name for the class,
# Example "Hydrolase (other)" must have "GO:0016787" (Hydrolase activity) but not one of the other specified
# hydrolases {"GO:0016810", "GO:0016817", "GO:0016788", "GO:0016798", "GO:0008233"}

default_functions = [({"GO:0030234"}, set(), "Enzyme regulator"),
                     ({"GO:0008233"}, set(), "Peptidase"),
                     ({"GO:0016810"}, set(), "Hydrolase (C-N bonds, no peptides)"),
                     ({"GO:0016817"}, set(), "Hydrolase (acid anhydrides)"),
                     ({"GO:0016788"}, set(), "Hydrolase (ester bonds)"),
                     ({"GO:0016798"}, set(), "Hydrolase (glycosyl bonds)"),
                     ({"GO:0016787"},
                      {"GO:0008233", "GO:0016810", "GO:0016817", "GO:0016788", "GO:0016798"},
                      "Hydrolase (other)"
                      ),

                     ({"GO:0016853"}, set(), "Isomerase"),
                     ({"GO:0016874"}, set(), "Ligase"),
                     ({"GO:0016829"}, set(), "Lyase"),
                     ({"GO:0016491"}, set(), "Oxidoreductase"),

                     # Signaling receptor
                     # Non-transmembrane signaling receptor
                     ({"GO:0038023"}, {"GO:0004888"}, "Signaling receptor (not transmembrane)"),

                     # Transmembrane signaling receptors:
                     ({"GO:0004930"}, set(), "G protein-coupled receptor"),
                     ({"GO:0019199"}, set(), "Transmembrane receptor protein kinase"),
                     ({"GO:0004888"},
                      {"GO:0019199", "GO:0004930"},
                      "Transmembrane signaling receptors (other)"
                      ),

                     ({"GO:0140110"}, set(), "Transcription regulator"),
                     ({"GO:0051378"}, set(), "Serotonin binding"),
                     ({"GO:0035240"}, set(), "Dopamine binding"),
                     ({"GO:0042166"}, set(), "Acetylcholine binding"),
                     ({"GO:0005542"}, set(), "Folic acid binding"),
                     ({"GO:0051379"}, set(), "Epinephrine binding"),
                     ({"GO:0016746"}, set(), "Transferase (acyl groups)"),
                     ({"GO:0016765"}, set(), "Transferase (alkyl or aryl groups, no methyl)"),
                     ({"GO:0016757"}, set(), "Transferase (glycosyl groups)"),
                     ({"GO:0016741"}, set(), "Transferase (one-carbon groups)"),
                     ({"GO:0016301"}, set(), "Kinase"),
                     ({"GO:0016772"}, {"GO:0016301"}, "Transferase (phosphorus-containing groups, non-kinase)"),
                     ({"GO:0016740"},
                      {"GO:0016746", "GO:0016765", "GO:0016757", "GO:0016741", "GO:0016772"},
                      "Transferase (other)"),
                     ({"GO:0005215"}, set(), "Transporter"),

                     ({"GO:0003677"}, set(), "DNA binding"),
                     ]


class DefaultAnnotation(SpecialFunctionAnnotation):
    def __init__(self):
        super(SpecialFunctionAnnotation, self).__init__()
        self._function_specifications = default_functions
