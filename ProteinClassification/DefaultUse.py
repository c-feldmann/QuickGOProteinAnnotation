from .GOAnotation import ProteinFactory


class PreconfiguredProteinClasses(ProteinFactory):
    def __init__(self):
        """A predefined classification scheme.

        Attention: transmembrane receptor protein kinases are ill defined and always listed as "Multiclass" since
        they are receptors and transferases...

        """
        super().__init__()

        # Define default classes:
        # First argument is the given name for the class,
        # second argument is a set of GO functions which a protein must have in order to be classified as this class.
        # The third argument defines all GO functions a protein must not have in order to be classified as this class.
        # Example "Hydrolase (other)" must have "GO:0016787" (Hydrolase activity) but not one of the other specified
        # hydrolases {"GO:0016810", "GO:0016817", "GO:0016788", "GO:0016798", "GO:0008233"} to avoid overlapping
        # definitions.
        self.define_protein_class("Enzyme regulator", {"GO:0030234"}, set())

        self.define_protein_class("Peptidase", {"GO:0008233"}, set())
        self.define_protein_class("Hydrolase (C-N bonds, no peptides)", {"GO:0016810"}, set())
        self.define_protein_class("Hydrolase (acid anhydrides)", {"GO:0016817"}, set())
        self.define_protein_class("Hydrolase (ester bonds)", {"GO:0016788"}, set())
        self.define_protein_class("Hydrolase (glycosyl bonds)", {"GO:0016798"}, set())
        self.define_protein_class("Hydrolase (other)",
                                  {"GO:0016787"},
                                  {"GO:0008233", "GO:0016810", "GO:0016817", "GO:0016788", "GO:0016798"})
        self.define_protein_class("Isomerase", {"GO:0016853"}, set())

        self.define_protein_class("Ligase", {"GO:0016874"}, set())

        self.define_protein_class("Lyase", {"GO:0016829"}, set())

        self.define_protein_class("Oxidoreductase", {"GO:0016491"}, set())

        # Singaling recptor without  transmembrane signaling receptor  (GO:0004888)
        self.define_protein_class("Signaling receptor (not transmembrane)", {"GO:0038023"}, {"GO:0004888"})

        # Transmembrane signaling receptors:
        # GPCR
        self.define_protein_class("G protein-coupled receptor", {"GO:0004930"}, set())

        # Transmembrane receptor protein kinase
        self.define_protein_class("Transmembrane receptor protein kinase", {"GO:0019199"}, set())

        # Other Transmembrane signaling receptors
        self.define_protein_class("Transmembrane signaling receptors (other)", {"GO:0004888"},
                                  {"GO:0019199", "GO:0004930"})

        self.define_protein_class("Transcription regulator", {"GO:0140110"}, set())

        self.define_protein_class("Serotonin binding", {"GO:0051378"}, set())
        self.define_protein_class("Dopamine binding", {"GO:0035240"}, set())
        self.define_protein_class("Acetylcholine binding", {"GO:0042166"}, set())
        self.define_protein_class("Folic acid binding", {"GO:0005542"}, set())

        self.define_protein_class("Transferase (acyl groups)", {"GO:0016746"}, set())
        self.define_protein_class("Transferase (alkyl or aryl groups, no methyl)", {"GO:0016765"}, set())
        self.define_protein_class("Transferase (glycosyl groups)", {"GO:0016757"}, set())
        self.define_protein_class("Transferase (one-carbon groups)", {"GO:0016741"}, set())
        self.define_protein_class("Transferase (phosphorus-containing groups)", {"GO:0016772"}, set(),
                                  check_definition_clash=False)
        self.define_protein_class("Transferase (other)",
                                  {"GO:0016740"},
                                  {"GO:0016746", "GO:0016765", "GO:0016757", "GO:0016741", "GO:0016772"})

        self.define_protein_class("Transporter", {"GO:0005215"}, set())

        self.define_protein_class("Other/Unclassified",
                                  set(),
                                  {"GO:0030234", "GO:0016787", "GO:0016853", "GO:0016874", "GO:0016829",
                                   "GO:0016491", "GO:0008233", "GO:0038023", "GO:0140110", "GO:0016740",
                                   "GO:0005215",
                                   })