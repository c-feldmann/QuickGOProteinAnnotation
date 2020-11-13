import pandas as pd
from ProteinClassification.DefaultUse import PreconfiguredProteinClasses
from tqdm import tqdm

if __name__ == "__main__":
    import argparse
    protein_mangager = PreconfiguredProteinClasses()

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infile", help="name of the input file")
    parser.add_argument("-o", "--outfile", help="name of the output file")
    parser.add_argument("-c", "--column", help="header of smiles column.")
    parser.add_argument("-s", "--seperator", help="symbol for delimiation.")
    args = parser.parse_args()

    # Checking parsed arguments for validity
    if not args.infile:
        raise ValueError("Please specify an input-file via the argument '-i'!")

    if args.outfile:
        out_file = args.outfile
    else:
        out_file = "go_function_annotation.tsv"

    if not args.column:
        raise ValueError("Please specify the name of the column containing the smiles '-c'!")

    if not args.seperator:
        sep_symbol = "\t"
    else:
        sep_symbol = str(args.seperator)

    protein_df = pd.read_csv(args.infile, sep=sep_symbol, low_memory=False)
    protein_df = protein_df[args.column].unique().tolist()

    out_dict_list = []
    for uniprot_id in tqdm(protein_df):
        protein_mangager.add_protein(uniprot_id)

    for uniprot_id in protein_df:
        functions = protein_mangager.get_matching_classes(uniprot_id)
        if functions:
            for func in functions:
                out_dict_list.append({"uniprot_id": uniprot_id,
                                      "protein_function": func})

    protein_class_df = pd.DataFrame(out_dict_list)
    protein_class_df.set_index("uniprot_id", inplace=True)
    protein_class_df.to_csv(out_file, sep="\t")
    print("Listing protein functions which lead by definition to multiclass behavior. This may be fixed later.")
    protein_mangager.check_classification_overlap()
