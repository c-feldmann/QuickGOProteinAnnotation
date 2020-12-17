import pandas as pd
from go_protein_annotation.default_use import DefaultAnnotation
from tqdm import tqdm

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infile", help="name of the input file")
    parser.add_argument("-o", "--outfile", help="name of the output file")
    parser.add_argument("-c", "--column", help="header of smiles column.")
    parser.add_argument("-s", "--separator", help="symbol for delimitation. Write 'tab' for tab-delimited files")
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
    if not args.separator or args.separator == "tab":
        sep_symbol = "\t"
    else:
        sep_symbol = str(args.separator)
        print(sep_symbol)
    protein_df = pd.read_csv(args.infile, sep=sep_symbol, low_memory=False)
    protein_df = protein_df[args.column].unique().tolist()
    default_annotation = DefaultAnnotation()
    out_dict_list = []
    for uniprot_id in tqdm(protein_df):
        out = default_annotation.get_protein_functions(uniprot_id, as_dataframe=False)
        out_dict_list.extend(out)

    protein_class_df = pd.DataFrame(out_dict_list)
    protein_class_df.set_index("uniprot_id", inplace=True)
    protein_class_df.to_csv(out_file, sep="\t")
