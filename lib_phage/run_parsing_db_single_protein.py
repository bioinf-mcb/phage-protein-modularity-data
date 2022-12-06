import sys
import pandas as pd

if __name__ == '__main__':
    protein_hhr_path = sys.argv[1]
    ecf = bool(sys.argv[2])

    sys.path.append(sys.argv[3])
    from lib_phage.utils import parse_hhr_single_file

    # parse single file
    protein_hhr_parsed = parse_hhr_single_file(protein_hhr_path, ecf)

    # save df after parsing
    protein_hhr_parsed.to_csv(protein_hhr_path.replace('.hhr', '.txt'),
                             index=False)
