# A simple scoring system for matching amino acid
def d(amino_acid_A, amino_acid_B):
  if amino_acid_A == amino_acid_B:
    return 3
    # If a(i) and aa(j) are different
  else:
    return 0
