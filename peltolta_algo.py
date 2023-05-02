import math
# Score system from matching amino acids
import score

# Codon to amino acid 
codonDict = {
    'UUU': 'F', 'UUC': 'F', #Phe
    'UUA': 'L', 'UUG':'L', 'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG':'L', #Leu
    'AUU': 'I', 'AUC': 'I', 'AUA': 'I', #Ile
    'AUG':'M', #Met
    'GUU':'V', 'GUC':'V', 'GUA':'V', 'GUG':'V', #Val
    'UCU':'S', 'UCC':'S', 'UCA': 'S', 'UCG': 'S','AGU':'S','AGC':'S',  #Ser
    'CCU': 'P', 'CCC': 'P', 'CCA':'P', 'CCG': 'P', #Pro
    'ACU': 'T', 'ACC':'T', 'ACA':'T', 'ACG':'T',  #Thr
    'GCU': 'A', 'GCC':'A','GCA':'A','GCG':'A',  #Ala
    'UAU': 'Y', 'UAC':'Y', #Tyr
    'CAU':'H', 'CAC':'H', #His
    'CAA':'Q', 'CAG':'Q', #Gln
    'AAA':'K', 'AAG':'K', #Lys
    'GAU':'D', 'GAC': 'D', #Asp
    'GAA':'E', 'GAG':'E', #Glu
    'UGU':'C', 'UGC':'C', #Cys
    'UGG': 'W', #Trp
    'CGU': 'R', 'GCG':'R','CGA':'R','CGG':'R','AGA':'R','AGG':'R', #Arg
    'GGU': 'G', 'GGC':'G', 'GGA':'G', 'GGG':'G', #Gly
    'UGA':'*', 'UAA':'*', 'UAG':'*' # Stop Codon
}

# Function to translate codon to amino acid
def aa(j):
  codon = b[j-2:j+1]
  # Check if it is a DNA Sequence
  if 'T' in codon:
    for i in range(len(codon)):
      if codon[i] == 'T':
        codon = codon.replace('T','U')
  return codonDict[codon]

def E(i,j):
  _del1 = 2 #del -1
  _del2 = 4 # del -2
  del1 = 2 # del 1
  del2 = 4 #del 2

  if(i > 0 and j > 4):
    matrix_E[i][j] = score.d(a[i-1], aa(j-3)) + min(
        E(i-1,j-1) + _del2,
        E(i-1,j-2) + _del1,
        E(i-1,j-3),
        E(i-1,j-4) + del1,
        E(i-1,j-5) + del2
    )
  return matrix_E[i][j]

def main():
  global a,b
  # Amino acid sequence
  a = 'LV'
  # Nucleotide sequence
  b = 'CTGGTT' 

  global n,m
  #Length of sequence a
  n = len(a)
  #Length of Sequence b
  m = len(b)

  global matrix_E
  matrix_E = [[0 for _ in range(m + 3)] for _ in range(n + 1)]

  global i,j
  
  j_index = -2

  # Evaluate column by column
  for j in range(m+3):
    for i in range(n+1):
      if j_index >= -2 and j_index <= -1 and i == 0:
        matrix_E[i][j] = math.inf
      elif j_index >= 0 and j_index <= m and i == 0:
        matrix_E[i][j] = 0
      elif i > 0 and i <= n and j_index >= -2 and j_index <= 2:
        matrix_E[i][j] = math.inf
      else:
        E(i,j)

    for rows in matrix_E:
        print(rows)
    print()
    j_index = j_index + 1

if __name__ == '__main__':
  main()

