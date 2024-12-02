'''Given: Two protein strings s and t in FASTA format (each of length at most 1000 aa).
Return: The maximum alignment score between s and t. Use:
The BLOSUM62 scoring matrix.
Linear gap penalty equal to 5 (i.e., a cost of -5 is assessed for each gap symbol).'''
#compute alignment score between 2 protein strings using blosum62 matrix and linear gap penalty of -5
import numpy as np
from Bio import SeqIO
blos_score_matrix= [
    [4, 0, -2, -1, -2, 0, -2, -1, -1, -1, -1, -2, -1, -1, -1, 1, 0, 0, -3, -2],
    [0, 9, -3, -4, -2, -3, -3, -1, -3, -1, -1, -3, -3, -3, -3, -1, -1, -1, -2, -2],
    [-2, -3, 6, 2, -3, -1, -1, -3, -1, -4, -3, 1, -1, 0, -2, 0, -1, -3, -4, -3],
    [-1, -4, 2, 5, -3, -2, 0, -3, 1, -3, -2, 0, -1, 2, 0, 0, -1, -2, -3, -2],
    [-2, -2, -3, -3, 6, -3, -1, 0, -3, 0, 0, -3, -4, -3, -3, -2, -2, -1, 1, 3],
    [0, -3, -1, -2, -3, 6, -2, -4, -2, -4, -3, 0, -2, -2, -2, 0, -2, -3, -2, -3],
    [-2, -3, -1, 0, -1, -2, 8, -3, -1, -3, -2, 1, -2, 0, 0, -1, -2, -3, -2, 2],
    [-1, -1, -3, -3, 0, -4, -3, 4, -3, 2, 1, -3, -3, -3, -3, -2, -1, 3, -3, -1],
    [-1, -3, -1, 1, -3, -2, -1, -3, 5, -2, -1, 0, -1, 1, 2, 0, -1, -2, -3, -2],
    [-1, -1, -4, -3, 0, -4, -3, 2, -2, 4, 2, -3, -3, -2, -2, -2, -1, 1, -2, -1],
    [-1, -1, -3, -2, 0, -3, -2, 1, -1, 2, 5, -2, -2, 0, -1, -1, -1, 1, -1, -1],
    [-2, -3, 1, 0, -3, 0, 1, -3, 0, -3, -2, 6, -2, 0, 0, 1, 0, -3, -4, -2],
    [-1, -3, -1, -1, -4, -2, -2, -3, -1, -3, -2, -2, 7, -1, -2, -1, -1, -2, -4, -3],
    [-1, -3, 0, 2, -3, -2, 0, -3, 1, -2, 0, 0, -1, 5, 1, 0, -1, -2, -2, -1],
    [-1, -3, -2, 0, -3, -2, 0, -3, 2, -2, -1, 0, -2, 1, 5, -1, -1, -3, -3, -2],
    [1, -1, 0, 0, -2, 0, -1, -2, 0, -2, -1, 1, -1, 0, -1, 4, 1, -2, -3, -2],
    [0, -1, -1, -1, -2, -2, -2, -1, -1, -1, -1, 0, -1, -1, -1, 1, 5, 0, -2, -2],
    [0, -1, -3, -2, -1, -3, -3, 3, -2, 1, 1, -3, -2, -2, -3, -2, 0, 4, -3, -1],
    [-3, -2, -4, -3, 1, -2, -2, -3, -3, -2, -1, -4, -4, -2, -3, -3, -2, -3, 11, 2],
    [-2, -2, -3, -2, 3, -3, 2, -1, -2, -1, -1, -2, -3, -1, -2, -2, -2, -1, 2, 7]
]

gap_penalty=-5
def read_fasta(filename):
    sequences=[]
    for record in SeqIO.parse(filename,"fasta"):#parse fasta using biopython function(record refers to header and seq)
        sequences.append(str(record.seq))#add sequences of dna to empty list
    return sequences#return list of sequences

#map aminoacids to their indices in the score matrix
amino_a=['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
amino_a_index={aa:idx for idx, aa in enumerate(amino_a)}

def get_score(a,b):#returns score for aligning characters 'a'and 'b'based on blosum matrix
    return blos_score_matrix[amino_a_index[a]][amino_a_index[b]]
    
def max_alignment_score(s,t):#compute max alignment between s and t, starting from blosum and including gap penalty
    m,n=len(s), len(t)
    dp=[[0]*(n+1) for _ in range(m+1)]
    #to initialize dp
    for i in range(1, m+1):
        dp[i][0]=dp[i-1][0]+ gap_penalty
    for j in range(1, n+1):
        dp[0][j]=dp[0][j-1] + gap_penalty
        #fill dp
    for i in range(1,m+1):#loop through s
        for j in range(1,n+1):#loop thorugh t
            match=dp[i-1][j-1]+ get_score(s[i-1],t[j-1])#match/mismatch
            delete=dp[i-1][j]+gap_penalty#deletion
            insert=dp[i][j-1]+gap_penalty #insertion
            dp[i][j]=max(match,delete,insert)#this  time we take the maximum score
    return dp[m][n]

def main():
    sequences=read_fasta("glob.fasta")
    s,t=sequences[0], sequences[1]
    score=max_alignment_score(s,t)
    print(score)
    
if __name__=="__main__":
    main()
