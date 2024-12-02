'''Given: Two DNA strings s and t (each of length at most 1 kbp) in FASTA format.
Return: One collection of indices of s in which the symbols of t appear as a subsequence of s. 
If multiple solutions exist, you may return any one.'''
#biopython to read fasta files
from Bio import SeqIO
def read_fasta(filename):
    sequences=[]
    for record in SeqIO.parse(filename,"fasta"):#parse fasta
        sequences.append(str(record.seq))
    return sequences

def subseq_indices(s,t):
    indices=[]#store indices
    s_index=0
    for char in t:#iterate in t to find matching subsequences
        while s_index<len(s)and s[s_index] != char:#find occurrence in s-index
            s_index+=1#if char found, record the index
        if s_index<len(s):
            indices.append(s_index+1)#comvert to 1-based index
            s_index +=1#move in s
        else:
            break#exit if char not found
    return indices
#main function that reads input and solves problem
def main():
    sequences=read_fasta("sseq.fasta")
    s=sequences[0]#first sequence
    t=sequences[1]#second sequence
    indices=subseq_indices(s,t)#indices where t is a subsequence of s
    print(" ".join(map(str, indices)))

#run
if __name__=="__main__":
    main()