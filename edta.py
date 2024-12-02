'''Given: Two protein strings s and t in FASTA format (with each string having length at most 1000 aa).
Return: The edit distance dE(s,t) followed by two augmented strings s′ and t′ representing an optimal 
alignment of s and t.'''
from Bio import SeqIO
import numpy as np
def read_fasta(filename):
    sequences=[]
    for record in SeqIO.parse(filename,"fasta"):#parse fasta using biopython function(record refers to header and seq)
        sequences.append(str(record.seq))#add sequences of dna to empty list
    return sequences#return list of sequences
def edit_distance_alignment(s,t):#using dp array we first calculate edit distance
    m,n=len(s), len(t)
    dp=[[0]*(n+1) for _ in range(m+1)]
    #to initialize dp
    for i in range(m+1):
        dp[i][0]=i
    for j in range(n+1):
        dp[0][j]=j
        #fill dp
    for i in range(1,m+1):#loop through s
        for j in range(1,n+1):#loop thorugh t
            if s[i-1]==t[j-1]:#the 2 characters match, we don't have to do nothing 
                dp[i][j]=dp[i-1][j-1]
            else:#if the characters are different, is because of 3 possible operations, so deletion, inserion and substituition:for any of thse we add 1
                dp[i][j]=min(dp[i-1][j]+1,dp[i][j-1]+1,dp[i-1][j-1]+1)#backtrack to construct alignment
    i,j=m,n#i current position of s during backtracking and j current position of t during backtracking
    aligned_s=[]
    aligned_t=[]
    while i>0 or j>0:
        if i>0 and j>0 and dp[i][j]==dp[i-1][j-1]+(0 if s[i-1]==t[j-1] else 1):
            aligned_s.append(s[i-1])
            aligned_t.append(t[j-1])
            i-=1
            j-=1
        elif i>0 and dp[i][j]==dp[i-1][j]+1:
            aligned_s.append(s[i-1])
            aligned_t.append('-')
            i-=1
        else:
            aligned_s.append('-')
            aligned_t.append(t[j-1])
            j-=1
    aligned_s=''.join(reversed(aligned_s))#we backtrack in dynamic programming so both are reversed
    aligned_t=''.join(reversed(aligned_t))
    return dp[m][n], aligned_s, aligned_t

def main():
    s,t=read_fasta("edta.fasta")
    edit_dist, aligned_s, aligned_t=edit_distance_alignment(s,t)
    print(edit_dist)
    print(aligned_s)
    print(aligned_t)
    
if __name__=="__main__":
    main()

