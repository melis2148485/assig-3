'''Given: Two DNA strings s and t (each having length at most 1 kbp) in FASTA format.
Return: A longest common subsequence of s and t. (If more than one solution exists, you may return any one.)'''
from Bio import SeqIO
def read_fasta(filename):
    sequences=[]
    for record in SeqIO.parse(filename,"fasta"):#parse fasta using biopython function(record refers to header and seq)
        sequences.append(str(record.seq))#add sequences of dna to empty list
    return sequences#return list of dna sequences

def longest_common_subseq(s,t):
    m,n=len(s), len(t)
    dp=[[0]*(n+1) for _ in range(m+1)]#initialize dp
    for i in range(1,m+1):
        for j in range(1,n+1):
            if s[i-1]==t[j-1]:
                dp[i][j]=dp[i-1][j-1]+1
            else:
                dp[i][j]=max(dp[i-1][j],dp[i][j-1])
    #recunstruct lcs from dp array
    lcs=[]
    i,j=m,n
    while i>0 and j>0:
        if s[i-1]==t[j-1]:
            lcs.append(s[i-1])
            i-=1
            j-=1
        elif dp[i-1][j]>dp[i][j-1]:
            i-=1
        else:
            j-=1
    lcs.reverse()
    return "".join(lcs)
         
def main():
    sequences=read_fasta("lcqs.fasta")
    s=sequences[0]
    t=sequences[1]
    lcs=longest_common_subseq(s,t)
    print(lcs)

if __name__=="__main__":
    main()
