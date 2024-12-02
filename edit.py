'''Given: Two protein strings s and t in FASTA format (each of length at most 1000 aa).
Return: The edit distance dE(s,t).'''
from Bio import SeqIO
def read_fasta(filename):
    sequences=[]
    for record in SeqIO.parse(filename,"fasta"):#parse fasta using biopython function(record refers to header and seq)
        sequences.append(str(record.seq))#add sequences of dna to empty list
    return sequences#return list of dna sequences
def edit_dist(s,t):
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
            else:#if the chacaters are fifferent, is because of 3 possible operations, so deletion, inserion and substituition:for any of thse we add 1
                dp[i][j]=min(dp[i-1][j]+1,dp[i][j-1]+1,dp[i-1][j-1]+1)
    return dp[m][n]
def main():
    sequences=read_fasta("edit.fasta")
    s=sequences[0]
    t=sequences[1]
    distance=edit_dist(s,t)
    print(distance)
    
if __name__=="__main__":
    main()