'''Given: Two protein strings s and t in FASTA format, each of length at most 1000 aa.
Return: The total number of optimal alignments of s and t with respect to edit alignment score, modulo 134,217,727 (227-1).'''
from Bio import SeqIO
MOD=134217727
def read_fasta(filename):
    sequences=[]
    for record in SeqIO.parse(filename,"fasta"):#parse fasta using biopython function(record refers to header and seq)
        sequences.append(str(record.seq))#add sequences of dna to empty list
    return sequences#return list of sequences

def edit_dist_count(s,t):#using dp array we first calculate edit distance and count to number the optimal alignments
    m,n=len(s), len(t)
    dp=[[0]*(n+1) for _ in range(m+1)]
    count=[[0]*(n+1) for _ in range(m+1)]
    #to initialize dp
    for i in range(m+1):
        dp[i][0]=i
        count[i][0]=1
    for j in range(n+1):
        dp[0][j]=j
        count[j][0]=1
    for i in range(1,m+1):#loop through s
        for j in range(1,n+1):#loop thorugh t
            if s[i-1]==t[j-1]:#the 2 characters match, we don't have to do nothing 
                dp[i][j]=dp[i-1][j-1]
                count[i][j]=count[i-1][j-1]
            else:#if the characters are different, is because of 3 possible operations, so deletion, inserion and substituition:for any of thse we add 1
                dp[i][j]=min(dp[i-1][j]+1,dp[i][j-1]+1,dp[i-1][j-1]+1)#backtrack to construct alignment
            if dp[i][j]==dp[i-1][j]+1:
                count [i][j]+=count[i-1][j]
            if dp[i][j]==dp[i][j-1]+1:
                count [i][j]+=count[i][j-1]
            if dp[i][j]==dp[i-1][j-1]+1:
                count [i][j]+=count[i-1][j-1]
            count[i][j]%=MOD
    return dp[m][n], count[m][n]
def main():
    sequences=read_fasta("ctea.fasta")
    s,t=sequences[0],sequences[1]
    edit_dist, aligned_count=edit_dist_count(s,t)
    print(aligned_count)
   
if __name__=="__main__":
    main()


