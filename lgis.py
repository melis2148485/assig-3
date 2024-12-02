'''Given: A positive integer n≤10000mfollowed by a permutation π of length n.
Return: A longest increasing subsequence of π, followed by a longest decreasing subsequence of π.'''
#longest increasing sequence: longest subset in permutations including only the elements that increase
#longest decreasing seqeuncesequence:longest subset in permutations including only the elements that decrease


def longest_increasing_sequence(sequence):
    n=len(sequence)
    dp=[1]*n#initialize dp considered that each element is a subset array of length 1
    parent=[-1]*n#the earray takes all elements in reverse order, starting from last element of the
    #sequence and taking only those of the iteration
    for i in range(1,n):
        for j in range(i):
            if sequence[i]>sequence[j] and dp[i]<dp[j]+1:#extending only valid increasing sequence and 
                #update dp to include element at last position, so the element j bigger than i, so we update
                dp[i]=dp[j]+1
                parent[i]=j
            #find length of lis and its last element
    fin_len=max(dp)
    fin_ind=dp.index(fin_len)
    #lis
    lis=[]
    while fin_ind != -1:#we are checking in reverse
       #get current element and add to lis subsequence
        lis.append(sequence[fin_ind])
        fin_ind=parent[fin_ind]#restart with previous element
    lis.reverse()
    return lis
#lis=longest_increasing_sequence(sequence=)
def longest_decreasing_sequence(sequence):
    n=len(sequence)
    dp=[1]*n#initialize dp considered that each element is a subset array of length 1
    parent=[-1]*n#the earray takes all elements in reverse order, starting from last element of the
    #sequence and taking only those of the iteration
    for i in range(1,n):
        for j in range(i):
            if sequence[i]<sequence[j] and dp[i]<dp[j]+1:#extending only valid decreasing sequence and 
                #update dp to include element at last position, so the element j bigger than i, so we update
                dp[i]=dp[j]+1
                parent[i]=j
            #find length of lis and its last element
    fin_len=max(dp)
    fin_ind=dp.index(fin_len)
    #lis
    lds=[]
    while fin_ind != -1:#we are checking in reverse
#get current element and add to lis subsequence
        lds.append(sequence[fin_ind])
        fin_ind=parent[fin_ind]#restart with previous element
    lds.reverse()
    return lds
def main_funct():
    #read input from a file
    with open ("lgis.txt", "r") as file:
        lines=file.readlines()

    n=int(lines[0].strip())#first line is length of permutation
    π=lines[1].strip()#second line is the permutation
    permutation=list(map(int,π.split()))#convert space seprataed numbers into list of integers
    lis=longest_increasing_sequence(permutation)
    lds=longest_decreasing_sequence(permutation)
    print(' '.join(map(str, lis)))
    print(' '.join(map(str, lds)))
if __name__ == '__main__':
    main_funct()