import math

def readObeservationFile():
    arr = []
    file = open("observation.csv", "r")
    for line in file:
        arr.append(line.split(','))
    file.close()
    B_N={}
    B_G={}
    i=0
    for j in ["A","C","G","T"]:
        for k in ["A", "C", "G", "T"]:
            for l in ["A", "C", "G", "T"]:
                    val1 = float(arr[0][i])
                    val2 = float(arr[1][i])
                    B_N[j+k+l]=val1
                    B_G[j+k+l]=val2
                    i=i+1
    B = [B_N, B_G]
    return B

def readFile(filename):
    file = open(filename, "r")
    string = ""
    for line in file:
        string = string + line.strip("\n")
    array = []
    pos=0
    for i in range(0,len(string)):
        if i!=0 and i%3==2:
            array.append(string[pos:i+1])
            pos = i+1
    file.close()
    return array

def Viterbi(sequence, PI, B, A):
    T = len(sequence)

    #initialization
    type = sequence[0]
    delta=[]
    fai = []
    for l in range(0,len(PI)):

        delta.append([PI[l]*B[l][type]])
        #if PI[l]==0:
        #    delta.append([0.0])
        #else:
        #    delta.append([-math.log(PI[l])-math.log(B[l][type])])
        fai.append([0])

    #iteration for t
    for t in range(1,T):
        for i in range(0,len(PI)):
            tmp_max1 = float('-inf')
            pos1 = 0
            for j in range(0, len(PI)):
                tmp1 = delta[j][t-1]*A[j][i]
                #tmp1 = delta[j][t - 1]-math.log(A[j][i])
                if tmp1>tmp_max1:
                    tmp_max1 = tmp1
                    pos1 = j
            type_t = sequence[t]
            delta[i].append(tmp_max1 * B[i][type_t])
            #if B[i][type_t]<= 0:
            #   print(B[i][type_t])
            #delta[i].append(tmp_max1-math.log(B[i][type_t]))
            fai[i].append(pos1)

    #finalization
    i_star = []
    tmp_max2 = 0.0
    pos2 = 0
    for m in range(0,len(PI)):
        tmp2 = delta[m][T-1]
        if tmp2>tmp_max2:
            tmp_max2 = tmp2
            pos2 = m
    P_star = tmp_max2
    #P_star = math.exp(-tmp_max2)
    i_star.insert(0,pos2)

    #backward
    for t in range(0,T-1):
        i_star.insert(0,fai[i_star[0]][T-t-1])

    I = i_star

    str = ""
    for n in range(0,len(I)):
        if I[n]==0:
            str = str+"N"
        elif I[n]==1:
            str = str+"G"

    #calculate gene's number
    count = 0
    coordinates = []
    u=0
    while u < len(str):
        end = 1
        if str[u] == 'G':
            count = count+1
            for v in range(u+1,len(str)):
                if str[v] == 'G' and v == len(str)-1:
                    coordinates.append([u, v])
                    break
                elif str[v] == 'N':
                    coordinates.append([u,v-1])
                    end = v+1-u
                    break
        u = u+end


    return [str, P_star, count, coordinates]

def main():
    #PI contains the initial probabilities of states
    pi_N = 1
    pi_G = 0
    PI = [pi_N, pi_G]

    #A contains states transition probabilities
    print("Please input the transition probabilities p(N,N) and p(G,G):")
    p_NN = 0.93 #float(input("P(N,N): "))
    p_GG = 0.93 #float(input("P(G,G): "))
    p_NG = 1 - p_NN
    p_GN = 1 - p_GG
    A = [[p_NN,p_NG],[p_GN, p_GG]]

    #B contains observation probabilities
    B = readObeservationFile()

    #Calculate the most possible states sequence and joint probability
    print("--------------First document---------------")
    sequence1 = readFile("1_39675.txt")
    print("Observation sequence of 1_39675.txt is:")
    print(sequence1)
    res1 = Viterbi(sequence1,PI,B,A)
    print("The most possible states sequence is:")
    print(res1[0])
    print("The joint probability of path and sequence is:")
    print(res1[1])
    print("The number of the predicted genes in the sequence provided is:")
    print(res1[2])
    print("The start-end coordinates of the predicted genes in the sequence provided are:")
    print(res1[3])
    print()

    print("--------------Second document---------------")
    sequence2 = readFile("1000_2000_1.txt")
    print("Observation sequence of 1000_2000_1.txt is:")
    print(sequence2)
    res2 = Viterbi(sequence2,PI,B,A)
    print("The most possible states sequence is:")
    print(res2[0])
    print("The joint probability of path and sequence is:")
    print(res2[1])
    print("The number of the predicted genes in the sequence provided is:")
    print(res2[2])
    print("The start-end coordinates of the predicted genes in the sequence provided are:")
    print(res2[3])

    print()

    print("--------------Third document---------------")
    sequence3 = readFile("1000_5000.txt")
    print("Observation sequence of 1000_5000.txt is:")
    print(sequence3)
    res3 = Viterbi(sequence3,PI,B,A)
    print("The most possible states sequence is:")
    print(res3[0])
    print("The joint probability of path and sequence is:")
    print(res3[1])
    print("The number of the predicted genes in the sequence provided is:")
    print(res3[2])
    print("The start-end coordinates of the predicted genes in the sequence provided are:")
    print(res3[3])

main()