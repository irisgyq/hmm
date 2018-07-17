def readFile():
    sequence = []
    file = open("diceSequences.csv", "r")
    for line in file:
        sequence.append(line.strip("\n").split(','))
    file.close()
    return sequence

def Viterbi(sequence, PI, B, A):
    T = len(sequence)

    # initialization
    type = sequence[0]
    delta = []
    fai = []
    for l in range(0, len(PI)):
        delta.append([PI[l] * B[l][type]])
        fai.append([0])

    # iteration for t
    for t in range(1, T):
        for i in range(0, len(PI)):
            tmp_max1 = -1.0
            pos1 = 0
            for j in range(0, len(PI)):
                tmp1 = delta[j][t - 1] * A[j][i]
                if tmp1 >= tmp_max1:
                    tmp_max1 = tmp1
                    pos1 = j
            type_t = sequence[t]
            delta[i].append(tmp_max1 * B[i][type_t])
            fai[i].append(pos1)

    # finalization
    i_star = []
    tmp_max2 = -1.0
    pos2 = 0
    for m in range(0, len(PI)):
        tmp2 = delta[m][T - 1]
        if tmp2 >= tmp_max2:
            tmp_max2 = tmp2
            pos2 = m
    P_star = tmp_max2
    i_star.insert(0, pos2)

    # backward
    for t in range(0, T - 1):
        i_star.insert(0, fai[i_star[0]][T - t - 1])

    I = i_star
    string = ""
    for n in range(0,len(I)):
        string = string+"D"+str(I[n]+1)+" "
    return [string,P_star]

def main():
    pi_1 = 1/3
    pi_2 = 1/3
    pi_3 = 1/3
    PI = [pi_1,pi_2,pi_3]
    sequence = readFile()
    A = [[0.5, 0.125, 0.125], [0.125, 0.5, 0.125], [0.125, 0.125, 0.5]]
    B = [{"1": 0.6, "2": 0.2, "3": 0.2}, {"1": 0.2, "2": 0.6, "3": 0.2}, {"1": 0.2, "2": 0.2, "3":0.6}]
    for i in range(0,len(sequence)):
        res = Viterbi(sequence[i], PI, B, A)
        print("-------The number "+str(i+1)+" sequence----------")
        print("The possible states sequence is:")
        print(res[0])
        print("The probability of this sequence is:")
        print(res[1])
        print()

main()