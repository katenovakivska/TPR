import numpy as np
import math
condition = [[10, 3, 3, 10, 4, 1, 6, 10, 5, 10, 3, 1],
             [2, 2, 7, 3, 8, 7, 9, 7, 1, 4, 8, 3],
             [2, 5, 1, 2, 2, 1, 9, 5, 2, 6, 4, 4],
             [9, 10, 1, 4, 7, 6, 3, 9, 7, 3, 4, 2],
             [1, 2, 1, 1, 5, 4, 2, 6, 4, 1, 10, 6],
             [4, 5, 4, 10, 7, 2, 5, 2, 5, 3, 2, 10],
             [10, 9, 1, 7, 7, 5, 6, 3, 6, 3, 10, 4],
             [5, 3, 2, 6, 5, 9, 3, 8, 7, 8, 5, 9],
             [8, 7, 7, 5, 1, 7, 2, 2, 6, 7, 6, 6],
             [5, 3, 7, 3, 1, 1, 8, 4, 3, 1, 4, 3],
             [5, 6, 2, 4, 9, 6, 9, 1, 10, 6, 8, 2],
             [7, 6, 10, 6, 7, 6, 6, 5, 9, 7, 10, 2],
             [3, 5, 7, 10, 6, 1, 10, 5, 7, 10, 7, 6],
             [8, 8, 4, 9, 6, 1, 4, 4, 8, 3, 9, 3],
             [6, 5, 5, 1, 1, 6, 5, 6, 7, 9, 6, 6]]
weight = [1, 10, 7, 7, 7, 8, 5, 10, 1, 7, 2, 1]
optimization2 = ["max","max","max","max","max","max","max","min","min","min","min","min"]
optimization1 = ["max","max","max","max","max","max","max","max","max","max","max","max"]
# condition = [[5, 8, 4],
#              [7, 6, 8],
#              [8, 8, 6],
#              [7, 4, 6]]
#
# weight = [0.3, 0.4, 0.3]
# optimization = ["max","max","max"]
#optimization2 = ["max","max","min"]
v = 0.5
def check_optimization(optimization):
    isAllMax = True
    for i in range(0, len(optimization)):
        if optimization[i] != "max":
            isAllMax = False
            break

    return isAllMax

def normalization_of_estimations(optimization):
    normalized = condition
    transponate = np.array(condition).transpose()
    sum = []
    desired_values = []
    worst_values = []
    if check_optimization(optimization) == True:
        for i in range(0, len(transponate)):
            sum_of_column = 0.00
            for j in range(0, len(transponate[i])):
                sum_of_column += math.pow(condition[j][i], 2)
            sum.append(sum_of_column)
        for i in range(0, len(normalized)):
            for j in range(0, len(normalized[i])):
                normalized[i][j] = normalized[i][j] / math.sqrt(sum[j])

    elif check_optimization(optimization) == False:
        for i in range(0, len(transponate)):
            if optimization[i] == "max":
                worst = np.array(transponate[i]).min()
                best = np.array(transponate[i]).max()
            elif optimization[i] == "min":
                worst = np.array(transponate[i]).max()
                best = np.array(transponate[i]).min()
            desired_values.append(best)
            worst_values.append(worst)

        for i in range(0, len(normalized)):
            for j in range(0, len(normalized[i])):
                if optimization[j] == "max":
                    normalized[i][j] = (condition[i][j] - worst_values[j]) / (desired_values[j] - worst_values[j])
                elif optimization[j] == "min":
                    normalized[i][j] = (worst_values[j] - condition[i][j]) / (worst_values[j] - desired_values[j])

    return normalized

def weighted_estimations(optimization):
    estimations = normalization_of_estimations(optimization)
    weight_sum = 0
    for i in range(0, len(weight)):
        weight_sum += weight[i]
    for i in range(0, len(weight)):
        weight[i] /= weight_sum
    for i in range(0, len(estimations)):
        for j in range(0, len(estimations[i])):
            estimations[i][j] *= weight[j]

    # for i in range(0, len(estimations)):
    #     print(estimations[i])

    return estimations

def find_nis(estimations):
    weighted_estimations = np.array(estimations).transpose()
    nis = []
    for i in range(0, len(weighted_estimations)):
        min_value = np.array(weighted_estimations[i]).min()
        nis.append(min_value)

    return nis

def find_pis(estimations):
    weighted_estimations = np.array(estimations).transpose()
    pis = []
    for i in range(0, len(weighted_estimations)):
        max_value = np.array(weighted_estimations[i]).max()
        pis.append(max_value)
    return pis

def distance_to_pis(PIS, weighted_estimations):
    d_better = []
    for i in range(0, len(weighted_estimations)):
        difference = 0.00
        for j in range(0, len(weighted_estimations[i])):
            difference += math.pow(weighted_estimations[i][j] - PIS[j],2)
        d_better.append(math.sqrt(difference))
    return d_better

def distance_to_nis(NIS, weighted_estimations):
    d_worse = []
    for i in range(0, len(weighted_estimations)):
        difference = 0.00
        for j in range(0, len(weighted_estimations[i])):
            difference += math.pow(weighted_estimations[i][j] - NIS[j], 2)
        d_worse.append(math.sqrt(difference))

    return  d_worse

def proximity(d_pis, d_nis):
    ck = []
    for i in range(0, len(d_pis)):
        ck.append(d_nis[i] / (d_pis[i] + d_nis[i]))

    return  ck

def topsis(optimization):
    estimations = weighted_estimations(optimization)
    PIS = find_pis(estimations)
    NIS = find_nis(estimations)
    d_pis = distance_to_pis(PIS, estimations)
    d_nis = distance_to_nis(NIS, estimations)
    ck = proximity(d_pis, d_nis)
    alternatives = np.argsort(np.array(ck))
    alternatives = [x+1 for x in alternatives]
    print("Найкраща альтернатива: ", (alternatives[len(alternatives)-1]))
    print("Ранжування: ", alternatives[::-1])

def find_min():
    transpose = np.array(condition).transpose()
    min_values = []
    for i in range(0, len(transpose)):
        min_value = np.array(transpose[i]).min()
        min_values.append(min_value)

    return min_values

def find_max():
    transpose = np.array(condition).transpose()
    max_values= []
    for i in range(0, len(transpose)):
        max_value = np.array(transpose[i]).max()
        max_values.append(max_value)

    return max_values

def find_sk(best, worse):
    sk = condition

    if np.array(weight).sum() > 1:
        for i in range(0, len(weight)):
            weight[i] /= np.array(weight).sum()

    for i in range(0, len(condition)):
        for j in range(0, len(condition[i])):
            sk[i][j] = weight[j] * math.fabs(best[j] - condition[i][j]) / math.fabs(best[j] - worse[j])

    # for i in range(0, len(sk)):
    #     print(sk[i])
    return sk

def find_sk_sum(sk):
    sum_sk = []
    for i in range(0, len(sk)):
        sum_sk.append(np.array(sk[i]).sum())

    return  sum_sk

def find_rk(sk):
    rk = []
    for i in range(0, len(sk)):
        max_value = np.array(sk[i]).max()
        rk.append(max_value)

    return rk

def find_Q(sk_sum, rk):
    min_sk = np.array(sk_sum).min()
    max_sk = np.array(sk_sum).max()
    min_rk = np.array(rk).min()
    max_rk = np.array(rk).max()

    Qk = []
    for i in range(0, len(sk_sum)):
        Qk.append((v * (sk_sum[i] - min_sk) / (max_sk - min_sk)) + ((1 - v) * (rk[i] - min_rk) / (max_rk - min_rk)))

    return Qk

def c1_condition(Q, Q_sorted):
    c1 = []
    Q_s = Q_sorted[::-1]
    c1.append(Q_s[0])

    for i in range(0, len(Q_sorted)):
        if i != 0:
            if (Q[Q_s[i]] - Q[Q_s[0]]) < (1 / (len(Q) - 1)) :
                c1.append(Q_sorted[i])
            elif (Q[Q_s[i]] - Q[Q_s[0]]) >= (1 / (len(Q) - 1)) :
                break
    if len(c1) < 2:
        print("Умова С1 не виконується")
    elif len(c1) >= 2:
        print("Умова С1 виконується")
    return c1

def c2_condition(Q_sorted, S_sorted, R_sorted, c1):
    c2 = []
    indexes = []

    for i in range(0, len(c1)):
        indexes.append(np.where(Q_sorted == c1[i]))

    for i in range(0, len(indexes)):
        if Q_sorted[indexes[i]] == S_sorted[indexes[i]] or Q_sorted[indexes[i]] == R_sorted[indexes[i]]:
            c2.append(c1[i])

    if len(c2) < 1:
        print("Умова С2 не виконується")
    elif len(c2) >= 1:
        print("Умова С2 виконується")
    return  c2

def vikor(v):
    best = find_max()
    worse = find_min()
    sk = find_sk(best, worse)
    sk_sum = find_sk_sum(sk)
    rk = find_rk(sk)
    Q = find_Q(sk_sum, rk)

    Q_sorted = np.argsort(np.array(Q))
    Q_sorted = Q_sorted[::-1]
    S_sorted = np.argsort(np.array(sk_sum))
    S_sorted = S_sorted[::-1]
    R_sorted = np.argsort(np.array(rk))
    R_sorted = R_sorted[::-1]

    c1 = c1_condition(Q, Q_sorted)
    c2 = c2_condition(Q_sorted, S_sorted, R_sorted, c1)
    c2 = [x + 1 for x in c2]
    Q_sorted = Q_sorted[::-1]
    S_sorted = S_sorted[::-1]
    R_sorted = R_sorted[::-1]
    print("Ранжування Q, S, R: ")
    Q_sorted = [x + 1 for x in Q_sorted]
    S_sorted = [x + 1 for x in S_sorted]
    R_sorted = [x + 1 for x in R_sorted]
    print(Q_sorted)
    print(S_sorted)
    print(R_sorted)
    if len(c2) > 1:
        print("Компромісні розв'язки: ", c2)
    elif len(c2) == 1:
        print("Найкраща альтернатива: ", c2[len(c2) - 1])
        print("Q значення найкращої альтернативи", np.array(Q).min())
    print("Ранжування: ", Q_sorted)


print("Метод TOPSIS - всі критерії максимізуються")
topsis(optimization1)
print("--------------------------------------------------------------------")
print("Метод TOPSIS - критерії k1-k7 максимізуються, а k8-k12 мінімізуються")
topsis(optimization2)
print("--------------------------------------------------------------------")
print("Метод VIKOR - всі критерії максимізуються")
vikor(v)
print("--------------------------------------------------------------------")
print("ДОСЛІДЖЕННЯ")
print("--------------------------------------------------------------------")

for i in range(0, 11):
    V = round(0.1 * i, 1)
    print("Дослідження при V = ", V, " :")
    vikor(V)
    print("--------------------------------------------------------------------")