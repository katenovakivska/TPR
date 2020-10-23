import numpy as np
condition = [[4, 6, 6, 2, 9, 10, 3, 10, 6, 3, 9, 3],
             [ 4, 9, 6, 5, 9, 10, 10, 10, 8, 4, 9, 3],
             [ 5, 9, 6, 6, 9, 10, 10, 10, 8, 9, 9, 3],
             [ 5, 2, 6, 4, 9, 5, 5, 2, 5, 6, 9, 3],
             [ 5, 6, 8, 5, 10, 10, 7, 10, 6, 6, 9, 6],
             [ 5, 7, 8, 9, 10, 10, 7, 10, 7, 10, 9, 9],
             [ 6, 7, 8, 9, 10, 10, 10, 10, 7, 10, 9, 9],
             [ 6, 7, 8, 9, 10, 10, 10, 10, 7, 10, 9, 9],
             [ 4, 2, 5, 3, 7, 3, 4, 4, 7, 4, 6, 1],
             [ 4, 7, 5, 6, 8, 7, 7, 4, 7, 4, 6, 2],
             [ 4, 4, 5, 6, 3, 3, 4, 4, 5, 2, 6, 2],
             [ 6, 10, 10, 7, 3, 8, 4, 4, 9, 2, 6, 3],
             [ 4, 2, 6, 3, 3, 5, 4, 2, 2, 2, 3, 3],
             [ 8, 6, 9, 3, 9, 5, 8, 8, 9, 2, 9, 9],
             [ 8, 10, 10, 7, 10, 8, 8, 8, 9, 3, 10, 9],
             [ 9, 10, 10, 7, 10, 8, 8, 8, 9, 10, 10, 9],
             [ 5, 4, 6, 7, 1, 8, 3, 1, 3, 8, 9, 2],
             [ 5, 4, 6, 1, 1, 6, 3, 1, 3, 8, 2, 2],
             [ 4, 2, 5, 1, 1, 5, 3, 1, 2, 2, 2, 2],
             [ 1, 2, 5, 1, 1, 4, 3, 1, 2, 2, 1, 2]]

# пошук відношення Подиновського
def podynovskiy(matrix):
    # сортування за спаданням
    fi_a = fi_alternatives(matrix)
    # формування сигма побудованих для фі(аі)
    sv = sigma_matrix(fi_a)
    # відношення Парето для множини фі(аі)
    result = pareto(sv)

    f = open("task5.txt", "w")
    f.write(" 5" + '\n')
    for i in range(0, len(result)):
        for j in range(0, len(result)):
            f.write(" {} ".format(result[i][j]))
        f.write('\n')
    f.close()
    return result

# пошук відношення Березовського
def berezovskiy():
    # Pb1, Ib1, Nb1
    pb1 = [[0] * 20 for i in range(20)]
    ib1 = [[0] * 20 for i in range(20)]
    nb1 = [[0] * 20 for i in range(20)]
    result = [[0] * 20 for i in range(20)]
    # розбиваємо критерії на класи
    class1, class2, class3 = divide_on_classes(SV)
    pareto1 = np.array(pareto(class1))
    pareto2 = np.array(pareto(class2))
    pareto3 = np.array(pareto(class3))
    # виділення симетричної, асиметричної та непорівнюваної частини
    # частини першої ітерації
    i1 = find_I(pareto1)
    p1 = find_P(pareto1)
    n1 = find_N(pareto1)
    # частини другої ітерації
    i2 = find_I(pareto2)
    p2 = find_P(pareto2)
    # частини третьої ітерації
    i3 = find_I(pareto3)
    p3 = find_P(pareto3)

    for i in range(0, len(result)):
        for j in range(0, len(result)):
            if p2[i][j] == 1 and p2[i][j] == p1[i][j]:
                pb1[i][j] = 1
            if p2[i][j] == 1 and p2[i][j] == n1[i][j]:
                pb1[i][j] = 1
            if i2[i][j] == 1 and i2[i][j] == p1[i][j]:
                pb1[i][j] = 1
            if p2[i][j] == 1 and p2[i][j] == i1[i][j]:
                pb1[i][j] = 1
            if i2[i][j] == 1 and i2[i][j] == i1[i][j]:
                ib1[i][j] = 1
    for i in range(0, len(result)):
        for j in range(0, len(result)):
            if pb1[i][j] == 0 and ib1[i][j] == 0:
                nb1[i][j] == 1
    # порівняння частин
    for i in range(0, len(result)):
        for j in range(0, len(result)):
            if p3[i][j] == 1 and p3[i][j] == pb1[i][j]:
                result[i][j] = 1
            if p3[i][j] == 1 and p3[i][j] == nb1[i][j]:
                result[i][j] = 1
            if i3[i][j] == 1 and i3[i][j] == pb1[i][j]:
                result[i][j] = 1
            if p3[i][j] == 1 and p3[i][j] == ib1[i][j]:
                result[i][j] = 1
    f = open("task4.txt", "w")
    f.write(" 4" + '\n')
    for i in range(0, len(result)):
        for j in range(0, len(result)):
            f.write(" {} ".format(result[i][j]))
        f.write('\n')
    f.close()
    return result

#пошук відношення парето
def pareto(sigma_matrix):
    matrix = np.array(sigma_matrix)
    result = [[0]*20 for i in range(20)]
    for i in range(0, 20):
        result[i][i] = 1
    for i in range(0, len(matrix)):
        for j in range(i+1, len(matrix)):
            flag = pareto_check(matrix[i][j])
            if flag == 1:
                result[i][j] = 1
                result[j][i] = 0
                continue
            elif flag == 2:
                result[i][j] = 0
                result[j][i] = 1
                continue
            elif flag == 3:
                result[i][j] = 0
                result[j][i] = 0
    f = open("task1.txt", "w")
    f.write(" 1" + '\n')
    for i in range(0, len(result)):
        for j in range(0, len(result)):
            f.write(" {} ".format(result[i][j]))
        f.write('\n')
    f.close()
    return result

#пошук відношення мажоритарності
def majority(sigma_matrix):
    matrix = np.array(sigma_matrix)
    result = [[0]*20 for i in range(20)]
    for i in range(0,len(matrix)):
        for j in range (i,len(matrix)):
            if matrix[i][j].sum()>0:
                result[i][j]=1
                result[j][i]=0
                continue
            elif matrix[i][j].sum()<0:
                result[i][j]=0
                result[j][i]=1
                continue
            elif matrix[i][j].sum()==0:
                result[i][j]=0
                result[j][i]=0
    f = open("task2.txt", "w")
    f.write(" 2" + '\n')
    for i in range(0, len(result)):
        for j in range(0, len(result)):
            f.write(" {} ".format(result[i][j]))
        f.write('\n')
    f.close()
    return result

# пошук лексикографічного відношення
def lexicographic(matrix):
    sorted_matrix = np.array(matrix)
    result = [[0]*20 for i in range(20)]
    for i in range (0, len(result)):
        result[i][i]=0
    for i in range(0, len(sorted_matrix)):
        for j in range(i+1, len(sorted_matrix)):
            alternative = sorted_matrix[i][j]
            for k in range(0, len(alternative)):
                if alternative[k]==1:
                    result[i][j]=1
                    result[j][i]=0
                    break
                elif alternative[k]==-1:
                    result[i][j]=0
                    result[j][i]=1
                    break
                else:
                    continue
    f = open("task3.txt", "w")
    f.write(" 3" + '\n')
    for i in range(0, len(result)):
        for j in range(0, len(result)):
            f.write(" {} ".format(result[i][j]))
        f.write('\n')
    f.close()
    return result

# найбільші по Р за домінуванням
def domination_max_p(matrix):
    matrix = np.array(matrix)
    max_alternatives = []
    for i in range(0, len(matrix)):
        if matrix[i][i] == 0 and matrix[i].sum() == len(matrix) - 1:
            max_alternatives.append(i)
    return max_alternatives

# найбільші по R за домінуванням
def domination_max_r(matrix):
    matrix = np.array(matrix)
    max_alternatives = []
    strong_max = []
    for i in range(0, len(matrix)):
        if matrix[i].sum() == len(matrix):
            max_alternatives.append(i)
            if matrix[:, i].sum() == 1:
                strong_max.append(i)
    return max_alternatives, strong_max

# максимальні по Р за блокуванням
def block_max_p(matrix):
    matrix = np.array(matrix)
    max_alternatives = []
    for i in range(0, len(matrix)):
        if matrix[:, i].sum() == 0:
            max_alternatives.append(i)
    return max_alternatives

# максимальні по R за блокуванням
def block_max_r(matrix):
    matrix = np.array(matrix)
    symmetric = find_I(matrix)
    max_alternatives = []
    strong_max = []
    for i in range(0, len(matrix)):
        if np.any(np.array_equal(matrix[:, i], symmetric[:, i]) == False) == False:
            max_alternatives.append(i)
            if matrix[:, i].sum() == 1 and matrix[i][i] == 1:
                strong_max.append(i)
    return max_alternatives, strong_max

#порівняння двох альтернатив
def compare_alternatives (alternative1, alternative2):
    comparison_result = np.arange(len(alternative1))
    for i in range(0, len(alternative1)):
        if alternative1[i] > alternative2[i]:
            comparison_result[i] = 1
            continue
        elif alternative1[i] < alternative2[i]:
            comparison_result[i] = -1
            continue
        elif alternative1[i] == alternative2[i]:
            comparison_result[i] = 0
    return comparison_result

#пошук сигма, попарні порівняння альтернатив
def sigma_matrix(matrix):
    alternatives = np.array(matrix)
    result = [[0]*20 for i in range(20)]
    for i in range(0, 20):
        for j in range(0, 20):
            result[i][j] = compare_alternatives(alternatives[i], alternatives[j])
    return result

SV = sigma_matrix(condition)

# виділення симетричної частини
def find_I(r):
    return (r == r.T) * r

# виділення асиметричної частини
def find_P(r):
    return r - find_I(r)

# виділення непорівнюваної частини
def find_N(r):
    return (r == r.T) - find_I(r)

# перевірка частини на симетричність
def check_is_symetric(matrix):
    symmetric = find_I(np.array(matrix))
    is_sum_positive = False
    for i in range(0, len(symmetric)):
        if symmetric[i].sum() > 0:
            is_sum_positive = True
            break
    return is_sum_positive

#вивід максимальних та найбільших елементів
def print_optimal_alternative(matrix):
    if check_is_symetric(matrix):
        m, sm = domination_max_r(matrix)
        if len(m)>0:
            print("Найбільші по R: {}".format(m))
            print("Строго найбільші по R: {}".format(sm))
        else:
            m, sm = block_max_r(matrix)
            print("Максимальні по R: {}".format(m))
            print("Строго максимальні по R: {}".format(sm))
    else:
        if len(domination_max_p(matrix))>0:
            print("Найбільший по Р: {}".format(domination_max_p(matrix)))
        else:
            print("Максимальні по Р: {}".format(block_max_p(matrix)))

#перевірка елементів сигма-масиву
def pareto_check(array):
    check = np.array(array)
    all_elements_not_negative = True
    for i in range (0, len(check)):
        if check[i]>=0:
            pass
        else:
            all_elements_not_negative = False
    if all_elements_not_negative:
        return 1
    else:
        all_elements_not_negative = True
        for i in range (0, len(check)):
            if check[i]<=0:
                pass
            else:
                all_elements_not_negative = False
        if all_elements_not_negative:
            return 2
        else:
            return 3

#k6>k9>k1>k8>k7>k4>k2>k11>k10>k12>k5>k3

# строге впорядкування критеріїв
def sort_strong(alternative):
    sorted = np.arange(len(alternative))
    order = np.array([ 5, 8, 0, 7, 6, 3, 1, 10, 9, 11, 4, 2])
    for i in range (0, len(alternative)):
        sorted[i]=alternative[order[i]]
    return sorted

# упорядкування критеріїв
def sort_sigma(sigma_matrix):
    matrix = np.array(sigma_matrix)
    sorted_matrix = [[0]*20 for i in range(20)]
    for i in range(0,len(matrix)):
        for j in range(0,len(matrix)):
            sorted_matrix[i][j] = sort_strong(matrix[i][j])
    return sorted_matrix

# {k1,k8,k10,k11} < {k2,k12} < {k3,k4,k5,k6,k7,k9}

#виділення класів
def select_classes(alternative):
    select1 = np.array([alternative[0], 0, 0, 0, 0, 0, 0, alternative[7], 0, alternative[9], alternative[10], 0])
    select2 = np.array([0, alternative[1], 0, 0, 0, 0, 0, 0, 0, 0, 0, alternative[11]])
    select3 = np.array([0, 0, alternative[2], alternative[3], alternative[4], alternative[5], alternative[6], 0, alternative[8], 0, 0, 0])
    return select1, select2, select3

# розділення сигма на класи
def divide_on_classes(sigma_matrix):
    matrix = np.array(sigma_matrix)
    class1 = [[0] * 20 for i in range(20)]
    class2 = [[0] * 20 for i in range(20)]
    class3 = [[0] * 20 for i in range(20)]
    for i in range(0, 20):
        for j in range(0, 20):
            alternative = matrix[i][j]
            class1[i][j], class2[i][j], class3[i][j] = select_classes(alternative)

    return class1, class2, class3

# сортування критеріїв за спаданням
def fi_alternatives(matrix):
    alternatives = np.array(matrix)
    for i in range(0,20):
        alternative=alternatives[i]
        alternative[::-1].sort()
        alternatives[i] = alternative
    return alternatives

print_optimal_alternative(pareto(SV))
pareto(SV)
print_optimal_alternative(majority(SV))
majority(SV)
print_optimal_alternative(lexicographic(sort_sigma(SV)))
lexicographic(sort_sigma(SV))
print_optimal_alternative(berezovskiy())
berezovskiy()
print_optimal_alternative(podynovskiy(condition))
podynovskiy(condition)