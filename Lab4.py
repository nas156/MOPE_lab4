import math
import numpy as np
from scipy.stats import t,f
import random as r
from functools import partial
import prettytable as p


def make_exp(m):
    return [[r.randint(math.floor(Y_min), math.floor(Y_max)) for _ in range(m)] for i in range(N)]


def get_dispersion(y_aver, y):
    return sum([(i-y_aver)**2 for i in y])/len(y)


def get_average(y):
    return sum(y)/len(y)


def count_F(a, b):
    return max([a, b])/min([a, b])

def y_regr_norm(x1, x2, x3):
        return b0 + x1*b1 + x2*b2 + x3*b3

def y_regr_abs(x1, x2, x3):
        return a0 + a1*x1 + a2*x2 + a3*x3

def get_fisher_critical(prob,f3, f4):
    for i in [j*0.001 for j in range(int(10/0.001))]:
        if abs(f.cdf(i,f4,f3)-prob) < 0.0001:
            return i


def get_student_critical(prob, f3):
    for i in [j*0.0001 for j in range(int(5/0.0001))]:
        if abs(t.cdf(i,f3)-(0.5 + prob/0.1*0.05)) < 0.000005:
            return i

def get_cohren_critical(prob, f1, f2):
    f_crit = f.isf((1-prob)/f2, f1, (f2-1)*f1)
    return f_crit/(f_crit+f2-1)

    
def is_significant_coef(t, tkr):
    return t > tkr

def craete_table_norm():
        global table1, X_norm, Y_exp
        table1 = p.PrettyTable()
        for i in range(len(X_norm)):
            table1.add_column(X_types[i], X_norm[i])
        for i in range(m):
            table1.add_column("Y{0}".format(i+1), [j[i] for j in Y_exp])

def create_table_abs():
        global table2, X_abs, Y_exp
        table2 = p.PrettyTable()
        for i in range(len(X_abs)):
            table2.add_column(X_types[i], X_abs[i])
        for i in range(m):
            table2.add_column("Y{0}".format(i+1), [j[i] for j in Y_exp])

def get_aprox_value_abs(x1, x2, x3, x12, x13, x23, x123):
    return b_abs[0]+b_abs[1]*x1+b_abs[2]*x2+b_abs[3]*x3+b_abs[4]*x12+b_abs[5]*x13+b_abs[6]*x23+b_abs[7]*x123

def get_aprox_value_norm(x1, x2, x3, x12, x13, x23, x123):
    return b_norm[0]+b_norm[1]*x1+b_norm[2]*x2+b_norm[3]*x3+b_norm[4]*x12+b_norm[5]*x13+b_norm[6]*x23+b_norm[7]*x123

m = 3
prob = 0.95
X1_min = -20
X1_max = 30
X2_min = -25
X2_max = 10
X3_min = -25
X3_max = -20
k = 3
X_ranges = [[X1_min, X1_max], [X2_min, X2_max], [X3_min, X3_max]]
X_types = ["X1", "X2", "X3", "X12", "X13", "X23", "X123"]
Xcp_max = (X1_max + X2_max + X3_max) / 3
Xcp_min = (X1_min + X2_min + X3_min) / 3

Y_min = 200 + Xcp_min
Y_max = 200 + Xcp_max
flag0 = True
ct = 3 #кількість спроб
while(flag0):
    flag0 = False
    X0_norm = [1, 1, 1, 1]
    X1_norm = [-1, -1, 1, 1]
    X2_norm = [-1, 1, -1, 1]
    X3_norm = [-1, 1, 1, -1]
    X_norm = [X1_norm, X2_norm, X3_norm]
    N = len(X1_norm)
    
    X_abs = [[max(X_ranges[j]) if i == 1 else min(X_ranges[j]) for i in X_norm[j]] for j in range(k)]


    for i in range(len(X_abs)):
        print("Абсолютні Х{0}: {1}".format(i+1, X_abs[i]))


    Y_exp = make_exp(m)

    flag = True
    while(flag):
        table1 = p.PrettyTable()
        craete_table_norm()
        print("Нормалізована матриця:\n", table1)

        mx_norm_list = [get_average(i) for i in X_norm] 
        y_aver = [get_average(i) for i in Y_exp]
        my = get_average(y_aver)
        a1 = get_average([X_norm[0][i]*y_aver[i] for i in range(N)])
        a2 = get_average([X_norm[1][i]*y_aver[i] for i in range(N)])
        a3 = get_average([X_norm[2][i]*y_aver[i] for i in range(N)])
        a11 = get_average([X_norm[0][i]**2 for i in range(N)])
        a22 = get_average([X_norm[1][i]**2 for i in range(N)])
        a33 = get_average([X_norm[2][i]**2 for i in range(N)])
        a12 = get_average([X_norm[0][i]*X_norm[1][i] for i in range(N)])
        a13 = get_average([X_norm[0][i]*X_norm[2][i] for i in range(N)])
        a23 = get_average([X_norm[1][i]*X_norm[2][i] for i in range(N)])
        a21 = a12
        a31 = a13
        a32 = a23
        znam = np.array([[1, mx_norm_list[0], mx_norm_list[1],  mx_norm_list[2]],
                        [mx_norm_list[0], a11, a12, a13],
                        [mx_norm_list[1], a12, a22, a32],
                        [mx_norm_list[2], a13, a23, a33]])

        b0_matr = np.array([[my, mx_norm_list[0], mx_norm_list[1],  mx_norm_list[2]],
                        [a1, a11, a12, a13],
                        [a2, a12, a22, a32],
                        [a3, a13, a23, a33]])

        b1_matr = np.array([[1, my, mx_norm_list[1],  mx_norm_list[2]],
                        [mx_norm_list[0], a1, a12, a13],
                        [mx_norm_list[1], a2, a22, a32],
                        [mx_norm_list[2], a3, a23, a33]])

        b2_matr = np.array([[1, mx_norm_list[0], my,  mx_norm_list[2]],
                        [mx_norm_list[0], a11, a1, a13],
                        [mx_norm_list[1], a12, a2, a32],
                        [mx_norm_list[2], a13, a3, a33]])

        b3_matr = np.array([[1, mx_norm_list[0], mx_norm_list[1],  my],
                        [mx_norm_list[0], a11, a12, a1],
                        [mx_norm_list[1], a12, a22, a2],
                        [mx_norm_list[2], a13, a23, a3]])

        
        zanm_value = np.linalg.det(znam)
        b0 = np.linalg.det(b0_matr)/zanm_value
        b1 = np.linalg.det(b1_matr)/zanm_value
        b2 = np.linalg.det(b2_matr)/zanm_value
        b3 = np.linalg.det(b3_matr)/zanm_value

        
        
        print("Рівняння регресії для нормованих значень:\ny = {0} + {1}*x1 + {2}*x2 + {3}*x3".format(b0, b1, b2, b3))

        print("Підставимо нормованi значення Х в рівння регресії")
        print("Р-ня регресії для Х11, Х21, Х31 =",y_regr_norm(X_norm[0][0], X_norm[1][0], X_norm[2][0]))
        print("Середнє y1 =", y_aver[0])
        print("Р-ня регресії для Х12, Х22, Х32 =", y_regr_norm(X_norm[0][1], X_norm[1][1], X_norm[2][1]))
        print("Середнє y2 =", y_aver[1])
        print("Р-ня регресії для Х13, Х23, Х33 =", y_regr_norm(X_norm[0][2], X_norm[1][2], X_norm[2][2]))
        print("Середнє y3 =", y_aver[2])

        delt_x1 = (X1_max - X1_min)/2
        delt_x2 = (X2_max - X2_min)/2
        delt_x3 = (X3_max - X3_min)/2
        x10 = (X1_max + X1_min)/2
        x20 = (X2_max + X2_min)/2
        x30 = (X3_max + X3_min)/2
        a0 = b0 - b1*(x10/delt_x1) - b2*(x20/delt_x2) - b3*(x30/delt_x3)
        a1 = b1/delt_x1
        a2 = b2/delt_x2
        a3 = b3/delt_x3

        print("Рівняння регресії для абсолютних значень:\ny = {0} + {1}*x1 + {2}*x2 + {3}*x3".format(a0, a1, a2, a3))

        print("Підставимо абсолютні значення Х в рівння регресії")
        print("Р-ня регресії для Х11, Х21, Х31 =",y_regr_abs(X_abs[0][0], X_abs[1][0], X_abs[2][0]))
        print("Середнє y1 =", y_aver[0])
        print("Р-ня регресії для Х12, Х22, Х32 =", y_regr_abs(X_abs[0][1], X_abs[1][1], X_abs[2][1]))
        print("Середнє y2 =", y_aver[1])
        print("Р-ня регресії для Х13, Х23, Х33 =", y_regr_abs(X_abs[0][2], X_abs[1][2], X_abs[2][2]))
        print("Середнє y3 =", y_aver[2])
        print("Р-ня регресії для Х14, Х24, Х34 =", y_regr_abs(X_abs[0][3], X_abs[1][3], X_abs[2][3]))
        print("Середнє y3 =", y_aver[3])


        #Кохрен
        y_disps = [get_dispersion(y_aver[i], Y_exp[i]) for i in range(N)]
        f1 = m - 1 
        f2 = N
        f3 = f2*f1
        Gp = max(y_disps)/sum(y_disps)
        Gkr = get_cohren_critical(prob, f1, f2)
        print("--------------------------------------------------------")
        if(Gkr > Gp):
            print("Gkr = {0} > Gp = {1} ---> Дисперсії однорідні".format(Gkr, Gp))
            flag = False
        else:
            print("Gkr = {0} < Gp = {1} ---> Дисперсії неоднорідні, збільшимо m і проведемо розрахунки".format(Gkr, Gp))
            Y_exp[0].append(r.randint(math.floor(Y_min), math.floor(Y_max)))
            Y_exp[1].append(r.randint(math.floor(Y_min), math.floor(Y_max)))
            Y_exp[2].append(r.randint(math.floor(Y_min), math.floor(Y_max)))
            Y_exp[3].append(r.randint(math.floor(Y_min), math.floor(Y_max)))
            m += 1



    #Стьюдент
    S2B = sum(y_disps)/N
    S2b = S2B/(N*m)
    Sb = math.sqrt(S2b)
    beta0 = sum([y_aver[i]*X0_norm[i] for i in range(N)])/N
    beta1 = sum([y_aver[i]*X1_norm[i] for i in range(N)])/N
    beta2 = sum([y_aver[i]*X2_norm[i] for i in range(N)])/N
    beta3 = sum([y_aver[i]*X3_norm[i] for i in range(N)])/N
    t0 = abs(beta0)/Sb
    t1 = abs(beta1)/Sb
    t2 = abs(beta2)/Sb
    t3 = abs(beta3)/Sb
    tkr = get_student_critical(prob, f3)


    d = sum([1 if is_significant_coef(tkr, i) else 0 for i in [t0, t1, t2, t3]])

    a0 = a0 if is_significant_coef(t0, tkr) else 0
    a1 = a1 if is_significant_coef(t1, tkr) else 0
    a2 = a2 if is_significant_coef(t2, tkr) else 0
    a3 = a3 if is_significant_coef(t3, tkr) else 0

    y_new = [y_regr_abs(X_abs[0][i], X_abs[1][i], X_abs[2][i]) for i in range(N)]

    print("--------------------------------------------------------\nПісля перевірки значимості коефіцієнтів: ")
    print("Кількість значимих коефіцієнтів:", d)
    print("Рівняння регресії для абсолютних значень:\ny = {0} + {1}*x1 + {2}*x2 + {3}*x3".format(a0, a1, a2, a3))
    print("Підставимо абсолютні значення Х в рівння регресії")
    print("Р-ня регресії для Х11, Х21, Х31 =", y_new[0])
    print("Р-ня регресії для Х12, Х22, Х32 =", y_new[1])
    print("Р-ня регресії для Х13, Х23, Х33 =", y_new[2])
    print("Р-ня регресії для Х14, Х24, Х34 =", y_new[3])


    #Фішер
    print("--------------------------------------------------------")
    f4 = N - d
    S2ad = (m/(N-d))*sum([(y_new[i] - y_aver[i])**2 for i in range(N)])
    Fp = S2ad/S2b
    Fkr = get_fisher_critical(prob, f3, f4)
    if(Fkr > Fp):
        print("Fkr = {0} > Fp = {1} ---> Р-ня адекватне оригіналу".format(Fkr, Fp))
    else:
        print("Рівняння неадекватне переходимо до ефекту взаємодії")
        print("-"*45)
        X1_norm = [-1, -1, -1, -1, 1, 1, 1, 1]
        X2_norm = [-1, -1, 1, 1, -1, -1, 1, 1]
        X3_norm = [-1, 1, -1, 1, -1, 1, -1, 1]
        N = len(X1_norm)
        Xcp_max = (X1_max + X2_max + X3_max) / 3
        Xcp_min = (X1_min + X2_min + X3_min) / 3
        X_norm = [X1_norm, X2_norm, X3_norm]
        X_abs = [[max(X_ranges[j]) if i == 1 else min(X_ranges[j]) for i in X_norm[j]] for j in range(k)]
        

        Y_exp = make_exp(m)
        flag1 = True
        while(flag1):
            table1 = p.PrettyTable()

            X12_norm = [X1_norm[i]*X2_norm[i] for i in range(N)]
            X13_norm = [X1_norm[i]*X3_norm[i] for i in range(N)]
            X23_norm = [X2_norm[i]*X3_norm[i] for i in range(N)]
            X123_norm = [X1_norm[i]*X2_norm[i]*X3_norm[i] for i in range(N)]

            X_norm.append(X12_norm)
            X_norm.append(X13_norm)
            X_norm.append(X23_norm)
            X_norm.append(X123_norm)
            craete_table_norm()
            print("Нормалізована матриця:\n", table1)
            X_abs = [[max(X_ranges[j]) if i == 1 else min(X_ranges[j]) for i in X_norm[j]] for j in range(k)]
            X12_abs, X13_abs, X23_abs, X123_abs = [], [], [], []
            for i in range(N):
                if X1_norm[i] == 1:
                    x1 = max(X_ranges[0])
                else:
                    x1 = min(X_ranges[0])

                if X2_norm[i] == 1:
                    x2 = max(X_ranges[1])
                else:
                    x2 = min(X_ranges[1])
                X12_abs.append(x1*x2)

            X_abs.append(X12_abs)

            for i in range(N):
                if X1_norm[i] == 1:
                    x1 = max(X_ranges[0])
                else:
                    x1 = min(X_ranges[0])

                if X3_norm[i] == 1:
                    x3 = max(X_ranges[2])
                else:
                    x3 = min(X_ranges[2])
                X13_abs.append(x1*x3)
            X_abs.append(X13_abs)

            for i in range(N):
                if X3_norm[i] == 1:
                    x3 = max(X_ranges[2])
                else:
                    x3 = min(X_ranges[2])

                if X2_norm[i] == 1:
                    x2 = max(X_ranges[1])
                else:
                    x2 = min(X_ranges[1])
                X23_abs.append(x3*x2)   

            X_abs.append(X23_abs)

            for i in range(N):
                if X1_norm[i] == 1:
                    x1 = max(X_ranges[0])
                else:
                    x1 = min(X_ranges[0])

                if X2_norm[i] == 1:
                    x2 = max(X_ranges[1])
                else:
                    x2 = min(X_ranges[1])

                if X3_norm[i] == 1:
                    x3 = max(X_ranges[2])
                else:
                    x3 = min(X_ranges[2])
                X123_abs.append(x1*x2*x3)
            X_abs.append(X123_abs)
            table2 = p.PrettyTable()
            create_table_abs()
            print("Абсолютна матриця:\n", table2)

            X1_abs = X_abs[0]
            X2_abs = X_abs[1]
            X3_abs = X_abs[2]
            mj0 = [N]
            for i in [sum(i) for i in X_abs]:
                mj0.append(i)

            mj1 = [sum(X1_abs),
                    sum([i**2 for i in X1_abs]),
                    sum(X12_abs),
                    sum(X13_abs),
                    sum([X1_abs[i]**2*X2_abs[i] for i in range(N)]),
                    sum([X1_abs[i]**2*X3_abs[i] for i in range(N)]),
                    sum(X123_abs),
                    sum([X1_abs[i]**2*X2_abs[i]*X3_abs[i] for i in range(N)])]

            mj2 = [sum(X2_abs),
                    sum(X12_abs),
                    sum([i**2 for i in X2_abs]),
                    sum(X23_abs),
                    sum([X1_abs[i]*X2_abs[i]**2 for i in range(N)]),
                    sum(X123_abs),
                    sum([X3_abs[i]*X2_abs[i]**2 for i in range(N)]),
                    sum([X1_abs[i]*X2_abs[i]**2*X3_abs[i] for i in range(N)])]

            mj3 = [sum(X3_abs),
                    sum(X13_abs),
                    sum(X23_abs),
                    sum([i**2 for i in X3_abs]),
                    sum(X123_abs),
                    sum([X1_abs[i]*X3_abs[i]**2 for i in range(N)]),
                    sum([X3_abs[i]**2*X2_abs[i] for i in range(N)]),
                    sum([X1_abs[i]*X2_abs[i]*X3_abs[i]**2 for i in range(N)])
                    ]

            mj4 = [sum(X12_abs),
                    sum([X1_abs[i]**2*X2_abs[i] for i in range(N)]),
                    sum([X1_abs[i]*X2_abs[i]**2 for i in range(N)]),
                    sum(X123_abs),
                    sum([X1_abs[i]**2*X2_abs[i]**2 for i in range(N)]),
                    sum([X1_abs[i]**2*X2_abs[i]*X3_abs[i] for i in range(N)]),
                    sum([X1_abs[i]*X2_abs[i]**2*X3_abs[i] for i in range(N)]),
                    sum([X1_abs[i]**2*X2_abs[i]**2*X3_abs[i] for i in range(N)])
                    ]

            mj5 = [sum(X13_abs),
                    sum([X1_abs[i]**2*X3_abs[i] for i in range(N)]),
                    sum(X123_abs),
                    sum([X1_abs[i]*X3_abs[i]**2 for i in range(N)]),
                    sum([X1_abs[i]**2*X2_abs[i]*X3_abs[i] for i in range(N)]),
                    sum([X1_abs[i]**2*X3_abs[i]**2 for i in range(N)]),
                    sum([X1_abs[i]*X2_abs[i]*X3_abs[i]**2 for i in range(N)]),
                    sum([X1_abs[i]**2*X2_abs[i]*X3_abs[i]**2 for i in range(N)])
                    ]

            mj6 = [sum(X23_abs),
                    sum(X123_abs),
                    sum([X3_abs[i]*X2_abs[i]**2 for i in range(N)]),
                    sum([X3_abs[i]**2*X2_abs[i] for i in range(N)]),
                    sum([X1_abs[i]*X2_abs[i]**2*X3_abs[i] for i in range(N)]),
                    sum([X1_abs[i]*X2_abs[i]*X3_abs[i]**2 for i in range(N)]),
                    sum([X2_abs[i]**2*X3_abs[i]**2 for i in range(N)]),
                    sum([X1_abs[i]*X2_abs[i]**2*X3_abs[i]**2 for i in range(N)])
                    ]

            mj7 = [sum(X123_abs),
                    sum([X1_abs[i]**2*X2_abs[i]*X3_abs[i] for i in range(N)]),
                    sum([X1_abs[i]*X2_abs[i]**2*X3_abs[i] for i in range(N)]),
                    sum([X1_abs[i]*X2_abs[i]*X3_abs[i]**2 for i in range(N)]),
                    sum([X1_abs[i]**2*X2_abs[i]**2*X3_abs[i] for i in range(N)]),
                    sum([X1_abs[i]**2*X2_abs[i]*X3_abs[i]**2 for i in range(N)]),
                    sum([X1_abs[i]*X2_abs[i]**2*X3_abs[i]**2 for i in range(N)]),
                    sum([X1_abs[i]**2*X2_abs[i]**2*X3_abs[i]**2 for i in range(N)])
                    ]

            y_aver = [get_average(i) for i in Y_exp]

            m_list = [mj0, mj1, mj2, mj3, mj4, mj5, mj6, mj7]
            m0j = [i[0] for i in m_list]
            m1j = [i[1] for i in m_list]
            m2j = [i[2] for i in m_list]
            m3j = [i[3] for i in m_list]
            m4j = [i[4] for i in m_list]
            m5j = [i[5] for i in m_list]
            m6j = [i[6] for i in m_list]
            m7j = [i[7] for i in m_list]

            kf = [sum(y_aver),
                    sum([y_aver[i]*X1_abs[i] for i in range(N)]),
                    sum([y_aver[i]*X2_abs[i] for i in range(N)]),
                    sum([y_aver[i]*X3_abs[i] for i in range(N)]),
                    sum([y_aver[i]*X12_abs[i] for i in range(N)]),
                    sum([y_aver[i]*X13_abs[i] for i in range(N)]),
                    sum([y_aver[i]*X23_abs[i] for i in range(N)]),
                    sum([y_aver[i]*X123_abs[i] for i in range(N)])
                    ]

            b_abs = list(np.linalg.solve([mj0, mj1, mj2, mj3, mj4, mj5, mj6, mj7], kf))


            print("-"*45)
            print("Рівняння регресії для абсолютних значень: {0} + {1}*x1 + {2}*x2 + {3}*x3 + {4}*X12 + {5}*X13 + {6}*x23 + {7}*x123".format(b_abs[0],b_abs[1],b_abs[2],b_abs[3],b_abs[4],b_abs[5],b_abs[6],b_abs[7]))
            print("Перевіримо")
            print("Y1 середнє =", y_aver[0])
            a = [X_abs[i][0] for i in range(7)]
            print("Y1 з рівняння =", get_aprox_value_abs(*a))
            print("Y4 середнє =", y_aver[3])
            b = [X_abs[i][3] for i in range(7)]
            print("Y4 з рівняння =", get_aprox_value_abs(*b))

            b_norm = [get_average(y_aver),
                        get_average([y_aver[i]*X1_norm[i] for i in range(N)]),
                        get_average([y_aver[i]*X2_norm[i] for i in range(N)]),
                        get_average([y_aver[i]*X3_norm[i] for i in range(N)]),
                        get_average([y_aver[i]*X12_norm[i] for i in range(N)]),
                        get_average([y_aver[i]*X13_norm[i] for i in range(N)]),
                        get_average([y_aver[i]*X23_norm[i] for i in range(N)]),
                        get_average([y_aver[i]*X123_norm[i] for i in range(N)])
                        ]



            print("-"*45)
            print("Рівняння регресії для нормалізовних значень: {0} + {1}*x1 + {2}*x2 + {3}*x3 + {4}*X12 + {5}*X13 + {6}*x23 + {7}*x123".format(b_norm[0],b_norm[1],b_norm[2],b_norm[3],b_norm[4],b_norm[5],b_norm[6],b_norm[7]))
            print("Перевіримо")
            print("Y1 середнє =", y_aver[0])
            a = [X_norm[i][0] for i in range(7)]
            print("Y1 з рівняння =", get_aprox_value_norm(*a))
            print("Y4 середнє =", y_aver[3])
            b = [X_norm[i][3] for i in range(7)]
            print("Y4 з рівняння =", get_aprox_value_norm(*b))
            y_disps = [get_dispersion(y_aver[i], Y_exp[i]) for i in range(N)]

            f1 = m - 1 
            f2 = N
            f3 = f2*f1
            Gp = max(y_disps)/sum(y_disps)
            Gkr = get_cohren_critical(prob, f1, f2)
            print("--------------------------------------------------------")
            if(Gkr > Gp):
                print("Gkr = {0} > Gp = {1} ---> Дисперсії однорідні".format(Gkr, Gp))
                flag1 = False
            else:
                print("Gkr = {0} < Gp = {1} ---> Дисперсії неоднорідні, збільшимо m і проведемо розрахунки".format(Gkr, Gp))
                Y_exp[0].append(r.randint(math.floor(Y_min), math.floor(Y_max)))
                Y_exp[1].append(r.randint(math.floor(Y_min), math.floor(Y_max)))
                Y_exp[2].append(r.randint(math.floor(Y_min), math.floor(Y_max)))
                Y_exp[3].append(r.randint(math.floor(Y_min), math.floor(Y_max)))
                Y_exp[4].append(r.randint(math.floor(Y_min), math.floor(Y_max)))
                Y_exp[5].append(r.randint(math.floor(Y_min), math.floor(Y_max)))
                Y_exp[6].append(r.randint(math.floor(Y_min), math.floor(Y_max)))
                Y_exp[7].append(r.randint(math.floor(Y_min), math.floor(Y_max)))
                m += 1

        #Стьюдент
        S2B = sum(y_disps)/N
        S2b = S2B/(N*m)
        Sb = math.sqrt(S2b)
        beta0 = get_average(y_aver)
        beta1 = sum([y_aver[i]*X1_norm[i] for i in range(N)])/N
        beta2 = sum([y_aver[i]*X2_norm[i] for i in range(N)])/N
        beta3 = sum([y_aver[i]*X3_norm[i] for i in range(N)])/N
        beta4 = sum([y_aver[i]*X12_norm[i] for i in range(N)])/N
        beta5 = sum([y_aver[i]*X13_norm[i] for i in range(N)])/N
        beta6 = sum([y_aver[i]*X23_norm[i] for i in range(N)])/N
        beta7 = sum([y_aver[i]*X123_norm[i] for i in range(N)])/N
        t0 = abs(beta0)/Sb
        t1 = abs(beta1)/Sb
        t2 = abs(beta2)/Sb
        t3 = abs(beta3)/Sb
        t4 = abs(beta4)/Sb
        t5 = abs(beta5)/Sb
        t6 = abs(beta6)/Sb
        t7 = abs(beta7)/Sb
        tkr = get_student_critical(prob, f3)
    
        d = sum([1 if is_significant_coef(i, tkr) else 0 for i in [t0, t1, t2, t3, t4, t5, t6, t7]])
        b_abs[0] = b_abs[0] if is_significant_coef(t0, tkr) else 0
        b_abs[1] = b_abs[1] if is_significant_coef(t1, tkr) else 0
        b_abs[2] = b_abs[2] if is_significant_coef(t2, tkr) else 0
        b_abs[3] = b_abs[3] if is_significant_coef(t3, tkr) else 0
        b_abs[4] = b_abs[4] if is_significant_coef(t4, tkr) else 0
        b_abs[5] = b_abs[5] if is_significant_coef(t5, tkr) else 0
        b_abs[6] = b_abs[6] if is_significant_coef(t6, tkr) else 0
        b_abs[7] = b_abs[7] if is_significant_coef(t7, tkr) else 0

        y_new = [get_aprox_value_abs(X1_abs[i], X2_abs[i], X3_abs[i], X12_abs[i], X13_abs[i], X23_abs[i], X123_abs[i]) for i in range(N)]
        
        print("-"*45)
        print("Після перевірки значимості коефіцієнтів: ")
        print("Кількість значимих коефіцієнтів:", d)
        print("Рівняння регресії для абсолютних значень:\ny = {0} + {1}*x1 + {2}*x2 + {3}*x3 + {4}*x12 + {5}*x13 + {6}*x23 + {7}*x123".format(b_abs[0], b_abs[1], b_abs[2], b_abs[3], b_abs[4], b_abs[5], b_abs[6], b_abs[7]))
        print("Підставимо абсолютні значення Х в рівння регресії")
        print("Y1 середнє =", y_aver[0])
        a = [X_abs[i][0] for i in range(7)]
        print("Y1 з рівняння =", get_aprox_value_abs(*a))
        print("Y4 середнє =", y_aver[3])
        b = [X_abs[i][3] for i in range(7)]
        print("Y4 з рівняння =", get_aprox_value_abs(*b))
        
        #Фішер
        print("-"*45)
        f4 = N - d
        S2ad = abs((m/(N-d))*sum([(y_new[i] - y_aver[i])**2 for i in range(N)]))
        Fp = S2ad/S2b
        Fkr = get_fisher_critical(prob, f3, f4)
        
        if(Fkr > Fp):
            print("Fkr = {0} > Fp = {1} ---> Р-ня адекватне оригіналу".format(Fkr, Fp))
            flag0 = False
        else:
            print("Fkr = {0} < Fp = {1} ---> Р-ня неадекватне оригіналу".format(Fkr, Fp))
            flag0 = ct!=0
            ct-=1
