import numpy as np
import math as ma
def carry(a, b, c):
    return int((a & b) ^ (a & c) ^ (b & c))

def numtotwo(a, n):
    m = np.zeros(n, dtype=int)
    for i in range(0, n):
        m[n - 1 - i] = a % 2
        a = a >> 1
    return m

def M(out_1, in_1, in_2, out_2):  # 2模加旁枝，1模加输出枝
    num = 4
    M = np.mat(np.zeros((4, 4), dtype=float))
    for j in range(0, num):
        m = numtotwo(j, 2)
        for y in range(0, 4):
            y_1 = numtotwo(y, 2)
            c_1 = carry(y_1[0], y_1[1], m[1])
            c_2 = carry(y_1[0] ^ in_1, y_1[1] ^ in_2, m[0])
            value = (out_1 & (m[0] ^ m[1] ^ in_1 ^ in_2)) ^ (out_2 & in_2)
            t = c_1 + c_2 * 2
            M[t, j] = M[t, j] + ma.pow((-1), (value))
    return M * ma.pow(2, -2)


def generate(out_1, out_2, S, c_1, c_2, index):  # 2模加旁枝，1模加输出枝
    M = np.mat(np.zeros((4, 4), dtype=float))
    for i in range(0, 4):
        i_4 = int(i) & (0x1)
        i_3 = int(i) >> 1
        M = M + (S[i_4 + i_3 * 2 + out_2 * 4 + out_1 * 8] * 0.5 * (1 + ma.pow(-1, i_3) * c_1[index]) * 0.5 * (
                    1 + ma.pow(-1, i_4) * c_2[index]))
    return M;


def Mall_dl(out_1, out_2, S, c_1, c_2):
    len = 16
    L = np.mat(([1.0], [0], [0], [0]))
    C = np.mat(([1, 1, 1, 1]))
    for i in range(0, len):
        ii = int(i)
        O_1 = (int(out_1) >> ii) & (0x1)
        O_2 = (int(out_2) >> ii) & (0x1)
        M = generate(O_1, O_2, S, c_1, c_2, i)
        L = np.dot(M, L)
    return (np.dot(C, L));

def shif_left(a, t):  # 32位数
    b = 0xffff
    return (a >> int(16 - t)) | (a << int(t) & b);

def shif_right(a, t):
    b = 0xffff
    return (a >> int(t)) | (a << int(16 - t) & b);

def calculate_dl(c_1, c_2, S):  # 2模加旁枝，1模加输出枝
    cc_1 = np.ones((16))
    cc_2 = np.ones((16))
    for i in range(0, 16):
        b = int(0x1 << int(i))
        b_1 = shif_left(b, 7)
        b_2 = shif_right(b, 2)
        cc_1[i] = Mall_dl(b_1, 0x0000, S, c_1, c_2)
        cc_2[i] = Mall_dl(b, b_2, S, c_1, c_2)
    return [cc_1, cc_2]

def L(a_1, a_2):  # a_1 模加 a_2旁枝
    b = (shif_right(a_2, 2))
    a = a_1^a_2
    return [a, b]

S = []
for i_1 in range(0, 2):
    for i_2 in range(0, 2):
        for i_3 in range(0, 2):
            for i_4 in range(0, 2):
                S.append(M(i_1, i_3, i_4, i_2))

def Speck_32_dl(a_1,a_2,out_1,out_2,round,S):
    a_1 = shif_right(a_1, 7)
    c_1 = np.ones((16))
    c_2 = np.ones((16))
    for i in range(0, 16):
        O_1 = (int(a_1) >> i) & (0x1)
        O_2 = (int(a_2) >> i) & (0x1)
        c_1[i] = ma.pow(-1, O_1)
        c_2[i] = ma.pow(-1, O_2)
    for i in range(0,round-1):
        [c_1, c_2] = calculate_dl(c_1, c_2, S)
    [a, b] = L(out_1, out_2)
    value = Mall_dl(a, b, S, c_1, c_2)
    return value
a_1 = 0x8100  # 0x8040
a_2 = 0x8102  # 0x8140
out_1=0x0008
out_2=0x0008
print(Speck_32_dl(a_1,a_2,out_1,out_2,4,S))