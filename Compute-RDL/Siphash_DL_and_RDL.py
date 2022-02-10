import numpy as np
import math as ma
def carry(a,b,c):
    return int((a & b) ^ (a & c) ^ (b & c))
def numtotwo(a,n):
    m=np.zeros(n,dtype=int)
    for i in range(0,n):
        m[n-1-i]=a%2
        a=a>>1
    return m
def M(out_1, in_1, in_2,out_2): #2模加旁枝，1模加输出枝
    num=4
    M=np.mat(np.zeros((4,4),dtype=float))
    for j in range(0,num):
        m=numtotwo(j,2)
        for y in range(0,4):
            y_1=numtotwo(y,2)
            c_1=carry(y_1[0],y_1[1],m[1])
            c_2=carry(y_1[0]^in_1,y_1[1]^in_2,m[0])
            value=(out_1&(m[0]^m[1]^in_1^in_2))^(out_2&in_2)
            t=c_1+c_2*2
            M[t,j]= M[t,j]+ma.pow((-1),(value))
    return  M*ma.pow(2,-2)

def generate(out_1,out_2,S,c_1,c_2,index):#2模加旁枝，1模加输出枝
    M = np.mat(np.zeros((4, 4), dtype=float))
    for i in range(0,4):
         i_4=int(i)&(0x1)
         i_3=int(i)>>1
         M=M+(S[i_4+i_3*2+out_2*4+out_1*8]*0.5*(1+ma.pow(-1,i_3)*c_1[index])*0.5*(1+ma.pow(-1,i_4)*c_2[index]))
    return  M;

def Mall_dl(out_1,out_2,S,c_1,c_2):
    len=64
    L=np.mat(([1.0],[0],[0],[0]))
    C=np.mat(([1,1,1,1]))
    for i in range(0,len):
        ii=int(i)
        O_1=(int(out_1)>>ii)&(0x1)
        O_2=(int(out_2)>>ii)&(0x1)
        M=generate(O_1,O_2,S,c_1,c_2,i)
        L=np.dot(M,L)
    return  (np.dot(C,L));

def Mall_rdl(out_1,out_2,S,c_1,c_2,t):
    len=64
    L_1=np.mat(([1.0],[0],[0],[0]))
    L_2=np.mat(([0],[1.0],[0],[0]))
    C_1=np.mat(([1,0,1,0]))
    C_2=np.mat(([0,1,0,1]))
    MC=np.mat(([1.0,1.0,0,0],[0,0,0,0],[0,0,1.0,1.0],[0,0,0,0]))
    for i in range(t,len):
        ii=int(i)
        O_1=(int(out_1)>>ii)&(0x1)
        O_2=(int(out_2)>>ii)&(0x1)
        M=generate(O_1,O_2,S,c_1,c_2,i)
        L_1=np.dot(M,L_1)
        L_2 = np.dot(M, L_2)
    L_1=np.dot(MC,L_1)
    L_2 = np.dot(MC, L_2)
    for i in range(0,t):
        ii=int(i)
        O_1=(int(out_1)>>ii)&(0x1)
        O_2=(int(out_2)>>ii)&(0x1)
        M=generate(O_1,O_2,S,c_1,c_2,i)
        L_1=np.dot(M,L_1)
        L_2 = np.dot(M, L_2)
    return  (np.dot(C_1,L_1)+np.dot(C_1,L_2));

def shif_left(a,t):#32位数
    b=0xffffffffffffffff
    return (a>>int(64-t))|(a<<int(t)&b);
def shif_right(a,t):
    b = 0xffffffffffffffff
    return (a>>int(t))|(a<<int(64-t)&b);

def calculate_dl(a,b,c,d,S,r_1,r_2,r_3):#2模加旁枝，1模加输出枝
   cc_a=np.ones((64))
   cc_b=np.ones((64))
   cc_c=np.ones((64))
   cc_d=np.ones((64))
   for i in range(0,64):
        cc_c[i]=Mall_dl(shif_right(((0x1<<int(i))^0x0000000000000000),r_3),0x0000000000000000,S,a,b)
        cc_a[i] = Mall_dl(((0x1 << int(i)) ^ 0x0000000000000000), 0x0000000000000000, S, c, d)
        bb=int(0x1<<int(i))
        b_2=shif_right(bb,r_2)
        b_1=shif_right(bb,r_1)
        cc_b[i]=Mall_dl(bb,b_2,S,a,b)
        cc_d[i] = Mall_dl(bb, b_1, S, c, d)
   return [cc_a,cc_b,cc_c,cc_d]

def calculate_rdl(a,b,c,d,S,r_1,r_2,r_3,t):#2模加旁枝，1模加输出枝
   cc_a=np.ones((64))
   cc_b=np.ones((64))
   cc_c=np.ones((64))
   cc_d=np.ones((64))
   for i in range(0,64):
        cc_c[i]=Mall_rdl(shif_right(((0x1<<int(i))^0x0000000000000000),r_3),0x0000000000000000,S,a,b,t)
        cc_a[i] = Mall_rdl(((0x1 << int(i)) ^ 0x0000000000000000), 0x0000000000000000, S, c, d,t)
        bb=int(0x1<<int(i))
        b_2=shif_right(bb,r_2)
        b_1=shif_right(bb,r_1)
        cc_b[i]=Mall_rdl(bb,b_2,S,a,b,t)
        cc_d[i] = Mall_rdl(bb, b_1, S, c, d,t)
   return [cc_a,cc_b,cc_c,cc_d]

def L(a,b,c,d,r_2,s_2):
    a_1=b^shif_right(c,32)
    b_1=shif_right(b,s_2)
    c_1=d^a
    d_1=shif_right(d,r_2)
    return [a_1,b_1,c_1,d_1]

def Siphash_DL(a_1,b_1,c_1,d_1,out_a,out_b,out_c,out_d,round,flag):
    a = np.ones((64))
    b = np.ones((64))
    c = np.ones((64))
    d = np.ones((64))
    for i in range(0, 64):
        O_1 = (int(a_1) >> i) & (0x1)
        O_2 = (int(b_1) >> i) & (0x1)
        O_3 = (int(c_1) >> i) & (0x1)
        O_4 = (int(d_1) >> i) & (0x1)
        a[i] = ma.pow(-1, O_1)
        b[i] = ma.pow(-1, O_2)
        c[i] = ma.pow(-1, O_3)
        d[i] = ma.pow(-1, O_4)
    if flag==1:
        r_1=16
        s_1=13
        r_2=21
        s_2=17
    if flag==2:
        r_1=21
        s_1=17
        r_2=16
        s_2=13

    for i in range(0,round-1):
        if int(i%2)==0:
            [a, b, c, d] = calculate_dl(a, b, c, d, S, r_1, s_1, 32)
        if int(i%2)==1:
            [a, b, c, d] = calculate_dl(a, b, c, d, S, r_2, s_2, 32)
    if round%2==0:
        [Out_a, Out_b, Out_c, Out_d]=L(out_a,out_b,out_c,out_d,r_2,s_2)
    if round%2==1:
        [Out_a, Out_b, Out_c, Out_d]=L(out_a,out_b,out_c,out_d,r_1,s_1)
    value = Mall_dl(Out_a, Out_b, S, a, b) * Mall_dl(Out_c, Out_d, S, c, d)
    return  value

def Siphash_RDL(a_1,b_1,c_1,d_1,out_a,out_b,out_c,out_d,round,flag,t):
    a = np.ones((64))
    b = np.ones((64))
    c = np.ones((64))
    d = np.ones((64))
    for i in range(0, 64):
        O_1 = (int(a_1) >> i) & (0x1)
        O_2 = (int(b_1) >> i) & (0x1)
        O_3 = (int(c_1) >> i) & (0x1)
        O_4 = (int(d_1) >> i) & (0x1)
        a[i] = ma.pow(-1, O_1)
        b[i] = ma.pow(-1, O_2)
        c[i] = ma.pow(-1, O_3)
        d[i] = ma.pow(-1, O_4)
    if flag==1:
        r_1=16
        s_1=13
        r_2=21
        s_2=17
    if flag==2:
        r_1=21
        s_1=17
        r_2=16
        s_2=13

    for i in range(0,round-1):
        if int(i%2)==0:
            [a, b, c, d] = calculate_rdl(a, b, c, d, S, r_1, s_1, 32,t)
        if int(i%2)==1:
            [a, b, c, d] = calculate_rdl(a, b, c, d, S, r_2, s_2, 32,t)
    if round%2==0:
        [Out_a, Out_b, Out_c, Out_d]=L(out_a,out_b,out_c,out_d,r_2,s_2)
    if round%2==1:
        [Out_a, Out_b, Out_c, Out_d]=L(out_a,out_b,out_c,out_d,r_1,s_1)
    value = Mall_rdl(Out_a, Out_b, S, a, b,t) * Mall_rdl(Out_c, Out_d, S, c, d,t)
    return  value

S=[]
for i_1 in range(0,2):
    for i_2 in range(0,2):
        for i_3 in range(0,2):
            for i_4 in range(0,2):
                S.append(M(i_1, i_3, i_4,i_2))


a_1=0x8000000000000000
b_1=0x0000000000000000
c_1=0x0000000000000000
d_1=0x0000000000000000
out_a=0x400000004000
out_b=0x400000004000
out_c=0x400000004000
out_d=0x400000004000
print(Siphash_DL(a_1,b_1,c_1,d_1,out_a,out_b,out_c,out_d,6,1))

a_1=0x0000000000000000
b_1=0x0000000000000000
c_1=0x8000000000000000
d_1=0x0000000000000000
out_a=0x2000000020000000
out_b=0x2000000020000000
out_c=0x2000000020000000
out_d=0x2000000020000000
print(Siphash_DL(a_1,b_1,c_1,d_1,out_a,out_b,out_c,out_d,7,2))



