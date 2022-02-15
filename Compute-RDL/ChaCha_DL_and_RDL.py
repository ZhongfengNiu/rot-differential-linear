import numpy as np
import math as ma
import  random
def carry(a,b,c):
    return int((a & b) ^ (a & c) ^ (b & c))

def numtotwo(a,n):
    m=np.zeros(n,dtype=int)
    for i in range(0,n):
        m[n-1-i]=a%2
        a=a>>1
    return m

def M(out_1, in_1, in_2,out_2,in_3,out_3): #2模加旁枝，1模加输出枝
    num=4
    M=np.mat(np.zeros((4,4),dtype=float))
    for j in range(0,num):
        m=numtotwo(j,2)
        for y in range(0,4):
            y_1=numtotwo(y,2)
            c_1=carry(y_1[0],y_1[1],m[1])
            c_2=carry(y_1[0]^in_1,y_1[1]^in_2,m[0])
            value=(out_1&(m[0]^m[1]^in_1^in_2))^(out_2&in_2)^(out_3&in_3)
            t=c_1+c_2*2
            M[t,j]= M[t,j]+ma.pow((-1),(value))
    return  M*ma.pow(2,-2)

def generate(out_1,out_2,out_3,S,c_1,c_2,c_3,index):#2模加旁枝，1模加输出枝 #额外枝
    M = np.mat(np.zeros((4, 4), dtype=float))
    for i in range(0,8):
         i_5=int(i)&(0x1)
         i_4=(int(i)>>1)&(0x1)
         i_3=(int(i)>>2)&(0x1)
         M=M+(S[i_5+i_4*2+i_3*4+out_3*8+out_2*16+out_1*32]*0.5*(1+ma.pow(-1,i_5)*c_1[index])*0.5*(1+ma.pow(-1,i_4)*c_2[index])*0.5*(1+ma.pow(-1,i_3)*c_3[index]))
    return  M;

def Mall_dl(out_1,out_2,out_3,S,c_1,c_2,c_3):
    len=32
    L_1=np.mat(([1.0],[0],[0],[0]))
    C_1=np.mat(([1,1,1,1]))
    for i in range(0,len):
        ii=int(i)
        O_1=(int(out_1)>>ii)&(0x1)
        O_2=(int(out_2)>>ii)&(0x1)
        O_3=(int(out_3)>>ii)&(0x1)
        M=generate(O_1,O_2,O_3,S,c_1,c_2,c_3,i)
        L_1=np.dot(M,L_1)
    return  (np.dot(C_1,L_1));

def shif_left(a,t):#32位数
    b=0xffffffff
    return (a>>int(32-t))|(a<<int(t)&b);
def shif_right(a,t):
    b = 0xffffffff
    return (a>>int(t))|(a<<int(32-t)&b);

def calculate_dl(c_1,c_2,c_3,c_4,S,r_1,r_2):#2模加旁枝，1模加输出枝
   cc_1=np.ones((32))
   cc_2=np.ones((32))
   cc_3=np.ones((32))
   cc_4=np.ones((32))
   for i in range(0,32):
        bb=(0x1<<int(i))^0x00000000
        cc_1[i]=Mall_dl(bb,0x00000000,0x00000000,S,c_1,c_2,c_4)
        b_1 = shif_right(bb, r_1)
        cc_4[i]=Mall_dl(b_1,0x00000000,b_1,S,c_1,c_2,c_4)
   for i in range(0,32):
        bb=(0x1<<int(i))^0x00000000
        cc_3[i]=Mall_dl(bb,0x00000000,0x00000000,S,c_3,c_4,c_2)
        b_1 = shif_right(bb, r_2)
        cc_2[i]=Mall_dl(b_1,0x00000000,b_1,S,c_3,c_4,c_2)
   return [cc_1,cc_2,cc_3,cc_4]

def Calculate_dl_even(CC,S,r_1,s_1):# 16,12; 8,7;
    [CC[0],CC[5],CC[10],CC[15]]=calculate_dl(CC[0],CC[5],CC[10],CC[15], S, r_1, s_1)
    [CC[1], CC[6], CC[11], CC[12]] = calculate_dl(CC[1], CC[6], CC[11], CC[12], S, r_1, s_1)
    [CC[2], CC[7], CC[8], CC[13]] = calculate_dl(CC[2], CC[7], CC[8], CC[13], S, r_1, s_1)
    [CC[3], CC[4], CC[9], CC[14]] = calculate_dl(CC[3], CC[4], CC[9], CC[14], S, r_1, s_1)
    return  CC
def Calculate_dl_odd(CC,S,r_1,s_1):# 16,12; 8,7;
    [CC[0],CC[4],CC[8],CC[12]]=calculate_dl(CC[0],CC[4],CC[8],CC[12], S, r_1, s_1)
    [CC[1], CC[5], CC[9], CC[13]] = calculate_dl(CC[1], CC[5], CC[9], CC[13], S, r_1, s_1)
    [CC[2], CC[6], CC[10], CC[14]] = calculate_dl(CC[2], CC[6], CC[10], CC[14], S, r_1, s_1)
    [CC[3], CC[7], CC[11], CC[15]] = calculate_dl(CC[3], CC[7], CC[11], CC[15], S, r_1, s_1)
    return  CC
def Base_DL(CC,S,start,round,i_1,j_1):
    cc=(np.ones((16,32)))
    for i in range(0,16):
        for j in range(0,32):
            oo=(int(CC[i]) >> j) & (0x1)
            cc[i][j]=ma.pow(-1,oo)
    for m in range(start,start+round):
        if(m%4==0):
            cc=Calculate_dl_odd(cc, S, 16, 12)
            #print(cc)
        if(m%4==1):
            cc=Calculate_dl_odd(cc, S, 8, 7)
        if(m%4==2):
            cc=Calculate_dl_even(cc, S,16, 12)
        if(m%4==3):
            cc=Calculate_dl_even(cc, S, 8, 7)
    return cc[i_1][j_1]

def Mall_rdl(out_1,out_2,out_3,S,c_1,c_2,c_3,t):
    len=32
    L_1=np.mat(([1.0],[0],[0],[0]))
    L_2=np.mat(([0],[1.0],[0],[0]))
    C_1=np.mat(([1,0,1,0]))
    C_2=np.mat(([0,1,0,1]))
    MC=np.mat(([1.0,1.0,0,0],[0,0,0,0],[0,0,1.0,1.0],[0,0,0,0]))
    for i in range(t,len):
        ii=int(i)
        O_1=(int(out_1)>>ii)&(0x1)
        O_2=(int(out_2)>>ii)&(0x1)
        O_3=(int(out_3)>>ii)&(0x1)
        M=generate(O_1,O_2,O_3,S,c_1,c_2,c_3,i)
        L_1=np.dot(M,L_1)
        L_2 = np.dot(M, L_2)
    L_1=np.dot(MC,L_1)
    L_2 = np.dot(MC, L_2)
    for i in range(0,t):
        ii=int(i)
        O_1=(int(out_1)>>ii)&(0x1)
        O_2=(int(out_2)>>ii)&(0x1)
        O_3=(int(out_3)>>ii)&(0x1)
        M=generate(O_1,O_2,O_3,S,c_1,c_2,c_3,i)
        L_1=np.dot(M,L_1)
        L_2 = np.dot(M, L_2)
    return  (np.dot(C_1,L_1)+np.dot(C_2,L_2));

def calculate_rdl(c_1,c_2,c_3,c_4,S,r_1,r_2,t):
   cc_1=np.ones((32))
   cc_2=np.ones((32))
   cc_3=np.ones((32))
   cc_4=np.ones((32))
   for i in range(0,32):
        bb=(0x1<<int(i))^0x00000000
        cc_1[i]=Mall_rdl(bb,0x00000000,0x00000000,S,c_1,c_2,c_4,t)
        b_1 = shif_right(bb, r_1)
        cc_4[i]=Mall_rdl(b_1,0x00000000,b_1,S,c_1,c_2,c_4,t)
   for i in range(0,32):
        bb=(0x1<<int(i))^0x00000000
        cc_3[i]=Mall_rdl(bb,0x00000000,0x00000000,S,c_3,c_4,c_2,t)
        b_1 = shif_right(bb, r_2)
        cc_2[i]=Mall_rdl(b_1,0x00000000,b_1,S,c_3,c_4,c_2,t)
   return [cc_1,cc_2,cc_3,cc_4]

def Calculate_rdl_even(CC,S,r_1,s_1,t):# 16,12; 8,7;
    [CC[0],CC[5],CC[10],CC[15]]=calculate_rdl(CC[0],CC[5],CC[10],CC[15], S, r_1, s_1,t)
    [CC[1], CC[6], CC[11], CC[12]] = calculate_rdl(CC[1], CC[6], CC[11], CC[12], S, r_1, s_1,t)
    [CC[2], CC[7], CC[8], CC[13]] = calculate_rdl(CC[2], CC[7], CC[8], CC[13], S, r_1, s_1,t)
    [CC[3], CC[4], CC[9], CC[14]] = calculate_rdl(CC[3], CC[4], CC[9], CC[14], S, r_1, s_1,t)
    return  CC
def Calculate_rdl_odd(CC,S,r_1,s_1,t):# 16,12; 8,7;
    [CC[0],CC[4],CC[8],CC[12]]=calculate_rdl(CC[0],CC[4],CC[8],CC[12], S, r_1, s_1,t)
    [CC[1], CC[5], CC[9], CC[13]] = calculate_rdl(CC[1], CC[5], CC[9], CC[13], S, r_1, s_1,t)
    [CC[2], CC[6], CC[10], CC[14]] = calculate_rdl(CC[2], CC[6], CC[10], CC[14], S, r_1, s_1,t)
    [CC[3], CC[7], CC[11], CC[15]] = calculate_rdl(CC[3], CC[7], CC[11], CC[15], S, r_1, s_1,t)
    return  CC
def Base_RDL(CC,S,start,round,i_1,j_1,t):
    cc=(np.ones((16,32)))
    for i in range(0,16):
        for j in range(0,32):
            oo=(int(CC[i]) >> j) & (0x1)
            cc[i][j]=ma.pow(-1,oo)
    for m in range(start,start+round):
        if(m%4==0):
            cc=Calculate_rdl_odd(cc, S, 16, 12,t)
            #print(cc)
        if(m%4==1):
            cc=Calculate_rdl_odd(cc, S, 8, 7,t)
        if(m%4==2):
            cc=Calculate_rdl_even(cc, S,16, 12,t)
        if(m%4==3):
            cc=Calculate_rdl_even(cc, S, 8, 7,t)
    return cc[i_1][j_1]

S=[]
for i_0 in range(0,2):#out_1
    for i_1 in range(0,2):#out_2
        for i_2 in range(0,2):#out_3
            for i_3 in range(0,2):#in_3
                for i_4 in range(0,2):#in_2
                    for i_5 in range(0,2):#in_1
                        S.append(M(i_0, i_5, i_4,i_1,i_3,i_2))

CC=[0x00000000,0x00000000,0x00000000,0x00000000,
    0x00000000, 0x00000000, 0x00000000,0x00000000,
    0x00000000,0x00000000,0x00000000,0x00000000,
    0x00000000,0x00002000,0x00000000,0x00000000]
print(Base_DL(CC,S,0,6,11,0))

CC=[0x00000004,0x00000000,0x00000000,0x00000000,
    0x02202002,0x00000000, 0x00000000, 0x00000000,
    0x00400404,0x00000000,0x00000000,0x00000000,
    0x00400004,0x00000000,0x00000000,0x00000000]
print(Base_DL(CC,S,2,5,1,0))

CC=[0x00000000,0x00000004,0x00000000,0x00000000,
    0x00000000, 0x02202002,0x00000000, 0x00000000,
    0x00000000,0x00400404,0x00000000,0x00000000,
    0x00000000,0x00400004,0x00000000,0x00000000]
print(Base_DL(CC,S,2,5,2,0))

CC=[0x00000000,0x00000000,0x00000004,0x00000000,
    0x00000000, 0x00000000,0x02202002,0x00000000,
    0x00000000,0x00000000,0x00400404,0x00000000,
    0x00000000,0x00000000,0x00400004,0x00000000]
print(Base_DL(CC,S,2,5,3,0))

CC=[0x00000000,0x00000000,0x00000000,0x00000004,
    0x00000000, 0x00000000, 0x00000000,0x02202002,
    0x00000000,0x00000000,0x00000000,0x00400404,
    0x00000000,0x00000000,0x00000000,0x00400004]
print(Base_DL(CC,S,2,5,0,0))

CC=[0x00000000,0x00000000,0x00000000,0x00000000,
    0x00000000, 0x00000000, 0x00000000,0x00000000,
    0x00000000,0x00000000,0x00000000,0x00000000,
    0x00000000,0x00000000,0x00000000,0x00000000]
print(Base_RDL(CC,S,2,2,15,1,1))

CC=[0x00000000,0x00000000,0x00000000,0x00000000,
    0x00000000, 0x00000000,0x00000000,0x00000000,
    0x00000000,0x00000000,0x00000000,0x00000000,
    0x00000000,0x00000000,0x00000000,0x20000000]
print(Base_DL(CC,S,0,5,10,0))