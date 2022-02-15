#include<iostream>
#include<cstdlib>
#include<cmath>
#include<cstdint>
#include<random>
#include<array>
#define TYPE uint32_t   //define a new type for uint32_t
#define ROL(x,t) (((x<<t)|(x>>(32-t)))&0xffffffff)
using namespace std;

//////////////////////////////////////////////// random number generator
uint64_t p = 0x100000000;
random_device rd;
mt19937 gen(rd());
uniform_int_distribution<> dis(0,p-1);
uint64_t rand_64_bit(){
  TYPE left = dis(gen);
  TYPE right = dis(gen);
  return ((left));
}

void ChaCha(array<TYPE,16> &state,int start, int Round)
        {
         int  a_1[4][4]={{0,4,8,12},{1,5,9,13},{2,6,10,14},{3,7,11,15}};
         int  a_2[4][4]={{0,5,10,15},{1,6,11,12},{2,7,8,13},{3,4,9,14}};
           TYPE z1,z2,z3,z4,w1,w2;
            for(int i=start;i<start+Round;i++)
            {
                if (i%2==0)
                {
                    for(int j=0;j<4;j++)
                    {
                        z1 = (state[a_1[j][0]]+state[a_1[j][1]])&0xffffffff;
                        w1=state[a_1[j][3]]^z1;
                        z4=  ROL(w1,16);
                        z3 = (state[a_1[j][2]]+z4)&0xffffffff;
                        w2=state[a_1[j][1]]^z3;
                        z2=  ROL(w2,12);
                        state[a_1[j][0]] = (z1+z2)&0xffffffff;
                        w1=z4^state[a_1[j][0]];
                        state[a_1[j][3]]=ROL(w1,8);
                        state[a_1[j][2]] = (z3+state[a_1[j][3]])&0xffffffff;
                        w2=z2^state[a_1[j][2]];
                        state[a_1[j][1]]=ROL(w2,7);
                    }

                }
                if(i%2==1)
                 {
                    for(int j=0;j<4;j++)
                    {
                        z1 = (state[a_2[j][0]]+state[a_2[j][1]])&0xffffffff;
                        w1=state[a_2[j][3]]^z1;
                        z4=  ROL(w1,16);
                        z3 = (state[a_2[j][2]]+z4)&0xffffffff;
                        w2=state[a_2[j][1]]^z3;
                        z2=  ROL(w2,12);
                        state[a_2[j][0]] = (z1+z2)&0xffffffff;
                        w1=z4^state[a_2[j][0]];
                        state[a_2[j][3]]=ROL(w1,8);
                        state[a_2[j][2]] = (z3+state[a_2[j][3]])&0xffffffff;
                        w2=z2^state[a_2[j][2]];
                        state[a_2[j][1]]=ROL(w2,7);
                    }

                 }
            }
        }

int main(){
  int number_of_round = 4;
 TYPE diff_in[16]={0x00000000, 0x00000000, 0x00000000,0x00000000,
                   0x00000000, 0x00000000, 0x00000000,0x00000000,
                   0x00000000, 0x00000000, 0x00000000,0x00000000,
                   0x00000000, 0x00000000, 0x00000000,0x20000000};

 TYPE diff_out[16]={0x00000000, 0x00000000, 0x00000000,0x00000000,
                    0x00000000, 0x00000000, 0x00000000,0x00000000,
                    0x00000000, 0x00000000, 0x00000000,0x00000000,
                    0x00000000, 0x00000000, 0x00000000,0x00000000};
  long int counter = 0;
  long int counter_prf = 0;
  long int data_size = pow(2,25);
  array<TYPE,16> state;
   array<TYPE,16>state_diff;
  for(int exper = 0; exper <= (data_size - 1); exper++)
 {
      state[0]=0x61707865;
      state[1]=0x3320646e;
      state[2]=0x79622d32;
      state[3]=0x6b206574;
      for(int i=0;i<4;i++)
      {
          state_diff[i]=state[i]^diff_in[i];
      }
      for(int i=4;i<16;i++)
      {
          state[i]=rand_64_bit()&0xffffffff;
          state_diff[i]=state[i]^diff_in[i];
      }
        ChaCha(state,0,number_of_round);
        ChaCha(state_diff,0,number_of_round);
     for(int i=0;i<16;i++)
      {
          diff_out[i]=state[i]^state_diff[i];
      }
    TYPE output_linear_approximation =((diff_out[10])&0x1)^((diff_out[0])&0x1)^((diff_out[15])&0x1)^((diff_out[15]>>8)&0x1)
    ^((diff_out[14]>>24)&0x1)^((diff_out[3])&0x1)^((diff_out[3]>>16)&0x1)^((diff_out[4]>>7)&0x1)^((diff_out[9])&0x1);
    if(output_linear_approximation == 0)
    {
      counter++;
    }
  }
 float dd=(float)counter/data_size;
  cout << "permutation" << endl;
  cout << counter <<endl;
  cout << "The probability is " << (float)counter/data_size << endl;
  cout << "The correlation is  " << (2.0*dd-1) << endl;
  return 0;
}

