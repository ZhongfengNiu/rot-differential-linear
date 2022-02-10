#include<iostream>
#include<cstdlib>
#include<cmath>
#include<cstdint>
#include<random>
#include<ctime>
#include<iomanip>
#define TYPE uint16_t   //define a new type for uint64_t
#define ROL(x,t) (((x<<t)|(x>>(16-t)))&0xffff)
#define ROR(x,t) (((x>>t)|(x<<(16-t)))&0xffff)
#include<array>
using namespace std;

//////////////////////////////////////////////// random number generator
uint32_t p = 0x10000;
random_device rd;
mt19937 gen(rd());
uniform_int_distribution<> dis(0,p-1);
uint64_t rand_64_bit(){
  TYPE left = dis(gen);
  return (left);
}
///////////////////////////////////////////////////
void speck_32(TYPE a_in, TYPE b_in,  TYPE &a_out, TYPE &b_out,  int rd_num, array<TYPE,11> &kk)
{
  TYPE a = a_in; TYPE b = b_in;
  TYPE z1,z2,w1,w2,z3,z4;
  for(int i=0;i<rd_num;i++)
  {
      a=ROR(a,7);
      a=(a+b)&0xffff;
      a=a^kk[i];
      b=ROL(b,2)^a;
  }
  a_out=a;
  b_out=b;
}
bool maskdot(TYPE value, TYPE mask){
  bool prod = 0;
  for(int i=0; i<=15; i++){
    prod = prod ^ (((value>>i)&0x1)&((mask>>i)&0x1));
  }
  return prod;
}
int main()
{
 double counter,sum;
long int data_size = 100; // The number of master key.
 double message_szie=pow(2,32); // The the full plaintext space

//////////////////////////////////////////
  // 10 round D-L distinguishers
  int number_of_round = 10;
  TYPE diff_in_a =0x0a20;
  TYPE diff_in_b =0x4205;
  TYPE mask_a = 0x5820;
  TYPE mask_b = 0x4020;
//////////////////////////////////////////
  // 9 round D-L distinguishers
  //int number_of_round = 9;
  //TYPE diff_in_a=0x0211;
  //TYPE diff_in_b=0x0A04;
  // TYPE mask_a = 0x5820;
  //TYPE mask_b = 0x4020;
//////////////////////////////////////////
  // 8 round D-L distinguishers
  // int number_of_round = 8;
  // TYPE diff_in_a=0x0211;
  // TYPE diff_in_b=0x0A04;
  // TYPE mask_a = 0x0008;
  // TYPE mask_b = 0x0008;


  for(int exper = 0; exper <= (data_size - 1); exper++)
 {
    array<TYPE,4> master_key;
    array<TYPE,11>kk;
    array<TYPE,14>l;
    cout << "Speck_32 " <<number_of_round<<" round"<<endl;
    cout<<"The master key is"<<endl;
    for(int i=0;i<4;i++)
    {
        master_key[i]=rand_64_bit()&0xffff;
        cout<<hex<<master_key[i]<<endl;
    }
    kk[0]=master_key[0];
    TYPE mid;
    for (int i=0;i<3;i++)
    {
        l[i]=master_key[i+1];
    }
    for(int i=0;i<11;i++)
    {
        mid=((ROR(l[i],7)+kk[i])&0xffff)^(i&0xffff);
        kk[i+1]=mid^(ROL(kk[i],2));
        l[i+3]=mid;
    }
        counter=0;
        for(int x=0; x<pow(2,16);x++)
        {
            for(int y=0x0; y<pow(2,16);y++)
            {
                    TYPE input_a = x&0xffff;   //generate random inputs
                    TYPE input_b = y&0xffff;
                    TYPE input_a_prime =(input_a)^diff_in_a;
                    TYPE input_b_prime = (input_b)^diff_in_b;
                    TYPE output_a = 0x0; TYPE output_b = 0x0;
                    TYPE output_a_prime = 0x0; TYPE output_b_prime = 0x0;
                    speck_32(input_a,input_b,output_a,output_b,number_of_round,kk);
                    speck_32(input_a_prime,input_b_prime,output_a_prime,output_b_prime,number_of_round,kk);
                      //compute the output difference
                    TYPE diff_out_a = (output_a) ^ output_a_prime;
                    TYPE diff_out_b = (output_b) ^ output_b_prime;
                    TYPE output_linear_approximation =maskdot(mask_a,diff_out_a)^maskdot(mask_b,diff_out_b);
                if(output_linear_approximation == 0)
                {
                     counter++;
                }
            }
        }
        cout<<"Counter "<<endl;
        cout<<fixed<<setprecision(13)<<counter<<endl;
        sum=sum+counter;
        double dd=counter/(message_szie);
        cout << "The probability is " << (double)counter/(message_szie) << endl;
        cout << "The correlation is  " << (2.0*dd-1) << endl;

  }
        double dd_1=sum/(pow(2,32)*100.0);
        cout << "Speck_32 " <<number_of_round<<" round"<<endl;
        cout << "The mean correlation is  " << (2.0*dd_1-1) << endl;

    return 0;
}
