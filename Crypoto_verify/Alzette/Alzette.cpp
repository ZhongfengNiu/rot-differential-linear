#include <iostream>
#include<iostream>
#include<cstdlib>
#include<cmath>
#include<cstdint>
#include<random>
#include<array>
#define TYPE uint32_t   //define a new type for uint64_t
#define ROL(x,t) (((x<<t)|(x>>(32-t)))&0xffffffff)
#define ROR(x,t) (((x>>t)|(x<<(32-t)))&0xffffffff)
using namespace std;

//////////////////////////////////////////////// random number generator
uint64_t p = 0x100000000;
random_device rd;
mt19937 gen(rd());
uniform_int_distribution<> dis(0,p-1);
uint64_t rand_64_bit(){
  TYPE left = dis(gen);
  //TYPE right = dis(gen);
  return (left);
}
void Alzette(TYPE a_in, TYPE b_in,  TYPE &a_out, TYPE &b_out,  int rd_num)
{
       TYPE c_0=0xB7E15162;
       TYPE c_1=0xBF715880;
       TYPE c_2=0x38B4DA56;
       TYPE c_3=0x324E7738;
       TYPE c_4=0xBB1185EB;
       TYPE c_5=0x4F7C7B57;
       TYPE c_6=0xCFBFA1C8;
       TYPE c_7=0xC2B3293D;
      TYPE a = a_in; TYPE b = b_in;
      TYPE key;
      array<int,8> rr={31,17,0,24,31,17,0,24};
      array<int,8> ll={24,17,31,16,24,17,31,16};
      for(int i=0;i<rd_num;i++)
      {
            if(i<4)
            {
               key=c_0;//Chosen the round constant;
            }
            else
            {
               key=c_2;//Chosen the round constant;
            }
            a=(a+ROR(b,rr[i]))&0xffffffff;
            b=b^ROR(a,ll[i]);
            a=a^key;
      }
        a_out = a;
        b_out = b;
}
bool maskdot(TYPE value, TYPE mask){
  bool prod = 0;
  for(int i=0; i<=31; i++){
    prod = prod ^ (((value>>i)&0x1)&((mask>>i)&0x1));
  }
  return prod;
}
int main()
{

   int  t;
    t=0; //Set the rotational parameters is 0, the D-L.

///////////////////////////////////////////////////////

    // The 8 round Alzette D-L distinguishers:
    // int number_of_round = 8;
    // TYPE diff_in_a =0x80020100;
    // TYPE diff_in_b =0x00010080;
    // TYPE mask_a = 0x80000000;
    // TYPE mask_b = 0x8000;

///////////////////////////////////////////////////////

  // The 5 round Alzette D-L distinguishers:
   int number_of_round = 5;
   TYPE diff_in_a =0x80000000;
   TYPE diff_in_b =0x0;
   TYPE mask_a = 0x80;
   TYPE mask_b = 0x8000;

////////////////////////////////////////////////////////

 // The 6 round Alzette D-L distinguishers:
 // int number_of_round = 6;
 // TYPE diff_in_a =0x80000000;
 // TYPE diff_in_b =0x0;
 // TYPE mask_a = 0x40;
 // TYPE mask_b = 0x200000;

//////////////////////////////////////////////////////
 // The 4 round Alzette D-L 2 distinguishers :
  //int number_of_round = 4;
  //TYPE diff_in_a =0x80000000;
  //TYPE diff_in_b =0x0;
  //TYPE mask_a = 0x80000;
  //TYPE mask_b = 0x8;

  //int number_of_round = 4;
  //TYPE diff_in_a =0x80000000;
  //TYPE diff_in_b =0x0;
  //TYPE mask_a = 0xc0;
  //TYPE mask_b = 0xc00000;

////////////////////////////////////////////////////////////
    // t=30;   // Set the rotational parameters is 30, the RD-L.

   // The 4 round Alzette RD-L  distinguishers :
    //int number_of_round = 4;
    //TYPE diff_in_a =0x7ffffffc;
    //TYPE diff_in_b =0x3fffffff;
    //TYPE mask_a = 0x4000;
    //TYPE mask_b = 0x40000000;

  long int counter = 0;
  long int data_size = pow(2,26);

    for(int exper = 0; exper <= (data_size - 1); exper++)
 {
      //input value
    TYPE input_a = rand_64_bit()&0xffffffff;   //generate random inputs
    TYPE input_b = rand_64_bit()&0xffffffff;
    //input value with the above difference
    TYPE input_a_prime =ROL(input_a,t)^diff_in_a;
    TYPE input_b_prime =ROL(input_b,t)^diff_in_b;

    TYPE output_a = 0x0; TYPE output_b = 0x0;
    TYPE output_a_prime = 0x0; TYPE output_b_prime = 0x0;
    Alzette(input_a,input_b,output_a,output_b,number_of_round);
    Alzette(input_a_prime,input_b_prime,output_a_prime,output_b_prime,number_of_round);
    //compute the output difference
    TYPE diff_out_a = ROL(output_a,t) ^ output_a_prime;
    TYPE diff_out_b = ROL(output_b,t) ^ output_b_prime;

    TYPE output_linear_approximation = maskdot(diff_out_a,mask_a)^maskdot(diff_out_b,mask_b);
    if(output_linear_approximation == 0)
    {
      counter++;
    }
  }
   float dd=(float)counter/data_size;
  cout << "Alzette " <<number_of_round<<" round"<< endl;
  cout<<"Counter"<<endl;
  cout << counter <<endl;
  cout << "The probability is " << (float)counter/data_size << endl;
  cout << "The correlation is  " << (2.0*dd-1) << endl;
  return 0;
}

