#include<iostream>
#include<cstdlib>
#include<cmath>
#include<cstdint>
#include<random>
#define TYPE uint64_t   //define a new type for uint64_t
#define ROL(x,t) (((x<<t)|(x>>(64-t)))&0xffffffffffffffff)
using namespace std;

//////////////////////////////////////////////// random number generator
uint32_t p = 0x100000000;
random_device rd;
mt19937 gen(rd());
uniform_int_distribution<> dis(0,p-1);
uint64_t rand_64_bit(){
  TYPE left = dis(gen);
  TYPE right = dis(gen);
  return ((left<<32)|right);
}
void siphash_permutation(TYPE a_in, TYPE b_in, TYPE c_in, TYPE d_in, TYPE &a_out, TYPE &b_out, TYPE &c_out, TYPE &d_out, int rd_num)
{  //siphash permutation
  TYPE a = a_in; TYPE b = b_in; TYPE c = c_in; TYPE d = d_in;
  TYPE z1,z2,w1,w2,z3,z4;
  for(int i=0; i<rd_num; i++){
    z1 = (a+b)&0xffffffffffffffff;
    z2 = (c+d)&0xffffffffffffffff;

    w1 = z1 ^ ROL(a,13);
    w2 = z2 ^ ROL(d,16);

    z3 = (w1+z2)&0xffffffffffffffff;
    z4 = (w2+(ROL(z1,32)))&0xffffffffffffffff;

    a = z3 ^ ROL(w1,17);
    b = z4;
    c = ROL(z3,32);
    d = z4 ^ ROL(w2,21);
  }
  a_out = a;
  b_out = b;
  c_out = c;
  d_out = d;
}
bool maskdot(TYPE value, TYPE mask){
  bool prod = 0;
  for(int i=0; i<=63; i++){
    prod = prod ^ (((value>>i)&0x1)&((mask>>i)&0x1));
  }
  return prod;
}
int main()
{
  int t;

///////////////////////////////////////////////////
  // The 3 round  siphash finalization D-L distinguisher
    t=0;
  int number_of_round = 3;
  TYPE diff_in_a = 0x0;
  TYPE diff_in_b = 0x8000000000000000;
  TYPE diff_in_c = 0x0;
  TYPE diff_in_d = 0x0;
  //output masks
  TYPE mask =0x400000004000;
  TYPE mask_a =0x400000004000;
  TYPE mask_b = 0x400000004000;
  TYPE mask_c =0x400000004000;
  TYPE mask_d =0x400000004000;

///////////////////////////////////////////////////
 //The 4 round  siphash finalization D-L distinguisher
 //    t=0;
//  int number_of_round = 4;
//  TYPE diff_in_a = 0x40000;
//  TYPE diff_in_b = 0x80040000;
//  TYPE diff_in_c = 0x0;
//  TYPE diff_in_d = 0x0;

  //output masks
  //TYPE mask =0x400000004000;
//  TYPE mask_a =0x2000000020000000;
//  TYPE mask_b = 0x2000000020000000;
//  TYPE mask_c =0x2000000020000000;
//  TYPE mask_d =0x2000000020000000;

///////////////////////////////////////////////
  long int counter = 0;
  long int data_size = pow(2,26);

  for(int exper = 0; exper <= (data_size - 1); exper++){
      //input value
    TYPE input_a = rand_64_bit()&0xfffffffffffffff;   //generate random inputs
    TYPE input_b = rand_64_bit()&0xfffffffffffffff;
    TYPE input_c = rand_64_bit()&0xfffffffffffffff;
    TYPE input_d = rand_64_bit()&0xfffffffffffffff;
    //input value with the above difference
    TYPE input_a_prime = ROL(input_a,t) ^ diff_in_a;
    TYPE input_b_prime = ROL(input_b,t) ^ diff_in_b;
    TYPE input_c_prime = ROL(input_c,t) ^ diff_in_c;
    TYPE input_d_prime = ROL(input_d,t) ^ diff_in_d;

    TYPE output_a = 0x0; TYPE output_b = 0x0; TYPE output_c = 0x0; TYPE output_d = 0x0;
    TYPE output_a_prime = 0x0; TYPE output_b_prime = 0x0; TYPE output_c_prime = 0x0; TYPE output_d_prime = 0x0;
    siphash_permutation(input_a,input_b,input_c,input_d,output_a,output_b,output_c,output_d,number_of_round);
    siphash_permutation(input_a_prime,input_b_prime,input_c_prime,input_d_prime,output_a_prime,output_b_prime,output_c_prime,output_d_prime,number_of_round);
    TYPE diff_out_a = ROL(output_a,t) ^ output_a_prime;
    TYPE diff_out_b = ROL(output_b,t) ^ output_b_prime;
    TYPE diff_out_c = ROL(output_c,t) ^ output_c_prime;
    TYPE diff_out_d = ROL(output_d,t) ^ output_d_prime;  //differences on four branches

    TYPE output_linear_approximation = maskdot(mask_a,diff_out_a) ^ maskdot(mask_b,diff_out_b) ^ maskdot(mask_c,diff_out_c) ^ maskdot(mask_d,diff_out_d);
    if(output_linear_approximation == 0){
      counter++;
    }
  }
  cout << "Siphash "<<number_of_round<<" round"<< endl;
  cout<<"Counter"<<endl;
  cout << counter <<endl;
  cout << "The probability is " << (float)counter/data_size << endl;
  cout << "The correlation is  " << ((2*((float)counter/data_size)-1));
  return 0;

}
