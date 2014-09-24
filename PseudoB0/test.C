#include "TF1.h"
int* fit(int a, int b){
static int test[2];
test[0] = a;
test[1] = b;
int *rt;
rt = test;
return rt;
}
void test(){
int* a = fit(10, 20);
std::cout<<a[0]<<std::endl;
int* b = fit(30, 40);
std::cout<<b[0]<<std::endl;
}
