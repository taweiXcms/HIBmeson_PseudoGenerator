#include "iostream"
#include "algorithm"
#include <stdlib.h>     /* qsort */
#include <TRandom3.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
using namespace std;
int compare (const void * a, const void * b)
{
  return ( *(int*)a - *(int*)b );
}
void quickSort(int arr[], int left, int right) {
      int i = left, j = right;
      int tmp;
      int pivot = arr[(left + right) / 2];
      /* partition */
      while (i <= j) {
            while (arr[i] < pivot)
                  i++;
            while (arr[j] > pivot)
                  j--;
            if (i <= j) {
                  tmp = arr[i];
                  arr[i] = arr[j];
                  arr[j] = tmp;
                  i++;
                  j--;
            }
      };
      /* recursion */
      if (left < j)
            quickSort(arr, left, j);
      if (i < right)
            quickSort(arr, i, right);
}
void test(){
	int arr[9] = {1, 12, 5, 26, 7, 14, 3, 7, 2};
	quickSort(arr, 0, 8);
	for(int i=0;i<9;i++){
		std::cout<<arr[i]<<std::endl;
	}
}

void ran(){
	TRandom3 *evtrandom = new TRandom3();
	int ran = 0;
	int myran[50000];
	for(int i=0; i<50000; i++){
		myran[i] = int(evtrandom->Uniform(500000000));
//		std::sort(std::begin(myran), std::end(myran));
	}
	quickSort(myran,0,49999);
//	qsort (myran, 50000, sizeof(int), compare);
	for(int i=0; i<50000; i++){
		std::cout<<myran[i]<<std::endl;
	}

//	for(int i = 0; i < 10000; i++){
//		ran = int(evtrandom->Uniform(5000000));
//		std::cout<<ran<<std::endl;
//	}  

//	int count = 0;
//	do{
//	    ran = int(evtrandom->Uniform(500000000));
//		std::cout<<ran<<std::endl;
//		count++;
//	}while(count<50000);
}
