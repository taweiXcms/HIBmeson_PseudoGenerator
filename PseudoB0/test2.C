#include <iostream>                                                                                                                                                                                         
#include <ctime>

using namespace std;

// function to generate and retrun random numbers.
int * getRandom(int a)
{
  static int  r[10];

  // set the seed
  srand( (unsigned)time( NULL ) );
  for (int i = 0; i < 10; ++i)
  {
    r[i] = rand()+a;
    cout << r[i] << endl;
  }

  return r;
}

// main function to call above defined function.
int main ()
{
   // a pointer to an int.
   int *p;
   p = getRandom(3);
   for ( int i = 0; i < 10; i++ )
   {
       cout << "*(p + " << i << ") : ";
       cout << *(p + i) << endl;
   }
   int *q;
   q = getRandom(5);
   for ( int i = 0; i < 10; i++ )
   {
       cout << "*(q + " << i << ") : ";
       cout << *(q + i) << endl;
   }

   return 0;
}
