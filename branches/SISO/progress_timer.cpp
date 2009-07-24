#include <itpp/itbase.h>

using namespace itpp;

#define MAX 10000000

int main() {
  Progress_Timer timer;

  timer.set_max(MAX);
  timer.progress(0.0);
  for(int i=0;i<MAX;i++)
  {
  	timer.progress(i+1);
   	//timer.progress((double)(i+1)/MAX);
  }
  timer.toc_print();
}
