#include <stdlib.h>
#include <math.h>

/*--------------------------------------------------------------------------*/  
double gauss(void){
  double rnd(void);

  return rnd() + rnd() + rnd() + rnd() + rnd() + rnd() +
    rnd() + rnd() + rnd() + rnd() + rnd() + rnd() - 6.0;
}

/*--------------------------------------------------------------------------*/  
double rnd(void){
  return (double)random() / (double) RAND_MAX;
}

