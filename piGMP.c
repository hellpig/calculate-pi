/*

Will output pi rounded to any number of digits in various bases.
I will use GMP to do this.

I will use the BBP formula that divides by 16^k each iteration...
https://en.wikipedia.org/wiki/List_of_formulae_involving_%CF%80#Efficient_infinite_series
https://en.wikipedia.org/wiki/Pi#Spigot_algorithms


How many iterations do we need? Well, lets start with...
  4/(8*k+1) - 2/(8*k + 4) - 1/(8*k+5) - 1/(8*k+6) 
  = (120 * k^2 + 151 * k + 47)/((2 * k + 1) * (4 * k + 3) * (8 * k + 1) * (8 * k + 5))
has asymptotes at -6/8, -5/8, -4/8, and -1/8
and has 2 zeros at negative k values.
For k >= 0, the value is always positive, and slope is always negative.
For large k, it is approximately 15 / (64 k^2),
so kth term in the sum is approximately...
  15 / (64 * k^2 * 16^k)
so, for k large enough, the bit of the fractional part the kth term will modify is
  bit = - log2(15 / (64 * k^2 * 16^k))
  = - log2(15 / 64) - 2*log2(k) - 4*k
and the decimal digit is
  digit = - log10(15 / (64 * k^2 * 16^k))
For k=11, digit = 16.0
For k=23, digit = 31.0
For k=80, digit = 100.8
The more useful thing is to find which final k is needed given digit...
  solve k^2 * 16^k = 15.0 / 64 * 10^digit for k
So, for 1000 digits, you get k=825.1, so 826 should be safe enough.
Note that you can safely approximate the solution for k as...
  k = digit * log2(10) / 4
though I'd still add on 10 just to be safe...
  k = digit * log2(10) / 4 + 10
Since digit * log2(10) suffers from severe roundoff error when converting
to integer for very large digit, do...
  k = (digit * log2(10) / 4 + 10) * 1.0001


GMP manual...
  https://gmplib.org/manual/Floating_002dpoint-Functions
  https://gmplib.org/manual/Float-Internals

On macOS, I used MacPorts to do...
  sudo port install gmp
which let me compile via...
  gcc -O3 piGMP.c -I/opt/local/include/ -L/opt/local/lib/ -lgmp

On Linux, I needed to do...
  sudo apt install libgmp-dev
which let me compile via...
  gcc -O3 piGMP.c -lgmp -lm

On Windows (using Cygwin with gcc-core and libgmp-devel)...
  gcc -O3 piGMP.c -lgmp -lm
I would worry a bit more about when running this on Windows because GMP uses
"unsigned long int", which is only 32-bits on 64-bit Windows.


To run...
  ./a.out digits [BASE]
where digits will be int32_t and is the fractional digits of pi in base BASE.
If BASE is not supplied, it default to base 10.
BASE may vary from 2 to 62.
Note that the arguments are to be entered in base 10.

Run time is proportional to digits^2
A million base-10 digits will take less than half an hour, with the final output being...
  ...46460422090106105779458151
Keep in mind that the final 1 may be a 0 that was rounded up.
Running extra digits shows that this final 1 is, in fact, a 1.
To check this...
  https://www.piday.org/million/
Actually, this site and my code gives a million + 1 digits since people don't count
the first 3 as a digit (we have a million digits of the fractional part of pi).


(c) 2021 Bradley Knockel

*/


#include <stdlib.h>
#include <stdio.h>

// do I need these?
#include <stdint.h>
#include <inttypes.h>

#include <string.h>

#include <math.h>   // for log()

#include <gmp.h>

#include <sys/time.h>
struct timeval tv1, tv2;



int main(int argc, char *argv[]) {

  if(argc < 2) {
    printf(" Error: argument is required\n");
    return -1;
  }

  int32_t digits = strtol(argv[1], NULL, 10);
  if(digits <= 0) {
    printf(" Error: positive digits is required\n");
    return -1;
  }

  int32_t BASE = 10;
  if(argc > 2) {
    BASE = strtol(argv[2], NULL, 10);
    if(BASE < 2 || BASE > 62) {
      printf(" Error: BASE may only vary from 2 to 62\n");
      return -1;
    }
  }



  int32_t steps = ((double)digits * log((double)BASE) / (4.0 * log(2.0)) + 10.0)*1.0001;
  int32_t k;

  /* the "+10" is for erring on the side of safety */
  mpf_set_default_prec( (mp_bitcnt_t)(digits*log((double)BASE)/log(2.0) + 10.0) );

  mpf_t n;
  mpf_init(n);  // initialized to 0

  mpf_t temp;
  mpf_init(temp);

  mpf_t temp2;
  mpf_init(temp2);

  mpf_t one;
  mpf_init(one);
  mpf_set_ui(one, 1);



  // start timer
  gettimeofday(&tv1, NULL);



  for(k = steps - 1; k >= 0; k--){    // iterate backwards to reduce roundoff error

    /* n += ( 4/(8*k+1) - 2/(8*k + 4) - 1/(8*k+5) - 1/(8*k+6) ) / 16^k  */

    mpf_div_ui(temp, one, (uint64_t)8*k+1);
    mpf_mul_2exp(temp, temp, 2);             // multiply by 4

    mpf_div_ui(temp2, one, (uint64_t)8*k+4);
    mpf_mul_2exp(temp2, temp2, 1);
    mpf_sub(temp, temp, temp2);

    mpf_div_ui(temp2, one, (uint64_t)8*k+5);
    mpf_sub(temp, temp, temp2);

    mpf_div_ui(temp2, one, (uint64_t)8*k+6);
    mpf_sub(temp, temp, temp2);

    mpf_div_2exp(temp, temp, (uint64_t)4*k);   // divide by 16^k

    mpf_add(n, n, temp);

  }

  // clean up as much as possible
  mpf_clear(temp);
  mpf_clear(temp2);
  mpf_clear(one);



  /* print number */

  // subtract 2 to make the integer part equal to 1 (so leading zeros of fractional part in base 2 print)
  mpf_sub_ui(n, n, 2);

  mp_exp_t exp[10];      // this is never used
  char* string = (char*)malloc(digits+2);
  mpf_get_str(string, exp, BASE, digits+1, n);   // string has no radix point
  string[0] = '.';       // replace the leading 1 with a radix point
  mpf_clear(n);

  // add any trailing zeros
  int32_t len = strlen(string) - 1;  // length of fractional part
  if (len < digits) {
    char padding[digits - len + 1];
    padding[0] = '\0';
    for (k = 0; k < digits - len; k++)
      strcat(padding, "0");
    strcat(string, padding);
  }

  if (BASE == 2) 
    printf("11%s\n", string);   // 3 = 0b11
  else if (BASE == 3) 
    printf("10%s\n", string);
  else
    printf("3%s\n", string);
  free(string);



  // print elapsed wall time
  gettimeofday(&tv2, NULL);
  printf("    Running took %e seconds\n",
      (double)(tv2.tv_usec - tv1.tv_usec) / 1000000.0 + (double)(tv2.tv_sec - tv1.tv_sec)); 

  return 0;
}
