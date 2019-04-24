#include <stdio.h>
#include <gmp.h>

#define init_values {\
  mpz_init(a);\
  mpz_init(b);\
  mpz_init(x);\
  mpz_init(x2);\
  mpz_init(y);\
  mpz_init(y2);\
  mpz_init(xo);\
  mpz_init(yo);\
  mpz_init(tmp_exp);\
  mpz_init(n);\
  mpz_init(r);\
}

#define clear_values {\
  mpz_clear(a);\
  mpz_clear(b);\
  mpz_clear(x);\
  mpz_clear(x2);\
  mpz_clear(y);\
  mpz_clear(y2);\
  mpz_clear(xo);\
  mpz_clear(yo);\
  mpz_clear(tmp_exp);\
  mpz_clear(n);\
  mpz_init(r);\
}


mpz_t x,x2,y,y2;   // two points (x,y)
mpz_t a,b,n;       // Elliptic Curve y²=x³+ax+b mod n

// Temporal Variables
mpz_t xo,yo;       // x_old and y_old , which are temporal copies of x and y
mpz_t tmp_exp,r;   // genereal temporal expression (try not to avoid where possible)

// Class of Functions
// 1. Functions for signed integer arithmetic, with names beginning with  mpz_
// 4. Fast low-level functions that operate on natural numbers begin with mpn_

void Evaluate_function(mpz_t result,const mpz_t input_x){
  mpz_mul(result,a,x);
  mpz_powm_ui(tmp_exp,x,3,n);
  mpz_add(result,result,tmp_exp);
  mpz_add(result,result,b);
  mpz_mod(result,result,n);
}




void point_addition(mpz_t x,mpz_t y,mpz_t x2,mpz_t y2){

  printf("(x1,y1)=(%lu,%lu) and (x2,y2)=(%lu,%lu) -> ",mpz_get_ui(x),mpz_get_ui(y),mpz_get_ui(x2),mpz_get_ui(y2));

  mpz_sub(xo,x,x2);
  mpz_sub(yo,y,y2);
  mpz_cdiv_qr(tmp_exp,r,yo,xo); // tmp_exp = (y-y2)/(x-x2)

  mpz_set(xo,x);
  mpz_set(yo,y);

  if(mpz_get_ui(r)!=0) printf("! ! ! Remainder in point addition unequal zero ! ! ! \n");
  mpz_mul(r,tmp_exp,tmp_exp);
  mpz_add(yo,x,x2);
  mpz_sub(x,r,yo);               // x_new = tmp_exp²-(x+x2)
  mpz_mod(x,x,n);

  mpz_sub(r,x2,x);
  mpz_mod(r,r,n); // spped enhancement due to calculation iwth small postive numbers



  mpz_mul(tmp_exp,tmp_exp,r);
  mpz_sub(tmp_exp,tmp_exp,y2);                // y_new = tmp_exp * (x-x_new) - y
  mpz_mod(y,tmp_exp,n);
}

void point_doubling(mpz_t x,mpz_t y){ // x=[(3x²+a)/(2y)]² - 2 x
  mpz_set(xo,x);
  mpz_set(yo,y);
  mpz_mul(tmp_exp,x,x);            // x²
  mpz_mul_ui(tmp_exp,tmp_exp,3);   //3x²
  mpz_add(tmp_exp,tmp_exp,a);      //3x²+a
  mpz_mul_2exp(y,yo,1);            // ## NOTE: y will be used as a tmp var here
  mpz_invert(y,y,n);               // ## (2y)⁻¹
  if(mpz_get_ui(y)==0) printf("Error, inverse doesn't exist \n"); // Check might be neglected
  mpz_mul(tmp_exp,tmp_exp,y);      // ##
  mpz_mod(tmp_exp,tmp_exp,n);      // s=[(3x²+a)/(2y)]
  mpz_mul(x,tmp_exp,tmp_exp);      // x=s²
  mpz_mul_2exp(y,xo,1);            // ## NOTE: y will be used as a tmp var here
  mpz_sub(x,x,y);                  // ## -2y
  mpz_mod(x,x,n);
  mpz_sub(y,xo,x);
  mpz_mul(y,y,tmp_exp);
  mpz_sub(y,y,yo);
  mpz_mod(y,y,n);                  // y_new = s * (x_old-x_new)-y_old
}

int main()
{
  printf("limb size: %i bit \n",mp_bits_per_limb);  // A limb means the part of a multi-precision number that fits in a single machine word. Normally a limb is 32 or 64 bits.
  printf("gmp version: %s \n",gmp_version);

  init_values // init values and sets values to zero

  mpz_init_set_str (a, "2", 10);
  mpz_init_set_str (n, "17", 10);
	mpz_set_ui(b,2);

  mpz_init_set_str (x, "5", 10);

  Evaluate_function(y,x); // y=f(x)=x

  mpz_set(x2,x);
  mpz_set(y2,y);


  unsigned long int y_ui,x_ui;

  x_ui = mpz_get_ui (x); // cast back to unsigned int
	y_ui = mpz_get_ui (y); // cast back to unsigned int
	printf("point= (%lu,%lu) \n",x_ui,y_ui);

  point_doubling(x,y);

  x_ui = mpz_get_ui (x); // cast back to unsigned int
  y_ui = mpz_get_ui (y); // cast back to unsigned int
	printf("point= (%lu,%lu) \n",x_ui,y_ui);

  for(int i=0;i!=5;i++)
    {
    point_addition(x,y,x2,y2);
    //point_doubling(x,y);

    x_ui = mpz_get_ui (x); // cast back to unsigned int
    y_ui = mpz_get_ui (y); // cast back to unsigned int
	   printf("point= (%lu,%lu) \n",x_ui,y_ui);
   }


  clear_values

  return 0;

}
