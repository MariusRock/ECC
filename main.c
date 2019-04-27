#include <stdio.h>
#include <gmp.h>

#include <string.h> // for strlen

#define init_values {\
  mpz_init(k);\
  mpz_init(Gx);\
  mpz_init(Gy);\
  mpz_init(a);\
  mpz_init(b);\
  mpz_init(x);\
  mpz_init(y);\
  mpz_init(xo);\
  mpz_init(yo);\
  mpz_init(tmp_exp);\
  mpz_init(n);\
  mpz_init(r);\
}

#define clear_values {\
  mpz_clear(k);\
  mpz_clear(Gx);\
  mpz_clear(Gy);\
  mpz_clear(a);\
  mpz_clear(b);\
  mpz_clear(x);\
  mpz_clear(y);\
  mpz_clear(xo);\
  mpz_clear(yo);\
  mpz_clear(tmp_exp);\
  mpz_clear(n);\
  mpz_init(r);\
}

mpz_t x,y;         // point (x,y)
mpz_t a,b,n;       // Elliptic Curve y²=x³+ax+b mod n
mpz_t Gx,Gy;       // start point on curve
mpz_t k;           // private key

// Temporal Variables
mpz_t xo,yo;       // x_old and y_old , which are temporal copies of x and y (used in Addition and Doubling)
mpz_t tmp_exp,r;   // genereal temporal expression (try not to avoid where possible)


void Evaluate_function(mpz_t result,const mpz_t input_x){
  mpz_mul(result,a,x);
  mpz_powm_ui(tmp_exp,x,3,n);
  mpz_add(result,result,tmp_exp);
  mpz_add(result,result,b);
  mpz_mod(result,result,n);
}

void point_addition(mpz_t x,mpz_t y,mpz_t x2,mpz_t y2){

  if (mpz_cmp_ui(x, 0UL) == 0 && mpz_cmp_ui(y, 0UL) == 0){
    mpz_set(x,x2);
    mpz_set(y,y2);
    return;
  }

	if (mpz_cmp_ui(x2, 0UL) == 0 && mpz_cmp_ui(y2, 0UL) == 0)
  return;

  if (mpz_cmp(x,x2) == 0){
    mpz_set_ui(x,0UL);
    mpz_set_ui(y,0UL);
    printf("Warning, addition with same x coordinate ");
    return;
  }

  mpz_sub(xo,x,x2);
  mpz_sub(yo,y,y2);
  mpz_invert(tmp_exp,xo,n);
  mpz_mul(tmp_exp,tmp_exp,yo);  // tmp_exp = ( (y-y2) * (x-x2)⁻¹ )
  mpz_mod(tmp_exp,tmp_exp,n);
  mpz_set(xo,x);
  mpz_set(yo,y);
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
  if(mpz_get_ui(y)==0) printf("Warning, inverse doesn't exist "); // Check might be neglected
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


void scalar_multiplication(mpz_t x,mpz_t y,mpz_t k){
  mpz_set(Gx,x);
  mpz_set(Gy,y);
  char *val_str = mpz_get_str(NULL, 2, k);
  size_t len = strlen(val_str);
  for (int i = 1; i != len;i++)  {
  point_doubling(x,y); // This leads to Doubling of point infinity in first round, which is uncritical
   if (val_str[i] == '1'){
    point_addition(x,y,Gx,Gy);
  }
  }
}

int main()
{
  printf("limb size: %i bit \n",mp_bits_per_limb);  // A limb means the part of a multi-precision number that fits in a single machine word. Normally a limb is 32 or 64 bits.
  printf("gmp version: %s \n",gmp_version);

  init_values // init values and sets values to zero

  mpz_set_ui(n,17);
  mpz_set_ui(a,2);
	mpz_set_ui(b,2);
  mpz_set_ui(k,97);

  mpz_init_set_str (x, "5", 10);
  Evaluate_function(y,x); // y=f(x)=x

	printf("point G = (%lu,%lu) \n",mpz_get_ui (x),mpz_get_ui (y));


  for(int i=2;i!=25;i++){
    mpz_set_ui(k,i);
    mpz_init_set_str (x, "5", 10);
    Evaluate_function(y,x); // y=f(x)=x
  scalar_multiplication(x,y,k);
  printf("point %lu G = (%lu,%lu) \n",mpz_get_ui(k),mpz_get_ui(x),mpz_get_ui(y));
  }

  clear_values

  return 0;
}
