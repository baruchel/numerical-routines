/* to be compiled with    gcc -lpari test.c        */

#include <pari/pari.h>
#include <time.h>

/*
 * Beginning of the function is hacked from div_ser in src/basemath/gen1.c
 */
static GEN
inv_ser2(GEN y)
{
  long i, j, n, l = -valp(y), ly = lg(y);
  GEN y_lead, p1, p2, z;
  long vy = varn(y);
  long k = 1;
  while(k<ly-2) k <<= 1;
  k >>= 1;

  if (!signe(y)) pari_err_INV("inv_ser2", y);

  pari_sp av_base = avma;

  y_lead = gel(y,2);
  if (gequal0(y_lead)) /* normalize denominator if leading term is 0 */
  {
    pari_warn(warner,"normalizing a series with 0 leading term");
    for (l--, ly--,y++; ly > 2; l--, ly--, y++)
    {
      y_lead = gel(y,2);
      if (!gequal0(y_lead)) break;
    }
    if (ly <= 2) pari_err_INV("inv_ser2", y);
  }

  z = cgetg(ly,t_SER);
  z[1] = evalvalp(l) | evalvarn(vy) | evalsigne(1);
  gel(z,2) = gdiv(gen_1, y_lead);
  
  p2 = cgetg(ly, t_VECSMALL);
  for (i=3; i<ly; i++)
  {
    p1 = gel(y,i);
    if (isrationalzero(p1)) p1 = NULL;
    gel(p2,i) = p1;
  }

  for (n=1; n < ly-2;) {
    /* use uncomputed coefficients in z as a temporary "buffer" */
    pari_sp av0 = avma;
    for (i=2; i<n+2; i++) {
      p1 = gen_0;
      for (j=2, l=n+i; j<n+2; j++, l--) {
        if (gel(p2,l)) p1 = gsub(p1, gmul(gel(z,j), gel(p2, l)));
      }
      gel(z, n+i) = p1;
    }
    pari_sp av1 = avma;
    for (i=n+1; i>1; i--) {
      p1 = gen_0;
      for (j=2, l=n+i; j<i+1; j++, l--) {
        p1 = gadd(p1, gmul(gel(z,j), gel(z,l)));
      }
      gel(z, n+i) = p1;
    }
    stackdummy(av0, av1);
    /* compute an additional term if needed */
    n <<= 1; k >>= 1;
    if(k&(ly-2)) {
      i = (++n)+1; p1 = gen_0;
      for (j=2, l=i; j<i; j++, l--)
        if (gel(p2,l)) p1 = gsub(p1, gmul(gel(z,j), gel(p2,l)));
      gel(z,i) = gdiv(p1, y_lead);
    }
  }
  return gerepileupto(av_base, normalize(z));
}

int main(void) {
    int i, j, n;
    clock_t start;
    static GEN t;
    static GEN a;

    pari_init(50000000,1000);

    t = gp_read_str("2/3+2*x+3*x^2+4*x^3 + O(x^8)");
    pari_printf("\nReciprocal of %Ps\n", t);

    a = ginv(t); pari_printf("(ginv) %Ps\n", a); cgiv(a);
    a = inv_ser2(t); pari_printf("(test) %Ps\n", a); cgiv(a);


    t = gp_read_str("2/3+2*x+3*x^2+4*x^3 + O(x^7)");
    pari_printf("\nReciprocal of %Ps\n", t);

    a = ginv(t); pari_printf("(ginv) %Ps\n", a); cgiv(a);
    a = inv_ser2(t); pari_printf("(test) %Ps\n", a); cgiv(a);


    /* Benchmark (1) */
    t = gp_read_str("2/3+2*x+3*x^2+4*x^3");

    pari_sp av = avma;
    for(n=8; n<=32; n<<=1) {
      a = gerepileupto(av, gadd(t, ggrando(pol_x(varn(t)), n)));
      pari_printf("\nTime comparison for %Ps\n", a);

      printf("Using ginv from libpari...");
      start = clock(); for (i=0;i<10000;i++) { cgiv(ginv(a)); }
      printf(" %f millisec.\n", ((double) (clock()-start))*1000/CLOCKS_PER_SEC);

      printf("Using inv_ser2...         ");
      start = clock(); for (i=0;i<10000;i++) { cgiv(inv_ser2(a)); }
      printf(" %f millisec.\n", ((double) (clock()-start))*1000/CLOCKS_PER_SEC);
    }

    /* Benchmark (2) */

    for(n=24; n<=36; n += 12) {
      a = gerepileupto(av, gadd(t, ggrando(pol_x(varn(t)), n)));
      pari_printf("\nTime comparison for %Ps\n", a);

      printf("Using ginv from libpari...");
      start = clock(); for (i=0;i<10000;i++) { cgiv(ginv(a)); }
      printf(" %f millisec.\n", ((double) (clock()-start))*1000/CLOCKS_PER_SEC);

      printf("Using inv_ser2...         ");
      start = clock(); for (i=0;i<10000;i++) { cgiv(inv_ser2(a)); }
      printf(" %f millisec.\n", ((double) (clock()-start))*1000/CLOCKS_PER_SEC);
    }
}
