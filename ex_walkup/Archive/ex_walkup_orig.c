//-----------------------------------------------
// example of vector permute operations
//-----------------------------------------------
#include <stdio.h>

#define NPTS 16
static double __attribute__((aligned(32))) x[NPTS];

int main(int argc, char * argv[])
{
   int i;
   vector4double v0, v1, v2, v3;
   vector4double p0, p1, p2, p3;
   vector4double q0, q1, q2, q3;

   for (i=0; i<NPTS; i++) x[i] = (double) i;

   v0 = vec_ld(0L, &x[0]);    // v0 has values x[0-3]
   v1 = vec_ld(0L, &x[4]);    // v1 has values x[4-7]
   v2 = vec_ld(0L, &x[8]);    // v2 has values x[8-11]
   v3 = vec_ld(0L, &x[12]);   // v3 has values x[12-15]

   p0 = vec_perm(v0, v1, vec_gpci(00415));
   p1 = vec_perm(v0, v1, vec_gpci(02637));
   p2 = vec_perm(v2, v3, vec_gpci(00415));
   p3 = vec_perm(v2, v3, vec_gpci(02637));

   q0 = vec_perm(p0, p2, vec_gpci(00145));
   q1 = vec_perm(p0, p2, vec_gpci(02367));
   q2 = vec_perm(p1, p3, vec_gpci(00145));
   q3 = vec_perm(p1, p3, vec_gpci(02367));

   printf("\n");

   printf("q0 = %4.1lf %4.1lf %4.1lf %4.1lf\n", q0[0], q0[1], q0[2], q0[3]);
   printf("q1 = %4.1lf %4.1lf %4.1lf %4.1lf\n", q1[0], q1[1], q1[2], q1[3]);
   printf("q2 = %4.1lf %4.1lf %4.1lf %4.1lf\n", q2[0], q2[1], q2[2], q2[3]);
   printf("q3 = %4.1lf %4.1lf %4.1lf %4.1lf\n", q3[0], q3[1], q3[2], q3[3]);

   return 0;
}
