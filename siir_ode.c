#include <R.h>
static double parms[7];

#define bA parms[0]
#define bB parms[1]
#define fA parms[2]
#define fB parms[3]
#define gA parms[4]
#define gB parms[5]
#define pop parms[6]
  /* initializer  */
  void initmod(void (* odeparms)(int *, double *))
  {
    int N=7;
    odeparms(&N, parms);
  }
/* Derivatives and 1 output variable */
void derivs (int *neq, double *t, double *y, double *ydot,
             double *yout, int *ip)
{
  ydot[0] = -(bA*y[3]+bB*y[4])*y[0]/pop;
  ydot[1] = bA*y[3]*y[0]/pop-fA*y[1];
  ydot[2] = bB*y[4]*y[0]/pop-fB*y[2];
  ydot[3] = fA*y[1]-gA*y[3];
  ydot[4] = fB*y[2]-gB*y[4];
  ydot[5] = gA*y[3];
  ydot[6] = gB*y[4];
}
/* The Jacobian matrix */
void jac(int *neq, double *t, double *y, int *ml, int *mu,
  double *pd, int *nrowpd, double *yout, int *ip)
{
  pd[0] = -(bA*y[3]+bB*y[4])/pop;
  pd[1] = bA*y[3]/pop;
  pd[2] = bB*y[4]/pop;
  pd[3] = 0;
  pd[4] = 0;
  pd[5] = 0;
  pd[6] = 0;
  
  pd[(*nrowpd)] = 0;
  pd[(*nrowpd) + 1] = -fA;
  pd[(*nrowpd) + 2] = 0;
  pd[(*nrowpd) + 3] = fA;
  pd[(*nrowpd) + 4] = 0;
  pd[(*nrowpd) + 5] = 0;
  pd[(*nrowpd) + 6] = 0;
  
  pd[(*nrowpd)*2] = 0;
  pd[(*nrowpd)*2 + 1] = 0;
  pd[(*nrowpd)*2 + 2] = -fB;
  pd[(*nrowpd)*2 + 3] = 0;
  pd[(*nrowpd)*2 + 4] = fB;
  pd[(*nrowpd)*2 + 5] = 0;
  pd[(*nrowpd)*2 + 6] = 0;
  
  pd[(*nrowpd)*3] = -bA*y[0]/pop;
  pd[(*nrowpd)*3 + 1] = bA*y[0]/pop;
  pd[(*nrowpd)*3 + 2] = 0;
  pd[(*nrowpd)*3 + 3] = -gA;
  pd[(*nrowpd)*3 + 4] = 0;
  pd[(*nrowpd)*3 + 5] = gA;
  pd[(*nrowpd)*3 + 6] = 0;
  
  pd[(*nrowpd)*4] = -bB*y[0]/pop;
  pd[(*nrowpd)*4 + 1] = 0;
  pd[(*nrowpd)*4 + 2] = bB*y[0]/pop;
  pd[(*nrowpd)*4 + 3] = 0;
  pd[(*nrowpd)*4 + 4] = -gB;
  pd[(*nrowpd)*4 + 5] = 0;
  pd[(*nrowpd)*4 + 6] = gB;
  
  pd[(*nrowpd)*5] = 0;
  pd[(*nrowpd)*5 + 1] = 0;
  pd[(*nrowpd)*5 + 2] = 0;
  pd[(*nrowpd)*5 + 3] = 0;
  pd[(*nrowpd)*5 + 4] = 0;
  pd[(*nrowpd)*5 + 5] = 0;
  pd[(*nrowpd)*5 + 6] = 0;
  
  pd[(*nrowpd)*6] = 0;
  pd[(*nrowpd)*6 + 1] = 0;
  pd[(*nrowpd)*6 + 2] = 0;
  pd[(*nrowpd)*6 + 3] = 0;
  pd[(*nrowpd)*6 + 4] = 0;
  pd[(*nrowpd)*6 + 5] = 0;
  pd[(*nrowpd)*6 + 6] = 0;
}