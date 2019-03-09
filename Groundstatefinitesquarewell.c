#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define tol 1E-9
#define pi 3.141592653
double func(double E);
int main()
{printf("This program calculates the ground state energy for the finite square well with using muller method\n");
printf("Mass of electron=9.11E-31(kg)\n");
printf("Depth of well=10<eV> and width of well=6 angstrom\n");
printf("V0=10<eV> and a=3<angstrom>\n");
double E,V0,d,Mass;
double hbarsq=7.6199682; /*  hbarsq=(hbar)^2=1.055E-34(j) */
double alpha,beta;
d=3; /*width of well*/
V0=16; /* Depth of well=10eV*/
Mass=9.11; /*mass of electron =9.11*10E-31*/;
alpha=sqrt(2*Mass*E/hbarsq);
beta=sqrt(2*Mass*(V0-E)/hbarsq);


double x0,x1,delta,x2,x3,lastguess;
x0=0;
x1=((pi)*(pi)*hbarsq)/(8*Mass*d*d);
int iterates=0;
double fx1,fx0,fx2,error,realerror;
double a,b,c;

do{
	
fx0=func(x0);
fx1=func(x1);
x2=(x1*fx0-x0*fx1)/(fx0-fx1);

fx2=func(x2);
a=((x1-x2)*(fx0-fx2)-(x0-x2)*(fx1-fx2))/((x0-x1)*(x0-x2)*(x1-x2));
b=((x0-x2)*(x0-x2)*(fx1-fx2)-(x1-x2)*(x1-x2)*(fx0-fx2))/((x0-x1)*(x0-x2)*(x1-x2));
c=fx2;
if(b>0)
{
x3=-2*c/(b+sqrt(b*b-4*a*c))+x2;
}
else if(b<0)
{
x3=2*c/(-b+sqrt(b*b-4*a*c))+x2;
}
error=(x3-x2)/x2;
if(error<0)
{error=-error;
}

iterates++;
if(error<tol)
break;
	realerror=error;
lastguess=x3;
if(fx0*fx1<0)
{
x1=x2;
}
else if(fx1*fx0>0)
{
x0=x2;
}
}while(1);
printf("The ground state Energy=%.9f\(eV)\n",lastguess);
printf("iterats=%d\n",iterates);
printf("error=%.9f\n",realerror);


return 0;

}
double func(double E)
{double d=3;
double V0=16;
double Mass=9.11;
double alpha,beta;
double hbarsq=7.6199682;
alpha=sqrt(2*Mass*E/hbarsq);
beta=sqrt(2*Mass*(V0-E)/hbarsq);
return beta*cos(alpha*d)-alpha*sin(alpha*d);
}
