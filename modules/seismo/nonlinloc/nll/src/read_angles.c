#include <stdio.h>
#include <stdlib.h>
typedef union
{
    float fval;               /* float value (dummy) */
    unsigned short ival[2];   /* unsigned short values:
                                   ival[0] bits 0-3 = quality
                                           bits 4-15 = dip
                                   ival[1] = azimuth */
} TakeOffAngles;
TakeOffAngles SetTakeOffAngles(double , double , int );
int GetTakeOffAngles(TakeOffAngles *,
		     double *, double *, int *);

int main(void)
{
  TakeOffAngles aun;
  double azi,dip;
  unsigned short ival[2];
  int qual;

  while(fread(&aun,sizeof(TakeOffAngles), 1, stdin)==1){
    GetTakeOffAngles(&aun,&azi,&dip,&qual);
    printf("%11g %11g %i\n",azi,dip,qual);

  }
  return 0;
}


/*** function to set angle values in take-off angles union */
TakeOffAngles SetTakeOffAngles(double azim, double dip, int iqual)
{
	TakeOffAngles angles;

	angles.ival[1] = (unsigned short) (0.5 + 10.0 * azim);
	angles.ival[0] = (unsigned short) iqual
		+ (unsigned short) 16
		* (unsigned short) (0.5 + 10.0 * dip);

	return(angles);
}

/*** function to get values in take-off angles union */
int GetTakeOffAngles(TakeOffAngles *pangles,
			double *pazim, double *pdip, int *piqual)
{
	*pazim = ((double) pangles->ival[1]) / 10.0;
	*pdip = ((double) (pangles->ival[0] / (int) 16)) / 10.0;
	*piqual = (int) pangles->ival[0] % (int) 16;

	return(*piqual);
}
