#include <stdio.h>
#include <math.h>
#include "nrmatrix.h"
#include "nrutil.h"


/*  (C) Copr. 1986-92 Numerical Recipes Software */

/** function gaussj & dgaussj */

#define SWAP(a,b) {float temp=(a);(a)=(b);(b)=temp;}

/* NOTE!!!! */
/* 30SEP1997  AJL  changed from 1->N to 0->(N-1) indexing */

int gaussj(a,n,b,m)
float **a,**b;
int n,m;
{
	int *indxc,*indxr,*ipiv;
	int i,icol=-1,irow=-1,j,k,l,ll,*ivector();
	float big,dum,pivinv;
	void nrerror(),free_ivector();

	/*indxc=ivector(1,n);*/
	indxc=ivector(0,n-1);
	/*indxr=ivector(1,n);*/
	indxr=ivector(0,n-1);
	/*ipiv=ivector(1,n);*/
	ipiv=ivector(0,n-1);
	for (j=0;j<n;j++) ipiv[j]=0;
	for (i=0;i<n;i++) {
		big=0.0;
		for (j=0;j<n;j++)
			if (ipiv[j] != 1)
				for (k=0;k<n;k++) {
					if (ipiv[k] == 0) {
						if (fabs(a[j][k]) >= big) {
							big=fabs(a[j][k]);
							irow=j;
							icol=k;
						}
					} else if (ipiv[k] > 1)
			 			return(nrerror_return("GAUSSJ: Singular Matrix-1"));
				}
		++(ipiv[icol]);
		if (irow != icol) {
			for (l=0;l<n;l++) SWAP(a[irow][l],a[icol][l])
			for (l=0;l<m;l++) SWAP(b[irow][l],b[icol][l])
		}
		indxr[i]=irow;
		indxc[i]=icol;
		if (a[icol][icol] == 0.0)
			return(nrerror_return("GAUSSJ: Singular Matrix-2"));
		pivinv=1.0/a[icol][icol];
		a[icol][icol]=1.0;
		for (l=0;l<n;l++) a[icol][l] *= pivinv;
		for (l=0;l<m;l++) b[icol][l] *= pivinv;
		for (ll=0;ll<n;ll++)
			if (ll != icol) {
				dum=a[ll][icol];
				a[ll][icol]=0.0;
				for (l=0;l<n;l++) a[ll][l] -= a[icol][l]*dum;
				for (l=0;l<m;l++) b[ll][l] -= b[icol][l]*dum;
			}
	}
	for (l=n-1;l>=0;l--) {
		if (indxr[l] != indxc[l])
			for (k=0;k<n;k++)
				SWAP(a[k][indxr[l]],a[k][indxc[l]]);
	}
	/*free_ivector(ipiv,1,n);*/
	free_ivector(ipiv,0,n-1);
	/*free_ivector(indxr,1,n);*/
	free_ivector(indxr,0,n-1);
	/*free_ivector(indxc,1,n);*/
	free_ivector(indxc,0,n-1);

	return(0);
}


int dgaussj(a,n,b,m)
double **a,**b;
int n,m;
{
	int *indxc,*indxr,*ipiv;
	int i,icol=-1,irow=-1,j,k,l,ll,*ivector();
	double big,dum,pivinv;
	void nrerror(),free_ivector();

	/*indxc=ivector(1,n);*/
	indxc=ivector(0,n-1);
	/*indxr=ivector(1,n);*/
	indxr=ivector(0,n-1);
	/*ipiv=ivector(1,n);*/
	ipiv=ivector(0,n-1);
	for (j=0;j<n;j++) ipiv[j]=0;
	for (i=0;i<n;i++) {
		big=0.0;
		for (j=0;j<n;j++)
			if (ipiv[j] != 1)
				for (k=0;k<n;k++) {
					if (ipiv[k] == 0) {
						if (fabs(a[j][k]) >= big) {
							big=fabs(a[j][k]);
							irow=j;
							icol=k;
						}
					} else if (ipiv[k] > 1)
			return(nrerror_return("GAUSSJ: Singular Matrix-1"));
				}
		++(ipiv[icol]);
		if (irow != icol) {
			for (l=0;l<n;l++) SWAP(a[irow][l],a[icol][l])
			for (l=0;l<m;l++) SWAP(b[irow][l],b[icol][l])
		}
		indxr[i]=irow;
		indxc[i]=icol;
		if (a[icol][icol] == 0.0) 
			return(nrerror_return("GAUSSJ: Singular Matrix-2"));
		pivinv=1.0/a[icol][icol];
		a[icol][icol]=1.0;
		for (l=0;l<n;l++) a[icol][l] *= pivinv;
		for (l=0;l<m;l++) b[icol][l] *= pivinv;
		for (ll=0;ll<n;ll++)
			if (ll != icol) {
				dum=a[ll][icol];
				a[ll][icol]=0.0;
				for (l=0;l<n;l++) a[ll][l] -= a[icol][l]*dum;
				for (l=0;l<m;l++) b[ll][l] -= b[icol][l]*dum;
			}
	}
	for (l=n-1;l>=0;l--) {
		if (indxr[l] != indxc[l])
			for (k=0;k<n;k++)
				SWAP(a[k][indxr[l]],a[k][indxc[l]]);
	}
	/*free_ivector(ipiv,1,n);*/
	free_ivector(ipiv,0,n-1);
	/*free_ivector(indxr,1,n);*/
	free_ivector(indxr,0,n-1);
	/*free_ivector(indxc,1,n);*/
	free_ivector(indxc,0,n-1);

	return(0);

}

#undef SWAP




/** function  to display float matrix */

void 	DisplayMatrix(char* name, float** matrix, int num_rows, int num_cols)
{
	int nrow, ncol;

	fprintf(stdout, "\n%s Matrix: %d rows X %d columns\n", 
		name, num_rows, num_cols);
	for (nrow = 0; nrow < num_rows; nrow++) {
		for (ncol = 0; ncol < num_cols; ncol++)
			fprintf(stdout, "%.1e ", matrix[nrow][ncol]);
		fprintf(stdout, "\n");
	}
	fprintf(stdout, "\n");

}




/** function  to display double matrix */

void 	DisplayDMatrix(char* name, double** matrix, int num_rows, int num_cols)
{
	int nrow, ncol;

	fprintf(stdout, "\n%s Matrix: %d rows X %d columns\n", 
		name, num_rows, num_cols);
	for (nrow = 0; nrow < num_rows; nrow++) {
		for (ncol = 0; ncol < num_cols; ncol++) {
			if (ncol == nrow)
				fprintf(stdout, "\\ ");
			fprintf(stdout, "%.1e ", matrix[nrow][ncol]);
			if (ncol == nrow)
				fprintf(stdout, "\\ ");
		}
		fprintf(stdout, "\n");
	}
	fprintf(stdout, "\n");

}



/** function  svdcmp0 */

/* NOTE!!!! */
/* 22JUN1998  AJL  changed from 1->N to 0->(N-1) indexing */

static float at,bt,ct;
#define PYTHAG(a,b) ((at=fabs(a)) > (bt=fabs(b)) ? \
(ct=bt/at,at*sqrt(1.0+ct*ct)) : (bt ? (ct=at/bt,bt*sqrt(1.0+ct*ct)): 0.0))

static float maxarg1,maxarg2;
#define MAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
	(maxarg1) : (maxarg2))
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

int svdcmp0(a,m,n,w,v)
float **a,*w,**v;
int m,n;
{
	int flag,i,its,j,jj,k,l=-1,nm=-1;
	float c,f,h,s,x,y,z;
	float anorm=0.0,g=0.0,scale=0.0;
	float *rv1,*vector();
	void nrerror(),free_vector();

	if (m < n) return(nrerror_return(
		"SVDCMP: You must augment A with extra zero rows"));
	/*rv1=vector(1,n);*/
	rv1=vector(0,n-1);
	/*for (i=1;i<=n;i++) {*/
	for (i=0;i<n;i++) {
		l=i+1;
		rv1[i]=scale*g;
		g=s=scale=0.0;
		/*if (i <= m) {*/
		if (i < m) {
			/*for (k=i;k<=m;k++)*/
			for (k=i;k<m;k++)
				scale += fabs(a[k][i]);
			if (scale) {
				/*for (k=i;k<=m;k++) {*/
				for (k=i;k<m;k++) {
					a[k][i] /= scale;
					s += a[k][i]*a[k][i];
				}
				f=a[i][i];
				g = -SIGN(sqrt(s),f);
				h=f*g-s;
				a[i][i]=f-g;
				/*if (i != n) {*/
				if (i != n-1) {
					/*for (j=l;j<=n;j++) {*/
					for (j=l;j<n;j++) {
						/*for (s=0.0,k=i;k<=m;k++)*/
						for (s=0.0,k=i;k<m;k++) 
							s += a[k][i]*a[k][j];
						f=s/h;
						/*for (k=i;k<=m;k++)*/
						for (k=i;k<m;k++)
							a[k][j] += f*a[k][i];
					}
				}
				/*for (k=i;k<=m;k++)*/
				for (k=i;k<m;k++)
					a[k][i] *= scale;
			}
		}
		w[i]=scale*g;
		g=s=scale=0.0;
		/*if (i <= m && i != n) {*/
		if (i < m && i != n-1) {
			/*for (k=l;k<=n;k++)*/
			for (k=l;k<n;k++)
				scale += fabs(a[i][k]);
			if (scale) {
				/*for (k=l;k<=n;k++) {*/
				for (k=l;k<n;k++) {
					a[i][k] /= scale;
					s += a[i][k]*a[i][k];
				}
				f=a[i][l];
				g = -SIGN(sqrt(s),f);
				h=f*g-s;
				a[i][l]=f-g;
				/*for (k=l;k<=n;k++)*/
				for (k=l;k<n;k++)
					rv1[k]=a[i][k]/h;
				/*if (i != m) {*/
				if (i != m-1) {
					/*for (j=l;j<=m;j++) {*/
					for (j=l;j<m;j++) {
						/*for (s=0.0,k=l;k<=n;k++)*/
						for (s=0.0,k=l;k<n;k++)
							s += a[j][k]*a[i][k];
						/*for (k=l;k<=n;k++)*/
						for (k=l;k<n;k++)
							a[j][k] += s*rv1[k];
					}
				}
				/*for (k=l;k<=n;k++)*/
				for (k=l;k<n;k++)
					a[i][k] *= scale;
			}
		}
		anorm=MAX(anorm,(fabs(w[i])+fabs(rv1[i])));
	}
	/*for (i=n;i>=1;i--) {*/
	for (i=n-1;i>=0;i--) {
		/*if (i < n) {*/
		if (i < n-1) {
			if (g) {
				/*for (j=l;j<=n;j++)*/
				for (j=l;j<n;j++)
					v[j][i]=(a[i][j]/a[i][l])/g;
				/*for (j=l;j<=n;j++) {*/
				for (j=l;j<n;j++) {
					/*for (s=0.0,k=l;k<=n;k++)*/
					for (s=0.0,k=l;k<n;k++)
						s += a[i][k]*v[k][j];
					/*for (k=l;k<=n;k++)*/
					for (k=l;k<n;k++)
						v[k][j] += s*v[k][i];
				}
			}
			/*for (j=l;j<=n;j++)*/
			for (j=l;j<n;j++)
				v[i][j]=v[j][i]=0.0;
		}
		v[i][i]=1.0;
		g=rv1[i];
		l=i;
	}
	/*for (i=n;i>=1;i--) {*/
	for (i=n-1;i>=0;i--) {
		l=i+1;
		g=w[i];
		/*if (i < n)*/
		if (i < n-1)
			/*for (j=l;j<=n;j++)*/
			for (j=l;j<n;j++)
				a[i][j]=0.0;
		if (g) {
			g=1.0/g;
			/*if (i != n) {*/
			if (i != n-1) {
				/*for (j=l;j<=n;j++) {*/
				for (j=l;j<n;j++) {
					/*for (s=0.0,k=l;k<=m;k++)*/
					for (s=0.0,k=l;k<m;k++)
						s += a[k][i]*a[k][j];
					f=(s/a[i][i])*g;
					/*for (k=i;k<=m;k++)*/
					for (k=i;k<m;k++)
						a[k][j] += f*a[k][i];
				}
			}
			/*for (j=i;j<=m;j++)*/
			for (j=i;j<m;j++)
				a[j][i] *= g;
		} else {
			/*for (j=i;j<=m;j++)*/
			for (j=i;j<m;j++)
				a[j][i]=0.0;
		}
		++a[i][i];
	}
	/*for (k=n;k>=1;k--) {*/
	for (k=n-1;k>=0;k--) {
		for (its=1;its<=30;its++) {
			flag=1;
			/*for (l=k;l>=1;l--) {*/
			for (l=k;l>=0;l--) {
				nm=l-1;
				if (fabs(rv1[l])+anorm == anorm) {
					flag=0;
					break;
				}
				if (fabs(w[nm])+anorm == anorm) break;
			}
			if (flag) {
				c=0.0;
				s=1.0;
				for (i=l;i<=k;i++) {
					f=s*rv1[i];
					if (fabs(f)+anorm != anorm) {
						g=w[i];
						h=PYTHAG(f,g);
						w[i]=h;
						h=1.0/h;
						c=g*h;
						s=(-f*h);
						/*for (j=1;j<=m;j++) {*/
						for (j=0;j<m;j++) {
							y=a[j][nm];
							z=a[j][i];
							a[j][nm]=y*c+z*s;
							a[j][i]=z*c-y*s;
						}
					}
				}
			}
			z=w[k];
			if (l == k) {
				if (z < 0.0) {
					w[k] = -z;
					/*for (j=1;j<=n;j++)*/
					for (j=0;j<n;j++)
						v[j][k]=(-v[j][k]);
				}
				break;
			}
			if (its == 30) return(nrerror_return(
				"No convergence in 30 SVDCMP iterations"));
			x=w[l];
			nm=k-1;
			y=w[nm];
			g=rv1[nm];
			h=rv1[k];
			f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
			g=PYTHAG(f,1.0);
			f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
			c=s=1.0;
			for (j=l;j<=nm;j++) {
				i=j+1;
				g=rv1[i];
				y=w[i];
				h=s*g;
				g=c*g;
				z=PYTHAG(f,h);
				rv1[j]=z;
				c=f/z;
				s=h/z;
				f=x*c+g*s;
				g=g*c-x*s;
				h=y*s;
				y=y*c;
				/*for (jj=1;jj<=n;jj++) {*/
				for (jj=0;jj<n;jj++) {
					x=v[jj][j];
					z=v[jj][i];
					v[jj][j]=x*c+z*s;
					v[jj][i]=z*c-x*s;
				}
				z=PYTHAG(f,h);
				w[j]=z;
				if (z) {
					z=1.0/z;
					c=f*z;
					s=h*z;
				}
				f=(c*g)+(s*y);
				x=(c*y)-(s*g);
				/*for (jj=1;jj<=m;jj++) {*/
				for (jj=0;jj<m;jj++) {
					y=a[jj][j];
					z=a[jj][i];
					a[jj][j]=y*c+z*s;
					a[jj][i]=z*c-y*s;
				}
			}
			rv1[l]=0.0;
			rv1[k]=f;
			w[k]=x;
		}
	}
	/*free_vector(rv1,1,n);*/
	free_vector(rv1,0,n-1);

	return(0);
}

#undef SIGN
#undef MAX
#undef PYTHAG




/** function  svdcmp0 */


static float at,bt,ct;
#define PYTHAG(a,b) ((at=fabs(a)) > (bt=fabs(b)) ? \
(ct=bt/at,at*sqrt(1.0+ct*ct)) : (bt ? (ct=at/bt,bt*sqrt(1.0+ct*ct)): 0.0))

static float maxarg1,maxarg2;
#define MAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
	(maxarg1) : (maxarg2))
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

int svdcmp(a,m,n,w,v)
float **a,*w,**v;
int m,n;
{
	int flag,i,its,j,jj,k,l=-1,nm=-1;
	float c,f,h,s,x,y,z;
	float anorm=0.0,g=0.0,scale=0.0;
	float *rv1,*vector();
	void nrerror(),free_vector();

	if (m < n) return(nrerror_return(
		"SVDCMP: You must augment A with extra zero rows"));
	rv1=vector(1,n);
	for (i=1;i<=n;i++) {
		l=i+1;
		rv1[i]=scale*g;
		g=s=scale=0.0;
		if (i <= m) {
			for (k=i;k<=m;k++) scale += fabs(a[k][i]);
			if (scale) {
				for (k=i;k<=m;k++) {
					a[k][i] /= scale;
					s += a[k][i]*a[k][i];
				}
				f=a[i][i];
				g = -SIGN(sqrt(s),f);
				h=f*g-s;
				a[i][i]=f-g;
				if (i != n) {
					for (j=l;j<=n;j++) {
						for (s=0.0,k=i;k<=m;k++) s += a[k][i]*a[k][j];
						f=s/h;
						for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
					}
				}
				for (k=i;k<=m;k++) a[k][i] *= scale;
			}
		}
		w[i]=scale*g;
		g=s=scale=0.0;
		if (i <= m && i != n) {
			for (k=l;k<=n;k++) scale += fabs(a[i][k]);
			if (scale) {
				for (k=l;k<=n;k++) {
					a[i][k] /= scale;
					s += a[i][k]*a[i][k];
				}
				f=a[i][l];
				g = -SIGN(sqrt(s),f);
				h=f*g-s;
				a[i][l]=f-g;
				for (k=l;k<=n;k++) rv1[k]=a[i][k]/h;
				if (i != m) {
					for (j=l;j<=m;j++) {
						for (s=0.0,k=l;k<=n;k++) s += a[j][k]*a[i][k];
						for (k=l;k<=n;k++) a[j][k] += s*rv1[k];
					}
				}
				for (k=l;k<=n;k++) a[i][k] *= scale;
			}
		}
		anorm=MAX(anorm,(fabs(w[i])+fabs(rv1[i])));
	}
	for (i=n;i>=1;i--) {
		if (i < n) {
			if (g) {
				for (j=l;j<=n;j++)
					v[j][i]=(a[i][j]/a[i][l])/g;
				for (j=l;j<=n;j++) {
					for (s=0.0,k=l;k<=n;k++) s += a[i][k]*v[k][j];
					for (k=l;k<=n;k++) v[k][j] += s*v[k][i];
				}
			}
			for (j=l;j<=n;j++) v[i][j]=v[j][i]=0.0;
		}
		v[i][i]=1.0;
		g=rv1[i];
		l=i;
	}
	for (i=n;i>=1;i--) {
		l=i+1;
		g=w[i];
		if (i < n)
			for (j=l;j<=n;j++) a[i][j]=0.0;
		if (g) {
			g=1.0/g;
			if (i != n) {
				for (j=l;j<=n;j++) {
					for (s=0.0,k=l;k<=m;k++) s += a[k][i]*a[k][j];
					f=(s/a[i][i])*g;
					for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
				}
			}
			for (j=i;j<=m;j++) a[j][i] *= g;
		} else {
			for (j=i;j<=m;j++) a[j][i]=0.0;
		}
		++a[i][i];
	}
	for (k=n;k>=1;k--) {
		for (its=1;its<=30;its++) {
			flag=1;
			for (l=k;l>=1;l--) {
				nm=l-1;
				if (fabs(rv1[l])+anorm == anorm) {
					flag=0;
					break;
				}
				if (fabs(w[nm])+anorm == anorm) break;
			}
			if (flag) {
				c=0.0;
				s=1.0;
				for (i=l;i<=k;i++) {
					f=s*rv1[i];
					if (fabs(f)+anorm != anorm) {
						g=w[i];
						h=PYTHAG(f,g);
						w[i]=h;
						h=1.0/h;
						c=g*h;
						s=(-f*h);
						for (j=1;j<=m;j++) {
							y=a[j][nm];
							z=a[j][i];
							a[j][nm]=y*c+z*s;
							a[j][i]=z*c-y*s;
						}
					}
				}
			}
			z=w[k];
			if (l == k) {
				if (z < 0.0) {
					w[k] = -z;
					for (j=1;j<=n;j++) v[j][k]=(-v[j][k]);
				}
				break;
			}
			if (its == 30) return(nrerror_return(
				"No convergence in 30 SVDCMP iterations"));
			x=w[l];
			nm=k-1;
			y=w[nm];
			g=rv1[nm];
			h=rv1[k];
			f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
			g=PYTHAG(f,1.0);
			f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
			c=s=1.0;
			for (j=l;j<=nm;j++) {
				i=j+1;
				g=rv1[i];
				y=w[i];
				h=s*g;
				g=c*g;
				z=PYTHAG(f,h);
				rv1[j]=z;
				c=f/z;
				s=h/z;
				f=x*c+g*s;
				g=g*c-x*s;
				h=y*s;
				y=y*c;
				for (jj=1;jj<=n;jj++) {
					x=v[jj][j];
					z=v[jj][i];
					v[jj][j]=x*c+z*s;
					v[jj][i]=z*c-x*s;
				}
				z=PYTHAG(f,h);
				w[j]=z;
				if (z) {
					z=1.0/z;
					c=f*z;
					s=h*z;
				}
				f=(c*g)+(s*y);
				x=(c*y)-(s*g);
				for (jj=1;jj<=m;jj++) {
					y=a[jj][j];
					z=a[jj][i];
					a[jj][j]=y*c+z*s;
					a[jj][i]=z*c-y*s;
				}
			}
			rv1[l]=0.0;
			rv1[k]=f;
			w[k]=x;
		}
	}
	free_vector(rv1,1,n);

	return(0);
}

#undef SIGN
#undef MAX
#undef PYTHAG
