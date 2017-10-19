#include <math.h>
#include <stdlib.h>
#include <stdio.h>
/* 

   make random sources and receivers and save pairs in paths.txt file

   mode

   source               receiver

   random location modes

   1     uniform              uniform
   2     uniform              gaussian center
   3     Gaussian center      around sides
   4     sides                around sides

   deterministic location mode

   -1: read in sources.txt and receivers.txt and try to compute ndata
   data


*/
struct loc{
  float x,y;
  int use;
};

struct model{
  struct loc *source;
  struct loc *receiver;
  int nsource, nreceiver;
  float xtot,ytot,distmin,xcen,ycen;
  int ndata,ipick;
  FILE *os,*or,*op;  
  long iseed;
  int mode;
  
};

/* f routines */

float ran2(long *);
float gasdev(long *);


/* C */
void make_new_source(struct model *);
void make_new_receiver(struct model *);
float dfunc(float, float, float, float);
void init_model(struct model *);
void calc_geo(struct model *);



int main(int argc, char **argv)
{
  struct model m[1];
  int ndata,i,j,nuse,this_source;
  float dist;

  
  init_model(m);
 
  fprintf(stdout,"horizontal length of gridded region?\n");
  fscanf(stdin,"%f",&m->xtot);

  fprintf(stdout,"vertical length of gridded region?\n");
  fscanf(stdin,"%f",&m->ytot);

  fprintf(stdout,"how many data to build?\n");
  fscanf(stdin,"%i",&m->ndata);
  
  fprintf(stdout,"minimum acceptable epicentral distance?\n");
  fscanf(stdin,"%f",&m->distmin);
  
  fprintf(stdout,"station distribution mode?\n");
  fscanf(stdin,"%i",&m->mode);

  fprintf(stdout,"events per station?\n");
  fscanf(stdin,"%i",&m->ipick);
  /*  */
  calc_geo(m);
  if(m->mode > 0){		/* make new */
    m->or=fopen("receivers.txt","w");
    m->os=fopen("sources.txt","w");
  }else{			/* read old */
    m->or=fopen("receivers.txt","r");
    m->os=fopen("sources.txt","r");
  }
  m->op=fopen("paths.txt","w");

  fprintf(stderr,"%s: using xt/yt %g/%g, ndata %i , distmin %g, ipick %i, mode %i\n",
	  argv[0],m->xtot,m->ytot,m->ndata,m->distmin,m->ipick,m->mode);

  if(m->mode > 0){
    /* random mode, make new sources and receivers */
    ndata = 0;
    while(ndata < m->ndata){
      
      
      /* add a new source */
      make_new_source(m);
      this_source = m->nsource-1;

      for(j=0;j<m->nreceiver;j++)
	m->receiver[j].use =0;
      /* how many stations are available? */
      for(j=0,nuse=0;j<m->nreceiver;j++){
	if(dfunc(m->source[this_source].x,m->receiver[j].x,
		 m->source[this_source].y,m->receiver[j].y) > m->distmin){
	  m->receiver[j].use = 1;
	  nuse++;
	}
      }
      if(nuse > m->ipick){
	/* pick random stations */
	for(j=0;j<m->nreceiver;j++)
	  m->receiver[j].use =0;
	nuse = 0;
	while(nuse < m->ipick){
	  i = (int)(ran2(&m->iseed)*m->nreceiver);
	  if(!m->receiver[i].use){
	    m->receiver[i].use = 1;
	    nuse++;
	  }
	}
      }else if(nuse < m->ipick){
	while(nuse < m->ipick){
	  dist = 0;
	  while(dist < m->distmin){
	    /* make new receiver */
	    make_new_receiver(m);
	    j = m->nreceiver - 1;
	    dist = dfunc(m->source[this_source].x,m->receiver[j].x,
			 m->source[this_source].y,m->receiver[j].y);
	  }
	  m->receiver[j].use = 1;
	  nuse++;
	}
      }
      for(i=0;(i<m->nreceiver)&&(ndata < m->ndata);i++)
	if(m->receiver[i].use){	/* path from source to receiver */
	  fprintf(m->op,"%12.6e %12.6e %12.6e %12.6e \n",
		  m->source[this_source].x,m->source[this_source].y,
		  m->receiver[i].x,m->receiver[i].y);
	  ndata++;
	}
    } /*  */
    fprintf(stderr,"%s: done. %i data (paths) from %i new sources %i new receivers\n",
	    argv[0],m->ndata,m->nsource,m->nreceiver);

  }else{
    /* sources and receivers precomputed */
    make_new_source(m);
    make_new_receiver(m);
    fprintf(stderr,"%s: read %i sources and %i receivers\n",argv[0],m->nsource,m->nreceiver);
    this_source = 0;
    ndata = 0;
    while((ndata < m->ndata)&&(this_source < m->nsource)){

      for(j=0;j<m->nreceiver;j++) /* reset */
	m->receiver[j].use =0;
      /* check which receivers we can use */
      for(j=0,nuse=0;j<m->nreceiver;j++){
	if(dfunc(m->source[this_source].x,m->receiver[j].x,
		 m->source[this_source].y,m->receiver[j].y) > m->distmin){
	  m->receiver[j].use = 1;
	  nuse++;
	}
      }
      if(nuse > m->ipick){	/* too many */
	/* pick random receivers */
	for(j=0;j<m->nreceiver;j++)
	  m->receiver[j].use =0;
	nuse = 0;
	while(nuse < m->ipick){
	  i = (int)(ran2(&m->iseed)*m->nreceiver);
	  if(!m->receiver[i].use){
	    m->receiver[i].use = 1;
	    nuse++;
	  }
	}
      }else if(nuse < m->ipick){
	fprintf(stderr,"%s: WARNING: could only find %i out of picks for source %i, add more receivers\n",
		argv[0],nuse,m->ipick);
      }
      for(i=0;(i<m->nreceiver)&&(ndata < m->ndata);i++)
	if(m->receiver[i].use){	/* path from source to receiver */
	  fprintf(m->op,"%12.6e %12.6e %12.6e %12.6e \n",
		  m->source[this_source].x,m->source[this_source].y,
		  m->receiver[i].x,m->receiver[i].y);

	  ndata++;
	}
      this_source++;		/* next source */
    }
    if(ndata != m->ndata)
      fprintf(stderr,"%s: WARNING: could not produce all data, %i out of %i, add more sources\n",
	      argv[0],ndata,m->ndata);
    else
      fprintf(stderr,"%s: done, written %i data\n",argv[0],m->ndata);

  }
  fclose(m->os);fclose(m->or);fclose(m->op);
  return 0;

}


void make_new_source(struct model *m)
{
  float xsran,ysran,dist;
  m->source = (struct loc *)realloc(m->source, sizeof(struct loc)*(m->nsource+1));
  if(m->mode > 0){
    xsran = -1.;
    ysran = -1.;
    while((xsran < 0)||(xsran > m->xtot)||(ysran < 0)||
	  (ysran > m->xtot)){
      switch(m->mode){
      case 3:
	// Gaussian in center
	xsran=m->xcen + m->xtot/3.*gasdev(&m->iseed);
	ysran=m->ycen + m->ytot/3.*gasdev(&m->iseed);
	break;
      case 4:
	/* around sided */
	dist = 0.;
	while(dist < m->xtot/2.5){
	  xsran=m->xcen + m->xtot/7*gasdev(&m->iseed);
	  ysran=m->ycen + m->ytot/7*gasdev(&m->iseed);
	  dist=dfunc(xsran,m->xcen,ysran,m->ycen);
	}
	break;
      default:
	xsran=m->xtot*ran2(&m->iseed);
	ysran=m->ytot*ran2(&m->iseed);
	break;
      }
    }
    m->source[m->nsource].x = xsran;
    m->source[m->nsource].y = ysran;
    fprintf(m->os,"%12.6e %12.6e\n",xsran,ysran);
    m->nsource++;
  }else{			/* read all sources */
     while(fscanf(m->os,"%f %f",&(m->source[m->nsource].x),&(m->source[m->nsource].y))==2){
      m->nsource++;
      m->source = (struct loc *)realloc(m->source, sizeof(struct loc)*(m->nsource+1));
    }
  }
}

void make_new_receiver(struct model *m)
{
  float xrran,yrran;
  float dist;
  m->receiver = (struct loc *)realloc(m->receiver, 
				      sizeof(struct loc)*(m->nreceiver+1));
  if(m->mode > 0){		/* one more randome */
    switch(m->mode){
    case 1:
      //     uniformly distributed
      xrran=m->xtot*ran2(&m->iseed);
      yrran=m->ytot*ran2(&m->iseed);
      break;
    case 2:
      //     gaussian distributed around center
      xrran = -1.;
      while((xrran < 0)||(xrran > m->xtot))
	xrran = m->xcen + m->xtot/7*gasdev(&m->iseed);
      
      yrran = -1.;
      while((yrran < 0)||(yrran > m->ytot))
	yrran=m->ycen + m->ytot/7*gasdev(&m->iseed);
      break;
    case 3:
    case 4:
      //     distributed around sides
      dist = 0.;
      xrran=-1.;
      yrran=-1.;
      while((dist < m->xtot/2.5) || (xrran < 0) || (xrran > m->xtot) || 
	    (yrran < 0) || (yrran > m->ytot)){
	xrran=m->xcen + m->xtot/7*gasdev(&m->iseed);
	yrran=m->ycen + m->ytot/7*gasdev(&m->iseed);
	dist = dfunc(xrran,m->xcen,yrran,m->ycen);
      }
      break;
    default:
      fprintf(stderr,"make receiver error, mode %i undefined\n",m->mode);
      exit(-1);
      break;
    }
    m->receiver[m->nreceiver].x = xrran;
    m->receiver[m->nreceiver].y = yrran;
    fprintf(m->or,"%12.6e %12.6e\n",xrran,yrran);
    m->nreceiver++;
  }else{			/* read all receivers */
    m->nreceiver=0;
    m->receiver = (struct loc *)realloc(m->receiver, sizeof(struct loc)*(m->nreceiver+1));
    while(fscanf(m->or,"%f %f",&(m->receiver[m->nreceiver].x),&(m->receiver[m->nreceiver].y))==2){
      m->nreceiver++;
      m->receiver = (struct loc *)realloc(m->receiver, sizeof(struct loc)*(m->nreceiver+1));
    }
  }
}

float dfunc(float x1,float x2,float y1,float y2)
{
  float t1,t2;
  t1 = x1-x2;
  t2=y1-y2;
  return hypot(t1,t2);
}
void init_model(struct model *m)
{
  m->iseed = -1;
  ran2(&m->iseed);
  
  m->nsource=m->nreceiver=0;
  m->source = m->receiver = NULL;

}
void calc_geo(struct model *m)
{
  m->xcen=m->xtot/2.;
  m->ycen=m->ytot/2.;


}
