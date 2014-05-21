/*
 * p3.c: This file is part of the IEC project.
 *
 * IEC : Inspector/Executor Compiler
 *
 * Copyright (C) 2011 the Ohio State University
 *
 * This program can be redistributed and/or modified under the terms
 * of the license specified in the LICENSE.txt file at the root of the
 * project.
 *
 * Contact: P. Sadayappan <saday@cse.ohio-state.edu>
 *
 */
/**
 * @file: p3.c
 * @author: Mahesh Ravishankar <ravishan@cse.ohio-state.edu>
 * @acknowledgements : Sandip Mazumder <mazumder.2@osu.edu>
 */
#include "stdio.h"
#include "stdlib.h"
#include "assert.h"
#include "math.h"
#include "alloc.h"
#include "ie.h"

int ncells,nfaces,nnodes,nbcfaces;
double *xc,*yc,*volcell,*xf,*yf,*areaf,*vecfx,*vecfy,*xv,*yv;
int *ia_cf,*ia_cv,*lcf,*lcv,*lfc0,*lfc1,*lfv0,*lfv1;
int *bface,*f_to_bf,*bf_to_f,*bf_to_c,*bnode,*bctype;

long int memcount;
double *wcv,**wbfv;
double *t,*tb,**sc,**sb;
double **phic_new,**phic_old,**phib,**phiv,**dphi;

double **fclink,**ftlink;
double **vol_vec; 
double **ba,**bb,**bc;



enum bdy_type{
  BOTTOM,
  LEFT,
  TOP,
  RIGHT
};

enum bdy_type* bcname;

void inverse_multiply(double y[4], double A[4][4])
{
  int i,j,k,maxpivot,m,n;
  double y_old[4];
  for( i = 0 ; i < 4 ; i++ )
    y_old[i] = y[i];
  for( i = 0 ; i < 4 ; i++ ){
    double pivot = abs(A[i][i]); maxpivot = i;
    for( j = i + 1 ; j <  4 ; j++ )
      if( abs(A[j][i]) > pivot )
	maxpivot = j;
    if( maxpivot != i ){
      double temp;
      for( j = i ; j < 4 ; j++ ){
	temp = A[i][j];
	A[i][j] = A[maxpivot][j];
	A[maxpivot][j] = temp;
      }
      temp = y_old[i];
      y_old[i] = y_old[maxpivot];
      y_old[maxpivot] = temp;
    }
    pivot = A[i][i];
    for( j = i + 1 ; j < 4 ; j++ ){
      double factor = A[j][i]/pivot;
      for( k = i + 1 ; k < 4 ; k++ )
	A[j][k] -= factor * A[i][k];
      y_old[j] -= factor * y_old[i];
      A[j][i] = 0.0;
    }
  }
  for( i = 3 ; i >= 0 ; i-- ){
    y[i] = y_old[i];
    for( j = i+1 ; j < 4 ; j++ )
      y[i] -= A[i][j] * y[j];
    y[i] /= A[i][i];
  }
}


void print_output()
{
  int i;
  FILE* outfile = fopen("orig_output.dat","w");
  for( i = 0 ; i < ncells ; i++ )
    fprintf(outfile,"%d %lf %lf %lf %lf\n",i,phic_new[i][0],phic_new[i][1],phic_new[i][2],phic_new[i][3]);
  fclose(outfile);
}


void buildcoeffs(double alpha0, double beta)
{
  double c83,c11,c7_3,c61,c1,c2,c3,c4,c5,c6,sinx,cosx;
  const double alpha1 = 3.0, alpha2 = 5.0, alpha3 = 7.0;
  int i;

  c83 = 8.00 / alpha3 + 3.00 / alpha1;
  c11 = 1.00 / alpha3 + 1.00 / alpha1;
  c7_3 = 7.00 / alpha3 - 3.00 / alpha1;
  c61 = 6.00 / alpha3 + 1.00 / alpha1;
  
  double gamma[4][5];

  gamma[0][0] = 2.00 * c83 / beta;
  gamma[0][1] = - c11 / beta;
  gamma[0][2] = 2.00 * c7_3 / beta;
  gamma[0][3] = 5.00 / alpha1 / beta;
  gamma[0][4] = -2.00 * alpha2 * beta;

  gamma[1][0] = -6.00 * c11 / beta;
  gamma[1][1] = c61 / beta;
  gamma[1][2] = -6.00 * c11 / beta;
  gamma[1][3] = -5.00 / alpha1 / beta;
  gamma[1][4] = - alpha2 * beta;

  gamma[2][0] = -2.00 * c7_3 / beta;
  gamma[2][1] = - c11 / beta;
  gamma[2][2] = 2.00 * c83 / beta;
  gamma[2][3] = 5.00 / alpha1 / beta;
  gamma[2][4] = -2.00 * alpha2 * beta;

  gamma[3][0] = 6.00 / alpha1 / beta;
  gamma[3][1] = - 1.00 / alpha1 / beta;
  gamma[3][2] = 6.00 / alpha1 / beta;
  gamma[3][3] = 5.00 / alpha1 / beta;
  gamma[3][4] = -5.00 * alpha0 * beta;

  vol_vec = malloc_2d(ncells,4);
  memcount += sizeof(double)*4*ncells + sizeof(double*)*ncells ;
  for( i = 0 ; i < ncells ; i++ ){
    vol_vec[i][0] = volcell[i] * gamma[0][4];
    vol_vec[i][1] = volcell[i] * gamma[1][4];
    vol_vec[i][2] = volcell[i] * gamma[2][4];
    vol_vec[i][3] = volcell[i] * gamma[3][4];
  }
  
  ftlink = malloc_2d(nfaces,16);
  fclink = malloc_2d(nfaces,16);
  memcount += sizeof(double)*16*nfaces*2 + sizeof(double*)*nfaces*2 ;

  ba = malloc_2d(nbcfaces,16);
  bb = malloc_2d(nbcfaces,16);
  bc = malloc_2d(nbcfaces,16);
  memcount += sizeof(double)*16*nbcfaces*3 + sizeof(double*)*nbcfaces*3 ;

  for( i = 0 ; i < nfaces ; i++ ){
    double x1,y1,x2,y2,xv1,xv2,yv1,yv2;
    if( bface[i] == 0 ){
      x1 = xc[lfc0[i]];
      y1 = yc[lfc0[i]];
      x2 = xc[lfc1[i]];
      y2 = yc[lfc1[i]];
    }
    else{
      x1 = xc[lfc0[i]];
      y1 = yc[lfc0[i]];
      x2 = xf[i];
      y2 = yf[i];
    }
    xv1 = xv[lfv0[i]];
    yv1 = yv[lfv0[i]];
    xv2 = xv[lfv1[i]];
    yv2 = yv[lfv1[i]];
    
    double delta = vecfx[i] * ( x2 - x1 ) + vecfy[i] * ( y2 - y1 );
    assert(delta >= 0.0 );
    double tgtx = ( xv2 - xv1 ) / sqrt( ( xv2 - xv1 ) * ( xv2 - xv1 ) + ( yv2 - yv1 ) * ( yv2 - yv1 ) );
    double tgty =  ( yv2 - yv1 ) / sqrt( ( xv2 - xv1 ) * ( xv2 - xv1 ) + ( yv2 - yv1 ) * ( yv2 - yv1 ) );
    double tdotl = ( x2 - x1 ) * tgtx  + ( y2 - y1 ) * tgty;
    
    c1 = areaf[i] / delta;
    c2 = 2.0 * vecfx[i] * vecfy[i] * c1;
    c3 = ( vecfx[i] * tgty + vecfy[i] * tgtx ) * areaf[i];
    c4 = ( vecfx[i] * tgty - vecfy[i] * tgtx ) * areaf[i];
    c5 = ( vecfx[i] * vecfx[i] - vecfy[i] * vecfy[i] ) * c1;
    c6 = ( vecfx[i] * tgtx - vecfy[i] * tgty ) * areaf[i];
    
    fclink[i][0] = gamma[0][0]*c1; ftlink[i][0] = - gamma[0][0] * tdotl * c1;
    fclink[i][1] = gamma[0][1]*c2; ftlink[i][1] = gamma[0][1]*(c3 - tdotl * c2);
    fclink[i][2] = 0.0; ftlink[i][2] = gamma[0][2]*c4;
    fclink[i][3] = gamma[0][3]*c2; ftlink[i][3] = gamma[0][3]*(c3 - tdotl * c2);
    
    fclink[i][4] = gamma[1][0] * c2 ;
    fclink[i][5] = gamma[1][1] * c1;
    fclink[i][6] = gamma[1][2] * c5 ;
    fclink[i][7] = gamma[1][3] * c1;

    ftlink[i][4] = gamma[1][0] * (c3 - tdotl * c2 );
    ftlink[i][5] = - gamma[1][1] * tdotl * c1;
    ftlink[i][6] = gamma[1][2] * (c6 - tdotl * c5 );
    ftlink[i][7] = - gamma[1][3] * tdotl * c1;

    fclink[i][8] = 0.0;
    fclink[i][9] = gamma[2][1] * c5;
    fclink[i][10] = gamma[2][2] * c1;
    fclink[i][11] = gamma[2][3] * c5 ;

    ftlink[i][8] = gamma[2][0] * c4 ;
    ftlink[i][9] = gamma[2][1] * (c6 - tdotl * c5);
    ftlink[i][10] = - gamma[2][1] * tdotl * c1;
    ftlink[i][11] = gamma[2][3] * (c6 - tdotl * c5);

    fclink[i][12] = gamma[3][0] * c2;
    fclink[i][13] = gamma[3][1] * c1;
    fclink[i][14] = gamma[3][2] * c5 ;
    fclink[i][15] = gamma[3][3] * c1;

    ftlink[i][12] = gamma[3][0] * (c3 - tdotl * c2 );
    ftlink[i][13] = - gamma[3][1] * tdotl * c1;
    ftlink[i][14] = gamma[3][2] * (c6 - tdotl * c5);
    ftlink[i][15] = - gamma[3][3] * tdotl * c1;
    

    if( bface[i] != 0 ){
      double cosz,sinx;
      int currbf = f_to_bf[i];
      switch (bcname[currbf]){
      case BOTTOM:
	cosx = 1.0;
	sinx = 0.0;
	break;
      case RIGHT:
	cosx = -1.0;
	sinx = 0.0;
	break;
      case TOP:
	cosx = 1.0;
	sinx = 0.0;
	break;
      case LEFT:
	cosx = -1.0;
	sinx = 0.0;
	break;
      default:
	assert(0);
      }
      ba[currbf][0] = -3.0 * sinx / 4.0 - 12.0 * sinx / ( 5.0 * alpha1 * beta * delta );
      ba[currbf][1] = -1.0 / 8.0 - 2.0 / ( 5.0 * alpha1 * beta * delta );
      ba[currbf][2] = -3.0 * cosx / 4.0  - 12.0 * cosx / (5.0 * alpha1 * beta * delta );
      ba[currbf][3] = 1.0  + 2.0 / ( alpha1 * beta * delta );

      bb[currbf][0] = -12.0 * sinx / ( 5.0 * alpha1 * beta * delta );
      bb[currbf][1] = -2.0 / ( 5.0 * alpha1 * beta * delta );
      bb[currbf][2] = -12.0 * cosx / (5.0 * alpha1 * beta * delta );
      bb[currbf][3] = + 2.0 / ( alpha1 * beta * delta );

      bc[currbf][0] = 12.0 * cosx / ( 5.0 * alpha1 * beta ) - 12.0 * sinx * tdotl / ( 5.0 * alpha1 * beta * delta );
      bc[currbf][1] = -2.0 * tdotl / ( 5.0 * alpha1 * beta * delta );
      bc[currbf][2] = -12.0 * sinx / ( 5.0 * alpha1 * beta ) - 12.0 * cosx * tdotl / (5.0 * alpha1 * beta * delta );
      bc[currbf][3] = 2.0 * tdotl / ( alpha1 * delta * beta );

      ba[currbf][4] = -3.0 * cosx / 2.0 - 12.0 * cosx / ( 5.0 * alpha1 * delta * beta );
      ba[currbf][5] = 0.0;
      ba[currbf][6] = 3.0 * sinx / 2.0  + 12.0 * sinx / (5.0 * alpha1 * delta * beta);
      ba[currbf][7] = 0.0;

      bb[currbf][4] = -12.0 * cosx / ( 5.0 * alpha1 * delta * beta );
      bb[currbf][5] = 0.0;
      bb[currbf][6] = +12.0 * sinx / (5.0 * alpha1 * delta * beta);
      bb[currbf][7] = 0.0;

      bc[currbf][4] = -12.0 * sinx / ( 5.0 * alpha1 * beta ) - 12.0 * cosx * tdotl / ( 5.0 * alpha1 * delta * beta );
      bc[currbf][5] = 2.0 / ( 5.0 * alpha1 * beta );
      bc[currbf][6] = -12.0 * cosx / ( 5.0 * alpha1 * beta ) + 12.0 * sinx * tdotl / (5.0 * alpha1 * delta * beta );
      bc[currbf][7] = -2.0 / ( alpha1 * beta );

      ba[currbf][8] = 3.0 * sinx + 72.0 * sinx / ( 5.0 * alpha3 * beta * delta );
      ba[currbf][9] = 1.0 / 2.0 + 12.0 / ( 5.0 * alpha3 * beta * delta );
      ba[currbf][10] = 3.0 * cosx + 72.0 * cosx / (5.0 * alpha3 * beta * delta );
      ba[currbf][11] = 1.0 ;

      bb[currbf][8] = +72.0 * sinx / ( 5.0 * alpha3 * delta * beta );
      bb[currbf][9] = +12.0 / ( 5.0 * alpha3 * delta * beta );
      bb[currbf][10] = +72.0 * cosx / (5.0 * alpha3 * delta * beta );
      bb[currbf][11] = 0.0;

      bc[currbf][8] = 48.0 * cosx / ( 5.0 * alpha3 * beta ) + 72.0 * sinx * tdotl / ( 5.0 * alpha3 * beta * delta );
      bc[currbf][9] = 12.0 * tdotl / ( 5.0 * alpha3 * delta * beta );
      bc[currbf][10] = -48.0 * sinx / ( 5.0 * alpha3 * beta ) + 72.0 * cosx * tdotl / (5.0 * alpha3 * delta * beta );
      bc[currbf][11] = 0.0;

      ba[currbf][12] = sinx / 2.0 + 8.0 * sinx / ( 5.0 * alpha3 * delta * beta );
      ba[currbf][13] = -1.0 / 4.0 - 4.0 / ( 5.0 * alpha3 * beta * delta );
      ba[currbf][14] = cosx / 2.0 + 8.0 * cosx / (5.0 * alpha3 * delta * beta );
      ba[currbf][15] = 0.0 ;

      bb[currbf][12] = 8.0 * sinx / ( 5.0 * alpha3 * delta * beta );
      bb[currbf][13] = -4.0 / ( 5.0 * alpha3 * delta * beta );
      bb[currbf][14] = 8.0 * cosx / (5.0 * alpha3 * delta * beta );
      bb[currbf][15] = 0.0;

      bc[currbf][12] = 16.0 * cosx / ( 5.0 * alpha3 * beta ) + 8.0 * sinx * tdotl / ( 5.0 * alpha3 * beta * delta );
      bc[currbf][13] = -4.0 * tdotl / ( 5.0 * alpha3 * delta  * beta );
      bc[currbf][14] = -16.0 * sinx / ( 5.0 * alpha3 * beta ) + 8.0 * cosx * tdotl / (5.0 * alpha3 * delta * beta );
      bc[currbf][15] = 0.0;
    }
  }
}


void readgrid()
{
  FILE* one_file = fopen("gridfolder/one.out","r");
  int i,j,k;
  fscanf(one_file,"%d %d %d %d",&ncells,&nfaces,&nnodes,&nbcfaces);
  printf("one.out : %d %d %d %d\n",ncells,nfaces,nnodes,nbcfaces);
  fclose(one_file);
  
  xc = (double*) malloc(ncells*sizeof(double));
  yc = (double*) malloc(ncells*sizeof(double));
  volcell = (double*) malloc(ncells*sizeof(double));
  memcount += sizeof(double)*ncells*3;

  FILE* celldata_file = fopen("gridfolder/celldata.out","r");
  for( i = 0 ; i < ncells ; i++ ){
    int temp;
    fscanf(celldata_file,"%d %lf %lf %lf",&temp,xc+i,yc+i,volcell+i);
  }
  printf("celldata.out : %lf %lf %lf\n",xc[ncells-1],yc[ncells-1],volcell[ncells-1]);
  fclose(celldata_file);
  
  xf = (double*) malloc(nfaces*sizeof(double));
  yf = (double*) malloc(nfaces*sizeof(double));
  areaf = (double*) malloc(nfaces*sizeof(double));
  vecfx = (double*) malloc(nfaces*sizeof(double));
  vecfy = (double*) malloc(nfaces*sizeof(double));
  memcount += sizeof(double)*nfaces*5;
  
  FILE* facedata_file = fopen("gridfolder/facedata.out","r");
  for( i = 0 ; i < nfaces ; i++ ){
    int temp;
    fscanf(facedata_file,"%d %lf %lf %lf %lf %lf",&temp,xf+i,yf+i,areaf+i,vecfx+i,vecfy+i);
  }
  printf("facedata.out : %lf %lf %lf %lf %lf\n",xf[nfaces-1],yf[nfaces-1],areaf[nfaces-1],vecfx[nfaces-1],vecfy[nfaces-1]);
  fclose(facedata_file);

  xv = (double*) malloc(nnodes*sizeof(double));
  yv = (double*) malloc(nnodes*sizeof(double));
  memcount += sizeof(double)*nnodes*2;
  FILE* nodes_file = fopen("gridfolder/nodesdata.out","r");
  for( i = 0 ; i < nnodes ; i++ ){
    int temp;
    fscanf(nodes_file,"%d %lf %lf",&temp,xv+i,yv+i);
  }
  printf("nodesdata.out : %lf %lf\n",xv[nnodes-1],yv[nnodes-1]);
  fclose(nodes_file); 

  ia_cf = (int*)malloc(sizeof(int)*(ncells+1));
  ia_cv = (int*)malloc(sizeof(int)*(ncells+1));
  memcount += sizeof(int)*(ncells+1)*3;
  ia_cf[0] = 0;
  ia_cv[0] = 0;
  FILE* nfcells_file = fopen("gridfolder/nfcells.out","r");
  for( i = 0 ; i < ncells ; i++ ){
    int temp;
    fscanf(nfcells_file,"%d",&temp);
    ia_cv[i+1] = ia_cv[i] + temp;
    ia_cf[i+1] = ia_cv[i+1];
  }
  printf("nfcells.out : %d %d\n",ia_cv[ncells],ia_cf[ncells]);
  fclose(nfcells_file);

  lcv = (int*) malloc(ia_cv[ncells]*sizeof(int));
  lcf= (int*) malloc(ia_cf[ncells]*sizeof(int));
  memcount += sizeof(int)*(ia_cv[ncells]+ia_cf[ncells]);
  FILE* celldata2_file = fopen("gridfolder/celldata2.out","r");
  for( i = 0 ; i < ncells ; i++ ){
    int temp;
    fscanf(celldata2_file,"%d",&temp);
    assert(temp == i+1);
    for( j = ia_cf[i] ; j < ia_cf[i+1] ; j++ ){
      int temp2;
      fscanf(celldata2_file,"%d",&temp2);
      lcf[j] = temp2 -1 ;
    }
    fscanf(celldata2_file,"%d",&temp);
    assert(temp == i+1);
    for( j = ia_cv[i] ; j < ia_cv[i+1] ; j++ ){
      int temp2;
      fscanf(celldata2_file,"%d",&temp2);
      lcv[j] = temp2 -1;
    }
  }
  fclose(celldata2_file);
  printf("celldata2.out: %d %d\n",lcv[ia_cv[ncells]-1],lcf[ia_cf[ncells]-1]);

  lfc0 = (int*)malloc(sizeof(int)*nfaces);
  lfc1 = (int*)malloc(sizeof(int)*nfaces);
  memcount += sizeof(int)*2*nfaces + sizeof(int*)*nfaces;
  lfv0 = (int*)malloc(sizeof(int)*nfaces);
  lfv1 = (int*)malloc(sizeof(int)*nfaces);
  memcount += sizeof(int)*2*nfaces + sizeof(int*)*nfaces;

  FILE* facedata2_file = fopen("gridfolder/facedata2.out","r");
  for( i = 0 ; i < nfaces ; i++ ){
    int temp,temp1,temp2;
    fscanf(facedata2_file,"%d %d %d",&temp,&temp1,&temp2);
    lfc0[i] = temp1 - 1;
    lfc1[i] = temp2 - 1;
    fscanf(facedata2_file,"%d %d %d",&temp,&temp1,&temp2);
    lfv0[i] = temp1 - 1;
    lfv1[i] = temp2 - 1;
  }
  fclose(facedata2_file);
  printf("facedat2.out: %d %d\n",lfc0[nfaces-1],lfv1[nfaces-1]);

  bface = (int*) malloc(sizeof(int)*nfaces);
  f_to_bf = (int*) malloc(sizeof(int)*nfaces);
  memcount += sizeof(int)*2*nfaces;
  FILE* facedata3_file = fopen("gridfolder/facedata3.out","r");
  for( i = 0 ; i < nfaces ; i++ ){
    int temp;
    fscanf(facedata3_file,"%d %d %d",&temp,bface+i,f_to_bf+i);
    f_to_bf[i]--;
  }
  fclose(facedata3_file);
  printf("facedata3.out: %d %d\n",bface[nfaces-1],f_to_bf[nfaces-1]);

  bcname = (enum bdy_type*) malloc(sizeof(enum bdy_type)*nbcfaces);
  bctype = (int*) malloc(sizeof(int)*nbcfaces);
  bf_to_f = (int*) malloc(sizeof(int)*nbcfaces);
  memcount += sizeof(int)*2*nbcfaces + sizeof(enum bdy_type)*nbcfaces;

  FILE* boundarydata_file = fopen("gridfolder/boundarydata.out","r");
  char btype_names[][7] = {"bottom","left","top","right"};
  for( i = 0 ; i < nbcfaces ; i++ ){
    int temp;
    char name[10];
    fscanf(boundarydata_file,"%d %d %d %s",&temp,bf_to_f+i,bctype+i,name);
    /* printf("Name : %s\n",name); */
    if( strcmp(name,btype_names[0]) == 0 )
      bcname[i] = BOTTOM;
    else if( strcmp(name,btype_names[1]) == 0 )
      bcname[i] = LEFT;
    else if( strcmp(name,btype_names[2]) == 0 )
      bcname[i] = TOP;
    else if( strcmp(name,btype_names[3]) == 0 )
      bcname[i] = RIGHT;
    else
      assert(0);

    bf_to_f[i]--;
  }
  fclose(boundarydata_file);
  printf("boundarydata.out: %d %d\n",bctype[nbcfaces-1],bf_to_f[nbcfaces-1]);
  
  bf_to_c = (int*) malloc(sizeof(int)*nbcfaces);
  memcount += sizeof(int)*nbcfaces ;

  for( i = 0 ; i < nbcfaces ; i++ )
    bf_to_c[i] = lfc0[bf_to_f[i]];
   

  bnode = (int*) malloc(sizeof(int)*nnodes);
  memcount += sizeof(int)*nnodes ;
  for( i = 0 ; i < nnodes ; i++ )
    bnode[i] = 0;
  for( i = 0 ; i < nbcfaces ; i++){
    bnode[ lfv0[bf_to_f[i]] ] = 1;
    bnode[ lfv1[bf_to_f[i]] ] = 1;
  }

  double *wv = (double*) malloc(sizeof(double)*nnodes);
  wcv = (double*) malloc(sizeof(double)*ia_cv[ncells]);
  memcount += sizeof(double)*nnodes + sizeof(double)*lcv[ncells];
  for( i = 0 ; i < nnodes ; i++ )
    wv[i] = 0.0;

  for( i = 0 ; i < ncells ; i++){
    double xcell = xc[i];
    double ycell = yc[i];
    for( j = ia_cv[i] ; j < ia_cv[i+1] ; j++ )
      if( bnode[lcv[j]] == 0 ){
	double xnode = xv[lcv[j]];
	double ynode = yv[lcv[j]];
	wcv[j] = 1.0 / sqrt( (xcell - xnode)*(xcell - xnode) + (ycell - ynode)*(ycell-ynode));
	wv[lcv[j]] += wcv[j];
      }
  }
  
  for( i = 0 ; i < ncells ; i++ )
    for( j = ia_cv[i]  ; j < ia_cv[i+1] ; j++ )
      if( bnode[lcv[j]] == 0 )
	wcv[j] /= wv[lcv[j]];

  wbfv = malloc_2d(nbcfaces,2);
  memcount += sizeof(double)*2*nbcfaces + sizeof(double*)*nbcfaces ;
  
  for( i = 0 ; i < nbcfaces ; i++ ){
    int gf = bf_to_f[i];
    double xbdy = xf[gf];
    double ybdy = yf[gf];
    int n1 = lfv0[gf];
    assert(bnode[n1] == 1 );
    double xn1 = xv[n1];
    double yn1 = yv[n1];
    double inv_dist = sqrt( (xbdy - xn1)*(xbdy-xn1) + (ybdy-yn1)*(ybdy-yn1) );
    wbfv[i][0] = 1.0 / inv_dist;
    int n2 = lfv1[gf];
    assert(bnode[n2] == 1 );
    double xn2 = xv[n2];
    double yn2 = yv[n2];
    wbfv[i][1] = 1.0 / sqrt((xbdy - xn2)*(xbdy-xn2) + (ybdy-yn2)*(ybdy-yn2) );
    wv[n1] += wbfv[i][0];
    wv[n2] += wbfv[i][1];
  }
  for( i = 0 ; i < nbcfaces ; i++ ){
    wbfv[i][0] /= wv[ lfv0[bf_to_f[i]] ];
    wbfv[i][1] /= wv[ lfv1[bf_to_f[i]] ];
  }
  free(wv);
  memcount -= sizeof(double)*nnodes;
}



int main(int argc, char** argv)
{
  int i,j,k,l;
  memcount = 0;
  readgrid();
  printf("Grid reading finished, memcount = %ld\n",memcount);
  
  double tol = 1.0e-6;
  double kappa = 0.0, beta = 0.5, omega = 1 - kappa/beta, alpha0 = 1-omega;
  double tin = 1000.0, tw = 0.0;
  double dh1 = 0.0, dh2 = 100.0;
  double sigma = 5.67e-8;

  buildcoeffs(alpha0,beta);
  printf("Built co-efficients, memcount = %ld\n",memcount);

  tb = (double*) malloc(sizeof(double)*nbcfaces);
  t = (double*) malloc(sizeof(double)*ncells);
  memcount += sizeof(double)*nbcfaces + sizeof(double)*ncells;
  printf("Memcount = %ld\n",memcount);

  for( i = 0 ; i < ncells ; i++ )
    t[i] = 0.0;
  
  for( i = 0 ; i < nbcfaces ; i++){
    if( bcname[i] == BOTTOM ){
      int face_num = bf_to_f[i];
      if( xf[face_num] >= dh1 && xf[face_num] <= dh2 )
	tb[i] = tin;
      else
	tb[i] = tw;
    }
    else
      tb[i] = tw;
  }
  /* for( i = 0 ; i < nbcfaces ; i++ ) */
  /*   printf("%lf\n",tb[i]); */
  
  phic_old = malloc_2d(ncells,4);
  phic_new = malloc_2d(ncells,4);
  sc = malloc_2d(ncells,4);
  memcount += sizeof(double)*4*ncells*3 + sizeof(double*)*3*ncells;
  printf("Memcount = %ld\n",memcount);
  for( i = 0 ; i < ncells ; i++ ){
    for( j = 0 ; j < 4  ;j++){
      phic_old[i][j] = 0.0;
      phic_new[i][j] = 0.0;
    }
  }
  
  sb = malloc_2d(nbcfaces,4);
  phib = malloc_2d(nbcfaces,4);
  memcount += sizeof(double)*4*nbcfaces*2 + sizeof(double*)*2*nbcfaces;
  printf("Memcount = %ld\n",memcount);
  for( i = 0 ; i < nbcfaces ; i++ ){
    for( j = 0 ; j < 4 ; j++ )
      phib[i][j] = 0.0;
  }

  phiv = malloc_2d(nnodes,4);
  memcount += sizeof(double)*4*nnodes + sizeof(double*)*nnodes;
  printf("Memcount = %ld\n",memcount);
  
  dphi = malloc_2d(nfaces,4);
  memcount += sizeof(double)*4*nfaces + sizeof(double*)*nfaces;
  printf("Memcount = %ld\n",memcount);

  double resid0,resid1,resid2,resid3,rtemp;

  resid0 = tol * 2; 
  resid1 = tol * 2; 
  resid2 = tol * 2; 
  resid3 = tol * 2; 

  for( i = 0 ; i < nbcfaces ; i++ ){
    sb[i][0] = 4.0 * sigma * tb[i] * tb[i] * tb[i] * tb[i];
    sb[i][1] = 0.0;
    sb[i][2] = 4.0 * sigma * tb[i] * tb[i] * tb[i] * tb[i];
    sb[i][3] = 0.0;
  }

  for( i = 0 ; i < ncells ; i++ ){
    sc[i][0] = 0.0;
    sc[i][1] = 0.0;
    sc[i][2] = 0.0;
    sc[i][3] = -20.0 * alpha0 * sigma * t[i]*t[i]*t[i]*t[i] * beta * volcell[i];
  }

  double new_phib[4],lhs_matrix[4][4];
  int currf, currbf, n1, n2, currv;
  double alpha;
  int currcell;
  int lcv_len = ia_cv[ncells];
  double xv1,yv1,yv2,xv2;
  double diag_matrix[4][4];
  double new_phi[4];
  
  int count = 0;
//#pragma orig_loop
#pragma arrays phib[nbcfaces][4],dphi[nfaces][4],phiv[nnodes][4],phic_new[ncells][4],phic_old[ncells][4],sb[nbcfaces][4],sc[ncells][4],fclink[nfaces][16],ftlink[nfaces][16],ba[nbcfaces][16],bb[nbcfaces][16],bc[nbcfaces][16],xv[nnodes][1],yv[nnodes][1],vol_vec[ncells][4], wcv[lcv_len][1], wbfv[nbcfaces][2], bf_to_f[nbcfaces][1],bf_to_c[nbcfaces][1],lfc0[nfaces][1],lfc1[nfaces][1],lfv0[nfaces][1],lfv1[nfaces][1],bnode[nnodes][1],ia_cv[ncells+1][1],lcv[lcv_len][1],ia_cf[ncells+1][1],lcf[lcv_len][1],bface[nfaces][1],f_to_bf[nfaces][1]

//#pragma inspector_begin

  while ( resid0 > tol || resid1 > tol || resid2 > tol || resid3 > tol ){

//#pragma part_loop
    for(i = 0 ; i < nbcfaces ; i++ ){
      currf = bf_to_f[i];
      for( j = 0 ; j < 4 ; j++ )
	new_phib[j] = sb[i][j];
      for( k = 0 ; k < 4 ; k++ )
	for( l = 0 ; l < 4 ; l++ ){
	  new_phib[k] += bc[i][4*k+l] * dphi[currf][l] + bb[i][4*k+l]*phic_new[bf_to_c[i]][l];
	  lhs_matrix[k][l] = ba[i][4*k+l];
	}
      inverse_multiply(new_phib,lhs_matrix);
      
      for( j = 0 ; j < 4 ; j++ )
	phib[i][j] = new_phib[j];
    }

//#pragma part_loop
    for( i = 0 ; i < nnodes ; i++ )
      for( j = 0 ; j < 4 ; j++ )
	phiv[i][j] = 0.0;
  
//#pragma part_loop
    for( i = 0 ; i < nbcfaces ; i++ ){
      currf = bf_to_f[i];
      n1 = lfv0[currf];
      n2 = lfv1[currf];
      for( j = 0 ; j < 4 ; j++ ){
	phiv[n1][j] += phib[i][j] * wbfv[i][0];
	phiv[n2][j] += phib[i][j] * wbfv[i][1];
      }
    }

//#pragma part_loop
    for( i = 0 ; i < ncells ; i++ ){
      for( j = ia_cv[i] ; j < ia_cv[i+1]  ; j++ ){
	currv = lcv[j];
	if( bnode[currv] == 0 )
	  for( k = 0 ; k < 4 ; k++ )
	    phiv[currv][k] += phic_new[i][k] * wcv[j];
      }
      for( j = 0 ; j < 4 ; j++ )
	phic_old[i][j] = phic_new[i][j];
    }
    

//#pragma part_loop
    for( i = 0 ; i < nfaces ; i++ ){
      n1 = lfv0[i];
      n2 = lfv1[i];
      xv1 = xv[n1];
      xv2 = xv[n2];
      yv1 = yv[n1];
      yv2 = yv[n2];
      for( j = 0 ; j < 4 ; j++ )
	dphi[i][j] = ( phiv[n2][j] - phiv[n1][j] ) / sqrt( ( xv2 - xv1 ) * ( xv2 - xv1 ) + ( yv2 - yv1 ) * ( yv2 - yv1 ) );
    }
    
    
//#pragma part_loop
    for( i = 0 ; i < ncells ; i++ ){
      for( j = 0 ; j < 4 ; j++ ){
	for( k = 0 ; k < 4 ; k++ )
	  diag_matrix[j][k] = 0.0;
	new_phi[j] = sc[i][j];
	diag_matrix[j][j] = vol_vec[i][j];
      }
    
      for( j = ia_cf[i] ; j < ia_cf[i+1] ; j++ ){
	currf = lcf[j];

	if( bface[currf] == 0 ){
	  if( lfc0[currf] == i ){
	    currcell = lfc1[currf];
	    alpha = 1;
	  }
	  else{
	    alpha = -1;
	    currcell = lfc0[currf];
	  }
	  
	  for( k = 0 ; k < 4 ; k++ )
	    for( l = 0 ; l < 4 ; l++ ){
	      new_phi[k] += -(fclink[currf][k*4+l] * phic_old[currcell][l] + alpha * ftlink[currf][k*4+l] * dphi[currf][l]) ;
	      diag_matrix[k][l] += -fclink[currf][k*4+l];
	    }
	}
	else{
	  currbf = f_to_bf[currf];
	  for( k = 0 ; k < 4 ; k++ )
	    for( l = 0 ; l < 4 ; l++ ){
	      new_phi[k] += -(fclink[currf][4*k+l] * phib[currbf][l] + ftlink[currf][4*k+l] * dphi[currf][l]);
	      diag_matrix[k][l] += -fclink[currf][4*k+l] ;
	    }
	}
      }
      inverse_multiply(new_phi,diag_matrix);
      
      for( j = 0 ; j < 4 ; j++ )
	phic_new[i][j] = new_phi[j];
    }
    
    resid0 = 0.0;
    resid1 = 0.0;
    resid2 = 0.0;
    resid3 = 0.0;
    
//#pragma part_loop
    for( i = 0 ; i < ncells ; i++ ){
      rtemp = phic_new[i][0] - phic_old[i][0];
      resid0 += rtemp * rtemp;
      rtemp = phic_new[i][1] - phic_old[i][1];
      resid1 += rtemp * rtemp;
      rtemp = phic_new[i][2] - phic_old[i][2];
      resid2 += rtemp * rtemp;
      rtemp = phic_new[i][3] - phic_old[i][3];
      resid3 += rtemp * rtemp;
    }
    resid0 = sqrt(resid0);
    resid1 = sqrt(resid1);
    resid2 = sqrt(resid2);
    resid3 = sqrt(resid3);

    count++;
    //printf("%d %lf %lf %lf %lf\n",count,resid0,resid1,resid2,resid3);
  }  

//#pragma inspector_end

  print_output();

  free(xc);
  free(yc);
  free(xf);
  free(yf);
  free(volcell);
  free(areaf);
  free(vecfx);
  free(vecfy);
  free(xv);
  free(yv);
  free(ia_cf);
  free(ia_cv);
  free(lcf);
  free(lcv);
  free(lfc0);
  free(lfc1);
  free(lfv0);
  free(lfv1);
  free(bface);
  free(f_to_bf);
  free(bf_to_f);
  free(bf_to_c);
  free(bnode);
  free(bctype);
  free(bcname);
  free(wcv);
  free_2d(wbfv);
  free(t);
  free(tb);
  free_2d(sc);
  free_2d(sb);
  free_2d(phic_new);
  free_2d(phic_old);
  free_2d(phib);
  free_2d(phiv);
  free_2d(dphi);
  free_2d(fclink);
  free_2d(ftlink);
  free_2d(vol_vec);
  free_2d(ba);
  free_2d(bb);
  free_2d(bc); 

  return 0;
}
