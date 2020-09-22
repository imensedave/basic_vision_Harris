
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <set>

#include <stdio.h>

#include "Imgx.h"

int round2(float x);
  
extern void BKS(float **a, int n, int *indx, float *b);

void rescaleBW_fast( Matrix<short> &im, float scal, Matrix<short> &im2) 
{
  int nr= im.rows(), nc = im.cols();
  float Nr2 = nr*scal+0.9999, Nc2 = nc*scal+0.9999;
  int nr2 = (int) Nr2, nc2 = (int) Nc2; 
  
  im2.reshape(nr2,nc2);
  im2.clear();
  
  int nr2x = nr2-1, nc2x = nc2-1;
  float R, C, xr, xc;
  int r0, c0;
  
  if( nr2 > nr ){
    const float dr = ((float) nr)/((float) nr2), dc = ((float) nc)/((float) nc2);
    for( int r=0; r<nr2x; r++){ 
      R = r * dr;
      r0 = (int) R;
      xr = R-r0;
      for( int c=0; c<nc2x; c++){ 
	C = ((float) c)* dc;
	c0 = (int) C;
	xc = C-c0;
	im2[r][c] = im[r0][c0]* (1-xr)*(1-xc) + im[r0+1][c0]* (xr)*(1-xc) +
	  im[r0][c0+1]* (1-xr)*(xc) +im[r0+1][c0+1]* (xr)*(xc);
      }
    }
    
    for( int r=0; r<nr2; r++) { 
      R = dr * r;
      //r0 = (int) (round2(R));
      r0 = (int)R;
      im2[r][nc2x] = im[r0][nc-1];
    }
    
    for( int c=0; c<nc2; c++) { 
      C = ((float) c)* dc;
      c0 = (int)C;
      im2[nr2x][c] = im[nr-1][c0];
    }
    
  } //if( scal > 1 ){
  else{
    
    Matrix<float> A(nr2,nc2);
    A.clear();
    
    for( int r=0; r<nr; r++) { 
      R = ((float) r)*scal;
      r0 = (int)R;
      if( r0 >= nr2) r0 = nr2x;
      xr = 1.0 - fabs(R-r0)*0.5;  // normalised distance from center r0
      
      for( int c=0; c<nc; c++) {
	C = ((float) c)*scal;
	c0 = (int)C;
	if( c0 >= nc2) c0 = nc2x;
	xc = 1.0 - fabs(C-c0)*0.5;  // normalised distance from center r0
	//**/if (debug && r==0) printf("c=%d: [%d][%d] += %d*%.2f*%.2f  ", c, r0, c0, im[r][c], xr, xc); 
	im2[r0][c0] += im[r][c]*xr*xc;
	A[r0][c0] += xr*xc;
      }
		
    } //for( int r=0; r<nr; r++) { 
    
    for( int r=0; r<nr2; r++) {
      for( int c=0; c<nc2; c++) {
	im2[r][c] = (int)((float)im2[r][c] / A[r][c]);
      }
    }
    
  } //else //if( scal > 1 ){
  
  return;

}

void rescaleBW( Matrix<int> &im, Matrix<int> &im2) 
{
  int nr= im.rows(), nc = im.cols();
  int nr2 = im2.rows(), nc2 = im2.cols();
  float scal = (float)nr2 / (float)nr;

  //im2.reshape(nr2,nc2);
  im2.clear();
  
  float dr = ((float) nr)/((float) nr2), dc = ((float) nc)/((float) nc2);
  int nr2x = nr2-1, nc2x = nc2-1;
  float R,C, xr, xc;
  int r0, c0;
  
  if( scal > 1 ){
		for( int r=0; r<nr2x; r++){ 
		  R = ((float) r)* dr;
		  r0 = (int) R;
		  xr = R-r0;
		  for( int c=0; c<nc2x; c++){ 
			C = ((float) c)* dc;
			c0 = (int) C;
			xc = C-c0;
			im2[r][c] = im[r0][c0]* (1-xr)*(1-xc) + im[r0+1][c0]* (xr)*(1-xc) +
						im[r0][c0+1]* (1-xr)*(xc) +im[r0+1][c0+1]* (xr)*(xc);
		  }
		}
		
		for( int r=0; r<nr2; r++) { 
		  R = ((float) r)* dr;
		  r0 = (int) (round2(R));
		  im2[r][nc2x] = im[r0][nc-1];
		}

		for( int c=0; c<nc2; c++) { 
		  C = ((float) c)* dc;
		  c0 = (int) (round2(C));
		  im2[nr2x][c] = im[nr-1][c0];
		}

  } //if( scal > 1 ){
  else{
		
		Matrix<float> A(nr2,nc2);
		A.clear();
		
		for( int r=0; r<nr; r++) { 
		  R = ((float) r)*scal;
		  r0 = (int) round2( R);
		  if( r0 >= nr2) r0 = nr2x;
		  xr = 1.0 - fabs(R-r0)*0.5;  // normalised distance from center r0
		  
		  for( int c=0; c<nc; c++) {
			C = ((float) c)*scal;
			c0 = (int) ( round2( C));
			if( c0 >= nc2) c0 = nc2x;
			xc = 1.0 - fabs(C-c0)*0.5;  // normalised distance from center r0
			//**/if (debug && r==0) printf("c=%d: [%d][%d] += %d*%.2f*%.2f  ", c, r0, c0, im[r][c], xr, xc); 
			im2[r0][c0] += im[r][c]*xr*xc;
			A[r0][c0] += xr*xc;
		  }
		
		} //for( int r=0; r<nr; r++) { 
		
		
		for( int r=0; r<nr2; r++) { 
		  for( int c=0; c<nc2; c++) {
			
			im2[r][c] = (int)round2((float)im2[r][c] / A[r][c]);
		  }
		}
		
  } //else //if( scal > 1 ){

  return;

}


/* bilinear reascaling 
   performs a little smothing for scal < 1.
*/
int round2(float x) {
	int top = ceil(x);
	int bottom = floor(x);
	if (fabs(top - x) < fabs(x - bottom))
		return top;
	return bottom;
}
void rescale( Matrix<float> &im, Matrix<float> &im2, float scal) 
{
  int nr= im.rows(), nc = im.cols();
  int nr2, nc2;

  nr2 = (int) (round2( ((float) nr)*scal));
  nc2 = (int) (round2( ((float) nc)*scal));
  im2.reshape(nr2,nc2);
  im2.clear();
  
  float dr = ((float) nr)/((float) nr2), dc = ((float) nc)/((float) nc2);
  int nr2x = nr2-1, nc2x = nc2-1;
  float R,C, xr, xc;
  int r0, c0;
  
  if( scal > 1 ){
    for( int r=0; r<nr2x; r++){ 
      R = ((float) r)* dr;
      r0 = (int) R;
      xr = R-r0;
      for( int c=0; c<nc2x; c++){ 
	C = ((float) c)* dc;
	c0 = (int) C;
	xc = C-c0;
	im2[r][c] = im[r0][c0]* (1-xr)*(1-xc) + im[r0+1][c0]* (xr)*(1-xc) +
	            im[r0][c0+1]* (1-xr)*(xc) +im[r0+1][c0+1]* (xr)*(xc);
      }
    }
    
    for( int r=0; r<nr2; r++){ 
      R = ((float) r)* dr;
      r0 = (int) (round2(R));
      im2[r][nc2x] = im[r0][nc-1];
    }

    for( int c=0; c<nc2; c++){ 
      C = ((float) c)* dc;
      c0 = (int) (round2(R));
      im2[nr2x][c] = im[nr-1][c0];
    }
  }
  else{
      Matrix<float> A(nr2,nc2);
    A.clear();
    for( int r=0; r<nr; r++){ 
      R = ((float) r)*scal;
      r0 = (int) round2( R);
      if( r0 >= nr2) r0 = nr2x;
      xr = 1.0 - fabs(R-r0)*0.5;  // normalised distance from center r0
      for( int c=0; c<nc; c++){
	C = ((float) c)*scal;
	c0 = (int) ( round2( C));
	if( c0 >= nc2) c0 = nc2x;
	xc = 1.0 - fabs(C-c0)*0.5;  // normalised distance from center r0
	im2[r0][c0] += im[r][c]*xr*xc;
	A[r0][c0] += xr*xc;
      }
    }
    
    for( int r=0; r<nr2; r++){ 
      for( int c=0; c<nc2; c++){
	im2[r][c] /= A[r][c];
      }
    }
    


  }
  return;

}


// bad for it temporarily allocates more memory
void  smoth(Matrix<int> &im, Matrix<int> &out, float sigma)
{
	
  /* make the apropriate sized gaussian mask */
	
	register int i,j,k;
	float *mask, tmpx, tmp;
	int   gauss_size, bx1,bx2,by1,by2;
	
	int nr = im.rows();
	int nc = im.cols();
	
	Matrix<int> scr_1(nr,nc);
	gauss_size = (int) ( 2.0*sigma );
	
	mask = new float [(2*gauss_size + 1)];
	mask += gauss_size;
	
	if ( sigma < 0.5 ) /* treat as zero: no smoothing */
	{
		mask[0] = 1.0;
		for ( i = 1; i <= gauss_size; i++ )
			mask[i] = mask[-i] = 0.0;
		
	}
	else  {
		float total;
		float sigma2 = 2*sigma*sigma;
		
		
		for ( i = 0; i <= gauss_size; i++ ){
			mask[i] = 0.0;
			tmpx = i - 0.5;
			for (j=0; j<6; j++){
				mask[i] += exp( - tmpx * tmpx/ (sigma2) ) ;
				tmpx += 0.2;
			}
			mask[i] = mask[i]/6.0;
		}
		total = mask[0];
		for ( i = 1; i <= gauss_size; i++ )
			total += 2.0*mask[i];
		
		mask[0] /= total;
		for ( i = 1; i <= gauss_size; i++ ) {
			mask[i] /= total;
			mask[-i] = mask[i];
		}
	}
	
  /* convelve the mask with the image */
	by1 = gauss_size; bx1 = gauss_size;
	by2 = nr-gauss_size; bx2 = nc-gauss_size;
	
	int *pi, *pend, *scr;
	
	
  /* smooth in x */
	for (i=0; i<nr;i++){
		pi = im[i] + bx1;
		pend = im[i] + bx2;
		scr = scr_1[i] + bx1;
		for( ; pi<pend; pi++, scr++){
			tmp = 0.0;
			for (k= -gauss_size; k<= gauss_size; k++){
				tmp += mask[k]*( (float) *(pi+k) );
			}
			*scr = (int) tmp;
		}
	 /* fix the boundary pts */
		for (j=0; j<bx1;j++){
			tmp = 0.0;
			tmpx = 0.0;
			for (k= -gauss_size; k<= gauss_size; k++){
				if (j+k >=0 ){
					tmp += mask[k]*( (float) im[i][j+k] );
					tmpx += mask[k];
				}
			}
			scr_1[i][j] = (int) (tmp/tmpx);
		}
		for (j=bx2; j<nc;j++){
			tmp = 0.0;
			tmpx = 0.0;
			for (k= -gauss_size; k<= gauss_size; k++){
				if (j+k <nc ){
					tmp += mask[k]*( (float) im[i][j+k] );
					tmpx += mask[k];
				}
			}
			scr_1[i][j] = (int) (tmp/tmpx);
		}
	}
	
  /* smooth in y */
	for (j=0; j<nc;j++){
		for (i=by1; i<by2;i++){
			tmp = 0.0;
			for (k= -gauss_size; k<= gauss_size; k++){
				tmp += mask[k]*(  scr_1[i+k][j] );
			}
			out[i][j] = (int) tmp;
		}
	 /* fix the boundary pts */
		for (i=0; i<by1;i++){
			tmp = 0.0;
			tmpx = 0.0;
			for (k= -gauss_size; k<= gauss_size; k++){
				if (i+k >=0 ){
					tmp += mask[k]*(  scr_1[i+k][j] );
					tmpx += mask[k];
				}
			}
			out[i][j] = (int) (tmp/tmpx);
		}
		for (i=bx2; i<nr;i++){
			tmp = 0.0;
			tmpx = 0.0;
			for (k= -gauss_size; k<= gauss_size; k++){
				if (i+k <nr ){
					tmp += mask[k]*( scr_1[i+k][j] );
					tmpx += mask[k];
				}
			}
			out[i][j] = (int) (tmp/tmpx);
		}
	}
	mask -= gauss_size;
	
	delete [] mask;
	
	return;
}



extern float  **new_fl_mat(int nr, int nc)
{
  int i;

  /* alocate ptrs to rows */
  float **mat = new float * [nr];

  /* whole shebang */
  float *m = new float [ nr * nc ];
  mat[0] = m;

  /* line up the pointers */
  for (i=1;i<nr;i++) mat[i] = mat[i-1] + nc;
  return mat;
}

extern double **newdmat(int nr, int nc)
{
  int i;

  /* alocate ptrs to rows */
  double **mat = new double * [nr];
  /* whole shebang */
  double *m = new double [ nr * nc ];
  mat[0] = m;
  /* line up the pointers */
  for (i=1;i<nr;i++) mat[i] = mat[i-1] + nc;
  return mat;
}

extern int  **new_int_mat(int nr, int nc)
{
  int i;

  /* alocate ptrs to rows */
  int **mat = new int * [nr];

  /* whole shebang */
  int *m = new int [ nr * nc ];

  /* line up the pointers */
  mat[0] = m;

  for (i=1;i<nr;i++) mat[i] = mat[i-1] + nc;
  return mat;
}

/*
void write_f_mat(float **mat, int nr, int nc, char * fname){
   std::ofstream out(fname, ios::out);
   
   out << nr << " " << nc << std::endl;
   
   for (int r = 0; r < nr; r++){
      out << mat[r][0];
      for (int c=1; c < nc; c++){
         out <<  " " << mat[r][c];
      }
      out <<  endl;
   }
   out.close();

   return;
};


void write_f_mat2(Matrix<float> &mat, char * fname){
   FILE *fp = fopen(fname, "w");
   
   int nr = mat.rows();
   int nc = mat.cols();
 
   fprintf(fp,"%d %d\n", nr,nc);
   
   for (int r = 0; r < nr; r++){
      fprintf(fp,"%.1f", mat[r][0]);
      for (int c=1; c < nc; c++){
         fprintf(fp," %.1f", mat[r][c]);
      }
      fprintf(fp,"\n");
   }
   fclose(fp);

   return;
};


// standard conversion from colour to monochrome
  // images assumed preallocated
*/

void C_to_BW( CxImage &cim, Matrix<short> &out )
{

  int nr=cim.rows();
  int nc=cim.cols();
  if( ( nr != (int) out.rows()) || ( nc != (int) out.cols() ) ){
    Matrix<short> df(nr, nc );
    out = df;
  }
  
  short *R, *G, *B;
  short *ptr, *pend;
  float tmp;
  ptr = out[0]; pend = ptr + nr*nc;
  R = cim.red()[0];
  G = cim.green()[0];
  B = cim.blue()[0];
  

  for (; ptr<pend; ptr++, R++, G++, B++){
    //tmp = (float) 0.34375* (*R) +  0.5* (*G) +  0.15625* (*B); 
    // tmp = (float) 0.2989* (*R) +  0.587* (*G);// +  0.114* (*B) + 0.5; // matlabs version.
    //    tmp = (float) 0.4* ( max( *R, *B)) +  0.6* (*G);// +  0.114* (*B) + 0.5; // matlabs version.
    tmp = (float) 0.35* (*R) +  0.56* (*G) +  0.114* (*B) ; // +  0.114* (*B) + 0.5; // matlabs version.
    
    
    *ptr = (int) tmp;
  }
}

void C_to_BW2( CxImage &cim, Matrix<short> &out )
{

  int nr=cim.rows();
  int nc=cim.cols();
  if( ( nr != (int) out.rows()) || ( nc != (int) out.cols() ) ){
    Matrix<short> df(nr, nc );
    out = df;
  }
  
  short *R, *G, *B;
  short *ptr, *pend;
  float tmp;
  ptr = out[0]; pend = ptr + nr*nc;
  R = cim.red()[0];
  G = cim.green()[0];
  B = cim.blue()[0];
  

  for (; ptr<pend; ptr++, R++, G++, B++){
    //tmp = (float) 0.34375* (*R) +  0.5* (*G) +  0.15625* (*B); 
    tmp = (float) (*G);
    *ptr = (int) tmp;
  }
}




#define TINY 1.0e-40;

int DLU(float **a,int n,int *indx,float *d)
{
  int i,imax=-1,j,k;
  float big,dum,sum,temp;
  float *vv;
	
  //        vv=vector(1,n);
	
  vv = new float [ n + 1];
	
  *d=1.0;
  for (i=1;i<=n;i++) {
    big=0.0;
    for (j=1;j<=n;j++)
      if ((temp=fabs(a[i][j])) > big) big=temp;
    if (big == 0.0){
      //if (debugmode>1)  std::cerr << "Singular matrix in routine LUDCMP" << std::endl;
      return(1);
    }
    vv[i]=1.0/big;
  }
  for (j=1;j<=n;j++) {
    for (i=1;i<j;i++) {
      sum=a[i][j];
      for (k=1;k<i;k++) sum -= a[i][k]*a[k][j];
      a[i][j]=sum;
    }
    big=0.0;
    for (i=j;i<=n;i++) {
      sum=a[i][j];
      for (k=1;k<j;k++)
	sum -= a[i][k]*a[k][j];
      a[i][j]=sum;
      if ( (dum=vv[i]*fabs(sum)) >= big) {
	big=dum;
	imax=i;
      }
    }
    if (j != imax) {
      for (k=1;k<=n;k++) {
	dum=a[imax][k];
	a[imax][k]=a[j][k];
	a[j][k]=dum;
      }
      *d = -(*d);
      vv[imax]=vv[j];
    }
    indx[j]=imax;
    if (a[j][j] == 0.0) a[j][j]= (float) TINY;
    if (j != n) {
      dum=1.0/(a[j][j]);
      for (i=j+1;i<=n;i++) a[i][j] *= dum;
    }
  }
  delete [] vv;
  return(0);
}


int inv(float **a, float **y,int dim)
{
  float d, *col;
  int i,j,*indx;
  int NNN;
  //  int *ivector();
  int sflag = 0;
  NNN =  dim;
	
  col = new float [dim +1];
  indx = new int [dim +1];
	
  //  indx = ivector(1,NNN);
	
  sflag = DLU(a,NNN,indx,&d);
  if( sflag > 0 ) return(1);

  for (j=1;j<=NNN;j++)
    {
      for (i=1;i<=NNN;i++)  col[i] = 0.0;
      col[j]=1.0;
      BKS(a,NNN,indx,col);
      for (i=1;i<=NNN;i++) y[i][j] = col[i];
    }
	
  delete [] col;
  delete [] indx;

  return(sflag);

}



extern void BKS(float **a, int n, int *indx, float *b)
{
  int i,ii=0,ip,j;
  float sum;
	
  for (i=1;i<=n;i++) {
    ip=indx[i];
    sum=b[ip];
    b[ip]=b[i];
    if (ii)
      for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
    else if (sum) ii=i;
    b[i]=sum;
  }
  for (i=n;i>=1;i--) {
    sum=b[i];
    for (j=i+1;j<=n;j++) sum -= a[i][j]*b[j];
    b[i]=sum/a[i][i];
  }
}



extern void print_mat( Matrix <float> m)
{

  int nr = m.rows();
  int nc = m.cols();

  for (int i=0; i<nr; i++){
    for (int j=0; j<nc; j++){
      fprintf( stdout,"%.3f ", m[i][j]);
    }
    fprintf( stdout,"\n");
  }
  fprintf( stdout,"\n");
}




int read_pgm( Matrix<short> &im, char *fname)
{

  FILE *fp = fopen(fname, "r");
  char strang[256];
  int nr,nc, goof=0;
  unsigned char *tmp;


  if (fp == NULL){
    printf("Error: Can't open file: %s\n",fname);
    return( 0 );
  }
  else{
    fscanf(fp,"%s", strang);
    fscanf(fp,"%d %d", &nc, &nr);

    im.reshape(nr,nc);
    tmp = new unsigned char [nr*nc];
    goof = (int) fread(tmp,sizeof(unsigned char ),nr*nc,fp);
    short *p = im[0], *pend = p+nr*nc;
    unsigned char *pu = tmp;
    for( ; p<pend; p++, pu++){
      *p = (short) *pu;
    }

    delete [] tmp;
  };

  return(goof);
}



void smo4th( Matrix<short> &im, Matrix<short> &out){

  int nr = im.rows(), nc = im.cols();
  int br = nr-1, bc = nc-1;

  short *ps, *pend, *po;
  out = im;

  for(int r=1; r<br; r++){
    ps = im[r];    po = out[r];
    pend = ps+bc;
    ps++, po++;
    for( ; ps <pend; ps++, po++){
      *ps = (*po) <<2;
      *ps += *(po-1) + *(po+1) + *(po+nc) + *(po-nc);
      *ps = (*ps) >>3;
    }
  }
  
  
  return;
};
