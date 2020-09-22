#ifndef _Imgx_h
#define _Imgx_h

#include <iostream>
#include <fstream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <set>
#include <vector>



using namespace std;

/*  copyright David Sinclair 2020.

    released under the Apache2.0 license.
    https://www.apache.org/licenses/LICENSE-2.0

*/


template <class T>
class Matrix
{
 protected:
  unsigned int   nr;
  unsigned int   nc;
  T            **bin;
  T             *bulk;
  
 public:
  /**/T** getData() { return bin; }
  
 Matrix(void) : nr(1), nc(1){
    bin      = new T* [nr];
    bulk     = new T[nr * nc];
    
    T   *row = bulk;
    for (T   **b = bin; b < bin + nr; row += nc)
      *b++ = row;
    clear();
  };
  
 Matrix(unsigned int r, unsigned int c) : nr(r), nc(c){
    bin      = new T* [nr];
    bulk     = new T[nr * nc];
    
    T   *row = bulk;
    for (T   **b = bin; b < bin + nr; row += nc)
      *b++ = row;
    clear();
  };
 Matrix(const Matrix &m) : nr(m.nr), nc(m.nc){
    bin      = new T* [nr];
    bulk     = new T[nr * nc];
    
    T   *row = bulk;
    for (T   **b = bin; b < bin + nr; row += nc)
      *b++ = row;
    
    memcpy((void *) bulk, (void *) m.bulk, sizeof(T) * nr * nc);
  };
  
  Matrix &operator=(const Matrix &m) {
    if (&m != this) {
      if (nr != m.nr || nc != m.nc) {
	nr = m.nr;
	nc = m.nc;
	delete [] bin;
	delete [] bulk;
	bin  = new T* [nr];
	bulk = new T[nr * nc];
	
	T *row = bulk;
	for (T   **b = bin; b < bin + nr; row += nc)
	  *b++ = row;
      }
      
      memcpy((void *) bulk, (void *) m.bulk, sizeof(T) * nr * nc);
    }
    return *this;
  };
  
  ~Matrix(void){
    delete [] bin;
    delete [] bulk;
  };
  
  void  clear(void){
    memset((void *) bulk, 0, sizeof(T) * nr * nc);
  };
  
  inline T            *operator[](unsigned int i) { return bin[i]; };
  inline const T      *operator[](unsigned int i) const { return bin[i]; };
  inline unsigned int  rows() const     { return nr; };
  inline unsigned int  cols() const     { return nc; };
  
  void reshape(int d1, int d2) {
    if(( nr != d1) || (nc != d2)){
      nr = d1;
      nc = d2;
      delete [] bin;
      delete [] bulk;
      bin  = new T* [nr];
      bulk = new T[nr * nc];
      
      T *row = bulk;
      for (T   **b = bin; b < bin + nr; row += nc)
	*b++ = row;
    }
  }
  
  void write_rlm_ascii(const char *fname) {
    std::ofstream out(fname, ios::out);
    
    out << nr << " " << nc <<" run length map(r,c) in columns" << std::endl;
    
    for (unsigned int i1 = 0; i1 < nc; ++i1){
      int val, count;
      val = bin[0][i1]; // indexed as m(r,c)
      count = 0;
      for (unsigned int i2 = 0; i2 < nr; ++i2){
	if( val == bin[i2][i1] ){
	  count ++;
	}
	else{
	  out<<val <<" " << count << " ";
	  val =  bin[i2][i1];
	  count = 1;
	}
	//				out.write( &value, sizeof(int));
      }
      out << val << " " << count << std::endl;
    }
    out.close();
    
    return;
  };
  void cerr_rlm_ascii() {
    cerr << nr << " " << nc <<" run length map(r,c) in columns" << std::endl;
    
    for (unsigned int i1 = 0; i1 < nc; ++i1){
      int val, count;
      val = bin[0][i1]; // indexed as m(r,c)
      count = 0;
      for (unsigned int i2 = 0; i2 < nr; ++i2){
	if( val == bin[i2][i1] ){
	  count ++;
	}
	else{
	  cerr<<val <<" " << count << " ";
	  val =  bin[i2][i1];
	  count = 1;
	}
	//				out.write( &value, sizeof(int));
      }
      cerr << val << " " << count << std::endl;
    }
    return;
  };
  
  int write_PGM(char *fname) {
    FILE *file_id;
    file_id = fopen(fname,"w");
    if (file_id == NULL){
      printf("Error: Can't open file: %s\n",fname);
      return( 0 );
    }
    fprintf(file_id, "P5\n%d %d\n255\n", nc, nr);
    unsigned char *tmpi = new unsigned char [nr*nc];
    T *pb = bulk;
    unsigned char *pi = tmpi, *pend = pi+nr*nc;
    for( ; pi<pend; pi++, pb++){
      *pi = (unsigned char) *pb;
    }
    fwrite(tmpi,sizeof(char), nr*nc, file_id);
    delete [] tmpi;
    
    fclose(file_id);
    return(1);
  };
    
};




class CxImage
{
 public:
  Matrix<short> r;
  Matrix<short> g;
  Matrix<short> b;
  

 CxImage(unsigned int row, unsigned int col)
   : r(row, col), g(row, col), b(row, col)
    {
    };
  ~CxImage(void)
    {
    };

  inline Matrix<short>     &red   (void) { return r; };
  inline Matrix<short>     &green (void) { return g; };
  inline Matrix<short>     &blue  (void) { return b; };
  inline unsigned int  rows  (void) const { return r.rows();   };
  inline unsigned int  cols(void) const { return r.cols(); };

  void write_C_image(const char *filename){
    unsigned int height = rows();
    unsigned int width  = cols();

    std::ofstream out(filename, ios::out);

    out << "P6 " <<  width <<" " << height <<" " <<  "255" << std::endl;

    for (unsigned int row = 0; row < height; ++row)
      for (unsigned int column = 0; column < width; ++column){
	char dog = (char) r[row][column];
	out.write( (char *) &dog, sizeof(char) );
	dog = (char) g[row][column];
	out.write( (char *) &dog, sizeof(char) );
	dog = (char) b[row][column];
	out.write( (char *) &dog, sizeof(char) );
      }

    out.close();

    return;
  };
};


extern int jpeg_read( const char * filename, CxImage &im);
extern int jpeg_write( const char * filename, int quality, CxImage &cim);
extern int jpeg_read_grey( const char * filename, Matrix<short> &im);

   
int readPPMx(char *filename, Matrix<short> &grau);
int decodeTwoSteps( char* filename, CxImage &im);
int read_pgm( Matrix<short> &im, char *fname);




#endif
