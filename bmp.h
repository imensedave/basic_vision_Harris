#ifndef bmpx_h
#define bmpx_h


typedef struct tagBITMAPFILEHEADER2 {
  short  bfType;
  unsigned int   bfSize;
  short  bfReserved1;
  short  bfReserved2;
  unsigned int   bfOffBits;
} BITMAPFILEHEADER2, *PBITMAPFILEHEADER2; 

#define LONGx int

typedef struct tagBITMAPINFOHEADER2
{
unsigned int biSize;  //specifies the number of bytes required by the struct
int biWidth;  //specifies width in pixels
int biHeight;  //species height in pixels
short biPlanes; //specifies the number of color planes, must be 1
short biBitCount; //specifies the number of bit per pixel
int biCompression;//spcifies the type of compression
int biSizeImage;  //size of image in bytes
int biXPelsPerMeter;  //number of pixels per meter in x axis
int biYPelsPerMeter;  //number of pixels per meter in y axis
int biClrUsed;  //number of colors used by th ebitmap
int biClrImportant;  //number of colors that are important
}BITMAPINFOHEADER2;


typedef struct bmpfile_magic {
  unsigned char magic[2];
} BMPX;
 

#define uint32_t int
#define uint16_t short

typedef struct bmpfile_header {
  uint32_t filesz;
  uint16_t creator1;
  uint16_t creator2;
  uint32_t bmp_offset;
} BMPY;


unsigned char *LoadBitmapFile(char *filename, BITMAPINFOHEADER2 *bitmapInfoHeader);


#endif
