# Makefile  for line models

LIBDIR = /usr/lib
INCDIR = /usr/include
STLDIR = /usr/include




#CFLAGS = -g -I$(INCDIR) -I$(STLDIR) -L$(LIBDIR) -Wno-deprecated -I/home/dasx/source/jpeg-6b
#
CFLAGS = -O3 -I$(INCDIR) -I$(STLDIR) -L$(LIBDIR) -Wno-deprecated -I/home/dasx/source/jpeg-6b
#CFLAGS = -g -I$(INCDIR) -I$(STLDIR) -L$(LIBDIR) -Wno-deprecated -I/Users/das/source/jpeg-6b
#CFLAGS = -O3 -I$(INCDIR) -I$(STLDIR) -L$(LIBDIR) -Wno-deprecated -I/Users/das/source/jpeg-6b

# for raspberry Pi build.
#CFLAGS = -g -I$(INCDIR) -I$(STLDIR) -L$(LIBDIR) -Wno-deprecated -I/home/pi/source/jpeg-6b -mfpu=neon
#CFLAGS = -O3 -I$(INCDIR) -I$(STLDIR) -L$(LIBDIR) -Wno-deprecated -I/home/pi/source/jpeg-6b -mfpu=neon -std=c++11 -ffast-math -ftree-vectorize 


#LDFLAGS = -L/Users/das/source/jpeg-6b
#
LDFLAGS = -L/home/dasx/source/jpeg-6b
#LDFLAGS = -L/home/pi/source/jpeg-6b
# LIBS   = -lm -llibjpeg -lpthread
LIBS   = -lm -ljpeg -lpthread
OBJS  = main.o jpegread.o \
example_encode.o example_decode.o lodepng.o \
read_ppm.o read_bmp.o \
IM_utilities.o \
br_Harris.o \

HEADS = bmp.h  IM_func.h Imgx.h lodepng.h  IM_struct.h


CC = g++


# Building rules for Makefile

stx: $(OBJS) $(HEADS)
	g++ $(CFLAGS) $(LDFLAGS) -o $@ $(OBJS)  $(LIBS) 

%.o: %.C
	g++ $(CFLAGS) -c  $< -o $@

%.o: %.cpp
	g++ $(CFLAGS) -c  $< -o $@

clean:
	-rm *.o stx

main.o:			main.cpp 		$(HEADS)
jpegread.o:		jpegread.C 		$(HEADS)

example_encode.o:	example_encode.cpp	$(HEADS)
example_decode.o:	example_decode.cpp 	$(HEADS)
lodepng.o:		lodepng.cpp 		$(HEADS)
read_ppm.o:		read_ppm.cpp 		$(HEADS)
read_bmp.o:		read_bmp.cpp 		$(HEADS)

IM_utilities.o:		IM_utilities.cpp 	$(HEADS)
br_canny.o:		br_canny.cpp		$(HEADS)
