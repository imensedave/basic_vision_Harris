#ifndef _IM_struct_h
#define _IM_struct_h

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

#include "Imgx.h"



struct Link;
struct Snort;
struct UDL;
struct Snork_int;
struct PT2;
 


struct Exline{
  int r0,c0, r1,c1, id, num, e0, e1, up, right;
  float mr,mc;
  float dr,dc, cvtr, mag;
  Exline()  : r0(0), c0(0), r1(0), c1(0), e0(0), e1(0), up(0), right(0),
	      id(0), num(0), dr(0.0), dc(0.0), mr(0.0), mc(0.0), cvtr(0.0) ,mag(0.0) 
  {};
  
  Exline(const Exline &p) : r0(p.r0), c0(p.c0), r1(p.r1), c1(p.c1), 
			    e0(p.e0), e1(p.e1), up(p.up), right(p.right),
			    id(p.id), num(p.num),  mr(p.mr), mc(p.mc), dr(p.dr), dc(p.dc), 
			    cvtr(p.cvtr), mag(p.mag) {};
  
  Exline &operator=(const Exline &p)
  {
    r0 = p.r0;    c0 = p.c0;
    r1 = p.r1;    c1 = p.c1;
    e0 = p.e0;    e1 = p.e1;
    up = p.up;    right = p.right;
    id = p.id;    num = p.num;

    mr = p.mr;    mc = p.mc;
    dr = p.dr;    dc = p.dc;
       
    cvtr = p.cvtr;
    mag = p.mag;
    return *this;
  };
};


// sort by increasing size of .num
inline bool operator<(const Exline &a, const Exline &b) 
{ 
	return a.mag < b.mag;
};




/* for use in fortune.C */
struct Link
{
	int a,b, c;
	float ro;
	Link(): ro(0.0){};
	Link(int x, int y) : a(x), b(y),c(0), ro(0.0){};
	Link(int x, int y, int z) : a(x), b(y), c(z), ro(0.0) {};
	Link(int x, int y, int z, float q) : a(x), b(y), c(z), ro(q) {};
	Link(const Link &p) : a(p.a), b(p.b), c(p.c), ro(p.ro) {};

	Link &operator=(const Link &p)
	{
		a = p.a;
		b = p.b;
		c = p.c;
		ro = p.ro;
		return *this;
	};
};

inline bool operator<(const Link &x, const Link &y)
{
	return x.ro < y.ro;
};


inline bool operator==(const Link &x, const Link &y)
{
	return (x.a == y.a && x.b == y.b && x.c == y.c);
};




struct Grot
{
	int a,b, c;
	std::vector<int> idx;
	Grot() {};
	Grot(int x, int y) : a(x), b(y),c(0), idx(){};
	Grot(int x, int y, int z) : a(x), b(y), c(z), idx() {};
	Grot(const  Grot &p) : a(p.a), b(p.b), c(p.c), idx(p.idx) {};

	Grot &operator=(const  Grot &p)
	{
		a = p.a;
		b = p.b;
		c = p.c;
		idx = p.idx;
		return *this;
	};
};

// memory kind clustering container 
// groups a max of 303 things together. (one line of text)


struct Grot_mem
{
  int a,b, c;
  int idx[156];
  Grot_mem() {
    memset((void *) idx, 0, sizeof(int) * 148);
  };
  Grot_mem(int x, int y) : a(x), b(y),c(0){
    memset((void *) idx, 0, sizeof(int) * 148);
  };
  Grot_mem(int x, int y, int z) : a(x), b(y), c(z) {
    memset((void *) idx, 0, sizeof(int) * 148);
  };
  Grot_mem(const  Grot_mem &p) : a(p.a), b(p.b), c(p.c){
    memcpy((void *) idx, (void *) p.idx, sizeof(int) *148);
  };

  Grot_mem &operator=(const  Grot_mem &p)
  {
    a = p.a;
    b = p.b;
    c = p.c;
    memcpy((void *) idx, (void *) p.idx, sizeof(int) * 148);
    //idx = p.idx;
    return *this;
  };
};







struct Snort
{
  int a,b, c;
  float mag;
  Snort(): a(0), b(0), c(0), mag(0) {};
  Snort(int x, int y) : a(x), b(y),c(0), mag(0){};
  Snort(int x, int y, int z) : a(x), b(y), c(z), mag(0) {};
  Snort(int x, int y, int z, float g) : a(x), b(y), c(z), mag(g) {};
  Snort(const Snort &p) : a(p.a), b(p.b), c(p.c), mag(p.mag) {};
  
  Snort &operator=(const Snort &p)
  {
		a = p.a;
		b = p.b;
		c = p.c;
		mag = p.mag;
		return *this;
	};
};

inline bool operator<(const Snort &x, const Snort &y)
{
	return x.mag < y.mag;
};

struct Snork
{
	int a,b;
	float mag;
    Snork(): a(0), b(0), mag(0){};

	Snork(const Snork &p) : a(p.a),b(p.b),  mag(p.mag) {};

	Snork &operator=(const Snork &p)
	{
		a = p.a;
		b = p.b;
		mag = p.mag;
		return *this;
	};
};

inline bool operator<(const Snork &x, const Snork &y)
{
	return x.mag < y.mag;
};


inline char safeChar(char c) {
	if (c<=20) return ' ';
	else return c;
}




struct Snork_int
{
	int a,b;
  int mag, val;
	Snork_int() {};

Snork_int(const Snork_int &p) : a(p.a),b(p.b),  mag(p.mag), val(p.val) {};

	Snork_int &operator=(const Snork_int &p)
	{
		a = p.a;
		b = p.b;
		mag = p.mag;
		val = p.val;
		return *this;
	};
};

inline bool operator<(const Snork_int &x, const Snork_int &y)
{
	return x.mag < y.mag;
};

struct Snork2
{
  int a,b;
  float mag, val;
	Snork2() {};

  Snork2(const Snork2 &p) : a(p.a),b(p.b), val(p.val),  mag(p.mag) {};

	Snork2 &operator=(const Snork2 &p)
	{
		a = p.a;
		b = p.b;
		val = p.val;
		mag = p.mag;
		return *this;
	};
};

inline bool operator<(const Snork2 &x, const Snork2 &y)
{
	return x.mag < y.mag;
};


struct Line
{
  float r1,r2,c1,c2;
  int id, v1,v2;   // chain_id, start and end point ids
  float s,w;   // stripe percentage for line and width
  float pr,pc;  // perpendicular to the line (pointing into region...)
  float length;

  std::set<int> ids;
  Line(float e1, float e2, float e3, float e4, int x1) :
     r1(e1), r2(e2),c1(e3), c2(e4),id(x1),
     v1(0), v2(0), s(0.0), w(0.0), pr(0.0), pc(0.0), length(0.0), ids()
    {};
  
  Line (void): r1(0.0), r2(0.0),c1(0.0), c2(0.0), id(-1),
               v1(0), v2(0), s(0.0), w(0.0), pr(0.0), pc(0.0), length(0.0), ids() 
  {};
  Line (const Line &m): r1(m.r1), r2(m.r2),c1(m.c1), c2(m.c2),id(m.id), 
                        v1(m.v1), v2(m.v2), s(m.s), w(m.w),  pr(m.pr),pc(m.pc),
      length(m.length), ids(m.ids)
  {};
  Line &operator=(const Line &m)  {
    if(&m != this){
		 r1 = m.r1; r2 = m.r2;
		 c1 = m.c1; c2 = m.c2; 
       id = m.id;
		 v1 = m.v1; v2 = m.v2; 
		 s = m.s; w = m.w;
		 pr = m.pr; pc = m.pc;
		 length = m.length;
       ids = m.ids;
	 }
	 return( *this );
  }
};


// standard boundary follower structure.
struct UDL
{
	int ui, uj, di, dj, lab;

	UDL (void): ui(0),uj(0),di(0),dj(0), lab(0) {};
	UDL (const UDL &m): ui(m.ui),  uj(m.uj), di(m.di), dj(m.dj),  lab(m.lab)
	{ };
	UDL &operator=(const UDL &m)  {
		if(&m != this){
			ui = m.ui;      uj = m.uj;
			di = m.di;      dj = m.dj;
			lab = m.lab;
		}
		return( *this );
	}
};


struct UDL_chain
{
  int           num,    label;          // number of pixels in chain, region label
  std::vector<UDL>           udl;           // pixels
  std::vector<Line>          lines;        // lines found on the region boundary
  std::set<int>              bpts;        // set of break points, relative to udl vect
   

  UDL_chain (void): num(0), label(0), udl(), lines(), bpts()
  {};
  
  UDL_chain (const UDL_chain &m): num(m.num), label(m.label), udl( m.udl), 
    lines(m.lines), bpts(m.bpts)
  {};
  
  UDL_chain &operator=( const UDL_chain &m){
    if(&m != this) {
      num   = m.num; 
      label = m.label; 
      udl   = m.udl;
      lines = m.lines;
      bpts  = m.bpts; 
    }
    return( *this );
  }

};

struct PT 
{
  int r,c, v;
  PT() {};
  PT(int a, int b) : r(a), c(b), v(0) {};
  PT(int a, int b, int x) : r(a), c(b), v(x) {};
  PT(const PT &p) : r(p.r), c(p.c), v(p.v) {};

  PT &operator=(const PT &p)
  {
    r = p.r;
    c = p.c;
    v = p.v;
    return *this;
  };

};

      
inline bool operator<(const PT &a, const PT &b)
{
  if( a.r == b.r)
    return a.c < b.c;
  return a.r < b.r;
};

inline bool operator==(const PT &x, const PT &y)
{
  return (x.r == y.r && x.c == y.c);
};


inline bool operator!=(const PT &x, const PT &y)
{
  return (x.r != y.r || x.c != y.c);
};


struct Pf2 
{
  float r,c;
  Pf2(): r(0), c(0) {};
  Pf2(float a, float b) : r(a), c(b) {};
  Pf2(float a, float b, int x) : r(a), c(b) {};
  Pf2(const Pf2 &p) : r(p.r), c(p.c) {};
  
  Pf2 &operator=(const Pf2 &p)
  {
    r = p.r;
    c = p.c;
    return *this;
  };
  
};


inline bool operator<(const Pf2 &a, const Pf2 &b)
{
  if( a.r == b.r)
    return a.c < b.c;
  return a.r < b.r;
};


struct PTr 
{
  int *r,*c, v,w;
  PTr() {};
  PTr(int *a, int *b) : r(a), c(b), v(0), w(0) {};
  PTr(int *a, int *b, int x) : r(a), c(b), v(x), w(0) {};
  PTr(int *a, int *b, int x, int y) : r(a), c(b), v(x), w(y) {};
  PTr(const PTr &p) : r(p.r), c(p.c), v(p.v), w(p.w) {};

  PTr &operator=(const PTr &p)
  {
    r = p.r;
    c = p.c;
    v = p.v;
    w = p.w;
    return *this;
  };
};


struct PT_chain
{
  int           num,    label;          // number of pixels in chain, region label
  std::vector<PT>         p;           // pixels
 

  PT_chain (void): num(0), label(0), p()
  {};
  
  PT_chain (const PT_chain &m): num(m.num), label(m.label), p( m.p)
  
  {};
  
  PT_chain &operator=( const PT_chain &m){
    if(&m != this) {
      num   = m.num; 
      label = m.label; 
      p   = m.p;
    }
    return( *this );
  }

};



struct PX
{
  int r,c, v, mag;
  PX() {};
PX(int a, int b) : r(a), c(b), v(0), mag(0) {};
PX(int a, int b, int x, int y) : r(a), c(b), v(x), mag(y) {};
PX(const PX &p) : r(p.r), c(p.c), v(p.v), mag(p.mag) {};

  PX &operator=(const PX &p)
  {
    r = p.r;
    c = p.c;
    v = p.v;
    mag = p.mag;
    return *this;
  };

};

      
inline bool operator<(const PX &a, const PX &b)
{
  return a.mag > b.mag;
};








struct Squat{
  float r0,c0, r1,c1;
  float mr,mc;
  float dr,dc;
  Squat()  : r0(0), c0(0),r1(0), c1(0),  dr(0.0), dc(0.0), mr(0.0), mc(0.0) 
  {};

  Squat(const Squat &p) : r0(p.r0), c0(p.c0), r1(p.r1), c1(p.c1),  mr(p.mr), mc(p.mc), dr(p.dr), dc(p.dc) {};
  
  Squat &operator=(const Squat &p)
  {
    r0 = p.r0;
    c0 = p.c0;
    r1 = p.r1;
    c1 = p.c1;

    mr = p.dr;
    mc = p.dc;
    dr = p.dr;
    dc = p.dc;
	
    return *this;
  };
};




struct PTf 
{
  float r,c, mag; 
  int v;
  PTf() {};
  PTf(float a, float b) : r(a), c(b), v(0), mag(0.0) {};
  PTf(float a, float b, int x) : r(a), c(b), v(x), mag(0.0) {};
  PTf(const PTf &p) : r(p.r), c(p.c), v(p.v), mag(p.mag) {};

  PTf &operator=(const PTf &p)
  {
    r = p.r;
    c = p.c;
    v = p.v;
    mag = p.mag;
    return *this;
  };

};

      
inline bool operator<(const PTf &a, const PTf &b)
{
  //  if( a.r == b.r)
  //return a.c < b.c;
  return a.mag < b.mag;
};



struct PT2 
{
  int r,c;
  PT2() {};
  PT2(int a, int b) : r(a), c(b) {};
  PT2(int a, int b, int x) : r(a), c(b) {};
  PT2(const PT2 &p) : r(p.r), c(p.c) {};
  
  PT2 &operator=(const PT2 &p)
  {
    r = p.r;
    c = p.c;
    return *this;
  };
  
};


inline bool operator<(const PT2 &a, const PT2 &b)
{
  if( a.c == b.c)
    return a.r < b.r;
  return a.c < b.c;
};

struct PT4 
{
  int id,id2,r,c,r2,c2,up;
  PT4(): id(0), id2(0), r(0), r2(0), c(0), c2(0), up(0) {};
  PT4(int a, int b, int x1, int x2) : r(a), c(b), id(x1), up(x2)  {};
  PT4(const PT4 &p) : r(p.r),c(p.c),id(p.id),r2(p.r2),c2(p.c2),id2(p.id2),up(p.up) {};
  
  PT4 &operator=(const PT4 &p)
  {
    r = p.r;
    c = p.c;
    id = p.id;

    r2 = p.r2;
    c2 = p.c2;
    id2 = p.id2;

    up = p.up;
    return *this;
  };
  
};


inline bool operator<(const PT4 &a, const PT4 &b)
{
  if( a.c == b.c)
    return a.r < b.r;
  return a.c < b.c;
};




struct Grotid
{
  int id, num;
  int a,b, c;
  float mag;
  std::vector<int> idx;
  Grotid() : id(0), num(0), a(0), b(0),c(0), mag(0), idx() {};

  Grotid(const  Grotid &p) : id(p.id), num(p.num), a(p.a), b(p.b), c(p.c), mag(p.mag), idx(p.idx) {};
  
  Grotid &operator=(const  Grotid &p)
  {
    a = p.a;
    b = p.b;
    c = p.c;
    idx = p.idx;
    num = p.num;
    mag = p.mag;
    return *this;
  };
  
};

inline bool operator<(const Grotid &x, const Grotid &y)
{
  return x.mag < y.mag;
};



struct PT3
{
  int r,c, rx, cx, id, h;
  PT3() : r(0), c(0), rx(0), cx(0), id(0), h(0){  };

  PT3(int a, int b, int g,  int j, int i, int k) : r(a), c(b), rx(g), cx(j), id(i), h(k) {}; 
  PT3(const PT3 &p) : r(p.r), c(p.c), rx(p.rx), cx(p.cx), id(p.id), h(p.h) {};

  PT3 &operator=(const PT3 &p)
  {
    r = p.r;
    c = p.c;
    rx = p.rx;
    cx = p.cx;
    id = p.id;
    h  = p.h;
    return *this;
  };

};


inline bool operator<(const PT3 &a, const PT3 &b)
{
  if( a.cx == b.cx)
    return a.rx < b.rx;
  return a.cx < b.cx;
};


struct SLUG
{
  int id,id2, n1, n2,  r,c,r2,c2;
  float m1, m2;
  SLUG(): n1(-1), n2(-1){};
  SLUG(int a, int b, int x1, int x2) : r(a), c(b), id(x1)  {};
  SLUG(const SLUG &p) : r(p.r),c(p.c),
    id(p.id),id2(p.id2),
    n1(p.n1),n2(p.n2),
    r2(p.r2),c2(p.c2),
    m1(p.m1), m2(p.m2) {};
  
  SLUG &operator=(const SLUG &p)
  {
    r = p.r;
    c = p.c;
    id = p.id;
    id2 = p.id2;

    n1 = p.n1; n2 = p.n2;
    
    r2 = p.r2;
    c2 = p.c2;

    m1 = p.m1;
    m2 = p.m2;
    
    return *this;
  };
  
};


inline bool operator<(const SLUG &a, const SLUG &b)
{
  return a.c < b.c;
};



struct ExlineV{
  int r0,c0, r1,c1, id, num, e0, e1, up, right;
  float mr,mc;
  float dr,dc, cvtr, mag;
  std::vector<PT2> pts; 
  ExlineV()  : r0(0), c0(0), r1(0), c1(0), e0(0), e1(0), up(0), right(0),
    id(0), num(0), dr(0.0), dc(0.0), mr(0.0), mc(0.0), cvtr(0.0) ,mag(0.0),
    pts()
  {};
  
  ExlineV(const ExlineV &p) : r0(p.r0), c0(p.c0), r1(p.r1), c1(p.c1), 
    e0(p.e0), e1(p.e1), up(p.up), right(p.right),
    id(p.id), num(p.num),  mr(p.mr), mc(p.mc), dr(p.dr), dc(p.dc), 
    cvtr(p.cvtr), mag(p.mag),
    pts(p.pts) {};
  
  ExlineV &operator=(const ExlineV &p)
  {
    r0 = p.r0;    c0 = p.c0;
    r1 = p.r1;    c1 = p.c1;
    e0 = p.e0;    e1 = p.e1;
    up = p.up;    right = p.right;
    id = p.id;    num = p.num;

    mr = p.mr;    mc = p.mc;
    dr = p.dr;    dc = p.dc;
       
    cvtr = p.cvtr;
    mag = p.mag;
    pts = p.pts;
    return *this;
  };


  void clear() {
    pts.clear();
  }
  
  ~ExlineV(void){
    pts.clear();
  };
  
};



inline bool operator<(const ExlineV &a, const ExlineV &b) 
{ 
	return a.mag < b.mag;
};












#endif

