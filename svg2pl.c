/*
+----------+
| svg2pl.c |
+----------+

 convert (reasonable subset of) svg to polylines

- uses very little memory, 0 dynamic allocations
- single file, no dependencies other than C
- extremely fast
- no error checking -- assumes svg is perfectly valid, otherwise UB
- intended for compact projects or embeded systems

supported svg features:
- transform attribute on all elements
- all path d= commands (beziers, arcs...)
- elements: g, path, line, polyline, polygon, circle, ellipse, rect, svg
- you can register custom parser functions for other tags (see usage)

limitations:
- curves and ellipses are discretized
- fills and other styles are ignored,
  you can choose whether to skip all elements without "stoke" attribute
- viewBox and nested svg are supported, 
  but preserveAspectRatio can only be "xMidYMid meet"(default) or "none"
- percentage(%) and physical units are not supported
- JS, CSS and other crazy stuff are obviously not supported

usage:
1. first define macros S2P_MOVETO(x,y) and S2P_LINETO(x,y)
  - these define what to do when the library has extracted polyline vertices
  - you can write to your own data structure, print, or control a machine, etc.
2. (optional) define macro S2P_SETDIM(w,h)
  - this define what to do when the width/height of the svg are extracted
  - e.g. you might scale all subsequent vertices
3. #include this file
4. (optional) call s2p_def_tag to register parsers for custom tags (e.g. <text)
5. call s2p_parse_from_file OR s2p_parse_from_str OR s2p_parse (FILE*)

(c) lingdong huang 2023, at MIT. MIT License
*/

#ifndef SVG2PL_C
#define SVG2PL_C

#include <stdio.h>
#include <string.h>
#include <math.h>

#define S2P_TAG_OPEN   0
#define S2P_TAG_SELFCL 1
#define S2P_TAG_CLOSE  2

#define S2P__ST_OUT  0
#define S2P__ST_OPEN 1
#define S2P__ST_TAG  2
#define S2P__ST_IN   3
#define S2P__ST_AK   4
#define S2P__ST_AKE  5
#define S2P__ST_AVQ  6
#define S2P__ST_AV   7
#define S2P__ST_EXCL 8
#define S2P__ST_CMMT 9

#define S2P_MAX_ATTR    32
#define S2P_MAX_AK_LEN  7
#define S2P_MAX_TAG_LEN  7

#define S2P_MAX_MAT 32

#define S2P_MAX_CUSTFUN 16

#define S2P__MAX(a,b) ((a) > (b) ? (a) : (b))
#define S2P__MIN(a,b) ((a) < (b) ? (a) : (b))

#ifndef S2P_MOVETO
#define S2P_MOVETO(x,y) {printf("M %f %f\n",(x),(y));}
#endif

#ifndef S2P_LINETO
#define S2P_LINETO(x,y) {printf("L %f %f\n",x,y);}
#endif

#ifndef S2P_SETDIM
#define S2P_SETDIM(w,h) {printf("M 0 0 h %f v %f h %f z\n",(w),(h),-(w));}
#endif

int s2p_reso_curv = 24;
int s2p_reso_circ = 24;
int s2p_skip_nostroke = 0;

typedef struct s2p_attr_st {
  int k_len;
  int v_len;
  long k;
  long v;
} s2p_attr_t;

typedef void (* s2p_tagdef_funptr_t) (FILE* fd, s2p_attr_t* attrs, int n_attr);

typedef struct s2p_tagdef_st{
  char type;
  char name[S2P_MAX_TAG_LEN];
  s2p_tagdef_funptr_t func;
} s2p_tagdef_t;

s2p_tagdef_t s2p__tagdefs[S2P_MAX_CUSTFUN];
int s2p__ntd = 0;

float s2p__mats[S2P_MAX_MAT][9];
int s2p__mati = 0;

#define S2P__MATSET(A,a,b,c,d,e,f) {A[0]=a;A[1]=b;A[2]=c;A[3]=d;A[4]=e;A[5]=f;A[6]=0;A[7]=0;A[8]=1;}

void s2p__matmul(float* A, float* B, float* C){
  float a = A[0] * B[0] + A[1] * B[3] + A[2] * B[6];
  float b = A[0] * B[1] + A[1] * B[4] + A[2] * B[7];
  float c = A[0] * B[2] + A[1] * B[5] + A[2] * B[8];

  float d = A[3] * B[0] + A[4] * B[3] + A[5] * B[6];
  float e = A[3] * B[1] + A[4] * B[4] + A[5] * B[7];
  float f = A[3] * B[2] + A[4] * B[5] + A[5] * B[8];

  float g = A[6] * B[0] + A[7] * B[3] + A[8] * B[6];
  float h = A[6] * B[1] + A[7] * B[4] + A[8] * B[7];
  float i = A[6] * B[2] + A[7] * B[5] + A[8] * B[8];

  C[0] = a; C[1] = b; C[2] = c;
  C[3] = d; C[4] = e; C[5] = f;
  C[6] = g; C[7] = h; C[8] = i;
}

void s2p__m_moveto(float x, float y){
  float* m = s2p__mats[s2p__mati];
  float x1 = m[0] * x + m[1] * y + m[2];
  float y1 = m[3] * x + m[4] * y + m[5];
  float z1 = m[6] * x + m[7] * y + m[8];
  S2P_MOVETO(x1/z1,y1/z1);
}

void s2p__m_lineto(float x, float y){
  float* m = s2p__mats[s2p__mati];
  float x1 = m[0] * x + m[1] * y + m[2];
  float y1 = m[3] * x + m[4] * y + m[5];
  float z1 = m[6] * x + m[7] * y + m[8];
  S2P_LINETO(x1/z1,y1/z1);
}


void s2p_quadratic_bezier(
  float x0,float y0,float x1,float y1,
  float x2,float y2,
  float t, float* xo, float* yo){
  float s = 1-t;
  float s2 = s*s;
  float t2 = t*t;
  (*xo) = s2*x0+2*s*t*x1+t2*x2;
  (*yo) = s2*y0+2*s*t*y1+t2*y2;
}

void s2p_cubic_bezier(
  float x0,float y0,float x1,float y1,
  float x2,float y2,float x3,float y3,
  float t, float* xo, float* yo){
  float s = 1-t;
  float s2 = s*s;
  float s3 = s*s2;
  float t2 = t*t;
  float t3 = t2*t;
  (*xo) = s3*x0+3*s2*t*x1+3*s*t2*x2+t3*x3;
  (*yo) = s3*y0+3*s2*t*y1+3*s*t2*y2+t3*y3;
}


float s2p__aang(float ux, float uy, float vx, float vy){
  float dot = ux*vx+uy*vy;
  float den = sqrt(ux*ux+uy*uy)*sqrt(vx*vx+vy*vy);
  float a = acos(dot/den);
  float c = ux*vy-uy*vx;
  if (c < 0){
    a *= -1;
  }
  return a;
}

void s2p_arc(
  float x1, float y1, float x2, float y2,
  int fa, int fs,
  float rx, float ry, float phi){
  //https://www.w3.org/TR/SVG/implnote.html#ArcImplementationNotes

  float x125 = (x1-x2)*0.5;
  float y125 = (y1-y2)*0.5;
  float cosph = cos(phi);
  float sinph = sin(phi);

  float x1p = cosph*x125+sinph*y125;
  float y1p =-sinph*x125+cosph*y125;
  float x1p2 = x1p*x1p;
  float y1p2 = y1p*y1p;

  if (rx == 0 || ry == 0){
    return s2p__m_lineto(x2,y2);
  }
  if (rx < 0) rx = -rx;
  if (ry < 0) ry = -ry;

  float lam = x1p2/(rx*rx) + y1p2/(ry*ry);
  if (lam >= 1.0){
    float sql = sqrt(lam);
    rx = sql*rx;
    ry = sql*ry;
  }

  float rx2 = rx*rx;
  float ry2 = ry*ry;

  float rr = (rx2*y1p2+ry2*x1p2);
  float rrr = fabs(rx2*ry2-rr); //or fmax(rx2*ry2-rr,0)??
  float sqrrr = sqrt(rrr/rr);

  float cxp = sqrrr * (rx*y1p) / ry;
  float cyp =-sqrrr * (ry*x1p) / rx;
  if (fa == fs){
    cxp *= -1;
    cyp *= -1;
  }

  float cx = cosph*cxp - sinph*cyp + (x1+x2)*0.5;
  float cy = sinph*cxp + cosph*cyp + (y1+y2)*0.5;
  float xpcr = (x1p-cxp)/rx;
  float ypcr = (y1p-cyp)/ry;

  float th1 = s2p__aang(1,0, xpcr, ypcr );
  float dth = s2p__aang(xpcr, ypcr,  (-x1p-cxp)/rx, (-y1p-cyp)/ry);
  dth = fmod(dth+M_PI*4, M_PI*2);
  if (fs == 0){
    while (dth > 0){
      dth -= M_PI*2;
    }
  }else{
    while (dth < 0){
      dth += M_PI*2;
    }
  }

  for (int i = 0; i < s2p_reso_curv; i++){
    float t = (float)(i+1)/(float)s2p_reso_curv;
    float th = th1 + dth * t;
    float x0 = cos(th)*rx;
    float y0 = sin(th)*ry;
    float xx = cx+cosph*x0-sinph*y0;
    float yy = cy+sinph*x0+cosph*y0;
    s2p__m_lineto(xx,yy);
  }

}

int s2p__d_next(FILE* fd){
  char c;
  while ( c=fgetc(fd), c==','||c<=' '){}
  int is_n = '+' <= c && c <= '9';
  ungetc(c,fd);
  return is_n;
}


void s2p__parse_transf(FILE* fd, long v_pos, int v_len, float* mat){
  fseek(fd, v_pos, SEEK_SET);

  char c;
  float x0,x1,x2;
  float m [9] = {1,0,0,0,1,0,0,0,1};
  float m1[9] = {1,0,0,0,1,0,0,0,1};

  while (ftell(fd) < v_pos + v_len){
    s2p__d_next(fd);
    c = fgetc(fd);
    if (c == 't'){//ranslate
      fseek(fd,8,SEEK_CUR);
      s2p__d_next(fd);
      fgetc(fd);
      s2p__d_next(fd);
      fscanf(fd, "%f", &x0);
      s2p__d_next(fd);
      fscanf(fd, "%f", &x1);
      s2p__d_next(fd);
      fgetc(fd);

      S2P__MATSET(m, 1,0,x0, 0,1,x1);
      s2p__matmul(mat,m,mat);
    }else if (c == 'r'){//otate
      fseek(fd,5,SEEK_CUR);
      s2p__d_next(fd);
      fgetc(fd);
      s2p__d_next(fd);
      fscanf(fd, "%f", &x0);
      float th = x0 * M_PI / 180.0;
      S2P__MATSET(m, cos(th), -sin(th), 0, sin(th), cos(th), 0);

      if (s2p__d_next(fd)){
        fscanf(fd, "%f", &x1);
        s2p__d_next(fd);
        fscanf(fd, "%f", &x2);
        s2p__d_next(fd);
        m[2] = x1;
        m[5] = x2;
        S2P__MATSET(m1, 1,0,-x1, 0,1,-x2);
        s2p__matmul(m,m1,m);
      }
      s2p__matmul(mat,m,mat);
      fgetc(fd);
    }else if (c == 'm'){//atrix
      fseek(fd,5,SEEK_CUR);
      s2p__d_next(fd);
      fgetc(fd);
      s2p__d_next(fd);

      for (int i = 0; i < 6; i++){
        fscanf(fd, "%f", &x0);
        m[(i&1)*3+(i>>1)] = x0; //col maj
        s2p__d_next(fd);
      }
      m[6] = 0; m[7] = 0; m[8] = 1;
      fgetc(fd);
      s2p__matmul(mat,m,mat);

    }else if (c == 's'){
      fseek(fd,3,SEEK_CUR);
      c = fgetc(fd);
      if (c == 'e'){//scale
        s2p__d_next(fd);
        fgetc(fd);
        s2p__d_next(fd);
        fscanf(fd, "%f", &x0);

        S2P__MATSET(m, x0,0,0, 0,x0,0);

        if (s2p__d_next(fd)){
          fscanf(fd, "%f", &x1);
          s2p__d_next(fd);

          m[4] = x1;
        }
        s2p__matmul(mat,m,mat);
        fgetc(fd);
      }else if (c == 'X'){//skewX
        s2p__d_next(fd);
        fgetc(fd);
        s2p__d_next(fd);
        fscanf(fd, "%f", &x0);
        s2p__d_next(fd);
        fgetc(fd);

        S2P__MATSET(m, 1,tan(x0*M_PI/180.0),0, 0,1,0);
        s2p__matmul(mat,m,mat);

      }else if (c == 'Y'){//skewY
        s2p__d_next(fd);
        fgetc(fd);
        s2p__d_next(fd);
        fscanf(fd, "%f", &x0);
        s2p__d_next(fd);
        fgetc(fd);

        S2P__MATSET(m, 1,0,0, tan(x0*M_PI/180.0),1,0);
        s2p__matmul(mat,m,mat);
      }

    }
  }  
}




void s2p__get_transf(FILE* fd, s2p_attr_t* attrs, int n_attr, float* mat){
  for (int i = 0; i < n_attr; i++){
    char keybuf[S2P_MAX_AK_LEN+1] = {0};
    fseek(fd, attrs[i].k, SEEK_SET);
    fread(keybuf,1,S2P__MIN(attrs[i].k_len,S2P_MAX_AK_LEN),fd);
    if (!strncmp(keybuf,"tra",3)){
      s2p__parse_transf(fd,attrs[i].v,attrs[i].v_len,mat);
      return;
    }
  }
}

int s2p__get_stroke(FILE* fd, s2p_attr_t* attrs, int n_attr){
  for (int i = 0; i < n_attr; i++){
    char keybuf[S2P_MAX_AK_LEN+1] = {0};
    fseek(fd, attrs[i].k, SEEK_SET);
    fread(keybuf,1,S2P__MIN(attrs[i].k_len,S2P_MAX_AK_LEN),fd);
    
    if (!strcmp(keybuf,"stroke")){
      fseek(fd, attrs[i].v, SEEK_SET);
      if (fgetc(fd) != 'n' || fgetc(fd) != 'o'){
        return 1;
      }else{
        return 0;
      }
    }
  }
  return 0;
}

void s2p__parse_line(FILE* fd, s2p_attr_t* attrs, int n_attr){
  
  float x1, y1, x2, y2;
  for (int i = 0; i < n_attr; i++){
    char keybuf[S2P_MAX_AK_LEN+1] = {0};
    fseek(fd, attrs[i].k, SEEK_SET);
    fread(keybuf,1,S2P__MIN(attrs[i].k_len,S2P_MAX_AK_LEN),fd);
    if (!strcmp(keybuf,"x1")){
      fseek(fd, attrs[i].v, SEEK_SET);
      fscanf(fd, "%f", &x1);
    }else if (!strcmp(keybuf,"y1")){
      fseek(fd, attrs[i].v, SEEK_SET);
      fscanf(fd, "%f", &y1);
    }else if (!strcmp(keybuf,"x2")){
      fseek(fd, attrs[i].v, SEEK_SET);
      fscanf(fd, "%f", &x2);
    }else if (!strcmp(keybuf,"y2")){
      fseek(fd, attrs[i].v, SEEK_SET);
      fscanf(fd, "%f", &y2);
    }
  }
  s2p__m_moveto(x1,y1);
  s2p__m_lineto(x2,y2);
}


void s2p__parse_rect(FILE* fd, s2p_attr_t* attrs, int n_attr){
  float x=0,y=0,w,h,rx=0,ry=0;
  for (int i = 0; i < n_attr; i++){
    char keybuf[S2P_MAX_AK_LEN+1] = {0};
    fseek(fd, attrs[i].k, SEEK_SET);
    fread(keybuf,1,S2P__MIN(attrs[i].k_len,S2P_MAX_AK_LEN),fd);
    if (!strcmp(keybuf,"x")){
      fseek(fd, attrs[i].v, SEEK_SET);
      fscanf(fd, "%f", &x);
    }else if (!strcmp(keybuf,"y")){
      fseek(fd, attrs[i].v, SEEK_SET);
      fscanf(fd, "%f", &y);
    }else if (!strcmp(keybuf,"width")){
      fseek(fd, attrs[i].v, SEEK_SET);
      fscanf(fd, "%f", &w);
    }else if (!strcmp(keybuf,"height")){
      fseek(fd, attrs[i].v, SEEK_SET);
      fscanf(fd, "%f", &h);
    }else if (!strcmp(keybuf,"rx")){
      fseek(fd, attrs[i].v, SEEK_SET);
      fscanf(fd, "%f", &rx);
    }else if (!strcmp(keybuf,"ry")){
      fseek(fd, attrs[i].v, SEEK_SET);
      fscanf(fd, "%f", &ry);
    }
  }
  if (rx && !ry){
    ry = rx;
  }
  if (!rx || !ry){
    s2p__m_moveto(x,y);
    s2p__m_lineto(x+w,y);
    s2p__m_lineto(x+w,y+h);
    s2p__m_lineto(x,y+h);
    s2p__m_lineto(x,y);
  }else{
    rx = fabs(rx);
    ry = fabs(ry);
    if (rx > w/2) rx = w/2;
    if (ry > h/2) ry = h/2;
    s2p__m_moveto(x,y+ry);
    s2p_arc(x,y+ry, x+rx,y, 0,1, rx,ry,0);
    s2p__m_lineto(x+w-rx,y);
    s2p_arc(x+w-rx,y, x+w,y+ry, 0,1, rx,ry,0);
    s2p__m_lineto(x+w,y+h-ry);
    s2p_arc(x+w,y+h-ry, x+w-rx,y+h, 0,1, rx,ry,0);
    s2p__m_lineto(x+rx,y+h);
    s2p_arc(x+rx,y+h, x,y+h-ry, 0,1, rx,ry,0);
    s2p__m_lineto(x,y+ry);
  }
}


void s2p__parse_ellipse(FILE* fd, s2p_attr_t* attrs, int n_attr){
  float cx=0,cy=0,rx=0,ry=0;
  for (int i = 0; i < n_attr; i++){
    char keybuf[S2P_MAX_AK_LEN+1] = {0};
    fseek(fd, attrs[i].k, SEEK_SET);
    fread(keybuf,1,S2P__MIN(attrs[i].k_len,S2P_MAX_AK_LEN),fd);
    if (!strcmp(keybuf,"cx")){
      fseek(fd, attrs[i].v, SEEK_SET);
      fscanf(fd, "%f", &cx);
    }else if (!strcmp(keybuf,"cy")){
      fseek(fd, attrs[i].v, SEEK_SET);
      fscanf(fd, "%f", &cy);
    }else if (!strcmp(keybuf,"r")){
      fseek(fd, attrs[i].v, SEEK_SET);
      fscanf(fd, "%f", &rx);
    }else if (!strcmp(keybuf,"rx")){
      fseek(fd, attrs[i].v, SEEK_SET);
      fscanf(fd, "%f", &rx);
    }else if (!strcmp(keybuf,"ry")){
      fseek(fd, attrs[i].v, SEEK_SET);
      fscanf(fd, "%f", &ry);
    }
  }
  if (rx && !ry) ry = rx;
  if (ry && !rx) rx = ry;
  rx = fabs(rx);
  ry = fabs(ry);
  for (int i = 0; i < s2p_reso_circ; i++){
    float t = (float)i/(float)(s2p_reso_circ-1);
    float a = t * M_PI * 2;
    float x = cx+cos(a)*rx;
    float y = cy+sin(a)*ry;
    if (i){
      s2p__m_lineto(x,y);
    }else{
      s2p__m_moveto(x,y);
    }
  }
}

void s2p__parse_poly(FILE* fd, s2p_attr_t* attrs, int n_attr, int closed){
  for (int i = 0; i < n_attr; i++){
    char keybuf[S2P_MAX_AK_LEN+1] = {0};
    fseek(fd, attrs[i].k, SEEK_SET);
    fread(keybuf,1,S2P__MIN(attrs[i].k_len,S2P_MAX_AK_LEN),fd);
    if (!strcmp(keybuf,"points")){
      fseek(fd, attrs[i].v, SEEK_SET);
      float x,y,x0,y0;
      int j=0;
      while (ftell(fd) < attrs[i].v + attrs[i].v_len){
        fscanf(fd, "%f", &x);
        fgetc(fd);
        fscanf(fd, "%f", &y);
        fgetc(fd);
        if (j){
          s2p__m_lineto(x,y);
        }else{
          s2p__m_moveto(x,y);
          x0 = x;
          y0 = y;
        }
        j++;
      }
      if (closed){
        s2p__m_lineto(x0,y0);
      }
    }
  }
}




void s2p__parse_d(FILE* fd, long v_pos, int v_len){
  fseek(fd, v_pos, SEEK_SET);
  char c;
  char lc;
  float lx=0, ly=0;
  float ox=0, oy=0;
  float cx=0, cy=0;
  float qx=0, qy=0;
  float x0,y0;
  float x1,y1;
  float x2,y2;
  int f0,f1;
  while (ftell(fd) < v_pos + v_len){
    s2p__d_next(fd);
    c = fgetc(fd);
    if (c < 'A'){
      ungetc(c,fd);
      c = lc;
    }
    if (c == 'M'){
      fscanf(fd, "%f", &x0);
      s2p__d_next(fd);
      fscanf(fd, "%f", &y0);
      s2p__m_moveto(ox=lx=x0,oy=ly=y0);
      while (s2p__d_next(fd)){
        fscanf(fd, "%f", &x0);
        s2p__d_next(fd);
        fscanf(fd, "%f", &y0);
        s2p__m_lineto(lx=x0,ly=y0);
      }
    }else if (c == 'm'){
      fscanf(fd, "%f", &x0);
      s2p__d_next(fd);
      fscanf(fd, "%f", &y0);
      s2p__m_moveto(ox=(lx+=x0),oy=(ly+=y0));
      while (s2p__d_next(fd)){
        fscanf(fd, "%f", &x0);
        s2p__d_next(fd);
        fscanf(fd, "%f", &y0);
        s2p__m_lineto(lx+=x0,ly+=y0);
      }
    }else if (c == 'L'){
      fscanf(fd, "%f", &x0);
      s2p__d_next(fd);
      fscanf(fd, "%f", &y0);
      s2p__m_lineto(lx=x0,ly=y0);

    }else if (c == 'l'){
      fscanf(fd, "%f", &x0);
      s2p__d_next(fd);
      fscanf(fd, "%f", &y0);
      s2p__m_lineto(lx+=x0,ly+=y0);

    }else if (c == 'H'){
      fscanf(fd, "%f", &x0);
      s2p__m_lineto(lx=x0,ly);

    }else if (c == 'h'){
      fscanf(fd, "%f", &x0);
      s2p__m_lineto(lx+=x0,ly);

    }else if (c == 'V'){
      fscanf(fd, "%f", &y0);
      s2p__m_lineto(lx,ly=y0);

    }else if (c == 'v'){
      fscanf(fd, "%f", &y0);
      s2p__m_lineto(lx,ly+=y0);

    }else if (c == 'Z' || c == 'z'){
      s2p__m_lineto(lx=ox,ly=oy);
    }else if (c == 'C' || c == 'c'){
      float _x = 0;
      float _y = 0;
      if (c == 'c'){
        _x = lx;
        _y = ly;
      }
      fscanf(fd, "%f", &x0);
      s2p__d_next(fd);
      fscanf(fd, "%f", &y0);
      s2p__d_next(fd);
      fscanf(fd, "%f", &x1);
      s2p__d_next(fd);
      fscanf(fd, "%f", &y1);
      s2p__d_next(fd);
      fscanf(fd, "%f", &x2);
      s2p__d_next(fd);
      fscanf(fd, "%f", &y2);
      for (int i = 0; i < s2p_reso_curv; i++){
        float t = (float)(i+1)/(float)(s2p_reso_curv);
        float xt,yt;
        s2p_cubic_bezier(lx,ly,_x+x0,_y+y0,_x+x1,_y+y1,_x+x2,_y+y2,t,&xt,&yt);
        s2p__m_lineto(xt,yt);
      }
      lx = _x+x2;
      ly = _y+y2;
      cx = _x+x1;
      cy = _y+y1;
    }else if (c == 'Q' || c == 'q'){
      float _x = 0;
      float _y = 0;
      if (c == 'q'){
        _x = lx;
        _y = ly;
      }
      fscanf(fd, "%f", &x0);
      s2p__d_next(fd);
      fscanf(fd, "%f", &y0);
      s2p__d_next(fd);
      fscanf(fd, "%f", &x1);
      s2p__d_next(fd);
      fscanf(fd, "%f", &y1);
      for (int i = 0; i < s2p_reso_curv; i++){
        float t = (float)(i+1)/(float)(s2p_reso_curv);
        float xt,yt;
        s2p_quadratic_bezier(lx,ly,_x+x0,_y+y0,_x+x1,_y+y1,t,&xt,&yt);
        s2p__m_lineto(xt,yt);
      }
      lx = _x+x1;
      ly = _y+y1;
      qx = _x+x0;
      qy = _y+y0;
    }else if (c == 'S' || c == 's'){
      float _x = 0, _y = 0;
      if (c == 's'){
        _x = lx;
        _y = ly;
      }
      if (!(lc == 'C' || lc == 'c' || lc == 'S' || lc == 's')){
        cx = lx;
        cy = ly;
      }
      x0 = (lx-cx)+lx;
      y0 = (ly-cy)+ly;
      fscanf(fd, "%f", &x1);
      s2p__d_next(fd);
      fscanf(fd, "%f", &y1);
      s2p__d_next(fd);
      fscanf(fd, "%f", &x2);
      s2p__d_next(fd);
      fscanf(fd, "%f", &y2);
      for (int i = 0; i < s2p_reso_curv; i++){
        float t = (float)(i+1)/(float)(s2p_reso_curv);
        float xt,yt;
        s2p_cubic_bezier(lx,ly,x0,y0,_x+x1,_y+y1,_x+x2,_y+y2,t,&xt,&yt);
        s2p__m_lineto(xt,yt);
      }
      lx = _x+x2;
      ly = _y+y2;
      cx = _x+x1;
      cy = _y+y1;
    }else if (c == 'T' || c == 't'){
      float _x = 0, _y = 0;
      if (c == 't'){
        _x = lx;
        _y = ly;
      }
      if (!(lc == 'Q' || lc == 'q' || lc == 'T' || lc == 't')){
        qx = lx;
        qy = ly;
      }
      x0 = (lx-qx)+lx;
      y0 = (ly-qy)+ly;
      fscanf(fd, "%f", &x1);
      s2p__d_next(fd);
      fscanf(fd, "%f", &y1);
      for (int i = 0; i < s2p_reso_curv; i++){
        float t = (float)(i+1)/(float)(s2p_reso_curv);
        float xt,yt;
        s2p_quadratic_bezier(lx,ly,x0,y0,_x+x1,_y+y1,t,&xt,&yt);
        s2p__m_lineto(xt,yt);
      }
      lx = _x+x1;
      ly = _y+y1;
      qx = x0;
      qy = y0;
    }else if (c == 'A' || c == 'a'){
      float _x = 0, _y = 0;
      if (c == 'a'){
        _x = lx;
        _y = ly;
      }
      fscanf(fd, "%f", &x0);//rx
      s2p__d_next(fd);
      fscanf(fd, "%f", &y0);//ry
      s2p__d_next(fd);
      fscanf(fd, "%f", &x1);//rot
      s2p__d_next(fd);
      fscanf(fd, "%d", &f0);//fA
      s2p__d_next(fd);
      fscanf(fd, "%d", &f1);//fS
      s2p__d_next(fd);
      fscanf(fd, "%f", &x2);//nx
      s2p__d_next(fd);
      fscanf(fd, "%f", &y2);//ny
      
      s2p_arc(lx,ly,_x+x2,_y+y2,f0,f1,x0,y0,x1*M_PI/180.0);
      lx = _x+x2;
      ly = _y+y2;
    }
    lc = c;
  }
}

void s2p__parse_path(FILE* fd, s2p_attr_t* attrs, int n_attr){
  for (int i = 0; i < n_attr; i++){
    char keybuf[S2P_MAX_AK_LEN+1] = {0};
    fseek(fd, attrs[i].k, SEEK_SET);
    fread(keybuf,1,S2P__MIN(attrs[i].k_len,S2P_MAX_AK_LEN),fd);
    if (!strcmp(keybuf,"d")){
      
      s2p__parse_d(fd,attrs[i].v,attrs[i].v_len);
    }
  }
}


void s2p_def_tag(char* name, char tag_type, s2p_tagdef_funptr_t fun){
  strncpy(s2p__tagdefs[s2p__ntd].name, name, S2P_MAX_TAG_LEN);
  s2p__tagdefs[s2p__ntd].type = tag_type;
  s2p__tagdefs[s2p__ntd].func = fun;
  s2p__ntd++;
}

void s2p_parse(FILE* fd){
  int c;
  int state = S2P__ST_OUT;
  long tag;
  char tag_len = 0;
  char closing = 0;
  s2p_attr_t attrs[S2P_MAX_ATTR] = {0};
  int n_attr = 0;
  char quote = 0;
  char cmmttyp = 0;
  int endcmmt = 0;
  s2p__mati = 0;
  S2P__MATSET(s2p__mats[s2p__mati], 1,0,0, 0,1,0);

  long stroke_stack = 0;

  while ((c = fgetc(fd)) != EOF){
    // printf("%c %d\n",c,state);
    if (state == S2P__ST_OUT){
      if (c == '<'){
        state = S2P__ST_OPEN;
      }
    }else if (state == S2P__ST_OPEN){
      if (c == '>') goto lbl_close_tag;
      if (c == '/') {closing = S2P_TAG_CLOSE; continue;}
      if (c == '!'){ state = S2P__ST_EXCL;continue;}
      if (c > ' '){
        state = S2P__ST_TAG;
        tag = ftell(fd)-1;
        tag_len++;
      }
    }else if (state == S2P__ST_TAG){
      if (c == '>') goto lbl_close_tag;
      if (c == '/') goto lbl_slash;
      if (c > ' '){
        tag_len++;
      }else{
        state = S2P__ST_IN;
      }
    }else if (state == S2P__ST_IN){
      if (c == '>'){
        goto lbl_close_tag;
      }
      if (c == '/') goto lbl_slash;
      if (c > ' '){
        state = S2P__ST_AK;
        attrs[n_attr].k = ftell(fd)-1;
        attrs[n_attr].k_len++;
      }
    }else if (state == S2P__ST_AK){
      if (c == '>') goto lbl_close_tag;
      if (c == '/') goto lbl_slash;
      if (c == '='){
        state = S2P__ST_AVQ;
      }else if (c > ' '){
        attrs[n_attr].k_len++;
      }else{
        state = S2P__ST_AKE;
      }
    }else if (state == S2P__ST_AKE){
      if (c == '>') goto lbl_close_tag;
      if (c == '/') goto lbl_slash;
      if (c == '='){
        state = S2P__ST_AVQ;
      }else if (c > ' '){
        n_attr ++;
        state = S2P__ST_AK;
        attrs[n_attr].k = ftell(fd)-1;
        attrs[n_attr].k_len++;
      }
    }else if (state == S2P__ST_AVQ){
      if (c == '>') goto lbl_close_tag;
      if (c == '/') goto lbl_slash;
      if (c == '"' || c == '\''){
        quote = c;
        state = S2P__ST_AV;
        attrs[n_attr].v = ftell(fd);
      }
    }else if (state == S2P__ST_AV){
      if (c == quote){
        n_attr ++;
        state = S2P__ST_IN;
      }else{
        attrs[n_attr].v_len++;
      }
    }else if (state == S2P__ST_EXCL){
      cmmttyp = c;
      if (c == '[') cmmttyp += 2;
      state = S2P__ST_CMMT;
      endcmmt = -1;
    }else if (state == S2P__ST_CMMT){
      if (c == '>' && (endcmmt >= 2 || (cmmttyp != '-' && cmmttyp != ']'))) goto lbl_close_tag;
      if (c == cmmttyp) endcmmt ++; else endcmmt = 0;
    }
    continue;
lbl_slash:;
    closing ++;
    continue;
lbl_close_tag:;

    if (!tag_len){
      goto lbl_reset;
    }

    int pos = ftell(fd);

    char tagbuf[S2P_MAX_TAG_LEN+1] = {0};
    
    fseek(fd, tag, SEEK_SET);
    fread(tagbuf,1,S2P__MIN(tag_len,7),fd);

    // fprintf(stderr,"%s %d\n",tagbuf,pos);

    for (int i = 0; i < s2p__ntd; i++){
      if (closing != s2p__tagdefs[i].type){
        continue;
      }
      if (!strncmp(s2p__tagdefs[i].name, tagbuf, S2P_MAX_TAG_LEN)){
        fseek(fd, pos, SEEK_SET);
        s2p__tagdefs[i].func(fd,attrs,n_attr);
        goto lbl_done_tag;
      }
    }

    if (!strcmp(tagbuf,"svg")){
      if (closing == 0){
        float x=0,y=0,w=300,h=150,vx,vy,vw,vh,par=1;
        char set_whv = 0;
        for (int i = 0; i < n_attr; i++){
          char keybuf[S2P_MAX_AK_LEN+1] = {0};
          fseek(fd, attrs[i].k, SEEK_SET);
          fread(keybuf,1,S2P__MIN(attrs[i].k_len,S2P_MAX_AK_LEN),fd);
          if (!strcmp(keybuf,"x")){
            fseek(fd, attrs[i].v, SEEK_SET);
            fscanf(fd, "%f", &x);
          }else if (!strcmp(keybuf,"y")){
            fseek(fd, attrs[i].v, SEEK_SET);
            fscanf(fd, "%f", &y);
          }else if (!strcmp(keybuf,"width")){
            set_whv |= 4;
            fseek(fd, attrs[i].v, SEEK_SET);
            fscanf(fd, "%f", &w);
          }else if (!strcmp(keybuf,"height")){
            set_whv |= 2;
            fseek(fd, attrs[i].v, SEEK_SET);
            fscanf(fd, "%f", &h);
            
          }else if (!strcmp(keybuf,"viewBox")){
            set_whv |= 1;
            fseek(fd, attrs[i].v, SEEK_SET);
            fscanf(fd, "%f", &vx);
            s2p__d_next(fd);
            fscanf(fd, "%f", &vy);
            s2p__d_next(fd);
            fscanf(fd, "%f", &vw);
            s2p__d_next(fd);
            fscanf(fd, "%f", &vh);
          }else if (!strncmp(keybuf,"prese",5)){
            fseek(fd, attrs[i].v, SEEK_SET);
            s2p__d_next(fd);
            if (fgetc(fd) == 'n'){
              par = 0;
            }
          }
        }
        

        if (set_whv == 6 || set_whv == 0){
          vx = 0; vy = 0; vw = w; vh = h;
        }else if (set_whv == 1){
          w = vw; h = vh;
        }else if (set_whv == 5){
          h = w * vh / vw;
        }else if (set_whv == 3){
          w = h * vw / vh;
        }else if (set_whv == 2){
          w = h;
          vx = 0; vy = 0; vw = w; vh = h;
        }else if (set_whv == 4){
          h = w;
          vx = 0; vy = 0; vw = w; vh = h;
        }

        if (!s2p__mati){
          S2P_SETDIM(w,h);
        }
        
        s2p__mati ++;
        memcpy(s2p__mats+s2p__mati,s2p__mats+(s2p__mati-1),9*sizeof(float));

        float sx = w/vw;
        float sy = h/vh;
        float m[9] = {sx,0,-vx*sx+x, 0,sy,-vy*sy+y, 0,0,1};
        if (par){
          float s = fmin(sx,sy);
          float px = (w-vw*s)/2.0;
          float py = (h-vh*s)/2.0;
          m[0] = s; m[2] = -vx*s+x+px;
          m[4] = s; m[5] = -vy*s+y+py;
        }

        s2p__matmul(s2p__mats[s2p__mati],m,s2p__mats[s2p__mati]);
        s2p__get_transf(fd,attrs,n_attr,s2p__mats[s2p__mati]);

      }else if (closing == S2P_TAG_CLOSE){
        s2p__mati--;
      }

    }else if (!strcmp(tagbuf,"g")){
      
      if (closing == S2P_TAG_OPEN){
        s2p__mati ++;
        memcpy(s2p__mats+s2p__mati,s2p__mats+(s2p__mati-1),9*sizeof(float));
        s2p__get_transf(fd,attrs,n_attr,s2p__mats[s2p__mati]);

        if (s2p_skip_nostroke){
          stroke_stack <<= 1;
          stroke_stack |= s2p__get_stroke(fd,attrs,n_attr);
        }

      }else if (closing == S2P_TAG_CLOSE){
        s2p__mati--;
        stroke_stack >>= 1;
      }
    }else if (closing != S2P_TAG_CLOSE){

      if (s2p_skip_nostroke){
        if (!stroke_stack && !s2p__get_stroke(fd,attrs,n_attr)){
          goto lbl_done_tag;
        } 
      }

      s2p__mati ++;
      memcpy(s2p__mats+s2p__mati,s2p__mats+(s2p__mati-1),9*sizeof(float));
      s2p__get_transf(fd,attrs,n_attr,s2p__mats[s2p__mati]);

      if (!strcmp(tagbuf,"line")){

        s2p__parse_line(fd,attrs,n_attr);

      }else if (!strncmp(tagbuf,"poly",4)){

        s2p__parse_poly(fd,attrs,n_attr,tagbuf[4]=='g');

      }else if (!strcmp(tagbuf,"path")){
        
        s2p__parse_path(fd,attrs,n_attr);

      }else if (!strcmp(tagbuf,"rect")){
        
        s2p__parse_rect(fd,attrs,n_attr);

      }else if (!strcmp(tagbuf,"ellipse") || !strcmp(tagbuf,"circle")){

        s2p__parse_ellipse(fd,attrs,n_attr);

      }
      s2p__mati --;

    }
lbl_done_tag:;
    
    fseek(fd,pos,SEEK_SET);

lbl_reset:;
    tag_len = 0;
    tag = 0;
    n_attr = 0;
    memset(attrs,0,sizeof(attrs));
    state = S2P__ST_OUT;
    closing = 0;
  }
}

void s2p_parse_from_file(char* path){
  FILE* fd = fopen(path,"r");
  s2p_parse(fd);
  fclose(fd);
}

void s2p_parse_from_str(char* s, int n){
  if (n == -1)  n = strlen(s);
  FILE* fd = fmemopen(s,n,"r");
  s2p_parse(fd);
  fclose(fd);
}

#endif
