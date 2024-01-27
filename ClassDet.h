#ifndef ClassDet_h
#define ClassDet_h

#include <cstdio>

#define MAXMULTI 1000

class Det{
public:
  Det(): multi(0) {Clear(); }

  unsigned short       multi; // max 65535
  unsigned short       id[MAXMULTI];
  unsigned short       ch[MAXMULTI];
  unsigned short       e[MAXMULTI];
  unsigned long long   t[MAXMULTI];

  unsigned short  index[MAXMULTI]; // id * 12 + ch;

  void Clear(){
    multi = 0;
    for( int i = 0; i < MAXMULTI; i++){
      id[i] = 0;
      ch[i] = 0;
       e[i] = 0;
       t[i] = 0;
       index[i] = 0;
    }
  }

  void Print(){
    printf("=============================== multi : %u\n", multi);
    for( int i = 0; i < multi; i++) {
      printf(" %3d | %2d-%-2d(%5d) %5u %15llu \n", i, id[i], ch[i], index[i], e[i], t[i]);
    }
  }

  void SetDetDimension(unsigned short maxID, unsigned maxCh){
    nID = maxID;
    nCh = maxCh;
  }

  void CalIndex(){
    for( int i = 0; i < multi; i++){
      index[i] = id[i] * nCh + ch[i] ;
    }
  }

private:

  unsigned short nID;
  unsigned short nCh;

};


#endif 