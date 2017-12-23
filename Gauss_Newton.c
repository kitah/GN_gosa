#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <nd_malloc.h>
#include <imageio.h>
#include "bpf.h"

#define SGMC 0.85 //0.72
#define PSF_HS 11
#define BPF_HS  3
#define DELTA 0.001
#define BPF_N 0
//#define h(x,y) hh[(x)+(L-1)/2][(y)+(L-1)/2]

void gaussNewtonMethod(double ***fn, int nz, int BS, int center_x, int center_y, double *dk){
  int i, j, bn, bm, count, NZ = nz;
  FILE *fp;
  char flnm[256];
  double sgmc = SGMC;
  double tmp_e1, tmp_e2, tmp_e3;
  double tmp_a00_01, tmp_a00_12, tmp_a00_02, tmp_a01_01, tmp_a01_12, tmp_a01_02, tmp_a11_01, tmp_a11_12, tmp_a11_02;
  double a00, a01, a10, a11, b0, b1;
  double **f0p1, **f1p0, **f1p2, **f2p1, **f0p2, **f2p0;
  double **f0pd1, **f1pd0, **f1pd2, **f2pd1, **f0pd2, **f2pd0;
  double **f0pk1, **f1pk0, **f1pk2, **f2pk1, **f0pk2, **f2pk0;
  double **p0, **p1, **p2;
  double e[NZ], mz[NZ], d0k0[2], e1[NZ];
  double k0, d0, deltad, deltak, dprev, kprev, eprev, enow;
  
  void psf(double k, int z, double d, double sgmc, int psf_hs, double **p);
  void dpsfdd(double k, int z, double d, double sgmc, int psf_hs, double **p);
  void dpsfdk(double k, int z, double d, double sgmc, int psf_hs, double **p);
  void make_fnpn(int bs, int center_x, int center_y, double **fn, double **pn, double **fnpn);
  void make_fnpdn(int bs, int center_x, int center_y, double **fn, double **fnpdn, double k, int z, double d, double sgmc);
  void make_fnpkn(int bs, int center_x, int center_y, double **fn, double **fnpkn, double k, int z, double d, double sgmc);  
  void num_diff_method(int bs, double ** fnpndelta, double **fnpn, double **f);
  void filter(double ***fn, int nz, int BS, int center_x, int center_y, double *mz);
  void calculate_d_k(double *mz, int nz, double *d0k0);

  //-----------------------------------------------------------------------------
  // malloc
  p0   = malloc_double_2d(2 * PSF_HS + 1, PSF_HS,  2 * PSF_HS + 1, PSF_HS);
  p1   = malloc_double_2d(2 * PSF_HS + 1, PSF_HS,  2 * PSF_HS + 1, PSF_HS);
  p2   = malloc_double_2d(2 * PSF_HS + 1, PSF_HS,  2 * PSF_HS + 1, PSF_HS);  

  f0p1 = malloc_double_2d(BS, 0, BS, 0);
  f1p0 = malloc_double_2d(BS, 0, BS, 0);
  f1p2 = malloc_double_2d(BS, 0, BS, 0);
  f2p1 = malloc_double_2d(BS, 0, BS, 0);
  f0p2 = malloc_double_2d(BS, 0, BS, 0);
  f2p0 = malloc_double_2d(BS, 0, BS, 0);

  f0pd1 = malloc_double_2d(BS, 0, BS, 0);
  f1pd0 = malloc_double_2d(BS, 0, BS, 0);
  f1pd2 = malloc_double_2d(BS, 0, BS, 0);
  f2pd1 = malloc_double_2d(BS, 0, BS, 0);
  f0pd2 = malloc_double_2d(BS, 0, BS, 0);
  f2pd0 = malloc_double_2d(BS, 0, BS, 0);

  f0pk1 = malloc_double_2d(BS, 0, BS, 0);
  f1pk0 = malloc_double_2d(BS, 0, BS, 0);
  f1pk2 = malloc_double_2d(BS, 0, BS, 0);
  f2pk1 = malloc_double_2d(BS, 0, BS, 0);
  f0pk2 = malloc_double_2d(BS, 0, BS, 0);
  f2pk0 = malloc_double_2d(BS, 0, BS, 0);

  filter(fn, NZ, BS, center_x, center_y, mz);
  calculate_d_k(mz, NZ, d0k0);

  //printf("\nest_d  est_k for SFF: \n");
  //printf("[%f %f]\n\n", d0k0[0], d0k0[1]);
  
  d0 = d0k0[0];
  k0 = d0k0[1];
  
  count = 0;
  //d = 0.584408;
  //k = 0.5;
  /*
  double tmp1 = 0.0, tmp2 = 0.0, tmp3 = 0.0;
  for(j = 0 ; j < BS ; j++){
    for(i = 0 ; i < BS ; i++){
      tmp1 = (f0p1[i][j] - f1p0[i][j]);
      tmp2 = (f1p2[i][j] - f2p1[i][j]);
      tmp3 = (f0p2[i][j] - f2p0[i][j]);
      
      e1[0] += tmp1 * tmp1;
      e1[1] += tmp2 * tmp2;
      e1[2] += tmp3 * tmp3;
    }
    }*/
  
  eprev = enow = 0.0;
  
  // do-whileのループで対象のブロックのdepth"d"とレンズ定数"k"を推定する
  do{
    eprev = enow;
    psf(k0, 0, d0, sgmc, PSF_HS, p0);
    psf(k0, 1, d0, sgmc, PSF_HS, p1);
    psf(k0, 2, d0, sgmc, PSF_HS, p2);
    
    make_fnpn(BS, center_x, center_y, fn[0], p1, f0p1);
    make_fnpn(BS, center_x, center_y, fn[1], p0, f1p0);
    make_fnpn(BS, center_x, center_y, fn[1], p2, f1p2);
    make_fnpn(BS, center_x, center_y, fn[2], p1, f2p1);
    make_fnpn(BS, center_x, center_y, fn[0], p2, f0p2);
    make_fnpn(BS, center_x, center_y, fn[2], p0, f2p0);    
    
    make_fnpdn(BS, center_x, center_y, fn[0], f0pd1, k0, 1, d0, SGMC);
    make_fnpdn(BS, center_x, center_y, fn[1], f1pd0, k0, 0, d0, SGMC);
    make_fnpdn(BS, center_x, center_y, fn[1], f1pd2, k0, 2, d0, SGMC);
    make_fnpdn(BS, center_x, center_y, fn[2], f2pd1, k0, 1, d0, SGMC);
    make_fnpdn(BS, center_x, center_y, fn[0], f0pd2, k0, 2, d0, SGMC);
    make_fnpdn(BS, center_x, center_y, fn[2], f2pd0, k0, 0, d0, SGMC);
        
    make_fnpkn(BS, center_x, center_y, fn[0], f0pk1, k0, 1, d0, SGMC);
    make_fnpkn(BS, center_x, center_y, fn[1], f1pk0, k0, 0, d0, SGMC);
    make_fnpkn(BS, center_x, center_y, fn[1], f1pk2, k0, 2, d0, SGMC);
    make_fnpkn(BS, center_x, center_y, fn[2], f2pk1, k0, 1, d0, SGMC);
    make_fnpkn(BS, center_x, center_y, fn[0], f0pk2, k0, 2, d0, SGMC);
    make_fnpkn(BS, center_x, center_y, fn[2], f2pk0, k0, 0, d0, SGMC);
    
    tmp_e1 = tmp_e2 = tmp_e3 = tmp_a00_01 = tmp_a00_12 = tmp_a00_02 = 0.0;
    tmp_a01_01 = tmp_a01_12 = tmp_a01_02 = 0.0, tmp_a11_01 = tmp_a11_12 = tmp_a11_02 = 0.0;
    a00 = a01 = a10 = a11 = b0 = b1 = 0.0;
    e[0] = e[1] = e[2] = 0.0;
    for(j = 0 ; j < BS ; j++) {
      for(i = 0 ; i < BS ; i++) {
	tmp_e1 = (f0p1[i][j] - f1p0[i][j]);
	tmp_e2 = (f1p2[i][j] - f2p1[i][j]);
	tmp_e3 = (f0p2[i][j] - f2p0[i][j]);

	e[0] += tmp_e1 * tmp_e1;
	e[1] += tmp_e2 * tmp_e2;
	e[2] += tmp_e3 * tmp_e3;
	
	tmp_a00_01 += (f0pd1[i][j] - f1pd0[i][j]) * (f0pd1[i][j] - f1pd0[i][j]);
	tmp_a00_12 += (f1pd2[i][j] - f2pd1[i][j]) * (f1pd2[i][j] - f2pd1[i][j]);  
	tmp_a00_02 += (f0pd2[i][j] - f2pd0[i][j]) * (f0pd2[i][j] - f2pd0[i][j]);

	tmp_a01_01 += (f0pd1[i][j] - f1pd0[i][j]) * (f0pk1[i][j] - f1pk0[i][j]);
	tmp_a01_12 += (f1pd2[i][j] - f2pd1[i][j]) * (f1pk2[i][j] - f2pk1[i][j]);  
	tmp_a01_02 += (f0pd2[i][j] - f2pd0[i][j]) * (f0pk2[i][j] - f2pk0[i][j]);
	
	tmp_a11_01 += (f0pk1[i][j] - f1pk0[i][j]) * (f0pk1[i][j] - f1pk0[i][j]);
	tmp_a11_12 += (f1pk2[i][j] - f2pk1[i][j]) * (f1pk2[i][j] - f2pk1[i][j]);  
	tmp_a11_02 += (f0pk2[i][j] - f2pk0[i][j]) * (f0pk2[i][j] - f2pk0[i][j]);
	
	b0 += tmp_e1 * (f0pd1[i][j] - f1pd0[i][j]) + tmp_e2 * (f1pd2[i][j] - f2pd1[i][j]) + tmp_e3 * (f0pd2[i][j] - f2pd0[i][j]);
	b1 += tmp_e1 * (f0pk1[i][j] - f1pk0[i][j]) + tmp_e2 * (f1pk2[i][j] - f2pk1[i][j]) + tmp_e3 * (f0pk2[i][j] - f2pk0[i][j]);

	/*if((j = 0) || (j == (BS-1))){
	  e[0] = e[0] * 0.2;
	  e[1] = e[0] * 0.2;
	  e[2] = e[0] * 0.2;

	  tmp_a00_01 = tmp_a00_01 * 0.2;
	  tmp_a00_12 = tmp_a00_12 * 0.2;
	  tmp_a00_02 = tmp_a00_02 * 0.2;

	  tmp_a01_01 = tmp_a01_01 * 0.2;
	  tmp_a01_12 = tmp_a01_12 * 0.2;
	  tmp_a01_02 = tmp_a01_02 * 0.2;

	  tmp_a11_01 = tmp_a11_01 * 0.2;
	  tmp_a11_12 = tmp_a11_12 * 0.2;
	  tmp_a11_02 = tmp_a11_02 * 0.2;
	  }*/
      }
    }

    //printf("E[%d] = %f\n", count, e[0] + e[1] + e[2]);
    enow = e[0] + e[1] + e[2];

    //printf("eprev = %f enow = %f\n", eprev, enow);
    a00 = tmp_a00_01 + tmp_a00_12 + tmp_a00_02;
    a01 = tmp_a01_01 + tmp_a01_12 + tmp_a01_02;
    a10 = a01;
    a11 = tmp_a11_01 + tmp_a11_12 + tmp_a11_02;
    b0 = (-1)*b0; //b0
    b1 = (-1)*b1; //b1
    
    // deltadの計算をして、d0を更新
    deltad = (a11 * b0 - a01 * b1) / (a00 * a11 - a01 * a10);
    dprev = d0;
    d0 = dprev + deltad;
    
    //printf("deltad: %f\n", deltad);
    //printf("d0: %f\n", d0);
    
    // deltakの計算をして、k0を更新
    deltak = (a00 * b1 - a01 * b0) / (a00 * a11 - a01 * a10);
    kprev = k0;
    k0 = kprev + deltak;
    
    //printf("deltak: %f\n", deltak);
    //printf("k0: %f\n", k0);
    
    count++;
  }while(fabs(eprev - enow) > 0.01);
//} while(fabs(deltad) > 0.01);
//printf("\n");
//  } while( fabs( (dprev - model_d) - (d0 - model_d) ) / fabs(dprev - model_d) > 0.01);
/*
  printf("--------------------------------------------\n");
  printf("loop = %d\n", count);
  printf("d = %f, k = %f\n", d0, k0);
  printf("dprev = %f\n", dprev);
  printf("deltad = %f, deltak = %f\n", deltad, deltak);
  printf("gosa = %f\n", d0 - d);
*/
//free_double_2d(p, 2 * PSF_HS + 1, PSF_HS,  2 * PSF_HS + 1, PSF_HS);
//rintf("estE = %f\n", enow);
//printf("trueE = %f\n", e1[0]+e1[1]+e1[2]);

  dk[0] = d0;
  dk[1] = k0;
}

/*---------------------------------------------------------------------*/
void psf(double k, int z, double d, double sgmc, int psf_hs, double **p)
{
  int i, j;
  double c, K, K2;
 
  K = k * ((double)z - d);
  K2 = K * K + sgmc * sgmc;

  c = 0.0;
  for(j = -PSF_HS ; j <= PSF_HS ; j++) {
    for(i = -PSF_HS ; i <= PSF_HS ; i++) {
      p[i][j] = exp( -(double)(i * i + j * j) / 2.0 / K2);
      c += p[i][j];
    }
  }

  c = 2.0 * M_PI * (K * K + sgmc * sgmc);
  for(j = -PSF_HS ; j <= PSF_HS ; j++) {
    for(i = -PSF_HS ; i <= PSF_HS ; i++) {
      p[i][j] /= c;
    }
  }
}

/*---------------------------------------------------------------------*/
void dpsfdd(double k, int z, double d, double sgmc, int psf_hs, double **dpdd)
{
  int i, j;
  double K, K2;
  double P0, P1, P2, P3, i2j2;


  K = k * ((double)z - d);
  K2 = K * K + sgmc * sgmc;

  P0 = 1.0 / (M_PI * K2 * K2);
  P3 = -k * k * (z - d);
  
  for(j = -PSF_HS ; j <= PSF_HS ; j++) {
    for(i = -PSF_HS ; i <= PSF_HS ; i++) {
      i2j2 = (double)(i * i + j * j);
      P1 = exp(-i2j2 / 2.0 / K2);
      P2 = i2j2 / 2.0 / K2 - 1.0;
      dpdd[i][j] = P0 * P1 * P2 * P3;
    }
  }
}

/*---------------------------------------------------------------------*/
void dpsfdk(double k, int z, double d, double sgmc, int psf_hs, double **dpdk)
{
  int i, j;
  double K, K2;
  double P0, P1, P2, P3, i2j2;


  K = k * ((double)z - d);
  K2 = K * K + sgmc * sgmc;

  P0 = 1.0 / (M_PI * K2 * K2); 
  P3 = k * (z - d) * (z - d);
  
  for(j = -PSF_HS ; j <= PSF_HS ; j++) {
    for(i = -PSF_HS ; i <= PSF_HS ; i++) {
      i2j2 = (double)(i * i + j * j);
      P1 = exp(-i2j2 / 2.0 / K2);
      P2 = i2j2 / 2.0 / K2 - 1.0;
      dpdk[i][j] = P0 * P1 * P2 * P3;
    }
  }
}

/*---------------------------------------------------------------------*/
void make_fnpn(int bs, int center_x, int center_y, double **fn, double **pn, double **fnpn)
{
  int i, j, n, m;
  int hbs = (bs - 1) / 2; // 7
  
  for(j = 0 ; j < bs ; j++) {
    for(i = 0 ; i < bs ; i++) {
      fnpn[i][j] = 0.0;
      for(m = -PSF_HS ; m <= PSF_HS ; m++) {
	for(n = -PSF_HS ; n <= PSF_HS ; n++) {
	  fnpn[i][j] += fn[center_x - hbs + i - n][center_y - hbs + j - m] * pn[n][m];
	}
      }
    }
  }
  
}

/*---------------------------------------------------------------------*/
void make_fnpdn(int bs, int center_x, int center_y, double **fn, double **fnpdn, double k, int z, double d, double sgmc)
{
  int i, j, n, m;
  double **dpdd;
  int hbs = (bs - 1) / 2;
  
  dpdd = malloc_double_2d(2 * PSF_HS + 1, PSF_HS,  2 * PSF_HS + 1, PSF_HS);
  dpsfdd(k, z, d, sgmc, PSF_HS, dpdd);
   
  for(j = 0 ; j < bs ; j++) {
    for(i = 0 ; i < bs ; i++) {
      fnpdn[i][j] = 0.0;
      for(m = -PSF_HS ; m <= PSF_HS ; m++) {
	for(n = -PSF_HS ; n <= PSF_HS ; n++) {
	  fnpdn[i][j] += fn[center_x - hbs + i - n][center_y - hbs + j - m] * dpdd[n][m];
	}
      }
    }
  }
  
}

/*---------------------------------------------------------------------*/
void make_fnpkn(int bs, int center_x, int center_y, double **fn, double **fnpkn, double k, int z, double d, double sgmc)
{
  int i, j, n, m;
  double **dpdk;
  int hbs = (bs - 1) / 2;
  
  dpdk = malloc_double_2d(2 * PSF_HS + 1, PSF_HS,  2 * PSF_HS + 1, PSF_HS);
  dpsfdk(k, z, d, sgmc, PSF_HS, dpdk);
  
  for(j = 0 ; j < bs ; j++) {
    for(i = 0 ; i < bs ; i++) {
      fnpkn[i][j] = 0.0;
      for(m = -PSF_HS ; m <= PSF_HS ; m++) {
	for(n = -PSF_HS ; n <= PSF_HS ; n++) {	  
	  fnpkn[i][j] += fn[center_x  - hbs + i - n][center_y - hbs + j - m] * dpdk[n][m];	  
	}
      }
    }
  }
}

/*---------------------------------------------------------------------*/
void num_diff_method(int bs, double ** fnpndelta, double **fnpn, double **f){
  int i, j;

  for(j = 0 ; j < bs ; j++){
    for(i = 0 ; i < bs ; i++){
      f[i][j] = (fnpndelta[i][j] - fnpn[i][j]) / DELTA;
    }
  }
}

/*---------------------------------------------------------------------*/
void filter(double ***fn, int nz, int BS, int center_x, int center_y, double *mz){
   //BPF
  /*double hh[9][9] =
  {
    {
       -0.001014896851, 
        0.000655810401, 
       -0.003138266557, 
       -0.014109352684, 
       -0.018787171444, 
       -0.014109352684, 
       -0.003138266557, 
        0.000655810401, 
       -0.001014896851
    }, 
    {
        0.000655810401, 
       -0.009738774136, 
       -0.031145719861, 
       -0.034508428742, 
       -0.027133646128, 
       -0.034508428742, 
       -0.031145719861, 
       -0.009738774136, 
        0.000655810401
    }, 
    {
       -0.003138266557, 
       -0.031145719861, 
       -0.030032998719, 
        0.006440072228, 
        0.017091672591, 
        0.006440072228, 
       -0.030032998719, 
       -0.031145719861, 
       -0.003138266557
    }, 
    {
       -0.014109352684, 
       -0.034508428742, 
        0.006440072228, 
        0.078528220659, 
        0.103452025961, 
        0.078528220659, 
        0.006440072228, 
       -0.034508428742, 
       -0.014109352684
    }, 
    {
       -0.018787171444, 
       -0.027133646128, 
        0.017091672591, 
        0.103452025961, 
        0.156989353988, 
        0.103452025961, 
        0.017091672591, 
       -0.027133646128, 
       -0.018787171444
    }, 
    {
       -0.014109352684, 
       -0.034508428742, 
        0.006440072228, 
        0.078528220659, 
        0.103452025961, 
        0.078528220659, 
        0.006440072228, 
       -0.034508428742, 
       -0.014109352684
    }, 
    {
       -0.003138266557, 
       -0.031145719861, 
       -0.030032998719, 
        0.006440072228, 
        0.017091672591, 
        0.006440072228, 
       -0.030032998719, 
       -0.031145719861, 
       -0.003138266557
    }, 
    {
        0.000655810401, 
       -0.009738774136, 
       -0.031145719861, 
       -0.034508428742, 
       -0.027133646128, 
       -0.034508428742, 
       -0.031145719861, 
       -0.009738774136, 
        0.000655810401
    }, 
    {
       -0.001014896851, 
        0.000655810401, 
       -0.003138266557, 
       -0.014109352684, 
       -0.018787171444, 
       -0.014109352684, 
       -0.003138266557, 
        0.000655810401, 
       -0.001014896851
    } 
    };
  */
  
  int L = 7, L1 = (L-1)/2;
   int zn, i, j, ii, jj, xsize, ysize, hbs, NZ = nz;
   double **g, f_tmp, **gz, ***in_img;
   FILE *fp;
   char flnm[256];
   int bn, bm;
   
   xsize = BS;
   ysize = BS;
   hbs = (BS - 1) / 2;
   
   //malloc
   in_img = malloc_double_3d(NZ, 0, xsize, 0, ysize, 0);
   //g = malloc_double_2d(L1*2 + xsize, L1, L1*2 + ysize, L1);
   g = malloc_double_2d(PSF_HS*2 + xsize, PSF_HS, PSF_HS*2 + ysize, PSF_HS);
   gz = malloc_double_2d(xsize, 0, ysize, 0);
   
   //init
   for(zn = 0 ; zn < NZ ; zn++){
     mz[zn] = 0.0;
   }

   //1ブロック分の画像生成
   for(zn = 0 ; zn < NZ ; zn++){
     for(j = 0 ; j < ysize ; j++){
       for(i = 0 ; i < xsize ; i++){
	 in_img[zn][i][j] = fn[zn][center_x - hbs + i][center_y - hbs + j];
       }
     }   

     //折り返し
     //copy
     for(j = 0 ; j < ysize ; j++){
       for(i = 0 ; i < xsize ; i++){
	 g[i][j] = in_img[zn][i][j];
       }
     }
     
     for(j = 0 ; j < ysize ; j++){
       for(i = 1 ; i <= PSF_HS ; i++){
	 g[-i][j] = g[i][j];
	 g[xsize - 1 + i][j] = g[xsize - 1 - i][j];
       }
     }
     
     for(i = -PSF_HS ; i < xsize + PSF_HS ; i++){
       for(j = 1 ; j <= PSF_HS ; j++){
	 g[i][-j] = g[i][j];
	 g[i][ysize - 1 + j] = g[i][ysize - 1 - j];
       }
     }
     
     //二次元のフィルタリング処理
     for(j = 0 ; j < ysize ; j++){
       for(i = 0 ; i < xsize ; i++){
	 f_tmp = 0.0;
	 for(jj = -L1 ; jj <= L1 ; jj++){
	   for(ii = -L1 ; ii <= L1 ; ii++){
	     f_tmp += g[i - ii][j - jj] * bpf[BPF_N][ii + BPF_HS][jj + BPF_HS];
	     //f_tmp += g[i - ii][j - jj] * h(ii, jj);
	   }
	 }

	 gz[i][j] = f_tmp;
	 if(f_tmp >= 256) gz[i][j] = 255.0;
	 else if (f_tmp > 0) gz[i][j] = 0.0;
	 
       }
     }
     /*
     sprintf(flnm, "./107_%d.dat", zn);
     fp = fopen(flnm, "w");
     for(bn = 0 ; bn < BS ; bn++){
       for(bm = 0 ; bm < BS ; bm++){
	 fprintf(fp, "%4d %4d  %f\n", bn, bm, gz[bn][bm]);
       }
       fprintf(fp, "\n");
       }*/
     
     //合焦評価値
     for(j = 0 ; j < ysize ; j++){
       for(i = 0 ; i < xsize ; i++){	 
	 mz[zn] += (gz[i][j]) * (gz[i][j]);
       }
     }
     mz[zn] = log(mz[zn]);
     //printf("mz[%d] = %f\n", zn, mz[zn]);
   }//zn

  
 }

/*---------------------------------------------------------------------*/
//p, q = bs. do this in each block

void calculate_d_k(double *mz, int nz, double *d0k0){
  int x1, x2, x3, max_n = 0, n, NZ = nz;
  double max = 0, y1, y2, y3, a, b;//c;
  double a1, a2, a3, b1, b2, b3, c1, c2, c3, d1, d2, d3;
  
  for(n = 0 ; n < NZ ; n++){
    if(mz[n] > max) {
      max = mz[n];
      max_n = n;
    }
  }
  
  if(max_n == 0){
    x1 = max_n;
    x2 = max_n + 1;
    x3 = max_n + 2;
  }else if(max_n == NZ){
    x1 = max_n - 2;
    x2 = max_n - 1;
    x3 = max_n;
  }else{
    x1 = max_n - 1;
    x2 = max_n;
    x3 = max_n + 1;
  }

  //計算する為には式は3つで十分
  y1 = mz[x1];
  y2 = mz[x2];
  y3 = mz[x3];

  //y = ax2+bx+c に(x,y)を代入
  a1 = (double)x1 * (double)x1;
  b1 = (double)x1;
  c1 = 1.0;
  a2 = (double)x2 * (double)x2;
  b2 = (double)x2;
  c2 = 1.0;
  a3 = (double)x3 * (double)x3;
  b3 = (double)x3;
  c3 = 1.0;
  d1 = y1;
  d2 = y2;
  d3 = y3;
  
  a = (d1 * (b2*c3 - c2*b3) - b1 * (d2*c3 - c2*d3) + c1 * (d2*b3 - b2*d3)) / (a1 * (b2*c3 - c2*b3) - b1 * (a2*c3 - c2*a3) + c1 * (a2*b3 - b2*a3));
  b = (a1 * (d2*c3 - c2*d3) - d1 * (a2*c3 - c2*a3) + c1 * (a2*d3 - d2*a3)) / (a1 * (b2*c3 - c2*b3) - b1 * (a2*c3 - c2*a3) + c1 * (a2*b3 - b2*a3));
  //c = (a1 * (b2*d3 - d2*b3) - b1 * (a2*d3 - d2*a3) + d1 * (a2*b3 - b2*a3)) / (a1 * (b2*c3 - c2*b3) - b1 * (a2*c3 - c2*a3) + c1 * (a2*b3 - b2*a3));

  d0k0[0] = (-1.0 * b) / (2.0 * a); //d0
  d0k0[1] = 2.0 * sqrt(fabs(a)) / M_PI;     //k0
}
