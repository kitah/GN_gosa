#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <imageio.h>
#include <math.h>
#include <sys/stat.h>
#include <time.h>
//#include "lpf.h"
//#include "bpf.h"
//#include "func.h"

/**** Default値の設定 ****/
#define BLOCK_SIZE 13
#define PITCH     13
#define RR         1

int main(int argc, char *argv[]){
  /*FILE *fp;
  IMAGE img_go;
  char flnm[256];
  int zn, n, m, k, l, bn, bm, b0[2], b1[2], bc[2];
  int nr, nw, hbsn = (BSNN-1) / 2, hbsm = (BSMM-1) / 2;
  double ***go, ***g_sff, ***y, **dy01, **dy12, **dy20, ***p, ***go_t;
  double **d, **d_model, **d_est, **k_est, **d_sff, **k_sff, ***E;
  double **dm_est, **dp_est, ***ydm, ***yd, ***ydp, **dp_model, **dm_model;
  double data[(NHP+NX+NHP) * (NHP+NY+NHP)];
  double dmin, kmin, Emin, dtest, ktest, z, sgm, d_tmp;
  double EE[3], sum_EE = 0.0, tmp, RMSE;
  time_t t1, t2;
  t1 = time(NULL);

  double alpha = 0.0, dmax[NNB], d_min[NNB];
  */
  int pitch = PITCH, block_size = BLOCK_SIZE, r = RR;
  char c;
  void usage(char *);
  /*
#if SRAND_FLG == 1
  srandom(time(NULL));
#endif
  */
  //---------------------------------------------------------------------------
  // 引数の処理
  if (argc < 2) usage(argv[0]);

  while((c = getopt(argc, argv, "b:p:r:")) != -1){
    switch(c) {
    case 'b':
      block_size = atoi(optarg);
      break;
    case 'p':
      pitch = atoi(optarg);
      break;
    case 'r':
      r = atoi(optarg);
      break;
    default:
      usage(argv[0]);
    }
  }
  //  argc -= optind - 1;
  //  argv += (optind - 1);
  //  if(argc > 2) usage(argv[0]);

  printf("BS = %d\n", block_size);
  printf("pitch = %d\n", pitch);
  printf("R = %d\n", r);
}


/*--------------------------------------------------------------------------*/
void usage(char *com)
{
  fprintf(stderr, "\n");
  fprintf(stderr, "  Usage : %s [-bpr]\n\n", com);
  fprintf(stderr, "          -b <value> : BlockSize Odd   (int)\n");
  fprintf(stderr, "          -p <value> : Block Pitch     (int)\n");
  fprintf(stderr, "          -r <value> : Exclusion Block (int)\n");
  fprintf(stderr, "\n");

  exit(1);
}


  /*
  //---------------------------------------------------------------------------
  // 領域確保 
  // go[zn][NX][NY]
  if((go = (double ***)malloc(sizeof(double **) * NZ)) == NULL){
    fprintf(stderr, "\nError : malloc\n\n");exit(1);
  }
  for(zn = 0 ; zn < NZ ; zn++){
    if((go[zn] = (double **)malloc(sizeof(double *) * (NX+2*NHP))) == NULL){
      fprintf(stderr, "\nError : malloc\n\n");exit(1);    
    }
    go[zn] += NHP;
    for(n = -NHP ; n < NX + NHP ; n++){
      if((go[zn][n] = (double *)malloc(sizeof(double ) * (NY+2*NHP))) == NULL){
        fprintf(stderr, "\nError : malloc\n\n");exit(1);    
      }
      go[zn][n] += NHP;
    }
  }

  // go_t[zn][15*20][15*13]
  if((go_t = (double ***)malloc(sizeof(double **) * NZ)) == NULL){
    fprintf(stderr, "\nError : malloc\n\n");exit(1);
  }
  for(zn = 0 ; zn < NZ ; zn++){
    if((go_t[zn] = (double **)malloc(sizeof(double *) * (NNX+2*NHP))) == NULL){
      fprintf(stderr, "\nError : malloc\n\n");exit(1);    
    }
    go_t[zn] += NHP;
    for(n = -NHP ; n < NNX + NHP ; n++){
      if((go_t[zn][n] = (double *)malloc(sizeof(double ) * (NNY+2*NHP))) == NULL){
        fprintf(stderr, "\nError : malloc\n\n");exit(1);    
      }
      go_t[zn][n] += NHP;
    }
  }

  // g_sff[zn][NX][NY]
  if((g_sff = (double ***)malloc(sizeof(double **) * NZ)) == NULL){
    fprintf(stderr, "\nError : malloc\n\n");exit(1);
  }
  for(zn = 0 ; zn < NZ ; zn++){
    if((g_sff[zn] = (double **)malloc(sizeof(double *) * NX)) == NULL){
      fprintf(stderr, "\nError : malloc\n\n");exit(1);    
    }
    for(n = 0 ; n < NX ; n++){
      if((g_sff[zn][n] = (double *)malloc(sizeof(double ) * NY)) == NULL){
        fprintf(stderr, "\nError : malloc\n\n");exit(1);    
      }
    }
  }

  // y[2*NZ][BSNN+4*NHP][BSMM+4*NHP]
  if((y = (double ***)malloc(sizeof(double **) * NZ*2)) == NULL){
    fprintf(stderr, "\nError : malloc\n\n");exit(1);
  }
  for(zn = 0 ; zn < 2*NZ ; zn++){
    if((y[zn] = (double **)malloc(sizeof(double *) * (BSNN+4*NHP))) == NULL){
      fprintf(stderr, "\nError : malloc\n\n");exit(1);    
    }
    y[zn] += 2*NHP;
    for(n = -2*NHP ; n < BSNN + 2*NHP ; n++){
      if((y[zn][n] = (double *)malloc(sizeof(double ) * (BSMM+4*NHP))) == NULL){
        fprintf(stderr, "\nError : malloc\n\n");exit(1);    
      }
      y[zn][n] += 2*NHP;
    }
  }

  // ydm, yd, ydp[][][]
  if((ydm = (double ***)malloc(sizeof(double **) * NZ * 2)) == NULL ||
     (yd = (double ***)malloc(sizeof(double **) * NZ * 2)) == NULL ||
      (ydp = (double ***)malloc(sizeof(double **) * NZ * 2)) == NULL){
    fprintf(stderr, "\nError : malloc\n\n");exit(1);
  }
  for(zn = 0 ; zn < 2*NZ ; zn++){
    if((ydm[zn] = (double **)malloc(sizeof(double *) * (INBS + 4*NHP))) == NULL ||
       (yd[zn] = (double **)malloc(sizeof(double *) * (INBS + EXBS*2 + 4*NHP))) == NULL ||
       (ydp[zn] = (double **)malloc(sizeof(double *) * (INBS + 4*NHP))) == NULL){
      fprintf(stderr, "\nError : malloc\n\n");exit(1);
    }
    ydm[zn] += 2*NHP;
    yd[zn] += 2*NHP;
    ydp[zn] += 2*NHP;
    for(n = -2*NHP ; n < INBS + EXBS*2 + 2*NHP ; n++){
      if((ydm[zn][n] = (double *)malloc(sizeof(double ) * (BSMM + 4*NHP))) == NULL ||
	 (yd[zn][n] = (double *)malloc(sizeof(double ) * (INBS + EXBS*2 + 4*NHP))) == NULL ||
	 (ydp[zn][n] = (double *)malloc(sizeof(double ) * (BSMM + 4*NHP))) == NULL){
	fprintf(stderr, "\nError : malloc\n\n");exit(1);
      }
      ydm[zn][n] += 2*NHP;
      yd[zn][n] += 2*NHP;
      ydp[zn][n] += 2*NHP;
    }
  }
  

  // dy01, dy12, dy20[BSNN][BSMM]
  if ((dy01  = (double **)malloc(sizeof(double *) * (BSNN + 2*NHP)) ) == NULL ||
      (dy12  = (double **)malloc(sizeof(double *) * (BSNN + 2*NHP)) ) == NULL ||
      (dy20  = (double **)malloc(sizeof(double *) * (BSNN + 2*NHP)) ) == NULL) {
        fprintf(stderr, "\nError : malloc\n\n"); exit(1);
  }
  dy01 += NHP;
  dy12 += NHP;
  dy20 += NHP;
  for(n = -NHP ; n < BSNN  + NHP; n++) {
    if ((dy01[n]  = (double *)malloc(sizeof(double) * (BSMM + 2*NHP)) ) == NULL ||
        (dy12[n]  = (double *)malloc(sizeof(double) * (BSMM + 2*NHP)) ) == NULL ||
        (dy20[n]  = (double *)malloc(sizeof(double) * (BSMM + 2*NHP)) ) == NULL) {
          fprintf(stderr, "\nError : malloc\n\n"); exit(1);
    }
    dy01[n] += NHP;
    dy12[n] += NHP;
    dy20[n] += NHP;
  }

  // p[zn][2*NHP+1][2*NHP+1]
  if((p = (double ***)malloc(sizeof(double **) * NZ)) == NULL){
    fprintf(stderr, "\nError : malloc\n\n");exit(1);
  }
  for(zn = 0 ; zn < NZ ; zn++){
    if((p[zn] = (double **)malloc(sizeof(double *) * (2*NHP+1))) == NULL){
      fprintf(stderr, "\nError : malloc\n\n");exit(1);
    }
    p[zn] += NHP;
    for(n = -NHP ; n < 1+NHP ; n++){
      if((p[zn][n] = (double *)malloc(sizeof(double ) * (2*NHP+1))) == NULL){
        fprintf(stderr, "\nError : malloc\n\n");exit(1);    
      }
      p[zn][n] += NHP;
    }
  }

  // d
  if ( (d = (double **)malloc(sizeof(double *) * (NNX*K + 2*NHLPF + 4*NHP*K)) ) == NULL) {
    fprintf(stderr, "\nError : malloc\n\n"); exit(1);
  }
  d += NHLPF + 2*NHP*K;
  for(n = -(NHLPF + 2*NHP*K) ; n < NNX*K + NHLPF + 2*NHP*K ; n++) {
    if ( (d[n] = (double *)malloc(sizeof(double) * (NNY*K + 2*NHLPF + 4*NHP*K)) ) == NULL) {
      fprintf(stderr, "\nError : malloc\n\n"); exit(1);
    }
    d[n] += NHLPF + 2*NHP*K;
  }

  // d_model, d_est, k_est[NB][MB], dm,dp_est
  if ((d_model  = (double **)malloc(sizeof(double *) * NNB) ) == NULL ||
      (dm_model  = (double **)malloc(sizeof(double *) * NNB) ) == NULL ||
      (dp_model  = (double **)malloc(sizeof(double *) * NNB) ) == NULL ||
      (d_est  = (double **)malloc(sizeof(double *) * NNB) ) == NULL ||
      (dm_est  = (double **)malloc(sizeof(double *) * NNB) ) == NULL ||
      (dp_est  = (double **)malloc(sizeof(double *) * NNB) ) == NULL ||
      (k_est  = (double **)malloc(sizeof(double *) * NNB) ) == NULL) {
        fprintf(stderr, "\nError : malloc\n\n"); exit(1);
  }
  for(n = 0 ; n < NNB ; n++) {
    if ((d_model[n]  = (double *)malloc(sizeof(double) * MMB) ) == NULL ||
	(dm_model[n]  = (double *)malloc(sizeof(double) * MMB) ) == NULL ||
	(dp_model[n]  = (double *)malloc(sizeof(double) * MMB) ) == NULL ||
        (d_est[n]  = (double *)malloc(sizeof(double) * MMB) ) == NULL ||
	(dm_est[n]  = (double *)malloc(sizeof(double ) * MMB) ) == NULL ||
	(dp_est[n]  = (double *)malloc(sizeof(double ) * MMB) ) == NULL ||
        (k_est[n]  = (double *)malloc(sizeof(double) * MMB) ) == NULL) {
          fprintf(stderr, "\nError : malloc\n\n"); exit(1);
    }
  }

  // d_sff, k_sff[NB][MB]
  if ((d_sff  = (double **)malloc(sizeof(double *) * NB) ) == NULL ||
      (k_sff  = (double **)malloc(sizeof(double *) * NB) ) == NULL) {
       fprintf(stderr, "\nError : malloc\n\n"); exit(1); 
  }
  for(n = 0 ; n < NB ; n++) {
    if ((d_sff[n]  = (double *)malloc(sizeof(double) * MB) ) == NULL ||
        (k_sff[n]  = (double *)malloc(sizeof(double) * MB) ) == NULL) {
          fprintf(stderr, "\nError : malloc\n\n"); exit(1);
    }
  }

  // E
  if ( (E = (double ***)malloc(sizeof(double **) * NNB) ) == NULL) {
    fprintf(stderr, "\nError : malloc\n\n"); exit(1);
  }
  for(n = 0 ; n < NNB ; n++) {
    if ( (E[n] = (double **)malloc(sizeof(double *) * MMB) ) == NULL) {
      fprintf(stderr, "\nError : malloc\n\n"); exit(1);
    }
    for(m = 0 ; m < MMB ; m++){    
      if ( (E[n][m] = (double *)malloc(sizeof(double) * (int)(DZ / DD)) ) == NULL) {
	      fprintf(stderr, "\nError : malloc\n\n"); exit(1);
      }
    }
  }

  //---------------------------------------------------------------------------
  // d[]に対象物体のdepth形状を設定する関数を作りたいです
#if HEIMEN == 0

#endif

#if HEIMEN == 1

#endif

  // d[]の保存   
  sprintf(flnm, "./d.dat");
  fp = fopen(flnm, "w");
  //for(n = -(2*NHP*K + NHLPF) ; n < (NX*K + 2*NHP*K + NHLPF) ; n++){
  //for(m = -(2*NHP*K + NHLPF) ; m < (NY*K + 2*NHP*K + NHLPF) ; m++){}
  for(n = 0 ; n < NNX ; n++){
    for(m = 0 ; m < NNY ; m++){
      fprintf(fp, "%d %d %f\n", n, m, d[n*K][m*K]);
    }
    fprintf(fp, "\n");
  }
  fclose(fp);

  //---------------------------------------------------------------------------  
  // freedする data -> goへ
#if HEIMEN == 0
  for(zn = 0 ; zn < NZ ; zn++){
    //sprintf(flnm, "./zn3_0515/d0515_bi_g%d", zn);
    sprintf(flnm, "./../obs_img/k16/d0515/d0515_ra_0%d.bin", zn);
    fp = fopen(flnm, "r");
    //nr = fread(data, sizeof(float), (NHP+NX+NHP) * (NHP+NY+NHP), fp);
    nr = fread(data, sizeof(double), (NHP+NX+NHP) * (NHP+NY+NHP), fp);
    printf("nr = %f\n", nr);
    nw = (NHP+NX+NHP) * (NHP+NY+NHP);
    if(nr != nw){
      fprintf(stderr, "\nError : fread\n\n");exit(1);
    }

    for(m = -NHP ; m < NY + NHP ; m++){
      for(n = -NHP ; n < NX + NHP ; n++){
        go[zn][n][m] = (double)data[(m+NHP) * (NHP+NX+NHP) + (n+NHP)];
      }
    }
  }
#endif

  // 平面
#if HEIMEN == 1
  for(zn = 0 ; zn < NZ ; zn++){
    sprintf(flnm, "./d09/d09_bi_g%d", zn);
    fp = fopen(flnm, "r");
    nr = fread(data, sizeof(float), (NHP+NX+NHP) * (NHP+NY+NHP), fp);
    nw = (NHP+NX+NHP) * (NHP+NY+NHP);
    if(nr != nw){
      fprintf(stderr, "\nError : fread\n\n");exit(1);
    }
    for(m = -NHP ; m < NY + NHP ; m++){
      for(n = -NHP ; n < NX + NHP ; n++){
	go[zn][n][m] = (double)data[(m+NHP) * (NHP+NX+NHP) + (n+NHP)];
      }
    }
  }
#endif

  //---------------------------------------------------------------------------
  // go[zn][]に観測雑音を付与、量子化してgo[zn][]に
  for(zn = 0 ; zn < NZ ; zn++){
    for(n = -NHP ; n < NX + NHP ; n++){
      for(m = -NHP ; m < NY + NHP ; m++){
        
#if ONOISE_FLG == 1
        go[zn][n][m] = go[zn][n][m] + ON_SGM * gauss();
#endif

#if QUANT_FLG == 1
        go[zn][n][m] = floor(go[zn][n][m] + 0.5);
#endif

#if LIMIT_FLG == 1
        if (go[zn][n][m] > 255.0) {
          go[zn][n][m] = 255.0;
        } else if (go[zn][n][m] < 0.0) {
          go[zn][n][m] = 0.0;
        } else {
          go[zn][n][m] = floor(go[zn][n][m] + 0.5);
        }
#endif
      }
    }

    //---------------------------------------------------------------------------
    // go[zn][]の保存
    sprintf(flnm, "./model/go%d.dat", zn);
    cIMAGE(NX + 2*NHP, NY + 2*NHP, &img_go, MONO);
    fp = fopen(flnm, "w");
    for(n = -NHP ; n < NX + NHP ; n++){
      for(m = -NHP ; m < NY + NHP ; m++){
        fprintf(fp, "%4d %4d  %f\n", n, m, go[zn][n][m]);
        D(img_go, n+NHP, m+NHP) = (unsigned char)(go[zn][n][m] + 0.5);
      }
    }
    fclose(fp);
    sprintf(flnm, "./model/go%d.bmp", zn);
    wIMAGE(img_go, flnm);
  }//zn

  printf("---orikaesi start---\n");
  // 折り返し処理 折り返した部分は誤差計算には使わない
  for(zn = 0 ; zn < NZ ; zn++){
    // left side
    for(m = 0 ; m < NNY ; m++){
      for(n = 1 ; n <= NHP ; n++){
	go_t[zn][-n][m] = go_t[zn][n][m];
      }
    }
    
    // right side
    for(m = 0 ; m < NNY ; m++){
      for(n = 1 ; n <= NHP ; n++){
	go_t[zn][NNX - 1 + n][m] = go_t[zn][NNX - 1 - n][m];
      }
    }
    
    // up side
    for(n = -NHP ; n < NNX + NHP ; n++){
      for(m = 1 ; m <= NHP ; m++){
	go_t[zn][n][-m] = go_t[zn][n][m];
      }
    }
    
    // down side
    for(n = -NHP ; n < NNX + NHP ; n++){
      for(m = 1 ; m <= NHP ; m++){
	      go_t[zn][n][NNY - 1 + m] = go_t[zn][n][NNY - 1 - m];
      }
    }


  //ピッチ分だけ動かしながら、ブロックの中心の座標とブロックサイズ等をGause-Newton関数へ
  
}
*/
