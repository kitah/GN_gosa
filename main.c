#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <imageio.h>
#include <math.h>
#include <sys/stat.h>
#include <time.h>
#include <nd_malloc.h>
#include "pre_filter.h"

/**** Default値の設定 ****/
#define BLOCK_SIZE  13
#define PITCH       13
#define RR           1
#define SGMC      0.85

// フラグ
#define QUANT_FLG    0
#define ONOISE_FLG   0
#define LIMIT_FLG    0
#define SRAND_FLG    0
#define HEIMEN       0

#define ON_SGM      0.5

#define NX         320
#define NY         240
#define NZ           3

#define D0         0.5
#define D1         1.5

#define K            8

#define NHP         10
#define NHLPF       25
#define PRE_HS       3
#define PRE_N        0

int main(int argc, char *argv[]){
  FILE *fp;
  IMAGE img_go;
  char flnm[256];
  int nr, nw, zn, n, m, bn, bm, b0[2], b1[2];
  double ***go, ***go_t, ***go_pre, **model_depth, **model_d;
  double **d, **d_est, **k_est, tmp;// **d_model;
  double data[(NHP+NX+NHP) * (NHP+NY+NHP)];

  time_t t1, t2;
  t1 = time(NULL);
  
  int pitch = PITCH, block_size = BLOCK_SIZE, r = RR, quant_flg = QUANT_FLG, onoise_flg = ONOISE_FLG;
  char c;
  double sgmc = SGMC;
  void usage(char *);
  void gaussNewtonMethod(double ***fn, int nz, int BS, double sgmc, int center_x, int center_y, double *dk);
  double gauss(void);
  
#if SRAND_FLG == 1
  srandom(time(NULL));
#endif

  // gngosa -b 15 -p 10 -r 1
  //---------------------------------------------------------------------------
  // 引数の処理
  if (argc < 2) usage(argv[0]);

  while((c = getopt(argc, argv, "b:p:r:c:qn")) != -1){
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
    case 'c':
      sgmc = atof(optarg);
      break;
    case 'q':
      quant_flg = 1;
      break;
    case 'n':
      onoise_flg = 1;
      break;
    default:
      usage(argv[0]);
    }
  }
  //  argc -= optind - 1;
  //  argv += (optind - 1);
  //  if(argc > 2) usage(argv[0]);

  printf("\nBS = %d\n", block_size);
  printf("pitch = %d\n", pitch);
  printf("R = %d\n", r);
  printf("sgmc = %.2f\n", sgmc);
  if(quant_flg == 1) printf("quant: ON  ");
  else printf("quant: OFF  ");
  if(onoise_flg == 1) printf("noise: ON\n");
  else printf("noise: OFF\n");

  int num_bsx = 0, num_bsy = 0;
  int  hbsn = (block_size-1) / 2, hbsm = (block_size-1) / 2;
  n = m = 0;
  while(block_size + pitch*(n+1) < NX){
    n++;
    num_bsx++;
  }
  while(block_size + pitch*(m+1) < NY){
    m++;
    num_bsy++;
  }
  printf("block_num_x, block_num_y: %d %d\n", num_bsx+1, num_bsy+1);
  
  int NNX, NNY;
  NNX = block_size + num_bsx*pitch;
  NNY = block_size + num_bsy*pitch;
  printf("NNX : %d\n", NNX);
  printf("NNY : %d\n", NNY);
  
  //---------------------------------------------------------------------------
  // 領域確保 nd_malloc
  // go[zn][NX][NY]
  go = malloc_double_3d(NZ, 0, NX + 2*NHP, NHP, NY + 2*NHP, NHP);
  
  // go_t[zn][15*20][15*13]
  go_t = malloc_double_3d(NZ, 0, NNX + 2*NHP + 2*PRE_HS, NHP + PRE_HS, NNY + 2*NHP + 2*PRE_HS, NHP + PRE_HS);

  // go_pre[zn][][]
  go_pre = malloc_double_3d(NZ, 0, NNX + 2*NHP, NHP, NNY + 2*NHP, NHP);

  // model.binをfreedして格納するための配列
  model_depth = malloc_double_2d(NX, 0, NY, 0);
  model_d = malloc_double_2d(num_bsx+1, 0, num_bsy+1, 0);

  // d
  d = malloc_double_2d(NX*K + 2*NHLPF + 4*NHP*K, NHLPF + 2*NHP*K, NY*K + 2*NHLPF + 4*NHP*K, NHLPF + 2*NHP*K);

  // d_model, d_est, k_est [num_bsx][num_bsy]
  //d_model = malloc_double_2d(num_bsx+1, 0, num_bsy+1, 0);
  d_est = malloc_double_2d(num_bsx+1, 0, num_bsy+1, 0);
  k_est = malloc_double_2d(num_bsx+1, 0, num_bsy+1, 0);

  //---------------------------------------------------------------------------
  // d[]に対象物体のdepth形状を設定する
#if HEIMEN == 0
  for(n = -(2*NHP*K + NHLPF) ; n < (NX*K + 2*NHP*K + NHLPF) ; n++){
    for(m = -(2*NHP*K + NHLPF) ; m < (NY*K + 2*NHP*K + NHLPF) ; m++){    
      d[n][m] = D0 + (D1-D0) / (double)(NX*K - 1) * n;
    }
  }
#endif

#if HEIMEN == 1
  for(n = -(2*NHP*K + NHLPF) ; n < (NNX*K + 2*NHP*K + NHLPF) ; n++){
    for(m = -(2*NHP*K + NHLPF) ; m < (NNY*K + 2*NHP*K + NHLPF) ; m++){    
      d[n][m] = 0.9;
    }
  }
#endif
    
  // d[]の保存   
  sprintf(flnm, "./d.dat");
  fp = fopen(flnm, "w");
  for(n = 0 ; n < NX ; n++){
    for(m = 0 ; m < NY ; m++){
      fprintf(fp, "%d %d %f\n", n, m, d[n*K][m*K]);
    }
    fprintf(fp, "\n");
  }
  fclose(fp);
  
  //---------------------------------------------------------------------------  
  // freedする data -> goへ
#if HEIMEN == 0
  for(zn = 0 ; zn < NZ ; zn++){
    //sprintf(flnm, "./../obs_img/k05/d0515/d0515_ra_0%d.bin", zn);
    sprintf(flnm, "/home/v2unix/com_img/md_img/kmfg/md3/gauss20/pd/k03/eg/md_ra_0%d.bin", zn);
    fp = fopen(flnm, "r");
    nr = fread(data, sizeof(double), NX*NY, fp);
    nw = NX*NY;
    if(nr != nw){
      fprintf(stderr, "\nError : fread\n\n");exit(1);
    }

    for(m = 0 ; m < NY ; m++){
      for(n = 0 ; n < NX ; n++){
        go[zn][n][m] = (double)data[m*NX + n];
      }
    }
  }

  // model.binをfreedする
  sprintf(flnm, "/home/v2unix/com_img/md_img/kmfg/md3/gauss20/pd/k03/eg/md_model.bin");
  fp = fopen(flnm, "r");
  nr = fread(data, sizeof(double), NX*NY, fp);
  nw = NX*NY;
  if(nr != nw){
    fprintf(stderr, "\nError : fread\n\n");exit(1);
  }
  
  for(m = 0 ; m < NY ; m++){
    for(n = 0 ; n < NX ; n++){
      model_depth[n][m] = (double)data[m*NX + n];
    }
  }

  // modelでブロックごとに平均をとり、中心の点を求める
  for(bn = 0 ; bn <= num_bsx ; bn++){
    b0[0] = bn * pitch;    // 0, 10, 20, ... exp. pitch = 10 
    b1[0] = b0[0] + block_size; // 15, 25, ...
    for(bm = 0 ; bm <= num_bsy ; bm++){
      b0[1] = bm * pitch;
      b1[1] = b0[1] + block_size;

      tmp = 0.0;
      for(n = b0[0] ; n < b1[0] ; n++){
        for(m = b0[1] ; m < b1[1] ; m++){
          tmp += model_depth[n][m];
        }
      }
      model_d[bn][bm] = tmp / (double)(block_size * block_size);
    }
  }
  sprintf(flnm, "./model_d.dat");
  fp = fopen(flnm, "w");
  for(bn = r ; bn <= num_bsx - r ; bn++){
    for(bm = r ; bm <= num_bsy - r ; bm++){
      fprintf(fp, "%4d %4d  %f\n", bn*pitch + hbsn, bm*pitch + hbsm, model_d[bn][bm]);
    }
    fprintf(fp, "\n");
  }
  fclose(fp);
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
	// noise
	if(onoise_flg == 1) go[zn][n][m] = go[zn][n][m] + ON_SGM * gauss();

	// quant
	if(quant_flg == 1) go[zn][n][m] = floor(go[zn][n][m] + 0.5);

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

  // go_t <- go
  for(zn = 0 ; zn < NZ ; zn++){
    for(n = 0 ; n < NNX ; n++){
      for(m = 0 ; m < NNY ; m++){
	go_t[zn][n][m] = go[zn][n][m];
      }
    }
  }

  //---------------------------------------------------------------------------
  printf("------------------------\n");
  // 折り返し処理 折り返した部分は誤差計算には使わない
  for(zn = 0 ; zn < NZ ; zn++){
    // left side right side
    for(m = 0 ; m < NNY ; m++){
      for(n = 1 ; n <= NHP + PRE_HS ; n++){
	go_t[zn][-n][m] = go_t[zn][n][m];
	go_t[zn][NNX - 1 + n][m] = go_t[zn][NNX - 1 - n][m];
      }
    }
       
    // up side down side
    for(n = -(NHP + PRE_HS) ; n < NNX + (NHP + PRE_HS) ; n++){
      for(m = 1 ; m <= NHP + PRE_HS ; m++){
	go_t[zn][n][-m] = go_t[zn][n][m];
	go_t[zn][n][NNY - 1 + m] = go_t[zn][n][NNY - 1 - m];
      }
    }
    sprintf(flnm, "./model/go_t_%d.dat", zn);
    fp = fopen(flnm, "w");
    for(n = -NHP ; n < NNX + NHP ; n++){
      for(m = -NHP ; m < NNY + NHP ; m++){
        fprintf(fp, "%4d %4d  %f\n", n, m, go_t[zn][n][m]);
      }
    }
    fclose(fp);
    sprintf(flnm, "./model/go_t_%d.bmp", zn);
  }
  
  //---------------------------------------------------------------------------
  // HPFをかけて低い周波数成分を除去
  int ii, jj;
  for(zn = 0 ; zn < NZ ; zn++){
    for(bm = -NHP ; bm < NNY + NHP ; bm++){
      for(bn = -NHP ; bn < NNX + NHP ; bn++){
	tmp = 0.0;
	for(jj = -PRE_HS ; jj <= PRE_HS ; jj++){
	  for(ii = -PRE_HS ; ii <= PRE_HS ; ii++){
	    tmp += go_t[zn][bn - ii][bm - jj] * pre[PRE_N][ii + PRE_HS][jj + PRE_HS]; 
	  }
	}
	go_pre[zn][bn][bm] = tmp;
      }
    }
  }
  
  //---------------------------------------------------------------------------
  //ピッチ分だけ動かしながら、ブロックの中心の座標とブロックサイズ等をGause-Newton関数へ
  double dk[2];
  
  for(bm = 0 ; bm <= num_bsy ; bm++){
    for(bn = 0 ; bn <= num_bsx ; bn++){
      gaussNewtonMethod(go_pre, NZ, block_size, sgmc, hbsn + pitch*bn, hbsm + pitch*bm, dk);
      d_est[bn][bm] = dk[0];
      k_est[bn][bm] = dk[1];
    }
  }
  
  sprintf(flnm, "./d_est.dat");
  fp = fopen(flnm, "w");
  for(bn = r ; bn <= num_bsx - r ; bn++){
    for(bm = r ; bm <= num_bsy - r ; bm++){
      fprintf(fp, "%4d %4d  %f\n", bn*pitch + hbsn, bm*pitch + hbsm, d_est[bn][bm]);
    }
    fprintf(fp, "\n");
  }
  fclose(fp);

  //---------------------------------------------------------------------------  
  // ブロック単位のモデル形状 こっちの方がRMSE高い
  /*
  for(bn = 0 ; bn <= num_bsx ; bn++){
    b0[0] = bn * pitch;    // 0, 10, 20, ... exp. pitch = 10 
    b1[0] = b0[0] + block_size; // 15, 25, ...
    for(bm = 0 ; bm <= num_bsy ; bm++){
      b0[1] = bm * pitch;
      b1[1] = b0[1] + block_size;

      tmp = 0.0;
      for(n = b0[0] ; n < b1[0] ; n++){
        for(m = b0[1] ; m < b1[1] ; m++){
          tmp += d[n * K][m * K];
        }
      }

      d_model[bn][bm] = tmp / (double)(block_size * block_size);
    }
  }
  sprintf(flnm, "./d_model.dat");
  fp = fopen(flnm, "w");
  for(bn = r ; bn <= num_bsx - r ; bn++){
    for(bm = r ; bm <= num_bsy - r ; bm++){
      fprintf(fp, "%4d %4d  %f\n", bn*pitch + hbsn, bm*pitch + hbsm, d_model[bn][bm]);
    }
    fprintf(fp, "\n");
  }
  fclose(fp);
  */
  //---------------------------------------------------------------------------
  // RMSE -rで指定した分だけ外側のブロックは誤差計測に使用しない
  int sum = 0;
  double RMSE = 0.0;
  tmp = 0.0;
  RMSE = 0.0;
  for(bm = r ; bm <= num_bsy - r ; bm++){
    for(bn = r ; bn <= num_bsx - r ; bn++){
      printf("(x,y) : %d %d\n", bn, bm);
      tmp = model_d[bn][bm] - d_est[bn][bm];
      printf("tmp = %f\n", tmp);
      printf("tmp^2 = %f\n", tmp*tmp);
      RMSE += (tmp) * (tmp);
      sum++;
    }
  }
  printf("all block num = %d\n", sum);
  //printf("rmse = %f\n", RMSE);
  RMSE = sqrt(RMSE / (double)(sum));
  printf("RMSE = %f\n", RMSE);
  
  t2 = time(NULL);
  printf("\n");
  printf("[time] : %d[s]\n\n", (int)((t2 - t1)));
}
  
/*--------------------------------------------------------------------------*/
void usage(char *com){
  fprintf(stderr, "\n");
  fprintf(stderr, "  Usage : %s [-b:p:r:c:qn]\n\n", com);
  fprintf(stderr, "          -b <value> : BlockSize Odd   (int)\n");
  fprintf(stderr, "          -p <value> : Block Pitch     (int)\n");
  fprintf(stderr, "          -r <value> : Exclusion Block (int)\n");
  fprintf(stderr, "          -c <value> : SGMC (double)\n");
  fprintf(stderr, "          -q : quant flag (bool)\n");
  fprintf(stderr, "          -n : noise flag (bool)\n");
  fprintf(stderr, "\n");

  exit(1);
}

/*--------------------------------------------------------------------------*/  
double gauss(void){
  double rnd(void);

  return rnd() + rnd() + rnd() + rnd() + rnd() + rnd() +
    rnd() + rnd() + rnd() + rnd() + rnd() + rnd() - 6.0;
}

/*--------------------------------------------------------------------------*/  
double rnd(void){
  return (double)random() / (double) RAND_MAX;
}
