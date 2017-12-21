#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <imageio.h>
#include <math.h>
#include <sys/stat.h>
#include <time.h>
#include <nd_malloc.h>

/**** Default値の設定 ****/
#define BLOCK_SIZE  13
#define PITCH       13
#define RR           1

// フラグ
#define QUANT_FLG    0
#define ONOISE_FLG   0
#define LIMIT_FLG    0
#define SRAND_FLG    0
#define HEIMEN       0

#define NX         320
#define NY         240
#define NZ           3

#define D0         0.5
#define D1         1.5

#define K            8

#define NHP         11
#define NHLPF       25

int main(int argc, char *argv[]){
  FILE *fp;
  IMAGE img_go;
  char flnm[256];
  int nr, nw, zn, n, m, bn, bm, b0[2], b1[2];
  double ***go, ***go_t, **model_depth, **model_d;
  double **d, **d_est, **k_est, tmp;//**d_model
  double data[(NHP+NX+NHP) * (NHP+NY+NHP)];

  time_t t1, t2;
  t1 = time(NULL);
  
  int pitch = PITCH, block_size = BLOCK_SIZE, r = RR;
  char c;
  void usage(char *);
  void gaussNewtonMethod(double ***fn, int nz, int BS, int center_x, int center_y, double *dk);
  
#if SRAND_FLG == 1
  srandom(time(NULL));
#endif
  
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

  printf("\nBS = %d\n", block_size);
  printf("pitch = %d\n", pitch);
  printf("R = %d\n", r);

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
  go_t = malloc_double_3d(NZ, 0, NNX + 2*NHP, NHP, NNY + 2*NHP, NHP);

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
    sprintf(flnm, "/home/v4/tyoshida/imagep/focus/kmfg/md3/gauss20/pd/k05/py/md_ra_0%d.bin", zn);
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
  sprintf(flnm, "/home/v4/tyoshida/imagep/focus/kmfg/md3/gauss20/pd/k05/py/md_model.bin");
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
  }

  //---------------------------------------------------------------------------
  //ピッチ分だけ動かしながら、ブロックの中心の座標とブロックサイズ等をGause-Newton関数へ
  double dk[2];
  for(bm = 0 ; bm <= num_bsy ; bm++){
    for(bn = 0 ; bn <= num_bsx ; bn++){
      gaussNewtonMethod(go, NZ, block_size, hbsn + pitch*bn, hbsm + pitch*bm, dk);
      d_est[bn][bm] = dk[0];
      k_est[bn][bm] = dk[1];
    }
  }
  //gaussNewtonMethod(go, block_size, 27, 27, dk);
  
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
  double RMSE;
  tmp = 0.0;
  RMSE = 0.0;
  for(bm = r ; bm <= num_bsy - r ; bm++){
    for(bn = r ; bn <= num_bsx - r ; bn++){
      tmp = model_d[bn][bm] - d_est[bn][bm];
      RMSE += tmp * tmp;
      sum++;
    }
  }
  printf("all block num = %d\n", sum);
  RMSE = sqrt(RMSE / (double)(sum));
  printf("RMSE = %f\n", RMSE);
  
  t2 = time(NULL);
  printf("\n");
  printf("[time] : %d[minutes]\n\n", (int)((t2 - t1)/60));
}
  
/*--------------------------------------------------------------------------*/
void usage(char *com){
  fprintf(stderr, "\n");
  fprintf(stderr, "  Usage : %s [-bpr]\n\n", com);
  fprintf(stderr, "          -b <value> : BlockSize Odd   (int)\n");
  fprintf(stderr, "          -p <value> : Block Pitch     (int)\n");
  fprintf(stderr, "          -r <value> : Exclusion Block (int)\n");
  fprintf(stderr, "\n");

  exit(1);
}
