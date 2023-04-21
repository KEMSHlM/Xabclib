#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define FILENAME "../MatrixData/ecl32.dat"

#define OpenATI_INIT openati_init_
#define OpenATI_DURMV_11 openati_durmv_11_
#define Xabclib_BICGSTAB xabclib_bicgstab_
extern void OpenATI_INIT(int *iatparam, double *ratparam, int *info);
extern void OpenATI_DURMV_11(int *n, int *nnz, int *irp, int *icol, double *val,
                      double *x, double *y);
extern void Xabclib_BICGSTAB(int *n, int *nnz, int *irp, int *icol, double *val,
                   double *b, double *x, double *precond, int *npre, 
                   int *iatparam, double *ratparam, double *wk, int *lwk,
                   int *info);


int main(int argc, char *argv[])
{
  FILE *fp;
  int i, n, nnz, npre, lwk, info;
  int iatparam[50];
  int *irp, *icol;
  double beta0 = 0.0, beta = 0.0;
  double *x, *b, *wk, *wk2, *val, *precond, ratparam[50];

  // Read input matrix
  if( (fp = fopen(FILENAME, "r")) == NULL ){
    fprintf(stdout, "File not found.\n");
    exit(1);
  }
  fscanf(fp, "%d %d", &n, &nnz);
  printf("n=%d nnz=%d\n", n, nnz);
  fflush(stdout);

  x = (double *)malloc(n * sizeof(double));
  b = (double *)malloc(n * sizeof(double));
  wk2 = (double *)malloc(n * sizeof(double));
  // row_ptr
  irp = (int *)malloc((n + 1) * sizeof(int));
  // col_ind
  icol = (int *)malloc(nnz * sizeof(int));
  val = (double *)malloc(nnz * sizeof(double));
  if( x==NULL || b==NULL || wk2==NULL ||
      irp==NULL || icol==NULL || val==NULL ) exit(1);

  for(i = 0; i < n + 1; i++) fscanf(fp,"%d", &irp[i]);
  for(i = 0; i < nnz; i++) fscanf(fp,"%d", &icol[i]);
  for(i = 0; i < nnz; i++) fscanf(fp,"%lf", &val[i]);
  fclose(fp);

  // Set solution vecter "x"
  for(i = 0; i < n; i++) x[i] = 1.0;
  // Make RHS vecter "b" from A*x
  // xを1で初期化して，A*xをbに代入する．
  OpenATI_DURMV_11(&n, &nnz, irp, icol, val, x, b);
  // Set initial approximate solution vecter "x_0"
  for(i = 0; i < n; i++) x[i] = 0.0;

  // Initialize OpenATLib & Xabclib parameters
  // うーん，このOpenATI_DURMV_11で判定して，OpenATI_INITで値を代入してるのかな？
  // 規定値で初期化する．
  OpenATI_INIT(iatparam, ratparam, &info);

  // iatparamはソルバーの設定値
  // おそらく，50-1表記はfortranでいう50番目のパラメータ設定
  iatparam[50-1] = 1; // debug info (0:Off, 1:On)

  // Set flag of Auto-tuning for selecting SpMV implementation (0:off 3:on)
  iatparam[9-1]=0;
  iatparam[10-1]=12;

  // Set preconditioner type
  // 1:None, 2:Jacobi, 3:SSOR, 4:ILU0D, 5:ILU0, 6:ILUT
  // ILU0D: 0レベルのフィルイン.
  iatparam[25-1]=5;
  // 前処理行列にJacobi反復法に用いた場合，CG１回の反復に対して何回のJacobi反復を行うか
  iatparam[26-1]=5;
  // iatparam[25-1]が，4,5の場合に値を捨てるしきい値．
  ratparam[25-1]=1.0E-8;

  // Allocate work space for preconditioner
  // 前処理行列のサイズを決める. 公式ドキュメントに載ってる規定値を使う．
  if( iatparam[25-1]<=4){
    npre = n;
  } else if( iatparam[25-1]==5){
    npre = 3*nnz/2 + 2*n + 50;
  } else if( iatparam[25-1]==6){
    npre = 3*(2*iatparam[26-1]+1)*n/2 + 3*n + 50 + 1;
  }
  precond = (double *)malloc(npre * sizeof(double));

  // Allocate work space for bicgstab
  // 公式ドキュメントに載ってる規定値を使う．
  lwk = 9*n+(n-1)/2+1;
  wk = (double *)malloc(lwk * sizeof(double));
  if( wk==NULL || precond==NULL ) exit(1);

  // Call Xabclib_BICGSTAB
  // 先ほどxを１で求めたAx=bで求めたbよりxを求める．xの解が１になるまで反復する．
  Xabclib_BICGSTAB(&n,&nnz,irp,icol,val,b,x,precond,&npre,iatparam,ratparam,wk,&lwk,&info);
  if( info ) printf("**ERROR** info=%d\n", info);

  // Check true relative residual
  for(i = 0; i < n; i++) beta0 += b[i] *  b[i];
  OpenATI_DURMV_11(&n, &nnz, irp, icol, val, x, wk2);
  for(i = 0; i < n; i++) beta += (wk2[i] - b[i]) * (wk2[i] - b[i]);
  printf("||b-Ax_i||_2 / ||b||_2 = %e\n", sqrt(beta) / sqrt(beta0));

  free(val);
  free(icol);
  free(irp);
  free(wk);
  free(wk2);
  free(b);
  free(x);
  free(precond);
  return 0;
}

