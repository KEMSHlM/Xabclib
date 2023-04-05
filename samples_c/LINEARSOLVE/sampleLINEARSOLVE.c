#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define FILENAME "../MatrixData/ecl32.dat"

#define OpenATI_INIT openati_init_
#define OpenATI_DURMV_11 openati_durmv_11_
#define OpenATI_LINEARSOLVE openati_linearsolve_
void OpenATI_INIT(int *iatparam, double *ratparam, int *info);
void OpenATI_DURMV_11(int *n, int *nnz, int *irp, int *icol, double *val,
                      double *x, double *y);
void OpenATI_LINEARSOLVE(int *n, int *nnz, int *irp, int *icol, double *val,
                         double *b, double *x, int *iatparam, double *ratparam,
                         int *info);

int main(int argc, char *argv[])
{
  FILE *fp;
  int i, n, nnz, info;
  int iatparam[50];
  int *irp, *icol;
  double beta0 = 0.0, beta = 0.0;
  double *x, *b, *wk, *val, ratparam[50];

  if( (fp = fopen(FILENAME, "r")) == NULL ){
    fprintf(stdout, "File not found.\n");
    exit(1);
  }
  fscanf(fp, "%d %d", &n, &nnz);
  printf("n=%d nnz=%d\n", n, nnz);
  fflush(stdout);

  x = (double *)malloc(n * sizeof(double));
  b = (double *)malloc(n * sizeof(double));
  wk = (double *)malloc(n * sizeof(double));
  irp = (int *)malloc((n + 1) * sizeof(int));
  icol = (int *)malloc(nnz * sizeof(int));
  val = (double *)malloc(nnz * sizeof(double));
  if( x==NULL || b==NULL || wk==NULL ||
      irp==NULL || icol==NULL || val==NULL ) exit(1);

  for(i = 0; i < n + 1; i++) fscanf(fp,"%d", &irp[i]);
  for(i = 0; i < nnz; i++) fscanf(fp,"%d", &icol[i]);
  for(i = 0; i < nnz; i++) fscanf(fp,"%lf", &val[i]);
  fclose(fp);

  for(i = 0; i < n; i++) x[i] = 1.0;
  OpenATI_DURMV_11(&n, &nnz, irp, icol, val, x, b);
  for(i = 0; i < n; i++) x[i] = 0.0;

  OpenATI_INIT(iatparam, ratparam, &info);

  iatparam[50-1] = 1; // debug info (0:Off, 1:On)

  OpenATI_LINEARSOLVE(&n,&nnz,irp, icol, val, b, x, iatparam, ratparam, &info);
  if( info ) printf("**ERROR** info=%d\n", info);

  for(i = 0; i < n; i++) beta0 += b[i] *  b[i];
  OpenATI_DURMV_11(&n, &nnz, irp, icol, val, x, wk);
  for(i = 0; i < n; i++) beta += (wk[i] - b[i]) * (wk[i] - b[i]);
  printf("||b-Ax_i||_2 / ||b||_2 = %e\n", sqrt(beta) / sqrt(beta0));

  free(val);
  free(icol);
  free(irp);
  free(wk);
  free(b);
  free(x);
  return 0;
}

