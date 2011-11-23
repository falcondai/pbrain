/*     Version March 1994  */

/* Procrustes point registration - read in lists of homologous points, 
   form matrix of scalar products of all point pairs, do singular
   value decomposition to diagonalize, yielding transformation which
   matches points in least squares sense.

C. Pelizzari, Sept 91
   
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#ifdef THINK_C
#include <console.h>
#endif

#include "nrutil.h"

extern char *GetaString();

typedef float XYZPOINT[3];
typedef float **NRMATRIX;
typedef float *NRVECTOR;


static void transpose(a, nl, na)

NRMATRIX a;
int nl, na;

{
  int i,j;
  float temp;
  
  for (i = nl; i <= na; i++) 
  
     for (j = nl; j <= na; j++) {
        temp = a[i][j];
        a[i][j] = a[j][i];
        a[j][i] = temp;
     }

  
  return;
}

static void UnloadNRmatrix(a, nrl, nrh, ncl, nch, dest)

NRMATRIX a;
int nrl, nrh, ncl, nch;
float *dest;

{

  int i, j;
  float *fptr = dest;

  for (i = nrl; i <= nrh; i++)
    for (j = ncl; j <= nch; j++, fptr++) (*fptr) = a[i][j];


} /*  */

static NRMATRIX LoadNRmatrix(source, nr, nc)

float *source;
int nr, nc;

{

  int i, j;
  NRMATRIX temp;

  temp = matrix(1, nr, 1, nc);

  for (i = 1; i <= nr; i++)
    for (j = 1; j <= nc; j++, source++)  temp[i][j] = *source;

  return(temp);

} /*  */

static float NRtrace(a, nrl, nrh)

NRMATRIX a;
int nrl, nrh;

{

  int i;
  float sum = 0.0;

  for (i = nrl; i <= nrh; i++) sum += a[i][i];

  return(sum);

} /*  */

static NRMATRIX NRmatrixSubtract( a, b, nrl, nrh, ncl, nch)

/* subtracts b from a */


NRMATRIX a, b;
int nrl, nrh, ncl, nch;

{
   int i,j,k;
   NRMATRIX temp;

   temp  = matrix(nrl, nrh, ncl, nch);
   
  
   for (i = nrl; i <= nrh; i++) 
     for (j = ncl; j <= nch; j++) {
        temp[i][j] = a[i][j] - b[i][j];                                   
     }

   
   return (temp);
}



static NRMATRIX NRmatrixProduct( a, b, na, nb, nc)

/* multiplies a (na x nb) times b (nb x nc) returning c (na x nc) */


NRMATRIX a, b;
int na, nb, nc;

{
   int i,j,k;
   NRMATRIX temp;

   temp  = matrix(1, na, 1, nc);
   
  
   for (i = 1; i <= na; i++) 
     for (j = 1; j <= nc; j++) {
        temp[i][j] = 0;
        for (k = 1; k <= nb; temp[i][j] += a[i][k] * b[k][j], k++);
     }

   
   return (temp);
}

static NRMATRIX NRmatrixTranspose( a, nrl, nrh, ncl, nch)

/* allocates and fills transpose of general NR matrix */


NRMATRIX a;
int nrl, ncl, nrh, nch;

{
   int i,j,k;
   NRMATRIX temp;

   temp  = matrix(ncl, nch, nrl, nrh);
   
  
   for (i = nrl; i <= nrh; i++) 
     for (j = ncl; j <= nch; j++) temp[j][i] = a[i][j];
   
   return (temp);
}

static NRVECTOR NRrowColNorm( a, nrl, nrh, ncl, nch, rows)

/* calculates row or column norm of matrix */


NRMATRIX a;
int nrl, ncl, nrh, nch, rows;

{
   int i,j,k;
   NRVECTOR temp;

   if (rows) {
     temp  = vector(nrl, nrh);
   
  
     for (i = nrl; i <= nrh; i++){

      temp[i] = 0.0;
       for (j = ncl; j <= nch; j++) temp[i] += a[i][j] * a[i][j];
      temp[i] = (float) sqrt((double)temp[i]);
     }  


    } /*    if (row)  */

    else    {
     temp  = vector(ncl, nch);
   
  
     for (i = ncl; i <= nch; i++){

      temp[i] = 0.0;
       for (j = nrl; j <= nrh; j++) {
           //printf("temp[%d] += a[%d][%d]=%f * a[%d][%d]=%f\n", i, j,i,a[j][i], j,i,a[j][i]);
           temp[i] += a[j][i] * a[j][i];
       }
      temp[i] = (float) sqrt((double)temp[i]);
     }  


    } /*    if (row)  */
   return (temp);
}

static void NRrowShift(mat, vec, nrl, nrh, ncl, nch, subtract)

NRMATRIX mat;
NRVECTOR vec;
int nrl, nrh, ncl, nch, subtract;

{

  int i, j;
  float factor;

  factor   = (subtract ? -1.0 : 1.0);

  for (i = nrl; i <= nrh; i++)
    for (j = ncl; j <= nch; j++) mat[i][j] += factor * vec[j];

  return;

} /*  */

static void NRscale(mat, scale, nrl, nrh, ncl, nch)

NRMATRIX mat;
float scale;
int nrl, nrh, ncl, nch;

{

  int i, j;


  for (i = nrl; i <= nrh; i++)
    for (j = nrl; j <= nch; j++) mat[i][j] *= scale;

  return;

} /*  */

static void NRcolShift(mat, vec, nrl, nrh, ncl, nch, subtract)

NRMATRIX mat;
NRVECTOR vec;
int nrl, nrh, ncl, nch, subtract;

{

  int i, j;
  float factor;

  factor   = (subtract ? -1.0 : 1.0);

  for (i = ncl; i <= nch; i++)
    for (j = nrl; j <= nrh; j++) mat[j][i] += factor * vec[j];

  return;

} /*  */



static void PrintNRmatrix(title, amat, il, ih, jl, jh, byrows)

char *title;
NRMATRIX amat;
int il, ih, jl, jh, byrows;

{
  int i,j;
  
  printf("%s", title);

  if (byrows)   
  for (i = il; i <= ih; i++)  {
    for(j = jl; j <= jh; j++) printf("%.2f ", amat[i][j]);
    printf("\n");
    
  }
  
  else /*  by columns */
  
    for(j = jl; j <= jh; j++)  {
     for (i = il; i <= ih; i++) printf("%.2f ", amat[i][j]);
     printf("\n");
    
   }

}


static void GaussJ(theMatrix, rank, theInverse)
float theMatrix[];
float theInverse[];
int rank;
/*  just copies matrix into a NR matrix and calls the NR "gaussj" function */
{
  int i, j;
  NRMATRIX a, b;

  a = LoadNRmatrix(theMatrix, rank, rank);
  /* gaussj wants to solve equations, not just invert a matrix.  give it
     zero RHS vectors to work with */
  b = matrix(1, rank, 1, rank);
  for (i = 1; i <= rank; i++)
    for (j = 1; j <= rank; j++) b[i][j] = 0.0;
  gaussj(a, rank, b, rank);

  UnloadNRmatrix(a, 1, rank, 1, rank, theInverse);
  free_matrix(a, 1, rank, 1, rank);
  free_matrix(b, 1, rank, 1, rank);
  return;


}

static HMXprintMatrix(theStream, theTitle, theMatrix)
FILE *theStream;
char *theTitle;
float theMatrix[4][4];

{
  int i, j;
  fprintf(theStream, "%s\n", theTitle);
  for(i = 0; i < 4; i++)
  {
    for(j = 0; j < 4; j++) fprintf(theStream, "%.4f ", theMatrix[i][j]);
    fprintf(theStream, "\n");
  }
  return;
}

static HMXvectorRotate(theVector, theMatrix, theResult)
float theVector[], theResult[];
float theMatrix[4][4];

{
  int i, j;
  float tmp[4];
  int k,l;

  /*
  printf("theMatrix:\n");
  for (k =0; k<4; k++) {
      for (l =0; l<4; l++) {
          printf ("%f ", theMatrix[k][l]);
      }
      printf("\n");

  }

  printf("vec is ");for (i = 0; i <= 3; i++) printf("%.2f ", theVector[i]); printf("\n");
  */

  for (i = 0; i < 4; i++)
  {
    tmp[i] = 0;
    for (j = 0; j < 4; j++)
      tmp[i] += theVector[j] * theMatrix[j][i];

  }
  for(i = 0; i < 4; i++) theResult[i] = tmp[i];

}


Mult4x4(Matrix, Other)
NRMATRIX Matrix, Other;

{
  NRMATRIX tmp, product;
  int i,j;

  tmp = matrix(1,4,1,4);

  tmp = NRmatrixProduct(Matrix, Other, 4, 4, 4);

  for(i = 1; i <= 4; i++)
    for (j = 1; j <= 4; j++)
      Other[i][j] = tmp[i][j];

  free_matrix(tmp, 1, 4, 1, 4);

  return;


}

Translate4x4(Matrix, Translation, sense)
NRMATRIX Matrix;
NRVECTOR Translation;
int sense;

{
  NRMATRIX tmp, product;
  int i,j;

  tmp = matrix(1,4,1,4);
  for(i = 1; i <= 4; i++)
    for (j = 1; j <= 4; j++)
      tmp[i][j] = (i == j ? 1.0 : 0.0);

  for (i = 1; i <= 3; i++) tmp[4][i] = (sense == -1 ? -Translation[i]: Translation[i]);

  product = NRmatrixProduct(tmp, Matrix, 4, 4, 4);

  for(i = 1; i <= 4; i++)
    for (j = 1; j <= 4; j++)
      Matrix[i][j] = product[i][j];

  free_matrix(tmp, 1, 4, 1, 4);
  free_matrix(product, 1, 4, 1, 4);

  return;


}


Scale4x4(Matrix, Scale)
NRMATRIX Matrix;
NRVECTOR Scale;

{
  NRMATRIX tmp, product;
  int i,j;

  tmp = matrix(1,4,1,4);
  for(i = 1; i <= 4; i++)
    for (j = 1; j <= 4; j++)
      tmp[i][j] = (i == j ? 1.0 : 0.0);

  for (i = 1; i <= 3; i++) tmp[i][i] = Scale[i];

  product = NRmatrixProduct(tmp, Matrix, 4, 4, 4);

  for(i = 1; i <= 4; i++)
    for (j = 1; j <= 4; j++)
      Matrix[i][j] = product[i][j];

  free_matrix(tmp, 1, 4, 1, 4);
  free_matrix(product, 1, 4, 1, 4);

  return;


}

NRMATRIX BuildHomogeneousTransform(Center, Matrix, Translation, Scale, outmat)

XYZPOINT Center, Translation;
float Matrix[3][3], Scale, outmat[4][4];

{

  int i, j;
  NRMATRIX a, mat;
  NRVECTOR cent, trans, scale;

  scale = vector(1,3);
  cent = vector(1,3);
  trans = vector(1,3);

  for  (i = 1; i <= 3; i++) {

    scale[i] = Scale;
    cent[i] = Center[i-1];
    trans[i] = Translation[i-1];

  } /* for  (i = 1; i <= 3; i++)  */

  a = matrix(1,4,1,4);
  mat = matrix(1,4,1,4);

  for (i = 1; i <= 4; i++)
    for(j = 1; j <= 4; j++)
      mat[i][j] = a[i][j] = (i == j ? 1.0 : 0.0);
      
  
  for (i = 1; i < 4; i++)
    for(j = 1; j < 4; j++)
    {
      a[j][i] = Matrix[i-1][j-1];
    }
  printf ("cent= %f %f %f\n", cent[1],cent[2],cent[3]);
  PrintNRmatrix("mat= \n", mat, 1, 4, 1, 4, 1);
  PrintNRmatrix("a= \n", a, 1, 4, 1, 4, 1);
  Translate4x4(mat, cent, 1);

  PrintNRmatrix("after Translate4x4 , mat= \n", mat, 1, 4, 1, 4, 1);
  printf ("cent= %f %f %f\n", cent[1],cent[2],cent[3]);

  PrintNRmatrix("mat= \n", mat, 1, 4, 1, 4, 1);

  
  Mult4x4(a, mat);

  PrintNRmatrix("mat= \n", mat, 1, 4, 1, 4, 1);
  

  Scale4x4(mat, scale);

  PrintNRmatrix("mat after Scale4x4= \n", mat, 1, 4, 1, 4, 1);


  Translate4x4(mat, trans, 1);

  PrintNRmatrix("mat after Translate4x4= \n", mat, 1, 4, 1, 4, 1);


  Translate4x4(mat, cent, -1);
 
  PrintNRmatrix("mat after Translate4x4= \n", mat, 1, 4, 1, 4, 1);
  

  free_matrix(a, 1, 4, 1, 4);

  free_vector(scale, 1, 3);
  free_vector(cent, 1, 3);
  free_vector(trans, 1, 3);

  UnloadNRmatrix(mat, 1, 4, 1, 4, outmat);
  free_matrix(mat, 1, 4, 1, 4);

  return;

} /*  */





void ProcrustesMatch(npoints, PointList1, PointList2, Center, Matrix, 
    Translation, Scale, MeanError, MeanSquareError, RMSError, doscale, 
      ifprint)

int npoints;
XYZPOINT *PointList1, *PointList2;
XYZPOINT Center;
float Matrix[3][3], *Scale;
XYZPOINT Translation;
int doscale, ifprint;
float *MeanError, *MeanSquareError, *RMSError;

/*  Procrustes:  finds the rigid body transformation which matches a 3D
    set of points to a homologous set, minimizing the mean squared mismatch
    between the point pairs.  The transformation found fits PointList1
    to PointList2.  To apply it, the caller must do the following:

      1)  Add "Translation" to the PointList1 space coordinates
      2)  Rotate by "Matrix" about "Center", i.e. subtract Center, rotate
          by Matrix, add Center back again.

    This function optionally prints out an analysis of the point-by-point
    misfits after fitting.

    Algorithm is as given by Schonemann, Psychometrika vol 35, p. 245 (1970).
    (see eq. 2.11).

    Needs to link with svdcmp and nrutil from the Numerical Recipes in
    C package.

    C. Pelizzari, Sept 91.

*/



{

  int i,j,k,l,m;

  XYZPOINT Centroid1, Centroid2;
  
  NRMATRIX  a, b, u, v, r, p1, p2, t1, t2, r1, delta;
  NRMATRIX norm1, norm2;
  NRMATRIX  part1, part2, vt, ut;

  NRVECTOR w, distance;
  
  float mean, ms, rms, scale;

  if (ifprint) {


  printf("Procrustes: %d\n", npoints);

  for (i = 0; i < npoints; i++) {

    for (j = 0; j < 3; j++) printf("%6.2f ", PointList1[i][j]);
    for (j = 0; j < 3; j++) printf("%6.2f ", PointList2[i][j]);
    printf("\n");
  } /* for (i = 0; i < npoints; i++)  */
  } /*   if (ifprint)  */


  a = matrix(1, 3, 1, 3);
  w = vector(1, 3);
  v = matrix(1, 3, 1, 3);
  
  p1 = LoadNRmatrix(PointList1, npoints, 3);
  p2 = LoadNRmatrix(PointList2, npoints, 3);
 
    for (j = 0; j < 3; j++) {
      Centroid1[j] = 0.0;
      Centroid2[j] = 0.0;
    } 


  for (i = 0; i < npoints; i++) 
    for (j = 0; j < 3; j++) {
      Centroid1[j] += PointList1[i][j];
      Centroid2[j] += PointList2[i][j];
    } 

    for (j = 0; j < 3; j++) {
      Centroid1[j] /= npoints;
      Centroid2[j] /= npoints;
    } 

  printf("Centroid 1 is ");for (i = 0; i < 3; i++) printf("%.2f ", Centroid1[i]); printf("\n");
  printf("Centroid 2 is ");for (i = 0; i < 3; i++) printf("%.2f ", Centroid2[i]); printf("\n");

    

  for (i = 1; i <= npoints; i++) {

    for (j = 1; j <= 3; j++) {
      p1[i][j] -= Centroid1[j-1];
      p2[i][j] -= Centroid2[j-1];
    } /*     for (j = 1; j <= 3; j++)  */


  } /*   for (i = 1; i <= npoints; i++)  */

  PrintNRmatrix("p1= \n", p1, 1, npoints, 1, 3, 1);
  PrintNRmatrix("p2= \n", p2, 1, npoints, 1, 3, 1);


  for (i = 1; i <= 3; i++)
    for (j = 1; j <= 3; j++)  a[i][j] = 0.0;
  
  
     for (j = 1; j <= npoints; j++) {
     

       for (k = 1; k <= 3; k++)
         for (l = 1; l <= 3; l++) 
            a[k][l] += p1[j][k] * p2[j][l];
       
     }
  
  PrintNRmatrix("a= \n", a, 1, 3, 1, 3, 1);
 
  svdcmp(a, 3, 3, w, v);
  u = a;
  PrintNRmatrix("u= \n", u, 1, 3, 1, 3, 1);
  //PrintNRmatrix("w= \n", w, 1, 1, 1, 3, 1);
  printf("w is ");for (i = 1; i <= 3; i++) printf("%.2f ", w[i]); printf("\n");
  PrintNRmatrix("v= \n", v, 1, 3, 1, 3, 1);

  ut = NRmatrixTranspose(u, 1, 3, 1, 3);
  vt = NRmatrixTranspose(v, 1, 3, 1, 3);

  if (ifprint)
      PrintNRmatrix("This is v: \n", v, 1, 3, 1, 3, 1);
  if (ifprint)
      PrintNRmatrix("This is the ut: \n", ut, 1, 3, 1, 3, 1);

  r = NRmatrixProduct(  v, ut, 3, 3, 3);
  if (ifprint)
      PrintNRmatrix("This is the rotation matrix: \n", r, 1, 3, 1, 3, 1);

  t1 = NRmatrixTranspose(p1, 1, npoints, 1, 3);
  PrintNRmatrix("p1: \n", p1, 1, npoints, 1, 3, 1);
  PrintNRmatrix("t1: \n", t1, 1, 3, 1, npoints, 1);

  norm1 = NRmatrixProduct(p1, t1, npoints, 3, npoints);  /* A x A' */
  PrintNRmatrix("norm1: \n", norm1, 1, npoints, 1, npoints);
  

  t2 = NRmatrixTranspose(p2, 1, npoints, 1, 3);

  r1 = NRmatrixProduct(r, t1, 3, 3, npoints);

  PrintNRmatrix("r1: \n", r1, 1, 3, 1, npoints);

  norm2 = NRmatrixProduct(p2, r1,  npoints, 3, npoints); /* B x RA' */

  PrintNRmatrix("norm2: \n", norm2, 1, npoints, 1, npoints);

  *Scale = scale = NRtrace(norm2, 1, npoints) / NRtrace(norm1, 1, npoints);
    /* scale according to Schonemann & Carroll's recipe */
    /* eq. 2.12, p. 247. */

  if (ifprint) printf("scale is %.2f\n", scale);

  for (i = 1; i <= npoints; i++) {

    for (j = 1; j <= 3; j++) {
      printf("r1[%d][%d]=%f\n", j, i, r1[j][i]);
      if (doscale) r1[j][i] *= scale;
      printf("r1[%d][%d]=%f\n", j, i, r1[j][i]);
      r1[j][i] += Centroid2[j-1];
      t2[j][i] += Centroid2[j-1];
      printf("r1[%d][%d]=%f\n", j, i, r1[j][i]);
    } /*     for (j = 1; j <= 3; j++)  */

  } /*   for (i = 1; i <= npoints; i++)  */

  if (ifprint) {

  printf("Centroid 1 is ");for (i = 0; i < 3; i++) printf("%.2f ", Centroid1[i]); printf("\n");
  printf("Centroid 2 is ");for (i = 0; i < 3; i++) printf("%.2f ", Centroid2[i]); printf("\n");
  } /*   if (ifprint)  */

  mean = ms = 0.0;

  PrintNRmatrix("t2: \n", t2, 1, 3, 1, npoints);
  PrintNRmatrix("r1: \n", r1, 1, 3, 1, npoints);
  
  
  delta = NRmatrixSubtract(t2, r1, 1, 3, 1, npoints);
  PrintNRmatrix("delta: \n", delta, 1, 3, 1, npoints);

  distance = NRrowColNorm(delta, 1, 3, 1, npoints, 0);
  printf("distance is ");for (i = 1; i <= npoints; i++) printf("%.2f ", distance[i]); printf("\n");
  

  for (i = 1; i <= npoints; i++) {

  if (ifprint) 
    printf("(%6.2f, %6.2f, %6.2f)     %.2f\n", 
       delta[1][i], delta[2][i], delta[3][i], distance[i]);

    mean += distance[i];
    ms += distance[i] * distance[i];

  } /* for (i = 0; i <= npoints; i++)  */

  mean /= npoints;
  ms /= npoints;
  rms = (float) sqrt((double)ms);

  *MeanError = mean;    /* return residual values to caller */
  *MeanSquareError = ms;
  *RMSError = rms;

  if (ifprint) 
  printf("mean distance = %.2f; mean square distance = %.2f; rms distance = %.2f\n", mean, ms, rms);

  UnloadNRmatrix(r, 1, 3, 1, 3, Matrix);
  for (i = 0; i < 3; i++) {

    Center[i] = Centroid2[i];
    Translation[i] = Centroid2[i] - Centroid1[i];

  } /* for (i = 1; i <= 3; i++)  */
 
  free_matrix(p1, 1, npoints, 1, 3);
  free_matrix(p2, 1, npoints, 1, 3);

  free_matrix(t1, 1, 3, 1, npoints);
  free_matrix(t2, 1, 3, 1, npoints);

  free_matrix(r, 1, 3, 1, 3);
  free_matrix(r1, 1, 3, 1, npoints);

  free_matrix(norm1, 1, npoints, 1, npoints);
  free_matrix(norm2, 1, npoints, 1, npoints);

  free_matrix(a, 1, 3, 1, 3);
  free_matrix(ut, 1, 3, 1, 3);
  free_matrix(vt, 1, 3, 1, 3);
  free_matrix(v, 1, 3, 1, 3);

  free_matrix(delta, 1, 3, 1, npoints);
  free_vector(distance, 1, npoints);

  return;
}  

#ifdef TEST

#define DIST3D(a) sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2])

static int GetInt(fd, IntVal)

int *IntVal;
FILE *fd;
  
{ 
  
  *IntVal = 0;
  fscanf(fd, "%d", IntVal);

  return (*IntVal );
  
}

static int ReadPointData( filename,list1, list2)


 XYZPOINT **list1, **list2;

 char *filename; /* the file to be read from */

{
  int i, j, n, npts;
  FILE *fd;
  float *fptr;

  char inbuf[132], *ptr, *end;

  XYZPOINT *buf1;
   
  if (ptr = strchr(filename, '~')) *ptr = 0;
    
  if (  ( fd = fopen(filename, "r") ) != NULL)  {  /* found the file - read it  */
   
    if (  GetInt(fd, &npts) > 0 )  {  /*  if any points in this file  */

      buf1 = (XYZPOINT *) malloc (2 * npts * sizeof (XYZPOINT));  
      /* space for coordinates*/
 
      
      n = 0;

        fgets(inbuf, sizeof(inbuf) - 1, fd);
      fptr = (float *) buf1;

      while (n < 6 * npts) {  

        /* read in funny way to accomodate both blank and comma delimiters */

        fgets(inbuf, sizeof(inbuf) - 1, fd);
        end = inbuf + strlen(inbuf);
        *(end + 1) = ' ';
        ptr = inbuf; 
        while (n < 6 * npts && ptr < end)  {
          while (isspace(*ptr)) ptr++;
          sscanf(ptr, "%f", fptr);
/*        printf("%f\n", *fptr);*/
          fptr++;
          n++;
          ptr = strpbrk(ptr + 1, ", \n") + 1; 

        } /* while (i = sscanf(inbuf, "%f", fptr) )  */

/*        printf("%d\n", n);*/
      } /*       while (n < 3 * npts)  */

      *list1 = (XYZPOINT *) malloc (npts * sizeof (XYZPOINT));  
      *list2 = (XYZPOINT *) malloc (npts * sizeof (XYZPOINT));  

      for (i = 0; i < npts; i++) {  /* calculate centroids, load lists  */
      
        for (j = 0; j < 3; j++) {
        
          (*list1)[i][j] = buf1[2*i][j];
          (*list2)[i][j] = buf1[2*i+1][j];
          
        }
        
/*        for (j = 0; j < 3; j++) printf("%.2f,", (*list1)[i][j] );
        for (j = 0; j < 3; j++) printf("%.2f,", (*list2)[i][j] );
        printf("\n");
*/
      }   /*  loop over points - read in and calculate centroid  */
      
      

   }  /* if any points in this hat  */
   
   fclose(fd);
   free(buf1);
   
     return ( npts );

    
  }  /* if opened the file ok  */
  
  else 
    printf("could not open file %s\n", filename);
    return( 0 );
  
}   /* end ReadPointData  */

  
  
static int ReadPointFile(infile, list1, list2)

/*  prompts for file name, reads in  two lists of 3D points into 
two NRUTIL vectors 
(which it also allocates).  returns a pointer to  the file name or NULL if
none.  */

 XYZPOINT **list1, **list2;
 char *infile;

{
  char *filename;  /* buffer for the file name  */
  int i;
  int npoints;
  
  if (infile) return(ReadPointData(infile, list1, list2));

  i = 0;
  while (i == 0) {

    if (!GetaString(&filename, "Name of points file: ")) return 0;  
    i = ReadPointData(filename, list1, list2);
    free(filename);

  } /* while (i == 0)    */  

  return(i);
          
}   /* end ReadPointFile  */




main (argc, argv)

int argc;
char **argv;

{

  XYZPOINT *PointList1, *PointList2;
  XYZPOINT *outpoints;
  int npoints;
  float amat[3][3], scale;
  float hom_mat[4][4], inv_mat[4][4], point4[4], out_point4[4];

  XYZPOINT origin, shift;
  float deltamean[3], delta[3];
  float meandist, thisdist;
  int i, j, doscale = 0;
  int mni=0;
  float mean, ms, rms;
  char *answer;
  char *outfile = NULL;
  char *infile = NULL;
  FILE *fd;
#ifdef THINK_C
  argc = ccommand(&argv);
#endif  
    if ( argc == 1 || !strcmp(argv[1], "-help")) {

      printf("Usage: procrustes inputfile [-scale] [-mni] [-out outfile]\n");
      exit(0);

    } 
   for (i = 1; i < argc; i++) {

    if ( (i == 1) && (argv[i][0] != '-')) {

      infile = (char *)malloc((strlen(argv[i]) + 1) * sizeof(char));
      strcpy(infile, argv[i]);

    } /*   if (!strcmp(argv[i], "-out")  */
    else if (!strcmp(argv[i], "-MNI")  || !strcmp(argv[i], "-mni")) mni = 1;
    /* Alan Evans uses coords with some reflections - mni contr\ols this */
    else if (!strcmp(argv[i], "-scale")) doscale = 1;
    else if ( (i + 1 < argc) && !strcmp(argv[i], "-out")) {

      outfile = (char *)malloc((strlen(argv[i+1]) + 1) * sizeof(char));
      strcpy(outfile, argv[i+1]);

    } /*   if (!strcmp(argv[i], "-out")  */

   } /* for */


  if (! (npoints = ReadPointFile(infile, &PointList1, &PointList2 ))) return(1);
    
  printf("readpointfile returned %d points\n", npoints);

    if (mni) for (i = 0; i < npoints; i++){


       PointList1[i][1] *= -1;
       PointList2[i][1] *= -1;

    } /*     if (mni)  */

     
  ProcrustesMatch(npoints, PointList1, PointList2, origin, amat, shift, 
   &scale, &mean, &ms, &rms, doscale, 1);


  printf("Procrustes returned matrix:\n");
  for (i = 0; i < 3; i++)  {

    for (j = 0; j < 3; j++) printf("%6.2f ", amat[i][j]);
    printf("\n");

  } 
  printf("Procrustes returned origin:\n");

    for (j = 0; j < 3; j++) printf("%6.2f ", origin[j]);
    printf("\n");

  printf("Procrustes returned shift:\n");

    for (j = 0; j < 3; j++) printf("%6.2f ", shift[j]);
    printf("\n");
  
  printf("Procrustes returned scale: %6.3f\n", scale);

  printf("Procrustes returned mean, ms, rms errors: %6.2f, %6.2f, %6.2f\n",
   mean, ms, rms);
 
 outpoints = (XYZPOINT *) malloc (npoints * sizeof(XYZPOINT));

 if (!doscale) scale = 1.0;

  if (outfile) {
    fd = fopen(outfile,"w");
    fprintf(fd, "origin: %.3f %.3f %.3f\n", origin[0], origin[1], origin[2]); 
    fprintf(fd, "shift: %.3f %.3f %.3f\n", shift[0], shift[1], shift[2]); 
    fprintf(fd, "scale: %.3f\n", scale); 
    fprintf(fd, "matrix:\n");

    for (i = 0; i < 3; i++)  {

      for (j = 0; j < 3; j++) fprintf(fd,"%6.2f ", amat[i][j]);
      fprintf(fd, "\n");

    } 

  } /*   if (outfile)  */


  /* finally, here's the result:  a 4x4 matrix that will rotate one of the points
     from the first list (space "1") and match it with the corresponding point
     from the second list (space "2").  Using good computer graphics conventions,
     the vector PREmultiplies the 4x4 matrix:

     (vector in space 2) x (matrix) = (vector in space 1).

  */
   meandist=0.0;
   deltamean[0] = deltamean[1] = deltamean[2] = 0.0;
   
   BuildHomogeneousTransform(origin, amat, shift, scale, hom_mat);
   fprintf(stderr, "using 4x4: \n");
   if (outfile)
       fprintf(fd,"        Original 2           Transformed 1->2              Error\n");
   printf(
       "       Original 2          Transformed 1->2            Error         Dist\n");
   for (i = 0; i < npoints; i++) {
       for (j = 0; j < 4; j++) {
           point4[j] = PointList1[i][j];
       }
       point4[3] = 1.0;

       HMXvectorRotate(point4, hom_mat, out_point4);
    
       for (j = 0; j < 3; j++) {
           delta[j] = PointList2[i][j]- out_point4[j];
       }


       thisdist = DIST3D(delta);
       meandist += thisdist;
       
       printf("(%6.2f,%6.2f,%6.2f)  (%6.2f,%6.2f,%6.2f) (%5.2f,%5.2f,%5.2f)  %5.2f\n", 
              PointList2[i][0], PointList2[i][1], PointList2[i][2],
              out_point4[0], out_point4[1], out_point4[2],
              delta[0], delta[1], delta[2], thisdist);
       
       for(j = 0; j < 3; j++) {
           deltamean[j] += delta[j] * delta[j];
       }
     
       
       if (outfile)
           fprintf(fd,
                   "(%6.2f, %6.2f, %6.2f)  (%6.2f, %6.2f, %6.2f) (%6.2f, %6.2f, %6.2f) %5.2f\n",
                   PointList2[i][0], PointList2[i][1], PointList2[i][2],
                   out_point4[0], out_point4[1], out_point4[2],
                   delta[0], delta[1], delta[2], thisdist);
      
 
    
   }
   for(j = 0; j < 3; j++) {
       deltamean[j] /= npoints;
       deltamean[j] = sqrt(deltamean[j]);
   }
   meandist /= npoints;
  
   printf("  Mean dist: %5.2f                rms (x,y,z): (%5.2f,%5.2f,%5.2f)\n",
          meandist, deltamean[0], deltamean[1], deltamean[2]);

   /*
     and, as my wise friends in computer science have shown me, just by inverting this
     matrix we can transform a point from set 2 to its correct position in space "1".
     
     (vector in space 2) x (matrix) = (vector in space 1).

   */
  
   meandist=0.0;
   deltamean[0] = deltamean[1] = deltamean[2] = 0.0;

  GaussJ(hom_mat, 4, inv_mat);
  fprintf(stderr, "inverting: \n");
    if (outfile)
    fprintf(fd,"        Original 1           Transformed 2->1              Error\n");
     printf(
"       Original 1          Transformed 2->1            Error         Dist\n");
  for (i = 0; i < npoints; i++)
  {
    for (j = 0; j < 4; j++) point4[j] = PointList2[i][j];
    point4[3] = 1.0;


    HMXvectorRotate(point4, inv_mat, out_point4);

    for (j = 0; j < 3; j++) delta[j] = PointList1[i][j]- out_point4[j];

    thisdist = DIST3D(delta);
    meandist += thisdist;

    printf("(%6.2f,%6.2f,%6.2f)  (%6.2f,%6.2f,%6.2f) (%5.2f,%5.2f,%5.2f)  %5.2f\n", 
      PointList1[i][0], PointList1[i][1], PointList1[i][2],
      out_point4[0], out_point4[1], out_point4[2],
      delta[0], delta[1], delta[2], thisdist);
      
      
     for(j = 0; j < 3; j++) deltamean[j] += delta[j] * delta[j];
    if (outfile)
    fprintf(fd,
      "(%6.2f, %6.2f, %6.2f)  (%6.2f, %6.2f, %6.2f) (%6.2f, %6.2f, %6.2f) %5.2f\n",
      PointList1[i][0], PointList1[i][1], PointList1[i][2],
      out_point4[0], out_point4[1], out_point4[2],
      delta[0], delta[1], delta[2], thisdist);
  }
      for(j = 0; j < 3; j++) 
      {
        deltamean[j] /= npoints;
        deltamean[j] = sqrt(deltamean[j]);
      }
      meandist /= npoints;
  
  printf("  Mean dist: %5.2f                rms (x,y,z): (%5.2f,%5.2f,%5.2f)\n",
     meandist, deltamean[0], deltamean[1], deltamean[2]);

HMXprintMatrix(stdout, "Forward (rotates list 1 to match list 2):", hom_mat);
HMXprintMatrix(stdout, "Inverse (rotates list 2 to match list 1):", inv_mat);


  if (outfile)   
  {
    fprintf(fd, "mean, ms, rms errors: %6.2f, %6.2f, %6.2f\n",
   mean, ms, rms);
   HMXprintMatrix(fd, "Forward (rotates list 1 to match list 2):", hom_mat);
   HMXprintMatrix(fd, "Inverse (rotates list 2 to match list 1):", inv_mat);

  }

  if (outfile) fclose(fd);

} /*  */ 
#endif

