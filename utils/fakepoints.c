#include <stdio.h>
#include <stdlib.h>

typedef float XYZPOINT[3];
main()

{

  int i,j,k;

  float inpoints[][3] = { 10,0,0, 0,10,0, 0,0,10, -10,0,0, 0,-10,0, 0,0,-10,
                          10,10,10, -10,10,10, -10,-10,10, -10,-10,-10};

  int numpoints = sizeof(inpoints) / sizeof(inpoints[0]);
  
  XYZPOINT *outpoints;
  
  float amat[3][3];
  float rx, ry, rz;
  float *list;
  float origin[] = {5.0,2.0,-1.0};
  
  FILE *fd;
  
  outpoints = (XYZPOINT *) malloc(numpoints * sizeof(XYZPOINT));
  
  for (i = 0; i < numpoints; i++) {
    for (j = 0; j < 3; j++) printf("%6.2f ", inpoints[i][j]);
    printf("\n");
  }
  
  printf("enter Rz, Ry, Rx:"); scanf("%f %f %f", &rz, &ry, &rx);
 printf("%f %f %f\n", rz, ry, rx); 
  getmatrix(rz, ry, rx, amat); 
  
 fprintf(stderr, "rotation matrix:\n");
     for (i = 0; i < 3; i++)  {
    
      for (j = 0; j < 3; j++) fprintf(stderr,  "%6.3f ", amat[i][j]); 
         fprintf(stderr,  "\n");
   }
  rotate_about_center(numpoints, inpoints, amat, origin, outpoints, 1);
  printf("back from rotate\n");
  if (!(fd = fopen("points.fake", "w") ) ) exit(0);
    
  fprintf(fd, "%d\n", numpoints);
  
     for (i = 0; i < numpoints; i++)  {
    
      for (j = 0; j < 3; j++) fprintf(fd, "%6.3f,", outpoints[i][j]); 
       for (j = 0; j < 2; j++) fprintf(fd, "%6.3f,", inpoints[i][j]); 
     
          fprintf(fd, "%6.3f\n", inpoints[i][2]);
          
    }

 fclose(fd);
 
 }
