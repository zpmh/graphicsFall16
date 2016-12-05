#include <FPT.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <D3d_matrix.h>

int numobjects;
int numpoints[10];
double x[10][10000], y[10][10000], z[10][10000];
int numpolys[10000];
int psize[10][10000];
int con[10][10000][200];
double red[10][5000], grn[10][5000], blu[10][5000];
double cx,cy,cz;

double lightx, lighty,lightz;
double rinit, ginit, binit;
double ambient, maxdiff, specpow ; 

typedef
struct {
  int objnum ;
  int polynum ;
  double dist ;
}
  THING;

THING allpolys[10000];

void print_array(int n) {
  int i ;
  for (i = 0 ; i < n ; i++) {
    //printf("%d %d %lf\n",allpolys[i].objnum, allpolys[i].polynum, allpolys[i].dist) ;
  }
  //printf("\n") ;
}

int compare (const void *p, const void *q) {
  THING *a, *b ;

  a = (THING*)p ;
  b = (THING*)q ;

  if  (((*a).dist) < ((*b).dist)) return -1 ;
  else if (((*a).dist) > ((*b).dist)) return 1 ;
  else return 0;
}


int init_array() {  

  int allpolysN = 0;
  int i;
  int j;

  //printf("numobjects = %d\n",numobjects) ;

  for (i = 0; i < numobjects; i++) { //COUNT TOTAL NUMBER OF POLYGONS IN SORT and INITIALIZE 
    for (j = 0; j < numpolys[i]; j++) {
      allpolys[allpolysN].objnum = i;
      allpolys[allpolysN].polynum = j;
      allpolys[allpolysN].dist = z[i][con[i][j][0]];
      allpolysN = allpolysN + 1;
    }  
  }

  qsort (allpolys, allpolysN, sizeof(THING), compare);
  //printf("hello  %d\n",allpolysN) ;

  print_array(allpolysN);

  return allpolysN;
}



//double theta = 20*M_PI/180;
//double H;
//H = tan(theta);  // ???Data definition has no type or storage class -- not constant??



//COULD USE A STRUCT FOR ALL OF THESE 
int crossprod(double *v,double *u, double w[3]){

  //This computes the cross product of two input vectors

  w[0]= v[1]*u[2] - v[2]*u[1];
  w[1]= v[2]*u[0] - v[0]*u[2];
  w[2]= v[0]*u[1] - v[1]*u[0];
	
  return 1;

} 

int vecs(double *v,double *u,double w[3]){

  //Finds vector within each polygon to find the Normal vector etc - needs to be called twice to find two
	
  w[0]= v[0]-u[0];
  w[1]= v[1]-u[1];
  w[2]= v[2]-u[2];

  if(w[0] == w[1] && w[1] == w[2]){ 
	printf("scalar multiples");
	return -1;
  }
	
  return 1;
}

double vecsize(double *v){

  return sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);

}

double dotprod(double *v,double *u){

  return u[0]*v[0]+u[1]*v[1]+u[2]*v[2];

}





double xmax(int onum){
  int i;
  double max= x[onum][0];
  for(i=0;i<numpoints[onum];i++){
    if(x[onum][i] > max){
      max = x[onum][i];
    }
  }
  return max;
}

double xmin(int onum){
  int i;
  double min= x[onum][0];
  for(i=0;i<numpoints[onum];i++){
    if(x[onum][i] < min){
      min = x[onum][i];
    }
  }
  return min;
}
double ymax(int onum){
  int i;
  double max= y[onum][0];
  for(i=0;i<numpoints[onum];i++){
    if(y[onum][i] > max){
      max = y[onum][i];
    }
  }
  return max;
}

double ymin(int onum){
  int i;
  double min= y[onum][0];
  for(i=0;i<numpoints[onum];i++){
    if(y[onum][i] < min){
      min = y[onum][i];
    }
  }
  return min;
}

double zmax(int onum){
  int i;
  double max= z[onum][0];
  for(i=0;i<numpoints[onum];i++){
    if(z[onum][i] > max){
      max = z[onum][i];
    }
  }
  return max;
}

double zmin(int onum){
  int i;
  double min= z[onum][0];
  for(i=0;i<numpoints[onum];i++){
    if(z[onum][i] < min){
      min = z[onum][i];
    }
  }
  return min;
}

int lightcolor(THING allpolys,double *shade){

  double intensity;
  double cosa,cosb;
  double Nu[3], Lu[3], Ru[3], Eu[3]; 
  double v[3], u[3];

  //find v and u inside polygon:

  double p[3], q[3], r[3];

  p[0]= x[allpolys.objnum][con[allpolys.objnum][allpolys.polynum][0]];
  p[1]= y[allpolys.objnum][con[allpolys.objnum][allpolys.polynum][0]];
  p[2]= z[allpolys.objnum][con[allpolys.objnum][allpolys.polynum][0]];

  q[0]= x[allpolys.objnum][con[allpolys.objnum][allpolys.polynum][1]]; 
  q[1]= y[allpolys.objnum][con[allpolys.objnum][allpolys.polynum][1]];
  q[2]= z[allpolys.objnum][con[allpolys.objnum][allpolys.polynum][1]]; 

  r[0]= x[allpolys.objnum][con[allpolys.objnum][allpolys.polynum][2]]; 
  r[1]= y[allpolys.objnum][con[allpolys.objnum][allpolys.polynum][2]];
  r[2]= z[allpolys.objnum][con[allpolys.objnum][allpolys.polynum][2]]; 

  vecs(p,q,v);
  vecs(p,r,u);

	// check that not scalar multiples - if vecs = -1 do something? 

  crossprod(v, u, Nu);

  Nu[0] = Nu[0] / vecsize(Nu) ; 
  Nu[1] = Nu[1] / vecsize(Nu) ; 
  Nu[2] = Nu[2] / vecsize(Nu) ;

  Lu[0] = lightx - p[0] / vecsize(Nu);
  Lu[1] = lighty - p[1] / vecsize(Nu);
  Lu[2] = lightz - p[2] / vecsize(Nu);

  if(dotprod(Nu,Lu) < 0){	Nu[0] = -1*Nu[0]; Nu[1] = -1*Nu[1]; Nu[2] = -1*Nu[2]; }

  Ru[0] = (2*dotprod(Nu,Lu)*Nu[0] / vecsize(Nu));
  Ru[1] = (2*dotprod(Nu,Lu)*Nu[1] / vecsize(Nu));
  Ru[2] = (2*dotprod(Nu,Lu)*Nu[2] / vecsize(Nu));

  Eu[0] = p[0] / vecsize(Nu); 
  Eu[1] = p[1] / vecsize(Nu);
  Eu[2] = p[2] / vecsize(Nu);

  if(dotprod(Nu,Eu) < 0){ 
	shade[0] = ambient/(ambient+maxdiff)*rinit; 
	shade[1] = ambient/(ambient+maxdiff)*ginit; 
	shade[2] = ambient/(ambient+maxdiff)*binit; 
	return 1; } // skip the rest and go to the next polygon !!!

  cosa = dotprod(Nu,Lu);
  cosb = dotprod(Eu,Ru);

  if(cosb<0){
	cosb=0;
  }

  intensity = ambient + maxdiff*cosa + (1-ambient-maxdiff)*pow(cosb,specpow);

  shade[0] = (intensity/(ambient+maxdiff))*rinit;    
  shade[1] = (intensity/(ambient+maxdiff))*ginit;
  shade[2] = (intensity/(ambient+maxdiff))*binit;
	
  return 1;
}



void draw_file(int onum){

  int i,j,k;
  double x_current[10000], y_current[10000], z_current[10000];
  double light[3];

  int allpolysN = init_array();

  for(k=allpolysN; k>=0; k--){
    onum = allpolys[k].objnum;
    i = allpolys[k].polynum;	
	

    //for(i=0; i<numpolys[onum];i++){
    for(j=0; j<psize[onum][i]; j++){
      // finding the original values of x, y and z points
      x_current[j] = x[onum][con[onum][i][j]];
      y_current[j] = y[onum][con[onum][i][j]];
      z_current[j] = z[onum][con[onum][i][j]];			 

      // find the distance each point has to projection screen
      y_current[j] = y_current[j] / z_current[j];
      x_current[j] = x_current[j] / z_current[j];
      z_current[j] = 1;  
	
      //translate the H,H perspective to actual graphics window
      x_current[j] = (300/tan(20*M_PI/180))*x_current[j] + 300;
      y_current[j] = (300/tan(20*M_PI/180))*y_current[j] + 300;

      //you could do the two steps above in one step...
    } // end for j


    //Instead of changing the color for each object here, you should call the light function to change each polgons color!!!
    //give it a polygon

    lightcolor(allpolys[k], light);
    G_rgb(light[0],light[1],light[2]);		

    /*
      if(onum%2 == 0){
      G_rgb(0,0,1);
      } else {
      G_rgb(1,0,0);
      }
    */

    //G_rgb(red[onum][i],grn[onum][i],blu[onum][i]);

    G_fill_polygon(x_current,y_current,psize[onum][i]);
    G_rgb(96,96,96);
    G_polygon(x_current, y_current, psize[onum][i]);
  }
}

void center(onum){

  double m[4][4];
  double minv[4][4];

  D3d_make_identity(m);
  D3d_make_identity(minv);

  double xBig =  xmax(onum);
  double xSmall = xmin(onum);
  double yBig =  ymax(onum);
  double ySmall = ymin(onum);
  double zBig = zmax(onum);
  double zSmall = zmin(onum);


  double xCenter = (xBig + xSmall) / 2 ;
  double yCenter = (yBig + ySmall) / 2;
  double zCenter = (zBig + zSmall) / 2;

  D3d_translate(m,minv,-xCenter, -yCenter, -zCenter);

  double sf ;
  if(xBig-xSmall > yBig-ySmall && xBig-xSmall > zBig-zSmall){
    sf = 1/(xBig-xSmall);  // it's 1 instead of 300 bc draw is centering it at 300,300
  } else if(yBig-ySmall > xBig-xSmall && yBig-ySmall > zBig-zSmall) {  
    sf = 1/(yBig-ySmall) ; // fits it into 1 by 1 box
  } else {
    sf = 1/(zBig-zSmall);
  }

  D3d_scale(m,minv,sf,sf,sf); 
  D3d_translate(m,minv,0,0,20);
  D3d_mat_mult_points(x[onum],y[onum],z[onum],m,x[onum],y[onum],z[onum],numpoints[onum]);
}


void objCenter(onum){

  double m[4][4];
  double minv[4][4];
  double xSum, ySum, zSum;
  int i;

  D3d_make_identity(m);
  D3d_make_identity(minv);

  xSum = ySum = zSum = 0 ;

  for(i=0; i<numpoints[onum]; i++){
    xSum += x[onum][i];
    ySum += y[onum][i];
    zSum += z[onum][i];
  }

  cx = xSum / numpoints[onum];
  cy = ySum / numpoints[onum];
  cz = zSum / numpoints[onum];

} 

int main(int argc, char **argv){

  FILE *f ;
  int a,b,c,d,e,g,i,j,onum ;

  double m[4][4];
  double minv[4][4];
  onum = 0;
  i = 0;

  //let user input values for LIGHT etc
  printf("Please input x,y,z coordinates for the light source \n");
  printf("x: \n");
  scanf("%lf", &lightx);
  printf("y: \n");
  scanf("%lf", &lighty);
  printf("z: \n");
  scanf("%lf", &lightz);

  //input values for R,G,B
  printf("Please enter RGB values for color all in one line with a space inbetween\n");
  scanf("%lf %lf %lf", &rinit, &ginit, &binit);

  printf("Please input values for Ambient, MaxDiffuse and Specpow \n");
  printf("Ambient: \n");
  scanf("%lf", &ambient);
  printf("MaxDiffuse: \n");
  scanf("%lf", &maxdiff);
  printf("Specpow: \n");
  scanf("%lf", &specpow);



  //initialize Graphics

  G_init_graphics(600,600);
  G_rgb(0,0,0);
  G_clear();	

  numobjects = argc-1;
    
  //attempt to open file
  for(i=0; i<numobjects; i++){

    f= fopen(argv[i+1], "r");
    if(f == NULL) {
      printf("can't open file\n");
      exit(0);
    }

    //read file
	
    fscanf(f,"%d",&numpoints[i]);

    for(a=0; a<numpoints[i];a++){
      fscanf(f, "%lf %lf %lf", &x[i][a], &y[i][a], &z[i][a]);
    }

    fscanf(f,"%d",&numpolys[i]);

    for(b=0; b<numpolys[i]; b++){
      fscanf(f, "%d", &psize[i][b]);	
      for(c=0; c < psize[i][b]; c++){
	fscanf(f, "%d", &con[i][b][c]);	
      }
    }
    center(i);
  }



  //     draw_file(0) ;    G_wait_key() ;    exit(0) ;

  //draw stage

  char v = '0' ;
  int go = 0; //go is the onum when you press keys 0 1 2 so on == onum

     
  //q = quit

  while(v != 'q'){

    G_rgb(0,0,0);
    D3d_make_identity(m);
    D3d_make_identity(minv);
    G_clear();
    objCenter(go);
    //find the bounding box!!
    if( v == 'x'){   // press r to rotate
      D3d_translate(m,minv, -cx, -cy,-cz);
      D3d_rotate_x(m,minv, .05);
      D3d_translate(m,minv, cx, cy, cz);}
    else if( v == 'y'){   // press r to rotate
      D3d_translate(m,minv, -cx, -cy,-cz);
      D3d_rotate_y(m,minv, .05);
      D3d_translate(m,minv, cx, cy, cz);}
    else if( v == 'z'){   // press r to rotate
      D3d_translate(m,minv, -cx, -cy,-cz);
      D3d_rotate_z(m,minv, .05);
      D3d_translate(m,minv, cx, cy, cz);} 
    else if( v == 't'){ // press q to move object forward x direction
		  	// in draw z[onum][i] !=1   draw_file(go,0.5);
      D3d_translate(m,minv,1,0,0);}
    else if( v == 'w'){ // forward in y direction
      D3d_translate(m,minv,0,1,0);}
    else if( v == 'e'){ //backward in z direction
      D3d_translate(m,minv,0,0,1);  
      z[go][numpoints[go]] += c; // WHY??????? Need to keep track of c for the 20 in rotate
    }
    else if( v == 'g'){ // backward in x direction
      D3d_translate(m,minv,-1,0,0);}
    else if( v == 's'){ // backward in y direction
      D3d_translate(m,minv,0,-1,0);}
    else if( v == 'd'){ // backward in z direction
      D3d_translate(m,minv,0,0,-1);}
    else{
      go = v - 48;  //otherwise press 0 1 2 to get other pic and translate - 48 ASCII-code 
    }

    D3d_mat_mult_points(x[go],y[go],z[go],m,x[go],y[go],z[go],numpoints[go]);
    draw_file(go);
    //printf("%c \n", v); // seg fault for j
    v = G_wait_key();
  }
}


