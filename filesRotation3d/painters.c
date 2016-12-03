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



void draw_file(int onum){

	int i,j,k;
	double x_current[10000], y_current[10000], z_current[10000];

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

		//compute the distance of poly and then put them into an array at this point

		if(onum%2 == 0){
			G_rgb(0,0,1);
		} else {
			G_rgb(1,0,0);
		}

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



