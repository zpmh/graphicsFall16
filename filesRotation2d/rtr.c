#include <FPT.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <D2d_matrix.h>

int numobjects;
int numpoints[10];
double x[10][1000], y[10][1000];
int numpolys[10];
int psize[10][1000];
int con[10][5000][20];
double red[10][5000], grn[10][5000], blu[10][5000];


// clipping window
int nw ;
double xw[100],yw[100] ;

double X[1000],Y[1000];







int click_and_save(double *xclip, double *yclip)
{
  double xy[2] ;
  int cpoints ;

  G_rgb(1, 0, 0) ; // red
  G_fill_rectangle(0,0,20,10) ; //lower conrner rectangle

  G_wait_click(xy);

  cpoints = 0 ;
  while (0 == 0) {

     printf("cpoints = %d\n",cpoints) ;

     if ((xy[0] >= 0) && (xy[0] <= 20) &&
         (xy[1] >= 0) && (xy[1] <= 10))  { break ; }  

     xclip[cpoints] = xy[0] ; yclip[cpoints] = xy[1] ;

     G_rgb(1,0,0);
     G_circle(xclip[cpoints],yclip[cpoints],3);


     if (cpoints > 0) {
       G_line (xclip[cpoints-1],yclip[cpoints-1], xclip[cpoints],yclip[cpoints]) ;
     }


     cpoints++ ;
     G_wait_click(xy);

  }

  return cpoints ;

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


void translate(int onum, double dx, double dy){
	int i;
	for(i=0; i<numpoints[onum]; i++){
		x[onum][i] += dx;
		y[onum][i] += dy;
	}
}

void scale(int onum, double sx, double sy){
	int i;
	for(i=0;i<numpoints[onum];i++){
		x[onum][i] *=sx;
		y[onum][i] *=sy;	
	}
}

void rotate(int onum, double t){
    //in radians
	int i;
	double temp=0;
	double s= sin(t);
	double c= cos(t);
	for(i=0;i<numpoints[onum];i++){
		temp = x[onum][i] * c - y[onum][i] * s;
		y[onum][i] = x[onum][i]*s + y[onum][i]*c;
		x[onum][i] = temp;
	}
}


//onum = file number

void draw_file(int onum){

	int i,j;
	double x_current[1000], y_current[1000];
    	for(i=0; i<numpolys[onum];i++){
       		for(j=0; j<psize[onum][i]; j++){
         		x_current[j] = x[onum][con[onum][i][j]];
	 		y_current[j] = y[onum][con[onum][i][j]];
			} 
       	 	G_rgb(red[onum][i],grn[onum][i],blu[onum][i]);
       	 	G_fill_polygon(x_current,y_current,psize[onum][i]);
    	}
}

/*int in_out(int n, double XY[2]){

	int i, j, flag;
	double CM[2],dx,dy,a,b,c,u,v,sum1,sum2;

	sum1 = sum2 = 0;

	//find center of mass of convex clipping window
	CM[0] = CM [1] = 0;
	for(i=0; i<nw; i++){
		sum1 += xw[i];
		sum2 += yw[i];
	}
	CM[0] = sum1; CM[1] = sum2;
	CM[0] /= nw; CM[1] /= nw;

	flag = 1;
	
	j = (n+1) % nw;

	dx = xw[n] - xw[j];
	dy = yw[n] - yw[j];
	a = dy ;
  	b = -dx ;
  	c = xw[n]*yw[j] - xw[j]*yw[n] ;   

  	u = a*XY[0] + b*XY[1] + c ;
  	v = a*CM[0] + b*CM[1] + c ;

  	if (u*v < 0) { flag = 0 ; }
 

  	return flag ;

} */

int in_out( double x, double y, int iter, int iter1, int iter2){

    if((xw[iter1]-xw[iter])*yw[iter2] + (yw[iter] - yw[iter1])*xw[iter2] + (xw[iter]*yw[iter1] - yw[iter]*xw[iter1]) >= 0){

          if((xw[iter1]-xw[iter])*y + (yw[iter]-yw[iter1])*x+(xw[iter]*yw[iter1]-yw[iter]*xw[iter1]) >= 0){
      return 1;
    }

    } else {
    if((xw[iter1]-xw[iter])*y + (yw[iter]-yw[iter1])*x+(xw[iter]*yw[iter1]-yw[iter]*xw[iter1]) <= 0){
      return 1;
    }

    }

    return 0;

}


int clip(double *x_current, double *y_current, int ppts, double *Xnew, double *Ynew){

// initialize variables :
  int a, c, d, f, p, s, pt, n, n2, count = 0;
  int INOUT1, INOUT2 ;
  double Xint, Yint, aXY[2], dXY[2] ;
  double A, B, C, D, E, F ;
  double DX, DY, DWX, DWY ;
  double Xtemp[100], Ytemp[100] ;



  // Copy Polygon to Xnew, Ynew & Xtemp, Ytemp:
  for ( a=0 ; a<ppts ; a++ ) {
    Xnew[a] = x_current[a] ;
    Ynew[a] = y_current[a] ;
    Xtemp[a] = x_current[a] ;
    Ytemp[a] = y_current[a] ;
printf("X = %lf \t Y = %lf \n", x_current[a], y_current[a]);
  }


  // Cycle Clipping Window Lines :
  for ( s=0 ; s<nw ; s++ ) {
    n = (s+1) % nw;
    n2 = (s+2) % nw;

    // Intersection Prep :
    f = (s+1) % nw ;
    DWX = xw[s] - xw[f] ;
    DWY = yw[s] - yw[f] ;
    D =  DWY ;
    E = -DWX ;
    F =  xw[s]*yw[f] - xw[f]*yw[s] ;
    //printf("Window line %d to %d!\n", s, f) ;


    // Cycle Points (ppts is modified later) :
    for ( a=0, pt=0 ; a<ppts ; a++ ) { 
      
      // Initializing for later :
      d = (a+1) % ppts ;
      DX = Xnew[a] - Xnew[d] ;
      DY = Ynew[a] - Ynew[d] ;
      A =  DY ;
      B = -DX ;
      C =  Xnew[a]*Ynew[d] - Xnew[d]*Ynew[a] ;

      // IN OUT testing
      aXY[0] = Xnew[a] ;  aXY[1] = Ynew[a] ;
      dXY[0] = Xnew[d] ;  dXY[1] = Ynew[d] ;
      /*INOUT1 = in_out(s, aXY ) ;
      INOUT2 = in_out(s, dXY ) ; */
      INOUT1 = in_out(Xnew[a], Ynew[a], s, n, n2) ;
      INOUT2 = in_out(Xnew[d], Ynew[d], s, n, n2) ;
	
	printf("inout1 = %d \t inout2 = %d\n", INOUT1, INOUT2);
     count++;
      // IN-->IN
      if ( INOUT1+INOUT2 == 2 ) {
	Xtemp[pt] = Xnew[d] ;
	Ytemp[pt] = Ynew[d] ;
	pt++ ;
printf("in to in\n");
      }

      // IN-->OUT or OUT-->IN
      else if ( INOUT1+INOUT2 == 1) {
	
	// Intersection Math :
	Xint = (B*F - C*E) / (A*E - B*D) ;
	Yint = (C*D - A*F) / (A*E - B*D) ;
	// Note : AE - BD != 0 because lines ! parallel

	// Record intersection first in either case
        Xtemp[pt] = Xint ;
	Ytemp[pt] = Yint ;
	pt++ ;
printf("in to out or out to in\n");
	if ( INOUT2 == 1 ) {
	  // OUT-->IN additionally requires :
	  Xtemp[pt] = Xnew[d] ;
	  Ytemp[pt] = Ynew[d] ;
	  pt++ ;
printf("out to in\n");
	}
      }
      
    } // exit inner loop


    // Copy from temp to new poly info
    for(c=0 ; c<pt ; c++) {
      Xnew[c] = Xtemp[c] ;
      Ynew[c] = Ytemp[c] ;
    }

    // Modify number of points in new poly (used in loop)
    ppts = pt ;  

  }  // exit outer loop
 
  return ppts ;
}

void draw_clip(int onum){

	int i,j;
	int cpp;
	double x_current[1000], y_current[1000];
    	for(i=0; i<numpolys[onum];i++){
       		for(j=0; j<psize[onum][i]; j++){
         		x_current[j] = x[onum][con[onum][i][j]];
	 		y_current[j] = y[onum][con[onum][i][j]];
			} 
		cpp = clip(x_current,y_current,j,x_current,y_current);  // XW and YW
//printf("cpp = %d\n",cpp);
       	 	G_rgb(red[onum][i],grn[onum][i],blu[onum][i]);
//G_rgb(0,0,0);
//G_clear();
//G_rgb(0,0,1);
       	 	G_fill_polygon(x_current,y_current,cpp);
		printf(" %d \n", cpp);
    	}
}


void center(onum){

    double m[3][3];
    double minv[3][3];

    D2d_make_identity(m);
    D2d_make_identity(minv);

    double xBig =  xmax(onum);
    double xSmall = xmin(onum);
    double yBig =  ymax(onum);
    double ySmall = ymin(onum);

    double xCenter = (xBig + xSmall) / 2 ;
    double yCenter = (yBig + ySmall) / 2;

    D2d_translate(m,minv,-xCenter, -yCenter);

    double sf ;
    if(xBig-xSmall > yBig-ySmall){
        sf = 300/(xBig-xSmall);
    } else {  
        sf = 300/(yBig-ySmall) ;
    }

    D2d_scale(m,minv,sf,sf) ;
    D2d_translate(m,minv, 300, 300);
    D2d_mat_mult_points(x[onum],y[onum],m,x[onum],y[onum],numpoints[onum]);
}

int main(int argc, char **argv){

    FILE *f ;
    int a,b,c,d,e,g,i,j ;

    //initialize Graphics
	
    for(i=0; i<argc-1; i++){

		f=fopen (argv[i+1], "r");
		if(f == NULL) {
	  		printf("can't open file %s\n", argv[i+1] );
	  		exit(0);
		}

	//read file
	
		fscanf(f,"%d",&numpoints[i]);

		for(a=0; a<numpoints[i];a++){
			fscanf(f, "%lf %lf", &x[i][a], &y[i][a]);
		}

		fscanf(f,"%d",&numpolys[i]);
		//printf("%d", numpolys[i]);

		for(b=0; b<numpolys[i]; b++){
			fscanf(f, "%d", &psize[i][b]);	
			for(c=0; c < psize[i][b]; c++){
		  		fscanf(f, "%d", &con[i][b][c]);	
			}
		}
		for(d=0; d<numpolys[i];d++){
			fscanf(f, "%lf %lf %lf", &red[i][d], &grn[i][d], &blu[i][d]);
		}
		center(i);
    }

     //draw stage





    G_init_graphics (600,600);
    G_rgb(0,0,0);
    G_clear();	

  //  draw_file(0) ;


    nw = click_and_save(xw,yw) ;
 //   G_rgb(1,0,0) ;
 //   G_polygon(xw,yw,nw) ;
 //   G_wait_key() ;


   // draw_clip(0) ;

   // G_wait_key() ;

    // exit(1) ;


     char v = 'a';
     int go = 0;
     
     //q = quit

    while(v != 'q'){
       	G_rgb(0,0,0);
		G_clear();
		if( v == 'r'){
	  		translate(go, -300, -300);
	  		rotate(go, .5);
	  		translate(go, 300, 300);
		} else{
	  		go = v - 48;
		}
		draw_clip(go);
G_rgb(1,1,1);
G_polygon(xw, yw, nw);
		v = G_wait_key();
    }
}
