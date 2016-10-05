#include <FPT.h>

double gap = 20 ;

void sort(double *x, int n){

  int k,s,i;
  double t;

  for(k=0; k<n; k++){
    s=k;
    for(i = k+1; i<n; i++){
      if(x[i]<x[s]){
	s=i;
      }
    }
    t=x[k];x[k]=x[s];x[s]=t;
  }
}

void my_fill_polygon(double *x, double *y, int numpoints){

  int j; //index of y-point arrays
  double i;  //y value of horizontal line
  int h=0;  //index of sorted xavlue array
  // --double intersection[][2]; //array of points of intersection
  int k=0; //increasing j and necessary to set loop back to 0
  double xvalues[100];

  for(i=0.1; i<600; i++){
  //instead of going from 0 to 600 just find the min value and the max value of y array
	   double slope=0;
	   int t=0;  //index of xvalues array
     for(j=0; j<numpoints;j++){
	//k=(j+1)%numpoints; just another way to set k back to 0 at the end of loop
	     k=j+1; if(k==numpoints) {j=0;}
       if((i>=y[j] && i<=y[k]) || (i<=y[j] && i>=y[k])){
	        slope=((y[j]-y[k])/(x[j]-x[k]));
	        xvalues[t]=(i-y[j]+x[j]*slope)/slope;
	        // --intersection[t][0]=(i-y[j]+x[j]*slope)/slope;
	        // --intersection[t][1]=i; double array not necessary bc y always the same
	        t++;
	     }

     }
   sort(xvalues,t);
    for(h=0;h<t;h++){
      if((h+1)%2==1){
        G_rgb(0,0,1);
        G_line(xvalues[h],i,xvalues[h+1],i);
      }
    }
  }
}


void grid()
{
  int i ;
  for (i = 0 ; i < 600 ; i+= gap) {
    G_line(i,0, i,600) ;
    G_line(0,i, 600,i) ;
  }
}



int click_and_save(double *x, double *y)
{
  double xy[2] ;
  int numpoints ;
  double xc,yc ;

  G_rgb(1, 0, 0) ; // red
  G_fill_rectangle(0,0,  20,10) ;

  numpoints = 0 ;
  while (0 == 0) {

     G_wait_click(xy) ;

     if ((xy[0] >= 0) && (xy[0] <= 20) &&
         (xy[1] >= 0) && (xy[1] <= 10))  { break ; }

     G_rgb(1, 1, 0) ; // yellow


     xc = gap*floor((xy[0]+0.5*gap)/gap) ;
     yc = gap*floor((xy[1]+0.5*gap)/gap) ;

     G_circle(xc, yc, 3) ;


     x[numpoints] = xc ; y[numpoints] = yc ;

     if (numpoints > 0) {
       G_line (x[numpoints-1],y[numpoints-1], x[numpoints],y[numpoints]) ;
     }

     numpoints++ ;

  }

  return numpoints ;
}







int  main()
{
  int q ;
  double xp[1000],yp[1000] ;
  int np ;
  int k ;

  G_init_graphics (600,600) ;

  do {

     G_rgb(0,0,0) ;
     G_clear() ;
     G_rgb(1,0.5,0) ;
     grid() ;

     np = click_and_save(xp,yp) ;
     G_rgb(1,0,0) ;
     my_fill_polygon(xp,yp,np) ;

     G_rgb(0,1,0) ;
     for (k = 0 ; k < np ; k++) {
      G_circle(xp[k],yp[k],2) ;
     }

     q = G_wait_key() ;

  } while (q != 'q') ;


}
