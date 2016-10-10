#include <FPT.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

int numobjects;
int numpoints[10];
double x[10][1000], y[10][1000];
int numpolys[10];
int psize[10][1000];
int con[10][5000][20];
double red[10][5000], grn[10][5000], blu[10][5000];






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


void center(onum){

    double xBig =  xmax(onum);
    double xSmall = xmin(onum);
    double yBig =  ymax(onum);
    double ySmall = ymin(onum);

    double xCenter = (xBig + xSmall) / 2 ;
    double yCenter = (yBig + ySmall) / 2;

    translate(onum,-xCenter, -yCenter);
    if(xBig-xSmall > yBig-ySmall){
    	scale(onum,300/(yBig-ySmall),300/(yBig-ySmall));
    } else {  
    	scale(onum,300/(xBig-xSmall),300/(xBig-xSmall));
    }
    translate(onum, 300, 300);
}

int main(int argc, char **argv){

    FILE *f ;
    int a,b,c,d,e,g,i,j ;

    //initialize Graphics

    G_init_graphics (600,600);
    G_rgb(0,0,0);
    G_clear();	
	
    for(i=0; i<argc-1; i++){

		f=fopen (argv[i+1], "r");
		if(f == NULL) {
	  		printf("can't open file\n");
	  		exit(0);
		}

	//read file
	
		fscanf(f,"%d",&numpoints[i]);

		for(a=0; a<numpoints[i];a++){
			fscanf(f, "%lf %lf", &x[i][a], &y[i][a]);
		}

		fscanf(f,"%d",&numpolys[i]);
		printf("%d", numpolys[i]);

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
		draw_file(go);
		v = G_wait_key();
    }
}

