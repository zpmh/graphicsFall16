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



void draw_file(int onum){

	int i,j;
	double x_current[1000], y_current[1000];
    	for(i=0; i<numpolys[onum];i++){
       		for(j=0; j<psize[onum][i]; j++){
         		x_current[j] = x[onum][con[onum][i][j]];
	 			y_current[j] = y[onum][con[onum][i][j]];
			} 
       	 	G_rgb(red[onum][i],grn[onum][i],blu[onum][i]);
		printf("%lf %lf %lf\n", red[onum][i],grn[onum][i],blu[onum][i]);
       	 	G_fill_polygon(x_current,y_current,psize[onum][i]);
		printf(" %d\n", psize[onum][i]);
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
    printf("x range : %lf %lf\n",xSmall,xBig) ;
    printf("y range : %lf %lf\n",ySmall,yBig) ;

    double xCenter = (xBig + xSmall) / 2 ;
    double yCenter = (yBig + ySmall) / 2;

    D2d_translate(m,minv,-xCenter, -yCenter);

    double sf ;
    if(xBig-xSmall > yBig-ySmall){
        sf = 300/(xBig-xSmall);
    } else {  
        sf = 300/(yBig-ySmall) ;
    }
    printf("sf = %lf\n",sf) ;

    D2d_scale(m,minv,sf,sf) ;
    D2d_translate(m,minv, 300, 300);
    D2d_mat_mult_points(x[onum],y[onum],m,x[onum],y[onum],numpoints[onum]);
}

int main(int argc, char **argv){

    FILE *f ;
    int a,b,c,d,e,g,i,j,onum ;

    double m[3][3];
    double minv[3][3];
    onum = 0;
    i = 0;

    //initialize Graphics

    G_init_graphics(600,600);
    G_rgb(0,0,0);
    G_clear();	
    
    //attempt to open file
    for(i=0; i<argc-1; i++){

		f= fopen(argv[i+1], "r");
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



    //    draw_file(0) ;    G_wait_key() ;    exit(0) ;

     //draw stage

     char v = '0' ;
     int go = 0; //go is the onum when you press keys 0 1 2 so on == onum

     
     //q = quit

    while(v != 'q'){
       	G_rgb(0,0,0);
        D2d_make_identity(m);
	D2d_make_identity(minv);
	G_clear();
		if( v == 'r'){   // press r to rotate
	  		D2d_translate(m,minv, -300, -300);
	  		D2d_rotate(m,minv, .05);
	  		D2d_translate(m,minv, 300, 300);
			D2d_mat_mult_points(x[go],y[go],m,x[go],y[go],numpoints[go]);
		} else{
	  		go = v - 48;  //otherwise press 0 1 2 to get other pic and translate - 48 ASCII-code 
		}
		draw_file(go);
		printf("%c \n", v); // seg fault for j
		v = G_wait_key();
    }
}

