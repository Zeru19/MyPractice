#include<stdio.h>
#include<math.h>
#include<complex.h>
#include<malloc.h>
#include<string.h>
#include<stdlib.h>
#define PI 3.14159
#define INF 10000000 

/*double h_bar=1.0546E-34*/
double h_bar=1;
double m=1.,k=1;

_Complex double* new_array1(int length){
	_Complex double *array=(_Complex double*)malloc(length*sizeof(_Complex double));
	return array;
} 

_Complex double** new_array2(int row, int col)
{
	int i = 0;
	_Complex double **array = (_Complex double**)malloc(sizeof(_Complex double*) * row);
	for (i = 0; i < row; ++i) {
		array[i] = (_Complex double*)malloc(sizeof(_Complex double) * col);
		memset(array[i], 0, sizeof(_Complex double) * col);
	}
	return array;
}

void debug(){
	printf("has been done\n");
}

_Complex double power_function(double x){
	double y;
	_Complex double z;
	y=20.73629*pow(x,8)-216.04999*pow(x,7)+974.40772*pow(x,6)-2483.84505*pow(x,5)+
		3916.20660*pow(x,4)-3917.27332*pow(x,3)+2434.62330*pow(x,2)-863.07535*pow(x,1)-
		707.17883;
	z=y+0*I;
	return z;
}

double sum(_Complex double *a,int n){
	int i;
	double sum=0;
	for(i=0;i<n;i++){
		sum+=cabs(a[i])*cabs(a[i]);
	}
	return sum;
}

int main(){
	int i,j,z,c,N,run_time=1000;
	//time interval and the division of space
	double delta_t=0.0001,delta_h=0.005;
	double x_min,x_max;
	//the range of x
	x_min=0.7813;x_max=1.823;
	for(i=0;x_min+i*delta_h<=x_max;i++){
		N++;
	}
	
	//input the potential function
	//Demo potential function
	_Complex double *V=new_array1(N);
	for(i=0;i<N;i++){
		V[i]=power_function(x_min+delta_h*i);
	}
	//V[0]=INF;V[N-1]=INF;
	//initial wave function
	//Demo function
	_Complex double **Psi=new_array2(N,run_time);
	double x=x_min;
	for(i=0;i<N;i++){
		Psi[i][0]=1./(pow(PI*0.03,1./4))*exp(-1*(x-1.2)*(x-1.2)/(2*0.03))*(cos(k*x)+I*sin(k*x));
		x=x+delta_h;
	}
	for(j=0;j<run_time;j++){
		Psi[0][j]=0;
		Psi[N-1][j]=0;
	}
	//the coefficients 
	_Complex double A=I*delta_t/(4*pow(delta_h,2));
	_Complex double *B=new_array1(N);
	for(i=0;i<N;i++){
		B[i]=I*delta_t*V[i]/2.;
	}
	//to store Psi when time= 't and t+1'
	_Complex double *u=new_array1(N-1);
	_Complex double *v=new_array1(N-1);
	//triple diagonal matrix 0(main), 1(up), negetive 1(under)
	_Complex double *dia0=new_array1(N-1);
	_Complex double *dia1=new_array1(N-1);
	_Complex double *dian1=new_array1(N-1);
	
	//main circulation
	for(j=0;j<run_time;j++){
		for(i=1;i<N-1;i++){
			u[i]=A*Psi[i-1][j]+(1-2*A-B[i])*Psi[i][j]+A*Psi[i+1][j];
		}
		u[1]=u[1]-A*Psi[0][j+1];                              
		u[N-2]=u[N-2]-A*Psi[N-1][j+1];
		
		for(z=0;z<N-1;z++){
				dia0[z]=1+2*A+B[z];dia1[z]=-1*A;dian1[z]=-1*A;
		}
		
		//Thomas algorithm
		for(z=1;z<N-1;z++){
			dian1[z-1]=dian1[z-1]/dia0[z-1];
			dia0[z]=dia0[z]-dia1[z-1]*dian1[z-1];
			u[z]=u[z]-u[z-1]*dian1[z-1];
		}
		v[N-2]=u[N-2]/dia0[N-2];
		for(z=N-3;z>=0;z--){
			v[z]=(u[z]-dia1[z]*v[z+1])/dia0[z];
		}	
		
		for(z=1;z<N-1;z++){
			Psi[z][j+1]=v[z];
		}
	}
	//write
	FILE *fp=NULL;                                        
		if((fp=fopen("evolution.txt","w+"))==NULL){
			printf("Error on open\n");
		exit(1);
	}
	for(i=0;i<N;i++){
		for(j=0;j<run_time;j++){
			if(j%10==0){
				fprintf(fp,"%g\t",cabs(Psi[i][j])*cabs(Psi[i][j]));
			}
		}
		fprintf(fp,"\n");
	}
	fclose(fp);
	//the sum of probability in sampling time
	_Complex double *prob=new_array1(N);
	double initial_value;
	FILE *fp1=NULL;
	if((fp1=fopen("reducibility.txt","w+"))==NULL){
			printf("Error on open\n");
		exit(1);
	}
	fprintf(fp1,"Iteration\tSum\tDifference with the initial value");
	for(j=0;j<run_time;j++){
		debug();
		if(j%100==0){
			for(z=0;z<N;z++){
				prob[z]=Psi[z][j];
			}
			if(j==0)initial_value=sum(prob,N)*delta_h;
			fprintf(fp1,"%d\t%lf\t%lf\n",j,sum(prob,N)*delta_h,sum(prob,N)*delta_h-initial_value);
		}
	}
	fclose(fp1);
}
