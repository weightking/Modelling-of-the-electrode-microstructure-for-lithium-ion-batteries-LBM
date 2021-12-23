#include<stdio.h>
void shift0(double ***a, double ***b, int m, int n, int q){
	int i;
	int j;
	int k;
#pragma omp parallel private(i,j,k)
	{
#pragma omp for schedule(dynamic)
		for (i=0;i<m;i++){
			for (j=0;j<n;j++){
				for (k=0;k<q;k++){
					b[i][j][k]=a[i][j][k];
				}
			}
		}
	}
}

void shift1(double ***a, double ***b, int m, int n, int	q){
	int i;
	int j;
	int k;
#pragma omp parallel private(i,j,k)
	{
#pragma omp for schedule(dynamic)
		for (i=0;i<m-1;i++){
			for (j=0;j<n;j++){
				for (k=0;k<q;k++){
					b[i+1][j][k]=a[i][j][k];
				}
			}
		}
#pragma omp for schedule(dynamic)
		for (j=0;j<n;j++){
			for (k=0;k<q;k++){
				b[0][j][k]=a[m-1][j][k];
			}
		}
	}
}

double neightbour1(int m, int n, int q, int boundary, double ***a){
	double b;
	if (m==(boundary-1)){
		b=a[0][n][q];
	}
	else{
		b=a[m+1][n][q];
	}
	return b;
}

void shift2(double ***a, double ***b, int m, int n, int q){
	int i;
	int j;
	int k;
#pragma omp parallel private(i,j,k)
	{
#pragma omp for schedule(dynamic)
		for (i=0;i<m;i++){
			for (j=0;j<n-1;j++){
				for (k=0;k<q;k++){
					b[i][j+1][k]=a[i][j][k];
				}
			}
		}
#pragma omp for schedule(dynamic)
		for (i=0;i<m;i++){
			for (k=0;k<q;k++){
				b[i][0][k]=a[i][n-1][k];
			}
		}
	}
}

double neightbour2(int m, int n, int q, int boundary, double ***a){
	double b;
	if (n==(boundary-1)){
		b=a[m][0][q];
	}
	else{
		b=a[m][n+1][q];
	}
	return b;
}

void shift3(double ***a, double ***b, int m, int n, int q){
	int i;
	int j;
	int k;
#pragma omp parallel private(i,j,k)
	{
#pragma omp for schedule(dynamic)
		for (i=0;i<m;i++){
			for (j=0;j<n;j++){
				for (k=0;k<q-1;k++){
					b[i][j][k+1]=a[i][j][k];
				}
			}
		}
#pragma omp for schedule(dynamic)
		for (i=0;i<m;i++){
			for (j=0;j<n;j++){
				b[i][j][0]=a[i][j][q-1];
			}
		}
	}
}

double neightbour3(int m, int n, int q, int boundary, double ***a){
	double b;
	if (q==(boundary-1)){
		b=a[m][n][0];
	}
	else{
		b=a[m][n][q+1];
	}
	return b;
}

void shift4(double ***a, double ***b, int m, int n, int q){
	int i;
	int j;
	int k;
#pragma omp parallel private(i,j,k)
	{
#pragma omp for schedule(dynamic)
		for (i=1;i<m;i++){
			for (j=0;j<n;j++){
				for (k=0;k<q;k++){
					b[i-1][j][k]=a[i][j][k];
				}
			}
		}
#pragma omp for schedule(dynamic)
		for (j=0;j<n;j++){
			for (k=0;k<q;k++){
				b[m-1][j][k]=a[0][j][k];
			}
		}
	}
}

double neightbour4(int m, int n, int q, int boundary, double ***a){
	double b;
	if (m==0){
		b=a[boundary-1][n][q];
	}
	else{
		b=a[m-1][n][q];
	}
	return b;
}

void shift5(double ***a, double ***b, int m, int n, int q){
	int i;
	int j;
	int k;
#pragma omp parallel private(i,j,k)
	{
#pragma omp for schedule(dynamic)
		for (i=0;i<m;i++){
			for (j=1;j<n;j++){
				for (k=0;k<q;k++){
					b[i][j-1][k]=a[i][j][k];
				}
			}
		}
#pragma omp for schedule(dynamic)
		for (i=0;i<m;i++){
			for (k=0;k<q;k++){
				b[i][n-1][k]=a[i][0][k];
			}
		}
	}
}

double neightbour5(int m, int n, int q, int boundary, double ***a){
	double b;
	if (n==0){
		b=a[m][boundary-1][q];
	}
	else{
		b=a[m][n-1][q];
	}
	return b;
}

void shift6(double ***a, double ***b, int m, int n, int q){
	int i;
	int j;
	int k;
#pragma omp parallel private(i,j,k)
	{
#pragma omp for schedule(dynamic)
		for (i=0;i<m;i++){
			for (j=0;j<n;j++){
				for (k=1;k<q;k++){
					b[i][j][k-1]=a[i][j][k];
				}
			}
		}
#pragma omp for schedule(dynamic)
		for (i=0;i<m;i++){
			for (j=0;j<n;j++){
				b[i][j][q-1]=a[i][j][0];
			}
		}
	}
}

double neightbour6(int m, int n, int q, int boundary, double ***a){
	double b;
	if (q==0){
		b=a[m][n][boundary-1];
	}
	else{
		b=a[m][n][q-1];
	}
	return b;
}

void shift7(double ***a, double ***b, int m, int n, int q){
	int i;
	int j;
	int k;
#pragma omp parallel private(i,j,k)
	{
#pragma omp for schedule(dynamic)
		for (i=0;i<m-1;i++){
			for (j=0;j<n-1;j++){
				for (k=0;k<q;k++){
					b[i+1][j+1][k]=a[i][j][k];
				}
			}
		}
#pragma omp for schedule(dynamic)
		for (j=0;j<n-1;j++){
			for (k=0;k<q;k++){
				b[0][j+1][k]=a[m-1][j][k];
			}
		}
#pragma omp for schedule(dynamic)
		for (i=0;i<m-1;i++){
			for (k=0;k<q;k++){
				b[i+1][0][k]=a[i][n-1][k];
			}
		}
#pragma omp for schedule(dynamic)
		for (k=0;k<q;k++){
			b[0][0][k]=a[m-1][n-1][k];
		}
	}
}

double neightbour7(int m, int n, int q, int boundary1, int boundary2, double ***a){
	double b;
	if (m==(boundary1-1)&&n==(boundary2-1)){
		b=a[0][0][q];
	}
	if (m==(boundary1-1)&&n!=(boundary2-1)){
		b=a[0][n+1][q];
	}
	if (m!=(boundary1-1)&&n==(boundary2-1)){
		b=a[m+1][0][q];
	}
	if (m!=(boundary1-1)&&n!=(boundary2-1)){
		b=a[m+1][n+1][q];
	}
	return b;
}

void shift8(double ***a, double ***b, int m, int n, int q){
	int i;
	int j;
	int k;
#pragma omp parallel private(i,j,k)
	{
#pragma omp for schedule(dynamic)
		for (i=0;i<m-1;i++){
			for (j=0;j<n;j++){
				for (k=0;k<q-1;k++){
					b[i+1][j][k+1]=a[i][j][k];
				}
			}
		}
#pragma omp for schedule(dynamic)
		for (j=0;j<n;j++){
			for (k=0;k<q-1;k++){
				b[0][j][k+1]=a[m-1][j][k];
			}
		}
#pragma omp for schedule(dynamic)
		for (i=0;i<m-1;i++){
			for (j=0;j<n;j++){
				b[i+1][j][0]=a[i][j][q-1];
			}
		}
#pragma omp for schedule(dynamic)
		for (j=0;j<n;j++){
			b[0][j][0]=a[m-1][j][q-1];
		}
	}
}

double neightbour8(int m, int n, int q, int boundary1, int boundary2, double ***a){
	double b;
	if (m==(boundary1-1)&&q==(boundary2-1)){
		b=a[0][n][0];
	}
	if (m==(boundary1-1)&&q!=(boundary2-1)){
		b=a[0][n][q+1];
	}
	if (m!=(boundary1-1)&&q==(boundary2-1)){
		b=a[m+1][n][0];
	}
	if (m!=(boundary1-1)&&q!=(boundary2-1)){
		b=a[m+1][n][q+1];
	}
	return b;
}

void shift9(double ***a, double ***b, int m, int n, int q){
	int i;
	int j;
	int k;
#pragma omp parallel private(i,j,k)
	{
#pragma omp for schedule(dynamic)
		for (i=0;i<m;i++){
			for (j=0;j<n-1;j++){
				for (k=0;k<q-1;k++){
					b[i][j+1][k+1]=a[i][j][k];
				}
			}
		}
#pragma omp for schedule(dynamic)
		for (i=0;i<m;i++){
			for (j=0;j<n-1;j++){
				b[i][j+1][0]=a[i][j][q-1];
			}
		}
#pragma omp for schedule(dynamic)
		for (i=0;i<m;i++){
			for (k=0;k<q-1;k++){
				b[i][0][k+1]=a[i][n-1][k];
			}
		}
#pragma omp for schedule(dynamic)
		for (i=0;i<m;i++){
			b[i][0][0]=a[i][n-1][q-1];
		}
	}
}

double neightbour9(int m, int n, int q, int boundary1, int boundary2, double ***a){
	double b;
	if (n==(boundary1-1)&&q==(boundary2-1)){
		b=a[m][0][0];
	}
	if (n==(boundary1-1)&&q!=(boundary2-1)){
		b=a[m][0][q+1];
	}
	if (n!=(boundary1-1)&&q==(boundary2-1)){
		b=a[m][n+1][0];
	}
	if (n!=(boundary1-1)&&q!=(boundary2-1)){
		b=a[m][n+1][q+1];
	}
	return b;
}

void shift10(double ***a, double ***b, int m, int n, int q){
	int i;
	int j;
	int k;
#pragma omp parallel private(i,j,k)
	{
#pragma omp for schedule(dynamic)
		for (i=1;i<m;i++){
			for (j=1;j<n;j++){
				for (k=0;k<q;k++){
					b[i-1][j-1][k]=a[i][j][k];
				}
			}
		}
#pragma omp for schedule(dynamic)
		for (j=1;j<n;j++){
			for (k=0;k<q;k++){
				b[m-1][j-1][k]=a[0][j][k];
			}
		}
#pragma omp for schedule(dynamic)
		for (i=1;i<m;i++){
			for (k=0;k<q;k++){
				b[i-1][n-1][k]=a[i][0][k];
			}
		}
#pragma omp for schedule(dynamic)
		for (k=0;k<q;k++){
			b[m-1][n-1][k]=a[0][0][k];
		}
	}
}

double neightbour10(int m, int n, int q, int boundary1, int boundary2, double ***a){
	double b;
	if (m==0&&n==0){
		b=a[boundary1-1][boundary2-1][q];
	}
	if (m==0&&n!=0){
		b=a[boundary1-1][n-1][q];
	}
	if (m!=0&&n==0){
		b=a[m-1][boundary2-1][q];
	}
	if (m!=0&&n!=0){
		b=a[m-1][n-1][q];
	}
	return b;
}

void shift11(double ***a, double ***b, int m, int n, int q){
	int i;
	int j;
	int k;
#pragma omp parallel private(i,j,k)
	{
#pragma omp for schedule(dynamic)
		for (i=1;i<m;i++){
			for (j=0;j<n;j++){
				for (k=1;k<q;k++){
					b[i-1][j][k-1]=a[i][j][k];
				}
			}
		}
#pragma omp for schedule(dynamic)
		for (j=0;j<n;j++){
			for (k=1;k<q;k++){
				b[m-1][j][k-1]=a[0][j][k];
			}
		}
#pragma omp for schedule(dynamic)
		for (i=1;i<m;i++){
			for (j=0;j<n;j++){
				b[i-1][j][q-1]=a[i][j][0];
			}
		}
#pragma omp for schedule(dynamic)
		for (j=0;j<n;j++){
			b[m-1][j][q-1]=a[0][j][0];
		}
	}
}

double neightbour11(int m, int n, int q, int boundary1, int boundary2, double ***a){
	double b;
	if (m==0&&q==0){
		b=a[boundary1-1][n][boundary2-1];
	}
	if (m==0&&q!=0){
		b=a[boundary1-1][n][q-1];
	}
	if (m!=0&&q==0){
		b=a[m-1][n][boundary2-1];
	}
	if (m!=0&&q!=0){
		b=a[m-1][n][q-1];
	}
	return b;
}

void shift12(double ***a, double ***b, int m, int n, int q){
	int i;
	int j;
	int k;
#pragma omp parallel private(i,j,k)
	{
#pragma omp for schedule(dynamic)
		for (i=0;i<m;i++){
			for (j=1;j<n;j++){
				for (k=1;k<q;k++){
					b[i][j-1][k-1]=a[i][j][k];
				}
			}
		}
#pragma omp for schedule(dynamic)
		for (i=0;i<m;i++){
			for (j=1;j<n;j++){
				b[i][j-1][q-1]=a[i][j][0];
			}
		}
#pragma omp for schedule(dynamic)
		for (i=0;i<m;i++){
			for (k=1;k<q;k++){
				b[i][n-1][k-1]=a[i][0][k];
			}
		}
#pragma omp for schedule(dynamic)
		for (i=0;i<m;i++){
			b[i][n-1][q-1]=a[i][0][0];
		}
	}
}

double neightbour12(int m, int n, int q, int boundary1, int boundary2, double ***a){
	double b;
	if (n==0&&q==0){
		b=a[m][boundary1-1][boundary2-1];
	}
	if (n==0&&q!=0){
		b=a[m][boundary1-1][q-1];
	}
	if (n!=0&&q==0){
		b=a[m][n-1][boundary2-1];
	}
	if (n!=0&&q!=0){
		b=a[m][n-1][q-1];
	}
	return b;
}

void shift13(double ***a, double ***b, int m, int n, int q){
	int i;
	int j;
	int k;
#pragma omp parallel private(i,j,k)
	{
#pragma omp for schedule(dynamic)
		for (i=0;i<m-1;i++){
			for (j=1;j<n;j++){
				for (k=0;k<q;k++){
					b[i+1][j-1][k]=a[i][j][k];
				}
			}
		}
#pragma omp for schedule(dynamic)
		for (j=1;j<n;j++){
			for (k=0;k<q;k++){
				b[0][j-1][k]=a[m-1][j][k];
			}
		}
#pragma omp for schedule(dynamic)
		for (i=0;i<m-1;i++){
			for (k=0;k<q;k++){
				b[i+1][n-1][k]=a[i][0][k];
			}
		}
#pragma omp for schedule(dynamic)
		for (k=0;k<q;k++){
			b[0][n-1][k]=a[m-1][0][k];
		}
	}
}

double neightbour13(int m, int n, int q, int boundary1, int boundary2, double ***a){
	double b;
	if (m==(boundary1-1)&&n==0){
		b=a[0][boundary2-1][q];
	}
	if (m==(boundary1-1)&&n!=0){
		b=a[0][n-1][q];
	}
	if (m!=(boundary1-1)&&n==0){
		b=a[m+1][boundary2-1][q];
	}
	if (m!=(boundary1-1)&&n!=0){
		b=a[m+1][n-1][q];
	}
	return b;
}

void shift14(double ***a, double ***b, int m, int n, int q){
	int i;
	int j;
	int k;
#pragma omp parallel private(i,j,k)
	{
#pragma omp for schedule(dynamic)
		for (i=0;i<m-1;i++){
			for (j=0;j<n;j++){
				for (k=1;k<q;k++){
					b[i+1][j][k-1]=a[i][j][k];
				}
			}
		}
#pragma omp for schedule(dynamic)
		for (j=0;j<n;j++){
			for (k=1;k<q;k++){
				b[0][j][k-1]=a[m-1][j][k];
			}
		}
#pragma omp for schedule(dynamic)
		for (i=0;i<m-1;i++){
			for (j=0;j<n;j++){
				b[i+1][j][q-1]=a[i][j][0];
			}
		}
#pragma omp for schedule(dynamic)
		for (j=0;j<n;j++){
			b[0][j][q-1]=a[m-1][j][0];
		}
	}
}

double neightbour14(int m, int n, int q, int boundary1, int boundary2, double ***a){
	double b;
	if (m==(boundary1-1)&&q==0){
		b=a[0][n][boundary2-1];
	}
	if (m==(boundary1-1)&&q!=0){
		b=a[0][n][q-1];
	}
	if (m!=(boundary1-1)&&q==0){
		b=a[m+1][n][boundary2-1];
	}
	if (m!=(boundary1-1)&&q!=0){
		b=a[m+1][n][q-1];
	}
	return b;
}

void shift15(double ***a, double ***b, int m, int n, int q){
	int i;
	int j;
	int k;
#pragma omp parallel private(i,j,k)
	{
#pragma omp for schedule(dynamic)
		for (i=0;i<m;i++){
			for (j=0;j<n-1;j++){
				for (k=1;k<q;k++){
					b[i][j+1][k-1]=a[i][j][k];
				}
			}
		}
#pragma omp for schedule(dynamic)
		for (i=0;i<m;i++){
			for (j=0;j<n-1;j++){
				b[i][j+1][q-1]=a[i][j][0];
			}
		}
#pragma omp for schedule(dynamic)
		for (i=0;i<m;i++){
			for (k=1;k<q;k++){
				b[i][0][k-1]=a[i][n-1][k];
			}
		}
#pragma omp for schedule(dynamic)
		for (i=0;i<m;i++){
			b[i][0][q-1]=a[i][n-1][0];
		}
	}
}

double neightbour15(int m, int n, int q, int boundary1, int boundary2, double ***a){
	double b;
	if (n==(boundary1-1)&&q==0){
		b=a[m][0][boundary2-1];
	}
	if (n==(boundary1-1)&&q!=0){
		b=a[m][0][q-1];
	}
	if (n!=(boundary1-1)&&q==0){
		b=a[m][n+1][boundary2-1];
	}
	if (n!=(boundary1-1)&&q!=0){
		b=a[m][n+1][q-1];
	}
	return b;
}

void shift16(double ***a, double ***b, int m, int n, int q){
	int i;
	int j;
	int k;
#pragma omp parallel private(i,j,k)
	{
#pragma omp for schedule(dynamic)
		for (i=1;i<m;i++){
			for (j=0;j<n-1;j++){
				for (k=0;k<q;k++){
					b[i-1][j+1][k]=a[i][j][k];
				}
			}
		}
#pragma omp for schedule(dynamic)
		for (j=0;j<n-1;j++){
			for (k=0;k<q;k++){
				b[m-1][j+1][k]=a[0][j][k];
			}
		}
#pragma omp for schedule(dynamic)
		for (i=1;i<m;i++){
			for (k=0;k<q;k++){
				b[i-1][0][k]=a[i][n-1][k];
			}
		}
#pragma omp for schedule(dynamic)
		for (k=0;k<q;k++){
			b[m-1][0][k]=a[0][n-1][k];
		}
	}
}

double neightbour16(int m, int n, int q, int boundary1, int boundary2, double ***a){
	double b;
	if (m==0&&n==(boundary2-1)){
		b=a[boundary1-1][0][q];
	}
	if (m==0&&n!=(boundary2-1)){
		b=a[boundary1-1][n+1][q];
	}
	if (m!=0&&n==(boundary2-1)){
		b=a[m-1][0][q];
	}
	if (m!=0&&n!=(boundary2-1)){
		b=a[m-1][n+1][q];
	}
	return b;
}

void shift17(double ***a, double ***b, int m, int n, int q){
	int i;
	int j;
	int k;
#pragma omp parallel private(i,j,k)
	{
#pragma omp for schedule(dynamic)
		for (i=1;i<m;i++){
			for (j=0;j<n;j++){
				for (k=0;k<q-1;k++){
					b[i-1][j][k+1]=a[i][j][k];
				}
			}
		}
#pragma omp for schedule(dynamic)
		for (j=0;j<n;j++){
			for (k=0;k<q-1;k++){
				b[m-1][j][k+1]=a[0][j][k];
			}
		}
#pragma omp for schedule(dynamic)
		for (i=1;i<m;i++){
			for (j=0;j<n;j++){
				b[i-1][j][0]=a[i][j][q-1];
			}
		}
#pragma omp for schedule(dynamic)
		for (j=0;j<n;j++){
			b[m-1][j][0]=a[0][j][q-1];
		}
	}
}

double neightbour17(int m, int n, int q, int boundary1, int boundary2, double ***a){
	double b;
	if (m==0&&q==(boundary2-1)){
		b=a[boundary1-1][n][0];
	}
	if (m==0&&q!=(boundary2-1)){
		b=a[boundary1-1][n][q+1];
	}
	if (m!=0&&q==(boundary2-1)){
		b=a[m-1][n][0];
	}
	if (m!=0&&q!=(boundary2-1)){
		b=a[m-1][n][q+1];
	}
	return b;
}

void shift18(double ***a, double ***b, int m, int n, int q){
	int i;
	int j;
	int k;
#pragma omp parallel private(i,j,k)
	{
#pragma omp for schedule(dynamic)
		for (i=0;i<m;i++){
			for (j=1;j<n;j++){
				for (k=0;k<q-1;k++){
					b[i][j-1][k+1]=a[i][j][k];
				}
			}
		}
#pragma omp for schedule(dynamic)
		for (i=0;i<m;i++){
			for (k=0;k<q-1;k++){
				b[i][n-1][k+1]=a[i][0][k];
			}
		}
#pragma omp for schedule(dynamic)
		for (i=0;i<m;i++){
			for (j=1;j<n;j++){
				b[i][j-1][0]=a[i][j][q-1];
			}
		}
#pragma omp for schedule(dynamic)
		for (i=0;i<m;i++){
			b[i][n-1][0]=a[i][0][q-1];
		}
	}
}

double neightbour18(int m, int n, int q, int boundary1, int boundary2, double ***a){
	double b;
	if (n==0&&q==(boundary2-1)){
		b=a[m][boundary1-1][0];
	}
	if (n==0&&q!=(boundary2-1)){
		b=a[m][boundary1-1][q+1];
	}
	if (n!=0&&q==(boundary2-1)){
		b=a[m][n-1][0];
	}
	if (n!=0&&q!=(boundary2-1)){
		b=a[m][n-1][q+1];
	}
	return b;
}

void shift19(double ***a, double ***b, int m, int n, int q){
	int i;
	int j;
	int k;
#pragma omp parallel private(i,j,k)
	{
#pragma omp for schedule(dynamic)
		for (i = 0; i<m - 1; i++){
			for (j = 0; j<n - 1; j++){
				for (k = 0; k<q-1; k++){
					b[i + 1][j + 1][k+1] = a[i][j][k];
				}
			}
		}
#pragma omp for schedule(dynamic)
		for (j = 0; j<n - 1; j++){
			for (k = 0; k<q-1; k++){
				b[0][j + 1][k+1] = a[m - 1][j][k];
			}
		}
#pragma omp for schedule(dynamic)
		for (i = 0; i<m - 1; i++){
			for (k = 0; k<q-1; k++){
				b[i + 1][0][k+1] = a[i][n - 1][k];
			}
		}
#pragma omp for schedule(dynamic)
		for (i = 0; i<m - 1; i++){
			for (j = 0; j<n - 1; j++){
				b[i + 1][j+1][0] = a[i][j][q-1];
			}
		}
#pragma omp for schedule(dynamic)
		for (i = 0; i<m - 1; i++){
				b[i + 1][0][0] = a[i][n-1][q - 1];
		}
#pragma omp for schedule(dynamic)
		for (j = 0; j<n - 1; j++){
				b[0][j + 1][0] = a[m-1][j][q - 1];
		}
#pragma omp for schedule(dynamic)
		for (k = 0; k<q - 1; k++){
				b[0][0][k + 1] = a[m-1][n - 1][k];
		}
			b[0][0][0] = a[m - 1][n - 1][q-1];
	}
}

double neightbour19(int m, int n, int q, int boundary1, int boundary2, int boundary3, double ***a){
	double b;
	if (m == (boundary1 - 1) && n == (boundary2 - 1) && q==(boundary3-1)){
		b = a[0][0][0];
	}
	if (m == (boundary1 - 1) && n != (boundary2 - 1) && q!=(boundary3-1)){
		b = a[0][n + 1][q+1];
	}
	if (m != (boundary1 - 1) && n == (boundary2 - 1) && q != (boundary3 - 1)){
		b = a[m+1][0][q + 1];
	}
	if (m != (boundary1 - 1) && n != (boundary2 - 1) && q == (boundary3 - 1)){
		b = a[m+1][n + 1][0];
	}
	if (m != (boundary1 - 1) && n == (boundary2 - 1) && q == (boundary3 - 1)){
		b = a[m+1][0][0];
	}
	if (m == (boundary1 - 1) && n != (boundary2 - 1) && q == (boundary3 - 1)){
		b = a[0][n+1][0];
	}
	if (m == (boundary1 - 1) && n == (boundary2 - 1) && q != (boundary3 - 1)){
		b = a[0][0][q+1];
	}
	if (m != (boundary1 - 1) && n != (boundary2 - 1) && q != (boundary3 - 1)){
		b = a[m+1][n+1][q+1];
	}
	return b;
}

void shift20(double ***a, double ***b, int m, int n, int q){
	int i;
	int j;
	int k;
#pragma omp parallel private(i,j,k)
	{
#pragma omp for schedule(dynamic)
		for (i = 0; i<m - 1; i++){
			for (j = 0; j<n - 1; j++){
				for (k = 1; k<q; k++){
					b[i + 1][j + 1][k - 1] = a[i][j][k];
				}
			}
		}
#pragma omp for schedule(dynamic)
		for (j = 0; j<n - 1; j++){
			for (k = 1; k<q; k++){
				b[0][j + 1][k - 1] = a[m - 1][j][k];
			}
		}
#pragma omp for schedule(dynamic)
		for (i = 0; i<m - 1; i++){
			for (k = 1; k<q; k++){
				b[i + 1][0][k - 1] = a[i][n - 1][k];
			}
		}
#pragma omp for schedule(dynamic)
		for (i = 0; i<m - 1; i++){
			for (j = 0; j<n - 1; j++){
				b[i + 1][j + 1][q-1] = a[i][j][0];
			}
		}
#pragma omp for schedule(dynamic)
		for (i = 0; i<m - 1; i++){
			b[i + 1][0][q-1] = a[i][n - 1][0];
		}
#pragma omp for schedule(dynamic)
		for (j = 0; j<n - 1; j++){
			b[0][j + 1][q-1] = a[m - 1][j][0];
		}
#pragma omp for schedule(dynamic)
		for (k = 1; k<q; k++){
			b[0][0][k - 1] = a[m - 1][n - 1][k];
		}
		b[0][0][q-1] = a[m - 1][n - 1][0];
	}
}

double neightbour20(int m, int n, int q, int boundary1, int boundary2, int boundary3, double ***a){
	double b;
	if (m == (boundary1 - 1) && n == (boundary2 - 1) && q == 0){
		b = a[0][0][boundary3-1];
	}
	if (m == (boundary1 - 1) && n != (boundary2 - 1) && q != 0){
		b = a[0][n + 1][q - 1];
	}
	if (m != (boundary1 - 1) && n == (boundary2 - 1) && q != 0){
		b = a[m + 1][0][q - 1];
	}
	if (m != (boundary1 - 1) && n != (boundary2 - 1) && q == 0){
		b = a[m + 1][n + 1][boundary3-1];
	}
	if (m != (boundary1 - 1) && n == (boundary2 - 1) && q == 0){
		b = a[m + 1][0][boundary3-1];
	}
	if (m == (boundary1 - 1) && n != (boundary2 - 1) && q == 0){
		b = a[0][n + 1][boundary3-1];
	}
	if (m == (boundary1 - 1) && n == (boundary2 - 1) && q != 0){
		b = a[0][0][q - 1];
	}
	if (m != (boundary1 - 1) && n != (boundary2 - 1) && q != 0){
		b = a[m + 1][n + 1][q - 1];
	}
	return b;
}

void shift21(double ***a, double ***b, int m, int n, int q){
	int i;
	int j;
	int k;
#pragma omp parallel private(i,j,k)
	{
#pragma omp for schedule(dynamic)
		for (i = 0; i<m - 1; i++){
			for (j = 1; j<n; j++){
				for (k = 0; k<q - 1; k++){
					b[i + 1][j - 1][k + 1] = a[i][j][k];
				}
			}
		}
#pragma omp for schedule(dynamic)
		for (j = 1; j<n; j++){
			for (k = 0; k<q - 1; k++){
				b[0][j - 1][k + 1] = a[m - 1][j][k];
			}
		}
#pragma omp for schedule(dynamic)
		for (i = 0; i<m - 1; i++){
			for (k = 0; k<q - 1; k++){
				b[i + 1][n-1][k + 1] = a[i][0][k];
			}
		}
#pragma omp for schedule(dynamic)
		for (i = 0; i<m - 1; i++){
			for (j = 1; j<n; j++){
				b[i + 1][j - 1][0] = a[i][j][q - 1];
			}
		}
#pragma omp for schedule(dynamic)
		for (i = 0; i<m - 1; i++){
			b[i + 1][n-1][0] = a[i][0][q - 1];
		}
#pragma omp for schedule(dynamic)
		for (j = 1; j<n; j++){
			b[0][j - 1][0] = a[m - 1][j][q - 1];
		}
#pragma omp for schedule(dynamic)
		for (k = 0; k<q - 1; k++){
			b[0][n-1][k + 1] = a[m - 1][0][k];
		}
		b[0][n-1][0] = a[m - 1][0][q - 1];
	}
}

double neightbour21(int m, int n, int q, int boundary1, int boundary2, int boundary3, double ***a){
	double b;
	if (m == (boundary1 - 1) && n == 0 && q == (boundary3 - 1)){
		b = a[0][boundary2-1][0];
	}
	if (m == (boundary1 - 1) && n != 0 && q != (boundary3 - 1)){
		b = a[0][n - 1][q + 1];
	}
	if (m != (boundary1 - 1) && n == 0 && q != (boundary3 - 1)){
		b = a[m + 1][boundary2-1][q + 1];
	}
	if (m != (boundary1 - 1) && n != 0 && q == (boundary3 - 1)){
		b = a[m + 1][n - 1][0];
	}
	if (m != (boundary1 - 1) && n == 0 && q == (boundary3 - 1)){
		b = a[m + 1][boundary2-1][0];
	}
	if (m == (boundary1 - 1) && n != 0 && q == (boundary3 - 1)){
		b = a[0][n - 1][0];
	}
	if (m == (boundary1 - 1) && n == 0 && q != (boundary3 - 1)){
		b = a[0][boundary2-1][q + 1];
	}
	if (m != (boundary1 - 1) && n != 0 && q != (boundary3 - 1)){
		b = a[m + 1][n - 1][q + 1];
	}
	return b;
}

void shift22(double ***a, double ***b, int m, int n, int q){
	int i;
	int j;
	int k;
#pragma omp parallel private(i,j,k)
	{
#pragma omp for schedule(dynamic)
		for (i = 0; i<m - 1; i++){
			for (j = 1; j<n; j++){
				for (k = 1; k<q; k++){
					b[i + 1][j - 1][k - 1] = a[i][j][k];
				}
			}
		}
#pragma omp for schedule(dynamic)
		for (j = 1; j<n; j++){
			for (k = 1; k<q; k++){
				b[0][j - 1][k - 1] = a[m - 1][j][k];
			}
		}
#pragma omp for schedule(dynamic)
		for (i = 0; i<m - 1; i++){
			for (k = 1; k<q; k++){
				b[i + 1][n-1][k - 1] = a[i][0][k];
			}
		}
#pragma omp for schedule(dynamic)
		for (i = 0; i<m - 1; i++){
			for (j = 1; j<n; j++){
				b[i + 1][j - 1][q-1] = a[i][j][0];
			}
		}
#pragma omp for schedule(dynamic)
		for (i = 0; i<m - 1; i++){
			b[i + 1][n-1][q-1] = a[i][0][0];
		}
#pragma omp for schedule(dynamic)
		for (j = 1; j<n; j++){
			b[0][j - 1][q-1] = a[m - 1][j][0];
		}
#pragma omp for schedule(dynamic)
		for (k = 1; k<q; k++){
			b[0][n-1][k - 1] = a[m - 1][0][k];
		}
		b[0][n-1][q-1] = a[m - 1][0][0];
	}
}

double neightbour22(int m, int n, int q, int boundary1, int boundary2, int boundary3, double ***a){
	double b;
	if (m == (boundary1 - 1) && n == (0) && q == (0)){
		b = a[0][boundary2-1][boundary3-1];
	}
	if (m == (boundary1 - 1) && n != 0 && q != 0){
		b = a[0][n - 1][q - 1];
	}
	if (m != (boundary1 - 1) && n == 0 && q != 0){
		b = a[m + 1][boundary2-1][q - 1];
	}
	if (m != (boundary1 - 1) && n != 0 && q == 0){
		b = a[m + 1][n - 1][boundary3-1];
	}
	if (m != (boundary1 - 1) && n == 0 && q == 0){
		b = a[m + 1][boundary2-1][boundary3-1];
	}
	if (m == (boundary1 - 1) && n != 0 && q == 0){
		b = a[0][n - 1][boundary3-1];
	}
	if (m == (boundary1 - 1) && n == 0 && q != 0){
		b = a[0][boundary2-1][q - 1];
	}
	if (m != (boundary1 - 1) && n != 0 && q != 0){
		b = a[m + 1][n - 1][q - 1];
	}
	return b;
}

void shift23(double ***a, double ***b, int m, int n, int q){
	int i;
	int j;
	int k;
#pragma omp parallel private(i,j,k)
	{
#pragma omp for schedule(dynamic)
		for (i = 1; i<m; i++){
			for (j = 0; j<n - 1; j++){
				for (k = 0; k<q - 1; k++){
					b[i - 1][j + 1][k + 1] = a[i][j][k];
				}
			}
		}
#pragma omp for schedule(dynamic)
		for (j = 0; j<n - 1; j++){
			for (k = 0; k<q - 1; k++){
				b[m-1][j + 1][k + 1] = a[0][j][k];
			}
		}
#pragma omp for schedule(dynamic)
		for (i = 1; i<m; i++){
			for (k = 0; k<q - 1; k++){
				b[i - 1][0][k + 1] = a[i][n - 1][k];
			}
		}
#pragma omp for schedule(dynamic)
		for (i = 1; i<m; i++){
			for (j = 0; j<n - 1; j++){
				b[i - 1][j + 1][0] = a[i][j][q - 1];
			}
		}
#pragma omp for schedule(dynamic)
		for (i = 1; i<m; i++){
			b[i - 1][0][0] = a[i][n - 1][q - 1];
		}
#pragma omp for schedule(dynamic)
		for (j = 0; j<n - 1; j++){
			b[m-1][j + 1][0] = a[0][j][q - 1];
		}
#pragma omp for schedule(dynamic)
		for (k = 0; k<q - 1; k++){
			b[m-1][0][k + 1] = a[0][n - 1][k];
		}
		b[m-1][0][0] = a[0][n - 1][q - 1];
	}
}

double neightbour23(int m, int n, int q, int boundary1, int boundary2, int boundary3, double ***a){
	double b;
	if (m == 0 && n == (boundary2 - 1) && q == (boundary3 - 1)){
		b = a[boundary1-1][0][0];
	}
	if (m == 0 && n != (boundary2 - 1) && q != (boundary3 - 1)){
		b = a[boundary1-1][n + 1][q + 1];
	}
	if (m != 0 && n == (boundary2 - 1) && q != (boundary3 - 1)){
		b = a[m - 1][0][q + 1];
	}
	if (m != 0 && n != (boundary2 - 1) && q == (boundary3 - 1)){
		b = a[m - 1][n + 1][0];
	}
	if (m != 0 && n == (boundary2 - 1) && q == (boundary3 - 1)){
		b = a[m - 1][0][0];
	}
	if (m == 0 && n != (boundary2 - 1) && q == (boundary3 - 1)){
		b = a[boundary1-1][n + 1][0];
	}
	if (m == 0 && n == (boundary2 - 1) && q != (boundary3 - 1)){
		b = a[boundary1-1][0][q + 1];
	}
	if (m != 0 && n != (boundary2 - 1) && q != (boundary3 - 1)){
		b = a[m - 1][n + 1][q + 1];
	}
	return b;
}

void shift24(double ***a, double ***b, int m, int n, int q){
	int i;
	int j;
	int k;
#pragma omp parallel private(i,j,k)
	{
#pragma omp for schedule(dynamic)
		for (i = 1; i<m; i++){
			for (j = 0; j<n - 1; j++){
				for (k = 1; k<q; k++){
					b[i - 1][j + 1][k - 1] = a[i][j][k];
				}
			}
		}
#pragma omp for schedule(dynamic)
		for (j = 0; j<n - 1; j++){
			for (k = 1; k<q; k++){
				b[m-1][j + 1][k - 1] = a[0][j][k];
			}
		}
#pragma omp for schedule(dynamic)
		for (i = 1; i<m; i++){
			for (k = 1; k<q; k++){
				b[i - 1][0][k - 1] = a[i][n - 1][k];
			}
		}
#pragma omp for schedule(dynamic)
		for (i = 1; i<m; i++){
			for (j = 0; j<n - 1; j++){
				b[i - 1][j + 1][q-1] = a[i][j][0];
			}
		}
#pragma omp for schedule(dynamic)
		for (i = 1; i<m; i++){
			b[i - 1][0][q-1] = a[i][n - 1][0];
		}
#pragma omp for schedule(dynamic)
		for (j = 0; j<n - 1; j++){
			b[m-1][j + 1][q-1] = a[0][j][0];
		}
#pragma omp for schedule(dynamic)
		for (k = 1; k<q; k++){
			b[m-1][0][k - 1] = a[0][n - 1][k];
		}
		b[m-1][0][q-1] = a[0][n - 1][0];
	}
}

double neightbour24(int m, int n, int q, int boundary1, int boundary2, int boundary3, double ***a){
	double b;
	if (m == 0 && n == (boundary2 - 1) && q == 0){
		b = a[boundary1-1][0][boundary3-1];
	}
	if (m == 0 && n != (boundary2 - 1) && q != 0){
		b = a[boundary1-1][n + 1][q - 1];
	}
	if (m != 0 && n == (boundary2 - 1) && q != 0){
		b = a[m - 1][0][q - 1];
	}
	if (m != 0 && n != (boundary2 - 1) && q == 0){
		b = a[m - 1][n + 1][boundary3-1];
	}
	if (m != 0 && n == (boundary2 - 1) && q == 0){
		b = a[m - 1][0][boundary3-1];
	}
	if (m == 0 && n != (boundary2 - 1) && q == 0){
		b = a[boundary1-1][n + 1][boundary3-1];
	}
	if (m == 0 && n == (boundary2 - 1) && q != 0){
		b = a[boundary1-1][0][q - 1];
	}
	if (m != 0 && n != (boundary2 - 1) && q != 0){
		b = a[m - 1][n + 1][q - 1];
	}
	return b;
}

void shift25(double ***a, double ***b, int m, int n, int q){
	int i;
	int j;
	int k;
#pragma omp parallel private(i,j,k)
	{
#pragma omp for schedule(dynamic)
		for (i = 1; i<m; i++){
			for (j = 1; j<n; j++){
				for (k = 0; k<q - 1; k++){
					b[i - 1][j - 1][k + 1] = a[i][j][k];
				}
			}
		}
#pragma omp for schedule(dynamic)
		for (j = 1; j<n; j++){
			for (k = 0; k<q - 1; k++){
				b[m-1][j - 1][k + 1] = a[0][j][k];
			}
		}
#pragma omp for schedule(dynamic)
		for (i = 1; i<m; i++){
			for (k = 0; k<q - 1; k++){
				b[i - 1][n-1][k + 1] = a[i][0][k];
			}
		}
#pragma omp for schedule(dynamic)
		for (i = 1; i<m; i++){
			for (j = 1; j<n; j++){
				b[i - 1][j - 1][0] = a[i][j][q - 1];
			}
		}
#pragma omp for schedule(dynamic)
		for (i = 1; i<m; i++){
			b[i - 1][n-1][0] = a[i][0][q - 1];
		}
#pragma omp for schedule(dynamic)
		for (j = 1; j<n; j++){
			b[m-1][j - 1][0] = a[0][j][q - 1];
		}
#pragma omp for schedule(dynamic)
		for (k = 0; k<q - 1; k++){
			b[m-1][n-1][k + 1] = a[0][0][k];
		}
		b[m-1][n-1][0] = a[0][0][q - 1];
	}
}

double neightbour25(int m, int n, int q, int boundary1, int boundary2, int boundary3, double ***a){
	double b;
	if (m == 0 && n == 0 && q == (boundary3 - 1)){
		b = a[boundary1-1][boundary2-1][0];
	}
	if (m == 0 && n != 0 && q != (boundary3 - 1)){
		b = a[boundary1-1][n - 1][q + 1];
	}
	if (m != 0 && n == 0 && q != (boundary3 - 1)){
		b = a[m - 1][boundary2-1][q + 1];
	}
	if (m != 0 && n != 0 && q == (boundary3 - 1)){
		b = a[m - 1][n - 1][0];
	}
	if (m != 0 && n == 0 && q == (boundary3 - 1)){
		b = a[m - 1][boundary2-1][0];
	}
	if (m == 0 && n != 0 && q == (boundary3 - 1)){
		b = a[boundary1-1][n - 1][0];
	}
	if (m == 0 && n == 0 && q != (boundary3 - 1)){
		b = a[boundary1-1][boundary2-1][q + 1];
	}
	if (m != 0 && n != 0 && q != (boundary3 - 1)){
		b = a[m - 1][n - 1][q + 1];
	}
	return b;
}

void shift26(double ***a, double ***b, int m, int n, int q){
	int i;
	int j;
	int k;
#pragma omp parallel private(i,j,k)
	{
#pragma omp for schedule(dynamic)
		for (i = 1; i<m; i++){
			for (j = 1; j<n; j++){
				for (k = 1; k<q; k++){
					b[i - 1][j - 1][k - 1] = a[i][j][k];
				}
			}
		}
#pragma omp for schedule(dynamic)
		for (j = 1; j<n; j++){
			for (k = 1; k<q; k++){
				b[m - 1][j - 1][k - 1] = a[0][j][k];
			}
		}
#pragma omp for schedule(dynamic)
		for (i = 1; i<m; i++){
			for (k = 1; k<q; k++){
				b[i - 1][n - 1][k - 1] = a[i][0][k];
			}
		}
#pragma omp for schedule(dynamic)
		for (i = 1; i<m; i++){
			for (j = 1; j<n; j++){
				b[i - 1][j - 1][q-1] = a[i][j][0];
			}
		}
#pragma omp for schedule(dynamic)
		for (i = 1; i<m; i++){
			b[i - 1][n - 1][q-1] = a[i][0][0];
		}
#pragma omp for schedule(dynamic)
		for (j = 1; j<n; j++){
			b[m - 1][j - 1][q-1] = a[0][j][0];
		}
#pragma omp for schedule(dynamic)
		for (k = 1; k<q; k++){
			b[m - 1][n - 1][k - 1] = a[0][0][0];
		}
		b[m - 1][n - 1][q-1] = a[0][0][0];
	}
}

double neightbour26(int m, int n, int q, int boundary1, int boundary2, int boundary3, double ***a){
	double b;
	if (m == 0 && n == 0 && q == 0){
		b = a[boundary1 - 1][boundary2 - 1][boundary3-1];
	}
	if (m == 0 && n != 0 && q != 0){
		b = a[boundary1 - 1][n - 1][q - 1];
	}
	if (m != 0 && n == 0 && q != 0){
		b = a[m - 1][boundary2 - 1][q - 1];
	}
	if (m != 0 && n != 0 && q == 0){
		b = a[m - 1][n - 1][boundary3-1];
	}
	if (m != 0 && n == 0 && q == 0){
		b = a[m - 1][boundary2 - 1][boundary3-1];
	}
	if (m == 0 && n != 0 && q == 0){
		b = a[boundary1 - 1][n - 1][boundary3-1];
	}
	if (m == 0 && n == 0 && q != 0){
		b = a[boundary1 - 1][boundary2 - 1][q - 1];
	}
	if (m != 0 && n != 0 && q != 0){
		b = a[m - 1][n - 1][q - 1];
	}
	return b;
}