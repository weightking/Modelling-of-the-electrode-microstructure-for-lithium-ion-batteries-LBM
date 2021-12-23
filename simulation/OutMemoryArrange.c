#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include "MemoryArrange.h"
#include "matrixMove.h"
#include "twoPhaseField.h"
#include "memoryFree.h"

double ***Out0=NULL;
double ***Out1=NULL;
double ***Out2=NULL;
double ***Out3=NULL;
double ***Out4=NULL;
double ***Out5=NULL;
double ***Out6=NULL;
double ***Out7=NULL;
double ***Out8=NULL;
double ***Out9=NULL;
double ***Out10=NULL;
double ***Out11=NULL;
double ***Out12=NULL;
double ***Out13=NULL;
double ***Out14=NULL;
double ***Out15=NULL;
double ***Out16=NULL;
double ***Out17=NULL;
double ***Out18=NULL;

void OutMemoryArrange(int m, int n, int q){
	Out0=memoryarrange(m,n,q);
	Out1=memoryarrange(m,n,q);
	Out2=memoryarrange(m,n,q);
	Out3=memoryarrange(m,n,q);
	Out4=memoryarrange(m,n,q);
	Out5=memoryarrange(m,n,q);
	Out6=memoryarrange(m,n,q);
//	Out7=memoryarrange(m,n,q);
//	Out8=memoryarrange(m,n,q);
//	Out9=memoryarrange(m,n,q);
//	Out10=memoryarrange(m,n,q);
//	Out11=memoryarrange(m,n,q);
//	Out12=memoryarrange(m,n,q);
//	Out13=memoryarrange(m,n,q);
//	Out14=memoryarrange(m,n,q);
//	Out15=memoryarrange(m,n,q);
//	Out16=memoryarrange(m,n,q);
//	Out17=memoryarrange(m,n,q);
//	Out18=memoryarrange(m,n,q);
}

void OutMemoryFree(int m, int n, int q){
	freememory(Out0,m,n,q);
	freememory(Out1,m,n,q);
	freememory(Out2,m,n,q);
	freememory(Out3,m,n,q);
	freememory(Out4,m,n,q);
	freememory(Out5,m,n,q);
	freememory(Out6,m,n,q);
//	freememory(Out7,m,n,q);
//	freememory(Out8,m,n,q);
//	freememory(Out9,m,n,q);
//	freememory(Out10,m,n,q);
//	freememory(Out11,m,n,q);
//	freememory(Out12,m,n,q);
//	freememory(Out13,m,n,q);
//	freememory(Out14,m,n,q);
//	freememory(Out15,m,n,q);
//	freememory(Out16,m,n,q);
//	freememory(Out17,m,n,q);
//	freememory(Out18,m,n,q);
}