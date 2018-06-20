#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <limits>
#include <fstream>

using namespace std;

/********************************************************************/
/*****THIS IS COMPLETELY DEVOTED TO THE RANDOM NUMEBR GENERATOR******/
#define FNORM   (2.3283064365e-10)
#define RANDOM  ((ira[ip++] = ira[ip1++] + ira[ip2++]) ^ ira[ip3++])
#define FRANDOM (FNORM * RANDOM)
#define pm1 ((FRANDOM > 0.5) ? 1 : -1)

unsigned myrand, ira[256];
unsigned char ip, ip1, ip2, ip3;

unsigned rand4init(void)
{
  unsigned long long y;

  y = (myrand*16807LL);
  myrand = (y&0x7fffffff) + (y>>31);
  if (myrand&0x80000000)
    myrand = (myrand&0x7fffffff) + 1;
  return myrand;
}                                                                                                                        

void Init_Random(void)
{
  int i;

  ip = 128;
  ip1 = ip - 24;
  ip2 = ip - 55;
  ip3 = ip - 61;

  for (i=ip3; i<ip; i++)
    ira[i] = rand4init();
}

float gauss_ran(void)
{
  static int iset=0;
  static float gset;
  float fac, rsq, v1, v2;

  if (iset == 0) {
    do {
      v1 = 2.0 * FRANDOM - 1.0;
      v2 = 2.0 * FRANDOM - 1.0;
      rsq = v1 * v1 + v2 * v2;
    } while (rsq >= 1.0 || rsq == 0.0);
    fac = sqrt(-2.0 * log(rsq) / rsq);
    gset = v1 * fac;
    iset = 1;
    return v2 * fac;
  } else {
    iset = 0;
    return gset;
  }
}

/***********HERE WE FINISH THE RANDOM NUMBER GENERATOR********/
/*************************************************************/
//#define N 8

//VARIABLES
unsigned myrand_back; //seed
int dim, N, vol, Tin, steps;
int n;
struct table{
        int n;			//0,1,2 (spin or number of electrons in eg orbitals)
        int vec[6];			//sentido horario, empezando por las 12 hrs
        int frontera;		//1 si es frontera, 0 si no
		int x,y,z;

};

//FUNCIONES
//INICIALIZACION
void initsys(struct table *);


void deb(struct table *, int *);

//FICHEROS
ofstream info("info.dat");


void initsys(struct table *site){		//se generan los vecinos en sentido horario
	int i,j;
	if(dim==1){
		for(i=0; i<N; i++){
			if(i==0){		//vecino de izq
				site[i].vec[0]=-1;
				site[i].frontera=1;
			} else{
				site[i].vec[0]=i-1;
			}

			if(i==N-1){		//vecino de la derecha
				site[i].vec[1]=-1;
				site[i].frontera=1;
			} else{
				site[i].vec[1]=i+1;
			}

			site[i].vec[3]=-1;
			site[i].vec[2]=-1;
			site[i].vec[4]=-1;
			site[i].vec[5]=-1;
			site[i].x=i;
			site[i].y=0;
			site[i].z=0;

		}
	} else if(dim==2){
		for(i=0; i<N*N; i++){
			if(i<N){		//vecino de arriba
				site[i].vec[0]=-1;
				site[i].frontera=1;
			} else{
				site[i].vec[0]=i-N;
			}

			if((i+1)%N==0){		//vecino de la derecha
				site[i].vec[1]=-1;
				site[i].frontera=1;
			} else{
				site[i].vec[1]=i+1;
			}

			if(i<N*(N-1)){		//vecino de abajo
				site[i].vec[2]=i+N;
			} else{
				site[i].vec[2]=-1;
				site[i].frontera=1;
			}

			if(i%N==0){		//vecino de la izquierda
				site[i].vec[3]=-1;
				site[i].frontera=1;
			} else{
				site[i].vec[3]=i-1;
			}
			site[i].vec[4]=-1;
			site[i].vec[5]=-1;
			site[i].x=i%N;
			site[i].y=i/N;
			site[i].z=0;

		}
	} else if(dim==3){
		for(j=0; j<N;j++){
			for(i=0; i<N*N; i++){
				if(i<N){		//vecino de arriba
					site[j*N*N+i].vec[0]=-1;
					site[j*N*N+i].frontera=1;
				} else{
					site[j*N*N+i].vec[0]=j*N*N+i-N;
				}

				if((i+1)%N==0){		//vecino de la derecha
					site[j*N*N+i].vec[1]=-1;
					site[j*N*N+i].frontera=1;
				} else{
					site[j*N*N+i].vec[1]=j*N*N+i+1;
				}

				if(i<N*(N-1)){		//vecino de abajo
					site[j*N*N+i].vec[2]=j*N*N+i+N;
				} else{
					site[j*N*N+i].vec[2]=-1;
					site[j*N*N+i].frontera=1;
				}

				if(i%N==0){		//vecino de la izquierda
					site[j*N*N+i].vec[3]=-1;
					site[j*N*N+i].frontera=1;
				} else{
					site[j*N*N+i].vec[3]=j*N*N+i-1;
				}

				if(j==0){		//vecino enfrente
					site[j*N*N+i].vec[4]=-1;
					site[j*N*N+i].frontera=1;
				} else{
					site[j*N*N+i].vec[4]=j*N*N+i-N*N;
				}

				if(j==N-1){		//vecino atras
					site[j*N*N+i].vec[5]=-1;
					site[j*N*N+i].frontera=1;
				} else{
					site[j*N*N+i].vec[5]=j*N*N+i+N*N;
				}
				site[j*N*N+i].x=i%N;
				site[j*N*N+i].y=i/N;
				site[j*N*N+i].z=j;
			}
		}
	}

}

void deb(struct table *site){		//funcion para debuggear, i.e., imprime en terminal los valores guarados en los vecinos del arreglo site

	int i;
	for(i=0; i<vol; i++){

		cout<<i<<"\t"<<site[i].vec[0]<<"\t"<<site[i].vec[1]<<"\t"<<site[i].vec[2]<<"\t"<<site[i].vec[3]<<"\t"<<site[i].vec[4]<<"\t"<<site[i].vec[5]<<"\t"<<endl;
	}
	cout<<endl;
	for(i=0; i<vol; i++){

		cout<<i<<"\t"<<site[i].frontera<<"\t"<<site[i].n<<"\t"<<endl;
	}
	cout<<endl;
	for(i=0; i<vol; i++){

		cout<<i<<"\t"<<site[i].x<<"\t"<<site[i].y<<"\t"<<site[i].z<<"\t"<<endl;
	}
	cout<<endl;
}


int main(int argc, char *argv[])
{
        if(argc!=7) {
                printf("usage: %s <seed> <linear dimension> <spatial dimension> <Temperature> <n> <# of steps> \n", argv[0]);		//el programa espera 3 argumentos: 1- el nombre el programa,
                															//2- una semilla para el numero aleatorio, i.e., numero entero 3-la dimension
                exit(1);
        }
        myrand=(unsigned)atoi(argv[1]);
        N=(int)atoi(argv[2]);
        dim=(int)atoi(argv[3]);
        Tin=(double)atof(argv[4]);
        n=(int)atoi(argv[5]);
        steps=(int)atoi(argv[6]);

        myrand_back=myrand;
        Init_Random();

        struct table * site;
        if(dim==1){
        	vol=N;
        	site=(struct table *) calloc(N, sizeof(struct table));

        } else if(dim==2){
        	vol=N*N;
        	site=(struct table *) calloc(N*N, sizeof(struct table));

        } else if(dim==3){
        	vol=N*N*N;
        	site=(struct table *) calloc(N*N*N, sizeof(struct table));

        }

        initsys(site);
        deb(site);

        return 0;
}
