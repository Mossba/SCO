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
//Valor provisional para la constante de Boltzman hasta que se haga el analisis de unidades
#define k_Boltzman 1

//VARIABLES
unsigned myrand_back; //seed
int dim, N, vol, Tin, steps;
int n, neighbor=6;
double lambda, alpha, energy; 

struct table{
        int n;			//0,1,2 (spin or number of electrons in eg orbitals)
        int vec[6];			//sentido horario, empezando por las 12 hrs
        int frontera;		//1 si es frontera, 0 si no
		int x,y,z;
		double l_interaction;   //interaccion de larga distancia 
		double s_interaction;   //interaccion de corta distancia o de primeros vecinos
		int degeneration;       //degeneracion de cada estado

};

//FUNCIONES


//INICIALIZACION
void initsys(struct table *);
double energia(struct table *);
void actualizacion(struct table *);
void MonteCarlo(struct table *);
double metropolis(struct table *, double);
void concFunction(struct table *, int );
void writingFunction(struct table *);
double kappa_sq = sqrt((2 * (2 - lambda * neighbor))/(4 - lambda * neighbor));


void deb(struct table *, int *);

//FICHEROS
ofstream info("info.dat");
ofstream conc("SpinMedio.dat");
ofstream data("data.dat");


void initsys(struct table *site){		//se generan los vecinos en sentido horario
	int i,j;
	double R_ij;
	
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
			//Inicializacion de los valores del termino de interaccion en cada sitio y del spin de cada sitio
			for(int i = 0; i < vol; i++){
				site[i].l_interaction = 0;
				site[i].s_interaction = 0;
				site[i].n = 2;
				}

			//Inicializacion de la defeneracion de cada estado de spin
			for (int i = 0; i < vol; i++)
			{
				if (site[i].n == 2)
				{
					site[i].degeneration = 15;
				}

				if (site[i].n == 1)
				{
					site[i].degeneration = 9;
				}

				if (site[i].n == 0)
				{
					site[i].degeneration = 1;
				}
			}
			
			//Calculo de la interaccion de larga distancia del sitio i con todos los sitios del sistema
			for(int i = 0; i < vol; i++){
				for(int j = 0; j < vol && j != i; j++){
					R_ij = sqrt((site[i].x - site[j].x)*(site[i].x - site[j].x) + (site[i].y - site[j].y)*(site[i].y - site[j].y) + (site[i].z - site[j].z)*(site[i].z - site[j].z));
					site[i].l_interaction = site[i].l_interaction + exp(-kappa_sq * R_ij)/R_ij;
				}
			}

			//Calculo de la interaccion de corta distancia o tipo Ising 
			for (int i = 0; i < vol; i++)
			{
				if (site[i].vec[0] != -1) //Vecino de arriba
				{
					site[i].s_interaction = site[i].s_interaction + (alpha * alpha * lambda) * (site[i].n * site[site[i].vec[0]].n);
				}

				if (site[i].vec[1] != -1) //Vecino de la derecha
				{
					site[i].s_interaction = site[i].s_interaction + (alpha * alpha * lambda) * (site[i].n * site[site[i].vec[1]].n);
				}

				if (site[i].vec[2] != -1) //Vecino de abajo
				{
					site[i].s_interaction = site[i].s_interaction + (alpha * alpha * lambda) * (site[i].n * site[site[i].vec[2]].n);
				}

				if (site[i].vec[3] != -1) //Vecino de la izquierda
				{
					site[i].s_interaction = site[i].s_interaction + (alpha * alpha * lambda) * (site[i].n * site[site[i].vec[3]].n);
				}

				if (site[i].vec[4] != -1) //Vecino de enfrente
				{
					site[i].s_interaction = site[i].s_interaction + (alpha * alpha * lambda) * (site[i].n * site[site[i].vec[4]].n);
				}

				if (site[i].vec[5] != -1) //Vecino de arriba
				{
					site[i].s_interaction = site[i].s_interaction + (alpha * alpha * lambda) * (site[i].n * site[site[i].vec[5]].n);
				}
			}

}

//Inicializacion de la funcion para calcular la energia en cada sitio

double energia(struct table *site){
	double H = 0;

	for (int i = 0; i < vol; i++)
	{
		H = H + (energy * site[i].n - alpha * alpha * site[i].n * site[i].n) - site[i].s_interaction - site[i].l_interaction - (k_Boltzman * Tin) * log(site[i].degeneration);
	}

	return H;
}	

void actualizacion(struct table * site){ 	//Cada vez que el sistema se evoluciona, se debe calcular las energias de interaccion
											//Esta funcion hace eso, i.e., largo y corto alcance

			int i,j;
			double R_ij;

			//Calculo de la interaccion de larga distancia del sitio i con todos los sitios del sistema
			for(int i = 0; i < vol; i++){
				for(int j = 0; j < vol && j != i; j++){
					R_ij = sqrt((site[i].x - site[j].x)*(site[i].x - site[j].x) + (site[i].y - site[j].y)*(site[i].y - site[j].y) + (site[i].z - site[j].z)*(site[i].z - site[j].z));
					site[i].l_interaction = site[i].l_interaction + exp(-kappa_sq * R_ij)/R_ij;
				}
			}

			//Calculo de la interaccion de corta distancia o tipo Ising 
			for (int i = 0; i < vol; i++)
			{
				if (site[i].vec[0] != -1) //Vecino de arriba
				{
					site[i].s_interaction = site[i].s_interaction + (alpha * alpha * lambda) * (site[i].n * site[site[i].vec[0]].n);
				}

				if (site[i].vec[1] != -1) //Vecino de la derecha
				{
					site[i].s_interaction = site[i].s_interaction + (alpha * alpha * lambda) * (site[i].n * site[site[i].vec[1]].n);
				}

				if (site[i].vec[2] != -1) //Vecino de abajo
				{
					site[i].s_interaction = site[i].s_interaction + (alpha * alpha * lambda) * (site[i].n * site[site[i].vec[2]].n);
				}

				if (site[i].vec[3] != -1) //Vecino de la izquierda
				{
					site[i].s_interaction = site[i].s_interaction + (alpha * alpha * lambda) * (site[i].n * site[site[i].vec[3]].n);
				}

				if (site[i].vec[4] != -1) //Vecino de enfrente
				{
					site[i].s_interaction = site[i].s_interaction + (alpha * alpha * lambda) * (site[i].n * site[site[i].vec[4]].n);
				}

				if (site[i].vec[5] != -1) //Vecino de arriba
				{
					site[i].s_interaction = site[i].s_interaction + (alpha * alpha * lambda) * (site[i].n * site[site[i].vec[5]].n);
				}
			}

}	

void MonteCarlo(struct table * site){		/*Aqui hacemos Monte Carlo, es decir,
											  Medimos la energia
											  seleccionamos un sitio y le cambiamos el spin de s a t (ojo, hay 3 estados, entonces al variar el spin,
											  con probabilidad 0.5 se cambia de s a t_1 y con la prob complementaria se cambia de s a t_2 )
											  Actualizamos la estructura, o sea, las energia de interaccion
											  Medimos energia
											  Hacemos Metropolis y aceptamos o rechazamos esa variacion de spin*/

	int sitio=floor(FRANDOM*vol), valorSpinInicial;
	double c=FRANDOM, energy_before, energy_after, dE;
	energy_before=energia(site);
	valorSpinInicial=site[sitio].n;
	if(site[sitio].n==2){
		if(c<=0.5){
			site[sitio].n=0;
		} else{
			site[sitio].n=1;
		}
	}else if(site[sitio].n==1){
		if(c<=0.5){
			site[sitio].n=0;
		} else{
			site[sitio].n=2;
		}
	} else{
		if(c<=0.5){
			site[sitio].n=1;
		} else{
			site[sitio].n=2;
		}
	}
	actualizacion(site);
	energy_after=energia(site);
	dE=energy_after-energy_before;
	if(dE>0){					//YA SE HA CAMBIADO EL SPIN DEL SITE SITIO, POR TANTO, BASTA CORRER EL ALGORITMO DE METROPOLIS PARA VER SI DEBEMOS RECHAZAR LA NUEVA CONFIGURACION
		double met=metropolis(site, dE);
		if( met==0){
			site[sitio].n=valorSpinInicial;
		}
	}

}

double metropolis(struct table * site, double dE){

	double r, x=FRANDOM;
	if(x< exp(-dE/Tin)) r=1;		//SI X ES MENOR, SE ACEPTA LA NUEVA CONFIGURACION (se devuelve r=1). SI X ES MAYOR, SE RECHAZA Y SE VUELVE A LA CONFIGURACION INICIAL (se devuelve r=0)
	else r=0;

	return r;

}

///////////////////////////////////////////////////////DATA EXPORTING

void concFunction(struct table * site, int MCStep){			//Saca la concentracion a un archivo llamado SpinValorMedio.dat
	double s=0;
	for(int i=0; i<vol; i++){
		s=s+site[i].n;
	}
	s=s/vol;
	conc<<MCStep<<" "<<s<<endl;
}

void writingFunction(struct table *site){			//Saca fotos de la configuracion de spines. OJO CON ESTA FUNCION PARA MUCHOS PASOS DE MONTE CARLO Y SISTEMA GRANDES YA QUE SE AGOTA LA MEMORIA DE DISCO DROPBOX
	for(int i=0; i<vol; i++){

		data<<site[i].n<<" ";
	}
	data<<endl;
}

///////////////////////////////////////////////////////////////////////////

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

////////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[])
{
        if(argc!=9) {
                printf("usage: %s <seed> <linear dimension> <spatial dimension> <Temperature> <# of steps> <lambda> <alpha> <enrgy> \n" , argv[0]);		//el programa espera 3 argumentos: 1- el nombre el programa,
                															//2- una semilla para el numero aleatorio, i.e., numero entero 3-la dimension
                exit(1);
        }
        myrand=(unsigned)atoi(argv[1]);
        N=(int)atoi(argv[2]);
        dim=(int)atoi(argv[3]);
        Tin=(double)atof(argv[4]);
        steps=(int)atoi(argv[5]);
		lambda = (double)atof(argv[6]);
		alpha = (double)atof(argv[7]);
		energy = (double)atof(argv[8]);

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
        //deb(site);

        
        writingFunction(site);
        for(int i=0; i<steps; i++){
        	concFunction(site, i);		
        	MonteCarlo(site);
        	if(dim==2) writingFunction(site);
        }
		
	

        //Aqui se sacan los datos del input. Para ser organizados.
		info<<"seed \t linear dim \t spatial dim \t Temperature \t steps \t lambda \t alpha \t energy \t #of pictures \t snap"<<endl;
        info<<myrand<<"\t"<<N<<"\t"<<dim<<"\t"<<Tin<<"\t"<<steps<<"\t"<<lambda<<"\t"<<alpha<<"\t"<<energy<<endl;




        //Esto te genera un archivo que luego hay que convertirlo en un ejecutable y ejecutarlo. Lo que hace es organizar el archivo data.dat en matrices.
        //Cada matriz corresponde a la configuracion del sistema en un paso de Monte Carlo (ojo, solo se genera el archivo si estamos en 2D).
        // Luego, en gnuplot se hace la animacion escribiendo
        // set terminal gif animate delay 100
        // set output "animacion.gif"
        // do for [i=1:#dePasosDeMC] {
        // plot 'd2.txt' index i matrix with image
        // }
        //Mucho cuidado con poner muchos paso de MC, porque el archivo data.dat crece muy rapido

        if(dim==2){
        	ofstream mat("generamatriz");
      		mat<<"fold -s -w"<<2*N<<" data.dat > d.txt"<<endl;
        	mat<<"awk '{ print; } NR % "<<N<<"==0 {print \" \";}' d.txt > d2.txt"<<endl;
        	mat<<"rm data.dat"<<endl;
        	mat<<"rm d.txt"<<endl;
        	mat<<"mv d2.txt data.dat"<<endl;
        }
		
        return 0;
}
























