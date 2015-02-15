#include <math.h>       /* sqrt */
#include <vector>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <list>
#include <sys/time.h>
#include <bitset>
#include "lock.h"


using namespace std;

timeval start_step3, stop_step3,start_step4, stop_step4,start_step5, stop_step5, start_step6, stop_step6;
double elapsedTime_step3 = 0,elapsedTime_step4 = 0,elapsedTime_step5 = 0,elapsedTime_step6 = 0;

//FacorBase dichiarata qui cosi non bisogna calcolarsi ogni volta tutti i primi.
std::vector<int64_t> factorBase_vect;
long long primo = 0;

// Implementazione potenza per device
__device__ int64_t pow_cuda(int64_t base, long exponent) // ok testata
{
    int64_t result = 1;
    int64_t n = exponent;
    int64_t a = base;
    while(n>0){
      
      result *= a * (n%2==1)+ 1*(n%2 == 0);
      n /= 2;
      a *=a;
    }
     return result;
}

//Algoritmo di euclide per il GCD
//Fonte : http://en.wikipedia.org/wiki/Euclidean_algorithm
__device__ int64_t gcd_cuda ( int64_t a, int64_t b )
{
  int64_t c;
  while ( a != 0 ) {
     c = a; a = b%a;  b = c;
  }
  return b;
}

int64_t gcd ( int64_t a, int64_t b )
{
  int64_t c;
  while ( a != 0 ) {
     c = a; a = b%a;  b = c;
  }
  return b;
}

//calcola (a * b) % c senza mandare in overflow con l' operazione a*b/////
//Fonte = http://stackoverflow.com/questions/20971888/modular-multiplication-of-large-numbers-in-c
int64_t mulmod(int64_t  a, int64_t  b, int64_t  c)
{

     int64_t result = 0;
    a %= c;
    b %= c;
    while(b) {
        if(b & 0x1) {
            result += a;
            result %= c;
        }
        b >>= 1;
        a <<= 1;
        a %= c;
    }
    return result;
}

// Versione device
__device__ int64_t mulmod_cuda(int64_t  a, int64_t  b, int64_t  c)
{

     int64_t result = 0;
    a %= c;
    b %= c;
    while(b) {
        if(b & 0x1) {
            result += a;
            result %= c;
        }
        b >>= 1;
        a <<= 1;
        a %= c;
    }
    return result;
}

//modular exponentiation
// Info : http://en.wikipedia.org/wiki/Modular_exponentiation
int64_t modulo(int64_t base, int64_t exponent, int64_t mod)
{

    int64_t r = 1;

    while (true)
    {
        if (exponent % 2 == 1) r = mulmod(r, base, mod);
        exponent /= 2;
        if (exponent == 0) break;
        base = mulmod(base, base, mod);
    }
    return r;
}

// versione device
__device__ int64_t modulo_cuda(int64_t base, int64_t exponent, int64_t mod)
{

    int64_t r = 1;

    while (true)
    {
        if (exponent % 2 == 1) r = mulmod_cuda(r, base, mod);
        exponent /= 2;
        if (exponent == 0) break;
        base = mulmod_cuda(base, base, mod);
    }
    return r;
}

// modular inversion
// Fonte : http://rosettacode.org/wiki/Modular_inverse#C.2B.2B
int64_t mul_inv(int64_t a, int64_t b)
{

	int64_t b0 = b, t, q;
	int64_t x0 = 0, x1 = 1;
	if (b == 1) return 1;
	while (a > 1) {
		q = a / b;
		t = b;
        b = a % b;
        a = t;
		t = x0;
		int64_t temp = q * x0;
		x0 = x1 - temp;
		x1 = t;
	}
	if (x1 < 0) x1 += b0;
	return x1;
}

// versione Device
__device__ int64_t mul_inv_cuda(int64_t a, int64_t b)
{
	int64_t b0 = b, t, q;
	int64_t x0 = 0, x1 = 1;
	if (b == 1) return 1;
	while (a > 1) {
		q = a / b;
		t = b;
        b = a % b;
        a = t;
		t = x0;
		int64_t temp = q * x0;
		x0 = x1 - temp;
		x1 = t;
	}
	if (x1 < 0) x1 += b0;
	return x1;
}


//test Miller-Rabin
// Info : http://it.wikipedia.org/wiki/Test_di_Miller-Rabin
bool Miller(int64_t p,long iteration)
{
    if (p < 2){

        return false;
    }

    if (p != 2 && p % 2==0)
    {

        return false;
    }

    int64_t s = p - 1;
    while (s % 2 == 0)
    {
        s /= 2;
    }

    for (long i = 0; i < iteration; i++)
    {
        int64_t a = rand() % (p - 1) + 1, temp = s;
        int64_t mod = modulo(a, temp, p);
        while (temp != p - 1 && mod != 1 && mod != p - 1)
        {
            mod = mulmod(mod, mod, p);
            temp *= 2;
        }
        if (mod != p - 1 && temp % 2 == 0)
        {
            return false;
        }
    }
    return true;
}

// Calcola il successivo secondo il lemma di hensell
// Info : http://en.wikipedia.org/wiki/Hensel%27s_lemma
__device__ int64_t hensell_cuda(int64_t r1, int64_t  n, int64_t p)
{

    if (r1 == 0){ return 0; }

  
    int64_t temp;
    temp = mul_inv_cuda( r1*2, p );
    temp = (r1 - ( ( (r1*r1) - n) % p ) * temp ) % p;

    return temp;
}
// Test di Legendre
bool testLegendre(int64_t number, int64_t nRSA){

    if ( modulo(nRSA, (number-1)/2 , number) != 1 ){
        return false;
    }
    return true;

}
// Numero di legendre
long legendre(int64_t a, int64_t p)
{

  
  
      //if (p < 2)  // prime test is expensive.
        //throw new ArgumentOutOfRangeException("p", "p must not be < 2");
      if (a == 0)
      {
        return 0;
      }
      if (a == 1)
      {
        return 1;
      }
      int64_t result;
      if (a % 2 == 0)
      {
	result = legendre(a / 2, p);

	if (((p * p - 1) & 8) != 0) // instead of dividing by 8, shift the mask bit
        {
          result = -result;
        }
      }
      else
      {

        result = legendre(p % a, a);
        if (((a - 1) * (p - 1) & 4) != 0) // instead of dividing by 4, shift the mask bit
        {
          result = -result;
        }
      }
      return result;
    }
    
// versione CUDA    
__device__ long legendre_cuda(int64_t a, int64_t p)
{
      //if (p < 2)  // prime test is expensive.
        //throw new ArgumentOutOfRangeException("p", "p must not be < 2");
      if (a == 0)
      {
        return 0;
      }
      if (a == 1)
      {
        return 1;
      }
      int64_t result;
      if (a % 2 == 0)
      {
	result = legendre_cuda(a / 2, p);

	if (((p * p - 1) & 8) != 0) // instead of dividing by 8, shift the mask bit
        {
          result = -result;
        }
      }
      else
      {

        result = legendre_cuda(p % a, a);
        if (((a - 1) * (p - 1) & 4) != 0) // instead of dividing by 4, shift the mask bit
        {
          result = -result;
        }
      }
      return result;
    }

// Test bit = 0 o 1, di un certo numero.    
__device__ long test_bit_cuda(int64_t a, long pos_test) // ok testata
{
  
  long rest = 0;
  long result = a;
  
  for (long i = 0; i<= pos_test ; i++){
    
      rest = result % 2;
      result = result / 2;
      
      //if (i < *pos_test && result == 0) return -1;
  }
  return rest;
}

// Calcola il residuo QUadratico
// Info: http://it.wikipedia.org/wiki/Residuo_quadratico
//Fonte: Tradotta in c++ dalla versione mpz -> https://gmplib.org/list-archives/gmp-devel/2006-May/000633.html
__device__ long mpz_sqrtm_cuda(int64_t n, int64_t p) {  // ok testata

    int64_t q_int;
    int64_t w, n_inv, y;
    long i, s;

    if ( n % p == 0){
        q_int = 0;
        return q_int;
    }

     if(test_bit_cuda(p,1 )) {
	
	int64_t p_temp = (p+1) / 4;
      
        q_int = modulo_cuda(n, p_temp , p  ) ;
        return q_int;
    }

    q_int = p - 1;
    s = 0; 

    while(!test_bit_cuda(q_int, s )){
        s++;
    }
    
    q_int = q_int / pow_cuda(2,s);
    
    w = 2;

    while (legendre_cuda(w,p) != -1){
            w++;	    
    }

    w = modulo_cuda(w,q_int,p);
    q_int++;

    q_int = q_int / 2;
    
    q_int = modulo_cuda(n , q_int , p);
    
    n_inv = mul_inv_cuda(n,p);

    for(;;) {

      y = modulo_cuda(q_int,2,p);
       y = mulmod_cuda(y,n_inv,p);
	
        i = 0;
        
	while (y != 1){
            i++;
       	     y = modulo_cuda(y,2,p);
        }
        
        if (i == 0 ) {
        
	  return q_int;
        }
        
        if (s-i == 1){
            q_int = q_int * w;
        } else {
	    
            y = modulo_cuda(w,1 << (s-i-1), p);
            q_int = q_int * y;
        }

        q_int = q_int % p;
    }
    
}

// Parallelizzazione della Generazione delle relazioni
__global__ void genlistypsilon(int64_t *lyn_pun, int64_t *nRSA, int64_t *x_0,long *m){
  
    //*(C+threadIdx.x) = *(A+threadIdx.x) + *(B+threadIdx.x);
  long idx = blockIdx.x*blockDim.x + threadIdx.x;
 
    lyn_pun[idx] = (*x_0 + idx-*m) * (*x_0 + idx-*m) - *nRSA; 
}

// Parallelizzazione del processo del setaccio con il fattore 2 
__global__ void div2whileposs(int64_t *lyn_pun){
  
  long idx = blockIdx.x*blockDim.x + threadIdx.x;
  
  while( lyn_pun[idx] % 2 == 0){

     lyn_pun[idx] /=  2;

  }
  
}

// Parallelizzazione del processo del setaccio dal fattore 3 in su 
__global__ void div3up_while_poss(Lock lock, int64_t *factorBase_d, int64_t *lyn_pun,int64_t *nRSA, long *start ,int64_t *x_0,long *m, int64_t *limit_div){
  
  long idx = blockIdx.x*blockDim.x + threadIdx.x + *start;

  //cont_factor++;
  int64_t l_fact, p_origin;
  long exp;
  
  
  l_fact = factorBase_d[idx];
  p_origin =  factorBase_d[idx];
  exp = 1;

  int64_t s,h, s_calc, h_calc,s_calc_neg, h_calc_neg, s_calc_temp, h_calc_temp ;

  
  bool ok_sieve = true;
  while( ok_sieve ){
	      
            if (exp == 1){			// se il primo è p

	      
	      s = mpz_sqrtm_cuda(*nRSA, p_origin);
              h = l_fact - s;

              s_calc_temp = s - (*x_0 % l_fact);
              if (s_calc_temp < 0) {
                        s_calc_temp = l_fact + s_calc_temp;
               }

                h_calc_temp = (h - (*x_0 % l_fact)) % l_fact;
                if (h_calc_temp < 0) {
                        h_calc_temp = l_fact + h_calc_temp;
                }
         
	    }else{
	      				
	      
	      s= hensell_cuda(s, *nRSA, l_fact);
		
		h= l_fact - s;

                s_calc_temp = (s - (*x_0 % l_fact)) % l_fact;
                if (s_calc_temp < 0){
                    s_calc_temp = l_fact + s_calc_temp;
                }
                
                h_calc_temp = (h - (*x_0 % l_fact)) % l_fact;
                if (h_calc_temp < 0) {
                        h_calc_temp = l_fact + h_calc_temp;
                }
            }
            //cout << "yo6" << endl;
            int64_t k = 0;

            while(true){

                s_calc = s_calc_temp  + (l_fact * k );
                h_calc = h_calc_temp  + (l_fact * k);
                s_calc_neg = s_calc_temp + (l_fact * (-k));
                h_calc_neg = h_calc_temp + (l_fact * (-k));

               int64_t a = s_calc + *m, b = h_calc + *m, c = s_calc_neg + *m, d = h_calc_neg + *m;

	       
                if ((a < 0 or a >= 2*(*m)) and (b<0 or b >= 2 * (*m)) and (c < 0 or c >= 2*(*m)) and (d < 0 or d >= 2*(*m))){
                 //   cout <<"exitt whileeeeeeeee"<<endl;
                    break;
                }
                // applica divisione j positivo
                if (a >=0 and a < 2*(*m) ){

		    lock.lock();
		    lyn_pun[a] /= p_origin;;
		    lock.unlock();
                }
                if (a >= 0 and b < 2*(*m) and b != a){
                    
		  lock.lock();
		  lyn_pun[b] /= p_origin;
		  lock.unlock();
                }
                if (c >= 0 and c < 2*(*m) and c != b and c !=a){

                  lock.lock();  
		  lyn_pun[c] /= p_origin;
		  lock.unlock();

                }
                if (d >= 0 and d< 2*(*m) and d !=c and d!= b and d != a){

                  lock.lock();  
		  lyn_pun[d] /= p_origin;
		  lock.unlock();
                }

                k++;
            }
            exp++;                                                     // aumentiamo l esponente
    
            l_fact = pow_cuda(p_origin, exp);
	    if ( l_fact > *limit_div ) ok_sieve = false;
    
  }
}

// Parallelizzazione del processo di valutazione dei vettori soluzione
__global__ void evaluate_solutions(int64_t *nRSA_d, int64_t *jBuoni_d, long *jBuonicont_d, long *mat_sol_d, int64_t *listExponentOk_d, int64_t *result_d, int64_t* factorBase_d, long * factorBaseSize_d, size_t pitch, size_t pitch2, int64_t *x_0_d ){
  
	    long idx = blockIdx.x*blockDim.x + threadIdx.x;
  
            int64_t a = 1;
            int64_t b = 1;
	    int64_t *exponent;
	    exponent = (int64_t*)malloc(sizeof(int64_t)*(*factorBaseSize_d));

	    
            for (long i = 0; i <  *factorBaseSize_d; i++){
                exponent[i]=0;
            }
	    
	    for (long j = 0; j< *jBuonicont_d; j++){
		long* ms = (long*)((char*)mat_sol_d + idx*pitch2 );
		 if (ms[j] == 1){              // mat_sol[n_free_var*2][contRelationOk];
		  
            	     
		   a = mulmod_cuda(a , *x_0_d + jBuoni_d[j] , *nRSA_d);
		     
                    for (long i = 1; i < *factorBaseSize_d+1; i++){
		
		      int64_t* aa = (int64_t*)((char*)listExponentOk_d + j*pitch); 
                        exponent[i-1] = exponent[i-1] + aa[i];
                    
			
		    }
                }
            }

            for (long p = 0; p < *factorBaseSize_d; p++){

	       // b = (exponent[p]>0)*((b * modulo_cuda(factorBase_d[p], (exponent[p] / 2), *nRSA_d))  % *nRSA_d) + (exponent[p]<=0)*b;
	   
	       b = (exponent[p]>0)* mulmod_cuda(b, modulo_cuda(factorBase_d[p], (exponent[p] / 2), *nRSA_d),*nRSA_d)  + (exponent[p]<=0)*b;
	    }

            //STEP 8 = CALCOLO SOLUZIONE
            // Si calcola d=mcd(a-b,n) e se 1<d<n allora d è divisore non banale di n,
            // altrimenti si torna al passo 2) con una scelta di k più grande
         //   cout << "step 8" << endl;

            int64_t risultato_1, risultato_2, temp;
  
            temp = a+b;
	   
	    risultato_1 = gcd_cuda(temp,*nRSA_d);
	    
	    temp = a-b;
	    if (temp < 0) temp = -temp;
	    
	    risultato_2 = gcd_cuda(temp,*nRSA_d);

	   if ( risultato_1 != 1 and risultato_1 != *nRSA_d ){

                 *result_d = risultato_1;
		 
            }
           if ( risultato_2 != 1 and risultato_2 != *nRSA_d ){

                 *result_d = risultato_2;
	     
	  }
  
}




int64_t factorize(int64_t  nRSA, int64_t  x_0,  long k, long m){

    //a = numero rsa
    //x_0 = radice n rsa
    //k = numero factorBase
    //m = numero di relazioni

    /*puo restituire    -2 = trovati solo fattori banali -> aumentare k
                        -1 = meno relazione che fattori -> aumentare b
                        soluzione!
    */
    // STEP 1 IMMETTERE nRSA

    // STEP2 : SI SCEGLIE UN  K > 0, i numeri primi della factor base devono essere minori di k.
   cout << "step 2" << endl;

    // STEP3 : CREO FACTOR BASE,
    //   cout << "step 3" << endl;
    // Si esaminano tutti i primi p <= k e si eliminano tutti i primi dispari tali che  ( n/p ) != 1,  dove con ( n/p )
    // si intende il Simbolo di Legendre, e si ottiene così la base di fattori B = {p_1,p_2,....,p_t}
    gettimeofday(&start_step3, NULL);
    
    // ci servono per le allocazioni....
    long size_int64_t, size_long;
    size_long = sizeof(long);
    size_int64_t = sizeof(int64_t);
    
    
    long factorBaseSize = factorBase_vect.size();
    
    int64_t primo_temp = 0;
    if (factorBaseSize != 0){
        primo_temp =  factorBase_vect.at(factorBaseSize-1);
    }
    // aggiorno la factorbase
    while(true){

        while(true){

            primo_temp++;
            //cout << primo_temp << endl;
            if ( Miller(primo_temp, 20 ) ) break;

        }

        if (primo_temp > k) {
                break;
        }else if (testLegendre(primo_temp,nRSA)){
            primo = primo_temp;
            factorBase_vect.push_back(primo);
            //cout <<primo<< endl;
        }else{
            primo = primo_temp;
        }
    }
    
    factorBaseSize = factorBase_vect.size();
    int64_t *factorBase;
    factorBase=(int64_t*)malloc(size_int64_t*factorBaseSize);
    
    //cout << "FACTOR BASE" << endl;
    
    for (long p = 0; p< factorBaseSize; p++){
       // cout << "vect   "<< factorBase_vect.at(p) << endl;
	factorBase[p] = factorBase_vect.at(p);
	//cout << factorBase[p] << endl;
	
    }
    
    // copio factorbase su device
    int64_t *factorBase_d;
    cudaMalloc((void**)&factorBase_d, size_int64_t*factorBaseSize);
    cudaMemcpy(factorBase_d, factorBase, size_int64_t*factorBaseSize, cudaMemcpyHostToDevice);
    
    

    gettimeofday(&stop_step3, NULL);
    elapsedTime_step3 += (stop_step3.tv_sec - start_step3.tv_sec) * 1000.0;               // sec to ms
    elapsedTime_step3 += (stop_step3.tv_usec - start_step3.tv_usec) / 1000.0;

    // STEP4 : CREO LISTAYPSILONNUBMER E LI SETACCIO
    // Facendo assumere ad r valori interi successivi a sqrt(n), si trovano almeno t+1 valori y=r^2-n che abbiano tutti i loro fattori primi in B.
    cout << "step 4"<< endl;
    gettimeofday(&start_step4, NULL);
    
//////////// CUDA start/////////////////////////////////////////

    int64_t *listYpsilonNumber, *listYpsilonNumber_test,*listYpsilonNumber_d, *nRSA_punt ,*nRSA_d, *x_0_punt, *x_0_d, *limit_div, *limit_div_d;
    long *m_punt, *m_d;
    
    
    listYpsilonNumber = (int64_t*)malloc(size_int64_t*(2*m)); // alloco su host
    listYpsilonNumber_test = (int64_t*)malloc(size_int64_t*(2*m)); // alloco su host
    nRSA_punt = (int64_t*)malloc(size_int64_t);
    x_0_punt = (int64_t*)malloc(size_int64_t);
    limit_div = (int64_t*)malloc(size_int64_t);
    m_punt = (long*)malloc(size_long);
    
    cudaMalloc((void**)&listYpsilonNumber_d, size_int64_t*(2*m)); // alloco su device
    cudaMalloc((void**)&nRSA_d, size_int64_t);
    cudaMalloc((void**)&x_0_d, size_int64_t);
    cudaMalloc((void**)&limit_div_d, size_int64_t);
    cudaMalloc((void**)&m_d, size_long);
    
    *nRSA_punt = nRSA;
    *x_0_punt = x_0;
    *m_punt = m;
    
    
    cudaMemcpy(nRSA_d, nRSA_punt, size_int64_t, cudaMemcpyHostToDevice);
    cudaMemcpy(x_0_d, x_0_punt, size_int64_t, cudaMemcpyHostToDevice);
    cudaMemcpy(m_d, m_punt, size_long, cudaMemcpyHostToDevice);
 
   dim3 DimGrid(2*m,1); dim3 DimBlock(1,1,1); 
   genlistypsilon<<<DimGrid,DimBlock>>>(listYpsilonNumber_d,nRSA_d,x_0_d,m_d);
   
    cudaMemcpy(listYpsilonNumber, listYpsilonNumber_d, size_int64_t*(2*m), cudaMemcpyDeviceToHost);
    cudaMemcpy(listYpsilonNumber_test, listYpsilonNumber_d, size_int64_t*(2*m), cudaMemcpyDeviceToHost);
    
    cudaThreadSynchronize();
//////////// CUDA finish ///////////////////////////////////////////////

    *limit_div = (-(listYpsilonNumber[0]) /2  )+1;
    cudaMemcpy(limit_div_d, limit_div, size_int64_t, cudaMemcpyHostToDevice);

    long *start;
    start = (long*)malloc(sizeof(long));
    *start = 0;
    
    
    if (factorBase[0] == 2){
	cout << "caso2" << endl;
	
	//////////// CUDA start/////////////////////////////////////////
	
	dim3 DimGrid(2*m,1); dim3 DimBlock(1,1,1); 
	div2whileposs<<<DimGrid,DimBlock>>>(listYpsilonNumber_d);
	
        cudaMemcpy(listYpsilonNumber_test, listYpsilonNumber_d, size_int64_t*(2*m), cudaMemcpyDeviceToHost);
	
	//////////// CUDA finish////////////////////////////////////////
	
        *start = 1;
    }
    ///////////////// fin qui ok ////////////////////////
    cudaThreadSynchronize();

   
    //////////// CUDA start/////////////////////////////////////////
    
    Lock lock;
    
    long *start_d;
    cudaMalloc((void**)&start_d, sizeof(long));
    cudaMemcpy(start_d, start, sizeof(long), cudaMemcpyHostToDevice);
    
    cudaMemcpy(listYpsilonNumber_d,listYpsilonNumber_test,  size_int64_t*(2*m), cudaMemcpyHostToDevice); 

    if (*start == 0){
      dim3 DimGrid2(factorBaseSize-1,1); dim3 DimBlock2(1,1,1);
      div3up_while_poss<<<DimGrid2,DimBlock2>>>(lock,factorBase_d,listYpsilonNumber_d,nRSA_d,start_d,x_0_d,m_d,limit_div_d);
 
    }else{
      dim3 DimGrid2(factorBaseSize-1,1); dim3 DimBlock2(1,1,1);
      div3up_while_poss<<<DimGrid2,DimBlock2>>>(lock,factorBase_d,listYpsilonNumber_d,nRSA_d,start_d,x_0_d,m_d,limit_div_d);
 
    }  
   cudaThreadSynchronize();
    
   cudaMemcpy(listYpsilonNumber_test, listYpsilonNumber_d, size_int64_t*(2*m), cudaMemcpyDeviceToHost);

///////////////////////// CUDA finish ////////////////////////////////////////////   
 
  gettimeofday(&stop_step4, NULL);
    elapsedTime_step4 += (stop_step4.tv_sec - start_step4.tv_sec) * 1000.0;               // sec to ms
    elapsedTime_step4 += (stop_step4.tv_usec - start_step4.tv_usec) / 1000.0;

    // STEP 5 modulo A DUE DEGLI ESPONENTI!
    // Per ognuno dei valori y_1, y_2,...,y_{t+1} si calcola il vettore in Z2: di v_2(y_i)=(e_1,e_2,...,e_t)
    // dove e_i è la riduzione modulo 2 dell'esponente di p_i nella fattorizzazione di y_i.
   cout << "step 5" << endl;
    gettimeofday(&start_step5, NULL);

   long contRelationOk = 0;
    //#pragma omp parallel for
    for (long j = -m ; j < +m ; j++){
        if (listYpsilonNumber_test[j+m] == -1 or listYpsilonNumber_test[j+m] == 1){
     
            contRelationOk++; // mi serve per sapere la dimensione dell array ed essere comodo dopo..
      	}
    }
   
    // controllo che ci siano piu relazioni che numeri nella factorbase
    if (factorBaseSize+1 >  contRelationOk){
         
            cout << "ERROR! WE NEED MORE COLUMNS THAN ROWS" << endl;
        return -1;
    }
   
    int listExponentOk_mod2[contRelationOk][factorBaseSize+1];
    int64_t listExponentOk[contRelationOk][factorBaseSize+1];

    for (long x = 0; x <contRelationOk; x++){
        for (long y = 0; y <factorBaseSize+1; y++){
            listExponentOk_mod2[x][y] = 0;
            listExponentOk[x][y] = 0;

        }
    }
    // metto il segno nella prima linea, recuper gli esponenti, e mi segno i jnumber ok!!!!
    // creo un altra matrice con i soli valori che mi interessano, cioè coi jnumber che sono fattorizzabili con
    // la factor base, e già che ci sono gli faccio il modulo 2.

    // cout << "LISTA RELAZIONI ACCETTATE :" << endl;

    long jBuonicont = 0;
    std::vector<int64_t> jBuoni ;

    for (long j = -m ; j < +m ; j++){
        if (listYpsilonNumber_test[j+m] == -1 or listYpsilonNumber_test[j+m] == 1){
            if (listYpsilonNumber_test[j+m] == -1){
                listExponentOk_mod2[jBuonicont][0] = 1 ;// numero negativo
            }else{
                listExponentOk_mod2[jBuonicont][0] = 0 ;// numero positivo
            }
         
           for (long p = 0; p < factorBaseSize; p++){

                int64_t fact_temp = 1;
                bool ok = true;
                do{
                          //  cout << "5c" << endl;
                    int64_t pmpz = factorBase[p];
		    fact_temp = fact_temp * pmpz;
                    if( (listYpsilonNumber[j+m] % fact_temp ) == 0){
                        listExponentOk[jBuonicont][p+1]++;

                    }else{
                        ok = false;
                    }
                }while (ok);
                listExponentOk_mod2[jBuonicont][p+1] = listExponentOk[jBuonicont][p+1];
                listExponentOk_mod2[jBuonicont][p+1] = listExponentOk_mod2[jBuonicont][p+1] % 2;

                    //cont_factor++;
            }

            jBuonicont++;
            jBuoni.push_back(j);
	  }
    }

    int64_t *jBuoni_arr;
    jBuoni_arr = (int64_t*)malloc(sizeof(int64_t)*jBuoni.size());
    
    
    for (long i = 0; i < jBuoni.size(); i++){
      
     jBuoni_arr[i] = jBuoni.at(i);
    }
    
    
       // cout << endl;
    gettimeofday(&stop_step5, NULL);
    elapsedTime_step5 += (stop_step5.tv_sec - start_step5.tv_sec) * 1000.0;               // sec to ms
    elapsedTime_step5 += (stop_step5.tv_usec - start_step5.tv_usec) / 1000.0;


     // DA PARALLELIZZARE
    // STEP 6 METODO DI GAUSS
    // Con il metodo di eliminazione di Gauss si determinano alcuni dei vettori v_2(y_i) che danno somma uguale al vettore nullo
     gettimeofday(&start_step6, NULL);
 
    long vect_solution[contRelationOk]; // 0 var libere, -1 da valutare
    long vect_pivot[contRelationOk];  //   -1 no pivot, n pos riga
    long n_free_var = 0;

    for (long i = 0; i< contRelationOk; i++){
        vect_solution[i] = 0;
        vect_pivot[i] = -1;
    }

    long cont_rig = 0;
    for(long j = 0; j < contRelationOk; j++){ // colonne
        long row_i = -1;  // row i = j

        for(long i = cont_rig; i< factorBaseSize+1; i++ ){ // righe, cerco Aij = 1
            // cerco la riga con l' 1 sulla colonna di adesso
            if(listExponentOk_mod2[j][i] == 1){
                row_i = i;
                    //cout << "i =" <<  i ;
                break;
            }
        }

        if (row_i != -1){
                // salvo row i nella matrice risultato
                //cout << "cambio riga " << row_i <<" con riga "<< cont_rig << endl ;
            vect_solution[j] = -1;
            vect_pivot[j] = cont_rig;

            if ( cont_rig == row_i){

                cont_rig++;
            }else{
                // devo scambiarle

                for (long z = 0; z < contRelationOk; z++){
      
                    int64_t temp = listExponentOk_mod2[z][row_i];
                    listExponentOk_mod2[z][row_i] = listExponentOk_mod2[z][cont_rig];
                    listExponentOk_mod2[z][cont_rig] = temp;
                }
                cont_rig++;
            }

            // cerco le righe che hanno 1 nella colonna di adesso e le sommo alla riga selezionata
            for(long k = cont_rig; k <  factorBaseSize+1; k++){
                if (listExponentOk_mod2[j][k] == 1){
                    for (long cc = 0; cc< contRelationOk; cc++){

                        listExponentOk_mod2[cc][k] = (listExponentOk_mod2[cc][k] + listExponentOk_mod2[cc][cont_rig-1]) % 2 ;
                    }

                }
            }
         }else{

             n_free_var++;
        }
    }


    // TROVARE LA SOLUZIONE ( vettore spazio nullo) valutando le x dal basso!!!!
   long mat_sol[n_free_var*2][contRelationOk];

   for (long s = 0; s < n_free_var; s++){

    	for (long i = 0; i< contRelationOk; i++){

	  mat_sol[s][i] = vect_solution[i];
    	}


    long myvar = s;
    for (long i = 0; i < contRelationOk; i++){

        if (vect_pivot[i] == -1){

            if (myvar == 0){

                    mat_sol[s][i] = 1;
                    break;
            }else{
                myvar--;
            }
        }
    }


    for (long i = contRelationOk-1; i >= 0; i--){
            // estrago n pivot
            long pivot = vect_pivot[i];
            if (pivot != -1){ // non è var libera

                mat_sol[s][i] = 0;
                for (long j = i+1; j < contRelationOk; j++){
                    if (listExponentOk_mod2[j][pivot] == 1){
                        mat_sol[s][i] =  (mat_sol[s][i] + mat_sol[s][j]) % 2;

                }
            }
        }
    }
  }
  
  for (long s = n_free_var; s < n_free_var*2; s++){

    for (long i = 0; i< contRelationOk; i++){


        mat_sol[s][i] = (mat_sol[s-n_free_var][i]+1) %2;
    }
  }

  ////////////////// fino a qui identico////////////////////////

  //Puntatori adhoccuda
  long * factorBaseSize_punt;
  factorBaseSize_punt = (long*)malloc(sizeof(long));
  *factorBaseSize_punt = factorBaseSize;
  
  long * jBuonicont_punt;
  jBuonicont_punt = (long*)malloc(sizeof(long));
  *jBuonicont_punt = jBuonicont;
  
  
  // factorbase, jbuoni, mat_sol, listexponentOK, result.
  long *factorBaseSize_d ;
  int64_t *jbuoni_d, *listExponentOk_d, *result_d, *result; 
  long *jBuonicont_d,*mat_sol_d;
  
  size_t pitch, pitch2;  // pitch listexp, pitch2 mat_sol
  
  cudaMalloc((void**)&factorBaseSize_d, size_long);
  cudaMemcpy(factorBaseSize_d , factorBaseSize_punt , size_long, cudaMemcpyHostToDevice);
  
  cudaMalloc((void**)&jbuoni_d,jBuonicont*sizeof(int64_t));
  cudaMemcpy(jbuoni_d , jBuoni_arr , jBuonicont*sizeof(int64_t), cudaMemcpyHostToDevice);
  
  //long mat_sol[n_free_var*2][contRelationOk];
  cudaMallocPitch(&mat_sol_d, &pitch2, sizeof(long)*contRelationOk,n_free_var*2);
  cudaMemcpy2D(mat_sol_d,pitch2,mat_sol,contRelationOk*sizeof(long),contRelationOk*sizeof(long),n_free_var*2, cudaMemcpyHostToDevice);
  
  //int64_t listExponentOk[contRelationOk][factorBaseSize+1];
  cudaMallocPitch(&listExponentOk_d, &pitch, size_int64_t*(factorBaseSize+1),contRelationOk);
  cudaMemcpy2D(listExponentOk_d,pitch,listExponentOk,(factorBaseSize+1)*sizeof(int64_t),(factorBaseSize+1)*sizeof(int64_t),contRelationOk, cudaMemcpyHostToDevice);
  
  result = (int64_t*)malloc(sizeof(int64_t));
  *result = 0;
  cudaMalloc((void**)&result_d, size_int64_t);
  cudaMemcpy(result_d , result, sizeof(int64_t), cudaMemcpyHostToDevice);
  
  cudaMalloc((void**)&jBuonicont_d, sizeof(long));
  cudaMemcpy(jBuonicont_d , jBuonicont_punt, sizeof(long), cudaMemcpyHostToDevice);
  
  //dim3 DimGrid0(n_free_var*2,1); dim3 DimBlock0(1,1,1);
  evaluate_solutions<<<n_free_var*2,1>>>( nRSA_d, jbuoni_d, jBuonicont_d, mat_sol_d, listExponentOk_d, result_d, factorBase_d, factorBaseSize_d, pitch, pitch2, x_0_d);


  cudaMemcpy(result , result_d , size_int64_t, cudaMemcpyDeviceToHost);
  cudaDeviceReset();
  
  
    cout << "result = " << *result <<  endl;
    
    
//////////////////////////// cuda finish ////////////////////////////////////////
    
    gettimeofday(&stop_step6, NULL);
    elapsedTime_step6 += (stop_step6.tv_sec - start_step6.tv_sec) * 1000.0;               // sec to ms
    elapsedTime_step6 += (stop_step6.tv_usec - start_step6.tv_usec) / 1000.0;
    
    if (*result != 0 ){
        
        return *result;
    }

    return -2;
   


}
int main(int argc, char *argv[])
{
timeval start_tot, stop_tot;
double elapsedTime_tot;
gettimeofday(&start_tot, NULL);
    // STEP1 : Viene dato in input il numero naturale dispari n>1.
    //long nRSA =  argc[1];

      int64_t nRSA, result, x_0;

    if ( argc > 1 ) {

	nRSA = atoll(argv[1]);
	
	 x_0 = sqrt(nRSA);

        long k = 10;
        long m = 100;

        bool soluzione_trovata = false;

        long cont_meno1 = 0;
        long cont_meno2 = 0;
        while(!soluzione_trovata){
	    result = factorize(nRSA,x_0,k,m);
            if ( result == -2){
                k = k+10;
                m=  m+30;
          
                cont_meno2++;
                cout << "k =" << k << ", m = " << m << endl;
            }
            else if ( result == -1){
               
                k=k+10;
                m=m+10;
                cont_meno1++;
           
                cout << "k =" << k << ", m = " << m << endl;
            
            	
            }else{

                soluzione_trovata = true;
            }

        }
        gettimeofday(&stop_tot, NULL);
        elapsedTime_tot = (stop_tot.tv_sec - start_tot.tv_sec) * 1000.0;               // sec to ms
        elapsedTime_tot += (stop_tot.tv_usec - start_tot.tv_usec) / 1000.0;            // us to ms


        cout << "k = " << k << ", m = " << m << endl;
        cout << "-1 volte = " << cont_meno1 << ", -2 volte = " << cont_meno2 << endl;
        cout << endl << "nRSA = " << result << " x " << nRSA/result <<endl;
        cout << "nrsa = " << nRSA << endl;
	cout << "nrsa ris = "<< result * (nRSA / result) << endl;
        cout << " t-totale = " << elapsedTime_tot << " ms.\n" << endl;
        cout << " t-step3 = " << elapsedTime_step3 << " ms.\n" << endl;
        cout << " t-step4 = " << elapsedTime_step4 << " ms.\n" << endl;
        cout << " t-step5 = " << elapsedTime_step5<< " ms.\n" << endl;
        cout << " t-step6 = " << elapsedTime_step6 << " ms.\n" << endl;
	cout << "[QSieve, ver. CUDA]" << endl;
        cout << endl;


    }else{

        cout << "nessun argomento immesso!!!" << endl;
    }
    return 0;
}
