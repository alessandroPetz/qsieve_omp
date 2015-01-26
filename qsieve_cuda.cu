#include <math.h>       /* sqrt */
#include <vector>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <list>
#include <sys/time.h>
#include <bitset>
#include "lock.h"



#include "cuda_runtime.h"

#include "device_launch_parameters.h"


#define BigInt long long


using namespace std;

timeval start_step3, stop_step3,start_step4, stop_step4,start_step5, stop_step5, start_step6, stop_step6;
double elapsedTime_step3 = 0,elapsedTime_step4 = 0,elapsedTime_step5 = 0,elapsedTime_step6 = 0;

//FacorBase dichiarata qui cosi non bisogna calcolarsi ogni volta tutti i primi.
std::vector<BigInt> factorBase_vect ;
long long primo = 0;

/*
// modular moltiplication without overflow
__device__ BigInt modmup( BigInt a, BigInt b, BigInt m )
{
    if(a>m)
        a=a%m;
    if(b>m)
        b=b%m;
    BigInt ret = 0;
    BigInt temp = 0;
    if(a<b)
    {
      temp = a;
      a = b;
      b = temp;
    }
    while(b)
    {
        while(a<m)
        {
            if(b&1)
                ret += a;
            a<<=1;
            b>>=1;
        }
        a-=m;
        while(ret>=m)
            ret-=m;
        if(a<b)
        {
	  temp = a;
	  a = b;
	  b = temp;
	}
    }
    return ret;
}

BigInt ( BigInt a, BigInt b, BigInt m )
{
    if(a>m)
        a=a%m;
    if(b>m)
        b=b%m;
    BigInt ret = 0;
    BigInt temp = 0;
    if(a<b)
    {
      temp = a;
      a = b;
      b = temp;
    }
    while(b)
    {
        while(a<m)
        {
            if(b&1)
                ret += a;
            a<<=1;
            b>>=1;
        }
        a-=m;
        while(ret>=m)
            ret-=m;
        if(a<b)
        {
	  temp = a;
	  a = b;
	  b = temp;
	}
    }
    return ret;
}
*/

///////// power cuda /////////////
__device__ BigInt pow_cuda(BigInt base, int exponent) // ok testata
{
    BigInt result = base;
    for (int i = 1; i < exponent; i++ ){
      
      result = base * result;
      
    }  
     return result;
}

/////////////gcd cuda ///////////////////////

__device__ BigInt gcd_cuda2( BigInt a, BigInt b )
{
  int c;
  while ( a != 0 ) {
     c = a; a = b%a;  b = c;
  }
  return b;
}


__device__ BigInt gcd_cuda(BigInt u, BigInt v)
{
  int shift;
 
  /* GCD(0,v) == v; GCD(u,0) == u, GCD(0,0) == 0 */
  if (u == 0) return v;
  if (v == 0) return u;
 
  /* Let shift := lg K, where K is the greatest power of 2
        dividing both u and v. */
  for (shift = 0; ((u | v) & 1) == 0; ++shift) {
         u >>= 1;
         v >>= 1;
  }
 
  while ((u & 1) == 0)
    u >>= 1;
 
  /* From here on, u is always odd. */
  do {
       /* remove all factors of 2 in v -- they are not common */
       /*   note: v is not zero, so while will terminate */
       while ((v & 1) == 0)  /* Loop X */
           v >>= 1;
 
       /* Now u and v are both odd. Swap if necessary so u <= v,
          then set v = v - u (which is even). For bignums, the
          swapping is just pointer movement, and the subtraction
          can be done in-place. */
       if (u > v) {
         BigInt t = v; v = u; u = t;}  // Swap u and v.
       v = v - u;                       // Here v >= u.
     } while (v != 0);
 
  /* restore common factors of 2 */
  return u << shift;
}


//calculates (a * b) % c taking into account that a * b might overflow/////
BigInt mulmod(BigInt  a, BigInt  b, BigInt  c)
{

     BigInt result = 0;
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

__device__ BigInt mulmod_cuda(BigInt  a, BigInt  b, BigInt  c)
{

     BigInt result = 0;
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



//modular exponentiation///////////////////////////////////////////////////////
BigInt modulo(BigInt base, BigInt exponent, BigInt mod)
{

    BigInt x = 1;
    BigInt y = base;

    while (exponent > 0)
    {
        if (exponent % 2 == 1)
            x = mulmod(x, y, mod);

        y = mulmod(y,y,mod);
        exponent = exponent / 2;
    }
    return x % mod;
}

__device__ BigInt modulo_cuda(BigInt base, BigInt exponent, BigInt mod)
{

    BigInt x = 1;
    BigInt y = base;

    while (exponent > 0)
    {
        if (exponent % 2 == 1)
            x = mulmod_cuda(x, y, mod);

        y = mulmod_cuda(y,y,mod);
        exponent = exponent / 2;
    }
    return x % mod;
}

// modular inversion
BigInt mul_inv(BigInt a, BigInt b)
{

	BigInt b0 = b, t, q;
	BigInt x0 = 0, x1 = 1;
	if (b == 1) return 1;
	while (a > 1) {
		q = a / b;
		t = b;
        b = a % b;
        a = t;
		t = x0;
		BigInt temp = q * x0;
		x0 = x1 - temp;
		x1 = t;
	}
	if (x1 < 0) x1 += b0;
	return x1;
}

// modular inversion
__device__ BigInt mul_inv_cuda(BigInt a, BigInt b)
{
	BigInt b0 = b, t, q;
	BigInt x0 = 0, x1 = 1;
	if (b == 1) return 1;
	while (a > 1) {
		q = a / b;
		t = b;
        b = a % b;
        a = t;
		t = x0;
		BigInt temp = q * x0;
		x0 = x1 - temp;
		x1 = t;
	}
	if (x1 < 0) x1 += b0;
	return x1;
}


//test Miller-Rabin
bool Miller(BigInt p,int iteration)
{
    if (p < 2){

        return false;
    }

    if (p != 2 && p % 2==0)
    {

        return false;
    }

    BigInt s = p - 1;
    while (s % 2 == 0)
    {
        s /= 2;
    }

    for (int i = 0; i < iteration; i++)
    {
        BigInt a = rand() % (p - 1) + 1, temp = s;
        BigInt mod = modulo(a, temp, p);
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

__device__ BigInt hensell_cuda(BigInt r1, BigInt  n, BigInt p)
{

    if (r1 == 0){ return 0; }

    //return (r1 - ( ( (r1*r1-n) % p ) * mul_inv( r1*2, p )  % p )) % p;

    BigInt temp;
    //BigInt temp2 = r1*2;
    //BigInt expneg = -1;

    //mpz_powm (temp.get_mpz_t(), temp2.get_mpz_t() , expneg.get_mpz_t(), p.get_mpz_t() );
    //cout << "inv( " << r1*2 << ", "<< p << ") = " endl;
    temp = mul_inv_cuda( r1*2, p );
   //r2 = r1 - (r1^2 - n )/ (2r1);
    //temp = modulo(r1*2, -1, p);
    //cout << "temp = " << temp << endl;
    //return temp;
    //return (r1 - ( ( (r1*r1) - n) % p ) * temp ) % p;
    temp = (r1 - ( ( (r1*r1) - n) % p ) * temp ) % p;

    return temp;
}

bool testLegendre(BigInt number, BigInt nRSA){

    //mpz_class a;
    //mpz_class b = (number-1)/2;
    //mpz_class numbermpz = number;

    //mpz_powm(a.get_mpz_t(),nRSA.get_mpz_t(),b.get_mpz_t(), number.get_mpz_t());

    if ( modulo(nRSA, (number-1)/2 , number) != 1 ){
        return false;
    }
    return true;

}



int legendre(BigInt a, BigInt p)
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
      BigInt result;
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
__device__ int legendre_cuda(BigInt a, BigInt p)
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
      BigInt result;
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

    
__device__ int test_bit_cuda(BigInt a, int pos_test) // ok testata
{
  
  // 18
  // 10010
  // 18 / 2 = 9 % 0
  // 9 / 2 = 4 % 1
  // 4 % 2 = 2 % 0
  // 2 % 2 = 1 % 0
  // 1 % 2 = 0 % 1
  
  int rest = 0;
  int result = a;
  
  for (int i = 0; i<= pos_test ; i++){
    
      rest = result % 2;
      result = result / 2;
      
      //if (i < *pos_test && result == 0) return -1;
  }
  return rest;
}

__device__ int mpz_sqrtm_cuda(BigInt n, BigInt p) {  // ok testata

    BigInt q_int;
    BigInt w, n_inv, y;
    int i, s;

    if ( n % p == 0){
        q_int = 0;
        return q_int;
    }

     if(test_bit_cuda(p,1 )) {
	
	BigInt p_temp = (p+1) / 4;
      
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

        //y = y * n_inv;
	//y = y % p;

	 y = mulmod_cuda(y,n_inv,p);
	
        i = 0;
        
	while (y != 1){
            i++;
            //y = (y * y) % p;
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

  ///////okkkkkkkkkkkkkkkkkkkk///////////////
__global__ void genlistypsilon(BigInt *lyn_pun, BigInt *nRSA, BigInt *x_0,long *m){
  
    //*(C+threadIdx.x) = *(A+threadIdx.x) + *(B+threadIdx.x);
  int idx = blockIdx.x*blockDim.x + threadIdx.x;

    
    lyn_pun[idx] = (*x_0 + idx-*m) * (*x_0 + idx-*m) - *nRSA;
}

///////okkkkkkkkkkkkkkkkkkkk///////////////
__global__ void div2whileposs(BigInt *lyn_pun){
  
  int idx = blockIdx.x*blockDim.x + threadIdx.x;
  
  while( lyn_pun[idx] % 2 == 0){

     lyn_pun[idx] /=  2;

  }
  
}
//// okkkkkkkkkkkkkkkkk /////////////////////////////////////////
//__global__ void div3up_while_poss( BigInt *factorBase_d, BigInt *lyn_pun,BigInt *nRSA, int *start ,BigInt *x_0,long *m){
__global__ void div3up_while_poss(Lock lock, BigInt *factorBase_d, BigInt *lyn_pun,BigInt *nRSA, int *start ,BigInt *x_0,long *m){
  
  int idx = blockIdx.x*blockDim.x + threadIdx.x + *start;

  //cont_factor++;
  BigInt l_fact, p_origin;
  int exp;
  
  
  l_fact = factorBase_d[idx];
  p_origin =  factorBase_d[idx];
  exp = 1;

  BigInt s,h, s_calc, h_calc,s_calc_neg, h_calc_neg, s_calc_temp, h_calc_temp ;

  
  bool ok_sieve = true;
  while( ok_sieve ){

            l_fact = pow_cuda(p_origin, exp);
	 
	    if ( l_fact > 2 * *m ) ok_sieve = false;
	      
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
            BigInt k = 0;

            while(true){

                s_calc = s_calc_temp  + (l_fact * k );
                h_calc = h_calc_temp  + (l_fact * k);
                s_calc_neg = s_calc_temp + (l_fact * (-k));
                h_calc_neg = h_calc_temp + (l_fact * (-k));

               BigInt a = s_calc + *m, b = h_calc + *m, c = s_calc_neg + *m, d = h_calc_neg + *m;

	       
                if ((a < 0 or a >= 2*(*m)) and (b<0 or b >= 2 * (*m)) and (c < 0 or c >= 2*(*m)) and (d < 0 or d >= 2*(*m))){
                 //   cout <<"exitt whileeeeeeeee"<<endl;
                    break;
                }
                
                // applica divisione j positivo
                if (a >=0 and a < 2*(*m) ){
                    //#pragma omp critical
		    lock.lock();
		    lyn_pun[a] /= p_origin;;
		    lock.unlock();
                }
                if (a >= 0 and b < 2*(*m) and b != a){
                    //#pragma omp critical
                    
		  lock.lock();
		  lyn_pun[b] /= p_origin;
		  lock.unlock();
                }
                if (c >= 0 and c < 2*(*m) and c != b and c !=a){
                    //#pragma omp critical
                  lock.lock();  
		  lyn_pun[c] /= p_origin;
		  lock.unlock();

                }
                if (d >= 0 and d< 2*(*m) and d !=c and d!= b and d != a){
                    //#pragma omp critical
                  lock.lock();  
		  lyn_pun[d] /= p_origin;
		  lock.unlock();
                }

                k++;
            }
            exp++;                                                     // aumentiamo l esponente
    }
}

__global__ void evaluate_solutions( BigInt *nRSA_d, BigInt *factorBase_d, BigInt *x_0_d , BigInt *jBuoni_d ,BigInt *listExponentOk_d, size_t pitch , int* listExponentOk_mod2_d, size_t pitch2 ,int *vect_solution_d, int *vect_pivot_d, BigInt *result_d, long* contRelationOk_d,long* factorBaseSize_d,long* jBuonicont_d)
{
  
  //for(long s = 0; s < n_free_var; s++){
    
    int idx = blockIdx.x*blockDim.x + threadIdx.x;
  
        //cout << "numeri free var" <<n_free_var << endl;
        // copia locale di vect solution... ???? è una soluzione buona? pare di si..
	int *vect_solution_local;
	vect_solution_local = (int*)malloc(sizeof(int)*(*contRelationOk_d));
	
	for (long i = 0; i< *contRelationOk_d; i++){
            
	    vect_solution_local[i] = vect_solution_d[i];

        }

        long myvar = idx;
        for (long i = 0; i< *contRelationOk_d; i++){
            if ( vect_pivot_d[i] == -1){

                if (myvar == 0){

                    vect_solution_local[i] = 1;
                    break;
                }else{

                    myvar--;
                }
            }
        }
      /*  
        for (long i = contRelationOk-1; i >= 0; i--){
            // estrago n pivot
            long pivot = vect_pivot[i];
            if (pivot != -1){ // non è var libera

                vect_solution_local[i] = 0;
                for (long j = i+1; j < contRelationOk; j++){
                    if (listExponentOk_mod2[j][pivot] == 1){
                        vect_solution_local[i] =  (vect_solution_local[i] + vect_solution_local[j]) % 2;
                    }
                }
            }
        }
      */
        // qui valuto le x in base alla variabile libera assegnata!
        // da migliorare
        for (long i = (*contRelationOk_d-1); i >= 0; i--){
            // estrago n pivot
            long pivot = vect_pivot_d[i];
            if (pivot != -1){ // non è var libera

               *(vect_solution_local+i) = 0;
	      
	       for (long j = i+1; j < *contRelationOk_d; j++){
		    
		    int* aa = (int*)((char*)listExponentOk_mod2_d + j*pitch2);
		    //if (listExponentOk_mod2[j][pivot] == 1){
                      if( aa[pivot] == 1){
			
			vect_solution_local[i] =  (vect_solution_local[i] + vect_solution_local[j]) % 2;
                    }
                }
            }
        }

        for (long v = 0; v < 2; v++){  // lo facciamo 2 volte, uno per il vettore, e l altro per il suo inverso

            if (v == 1){

                 for (long i = 0; i< *contRelationOk_d; i++){
                    if ( vect_solution_local[i]== 0){
                       vect_solution_local[i] = 1;
                    }else{
                       vect_solution_local[i] = 0;
                    }
                }

            }

        // STEP 7 = calcolo di a e b
        // Si pone a uguale al prodotto degli r_i corrispondenti agli y_i trovati nel passo 6)
        // e si pone b uguale al prodotto delle potenze di p_1,p_2,...,p_t con esponenti uguali
        // alla semisomma degli esponenti della fattorizzazione degli stessi y_i
	    
            BigInt a = 1;
            BigInt b = 1;
            long cont = 0;
            //vector<BigInt>::iterator p2;
      //cout << "step 7" << endl;
	    
	    BigInt *exponent;
	    exponent = (BigInt*)malloc(sizeof(BigInt)*(*factorBaseSize_d));
	    
            //BigInt exponent[*factorBaseSize_d];
            for (long i = 0; i <  *factorBaseSize_d; i++){
                exponent[i]=0;
            }
	    
	    for (int i = 0; i< *jBuonicont_d; i++){
                if ( vect_solution_local[cont] == 1){

                    a = (a * (*x_0_d + jBuoni_d[i])) % *nRSA_d;
		    
                    for (long i = 1; i < *factorBaseSize_d+1; i++){
		      
		      //exponent[i-1] = exponent[i-1] + listExponentOk[cont][i];
			BigInt* aa = (BigInt*)((char*)listExponentOk_d + cont*pitch); 
                        exponent[i-1] = exponent[i-1] + aa[i];
                    
			
		    }
                }
                cont++;
            }
            cont=0;
            for (long p = 0; p < *factorBaseSize_d; p++){

                if (exponent[cont]>0){

                    BigInt p3mpz = factorBase_d[p];
                    BigInt temp = exponent[cont] * 2;
                   
		    temp = modulo_cuda(p3mpz, temp, *nRSA_d);
                    temp = temp * b;
                    b = temp % *nRSA_d;
                 
                }
                cont++;
            }

            //STEP 8 = CALCOLO SOLUZIONE
            // Si calcola d=mcd(a-b,n) e se 1<d<n allora d è divisore non banale di n,
            // altrimenti si torna al passo 2) con una scelta di k più grande
         //   cout << "step 8" << endl;

            BigInt risultato_1, risultato_2, temp;
  
            temp = abs(a+b);
	   
	    risultato_1 = gcd_cuda(temp,*nRSA_d);
	    
	    temp = abs(a-b);
            
	    risultato_2 = gcd_cuda(temp,*nRSA_d);

	   if ( risultato_1 != 1 and risultato_1 != *nRSA_d ){

                 *result_d = risultato_1;
		 
            }
           if ( risultato_2 != 1 and risultato_2 != *nRSA_d ){

                 *result_d = risultato_2;
	     
	  }
  
	}

	
}


BigInt factorize(BigInt  nRSA, BigInt  x_0,  long k, long m){

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
    int size_BigInt, size_long;
    size_long = sizeof(long);
    size_BigInt = sizeof(BigInt);
    
    long factorBaseSize = factorBase_vect.size();
    
    BigInt primo_temp = 0;
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
    BigInt *factorBase;
    factorBase=(BigInt*)malloc(size_BigInt*factorBaseSize);
    
    //cout << "FACTOR BASE" << endl;
    
    for (long p = 0; p< factorBaseSize; p++){
       // cout << "vect   "<< factorBase_vect.at(p) << endl;
	*(factorBase+p) = factorBase_vect.at(p);
//	cout << " array"<< *(factorBase+p) << endl;
	
    }
    
    // copio factorbase su device
    BigInt *factorBase_d;
    cudaMalloc((void**)&factorBase_d, size_BigInt*factorBaseSize);
    cudaMemcpy(factorBase_d, factorBase, size_BigInt*factorBaseSize, cudaMemcpyHostToDevice);
    
    

    gettimeofday(&stop_step3, NULL);
    elapsedTime_step3 += (stop_step3.tv_sec - start_step3.tv_sec) * 1000.0;               // sec to ms
    elapsedTime_step3 += (stop_step3.tv_usec - start_step3.tv_usec) / 1000.0;

    // STEP4 : CREO LISTAYPSILONNUBMER E LI SETACCIO
    // Facendo assumere ad r valori interi successivi a sqrt(n), si trovano almeno t+1 valori y=r^2-n che abbiano tutti i loro fattori primi in B.
    cout << "step 4"<< endl;
    gettimeofday(&start_step4, NULL);
    
//////////// CUDA start/////////////////////////////////////////

    BigInt *listYpsilonNumber, *listYpsilonNumber_test,*listYpsilonNumber_d, *nRSA_punt ,*nRSA_d, *x_0_punt, *x_0_d;
    long *m_punt, *m_d;
    
    listYpsilonNumber = (BigInt*)malloc(size_BigInt*(2*m)); // alloco su host
    listYpsilonNumber_test = (BigInt*)malloc(size_BigInt*(2*m)); // alloco su host
    nRSA_punt = (BigInt*)malloc(size_BigInt);
    x_0_punt = (BigInt*)malloc(size_BigInt);
    m_punt = (long*)malloc(size_long);
    
    cudaMalloc((void**)&listYpsilonNumber_d, size_BigInt*(2*m)); // alloco su device
    cudaMalloc((void**)&nRSA_d, size_BigInt);
    cudaMalloc((void**)&x_0_d, size_BigInt);
    cudaMalloc((void**)&m_d, size_long);
    
    *nRSA_punt = nRSA;
    *x_0_punt = x_0;
    *m_punt = m;
    
    cudaMemcpy(nRSA_d, nRSA_punt, size_BigInt, cudaMemcpyHostToDevice);
    cudaMemcpy(x_0_d, x_0_punt, size_BigInt, cudaMemcpyHostToDevice);
    cudaMemcpy(m_d, m_punt, size_long, cudaMemcpyHostToDevice);
    
    cout << "listyspilon number n thread = " << 2*m << endl;
    //int dim_a = 2*m/512 + 1;  
    //dim3 DimGrid(dim_a,1); dim3 DimBlock(512,1,1);  
  
    dim3 DimGrid(2*m,1); dim3 DimBlock(1,1,1); 
    
    genlistypsilon<<<DimGrid,DimBlock>>>(listYpsilonNumber_d,nRSA_d,x_0_d,m_d);
   
    cudaMemcpy(listYpsilonNumber, listYpsilonNumber_d, size_BigInt*(2*m), cudaMemcpyDeviceToHost);
    cudaMemcpy(listYpsilonNumber_test, listYpsilonNumber_d, size_BigInt*(2*m), cudaMemcpyDeviceToHost);
    
    cudaThreadSynchronize();

    
//////////// CUDA finish ///////////////////////////////////////////////



    int *start;
    start = (int*)malloc(sizeof(int));
    *start = 0;
    
    
    if (factorBase[0] == 2){
	cout << "caso2" << endl;
	
	//////////// CUDA start/////////////////////////////////////////
	
	//int dim_a = 2*m/512 + 1;  
	
	//dim3 DimGrid(dim_a,1); dim3 DimBlock(512,1,1);  
	
	dim3 DimGrid(2*m,1); dim3 DimBlock(1,1,1); 
	
	cout << "div2whileposs number n thread = " << 2*m << endl;
	
	div2whileposs<<<DimGrid,DimBlock>>>(listYpsilonNumber_d);
	
        cudaMemcpy(listYpsilonNumber_test, listYpsilonNumber_d, size_BigInt*(2*m), cudaMemcpyDeviceToHost);
	
	//////////// CUDA finish////////////////////////////////////////
	
        *start = 1;
    }
    ///////////////// fin qui ok ////////////////////////
    cudaThreadSynchronize();

   
    //////////// CUDA start/////////////////////////////////////////
    
    Lock lock;
    
    int *start_d;
    cudaMalloc((void**)&start_d, sizeof(int));
    cudaMemcpy(start_d, start, sizeof(int), cudaMemcpyHostToDevice);
    
    cudaMemcpy(listYpsilonNumber_d,listYpsilonNumber_test,  size_BigInt*(2*m), cudaMemcpyHostToDevice); 
    
    
    cout << "div3whileposs number n thread = " << factorBaseSize << endl;
    
    if (*start == 0){
      dim3 DimGrid2(factorBaseSize-1,1); dim3 DimBlock2(1,1,1);
      div3up_while_poss<<<DimGrid2,DimBlock2>>>(lock,factorBase_d,listYpsilonNumber_d,nRSA_d,start_d,x_0_d,m_d);
      //div3up_while_poss<<<DimGrid2,DimBlock2>>>(factorBase_d,listYpsilonNumber_d,nRSA_d,start_d,x_0_d,m_d);
    }else{
      dim3 DimGrid2(factorBaseSize-1,1); dim3 DimBlock2(1,1,1);
      //div3up_while_poss<<<DimGrid2,DimBlock2>>>(factorBase_d,listYpsilonNumber_d,nRSA_d,start_d,x_0_d,m_d);
      div3up_while_poss<<<DimGrid2,DimBlock2>>>(lock,factorBase_d,listYpsilonNumber_d,nRSA_d,start_d,x_0_d,m_d);
    }  
   
    
   
    
   //cudaMemcpy(listYpsilonNumber, listYpsilonNumber_d, size_BigInt*(2*m), cudaMemcpyDeviceToHost);
   cudaMemcpy(listYpsilonNumber_test, listYpsilonNumber_d, size_BigInt*(2*m), cudaMemcpyDeviceToHost);

///////////////////////// CUDA finish ////////////////////////////////////////////   
    cudaThreadSynchronize();
/*
    cout << "//////////dopo step 3//////////////////////////////////////////////////" << endl;
    for (long j = -m ; j < m ; j++){
        
            cout << "a(" << j << ") = " << listYpsilonNumber_test[j+m] << endl;
    }
  
*/
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
            //#pragma omp atomic
            contRelationOk++; // mi serve per sapere la dimensione dell array ed essere comodo dopo..
        //cout <<"ok = " << j+m << endl;
	  
	}
    }

    
    // controllo che ci siano piu relazioni che numeri nella factorbase
    if (factorBaseSize+1 >  contRelationOk){
          //  cout << endl;
            cout << "ERROR! WE NEED MORE COLUMNS THAN ROWS" << endl;
        return -1;
    }
   
    //cout << "m = "<< m << endl;
    
   
    //BigInt listYpsilonNumberOk[contRelationOk];
    int listExponentOk_mod2[contRelationOk][factorBaseSize+1];
    BigInt listExponentOk[contRelationOk][factorBaseSize+1];

    //#pragma omp parallel for
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
    std::vector<BigInt> jBuoni ;

    for (long j = -m ; j < +m ; j++){
        if (listYpsilonNumber_test[j+m] == -1 or listYpsilonNumber_test[j+m] == 1){
            if (listYpsilonNumber_test[j+m] == -1){
                listExponentOk_mod2[jBuonicont][0] = 1 ;// numero negativo
            }else{
                listExponentOk_mod2[jBuonicont][0] = 0 ;// numero positivo
            }
                //listYpsilonNumberOk[jBuonicont] = listYpsilonNumber[j+m];

           for (long p = 0; p < factorBaseSize; p++){

                BigInt fact_temp = 1;
                bool ok = true;
                do{
                          //  cout << "5c" << endl;
                    BigInt pmpz = factorBase[p];
		   // cout << pmpz <<", " ;
		   //cout << " lynumber" <<listYpsilonNumber[j+m] ;
		   //cout << "  || " << endl;
		    
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
    cout << endl;
    
    BigInt *jBuoni_arr;
    jBuoni_arr = (BigInt*)malloc(sizeof(BigInt)*jBuoni.size());
    
    
    for (int i = 0; i < jBuoni.size(); i++){
      
      //cout << "jBuoni ("<<i<<") = "<< jBuoni.at(i)<< endl; 
      *(jBuoni_arr+i) = jBuoni.at(i);
    }
    
    
       // cout << endl;
    gettimeofday(&stop_step5, NULL);
    elapsedTime_step5 += (stop_step5.tv_sec - start_step5.tv_sec) * 1000.0;               // sec to ms
    elapsedTime_step5 += (stop_step5.tv_usec - start_step5.tv_usec) / 1000.0;


     // DA PARALLELIZZARE
    // STEP 6 METODO DI GAUSS
    // Con il metodo di eliminazione di Gauss si determinano alcuni dei vettori v_2(y_i) che danno somma uguale al vettore nullo
     gettimeofday(&start_step6, NULL);
  /*  
    for (int i = 0; i <contRelationOk ; i++ ){
      
	for (int j = 0; j <factorBaseSize+1 ; j++ ){
      
	  cout << listExponentOk_mod2[i][j] << ", " ;
	}      
	cout << endl;
    }
    */
 
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
            if ( cont_rig == row_i){

                cont_rig++;
            }else{
                // devo scambiarle

                //#pragma omp parallel for
                for (int z = 0; z < contRelationOk; z++){
                    //cout << omp_get_thread_num() << endl;
                    BigInt temp = listExponentOk_mod2[z][row_i];
                    listExponentOk_mod2[z][row_i] = listExponentOk_mod2[z][cont_rig];
                    listExponentOk_mod2[z][cont_rig] = temp;
                }
                cont_rig++;
            }

            // cerco le righe che hanno 1 nella colonna di adesso e le sommo alla riga selezionata
            //#pragma omp parallel for schedule(dynamic)
            for(int k = cont_rig; k <  factorBaseSize+1; k++){
                if (listExponentOk_mod2[j][k] == 1){
                    for (int cc = 0; cc < contRelationOk; cc++){

                        listExponentOk_mod2[cc][k] = (listExponentOk_mod2[cc][k] + listExponentOk_mod2[cc][cont_rig-1]) % 2 ;
                    }

                }
            }
         }
    }
  
  /*for (int i = 0; i <contRelationOk ; i++ ){
      
	for (int j = 0; j <factorBaseSize+1 ; j++ ){
      
	  cout << listExponentOk_mod2[i][j] << ", " ;
	}      
	cout << endl;
    }
   */ 
    

    // TROVARE LA SOLUZIONE ( vettore spazio nullo) valutando le x dal basso!!!!
    int vect_solution[contRelationOk]; // 0 var libere, -1 da valutare
    int vect_pivot[contRelationOk];  // -1 no pivot, n pos riga

    //#pragma omp parallel for
    for (long i = 0; i< contRelationOk; i++){
        vect_solution[i] = 0;
        vect_pivot[i] = -1;
    }

    long start_to_check_pivot = 0;
    long n_free_var = contRelationOk;

    for (long i = 0; i< contRelationOk; i++){
        for (long j = start_to_check_pivot; j < factorBaseSize+1; j++){
            if (listExponentOk_mod2[i][j] == 1){

                vect_solution[i] = -1;
                vect_pivot[i] = j;
                start_to_check_pivot++;
                n_free_var--;
                break;
            }
        }
    }
    
    BigInt *result, *result_d, *jBuoni_d;
    long *contRelationOk_d, *contRelationOk_pu, *factorBaseSize_d, *factorBaseSize_pu, *jBuonicont_d, *jBuonicont_pu;
    
    int *vect_solution_d, *vect_pivot_d;
    
    result = (BigInt*)malloc(sizeof(BigInt));
    *result = 0;
    
    contRelationOk_pu = (long*)malloc(sizeof(long));
    *contRelationOk_pu = contRelationOk;
    
    factorBaseSize_pu = (long*)malloc(sizeof(long));
    *factorBaseSize_pu = factorBaseSize;
    
    jBuonicont_pu = (long*)malloc(sizeof(long));
    *jBuonicont_pu = jBuonicont;
    
    cudaMalloc((void**)&result_d, size_BigInt);
    cudaMalloc((void**)&contRelationOk_d, size_long);
    cudaMalloc((void**)&factorBaseSize_d, size_long);
    cudaMalloc((void**)&jBuonicont_d, size_long);
    cudaMalloc((void**)&vect_solution_d,sizeof(int)*contRelationOk);
    cudaMalloc((void**)&vect_pivot_d,sizeof(int)*contRelationOk);
    cudaMalloc((void**)&jBuoni_d,jBuonicont*sizeof(BigInt));
    
    // BigInt listExponentOk[contRelationOk][factorBaseSize+1];
    // BigInt listExponentOk[n][m];
    BigInt* listExponentOk_d;
    int* listExponentOk_mod2_d;
    size_t pitch, pitch2;
    
   cudaMallocPitch(&listExponentOk_d, &pitch, size_BigInt*(factorBaseSize+1),contRelationOk);
   cudaMemcpy2D(listExponentOk_d,pitch,listExponentOk,(factorBaseSize+1)*sizeof(BigInt),(factorBaseSize+1)*sizeof(BigInt),contRelationOk, cudaMemcpyHostToDevice);
  
   cudaMallocPitch(&listExponentOk_mod2_d, &pitch2, sizeof(int)*(factorBaseSize+1) , contRelationOk); 
   cudaMemcpy2D(listExponentOk_mod2_d,pitch2,listExponentOk_mod2,(factorBaseSize+1)*sizeof(int),(factorBaseSize+1)*sizeof(int),contRelationOk, cudaMemcpyHostToDevice);
   
    cudaMemcpy(result_d, result , size_BigInt, cudaMemcpyHostToDevice);
    cudaMemcpy(contRelationOk_d , contRelationOk_pu , size_long, cudaMemcpyHostToDevice);
    cudaMemcpy(factorBaseSize_d , factorBaseSize_pu , size_long, cudaMemcpyHostToDevice);
    cudaMemcpy(jBuonicont_d , jBuonicont_pu , size_long, cudaMemcpyHostToDevice);
    cudaMemcpy(vect_solution_d , vect_solution , sizeof(int)*contRelationOk, cudaMemcpyHostToDevice);
    cudaMemcpy(vect_pivot_d , vect_pivot , sizeof(int)*contRelationOk, cudaMemcpyHostToDevice);
    cudaMemcpy(jBuoni_d , jBuoni_arr , jBuonicont*sizeof(BigInt), cudaMemcpyHostToDevice);
    
    //cout << "evaluateexpression number n thread = " << n_free_var << endl;
    
    dim3 DimGrid2(n_free_var,1); dim3 DimBlock2(1,1,1);
    
       
    evaluate_solutions<<<DimGrid2,DimBlock2>>>( nRSA_d,factorBase_d,x_0_d ,jBuoni_d ,listExponentOk_d, pitch ,listExponentOk_mod2_d, pitch2, vect_solution_d, vect_pivot_d, result_d, contRelationOk_d, factorBaseSize_d, jBuonicont_d);
  
    cudaThreadSynchronize();
    
    int yo[contRelationOk][factorBaseSize+1];
    
    cudaMemcpy2D(yo,(factorBaseSize+1)*sizeof(int),listExponentOk_mod2_d,pitch2,(factorBaseSize+1)*sizeof(int),contRelationOk, cudaMemcpyDeviceToHost);

    cudaMemcpy(result , result_d , size_BigInt, cudaMemcpyDeviceToHost);
   
    cudaDeviceReset();
    
    cout << "result = " << *result <<  endl;
    
   // exit(1);
    
    
//////////////////////////// cuda finish ////////////////////////////////////////
    
    
    if (*result != 0 ){
        
        return *result;
    }

    gettimeofday(&stop_step6, NULL);
    elapsedTime_step6 += (stop_step6.tv_sec - start_step6.tv_sec) * 1000.0;               // sec to ms
    elapsedTime_step6 += (stop_step6.tv_usec - start_step6.tv_usec) / 1000.0;

    return -2;

}
int main(int argc, char *argv[])
{
timeval start_tot, stop_tot;
double elapsedTime_tot;
gettimeofday(&start_tot, NULL);
    // STEP1 : Viene dato in input il numero naturale dispari n>1.
    //long nRSA =  argc[1];
    




      BigInt nRSA, result, x_0;

    if ( argc > 1 ) {
	//nRSA = 1100017;
        //nRSA = 746294513;
	nRSA = atoll(argv[1]);
	//nRSA = 492700411219; // facile 3 cifre
	//nRSA = 6179743647887;
	//nRSA = 581677579783;
        //long k = pow(exp(sqrt(log(nRSA.get_ui())*log(log(nRSA.get_ui())))),0.35);
        //long m = k * k*k;

        //mpz_sqrt(x_0.get_mpz_t(), nRSA.get_mpz_t());

        x_0 = sqrt(nRSA);

        long k = 10;
        long m = 100;

	//long k = 790;
	//long m = 1920;
	
        bool soluzione_trovata = false;

        int cont_meno1 = 0;
        int cont_meno2 = 0;
        while(!soluzione_trovata){
	
         //   cout << "lancio con nRSA = "<< nRSA << ", k = "<< k << ", m = " << m << endl;
	 
            result = factorize(nRSA,x_0,k,m);
            if ( result == -2){
                k = k+10;
                m=  m+30;
                //k = k+20;
                //m = k*k*k;


                cont_meno2++;
               // cout << "cont mano 1 = " << cont_meno1 << endl;
               // cout << "cont mano 2 = " << cont_meno2 << endl;
                cout << "k =" << k << ", m = " << m << endl;
               // cout << endl;
            }
            else if ( result == -1){
                //k = k +5;
                //m = k*k*k;
                k=k+10;
                m=m+10;
                cont_meno1++;
            //    cout << "k = " << k << ", m = " << m<< endl;
              //  cout << "cont mano 1 = " << cont_meno1 << endl;
              //  cout << "cont mano 2 = " << cont_meno2 << endl;
                cout << "k =" << k << ", m = " << m << endl;
              //  cout << endl;
            }else{

                soluzione_trovata = true;
            }

        }
        gettimeofday(&stop_tot, NULL);
        elapsedTime_tot = (stop_tot.tv_sec - start_tot.tv_sec) * 1000.0;               // sec to ms
        elapsedTime_tot += (stop_tot.tv_usec - start_tot.tv_usec) / 1000.0;            // us to ms


        //cout << "k = " << k << ", m = " << m << endl;
        cout << "-1 volte = " << cont_meno1 << ", -2 volte = " << cont_meno2 << endl;
        cout << endl << "nRSA = " << result << " x " << nRSA/result <<endl;
        cout << "nrsa = " << nRSA << endl;
        cout << " t-totale = " << elapsedTime_tot << " ms.\n" << endl;
        cout << " t-step3 = " << elapsedTime_step3 << " ms.\n" << endl;
        cout << " t-step4 = " << elapsedTime_step4 << " ms.\n" << endl;
        cout << " t-step5 = " << elapsedTime_step5<< " ms.\n" << endl;
        cout << " t-step6 = " << elapsedTime_step6 << " ms.\n" << endl;


    }else{

        cout << "nessun argomento immesso!!!" << endl;
    }
    return 0;
}

