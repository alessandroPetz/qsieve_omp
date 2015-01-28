#include <math.h>       /* sqrt */
#include <vector>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <list>
#include <sys/time.h>
#include <bitset>
#include <omp.h>

#define BigInt long long


using namespace std;

timeval start_step3, stop_step3,start_step4, stop_step4,start_step5, stop_step5, start_step6, stop_step6;
double elapsedTime_step3 = 0,elapsedTime_step4 = 0,elapsedTime_step5 = 0,elapsedTime_step6 = 0;

//FacorBase dichiarata qui cosi non bisogna calcolarsi ogni volta tutti i primi.
std::vector<BigInt> factorBase ;
long long primo = 0;


//Fonte = http://en.wikipedia.org/wiki/Binary_GCD_algorithm
BigInt gcd(BigInt u, BigInt v)
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

//calculates (a * b) % c taking into account that a * b might overflow
//Fonte = http://stackoverflow.com/questions/20971888/modular-multiplication-of-large-numbers-in-c
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


//modular exponentiation
// Info : http://en.wikipedia.org/wiki/Modular_exponentiation
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

// modular inversion
// Fonte : http://rosettacode.org/wiki/Modular_inverse#C.2B.2B
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

//test Miller-Rabin
// Info : http://it.wikipedia.org/wiki/Test_di_Miller-Rabin
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

// Calcola il successivo secondo il lemma di hensell
// Info : http://en.wikipedia.org/wiki/Hensel%27s_lemma
BigInt  hensellemma(BigInt r1, BigInt  n, BigInt p)
{

    if (r1 == 0){ return 0; }

    //return (r1 - ( ( (r1*r1-n) % p ) * mul_inv( r1*2, p )  % p )) % p;

    BigInt temp;
    //BigInt temp2 = r1*2;
    //BigInt expneg = -1;

    //mpz_powm (temp.get_mpz_t(), temp2.get_mpz_t() , expneg.get_mpz_t(), p.get_mpz_t() );
    //cout << "inv( " << r1*2 << ", "<< p << ") = " endl;
    temp = mul_inv( r1*2, p );
   //r2 = r1 - (r1^2 - n )/ (2r1);
    //temp = modulo(r1*2, -1, p);
    //cout << "temp = " << temp << endl;
    //return temp;
    //return (r1 - ( ( (r1*r1) - n) % p ) * temp ) % p;
    temp = (r1 - ( ( (r1*r1) - n) % p ) * temp ) % p;

    return temp;
}

// effettua il test di Legendre
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

// Calcola il numero di Legendre
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

// Calcola il residuo QUadratico
// Info: http://it.wikipedia.org/wiki/Residuo_quadratico
//Fonte: Tradotta in c++ dalla versione mpz -> https://gmplib.org/list-archives/gmp-devel/2006-May/000633.html
int mpz_sqrtm(BigInt n, BigInt p) {
//cout << " mpz 11 " << endl;

    BigInt q_int;
    BigInt w, n_inv, y;
    int i, s;

    if ( n % p == 0){
        q_int = 0;
        return q_int;
    }

    std::bitset<(sizeof(p)*8)-1> foo2 (p);
    if(foo2.test(1)) {
        //cout << "uscita checkbit" << endl;
        q_int = modulo (n, (p+1)/4 , p  ) ;
        return q_int;
    }


    q_int = p - 1;
    s = 0; /* Factor out 2^s from q */

    std::bitset<(sizeof(q_int)*8)-1> foo (q_int);
    while (!foo.test(s)) {
        s++;
    }

    q_int = q_int / pow(2,s);

    w = 2;

    while (legendre(w,p) != -1){
            w++;
    }

    w = modulo(w,q_int,p);
    q_int++;
    q_int = q_int / 2;
    q_int = modulo(n , q_int , p);
    n_inv = mul_inv(n,p);


    for(;;) {

        y = modulo(q_int,2,p);
        y = mulmod(y,n_inv,p);

        i = 0;

        while (y != 1){
           // cout << "aaaaaaaaaaaa, y = " << y << "^2 mod"<< p << endl;
            i++;
            y = modulo(y,2,p);
        }

        if (i == 0 ) {

            return q_int;
        }

        if (s-i == 1){
            q_int = q_int * w;
        } else {
            y = modulo(w,1 << (s-i-1), p);
            q_int = q_int * y;
        }

        //mpz_mod(q, q, p); /* r = r * w^(2^(s-i-1)) (mod p) */
        q_int = q_int % p;
    }

    //mpz_clear(w); mpz_clear(n_inv); mpz_clear(y);
    //return 0;
    return q_int;
}

// Funzione che cerca di fattorizzare un numero RSA applicando
// l'algoritmo del Quadratic Sieve
// nRSA = numero RSA da fattorizzare
// x_0 = sqrt(RSA), numero da usare per il calcolo delle relazioni
// k = Numero max dentro factorBase
// m = intervallo min e max delle relazioni da considerare

/*puo restituire        -2 = trovati solo fattori banali -> aumentare k
                        -1 = meno relazione che fattori -> aumentare b
                        - soluzione!
*/

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

    BigInt primo_temp = 0;
    if (factorBase.size() != 0){
        primo_temp =  factorBase.at(factorBase.size()-1);
    }
    // aggiorno la factorbase
    while(true){

        while(true){

            primo_temp++;
            if ( Miller(primo_temp, 20 ) ) break;

        }

        if (primo_temp > k) {
                break;
        }else if (testLegendre(primo_temp,nRSA)){
            primo = primo_temp;
            factorBase.push_back(primo);
        }else{
            primo = primo_temp;
        }
    }

    gettimeofday(&stop_step3, NULL);
    elapsedTime_step3 += (stop_step3.tv_sec - start_step3.tv_sec) * 1000.0;               // sec to ms
    elapsedTime_step3 += (stop_step3.tv_usec - start_step3.tv_usec) / 1000.0;

    // STEP4 : CREO LISTAYPSILONNUBMER E LI SETACCIO
    // Facendo assumere ad r valori interi successivi a sqrt(n), si trovano almeno t+1 valori y=r^2-n che abbiano tutti i loro fattori primi in B.
    cout << "step 4"<< endl;
    gettimeofday(&start_step4, NULL);

    BigInt listYpsilonNumber[2*m];
    BigInt listYpsilonNumber_test[2*m];

    //cout << "step 4b parallel" << endl;
    //cout << "LISTA RELAZIONI :" << endl;
    #pragma omp parallel for
    for (long j = -m ; j < m ; j++){

        BigInt temp = (pow(x_0 + j, 2)) - nRSA;
        listYpsilonNumber[j+m] = temp;
        listYpsilonNumber_test[j+m] = temp;

    }
    //cout << "m = "<< m << endl;

    // FASE SETACCIO, CONTROLLARE QUALI DEI NUMERI GENERATI NELLA YLIST SONO FATTORIZZABILI
    // CON LA FACTORBASE
    // tirato fuori il caso 2 per vedere se diminuisce il tempo.. pare di si.
    int start = 0;
    if (factorBase.at(0) == 2){
    cout << "caso2" << endl;
	#pragma omp parallel for schedule(dynamic)
        for (long j = -m ; j < +m ; j++){

            while( listYpsilonNumber_test[j+m] % 2 == 0){

                listYpsilonNumber_test[j+m] /= 2;
            }
        }
        start = 1;
    }


    //cout << "caso3+" << endl;
    #pragma omp parallel for schedule(dynamic)
    for (long p = start; p< factorBase.size(); p++){

        //cont_factor++;
        BigInt l_fact = factorBase.at(p);


        BigInt p_origin = factorBase.at(p);
        // finche l_fact < radice( numero piu grande da checkare  )
        long exp = 1;

        //cout << "yo"<< endl;
        BigInt s,h, s_calc, h_calc,s_calc_neg, h_calc_neg, s_calc_temp, h_calc_temp ;

        bool ok_sieve = true;
        while( ok_sieve ){
            //cout << "yo" << endl;
            //mpz_pow_ui(l_fact.get_mpz_t(),p_origin.get_mpz_t(), exp);	// primo elevato a numero

            l_fact = pow(p_origin, exp);

            if (l_fact > 2*m ) ok_sieve = false;
            // tutte le j = s(h) - X_0 + k*l_fact per k < l_fact
            if (exp == 1){						// se il primo è p
//cout << "yo3" << endl;
            //   cout << "caso exp = 1"<< endl;
	        //cout << "nRSA= " << nRSA << "p_origin =" << p_origin << endl;
                s = mpz_sqrtm(nRSA,p_origin);
		//cout <<"s = "<< s << endl;

                h = l_fact - s;

                s_calc_temp = s - (x_0 % l_fact);
                if (s_calc_temp < 0) {
                        s_calc_temp = l_fact + s_calc_temp;
                }

                h_calc_temp = (h - (x_0 % l_fact)) % l_fact;
                if (h_calc_temp < 0) {
                        h_calc_temp = l_fact + h_calc_temp;
                }

            }else{
//cout << "yo5" << endl;
             // cout << "caso exp = 2"<< endl;							// se abbiamo un p^k
            //  cout << "prima..s = "<< s << " nRSA ="<< nRSA <<", l_fact = "<< l_fact << endl;
                s= hensellemma(s,nRSA, l_fact);				// ho r2 e r2a
              //cout << "ris = "<< s << endl;
//		cout << "yo 51" << endl;
                h= l_fact - s;

                s_calc_temp = (s - (x_0 % l_fact)) % l_fact;
                if (s_calc_temp < 0){
                    s_calc_temp = l_fact + s_calc_temp;
                }
                h_calc_temp = (h - (x_0 % l_fact)) % l_fact;

                if (h_calc_temp < 0) {
                        h_calc_temp = l_fact + h_calc_temp;
                }


            }
            //cout << "yo6" << endl;
            BigInt k = 0;
            while(true){

              //  cout << "whileeeeeee con k ="<< k <<endl;
                s_calc = s_calc_temp  + (l_fact * k );
                h_calc = h_calc_temp  + (l_fact * k);
                s_calc_neg = s_calc_temp + (l_fact * (-k));
                h_calc_neg = h_calc_temp + (l_fact * (-k));

                BigInt a = s_calc + m, b = h_calc + m, c = s_calc_neg + m, d = h_calc_neg + m;
               // cout << "a= "<< a << ", b = " << b <<" c= "<< c << ", d = " << d << endl;

                if ((a < 0 or a >= 2*m) and (b<0 or b >= 2 * m) and (c < 0 or c >= 2*m) and (d < 0 or d >= 2*m)){
                 //   cout <<"exitt whileeeeeeeee"<<endl;
                    break;
                }
                // applica divisione j positivo
                if (a >=0 and a < 2*m ){
                    #pragma omp critical
                    listYpsilonNumber_test[a] = listYpsilonNumber_test[a] / factorBase.at(p);

                }
                if (a >= 0 and b < 2*m and b != a){
                    #pragma omp critical
		    listYpsilonNumber_test[b] = listYpsilonNumber_test[b] / factorBase.at(p);
                }
                if (c >= 0 and c < 2*m and c != b and c !=a){
                    #pragma omp critical
                    listYpsilonNumber_test[c] = listYpsilonNumber_test[c] / factorBase.at(p);

                }
                if (d >= 0 and d< 2*m and d !=c and d!= b and d != a){
                    #pragma omp critical
                    listYpsilonNumber_test[d] = listYpsilonNumber_test[d] / factorBase.at(p);
                }

                k++;
            }
            exp++;                                                     // aumentiamo l esponente
        }
    }

    gettimeofday(&stop_step4, NULL);
    elapsedTime_step4 += (stop_step4.tv_sec - start_step4.tv_sec) * 1000.0;               // sec to ms
    elapsedTime_step4 += (stop_step4.tv_usec - start_step4.tv_usec) / 1000.0;

    // STEP 5 modulo A DUE DEGLI ESPONENTI!
    // Per ognuno dei valori y_1, y_2,...,y_{t+1} si calcola il vettore in Z2: di v_2(y_i)=(e_1,e_2,...,e_t)
    // dove e_i è la riduzione modulo 2 dell'esponente di p_i nella fattorizzazione di y_i.
    gettimeofday(&start_step5, NULL);

    long contRelationOk = 0;
    #pragma omp parallel for
    for (long j = -m ; j < +m ; j++){
        if (listYpsilonNumber_test[j+m] == -1 or listYpsilonNumber_test[j+m] == 1){
            #pragma omp atomic
            contRelationOk++; // mi serve per sapere la dimensione dell array ed essere comodo dopo..
        //cout <<"ok = " << j+m << endl;

	}
    }

    // controllo che ci siano piu relazioni che numeri nella factorbase
    if (factorBase.size()+1 >  contRelationOk){
          //  cout << endl;
            cout << "ERROR! WE NEED MORE COLUMNS THAN ROWS" << endl;
        return -1;
    }



    BigInt listYpsilonNumberOk[contRelationOk];
    BigInt listExponentOk_mod2[contRelationOk][factorBase.size()+1];
    BigInt listExponentOk[contRelationOk][factorBase.size()+1];

    #pragma omp parallel for
    for (long x = 0; x <contRelationOk; x++){
        for (long y = 0; y <factorBase.size()+1; y++){
            listExponentOk_mod2[x][y] = 0;
            listExponentOk[x][y] = 0;

        }
    }
    // metto il segno nella prima linea, recuper gli esponenti, e mi segno i jnumber ok!!!!
    // creo un altra matrice con i soli valori che mi interessano, cioè coi jnumber che sono fattorizzabili con
    // la factor base, e già che ci sono gli faccio il modulo 2.

    // cout << "LISTA RELAZIONI ACCETTATE :" << endl;

    long cont = 0;
    std::vector<BigInt> jBuoni ;

    for (long j = -m ; j < +m ; j++){
        if (listYpsilonNumber_test[j+m] == -1 or listYpsilonNumber_test[j+m] == 1){
            if (listYpsilonNumber_test[j+m] == -1){
                listExponentOk_mod2[cont][0] = 1 ;// numero negativo
            }else{
                listExponentOk_mod2[cont][0] = 0 ;// numero positivo
            }
                listYpsilonNumberOk[cont] = listYpsilonNumber[j+m];
	  #pragma omp parallel for schedule(dynamic)
           for (long p = 0; p < factorBase.size(); p++){

                BigInt fact_temp = 1;
                bool ok = true;
                do{
                          //  cout << "5c" << endl;
                    BigInt pmpz = factorBase.at(p);
		   // cout << pmpz <<", " ;
		   // cout << " lynumber" <<listYpsilonNumber[j+m] ;
		   // cout << "  || " << endl;
                    fact_temp = fact_temp * pmpz;
                    if( (listYpsilonNumber[j+m] % fact_temp )== 0){
                      #pragma omp critical
		      listExponentOk[cont][p+1]++;
			//cout<< "yo ";
                    }else{
                        ok = false;
                    }
                }while (ok);
                listExponentOk_mod2[cont][p+1] = listExponentOk[cont][p+1];
                listExponentOk_mod2[cont][p+1] = listExponentOk_mod2[cont][p+1] % 2;

                    //cont_factor++;
            }

            cont++;
            jBuoni.push_back(j);
        }
    }
    cout << endl;


       // cout << endl;
    gettimeofday(&stop_step5, NULL);
    elapsedTime_step5 += (stop_step5.tv_sec - start_step5.tv_sec) * 1000.0;               // sec to ms
    elapsedTime_step5 += (stop_step5.tv_usec - start_step5.tv_usec) / 1000.0;


     // DA PARALLELIZZARE
    // STEP 6 METODO DI GAUSS
    // Con il metodo di eliminazione di Gauss si determinano alcuni dei vettori v_2(y_i) che danno somma uguale al vettore nullo

   cout << "step 6" << endl;
    gettimeofday(&start_step6, NULL);
//cout << "mod2 primaaaaa" << endl;
    gettimeofday(&start_step6, NULL);

    long cont_rig = 0;
    for(long j = 0; j < contRelationOk; j++){ // colonne
        long row_i = -1;  // row i = j

        for(long i = cont_rig; i< factorBase.size()+1; i++ ){ // righe, cerco Aij = 1
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

                #pragma omp parallel for
                for (int z = 0; z < contRelationOk; z++){
                    //cout << omp_get_thread_num() << endl;
                    BigInt temp = listExponentOk_mod2[z][row_i];
                    listExponentOk_mod2[z][row_i] = listExponentOk_mod2[z][cont_rig];
                    listExponentOk_mod2[z][cont_rig] = temp;
                }
                cont_rig++;
            }

            // cerco le righe che hanno 1 nella colonna di adesso e le sommo alla riga selezionata
            #pragma omp parallel for schedule(dynamic)
            for(int k = cont_rig; k <  factorBase.size()+1; k++){
                if (listExponentOk_mod2[j][k] == 1){
                    for (int cc = 0; cc< contRelationOk; cc++){

                        listExponentOk_mod2[cc][k] = (listExponentOk_mod2[cc][k] + listExponentOk_mod2[cc][cont_rig-1]) % 2 ;
                    }

                }
            }
         }
    }



    // TROVARE LA SOLUZIONE ( vettore spazio nullo) valutando le x dal basso!!!!
    int vect_solution[contRelationOk]; // 0 var libere, -1 da valutare
    int vect_pivot[contRelationOk];  // -1 no pivot, n pos riga

    #pragma omp parallel for
    for (long i = 0; i< contRelationOk; i++){
        vect_solution[i] = 0;
        vect_pivot[i] = -1;
    }

    long start_to_check_pivot = 0;
    long n_free_var = contRelationOk;

    for (long i = 0; i< contRelationOk; i++){
        for (long j = start_to_check_pivot; j < factorBase.size()+1; j++){
            if (listExponentOk_mod2[i][j] == 1){

                vect_solution[i] = -1;
                vect_pivot[i] = j;
                start_to_check_pivot++;
                n_free_var--;
                break;
            }
        }
    }



    BigInt result = 0;
    // PARALLELIZZ0 IN MODO CHE TUTTI GLI ARRAI SOLUZIONI POSSIBILI
    // VEGONO VALUTATI IN PARALLELO
    cout << "step 7" << endl;
    #pragma omp parallel for schedule(dynamic)
    for(long s = 0; s < n_free_var; s++){
            //cout << "numeri free var" <<n_free_var << endl;
        // copia locale di vect solution... ???? è una soluzione buona? pare di si..
        int vect_solution_local[contRelationOk];
        for (long i = 0; i< contRelationOk; i++){
            vect_solution_local[i] = vect_solution[i];

        }

        long myvar = s;
        for (long i = 0; i< contRelationOk; i++){
            if (vect_pivot[i] == -1){

                if (myvar == 0){

                    vect_solution_local[i] = 1;
                    break;
                }else{

                    myvar--;
                }
            }
        }

        // qui valuto le x in base alla variabile libera assegnata!
        // da migliorare
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

        for (long v = 0; v < 2; v++){  // lo facciamo 2 volte, uno per il vettore, e l altro per il suo inverso

            if (v == 1){

                 for (long i = 0; i< contRelationOk; i++){
                    if (vect_solution_local[i] == 0){
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
            vector<BigInt>::iterator p2;
      //cout << "step 7" << endl;
            BigInt exponent[factorBase.size()];
            for (long i = 0; i < factorBase.size(); i++){
                exponent[i]=0;
            }
     //  cout << "step 7b" << endl;
            for (p2 = jBuoni.begin(); p2 != jBuoni.end(); p2++){

                if (vect_solution_local[cont] == 1){

                    a = (a * (x_0 + *p2)) % nRSA;
                    for (long i = 1; i < factorBase.size()+1; i++){

                        exponent[i-1] = exponent[i-1] + listExponentOk[cont][i];
                    }
                }
                cont++;
            }
            cont=0;
            for (long p = 0; p <  factorBase.size(); p++){

                if (exponent[cont]>0){

                   //  cout << "b = "<< b << "p = " << *p3 << "^" << exponent[cont] << endl;
                    BigInt p3mpz = factorBase.at(p);
                    BigInt temp = exponent[cont] * 2;
                    //mpz_powm(temp.get_mpz_t(),p3mpz.get_mpz_t(),temp.get_mpz_t(),nRSA.get_mpz_t());
                    temp = modulo(p3mpz, temp, nRSA);
                    temp = temp * b;
                    b = temp % nRSA;
                    //mpz_mod(b.get_mpz_t(),temp.get_mpz_t(),nRSA.get_mpz_t());

                    //b = (b * mod_pow(*p3,exponent[cont]/2,nRSA)) % nRSA;
                }
                cont++;
            }

            //STEP 8 = CALCOLO SOLUZIONE
            // Si calcola d=mcd(a-b,n) e se 1<d<n allora d è divisore non banale di n,
            // altrimenti si torna al passo 2) con una scelta di k più grande
         //   cout << "step 8" << endl;

            BigInt risultato_1, risultato_2, temp;

            temp = abs(a+b);
            risultato_1 = gcd(temp,nRSA);

	    // cout << "risultato 1 "<<risultato_1<< endl;

	    temp = abs(a-b);
        risultato_2 = gcd(temp,nRSA);

	    //cout << "risultato 2 "<<risultato_2<< endl;

	    //    cout << "step 8a" << endl;
           if ( risultato_1 != 1 and risultato_1 != nRSA ){

                 result = risultato_1;
            }
           if ( risultato_2 != 1 and risultato_2 != nRSA ){

                 result = risultato_2;
            }
        }
    }

    if (result != 0 ){
        // non si puù uscire subito una volta trovato il risultato??
        return result;
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
        nRSA = atoll(argv[1]);
         //long k = pow(exp(sqrt(log(nRSA.get_ui())*log(log(nRSA.get_ui())))),0.35);
        //long m = k * k*k;

        //mpz_sqrt(x_0.get_mpz_t(), nRSA.get_mpz_t());

        x_0 = sqrt(nRSA);

        long k = 10;
        long m = 100;

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


        cout << "k = " << k << ", m = " << m << endl;
        cout << "-1 volte = " << cont_meno1 << ", -2 volte = " << cont_meno2 << endl;
        cout << endl << "nRSA = " << result << " x " << nRSA/result <<endl;
        cout << "nrsa = " << nRSA << endl;
        cout << " t-totale = " << elapsedTime_tot << " ms.\n" << endl;
        cout << " t-step3 = " << elapsedTime_step3 << " ms.\n" << endl;
        cout << " t-step4 = " << elapsedTime_step4 << " ms.\n" << endl;
        cout << " t-step5 = " << elapsedTime_step5<< " ms.\n" << endl;
        cout << " t-step6 = " << elapsedTime_step6 << " ms.\n" << endl;
        cout << "[QSieve, ver. OMP senza GMP]" << endl;
        cout << endl;

    }else{

        cout << "nessun argomento immesso!!!" << endl;
    }
    return 0;
}

