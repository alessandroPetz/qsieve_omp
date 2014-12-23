#include <math.h>       /* sqrt */
#include <vector>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <list>
#include <sys/time.h>
#include <gmp.h>
#include <gmpxx.h>
#include <omp.h>

using namespace std;

timeval start_step3, stop_step3,start_step4, stop_step4,start_step5, stop_step5, start_step6, stop_step6;
double elapsedTime_step3 = 0,elapsedTime_step4 = 0,elapsedTime_step5 = 0,elapsedTime_step6 = 0;

//FacorBase dichiarata qui cosi non bisogna calcolarsi ogni volta tutti i primi.
std::vector<mpz_class> factorBase ;
mpz_class primo = 0;


int mpz_sqrtm(mpz_t q, const mpz_t n, const mpz_t p) {
    mpz_t w, n_inv, y;
    unsigned int i, s;
    //TMP_DECL;
    //TMP_MARK;
    if(mpz_divisible_p(n, p)) { /* Is n a multiple of p? */
        mpz_set_ui(q, 0); /* Yes, then the square root is 0. */
        return 1; /* (special case, not caught */
    } /* otherwise) */
    /*if(mpz_legendre(n, p) != 1) /* Not a quadratic residue? */
    /* return 0; /* No, so return error */
    if(mpz_tstbit(p, 1) == 1) { /* p = 3 (mod 4) ? */
        mpz_set(q, p);
        mpz_add_ui(q, q, 1);
        mpz_fdiv_q_2exp(q, q, 2);
        mpz_powm(q, n, q, p); /* q = n ^ ((p+1) / 4) (mod p) */
        return 1;
    }
    //MPZ_TMP_INIT(y, 2*SIZ(p));
    //MPZ_TMP_INIT(w, 2*SIZ(p));
    //MPZ_TMP_INIT(n_inv, 2*SIZ(p));
    mpz_init(y);
    mpz_init(w);
    mpz_init(n_inv);
    mpz_set(q, p);
    mpz_sub_ui(q, q, 1); /* q = p-1 */
    s = 0; /* Factor out 2^s from q */
    while(mpz_tstbit(q, s) == 0) s++;
    mpz_fdiv_q_2exp(q, q, s); /* q = q / 2^s */
    mpz_set_ui(w, 2); /* Search for a non-residue mod p */
    while(mpz_legendre(w, p) != -1) /* by picking the first w such that */
        mpz_add_ui(w, w, 1); /* (w/p) is -1 */
    mpz_powm(w, w, q, p); /* w = w^q (mod p) */
    mpz_add_ui(q, q, 1);
    mpz_fdiv_q_2exp(q, q, 1); /* q = (q+1) / 2 */
    mpz_powm(q, n, q, p); /* q = n^q (mod p) */
    mpz_invert(n_inv, n, p);
    for(;;) {
        mpz_powm_ui(y, q, 2, p); /* y = q^2 (mod p) */
        mpz_mul(y, y, n_inv);
        mpz_mod(y, y, p); /* y = y * n^-1 (mod p) */
        i = 0;
        while(mpz_cmp_ui(y, 1) != 0) {
            i++;
            mpz_powm_ui(y, y, 2, p); /* y = y ^ 2 (mod p) */
        }
        if(i == 0) { /* q^2 * n^-1 = 1 (mod p), return */
            //TMP_FREE;
            return 1;
        }
        if(s-i == 1) { /* In case the exponent to w is 1, */
            mpz_mul(q, q, w); /* Don't bother exponentiating */
        } else {
            mpz_powm_ui(y, w, 1 << (s-i-1), p);
            mpz_mul(q, q, y);
        }
        mpz_mod(q, q, p); /* r = r * w^(2^(s-i-1)) (mod p) */
    }
    mpz_clear(w); mpz_clear(n_inv); mpz_clear(y);
    return 0;
}

mpz_class  hensellemma(mpz_class r1, mpz_class  n, mpz_class p)
{

    if (r1 == 0){ return 0; }

    mpz_class temp;
    mpz_class temp2 = r1*2;
    mpz_class expneg = -1;

    mpz_powm (temp.get_mpz_t(), temp2.get_mpz_t() , expneg.get_mpz_t(), p.get_mpz_t() );

    return (r1 - ( ( (r1*r1) - n) % p ) * temp ) % p;
}


bool testLegendre(mpz_class number, mpz_class nRSA){

    mpz_class a;
    mpz_class b = (number-1)/2;
    mpz_class numbermpz = number;

    mpz_powm(a.get_mpz_t(),nRSA.get_mpz_t(),b.get_mpz_t(), number.get_mpz_t());

    if ( a != 1 ){
        return false;
    }
    return true;

}

mpz_class factorize(mpz_class nRSA, mpz_class x_0,  long k, long m){

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
    mpz_class primo_temp;
    // aggiorno la factorbase
    while(true){

        mpz_nextprime(primo_temp.get_mpz_t(), primo.get_mpz_t());
        if (primo_temp > k) {
                break;
        }else if (testLegendre(primo_temp,nRSA)){
            primo = primo_temp;
            factorBase.push_back(primo);
            //cout <<primo<< endl;
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

    mpz_class listYpsilonNumber[2*m];
    mpz_class listYpsilonNumber_test[2*m];

    cout << "step 4b parallel" << endl;
           //cout << "LISTA RELAZIONI :" << endl;
    #pragma omp parallel for
    for (long j = -m ; j < m ; j++){
        //cout << "a(" << j << ") = " << yNumb << endl;
        mpz_class temp = x_0 + j;
        mpz_pow_ui(temp.get_mpz_t(),temp.get_mpz_t(),2) ;
        mpz_sub(temp.get_mpz_t(),temp.get_mpz_t(),nRSA.get_mpz_t());
        listYpsilonNumber[j+m] = temp;
        listYpsilonNumber_test[j+m] = temp;
            //cout << "a(" << j << ") = " << listYpsilonNumber[j+m] << endl;
    }

    //cout << "m = "<< m << endl;

    // FASE SETACCIO, CONTROLLARE QUALI DEI NUMERI GENERATI NELLA YLIST SONO FATTORIZZABILI
    // CON LA FACTORBASE
    cout << "step 4c serie" << endl;
    mpz_class limiteFact;
    mpz_sqrt(limiteFact.get_mpz_t(),listYpsilonNumber[(2*m)-1].get_mpz_t());
    limiteFact = limiteFact;

    // tirato fuori il caso 2 per vedere se diminuisce il tempo.. pare di si.
    int start = 0;
    if (factorBase.at(0) == 2){
    cout << "caso2" << endl;

        #pragma omp parallel for schedule(dynamic)
        for (long j = -m ; j < +m ; j++){
            //cout <<"qui"<<endl;
            while( listYpsilonNumber_test[j+m] % 2 == 0){

                listYpsilonNumber_test[j+m] /= 2;
            }
        }
        start = 1;
    }


    cout << "caso2+" << endl;
    #pragma omp parallel for schedule(dynamic)
    for (long p = start; p< factorBase.size(); p++){

        //cont_factor++;
        mpz_class l_fact = factorBase.at(p);
        mpz_class p_origin = factorBase.at(p);
        // finche l_fact < radice( numero piu grande da checkare  )
        long exp = 1;

        //cout << "yo"<< endl;
        mpz_class s,h, s_calc, h_calc,s_calc_neg, h_calc_neg, s_calc_temp, h_calc_temp ;

        bool ok_sieve = true;
        while( ok_sieve ){

            mpz_pow_ui(l_fact.get_mpz_t(),p_origin.get_mpz_t(), exp);	// primo elevato a numero
            // cout << "l_fact = "<<l_fact<< " p_origin = "<<p_origin<< endl;
            if (l_fact > 2*m ) ok_sieve = false;
            // tutte le j = s(h) - X_0 + k*l_fact per k < l_fact
            if (exp == 1){						// se il primo è p

            //   cout << "caso exp = 1"<< endl;
                mpz_sqrtm(s.get_mpz_t(),nRSA.get_mpz_t(),p_origin.get_mpz_t()); // ho s e h
                h = l_fact - s;

             //  cout << "s = "<< s << ", h = "<< h << endl;
                s_calc_temp = s - (x_0 % l_fact);
                mpz_powm_ui(s_calc_temp.get_mpz_t(),s_calc_temp.get_mpz_t(),1,l_fact.get_mpz_t());
                h_calc_temp = (h - (x_0 % l_fact)) % l_fact;
                mpz_powm_ui(h_calc_temp.get_mpz_t(),h_calc_temp.get_mpz_t(),1,l_fact.get_mpz_t());

              // cout << "scalc_temp = "<< s_calc_temp << ", hcalc_temp = "<< h_calc_temp << endl;
            }else{

             // cout << "caso exp = 2"<< endl;							// se abbiamo un p^k
            //  cout << "prima..s = "<< s << ", h = "<< h << endl;
                s= hensellemma(s,nRSA, l_fact);				// ho r2 e r2a
                h= l_fact - s;
            //  cout << "dopo..s = "<< s << ", h = "<< h << endl;

                s_calc_temp = (s - (x_0 % l_fact)) % l_fact;
                mpz_powm_ui(s_calc_temp.get_mpz_t(),s_calc_temp.get_mpz_t(),1,l_fact.get_mpz_t());

                h_calc_temp = (h - (x_0 % l_fact)) % l_fact;
                mpz_powm_ui(h_calc_temp.get_mpz_t(),h_calc_temp.get_mpz_t(),1,l_fact.get_mpz_t());

             // cout << "scalc = "<< s_calc_temp << ", hcalc = "<< h_calc_temp << endl;
            }
            mpz_class k = 0;
            while(true){

              //  cout << "whileeeeeee con k ="<< k <<endl;
                s_calc = s_calc_temp  + (l_fact * k );
                h_calc = h_calc_temp  + (l_fact * k);
                s_calc_neg = s_calc_temp + (l_fact * (-k));
                h_calc_neg = h_calc_temp + (l_fact * (-k));

               // cout << "m= "<< m << endl;
              //  cout << "scalc_fin = "<< s_calc << ", hcalc_fin = "<< h_calc << endl;
              //  cout << "scalc_fin_neg = "<< s_calc_neg << ", hcalc_fin = "<< h_calc_neg << endl;
                mpz_class a = s_calc + m, b = h_calc + m, c = s_calc_neg + m, d = h_calc_neg + m;
               // cout << "a= "<< a << ", b = " << b <<" c= "<< c << ", d = " << d << endl;

                if ((a < 0 or a >= 2*m) and (b<0 or b >= 2 * m) and (c < 0 or c >= 2*m) and (d < 0 or d >= 2*m)){
                 //   cout <<"exitt whileeeeeeeee"<<endl;
                    break;
                }
                // applica divisione j positivo
                if (a >=0 and a < 2*m ){
                    #pragma omp critical
                    listYpsilonNumber_test[a.get_si()] = listYpsilonNumber_test[a.get_si()] / factorBase.at(p);
                }
                if (a >= 0 and b < 2*m and b != a){
                    #pragma omp critical
                    listYpsilonNumber_test[b.get_si()] = listYpsilonNumber_test[b.get_si()] / factorBase.at(p);
                }
                if (c >= 0 and c < 2*m and c != b and c !=a){
                    #pragma omp critical
                    listYpsilonNumber_test[c.get_si()] = listYpsilonNumber_test[c.get_si()] / factorBase.at(p);
                }
                if (d >= 0 and d< 2*m and d !=c and d!= b and d != a){
                    #pragma omp critical
                    listYpsilonNumber_test[d.get_si()] = listYpsilonNumber_test[d.get_si()] / factorBase.at(p);
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
   cout << "step 5" << endl;
    gettimeofday(&start_step5, NULL);

    long contRelationOk = 0;
    #pragma omp parallel for
    for (long j = -m ; j < +m ; j++){
        if (listYpsilonNumber_test[j+m] == -1 or listYpsilonNumber_test[j+m] == 1){
            #pragma omp atomic
            contRelationOk++; // mi serve per sapere la dimensione dell array ed essere comodo dopo..
        }
    }

    // controllo che ci siano piu relazioni che numeri nella factorbase
    if (factorBase.size()+1 >  contRelationOk){
          //  cout << endl;
            cout << "ERROR! WE NEED MORE COLUMNS THAN ROWS" << endl;
        return -1;
    }

    mpz_class listYpsilonNumberOk[contRelationOk];
    mpz_class listExponentOk_mod2[contRelationOk][factorBase.size()+1];
    mpz_class listExponentOk[contRelationOk][factorBase.size()+1];

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
    std::vector<mpz_class> jBuoni ;

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

                mpz_class fact_temp = 1;
                bool ok = true;
                do{
                          //  cout << "5c" << endl;
                    mpz_class pmpz = factorBase.at(p);
                    fact_temp = fact_temp * pmpz;
                         //   cout << "listynumber = "<< listYpsilonNumber[j+m] << "fact_temp = " << fact_temp << endl;
                    if( (listYpsilonNumber[j+m] % fact_temp )== 0){
                                //cout << "ok" <<endl;
                                //#pragma omp critical
                        listExponentOk[cont][p+1]++;

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
                //cout << "5e" << endl;
                // cout << "j = " << j;
        }
    }
       // cout << endl;
    gettimeofday(&stop_step5, NULL);
    elapsedTime_step5 += (stop_step5.tv_sec - start_step5.tv_sec) * 1000.0;               // sec to ms
    elapsedTime_step5 += (stop_step5.tv_usec - start_step5.tv_usec) / 1000.0;

  //  cout << "MATRICE MODULO 2"<<endl;
        // stampo a video la nuova matrice modulata!!
//   for (int x = 0 ; x< factorBase.size()+1; x++){
  //          for (int j = 0 ; j < contRelationOk ; j++){

//                cout <<listExponentOk_mod2[j][x] << ", ";
  //          }
    //        cout << endl;
      //   }

     // DA PARALLELIZZARE
    // STEP 6 METODO DI GAUSS
    // Con il metodo di eliminazione di Gauss si determinano alcuni dei vettori v_2(y_i) che danno somma uguale al vettore nullo

   cout << "step 6" << endl;
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
                    mpz_class temp = listExponentOk_mod2[z][row_i];
                    listExponentOk_mod2[z][row_i] = listExponentOk_mod2[z][cont_rig];
                    listExponentOk_mod2[z][cont_rig] = temp;

                       // cout << "z = " << z <<", cont_rig = "<< cont_rig << ", row_i = "<< row_i<<  endl;
                       // cout << "scambio " << listExponentOk_mod2[z][row_i] <<" con "<< listExponentOk_mod2[z][cont_rig] << endl;
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

    mpz_class result = 0;
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

        /*cout<< endl << "VETTORE SOLUZIONE " << cont_ass_1_tot << "= "<<endl;
       / for (long i = 0; i< contRelationOk; i++){

            cout << vect_solution[i] << ", ";
        }
        cout << endl;
        */
        for (long v = 0; v < 2; v++){  // lo facciamo 2 volte, uno per il vettore, e l altro per il suo inverso

            if (v == 1){
               //cout << "vettore opposto"<< endl;

                 for (long i = 0; i< contRelationOk; i++){
                    if (vect_solution_local[i] == 0){
                        vect_solution_local[i] = 1;
                    }else{
                        vect_solution_local[i] = 0;
                    }
                  //  cout << vect_solution[i] << ", ";
                }

            }

        // STEP 7 = calcolo di a e b
        // Si pone a uguale al prodotto degli r_i corrispondenti agli y_i trovati nel passo 6)
        // e si pone b uguale al prodotto delle potenze di p_1,p_2,...,p_t con esponenti uguali
        // alla semisomma degli esponenti della fattorizzazione degli stessi y_i

            mpz_class a = 1;
            mpz_class b = 1;
            long cont = 0;
            vector<mpz_class>::iterator p2;
      //cout << "step 7" << endl;
            mpz_class exponent[factorBase.size()];
            for (long i = 0; i < factorBase.size(); i++){
                exponent[i]=0;
            }
     //  cout << "step 7b" << endl;
            for (p2 = jBuoni.begin(); p2 != jBuoni.end(); p2++){

                if (vect_solution_local[cont] == 1){

                    a = (a * (x_0 + *p2)) % nRSA;
                    //cout << "a =" << a << endl;
                    for (long i = 1; i < factorBase.size()+1; i++){

                        exponent[i-1] = exponent[i-1] + listExponentOk[cont][i];
                    }
                }
                cont++;
            }
        //cout << "a = " << a << endl;
            cont=0;
            for (long p = 0; p <  factorBase.size(); p++){

                if (exponent[cont]>0){
                   //  cout << "b = "<< b << "p = " << *p3 << "^" << exponent[cont] << endl;
                    mpz_class p3mpz = factorBase.at(p);
                    mpz_class temp;
                    temp = exponent[cont] * 2;
                    mpz_powm(temp.get_mpz_t(),p3mpz.get_mpz_t(),temp.get_mpz_t(),nRSA.get_mpz_t());
                    temp = temp * b;
                    mpz_mod(b.get_mpz_t(),temp.get_mpz_t(),nRSA.get_mpz_t());

                    //b = (b * mod_pow(*p3,exponent[cont]/2,nRSA)) % nRSA;
                }
                cont++;
            }
           // cout << endl<< "a = "<< a << ", b = " << b << endl;

            //STEP 8 = CALCOLO SOLUZIONE
            // Si calcola d=mcd(a-b,n) e se 1<d<n allora d è divisore non banale di n,
            // altrimenti si torna al passo 2) con una scelta di k più grande
         //   cout << "step 8" << endl;

            mpz_class risultato_1, risultato_2, temp;

            temp = a+b;
            mpz_gcd(risultato_1.get_mpz_t(),temp.get_mpz_t(),nRSA.get_mpz_t());
            temp = a-b;
            mpz_gcd(risultato_2.get_mpz_t(),temp.get_mpz_t(),nRSA.get_mpz_t());

        //    cout << "step 8a" << endl;
         //   cout << "gcd(a+b,nRSA) = " << risultato_1 << endl;
         //   cout << "gcd(a-b,nRSA) = " << risultato_2 << endl;

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
    mpz_class nRSA, result, x_0;

    if ( argc > 1 ) {

        nRSA = argv[1];
        //long k = pow(exp(sqrt(log(nRSA.get_ui())*log(log(nRSA.get_ui())))),0.35);
        //long m = k * k*k;

        mpz_sqrt(x_0.get_mpz_t(), nRSA.get_mpz_t());
        long k = 10;
        long m = 100;

        bool soluzione_trovata = false;

        int cont_meno1 = 0;
        int cont_meno2 = 0;
        while(!soluzione_trovata){

         //   cout << "lancio con nRSA = "<< nRSA << ", k = "<< k << ", m = " << m << endl;
            result = factorize(nRSA,x_0,k,m);
            if ( result == -2){
                //k = k+10;
                //m=m+30;
                k = k+20;
                m = k*k*k;


                cont_meno2++;
               // cout << "cont mano 1 = " << cont_meno1 << endl;
               // cout << "cont mano 2 = " << cont_meno2 << endl;
                cout << "k =" << k << ", m = " << m << endl;
               // cout << endl;
            }
            else if ( result == -1){
                k = k +5;
                m = k*k*k;
                //k=k+10;
                //m=m+5;
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

