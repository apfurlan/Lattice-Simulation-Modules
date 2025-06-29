//  bond percolation on a square lattice of size WIDTH x WIDTH = N.  2N bonds
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <numeric>
#include <climits>

#define LX         64
#define LY         64
#define N          (LX*LY)   
#define NEDGES     (2*N)
#define nRealz     100000 
#define maxValIn   2  
#define maxValOut  2  

//====== MERSENNE TWISTER generator ===================
#define NN 624
#define MM 397
#define MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define UPPER_MASK 0x80000000UL /* most significant w-r bits */
#define LOWER_MASK 0x7fffffffUL /* least significant r bits */

static unsigned long mt[NN]; /* the array for the state vector  */
static int mti=NN+1; /* mti==N+1 means mt[N] is not initialized */

void init_genrand(unsigned long ss) ;
unsigned long genrand_int32(void)  ;
double genrand_real2(void)         ;
//==========================================================

using namespace std ; 

void nnBondList(int Lx, int Ly, long nEdges, 
		long vertex1[],long vertex2[]) ; 

int findroot(int i) ;
int findrootPlus(int j) ;
//int findrootMinus(int k) ;

long ptr[N] ; 
int findroot(int i){
  if (ptr[i] < 0) return i ;
  return ptr[i] = findroot(ptr[i]) ;
}


long ptrPlus[N] ; 
int findrootPlus(int j){
    if (ptrPlus[j] < 0) return j ;
    return ptrPlus[j] = findrootPlus(ptrPlus[j]) ;
}


//long ptrMinus[N] ; 
//int findrootMinus(int k){
//  if (ptrMinus[k] < 0) return k ;
//    return ptrMinus[k] = findrootMinus(ptrMinus[k]) ;
//}


int main(int argc, char * argv[]){

  static long vertex1[NEDGES] ;
  static long vertex2[NEDGES] ;
  static int order[NEDGES]    ;  
  
  static int valIn[N]  ;
  static int valOut[N] ;
 
  static double smax[NEDGES]   ;
  static double smax2[NEDGES]  ;
  static double smax3[NEDGES]  ;
  static double smax4[NEDGES]  ;

  static double smaxPlus[NEDGES]   ;
  static double smaxPlus2[NEDGES]  ;
  static double smaxPlus3[NEDGES]  ;
  static double smaxPlus4[NEDGES]  ;

  static double smaxMinus[NEDGES]   ;
  static double smaxMinus2[NEDGES]  ;
  static double smaxMinus3[NEDGES]  ;
  static double smaxMinus4[NEDGES]  ;

  static double M2TotalPlus[NEDGES] ;
  static double M2PrimePlus[NEDGES] ;

  static double M2TotalMinus[NEDGES] ;
  static double M2PrimeMinus[NEDGES] ;
    
  static double M2Total[NEDGES]  ;
  static double M2Prime[NEDGES]  ;
  
  static double nBondsPlus[NEDGES]  ;
  static double nBondsMinus[NEDGES] ;

  ofstream ofile2;
  
  int ibdm,ibdp ; // i-th bond minus/plus
  bool addBond  ; // if the rule of valencies is satisfied

  int r1,r2,r1p,r2p,r1m,r2m ;  // the root of cluster left/righth
  int sumbondsIn, sumbondsOut ; 
  int direction , ranDir ;
  //int sumBondsPlus, sumBondsMinus ; 
  //double nBondsPos, nBondsNeg; 
  int   i , v1,v2,n,runs,big,bigPlus,bigMinus ; //, n, runs, nsites;
  //int bondPos , bondNeg ; 
  //long 	big, r1, r2, v1, v2 ;
  
  unsigned long int seed = atoi(argv[1]); 


  nnBondList(LX,LY,NEDGES,vertex1,vertex2) ; 
  init_genrand(seed) ;

  string filename ="bondsPerRun_L"+to_string(LX)+".dat" ;
  ofile2.open(filename) ;
  ofile2 << "#   run          nBonds+         nBonds- " << endl ;
  ofile2 << "#" ;
  for(int i=0; i<50 ; i++) ofile2 << "=" ; 
  ofile2 << endl ; 

  
  
  for (i = 0; i < NEDGES ; i++) {
    M2Total[i] = M2Prime[i] = smax[i] = smax2[i] = smax3[i] = smax4[i] = 0. ;
    M2TotalMinus[i] = M2PrimeMinus[i] =  0 ;
    M2TotalPlus[i]  = M2PrimePlus[i]  =  0 ;
    smaxMinus[i] = smaxMinus2[i] =  smaxMinus3[i] = smaxMinus4[i] = 0. ;
    smaxPlus[i]  = smaxPlus2[i]  =  smaxPlus3[i]  = smaxPlus4[i]  = 0. ;
    order[i]=i ;
  }
  
  //nBondsPos = nBondsNeg = 0 ;
  for (runs = 0; runs < nRealz; runs++ ){ //loop over realizations

    for(int i=0; i< N; i++) {
      valIn[i] = valOut[i] = 0 ;
      ptr[i] = ptrPlus[i] = -1 ;
    } 

    srand(genrand_int32()); 
    random_shuffle(begin(order), end(order)); //shuffle order array
    

    double M2     = N ;
    double M2Plus = N ;
    double M2Minus= N ; 
    for (ibdp=ibdm=big=bigPlus=bigMinus=i=0; i<NEDGES; i++){

      v1 = vertex1[order[i]];
      v2 = vertex2[order[i]];

      ranDir = 2*round(genrand_real2())-1 ;
      addBond = false ;

      for(int nAttempts = 0; nAttempts < 2 ; nAttempts++){
	
	direction = ranDir - nAttempts*(2*ranDir) ;
	
	if(direction > 0) {
	  sumbondsIn  = valIn[v2]  ; 
	  sumbondsOut = valOut[v1] ;
	} else {
	  sumbondsIn  = valIn[v1]  ; 
	  sumbondsOut = valOut[v2] ;
	}

	if (sumbondsIn < maxValIn && sumbondsOut < maxValOut ) {
	  addBond = true ;
	  break ; 
	}

      }
      
      if(addBond){
		
	if(direction > 0){
	  valIn[v2]++    ; 
	  valOut[v1]++   ;
	  ibdp++         ; 
	} else {
	  valIn[v1]++    ; 
	  valOut[v2]++   ;
	  ibdm++         ; 
	}

	r1 = findroot(v1);
	r2 = findroot(v2); 

	// if bond added is positive ;
	if(direction > 0){
	  
	  r1p = findrootPlus(v1) ;
	  r2p = findrootPlus(v2) ;
	  
	  if(r2p != r1p){
	    
	    M2Plus += 2.*ptrPlus[r1p]*ptrPlus[r2p] ; 
	    
	    if(ptrPlus[r1p] > ptrPlus[r2p]){
	      
	      ptrPlus[r2p] += ptrPlus[r1p] ;
	      ptrPlus[r1p] = r2p ;
	      r1p = r2p ;
	      
	    } else {
		
	      ptrPlus[r1p] += ptrPlus[r2p] ;
	      ptrPlus[r2p] = r1p ;
	          
	    }
	  }
	  
	  if(-ptrPlus[r1p] > bigPlus) bigPlus = -ptrPlus[r1p] ;  

	  //if negative ; 
	} else {
	  r1m =0 ;
	  //r1m = findrootMinus(v1) ;
	  //r2m = findrootMinus(v2) ;

	  //if(r2m != r1m){

	  //M2Minus += 2.*ptrMinus[r1m]*ptrMinus[r2m] ;

	  //if(ptrMinus[r1m] > ptrMinus[r2m]){

	  //  ptrMinus[r2m] += ptrMinus[r1m] ;
	  //  ptrMinus[r1m] = r2m ;
	  //  r1m = r2m ; 
	      
	  //} else {

	  //  ptrMinus[r1m] += ptrMinus[r2m] ;
	  //   ptrMinus[r2p] = r1p ; 
	        
	  // }
	    
	  //}

	  //if(-ptrMinus[r1p] > bigMinus) bigMinus = -ptrMinus[r1] ;   
	}
	
	if(r2 != r1){

	  M2  += 2.*ptr[r1]*ptr[r2] ; 

	  if (ptr[r1] > ptr[r2]) {
	    
	    ptr[r2] += ptr[r1] ;
	    ptr[r1] = r2 ;
	    r1 = r2 ; 
	    
	 } else {

	    ptr[r1] += ptr[r2] ;
	    ptr[r2] = r1 ; 
	    
	  } // percolation - r2 != r1

	  if(-ptr[r1] > big) big = -ptr[r1] ;  
	}
		  
      } 
	
      nBondsPlus[i]  += ibdp ;
      nBondsMinus[i] += ibdm ;

      smaxPlus[i]  += bigPlus ;
      smaxPlus2[i] += bigPlus*1.*bigPlus ;
      smaxPlus3[i] += bigPlus*1.*bigPlus*1.*bigPlus ;
      smaxPlus4[i] += bigPlus*1.*bigPlus*1.*bigPlus*1.*bigPlus ;

      M2TotalPlus[i] += M2Plus ; 
      M2PrimePlus[i] += (M2Plus - bigPlus*1.*bigPlus) ; 
      
      smaxMinus[i]  += bigMinus ;
      smaxMinus2[i] += bigMinus*1.*bigMinus ;
      smaxMinus3[i] += bigMinus*1.*bigMinus*1.*bigMinus ;
      smaxMinus4[i] += bigMinus*1.*bigMinus*1.*bigMinus*1.*bigMinus ;

      M2TotalMinus[i] += M2Minus ;
      M2PrimeMinus[i] += (M2Minus - bigMinus*1.*bigMinus);

      
      smax[i]   += big ;
      smax2[i]  += big*1.*big ; 
      smax3[i]  += big*1.*big*1.*big ;
      smax4[i]  += big*1.*big*1.*big*1.*big ; 

      M2Total[i]+=  M2 ;
      M2Prime[i]+= (M2 - big*1.0*big) ; 
      
    } // ========= end loop over bonds ===================	
    
    ofile2 << fixed << setprecision(8) ; 
    ofile2 << setfill ('0') << std::setw (6) << runs ;
    ofile2 << "    " << scientific
	   << 1.*ibdp/(1.*NEDGES)  << "   "
	   << 1.*ibdm/(1.*NEDGES)  << "   " << endl ; 
  } // ======= end loop over replicas ================

  ofile2.close();
  

  ofstream ofile;
  filename ="icePlus_L"+to_string(LX)+".dat" ;
  ofile.open(filename) ;

  ofile << "#idx         prob        Smax      Smax2      Smax3     Smax4      M2T      M2M      nBonds+         nBonds- " << endl ;
  ofile << "#" ;
  for(int i=0; i<120 ; i++) ofile << "=" ; 
  ofile << endl ; 
  
  for (double prob = 1.0/NEDGES; prob < 1; prob += 1.0/NEDGES) {
    
    n = prob*NEDGES;
    
    ofile << fixed << setprecision(8) ; 
    ofile << setfill ('0') << std::setw (6) << n ;
    ofile << "    " << prob << "  "
	  << scientific
	  << smax[n]/(1.*runs)   << "  "
	  << smax2[n]/(runs*1.)  << "  "
	  << smax3[n]/(runs*1.)  << "  "
	  << smax4[n]/(runs*1.)  << "  " 
	  << M2Total[n]/(runs*1.*N)  << "  "
	  << M2Prime[n]/(runs*1.*N)  << "  "
	  << smaxPlus[n]/(1.*runs)   << "  "
	  << smaxPlus2[n]/(runs*1.)  << "  "
	  << smaxPlus3[n]/(runs*1.)  << "  "
	  << smaxPlus4[n]/(runs*1.)  << "  " 
	  << M2TotalPlus[n]/(runs*1.*N)  << "  "
	  << M2PrimePlus[n]/(runs*1.*N)  << "  "
	  << nBondsPlus[n]/(1.*runs*NEDGES)  << "    "
	  << nBondsMinus[n]/(1.*runs*NEDGES) << "   " << endl ; 
    
  }
  
  ofile.close();
    
  return 0 ;
}

//======================================================================
//================ FUNCTIONS AND ROUTINES ==============================
//======================================================================
void nnBondList(int Lx, int Ly, long nEdges, 
		long vertex1[],long vertex2[]){

  int v1,v2 ; 
  int index ; 
  
  ofstream outvertex ;
  
  string outname =  "vertexlist"+to_string(Lx)+"x"+to_string(Ly)+".dat" ; 
  outvertex.open(outname) ; 

  index = 0 ; 
  for (int x = 0; x < Lx; ++x) {
    for (int y = 0; y < Ly; ++y) {
      
      v1 = x + Lx*y;
      v2 = ((x + 1) % Lx) + Lx*y;  //horizontal bond

      vertex1[index]  = v1 ;
      vertex2[index]  = v2 ;
      
      outvertex << index << "  " << x  << "  " << y << "  "
		<<  x+1  << "  " << y  << endl ;
      
      index++ ;
      

      v2 = x + Lx*((y+1) % Ly);   //vertical bond

      vertex1[index]  = v1 ;
      vertex2[index]  = v2 ;

      outvertex << index << "  " <<  x  << "  " << y << "  "
		<<   x   << "  " << y+1 << endl ;
      
      index++ ; 
      
    }
  }
  outvertex.close(); 

  //cout << index << "  " << nEdges << endl ; 

  
  if(index != nEdges){
    std::cout<< "Error num. of  bonds " << endl ;
    exit(EXIT_FAILURE);
  }

  return ; 
}
  

/* initializes mt[N] with a seed */
void init_genrand(unsigned long ss){
    mt[0]= ss & 0xffffffffUL;
    for (mti=1; mti<NN; mti++) {
        mt[mti] =
      (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti);
        /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
        /* In the previous versions, MSBs of the seed affect   */
        /* only MSBs of the array mt[].                        */
        /* 2002/01/09 modified by Makoto Matsumoto             */
        mt[mti] &= 0xffffffffUL;
        /* for >32 bit machines */
    }
}


/* generates a random number on [0,0xffffffff]-interval */
unsigned long genrand_int32(void)
 {
    unsigned long y;
    static unsigned long mag01[2]={0x0UL, MATRIX_A};
    /* mag01[x] = x * MATRIX_A  for x=0,1 */

    if (mti >= NN) { /* generate N words at one time */
        int kk;

        if (mti == NN+1)   /* if init_genrand() has not been called, */
            init_genrand(5489UL); /* a default initial seed is used */

        for (kk=0;kk<NN-MM;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+MM] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        for (;kk<NN-1;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+(MM-NN)] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        y = (mt[NN-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
        mt[NN-1] = mt[MM-1] ^ (y >> 1) ^ mag01[y & 0x1UL];

        mti = 0;
    }

    y = mt[mti++];

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return y;
}

/* generates a random number on [0,1)-real-interval */
double genrand_real2(void)
{
    return genrand_int32()*(1.0/4294967296.0);
    /* divided by 2^32 */
}


