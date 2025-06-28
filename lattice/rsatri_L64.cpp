//  bond percolation on a square lattice of size WIDTH x WIDTH = N.  2N bonds
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <iomanip>
#include <numeric>
#include <random>
#include <cstring>


#define WIDTH 64

#define HEIGHT WIDTH
#define N (WIDTH*HEIGHT)


void nnbdlst(int ncoord, long vertex1[],
	     long vertex2[], long nn[]) ;

using namespace std ; 

struct Properties{
  string name  ; 
  double value ; 
} ;


int		findroot(int i);
long ptr[N] ;

int findroot(int i){
  if (ptr[i]<0) return i;           //returns the address of the root
  return ptr[i] = findroot(ptr[i]); //else returns this statement
}

int main(int argc, char * argv[])
{
  
  int minval = atoi(argv[1]) ; 
  int maxval = atoi(argv[2]) ;
  
  const short Nprops = 8 ; struct Properties prop[Nprops] ;
  
  prop[0].name = "idx"   ; prop[1].name = "   prob"  ; prop[2].name = "     Smax"  ;
  prop[3].name = "       chi"   ; prop[4].name = "    M2t"   ;
  prop[5].name = "      M2"    ;
  prop[6].name = "     nbonds"; prop[7].name = " nsites";

  
  const int printfreq  = 50000000 ;
  const int printfreq2 = 50000000 ;
  const int runsmax    = 50000001 ;
  const int ncoord     = 6 ;
  const int ncd2       = ncoord/2. ;

                         
  static long vertex1[ncd2*N] , vertex2[ncd2*N] , order[ncd2*N], avval[ncd2*N] ;
  static long neighval[ncd2*N], adsorblist[ncd2*N], adsorborderlist[ncd2*N];
  static long val[N], imax[N];
  
  static long nn[ncoord*N] ;
  
  static double smax[ncd2*N],smax2[ncd2*N],M2tot[3*N], M2minus[3*N], binom[3*N] ; 
  static double nbondstot[3*N], nsitestot[3*N];

  long    x, y, nsites,nbonds,i,j ;
  //xo, xp, y, i, j, n, dir, nocc, runs, index, nbonds, nsites;
  long 	 big, r1, r2, v1, v2;
  double  M1, M2, prob, sum, bondfluct, sitefluct;
  long    v10, v20, v1max, v2max, val0, valmax, valtot;
  
  double sumsqrbonds , sumsqrsites ; 
  long runs,n ; 

  nnbdlst(ncoord,vertex1,vertex2,nn) ; 

  long index = 3*N ; 
 
  
  const unsigned int seed = time(0) ;
  mt19937_64 rng(seed);
  uniform_int_distribution<int> dist6(minval,maxval);
  // here the random number are generated into [minval,maxval], maxval inclusive
  
  for(int i=0; i<N; i++)imax[i]= dist6(rng); 
  
  ofstream ofile1,ofile2;
  string filename1 = "rsatrg_min"+to_string(minval)+
    "_max"+to_string(maxval)+"_L"+to_string(WIDTH)+".dat" ;
  
  string filename2 = "rsatrgNN_min"+to_string(minval)+
    "_max"+to_string(maxval)+"_L"+to_string(WIDTH)+"tri.dat" ;

  ofile1.open(filename1) ;
  ofile1 << "#" ;
  for(int iname=0; iname < Nprops; iname++) ofile1 << prop[iname].name << "    " ;
  ofile1 << endl ; 
  ofile1 << "#" ;
  for(int i=0; i < 18*Nprops ; i++) ofile1 << "=" ; 
  ofile1 << endl ;

  
  // initialize arrays and order list
  for(i = 0; i < index; ++i){
    M2tot[i] = M2minus[i] = smax[i] = nbondstot[i] = nsitestot[i] = 0 ;
    avval[i] = adsorblist[i] = adsorborderlist[i] = neighval[i] = 0;
  }

  for (i = 0; i < index; ++i)  order[i] = i;
  
  sumsqrbonds = sumsqrsites = val0 = valmax = valtot = 0;
  for (runs = 1; runs <= runsmax; ++runs){
    
    /* permutation of the order list */
    shuffle(begin(order), end(order),rng);

    //initially each site is a cluster of size 1
    for (i = 0; i < N; i++) {ptr[i] = -1; val[i] = 0;}
    M1 = M2 = N;
        
    for (nsites=nbonds=big=i=0; i<index; i++){
      //i = number of occupied bonds minus 1
	  
      v1 = vertex1[order[i]];
      v2 = vertex2[order[i]];

      if(val[v2] < imax[v2] && val[v1] < imax[v1]){
	
	if (i==0) {v10=v1; v20=v2;}

	v1max = v1;
	v2max = v2;
	if(val[v1]==0) ++nsites;
	if(val[v2]==0) ++nsites;
	++val[v1];
	++val[v2];
	adsorborderlist[nbonds] = i;
	++adsorblist[nbonds];
	++nbonds;  //note, nbonds = number of bonds adsorbed - 1
	r1 = findroot(v1);
	r2 = findroot(v2);    //r1 and r2 are the addresses of the roots

	if (r2 != r1){
                
	  M2 += 2*ptr[r1]*1.0*ptr[r2];
	  // increment the second moment by
	  //(ptr[r1]+ptr[r2])^2-ptr[r1]^2-ptr[r2]^2
	  
	  if (ptr[r1] > ptr[r2]){
	    // size of cluster connected to vertex 1 is LESS THAN
	    //the size of cluster connected to vertex 2
	    
	    ptr[r2] += ptr[r1];
	    // root of cluster connected to vertex 2 remains the root
	    //of the combined cluster
	    
	    ptr[r1] = r2;
	    // link cluster at vertex 1 to the combined root
	    
	    r1 = r2;
	    // variable r1 is the address of the root of the combined cluster

	  } else {

	    ptr[r1] += ptr[r2];
	    // size of cluster connected to vertex 1 is BIGGER THAN
	    //the size of cluster connected to vertex 2
	    ptr[r2] = r1;
	  }

          if (-ptr[r1]>big) big = -ptr[r1];

	}
      }
      nbondstot[i] += nbonds;
      nsitestot[i] += nsites;
      smax[i] += big;
      smax2[i] += big*1.*big ; 
      M2tot[i] += M2;
      M2minus[i] += (M2 - big*1.0*big);
    }
  
    val0 += val[v10]+val[v20];
    valmax += val[v1max]+val[v2max];

    for (i=0;i<nbonds;++i){
      
      v1 = vertex1[order[adsorborderlist[i]]];
      v2 = vertex2[order[adsorborderlist[i]]];
      avval[i] += (val[v1]+val[v2]);
      for (j = 0; j < ncoord; ++j)
	neighval[i] += (val[nn[v1*ncoord + j]] + val[nn[v2*ncoord + j]]);
    }
        
    sumsqrbonds += nbonds*1.0*nbonds;
    sumsqrsites += nsites*1.0*nsites;
    
    //valency of the neighbors
    for (i = 0; i < N; ++i) valtot += val[i];
	
    if ((runs % printfreq2) == 0) {
      
      cout << runs << "  " << val0*0.5/runs << "  "
	   <<  valmax*0.5/runs << "  " << valtot/(runs*1.0*N) << endl ; 
      
      for (i = 0; i < 3*N; ++i) {
	if (adsorblist[i]) {
	      
	  cout << i << "  " << i/(3.*N) << "  "
	       << avval[i]/(adsorblist[i]*2.0) << "  "
	       << neighval[i]/(adsorblist[i]*(2.0*ncoord-2)) << "  "
	       << maxval-(neighval[i]-avval[i])/(adsorblist[i]*(2.0*ncoord-2))
	       << "  " << adsorblist[i]/(runs*1.0) << "  "
	       << adsorblist[i] << endl ; 
	  
	}
      }

      ofile2.open(filename2) ;
      ofile2 <<  fixed << setprecision(8) ;

      for (i = 0; i < 3*N; ++i)
	if (adsorblist[i]){
	  ofile2 << setfill ('0') << std::setw (6) << i ;
	  ofile2 << "   " << i*0.5/N << "   " << avval[i]/(adsorblist[i]*2.0)
		 << "   " << neighval[i]/(adsorblist[i]*(2.0*ncoord-2))
		 << "   "
		 << maxval-(neighval[i]-avval[i])/(adsorblist[i]*(2.0*ncoord-2))
		 << "   " << adsorblist[i]/(runs*1.0)
		 << "   " << adsorblist[i] << endl ;
	}
      ofile2.close();
    }   
    	  
    if ((runs % printfreq) == 0) {

      cout << 1.*val0/runs << "   " << valmax*1.0/runs << endl ; 
      
      //for finding the errors on jamming coverage
      bondfluct = sqrt(sumsqrbonds/(runs*1.0)-nbondstot[index-1]
		       /(runs*1.0)*nbondstot[index-1]/(runs*1.0))/2/N;
          
      sitefluct = sqrt(sumsqrsites/(runs*1.0)-nsitestot[index-1]
		       /(runs*1.0)*nsitestot[index-1]/(runs*1.0))/N;
	  
            
      //sitefluct = sqrt(sumsqrsites/(runs*1.0*N))/(runs*1.0*N);
      cout << WIDTH << "  " << HEIGHT << endl ;
      cout << runs << "  " << nbonds << "  " <<
	nbondstot[index-1]/(runs*3.0*N) << "  " << nsitestot[index-1]/(1.*runs*N) ; 
      cout << bondfluct << "  " << sitefluct << "  "
	   <<  sitefluct/sqrt(runs) << "  " << bondfluct/sqrt(runs) << endl ;
	  
      //ofile1.open(filename1) ;
      ofile1 << fixed << setprecision(8) ;
      
      for (prob = 1.0/index; prob < 1; prob += 1.0/index)  {
	
	n = prob*index;
	    
	ofile1 << setfill ('0') << std::setw (6) << n ;
	ofile1 << "   " <<    prob     << "   "  << smax[n]/(runs*1.0*N)
	       << "   " << sqrt(smax2[n]/(runs*1.0)
				 - smax[n]/(runs*1.0)*smax[n]/(runs*1.0)) << "  "
	       << M2tot[n]/(1.*runs*N) << "  "
	       << M2minus[n]/(runs*1.0*N) << "  "
	       << nbondstot[n]/(runs*1.0*N) << "  "
	       << nsitestot[n]/(runs*1.0*N) << "  " ; 
	
	/*  This is the convolution for the grand canonical (given p) vs.
	    the canonical (fixed number of bonds) ensemble.  Taking same p 
	    values as n/N. (Commented out because it is very slow and not 
	    really  //necessary for large systems!) */
	
	binom[n] = 1;
	for (i = n + 1; i < index; ++i)
	  binom[i] = binom[i-1]*(index+1-i)*1.0/i*prob/(1-prob);
	
	for (i = n - 1; i >= 0; --i)
	  binom[i] = binom[i+1]*(i+1)*1.0/(index-i)*(1-prob)/prob;
	
	for (i = sum = 0; i < index; ++i) sum += binom[i];
	
	for (i = 0; i < index; ++i) binom[i] /= sum;
        
	// un-commented this out, Nov 26 2018
	for (sum = i = 0; i < index-1; ++i)
	  sum += smax[i]*binom[i+1];
	
	ofile1 << sum/(runs*1.0*N) << "  " ;
	
	//for (sum = i = 0; i < index-1; ++i)
	//sum += M2tot[i]*binom[i+1]; fprintf(fp1,"%18.8e", sum/(runs*1.0*N));
	
	for (sum = i = 0; i < index-1; ++i)
	  sum += M2minus[i]*binom[i+1];
	
	ofile1 << sum/(runs*1.0*N) << "  " ;
	//first derivative  (n - N p)/(p - p^2)
	
	for (sum = i = 0; i < index-1; ++i)
	  sum += (i+1 - index*prob)/(prob - prob*prob)*M2minus[i]*binom[i+1];
	
	ofile1 << sum/(runs*1.0*N) << "  " ;
	//second derivative (n^2 + (-1 + N) N p^2 + n (-1 - 2 (-1 + N) p))/((-1 + p)^2 p^2)
	for (sum = i = 0; i < index-1; ++i){
	  sum += ((i+1)*(i+1) + (-1 + index)*index*prob*prob +
		  (i+1)*(-1-2*(-1+index)*prob))/(prob - prob*prob)*
	    M2minus[i]*binom[i+1];
	}
	
	ofile1 << sum/(runs*1.0*N) << endl ;
      }
      printf("done convolution\n");
      ofile1.close();
    }
  }
}


void nnbdlst( int ncoord ,long vertex1[],
	      long vertex2[], long nn[]){
  // define the bonds in the system for a triangular lattice
  
  //#define connect(A,B) {vertex1[index] = A; vertex2[index++] = B;}
  //#define address(X,Y) ((((X)+WIDTH)%WIDTH)+(((Y)+HEIGHT)%HEIGHT)*WIDTH)
  
  int x,y,index ; 
  for (index = x = 0; x < WIDTH; ++x){
    for (y = 0; y < HEIGHT; ++y) {

      int v1 = x + WIDTH*y;
      
      int v2 =  ((((x+1)+WIDTH)%WIDTH)+(((y)+HEIGHT)%HEIGHT)*WIDTH) ; 
            
      vertex1[index]   = v1 ;
      vertex2[index++] = v2 ;
      
      nn[v1*ncoord + 0]=v2 ;
      
      v2 = ((((x)+WIDTH)%WIDTH)+(((y+1)+HEIGHT)%HEIGHT)*WIDTH) ;
      
      vertex1[index]   = v1 ;
      vertex2[index++] = v2 ;
      
      nn[v1*ncoord + 1] = v2 ; 
            
      v2 = ((((x+1)+WIDTH)%WIDTH)+(((y+1)+HEIGHT)%HEIGHT)*WIDTH) ;
      
      vertex1[index]   = v1 ;
      vertex2[index++] = v2 ;
      
      nn[v1*ncoord + 2] = v2 ;  
      nn[v1*ncoord + 3] = ((((x)+WIDTH)%WIDTH)+(((y-1)+HEIGHT)%HEIGHT)*WIDTH) ;
      //address(x,y-1);
      nn[v1*ncoord + 4] = ((((x-1)+WIDTH)%WIDTH)+(((y)+HEIGHT)%HEIGHT)*WIDTH) ;
      //address(x-1,y);
      nn[v1*ncoord + 5] = ((((x-1)+WIDTH)%WIDTH)+(((y-1)+HEIGHT)%HEIGHT)*WIDTH) ;
      //address(x-1,y-1);
      
    }
  }

  if(index != (ncoord/2)*N) {std::cout<< "Error num. of  bonds "; exit(EXIT_FAILURE);}
  
  return ; 
}
