/* Cohesion tradeoff when all new nests are of equal quality.
   Simulation continues until all ants have emigrated to new nests.

Mersenne Twister pseudo-random number generator code mt19937ar.c is needed. 
One can download mt19937ar.c from various websites.
One needs to correct the folder's name in
    #include "../lib/stat/mt19937ar.c"
below to your local folder's name.

In the following, p={a,b,c} means that one needs to run the code with parameter p=a, p=b, and p=c, and collect results.

Produce the results shown in Figure 9 by running

    a.out 0.1 z={0.02,0.04,0.06,0.08,0.1,0.12,0.14,0.16,0.18,0.2,0.22,0.24,0.26,0.28,0.3,0.32,0.34,0.36,0.38,0.4} 100 2

    a.out 0.1 z={0.04,0.08,0.12,0.16,0.2,0.24,0.28,0.32,0.36,0.4} 100 4

    a.out 0.1 z={0.06,0.12,0.18,0.24,0.30,0.36} 100 6

    a.out {0.01,0.01778,0.03162,0.05623,0.1,0.1778,0.3162,0.5623,1} z=0.12 100 {2,4,6}

Produce the results shown in Figure S6 by running

    a.out {0.01,0.1,1} z={0.02,0.04,0.06,0.08,0.1,0.12,0.14,0.16,0.18,0.2,0.22,0.24,0.26,0.28,0.3,0.32,0.34,0.36,0.38,0.4} 100 2 with alpha_leak=0

    a.out {0.01,0.1,1} z={0.04,0.08,0.12,0.16,0.2,0.24,0.28,0.32,0.36,0.4} 100 4 with alpha_leak=0

    a.out {0.01,0.1,1} z={0.06,0.12,0.18,0.24,0.30,0.36} 100 6 with alpha_leak=0

    a.out {0.01,0.01778,0.03162,0.05623,0.1,0.1778,0.3162,0.5623,1} z={0.12,0.36} 100 {2,4,6} with alpha_leak=0

Produce the results shown in Figure S7 by running

    a.out {0.2,1} z={0.02,0.04,0.06,0.08,0.1,0.12,0.14,0.16,0.18,0.2,0.22,0.24,0.26,0.28,0.3,0.32,0.34,0.36,0.38,0.4} 100 2

    a.out {0.02,1} z={0.04,0.08,0.12,0.16,0.2,0.24,0.28,0.32,0.36,0.4} 100 4

    a.out {0.02,1} z={0.06,0.12,0.18,0.24,0.30,0.36} 100 6

    a.out {0.01,0.01778,0.03162,0.05623,0.1,0.1778,0.3162,0.5623,1} z=0.36 100 {2,4,6}


*/
#include <iostream>
using namespace std;
#include <fstream>
#include <cstdlib> // atoi
#include <cmath> // cos, sin
#include "../lib/stat/mt19937ar.c"
#include <ctime>

int main (int argc, char **argv) {

  if (argc != 5) {
    cerr << "cohesion-tradeoff.out alpha z Na Nnest" << endl;
    cerr << "alpha: rate at which committed ants convert to recruits, common to low-threshold and high-threshold ants" << endl;
     cerr << "z: initial fraction of recruiters" << endl;
     cerr << "Na: number of ants" << endl;
     cerr << "Nnest: number of new nests" << endl;
    exit(8);
  }

  init_genrand(time(NULL));
  double alpha = atof(argv[1]); // conversion rate from committed to recruit
  double alpha_leak = 0.05; // = 0.0 for fig. S6. = 0.05 otherwise
  cerr << "leak rate = " << alpha_leak << endl;
  double z = atof(argv[2]); // initial fraction of recruiters
  int Na = atoi(argv[3]);
  int tmp_sum;

  double t;
  double t_final_ave = 0.0;
  double t_final_std = 0.0;
  double cohesion_tmp;  
  double cohesion_ave = 0.0;
  double cohesion_std = 0.0;
  double Entropy; // entropy (working var)

  int i,j;
  int trials = 10000;
  int tr;
  int Nnest = atoi(argv[4]); // # new nests
  double rate[4*Nnest], rate_sum;
  double ra;

  int oldnest, com[Nnest], rec[Nnest]; // # ants in each category

  tr = 0;
  while (tr < trials) {
    // initialization
    t = 0.0;
    oldnest = (int)(Na*(1-z)+1e-8);
    tmp_sum = oldnest;
    for (i=0 ; i<Nnest ; i++) {
      com[i] = (int)((Na*z+1e-8)/Nnest);
      rec[i] = 0;
      tmp_sum += com[i];
    }
    if (tmp_sum < Na) {
      if (tr==0)
	cerr << "Na - tmp_sum = " << Na - tmp_sum << " " << oldnest << " " << com[0] << endl;
      oldnest += Na - tmp_sum;
    }

    while (oldnest < Na && oldnest > 0.1 * Na) {

	rate_sum = 0.0;
	for (i=0 ; i<Nnest ; i++) {
	  rate[4*i] = (double)rec[i] * oldnest / Na;
	  rate[4*i+1] = alpha * com[i];
	  rate[4*i+2] = alpha_leak * com[i];
	  rate[4*i+3] = alpha_leak * rec[i];
	  for (j=0 ; j<4 ; j++)
	    rate_sum += rate[4*i+j];
	}
	ra = (genrand_int32()+0.5)/4294967296.0 * rate_sum; // 0 < ra < rate_sum

	i=0;
	while (ra > rate[i]) {
	  ra -= rate[i];
	  i++;
	}
	if (i >= 4*Nnest) {
	  cerr << "error in ra" << endl;
	  exit(8);
	}

	if (i % 4 == 0) { // recruitment occurs
	  oldnest--;
	  com[i/4]++;
	} else if (i % 4 == 1) { // a committed ant turns to recruit
	  com[i/4]--;
	  rec[i/4]++;
	} else if (i % 4 == 2) { // leak from committed to oldnest
	  com[i/4]--;
	  oldnest++;
	} else { // leak from recruit to oldnest
	  rec[i/4]--;
	  oldnest++;
	}
	t += -1.0/rate_sum*log((genrand_int32()+0.5)/4294967296.0);
    } // dynamics done

    if (oldnest != Na) { // emigration done
      tr++;
      t_final_ave += t;
      t_final_std += t*t;
      tmp_sum = 0;
      for (i=0 ; i<Nnest ; i++)
	tmp_sum += com[i] + rec[i];

      Entropy = 0.0;
      for (i=0 ; i<Nnest ; i++) {
	if (com[i]+rec[i]>0)
	  Entropy += - (double)(com[i]+rec[i])/tmp_sum * log((double)(com[i]+rec[i])/tmp_sum);
#ifdef VERBOSE
	cerr << com[i]+rec[i] << " ";
#endif // VERBOSE
      } 
	
#ifdef VERBOSE
      cerr << tmp_sum << " " << 1 - Entropy/log(Nnest) << endl;
#endif // VERBOSE
      cohesion_tmp = 1 - Entropy/log(Nnest);
      cohesion_ave += cohesion_tmp;
      cohesion_std += cohesion_tmp*cohesion_tmp;
    }
  } // a single trial done

      t_final_ave /= trials;
      t_final_std = sqrt(t_final_std/trials - t_final_ave*t_final_ave);
      cohesion_ave /= trials;
      cohesion_std = sqrt(cohesion_std/trials - cohesion_ave*cohesion_ave);
      
      cout << alpha << " " << z << " " << t_final_ave << " " << t_final_std << " " << cohesion_ave << " " << cohesion_std << endl;

  return 0;
}
