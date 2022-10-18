/* nest-choice numerical simulations

Mersenne Twister pseudo-random number generator code mt19937ar.c is needed. 
One can download mt19937ar.c from various websites.
One needs to correct the folder's name in
    #include "../lib/stat/mt19937ar.c"
below to your local folder's name.

In the following, p={a,b,c} means that one needs to run the code with parameter p=a, p=b, and p=c, and collect results.

Produce the results shown in Figure 4 by running
        a.out alpha=0.1 alpha_s=0.1 H={0.06667,0.13334,0.2,0.26667,0.33334,0.4} z=0.3 threshold=0.5 Na=100

Produce the results shown in Figure 7 by running
        a.out alpha=0.1 alpha_s=0.1 H=0.2 z={0.1,0.2,0.3,0.4} threshold=0.5 Na=100

Produce the results shown in Figure S2 by running
        a.out alpha=0.1 alpha_s=0.1 H=0.2 z=0.3 threshold={0.25,0.3,0.35,0.4,0.45,0.5} Na=100

Produce the results shown in Figure S4 by running
        a.out alpha=0.1 alpha_s={0,0.1,0.2,0.3,0.4,0.5} H=0.2 z=0.3 threshold=0.5 Na=100 
*/


#include <iostream>
using namespace std;
#include <fstream>
#include <cstdlib> // atoi
#include <cmath> // cos, sin
#include "../lib/stat/mt19937ar.c"
#include <ctime>

int compare_finite_nestchoice (const void *a, const void *b) {
  return (*(double*)a-*(double*)b);
}

int main (int argc, char **argv) {

  if (argc != 7) {
    cerr << "finite-nestchoice.out alpha alpha_s H z threshold Na" << endl;
    cerr << "alpha: rate at which committed ants convert to recruits, common to low-threshold and high-threshold ants" << endl;
    cerr << "alpha_s: rate at which high-threshold ants visiting the poor nest moves to the good nest" << endl;
    cerr << "H: fraction of high-threshold ants" << endl;
     cerr << "z: initial fraction of recruiters" << endl;
     cerr << "threshold: quorum threshold, between 0 and 1" << endl;
     cerr << "Na: number of ants" << endl;
    exit(8);
  }

  init_genrand(time(NULL));
  double alpha = atof(argv[1]); // conversion rate from committed to recruit
  double alpha_g = alpha;
  double alpha_p = alpha;
  double alpha_s = atof(argv[2]); // rate at which high-threshold ants visiting the poor nest move to the good nest
  double alpha_leak = 0.05;
  cerr << "leak rate = " << alpha_leak << endl;
  double H = atof(argv[3]); // fraction of high-threshold ants
  double L = 1-H; // fraction of low-threshold ants
  double z = atof(argv[4]); // initial fraction of recruiters
  double th_quorum_frac = atof(argv[5]); // quorum threshold (normalized by N_ant)
  int Na = atoi(argv[6]);
  int th_quorum = (int)(th_quorum_frac*Na); // quorum threshold
  int tmp_sum;

  int trials = 10000;
  int tr;

  cerr << "Quorum threshold (fractional) = " << th_quorum_frac << endl;
  double precision = 0.0; // % correct collective decision

  double t;
  //  double t_quorum_list[trials]; // time to quorum for each trial
  double t_quorum_ave = 0.0;
  double t_quorum_std = 0.0;
  
  int i;
  double accum_rate[15], ra;


  int L_oldnest, H_oldnest, L_poor_com, L_poor_rec, H_poor_vis, L_good_com, H_good_com, L_good_rec, H_good_rec;

  tr = 0;
  while (tr < trials) {
    // initialization
    t = 0.0;
    H_oldnest = (int)(Na*H*(1-z)+1e-8);
    L_oldnest = (int)(Na*(1-z)+1e-8) - H_oldnest;
    H_poor_vis = H_good_com = (int)((Na*H*z+1e-8)/2);
    L_poor_com = (int)((Na*z+1e-8)/2) - H_poor_vis;
    L_good_com = (int)((Na*z+1e-8)/2) - H_good_com;
    L_poor_rec = L_good_rec = H_good_rec = 0;
    tmp_sum = H_oldnest + L_oldnest + H_poor_vis + L_poor_com + H_good_com + L_good_com;
    if (tmp_sum < Na)
      H_oldnest += Na - tmp_sum;
    if (tr==0)
      cerr << "Na - tmp_sum = " << Na - tmp_sum << "; " << H_oldnest << " " << L_oldnest << " " << H_poor_vis << " " << L_poor_com << " " << H_good_com << " " << L_good_com << endl;

    while (L_oldnest + H_oldnest < Na && L_good_com + H_good_com + L_good_rec + H_good_rec < th_quorum && L_poor_com + L_poor_rec + H_poor_vis < th_quorum) {

#ifdef VERBOSE
	if ((int)(t*100) > (int)((t-dt)*100))
	  cout << t << " " << (double)good_com/Na << " " << (double)(good_com + good_rec)/Na << " " << (double)L_poor_com/Na << " " << (double)(L_poor_com + L_poor_rec)/Na << endl;
#endif // VERBOSE

	accum_rate[0] = (double)L_poor_rec * L_oldnest / Na;
      	accum_rate[1] = accum_rate[0] + (double)L_poor_rec * H_oldnest / Na;
	accum_rate[2] = accum_rate[1] + (double)(L_good_rec + H_good_rec) * L_oldnest / Na;
	accum_rate[3] = accum_rate[2] + (double)(L_good_rec + H_good_rec) * H_oldnest / Na;
	accum_rate[4] = accum_rate[3] + alpha_p * L_poor_com; // commited -> recruiter
	accum_rate[5] = accum_rate[4] + alpha_g * L_good_com;
	accum_rate[6] = accum_rate[5] + alpha_g * H_good_com;
	accum_rate[7] = accum_rate[6] + alpha_s * H_poor_vis;
	accum_rate[8] = accum_rate[7] + alpha_leak * L_poor_com;
	accum_rate[9] = accum_rate[8] + alpha_leak * H_poor_vis;
	accum_rate[10] = accum_rate[9] + alpha_leak * L_poor_rec;
	accum_rate[11] = accum_rate[10] + alpha_leak * L_good_com;
	accum_rate[12] = accum_rate[11] + alpha_leak * H_good_com;
	accum_rate[13] = accum_rate[12] + alpha_leak * L_good_rec;
	accum_rate[14] = accum_rate[13] + alpha_leak * H_good_rec;

	ra = (genrand_int32()+0.5)/4294967296.0 * accum_rate[14];

	if (ra < accum_rate[0]) {
	  L_oldnest--;
	  L_poor_com++;
	} else if (ra < accum_rate[1]) {
	  H_oldnest--;
	  H_poor_vis++;
	} else if (ra < accum_rate[2]) {
	  L_oldnest--;
	  L_good_com++;
	} else if (ra < accum_rate[3]) {
	  H_oldnest--;
	  H_good_com++;
	} else if (ra < accum_rate[4]) {
	  L_poor_com--;
	  L_poor_rec++;
	} else if (ra < accum_rate[5]) {
	  L_good_com--;
	  L_good_rec++;
	} else if (ra < accum_rate[6]) {
	  H_good_com--;
	  H_good_rec++;
	} else if (ra < accum_rate[7]) { 
	  H_poor_vis--;
	  H_good_com++;
	} else if (ra < accum_rate[8]) { // leak
	  L_poor_com--;
	  L_oldnest++;
	} else if (ra < accum_rate[9]) {
	  H_poor_vis--;
	  H_oldnest++;
	} else if (ra < accum_rate[10]) {
	  L_poor_rec--;
	  L_oldnest++;
	} else if (ra < accum_rate[11]) {
	  L_good_com--;
	  L_oldnest++;
	} else if (ra < accum_rate[12]) {
	  H_good_com--;
	  H_oldnest++;
	} else if (ra < accum_rate[13]) {
	  L_good_rec--;
	  L_oldnest++;
	} else if (ra < accum_rate[14]) {
	  H_good_rec--;
	  H_oldnest++;
	}
	t += -1.0/accum_rate[14]*log((genrand_int32()+0.5)/4294967296.0);
    } // dynamics done


      if (L_oldnest + H_oldnest < Na) { // quorum reached in either new site
	//	t_quorum_list[tr] = t;
	tr++;
	t_quorum_ave += t;
	t_quorum_std += t*t;
	if (L_good_com + H_good_com + L_good_rec + H_good_rec == th_quorum)
	  precision += 1.0;
      }
  } // all trials done
      
  //  qsort(t_quorum_list,trials,sizeof(double),compare_finite_nestchoice);
  
  t_quorum_ave /= trials;
  t_quorum_std = sqrt(t_quorum_std/trials - t_quorum_ave*t_quorum_ave);
  precision /= trials;
  cout << H << " " << alpha_s << " " << z << " " << th_quorum_frac << " " << t_quorum_ave << " " << t_quorum_std << " " << precision << endl;
      // cout << " " << t_quorum_list[(int)(trials*0.05)] << " " << t_quorum_list[(int)(trials*0.25)] << " " << t_quorum_list[(int)(trials*0.5)] << " " << t_quorum_list[(int)(trials*0.75)] << " " << t_quorum_list[(int)(trials*0.95)] << endl;

  return 0;
}
