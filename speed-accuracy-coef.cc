/* Speed-accuracy trade-off: Pearson correlation coefficient between the time to quorum and the accuracy when one parameter is varied.

Mersenne Twister pseudo-random number generator code mt19937ar.c is needed. 
One can download mt19937ar.c from various websites.
One needs to correct the folder's name in
    #include "../lib/stat/mt19937ar.c"
below to your local folder's name.

In the following, p={a,b,c} means that one needs to run the code with parameter p=a, p=b, and p=c, and collect results.

Produce the results shown in Figure 6 by running
        a.out speed-accuracy-coef.out 0
        Note: We need to run this code three times and collect the results by setting alpha_s={0.01,0.1,1} within the first "if (to_vary==0)" loop

Produce the results shown in Figure 8 by running
        a.out speed-accuracy-coef.out 1
        Note: We need to run this code three times and collect the results by setting alpha_s={0.01,0.1,1} within the first "if (to_vary==0)" loop

Produce the results shown in Figure S1 by running

        a.out speed-accuracy-coef.out 0
        Note: We need to run this code three times and collect the results by setting alpha_s={0.01,0.1,1} within the first "if (to_vary==0)" loop and setting alpha_leak=0

        a.out speed-accuracy-coef.out 1
        Note: We need to run this code three times and collect the results by setting alpha_s={0.01,0.1,1} within the first "if (to_vary==1)" loop and setting alpha_leak=0

        a.out speed-accuracy-coef.out 3
        Note: We need to run this code three times and collect the results by setting alpha_s={0.01,0.1,1} within the first "if (to_vary==3)" loop and setting alpha_leak=0

        a.out speed-accuracy-coef.out 2
        Note: We need to run this code three times and collect the results by setting z={0.1,0.3} within the first "if (to_vary==2)" loop and setting alpha_leak=0

Produce the results shown in Figure S3 by running
        a.out speed-accuracy-coef.out 3
        Note: We need to run this code three times and collect the results by setting alpha_s={0.01,0.1,1} within the first "if (to_vary==3)" loop

Produce the results shown in Figure S5 by running
        a.out speed-accuracy-coef.out 2
        Note: We need to run this code three times and collect the results by setting z={0.1,0.3} within the first "if (to_vary==2)" loop

*/
#include <iostream>
using namespace std;
#include <fstream>
#include <cstdlib> // atoi
#include <cmath> // cos, sin
#include "../lib/stat/mt19937ar.c"
#include <ctime>

int main (int argc, char **argv) {

  if (argc != 2) {
    cerr << "speed-accuracy-coef.out to_vary" << endl;
    cerr << "to_vary is the parameter to vary. 0: H, 1: z, 2: alpha_s, 3: quorum threshold" << endl;
    exit(8);
  }

  init_genrand(time(NULL));
  int Na = 100; // # ants
  int to_vary = atoi(argv[1]); // vary H if to_vary=0, vary z if to_vary=1, vary alpha_s if to_vary=2
  double alpha_s; // rate at which high-threshold ants visiting the poor nest move to the good nest
  double z; // initial fraction of recruiters
  double H; // fraction of high-threshold ants
  int samples; // # values for the varied var
  if (to_vary==0) { // vary H
    z = 0.3;
    alpha_s = 1; // alpha_s = 0.01, 0.1, 1
    samples = 6;
    cerr << "alpha_s = " << alpha_s << endl;
  } else if (to_vary==1) { // vary z
    H = 0.2;
    alpha_s = 1; // alpha_s = 0.01, 0.1, 1
    samples = 4;
    cerr << "alpha_s = " << alpha_s << endl;
  } else if (to_vary==2) { // vary alpha_s
    H = 0.2;
    z = 0.3; // z = 0.1 or 0.3
    samples = 6;
  } else if (to_vary==3) { // vary quorum threshold
    H = 0.2;
    alpha_s = 1; // alpha_s = 0.01, 0.1, 1
    samples = 6; // quorum threshold = 0.25, 0.3, 0.35, 0.4, 0.45, 0.5
  } else {
    cerr << "to_vary must be 0, 1, or 2" << endl;
    exit(8);
  }
  double L = 1-H; // fraction of low-threshold ants

  int tmp_sum;
  double corr; // Pearson correlation coefficient
  double tmp[5];

  double th_quorum_frac;
  int th_quorum; // quorum threshold
  double t;
  
  int i;
  double accum_rate[15], ra;

  int trials = 10000;
  int tr;

  int L_oldnest, H_oldnest, L_poor_com, L_poor_rec, H_poor_vis, L_good_com, H_good_com, L_good_rec, H_good_rec;

  int samples_alpha = 15;
  int samples_y;
  if (to_vary==0 || to_vary==1 || to_vary==2)
    samples_y = 6; // vary quorum threshold on the y-axis
  else
    samples_y = 4; // vary z on the y-axis
  int ind, ind_alpha, ind_y;
  double alpha, alpha_g, alpha_p; // conversion rate from committed to recruit
  double alpha_leak = 0.05; // = 0.0 in Fig. S1
  double t_quorum[samples];
  double precision[samples]; // % correct collective decision

  for (ind_alpha = -1 ; ind_alpha < samples_alpha ; ind_alpha++) { // -1: dummy
    alpha = 0.1 * (ind_alpha+1);
    alpha_g = alpha;
    alpha_p = alpha;
    for (ind_y = -1 ; ind_y < samples_y ; ind_y++) { // -1: dummy
      if (to_vary==0 || to_vary==1 || to_vary==2) {
	th_quorum_frac = 0.25 + 0.05 * ind_y;
	th_quorum = (int)(th_quorum_frac*Na); // quorum threshold
      } else { // to_vary==3
	z = 0.1*(ind_y+1);
      }

      if (ind_alpha == -1 || ind_y == -1) { // dummy
	corr = -10.0;
      } else { // really calculate the corr coef

  for (ind = 0 ; ind < samples ; ind++) {
    if (to_vary==0)
      H = (double)(ind+1)/15+1e-8;
    else if (to_vary==1)
      z = 0.1*(ind+1);
    else if (to_vary==2)
      alpha_s = 0.1 * ind;
    else if (to_vary==3) {
      th_quorum_frac = 0.25 + 0.05 * ind;
      th_quorum = (int)(th_quorum_frac*Na); // quorum threshold
    }

    t_quorum[ind] = precision[ind] = 0.0;

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
    if (tr==0 && Na - tmp_sum > 0)
      cerr << "Na - tmp_sum = " << Na - tmp_sum << "; " << H_oldnest << " " << L_oldnest << " " << H_poor_vis << " " << L_poor_com << " " << H_good_com << " " << L_good_com << endl;

    while (L_oldnest + H_oldnest < Na && L_good_com + H_good_com + L_good_rec + H_good_rec < th_quorum && L_poor_com + L_poor_rec + H_poor_vis < th_quorum) {

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
      tr++;
      t_quorum[ind] += t;
      if (L_good_com + H_good_com + L_good_rec + H_good_rec == th_quorum)
	precision[ind] += 1.0;
    }
    } // all trials done
      
      t_quorum[ind] /= trials;
      precision[ind] /= trials;
      //            cerr << z << " " << th_quorum_frac << " " << t_quorum[ind] << " " << precision[ind] << endl;
  } // all samples necessary for calculating the corr coeff done

  // calculate Pearson's correlation coefficient
  for (i=0 ; i<5 ; i++) tmp[i] = 0.0;
  for (ind=0 ; ind < samples ; ind++) {
    tmp[0] += t_quorum[ind];
    tmp[1] += precision[ind];
    tmp[2] += t_quorum[ind] * t_quorum[ind];
    tmp[3] += precision[ind] * precision[ind];
    tmp[4] += t_quorum[ind] * precision[ind];
  }
  tmp[0] /= samples;
  tmp[1] /= samples;
  //  for (i=0 ; i<5 ; i++)
  //    cerr << " " << tmp[i];
  //  cerr << endl;
  corr = (tmp[4]/samples - tmp[0]*tmp[1]) / sqrt(tmp[2]/samples - tmp[0]*tmp[0]) / sqrt(tmp[3]/samples - tmp[1]*tmp[1]);
      }

  cout << alpha << " " << ((to_vary < 3)? th_quorum_frac : z) << " " << corr <<  " " << alpha + 0.05 << " " <<  ((to_vary < 3)? th_quorum_frac + 0.025 : z + 0.05) << endl;
    } // all y samples (qthresh if to_vary==0, 1, or 2; z if to_vary==3) done
    cout << endl;
  } // all alpha samples done

  return 0;
}
