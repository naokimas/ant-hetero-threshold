/* Differential equation version of the netst choice model

Produce the results shown in Figure 3 by running
        a.out of diffeqn-nestchoice.out alpha=0.1 alpha_s=0.1 H=0.2 z=0.3
*/

#include <iostream>
using namespace std;
#include <fstream>
#include <cstdlib> // atoi
#include <cmath> // cos, sin
#include <ctime>
#define VERBOSE

int main (int argc, char **argv) {

  if (argc != 5) {
    cerr << "diffeqn-nestchoice.out alpha alpha_s H z" << endl;
    cerr << "alpha: rate at which committed ants convert to recruits, common to low-threshold and high-threshold ants" << endl;
    cerr << "alpha_s: rate at which high-threshold ants visiting the poor nest moves to the good nest" << endl;
    cerr << "H: fraction of high-threshold ants" << endl;
     cerr << "z: initial fraction of recruiters" << endl;
    exit(8);
  }

  double alpha = atof(argv[1]); // conversion rate from committed to recruit
  double alpha_g = alpha;
  double alpha_p = alpha;
  double alpha_s = atof(argv[2]); // rate at which high-threshold ants visiting the poor nest move to the good nest
  double alpha_leak = 0.05;
  double H = atof(argv[3]); // fraction of high-threshold ants
  double L = 1-H; // fraction of low-threshold ants
  double z = atof(argv[4]); // initial fraction of recruiters

  double t=0;
  double dt = 0.001;
  double t_quorum;
  double sum; // for checking

  double L_oldnest, H_oldnest, L_poor_com, L_poor_rec, H_poor_vis, L_good_com, H_good_com, L_good_rec, H_good_rec;
  double L_oldnest_prev, H_oldnest_prev, L_poor_com_prev, L_poor_rec_prev, H_poor_vis_prev, L_good_com_prev, H_good_com_prev, L_good_rec_prev, H_good_rec_prev;

  // initialization
  L_oldnest = L*(1-z);
  H_oldnest = H*(1-z);
  H_poor_vis = H_good_com = H*z/2;
  L_poor_com = L_good_com = L*z/2;
  L_poor_rec = L_good_rec = H_good_rec = 0.0;

  double eps = 0.1; // dynamics stop when "# ants in the current nest <= eps" is reached for the first time

  while (L_oldnest + H_oldnest > eps) {

    if ((int)(t*100) > (int)((t-dt)*100))
      cout << t << " " << L_good_rec + H_good_rec << " " << L_good_com + H_good_com + L_good_rec + H_good_rec << " " << L_poor_rec << " " << L_poor_com + L_poor_rec + H_poor_vis << endl;

    if (L_good_com + H_good_com + L_good_rec + H_good_rec > 0.5 
	&& L_good_com_prev + H_good_com_prev + L_good_rec_prev + H_good_rec_prev < 0.5)
      t_quorum = t;
    L_oldnest_prev = L_oldnest;
    H_oldnest_prev = H_oldnest;
    L_poor_com_prev = L_poor_com;
    L_poor_rec_prev = L_poor_rec;
    H_poor_vis_prev = H_poor_vis;
    L_good_com_prev = L_good_com;
    H_good_com_prev = H_good_com;
    L_good_rec_prev = L_good_rec;
    H_good_rec_prev = H_good_rec;

    L_oldnest += dt * (alpha_leak * (L_poor_com_prev + L_poor_rec_prev + L_good_com_prev + L_good_rec_prev) -
       (L_poor_rec_prev + L_good_rec_prev + H_good_rec_prev) * L_oldnest_prev);
    H_oldnest += dt * (alpha_leak * (H_poor_vis_prev + H_good_com_prev + H_good_rec_prev) -
       (L_poor_rec_prev + L_good_rec_prev + H_good_rec_prev) * H_oldnest_prev);
    L_poor_com += dt * (L_poor_rec_prev * L_oldnest_prev - alpha_p * L_poor_com_prev - alpha_leak * L_poor_com_prev);
    H_poor_vis += dt * (L_poor_rec_prev * H_oldnest_prev - alpha_s * H_poor_vis_prev - alpha_leak * H_poor_vis_prev);
    L_good_com += dt * ((L_good_rec_prev + H_good_rec_prev) * L_oldnest_prev - alpha_g * L_good_com_prev - alpha_leak * L_good_com_prev);
    H_good_com += dt * ((L_good_rec_prev + H_good_rec_prev) * H_oldnest_prev + alpha_s * H_poor_vis_prev - alpha_g * H_good_com_prev - alpha_leak * H_good_com_prev);
    L_poor_rec += dt * (alpha_p * L_poor_com_prev - alpha_leak * L_poor_rec_prev);
    L_good_rec += dt * (alpha_g * L_good_com_prev - alpha_leak * L_good_rec_prev);
    H_good_rec += dt * (alpha_g * H_good_com_prev - alpha_leak * H_good_rec_prev);

    sum = L_oldnest + H_oldnest + L_poor_com + H_poor_vis + L_good_com + H_good_com + L_poor_rec + L_good_rec + H_good_rec;
    //    cerr << sum << endl;
    L_oldnest /= sum;
    H_oldnest /= sum;
    L_poor_com /= sum;
    H_poor_vis /= sum;
    L_good_com /= sum;
    H_good_com /= sum;
    L_poor_rec /= sum;
    L_good_rec /= sum;
    H_good_rec /= sum;
    t += dt;
  }

  return 0;
}
