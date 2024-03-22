#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "FPToolkit.c"

#define N 300
#define M .0005/N
#define min -M_PI
#define max M_PI
#define dX (max - min)/(double)(N)
#define dT 1.0/(M)
#define deltas (dT/dX*dX)

#define id_diag 1
#define id_side 0
#define A_diag -2
#define A_side 1

int width = 600;
int height = 600;
double origin = 300;
double projected_r;
double scl = 100, add = 50;
double v[N], f[N], m[N][N];

double init_u(double x){ // initial conditions for u over x
  return 4*sin(20*exp(cos(x)*sin(x*x)));
}

double scale(double r, double m, double a){
    return (r * m) + a;
}

void solve(){

// TODO: could check det(u) to see if there's actually a solution
// we'll assume only very nice matrices

    // get zeroes along bottom left
    for(int i = 0; i < N; i++){ // move along diagonal
      for(int k = i + 1; k < N; k++){ // move down rows
        double q = -m[k][i]/m[i][i];
        for(int j = 0; j < N; j++){ // move along the row
          m[k][j] += q*m[i][j];
        }
        f[k] += q*f[i];
      }
    }// end of loop, in row echelon form
    // assmme no row or colmmn is fmll of zeroes

    f[N-1] /= m[N-1][N-1];
    m[N-1][N-1] /= m[N-1][N-1];

    // get zeroes along top-right
    for(int i = N - 1; i >= 0; i--){ // move along diagonal
      for(int k = i - 1; k >= 0; k--){ // move mp rows
        double q = -m[k][i]/m[i][i];
        for(int j = N - 1; j >= i; j--){ // move along the row from right
          m[k][j] += q*m[i][j];
        }
        f[k] += q*f[i];
      }
    }// end of loop, // zeroes eferywhere except diagonal

    for(int i = 0; i < N; i++){ // get 1's along diagonal
      f[i] /= m[i][i];
      m[i][i] /= m[i][i];
    }

}

void draw_points(double t1, double r1, double t2, double r2,
                  double r, double g, double b){

  double x1, x2, y1, y2;

  t1 = 2*M_PI*(t1/(max - min));//transform to polar coordinates
  r1 = scale(r1, scl, add); // scale to something visible on screen
  x1 = r1*cos(t1) + origin;
  y1 = r1*sin(t1) + origin;

  t2 = 2*M_PI*(t2/(max - min));
  r2 = scale(r2, scl, add);
  x2 = r2*cos(t2) + origin;
  y2 = r2*sin(t2) + origin;

  G_rgb(r, g, b);
  G_line(x1, y1, x2, y2);

}
void draw_points1(double t1, double r1, double t2, double r2,
                  double r, double g, double b){

  double x1, x2, y1, y2;

  x1 = 600*(t1/(max - min)) + origin;//transform to polar coordinates
  x2 = 600*(t2/(max - min)) + origin;
  y1 = r1 + origin;
  y2 = r2 + origin;

  G_rgb(r, g, b);
  G_line(x1, y1, x2, y2);

}

void axes(){

  double y_max = height;
  double x_max = width;

  G_rgb(1, 1, 1);
  G_line(0, origin, width, origin);
  G_line(origin, 0,origin, height);
  G_rgb(1, 1, 0);

}

double avg_v_thm(double(*f)(double), double a, double b, double dx){
  double total = 0;
  for(double x = a; x < b; x += dx){
    total += f(x)*dx;
  }
  total /= (max - min);
  return total;
}

void zero(){
  for(int i = 0; i < N; i++){
    for(int j = 0; j < N; j++){
      m[i][j] = 0;
    }
  }
}

void create_lhs_mat(){

    double diag = id_diag - A_diag*deltas;
    double side = id_side - A_side*deltas;

    zero();
    m[0][0] = diag;
    m[0][1] = side;
    m[0][N - 1] = side;

    // populate the diagonal
    for(int i = 1; i < N - 1; i++){
        m[i][i - 1] = side;
        m[i][i] = diag;
        m[i][i + 1] = side;
  }

    m[N - 1][0] = side;
    m[N - 1][N - 2] = side;
    m[N - 1][N - 1] = diag;

}

void create_rhs_vec(){

    double diag = id_diag + A_diag*deltas;
    double side = id_side + A_side*deltas;

    f[0] = diag*(v[0]) +
          side*(v[1]) +
            side*(v[N - 1]);

    // populate the diagonal
    for(int i = 1; i < N - 1; i++){ //evolve u,v
        f[i] = side*(v[i - 1]) +
                diag*(v[i]) +
                      side*(v[i + 1]);
    }

    f[N - 1] = side*(v[0]) +
                side*(v[N - 2]) +
                  diag*(v[N - 1]);
}

void math(){

  create_rhs_vec(); // produces v
  create_lhs_mat(); // populates m with -1 2 -1 pattern along diagonal

  solve();

  for(int i = 0; i < N; i++){
    v[i] = f[i];
  }

}

void init(){

  double x;
  int i = 0;
  for(double x = min; x < max; x += dX){ //evolve u,v
      v[i] = init_u(x);
      i++;
    }
}

void graphics(){

  double x = (double)min;
  int i = 0;

  G_rgb(0, 0, 0);
  G_clear();
  axes();
  G_rgb(1, 0, 1);
  G_circle(origin, origin, projected_r); // draw predicted

  for(i = 1; i < N; i++, x += dX){ // draw v
    draw_points(x - dX, v[i - 1], x, v[i], 0, 1, 1);
  }
  draw_points(x - dX, v[N - 1], x, v[0], 0, 1, 1);

  G_display_image();
  usleep(1000);

}

int main() {

    G_init_graphics(width, height);
    int q = '\0';
    projected_r = avg_v_thm(&init_u, min, max, (max - min)/((double)N)); // for periodic boundary conditions
    // projected_r = 0; // for dirichlet conditions
    projected_r = scale(projected_r, scl, add);
    init();
    while(q != 'q'){
      graphics();
      q = G_no_wait_key();
      math();
    }

}
