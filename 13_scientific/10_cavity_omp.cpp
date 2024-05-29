#include <stdio.h>
#include <omp.h>

int main() {
    omp_set_num_threads(4);
    int nx = 41;
    int ny = 41;
    int nt = 500;
    int nit = 50;
    double dx;
    double dy;
    double dt = 0.01;
    double rho = 1.0;
    double nu = 0.02;
    int i = 0;
    int j = 0;
    int n;
    int it;
    dx = 2.0 / (nx - 1);  // dx, dyが0にならないように2を2.0とした．
    dy = 2.0 / (ny - 1);

    double u[ny][nx];
    double un[ny][nx];
    double v[ny][nx];
    double vn[ny][nx];
    double p[ny][nx];
    double pn[ny][nx];
    double b[ny][nx];

    // initializing
#pragma omp parallel for private(j,i)
    for(j = 0; j < ny; j++){
        for(i = 0; i < nx; i++){
            u[j][i] = 0;
            v[j][i] = 0;
            p[j][i] = 0;
            b[j][i] = 0;
        }
    }

    // time evolution
    for(n = 0; n < nt; n++){
        // solve b
#pragma omp parallel for private(j,i)
        for(j = 1; j < ny-1; j++){
            for(i = 1; i < nx-1; i++){   
                b[j][i] = rho * (1 / dt *
                    ((u[j][i+1] - u[j][i-1]) / (2 * dx) + (v[j+1][i] - v[j-1][i]) / (2 * dy)) -
                    ((u[j][i+1] - u[j][i-1]) / (2 * dx)) * ((u[j][i+1] - u[j][i-1]) / (2 * dx)) - 2 * ((u[j+1][i] - u[j-1][i]) / (2 * dy) *
                     (v[j][i+1] - v[j][i-1]) / (2 * dx)) - ((v[j+1][i] - v[j-1][i]) / (2 * dy)) * ((v[j+1][i] - v[j-1][i]) / (2 * dy)));
            }
        }

        // solve Poisson eq (p) iterately
        for(it = 0; it < nit; it++){
            // copy p to pn
#pragma omp parallel for private(j,i)
            for(j = 0; j < ny; j++){
                for(i = 0; i < nx; i++){
                    pn[j][i] = p[j][i];
                }
            }            

            // compute p
#pragma omp parallel for private(j,i)
            for(j = 1; j < ny-1; j++){
                for(i = 1; i < nx-1; i++){
                    p[j][i] = (dy*dy * (pn[j][i+1] + pn[j][i-1]) +
                           dx*dx * (pn[j+1][i] + pn[j-1][i]) -
                           b[j][i] * dx*dx * dy*dy)
                          / (2 * (dx*dx + dy*dy));
                }
            }

            // boundary condition of p
#pragma omp parallel for
            for(j = 0; j < ny; j++){
		p[j][nx-1] = p[j][nx-2];
                p[j][0] = p[j][1];
            }
#pragma omp parallel for
            for(i = 0; i < nx; i++){
                p[0][i] = p[1][i];
		p[ny-1][i] = 0;
            }
        }
        


        // solve NS wq (u, v)
        // copy u, v to un, vn
#pragma omp parallel for private(j,i)
        for(j = 0; j < ny; j++){
            for(i = 0; i < nx; i++){
                un[j][i] = u[j][i];
                vn[j][i] = v[j][i];
            }
        }

        // compute u, v
#pragma omp parallel for private(j,i)
        for(j = 1; j < ny-1; j++){
            for(i = 1; i < nx-1; i++){
                u[j][i] = un[j][i] - un[j][i] * dt / dx * (un[j][i] - un[j][i - 1])
                               - un[j][i] * dt / dy * (un[j][i] - un[j - 1][i])
                               - dt / (2 * rho * dx) * (p[j][i+1] - p[j][i-1])
                               + nu * dt / (dx*dx) * (un[j][i+1] - 2 * un[j][i] + un[j][i-1])
                               + nu * dt / (dy*dy) * (un[j+1][i] - 2 * un[j][i] + un[j-1][i]);

                v[j][i] = vn[j][i] - vn[j][i] * dt / dx * (vn[j][i] - vn[j][i - 1])
                               - vn[j][i] * dt / dy * (vn[j][i] - vn[j - 1][i])
                               - dt / (2 * rho * dx) * (p[j+1][i] - p[j-1][i])
                               + nu * dt / (dx*dx) * (vn[j][i+1] - 2 * vn[j][i] + vn[j][i-1])
                               + nu * dt / (dy*dy) * (vn[j+1][i] - 2 * vn[j][i] + vn[j-1][i]);
            }
        }

        // boundary condition of u, v
#pragma omp parallel for
        for(j = 0; j < ny; j++){
            u[j][0]  = 0;
	    u[j][nx-1] = 0;
            v[j][0]  = 0;
	    v[j][nx-1] = 0;
        }
#pragma omp parallel for
        for(i = 0; i < nx; i++){
            u[0][i]  = 0;
            u[ny-1][i] = 1;
            v[0][i]  = 0;
	    v[ny-1][i] = 0;
        }
         printf("%f\n", u[30][30]);
         printf("%f\n", v[30][30]);
    }

}
