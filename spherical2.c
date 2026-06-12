#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#define NX 96
#define NY 96
#define NZ 96
#define N (NX*NY*NZ)

static inline int id(int i,int j,int k){ return i*NY*NZ + j*NZ + k; }

/* ---------------- physics ---------------- */
const double gamma_g = 5.0/3.0;
const double kB = 1.380649e-23;
const double mp = 1.6726219e-27;

double dx = 0.01;
double dt = 2e-12;

/* ---------------- fields ---------------- */
float *rho, *ux, *uy, *uz, *E;
float *rho_n, *ux_n, *uy_n, *uz_n, *E_n;
float *fusionF, *alphaH, *radL;

/* ---------------- kernel for alpha deposition ---------------- */
#define K 7
double alpha_kernel[K][K][K];

/* ============================================================
   EOS
   ============================================================ */
static inline double pressure(double r,double E,double u2){
    return (gamma_g-1.0)*(E - 0.5*r*u2);
}

/* ============================================================
   fusion + radiation
   ============================================================ */
double reactivity(double T){
    double keV = T/1.16e7;
    if(keV < 1e-6) keV = 1e-6;
    return 1e-24 * keV*keV * exp(-20.0/keV);
}

double brems(double r,double T){
    double ne = r/mp;
    return 1e-40 * ne*ne * sqrt(T);
}

/* ============================================================
   alpha deposition kernel (smooth, energy-conserving)
   ============================================================ */
void init_kernel()
{
    double sigma = 1.5;
    double sum = 0.0;

    for(int i=0;i<K;i++)
    for(int j=0;j<K;j++)
    for(int k=0;k<K;k++){
        double x = i - (K-1)/2;
        double y = j - (K-1)/2;
        double z = k - (K-1)/2;
        double g = exp(-(x*x + y*y + z*z)/(2*sigma*sigma));
        alpha_kernel[i][j][k] = g;
        sum += g;
    }

    for(int i=0;i<K;i++)
    for(int j=0;j<K;j++)
    for(int k=0;k<K;k++)
        alpha_kernel[i][j][k] /= sum;

    printf("alpha kernel initialized (smooth, energy-conserving)\n");
}

/* ============================================================
   INIT: smooth spherical target + laser drive
   ============================================================ */
void init()
{
    rho = malloc(sizeof(float)*N);
    ux  = malloc(sizeof(float)*N);
    uy  = malloc(sizeof(float)*N);
    uz  = malloc(sizeof(float)*N);
    E   = malloc(sizeof(float)*N);

    rho_n = malloc(sizeof(float)*N);
    ux_n  = malloc(sizeof(float)*N);
    uy_n  = malloc(sizeof(float)*N);
    uz_n  = malloc(sizeof(float)*N);
    E_n   = malloc(sizeof(float)*N);

    fusionF = malloc(sizeof(float)*N);
    alphaH  = malloc(sizeof(float)*N);
    radL    = malloc(sizeof(float)*N);

    double R0 = 0.25;        // shell radius
    double thickness = 0.05; // shell width
    double rho0 = 200.0;     // peak density
    double T0 = 1e8;         // initial temperature

    #pragma omp parallel for collapse(3)
    for(int i=0;i<NX;i++)
    for(int j=0;j<NY;j++)
    for(int k=0;k<NZ;k++)
    {
        int c = id(i,j,k);
        double x = (i-(NX-1)/2.0)*dx;
        double y = (j-(NY-1)/2.0)*dx;
        double z = (k-(NZ-1)/2.0)*dx;
        double r = sqrt(x*x+y*y+z*z);

        /* smooth shell using tanh SDF */
        double sdf = r - R0;
        double h = 0.5*(1 - tanh(sdf/thickness));

        rho[c] = rho0 * h;
        ux[c] = uy[c] = uz[c] = 0.0;
        double p = (rho[c]/(2.5*mp))*kB*T0;
        E[c] = p/(gamma_g-1.0);

        fusionF[c] = alphaH[c] = radL[c] = 0.0;
    }
}

/* ============================================================
   Rusanov flux (simplified 3D Euler)
   ============================================================ */
void step()
{
    #pragma omp parallel for collapse(3)
    for(int i=1;i<NX-1;i++)
    for(int j=1;j<NY-1;j++)
    for(int k=1;k<NZ-1;k++)
    {
        int c = id(i,j,k);

        double r = rho[c];
        double u = ux[c];
        double v = uy[c];
        double w = uz[c];
        double E0 = E[c];

        double u2 = u*u+v*v+w*w;
        double P = pressure(r,E0,u2);
        double cs = sqrt(gamma_g*P/r);
        double a = fabs(u)+cs;

        /* density */
        rho_n[c] = r - dt * r *
        ((ux[id(i+1,j,k)]-ux[id(i-1,j,k)]) +
         (uy[id(i,j+1,k)]-uy[id(i,j-1,k)]) +
         (uz[id(i,j,k+1)]-uz[id(i,j,k-1)])) / (2*dx);

        /* momentum */
        ux_n[c] = u - dt * (P - pressure(rho[id(i-1,j,k)],E[id(i-1,j,k)],0)) / (r*dx);
        uy_n[c] = v;
        uz_n[c] = w;

        /* energy flux */
        E_n[c] = E0 - dt * a * (E[id(i+1,j,k)] - E[id(i-1,j,k)])/(2*dx);

        /* fusion + radiation */
        double T = P/(r/(2.5*mp))/kB;
        double R = reactivity(T);
        double n = r/mp;
        double fusion = n*n*R * 17.6e-19;
        double alpha  = n*n*R * 3.5e-19;
        double rad    = brems(r,T);
        fusionF[c]=fusion;
        radL[c]=rad;

        /* NONLOCAL alpha deposition */
        int ci=i,cj=j,ck=k;
        for(int di=0;di<K;di++)
        for(int dj=0;dj<K;dj++)
        for(int dk=0;dk<K;dk++){
            int ii = ci + di - (K-1)/2;
            int jj = cj + dj - (K-1)/2;
            int kk = ck + dk - (K-1)/2;
            if(ii<1||jj<1||kk<1||ii>=NX-1||jj>=NY-1||kk>=NZ-1) continue;
            int nidx = id(ii,jj,kk);
            E_n[nidx] += dt * alpha * alpha_kernel[di][dj][dk];
        }
    }

    /* swap */
    float *tmp;
    tmp=rho; rho=rho_n; rho_n=tmp;
    tmp=ux; ux=ux_n; ux_n=tmp;
    tmp=uy; uy=uy_n; uy_n=tmp;
    tmp=uz; uz=uz_n; uz_n=tmp;
    tmp=E;  E=E_n;   E_n=tmp;
}

/* ============================================================
   Swiss-cheesed VTI writer
   ============================================================ */
void write_vti(int stepnum)
{
    char fn[256];
    sprintf(fn,"icf_ignite_%05d.vti",stepnum);
    FILE *f=fopen(fn,"w");
    double ox=-((NX-1)/2.0)*dx;

    fprintf(f,"<?xml version=\"1.0\"?>\n");
    fprintf(f,"<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
    fprintf(f,"<ImageData WholeExtent=\"0 %d 0 %d 0 %d\" Origin=\"%g %g %g\" Spacing=\"%g %g %g\">\n",
            NX-1,NY-1,NZ-1,ox,ox,ox,dx,dx,dx);
    fprintf(f,"<Piece Extent=\"0 %d 0 %d 0 %d\">\n",NX-1,NY-1,NZ-1);
    fprintf(f,"<PointData Scalars=\"rho\">\n");

#define W(name,arr) \
    fprintf(f,"<DataArray type=\"Float32\" Name=\"" name "\" format=\"ascii\">\n"); \
    for(int k=0;k<NZ;k++) \
    for(int j=0;j<NY;j++) \
    for(int i=0;i<NX;i++) \
        fprintf(f,"%g ",arr[id(i,j,k)]); \
    fprintf(f,"\n</DataArray>\n")

    W("rho",rho);
    W("ux",ux);
    W("uy",uy);
    W("uz",uz);
    W("E",E);
    W("fusion",fusionF);
    W("alpha",alphaH);
    W("radiation",radL);

#undef W

    fprintf(f,"</PointData>\n");
    fprintf(f,"</Piece>\n");
    fprintf(f,"</ImageData>\n");
    fprintf(f,"</VTKFile>\n");
    fclose(f);

    printf("wrote %s\n",fn);
}

/* ============================================================
   MAIN
   ============================================================ */
int main()
{
    init_kernel();
    init();

    for(int t=0;t<200;t++){
        step();
        if(t%10==0) write_vti(t);
    }

    printf("done\n");
    return 0;
}
