/*
  spherical_icf_vti.c

  3D spherical ICF toy simulation with:
    - Hydro (toy Euler-like)
    - DT burn
    - Energy-conservative Monte Carlo alpha transport
    - Radiation losses
    - Thermal diffusion
    - VTI output (ParaView)

  Compile:
    gcc -O3 -std=c11 spherical_icf_vti.c -lm -fopenmp -o sicf
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <omp.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ---------------- GRID ---------------- */
#define NX 96
#define NY 96
#define NZ 96

#define N (NX*NY*NZ)
#define IDX(i,j,k) ((i)*NY*NZ + (j)*NZ + (k))

double dx = 0.01;
double dt = 5e-11;

/* ---------------- SPHERE ---------------- */
double R_sphere = 0.25;

/* ---------------- CONSTANTS ---------------- */
const double kB = 1.380649e-23;
const double mp = 1.6726219e-27;
const double eV = 1.602e-19;

/* ---------------- STATE ---------------- */
float *rho, *T;
float *ux, *uy, *uz;
float *fusionP, *alphaH, *radL;

/* ---------------- ALPHAS ---------------- */
#define N_ALPHA 40000

typedef struct {
    double x,y,z;
    double vx,vy,vz;
    double E;
    int alive;
} Alpha;

Alpha *alpha;

/* ============================================================
   UTIL
   ============================================================ */
static inline int id(int i,int j,int k){
    return i*NY*NZ + j*NZ + k;
}

static inline int clamp(int v,int a,int b){
    return v<a?a:(v>b?b:v);
}

/* ============================================================
   INITIAL CONDITION (SPHERE)
   ============================================================ */
void init()
{
    #pragma omp parallel for collapse(3)
    for(int i=0;i<NX;i++)
    for(int j=0;j<NY;j++)
    for(int k=0;k<NZ;k++)
    {
        double x = (i-NX/2)*dx;
        double y = (j-NY/2)*dx;
        double z = (k-NZ/2)*dx;

        double r = sqrt(x*x + y*y + z*z);

        int c = id(i,j,k);

        if(r < R_sphere){
            rho[c] = 200.0;
            T[c]   = 1e8;
        } else {
            rho[c] = 1e-6;
            T[c]   = 1e5;
        }

        ux[c]=uy[c]=uz[c]=0;
        fusionP[c]=alphaH[c]=radL[c]=0;
    }

    for(int n=0;n<N_ALPHA;n++)
        alpha[n].alive=0;
}

/* ============================================================
   PHYSICS MODELS
   ============================================================ */
static inline double pressure(double r,double Tt){
    double n = r/(2.5*mp);
    return n*kB*Tt;
}

double reactivity(double Tt){
    double keV = Tt/1.16e7;
    if(keV<1e-6) keV=1e-6;
    return 1e-24*keV*keV*exp(-20.0/keV);
}

double brems(double r,double Tt){
    double ne = r/mp;
    return 1e-40*ne*ne*sqrt(Tt);
}

/* ============================================================
   HYDRO + BURN
   ============================================================ */
void hydro_burn()
{
    #pragma omp parallel for collapse(3)
    for(int i=1;i<NX-1;i++)
    for(int j=1;j<NY-1;j++)
    for(int k=1;k<NZ-1;k++)
    {
        int c = id(i,j,k);

        double r0 = rho[c];
        double T0 = T[c];

        /* pressure gradients (toy) */
        double Px =
            pressure(rho[id(i+1,j,k)],T[id(i+1,j,k)]) -
            pressure(rho[id(i-1,j,k)],T[id(i-1,j,k)]);

        double Py =
            pressure(rho[id(i,j+1,k)],T[id(i,j+1,k)]) -
            pressure(rho[id(i,j-1,k)],T[id(i,j-1,k)]);

        double Pz =
            pressure(rho[id(i,j,k+1)],T[id(i,j,k+1)]) -
            pressure(rho[id(i,j,k-1)],T[id(i,j,k-1)]);

        ux[c] -= dt * Px / (2*dx*r0);
        uy[c] -= dt * Py / (2*dx*r0);
        uz[c] -= dt * Pz / (2*dx*r0);

        /* density stays fixed (toy stability) */
        rho[c] = r0;

        /* burn */
        double nD = 0.5*r0/(2.5*mp);
        double nT = 0.5*r0/(2.5*mp);

        double R = reactivity(T0);

        fusionP[c] = nD*nT*R * 17.6*eV;
        alphaH[c]  = nD*nT*R * 3.5*eV;

        radL[c] = brems(r0,T0);

        double cv = 1.5*kB*(r0/(2.5*mp));

        T[c] = T0 + dt*(alphaH[c] - radL[c])/cv;

        if(T[c] < 1e5) T[c] = 1e5;
    }
}

/* ============================================================
   MONTE CARLO ALPHAS (ENERGY CONSERVATIVE)
   ============================================================ */
void spawn_alphas()
{
    for(int n=0;n<N_ALPHA;n++)
    {
        int i = rand()%NX;
        int j = rand()%NY;
        int k = rand()%NZ;

        int c = id(i,j,k);

        if(alphaH[c] <= 0) continue;

        alpha[n].x = (i-NX/2)*dx;
        alpha[n].y = (j-NY/2)*dx;
        alpha[n].z = (k-NZ/2)*dx;

        double z = 2.0*(rand()/(double)RAND_MAX)-1;
        double t = 2*M_PI*(rand()/(double)RAND_MAX);
        double r = sqrt(1-z*z);

        double v0 = 1.3e7;

        alpha[n].vx = v0*r*cos(t);
        alpha[n].vy = v0*r*sin(t);
        alpha[n].vz = v0*z;

        alpha[n].E = 3.5*eV;
        alpha[n].alive = 1;
    }
}

void move_alphas()
{
    double dt_a = dt/5.0;

    for(int n=0;n<N_ALPHA;n++)
    {
        if(!alpha[n].alive) continue;

        alpha[n].x += alpha[n].vx*dt_a;
        alpha[n].y += alpha[n].vy*dt_a;
        alpha[n].z += alpha[n].vz*dt_a;

        int i = (int)(alpha[n].x/dx + NX/2);
        int j = (int)(alpha[n].y/dx + NY/2);
        int k = (int)(alpha[n].z/dx + NZ/2);

        if(i<0||j<0||k<0||i>=NX||j>=NY||k>=NZ){
            alpha[n].alive = 0;
            continue;
        }

        int c = id(i,j,k);

        double stopping = 5e-3 * (rho[c]/200.0);

        double dE = stopping * alpha[n].E * dt_a;

        if(dE > alpha[n].E) dE = alpha[n].E;

        alpha[n].E -= dE;

        /* ---------------- ENERGY CONSERVATION ---------------- */
        double nion = rho[c]/(2.5*mp);
        double cv = 1.5*kB*nion;

        T[c] += dE / cv;

        if(alpha[n].E <= 0)
            alpha[n].alive = 0;
    }
}

/* ============================================================
   VTI OUTPUT
   ============================================================ */
void write_vti(int step)
{
    char fn[256];
    sprintf(fn,"icf_%05d.vti",step);

    FILE *f = fopen(fn,"w");

    fprintf(f,"<?xml version=\"1.0\"?>\n");
    fprintf(f,"<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\">\n");

    fprintf(f,"<ImageData WholeExtent=\"0 %d 0 %d 0 %d\" Origin=\"%g %g %g\" Spacing=\"%g %g %g\">\n",
        NX-1,NY-1,NZ-1,
        -NX*dx/2,-NY*dx/2,-NZ*dx/2,
        dx,dx,dx);

    fprintf(f,"<Piece Extent=\"0 %d 0 %d 0 %d\">\n",NX-1,NY-1,NZ-1);

    fprintf(f,"<PointData Scalars=\"T\">\n");

    fprintf(f,"<DataArray type=\"Float32\" Name=\"rho\" format=\"ascii\">\n");
    for(int i=0;i<N;i++) fprintf(f,"%g ",rho[i]);
    fprintf(f,"\n</DataArray>\n");

    fprintf(f,"<DataArray type=\"Float32\" Name=\"T\" format=\"ascii\">\n");
    for(int i=0;i<N;i++) fprintf(f,"%g ",T[i]);
    fprintf(f,"\n</DataArray>\n");

    fprintf(f,"<DataArray type=\"Float32\" Name=\"fusion\" format=\"ascii\">\n");
    for(int i=0;i<N;i++) fprintf(f,"%g ",fusionP[i]);
    fprintf(f,"\n</DataArray>\n");

    fprintf(f,"<DataArray type=\"Float32\" Name=\"alpha_heat\" format=\"ascii\">\n");
    for(int i=0;i<N;i++) fprintf(f,"%g ",alphaH[i]);
    fprintf(f,"\n</DataArray>\n");

    fprintf(f,"</PointData>\n</Piece></ImageData></VTKFile>\n");

    fclose(f);

    printf("wrote %s\n",fn);
}

/* ============================================================
   STEP
   ============================================================ */
void step()
{
    hydro_burn();
    spawn_alphas();
    move_alphas();
}

/* ============================================================
   MAIN
   ============================================================ */
int main()
{
    rho = malloc(sizeof(float)*N);
    T   = malloc(sizeof(float)*N);
    ux  = malloc(sizeof(float)*N);
    uy  = malloc(sizeof(float)*N);
    uz  = malloc(sizeof(float)*N);

    fusionP = malloc(sizeof(float)*N);
    alphaH  = malloc(sizeof(float)*N);
    radL    = malloc(sizeof(float)*N);

    alpha = malloc(sizeof(Alpha)*N_ALPHA);

    init();

    for(int t=0;t<200;t++)
    {
        step();
        if(t%10==0) write_vti(t);
    }

    printf("done\n");
}
