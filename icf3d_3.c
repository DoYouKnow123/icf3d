/*
    icf3d_module1.c

    Educational 3D inertial confinement fusion framework.

    Build:
        gcc -O3 icf3d_module1.c -lm -o icf3d

    Output:
        snapshot_000.vti
        snapshot_001.vti
        snapshot_002.vti
        snapshot_003.vti
        snapshot_004.vti
        snapshot_005.vti
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define NX 128
#define NY 128
#define NZ 128

#define DX 1.0e-2

#define SNAPSHOT_INTERVAL 2.0e-9
#define TOTAL_TIME        4.0e-8
#define HEATING_FACTOR    1.0e-17
#define DT 1.0e-10
static double *newRho;
static double *newT;
static double *newVx;
static double *newVy;
static double *newVz;
typedef struct
{
    double rho;       /* density */
    double T;         /* temperature */
    double burn;      /* burn fraction */

    double vx;
    double vy;
    double vz;

    double pressure;
    double fusionPower;

} Cell;

static Cell *grid = NULL;

#define IDX(i,j,k) \
((k)*NY*NX + (j)*NX + (i))

static inline Cell *cell(int i,int j,int k)
{
    return &grid[IDX(i,j,k)];
}
double reactivity(double T)
{
    double TkeV = T / 1.160e7;

    if(TkeV < 0.1)
        return 0.0;

    return 1.0e-24 *
           pow(TkeV,2.0) *
           exp(-19.94/sqrt(TkeV));
}
double sigmav_DT(double T)
{
    const double Bg   = 34.3827;
    const double mr_c2 = 1124656.0;

    const double C1 = 1.17302e-9;
    const double C2 = 1.51361e-2;
    const double C3 = 7.51886e-2;
    const double C4 = 4.60643e-3;
    const double C5 = 1.35000e-2;
    const double C6 = -1.06750e-4;
    const double C7 = 1.36600e-5;

    double theta =
        T / (1.0
             - T * (C2 + T * (C4 + T * C6))
             / (1.0 + T * (C3 + T * (C5 + T * C7))));

    double xi = pow((Bg * Bg) / (4.0 * theta), 1.0 / 3.0);

    return C1 * theta *
           sqrt(xi / (mr_c2 * T * T * T)) *
           exp(-3.0 * xi);
}

/* rho from reaction rate density */
double rho_DT(double T, double R)
{
    const double mD = 3.34358372e-24; /* g */
    const double mT = 5.00735675e-24; /* g */

    double sigmav = sigmav_DT(T);

    return (mD + mT) * sqrt(R / sigmav);
}
double ideal_pressure(double rho,double T)
{
    const double R = 8.314e3;
    return rho * R * T;
}
double clamp(double x,double lo,double hi)
{
    if(x < lo) return lo;
    if(x > hi) return hi;
    return x;
}

double divergence_velocity(int i,int j,int k)
{
    int im = (i>0)    ? i-1 : i;
    int ip = (i<NX-1) ? i+1 : i;

    int jm = (j>0)    ? j-1 : j;
    int jp = (j<NY-1) ? j+1 : j;

    int km = (k>0)    ? k-1 : k;
    int kp = (k<NZ-1) ? k+1 : k;

    Cell *cxm = cell(im,j,k);
    Cell *cxp = cell(ip,j,k);

    Cell *cym = cell(i,jm,k);
    Cell *cyp = cell(i,jp,k);

    Cell *czm = cell(i,j,km);
    Cell *czp = cell(i,j,kp);

    double dudx =
        (cxp->vx - cxm->vx)/(2.0*DX);

    double dvdy =
        (cyp->vy - cym->vy)/(2.0*DX);

    double dwdz =
        (czp->vz - czm->vz)/(2.0*DX);

    return dudx + dvdy + dwdz;
}
void allocate_grid(void)
{
	size_t n = (size_t)NX*NY*NZ;

	newRho = calloc(n,sizeof(double));
	newT   = calloc(n,sizeof(double));

	newVx  = calloc(n,sizeof(double));
	newVy  = calloc(n,sizeof(double));
	newVz  = calloc(n,sizeof(double));

    grid = (Cell*)calloc(n,sizeof(Cell));

    if(!grid)
    {
        fprintf(stderr,"allocation failed\n");
        exit(EXIT_FAILURE);
    }
}

void initialize_capsule(void)
{

    double cx = NX*0.5;
    double cy = NY*0.5;
    double cz = NZ*0.5;

    double capsuleRadius = NX*0.18;
    double hotspotRadius = NX*0.04;

    for(int k=0;k<NZ;k++)
    {
        for(int j=0;j<NY;j++)
        {
            for(int i=0;i<NX;i++)
            {
                Cell *c = cell(i,j,k);

                double dx=i-cx;
                double dy=j-cy;
                double dz=k-cz;

                double r =
                    sqrt(dx*dx+dy*dy+dz*dz);

                if(r < capsuleRadius)
                {
                    c->rho = 250.0;
                    c->T   = 2.0e6;
                }
                else
                {
                    c->rho = 1.0;
                    c->T   = 1.0e5;
                }

                if(r < hotspotRadius)
                {
                    c->T = 8.0e6;
                }

                c->burn = 0.0;

                c->vx = 0.0;
                c->vy = 0.0;
                c->vz = 0.0;

	
                c->pressure =
                    ideal_pressure(
                        c->rho,
                        c->T);
		double rate =
			reactivity(c->T);

         //       c->rho = rho_DT(c->T,rate);

		double burnRate =
			c->rho * c->rho * rate;
		
		if(c->rho < 50.0)
		{
		burnRate *= c->rho/50.0;
		}
		c->burn += burnRate * DT;

		if(c->burn > 1.0)
    			c->burn = 1.0;
		c->fusionPower =
			burnRate * 2.82e-12;
		double alphaEnergy =
			burnRate * 5.6e-13;

		c->T += alphaEnergy * HEATING_FACTOR;
            }
        }
    }
}

void write_vti_snapshot(int frame)
{
    char filename[128];

    sprintf(filename,
            "snapshot_%03d.vti",
            frame);

    FILE *f = fopen(filename,"w");

    if(!f)
    {
        fprintf(stderr,
                "cannot write %s\n",
                filename);
        return;
    }

    fprintf(f,
"<?xml version=\"1.0\"?>\n"
"<VTKFile type=\"ImageData\" version=\"0.1\" "
"byte_order=\"LittleEndian\">\n");

    fprintf(f,
"<ImageData WholeExtent=\"0 %d 0 %d 0 %d\" "
"Origin=\"0 0 0\" "
"Spacing=\"%g %g %g\">\n",
            NX-1,
            NY-1,
            NZ-1,
            DX,DX,DX);

    fprintf(f,
"<Piece Extent=\"0 %d 0 %d 0 %d\">\n",
            NX-1,
            NY-1,
            NZ-1);

    fprintf(f,
"<PointData Scalars=\"density\">\n");

    fprintf(f,
"<DataArray "
"type=\"Float64\" "
"Name=\"density\" "
"format=\"ascii\">\n");

    for(int k=0;k<NZ;k++)
    for(int j=0;j<NY;j++)
    for(int i=0;i<NX;i++)
    {
        fprintf(f,"%e ",
                cell(i,j,k)->rho);
    }

    fprintf(f,"\n</DataArray>\n");

    fprintf(f,
"<DataArray "
"type=\"Float64\" "
"Name=\"temperature\" "
"format=\"ascii\">\n");

    for(int k=0;k<NZ;k++)
    for(int j=0;j<NY;j++)
    for(int i=0;i<NX;i++)
    {
        fprintf(f,"%e ",
                cell(i,j,k)->T);
    }

    fprintf(f,"\n</DataArray>\n");

    fprintf(f,
"<DataArray "
"type=\"Float64\" "
"Name=\"burnFraction\" "
"format=\"ascii\">\n");

    for(int k=0;k<NZ;k++)
    for(int j=0;j<NY;j++)
    for(int i=0;i<NX;i++)
    {
        fprintf(f,"%e ",
                cell(i,j,k)->burn);
    }

    fprintf(f,"\n</DataArray>\n");

    fprintf(f,
"<DataArray "
"type=\"Float64\" "
"Name=\"fusionPower\" "
"format=\"ascii\">\n");

    for(int k=0;k<NZ;k++)
    for(int j=0;j<NY;j++)
    for(int i=0;i<NX;i++)
    {
        fprintf(f,"%e ",
                cell(i,j,k)->fusionPower);
    }

    fprintf(f,"\n</DataArray>\n");

    fprintf(f,
"<DataArray "
"type=\"Float64\" "
"NumberOfComponents=\"3\" "
"Name=\"velocity\" "
"format=\"ascii\">\n");

    for(int k=0;k<NZ;k++)
    for(int j=0;j<NY;j++)
    for(int i=0;i<NX;i++)
    {
        Cell *c = cell(i,j,k);

        fprintf(f,
                "%e %e %e ",
                c->vx,
                c->vy,
                c->vz);
    }

    fprintf(f,
"\n</DataArray>\n");

    fprintf(f,
"</PointData>\n"
"<CellData/>\n"
"</Piece>\n"
"</ImageData>\n"
"</VTKFile>\n");

    fclose(f);

    printf("wrote %s\n",filename);
}

/*
    Placeholder update.

    Module 2 will replace this with:
        - implosion
        - Euler hydro
        - pressure gradients
        - compression
*/
void update_physics(double time)
{
    double cx = NX*0.5;
    double cy = NY*0.5;
    double cz = NZ*0.5;

    double driverStrength;

    if(time < 2.0e-6)
        driverStrength = 2.0e7;
    else
        driverStrength = 0.0;

    for(int k=1;k<NZ-1;k++)
    {
        for(int j=1;j<NY-1;j++)
        {
            for(int i=1;i<NX-1;i++)
            {
                Cell *c = cell(i,j,k);

                double dPdx =
                    (cell(i+1,j,k)->pressure -
                     cell(i-1,j,k)->pressure)
                    /(2.0*DX);

                double dPdy =
                    (cell(i,j+1,k)->pressure -
                     cell(i,j-1,k)->pressure)
                    /(2.0*DX);

                double dPdz =
                    (cell(i,j,k+1)->pressure -
                     cell(i,j,k-1)->pressure)
                    /(2.0*DX);

                double x=(i-cx)*DX;
                double y=(j-cy)*DX;
                double z=(k-cz)*DX;

                double r =
                    sqrt(x*x+y*y+z*z)
                    + 1e-20;

                double ex = x/r;
                double ey = y/r;
                double ez = z/r;

                double ax =
                    -(dPdx/c->rho)
                    - driverStrength*ex;

                double ay =
                    -(dPdy/c->rho)
                    - driverStrength*ey;

                double az =
                    -(dPdz/c->rho)
                    - driverStrength*ez;

                double vx =
                    c->vx + ax*DT;

                double vy =
                    c->vy + ay*DT;

                double vz =
                    c->vz + az*DT;

                double divU =
                    divergence_velocity(
                        i,j,k);
                double rho =
                    c->rho *
                    (1.0 - divU*DT);
		
                if(rho < 1e-3)
                    rho = 1e-3;

                double gamma = 5.0/3.0;

                double T =
                    c->T *
                    pow(rho/c->rho,
                        gamma-1.0);
		if(divU > 0.0)
                {
                double cooling = exp(-0.20 * divU * DT);
                T *= cooling;
                }


                if(time > 2.5e-6 && divU > 0.0)
                {
                T *= exp(-2.0 * divU * DT);
                }  
                double qvisc = 0.0;

                if(divU < 0.0)
                {
                    qvisc =
                        0.5 *
                        rho *
                        divU *
                        divU *
                        DX*DX;

                    T += qvisc*1e-4;
                }

                int idx =
                    IDX(i,j,k);

                newRho[idx] = rho;
                newT[idx]   = T;

                newVx[idx]  = vx;
                newVy[idx]  = vy;
                newVz[idx]  = vz;
            }
        }
    }

    for(int k=1;k<NZ-1;k++)
    {
        for(int j=1;j<NY-1;j++)
        {
            for(int i=1;i<NX-1;i++)
            {
                int idx =
                    IDX(i,j,k);

                Cell *c =
                    cell(i,j,k);

                c->rho =
                    newRho[idx];

                c->T =
                    newT[idx];

                c->vx =
                    newVx[idx];

                c->vy =
                    newVy[idx];

                c->vz =
                    newVz[idx];

                c->pressure =
                    ideal_pressure(
                        c->rho,
                        c->T);
            }
        }
    }
}
int main(void)
{
    allocate_grid();

    initialize_capsule();

    double time = 0.0;
    double nextSnapshot = 0.0;

    int frame = 0;

    while(time <= TOTAL_TIME)
    {
        update_physics(time);

        if(time >= nextSnapshot)
        {
            write_vti_snapshot(frame);

            frame++;

            nextSnapshot +=
                SNAPSHOT_INTERVAL;
        }

        time += DT;
    }

    free(grid);

    return 0;
}
