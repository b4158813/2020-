#include<bits/stdc++.h>
using namespace std;

const int Q = 9;
const int NX = 200;
const int NY = 100;
int e[Q][2] = {{0,0},{1,0},{0,1},{-1,0},{0,-1},{1,1},{-1,1},{-1,-1},{1,-1}};
double w[Q] = {4.0/9,1.0/9,1.0/9,1.0/9,1.0/9,1.0/36,1.0/36,1.0/36,1.0/36};
double rho[NX+1][NY+1], u[NX+1][NY+1][2], u0[NX+1][NY+1][2], f[NX+1][NY+1][Q], F[NX+1][NY+1][Q];
double U,Re,dx,dy,dt,Lx,Ly,rho0,tau,niu,c,cs;

inline double feq(int k,double rho,double u[2]){
    double res;
    double eu = e[k][0]*u[0]+e[k][1]*u[1];
    double uu = u[0]*u[0]+u[1]*u[1];
    res = w[k] * rho * (1 + eu/(cs*cs) + (eu*eu)/(2*cs*cs*cs*cs) - uu/(2*cs*cs));
    return res;
}

inline void init(){
    dx = 1.0;
    dy = 1.0;
    Lx = 1.0*dx*NX;
    Ly = 1.0*dy*NY;
    dt = 1.0;
    U = 0.1;
    rho0 = 1.0;
    Re = 100;
    niu = U*Lx/Re;
    c = dx/dt;
    cs = c/sqrt(3.0);
    tau = niu/(cs*cs) + 0.5*dt;
    memset(u,0,sizeof(u));
    for(int i=0;i<=NX;i++){
        u[i][NY][0] = U; // lid speed = U
        for(int j=0;j<=NY;j++){
            rho[i][j] = rho0;
            for(int k=0;k<Q;k++){
                f[i][j][k] = feq(k,rho[i][j],u[i][j]);
            }
        }
    }
}

inline void evolution(){
    // collision
    for(int i=1;i<NX;i++){
        for(int j=1;j<NY;j++){
            for(int k=0;k<Q;k++){
                int ip = i - e[k][0];
                int jp = j - e[k][1];
                F[i][j][k] = f[ip][jp][k] + (feq(k,rho[ip][jp],u[ip][jp]) - f[ip][jp][k])/tau;
            }
        }
    }

    memcpy(u0,u,sizeof(u));
    memset(rho,0,sizeof(rho));
    memset(u,0,sizeof(u));

    // streaming
    for(int i=1;i<NX;i++){
        for(int j=1;j<NY;j++){
            for(int k=0;k<Q;k++){
                f[i][j][k] = F[i][j][k];
                rho[i][j] += f[i][j][k];
                u[i][j][0] += e[k][0] * f[i][j][k];
                u[i][j][1] += e[k][1] * f[i][j][k];
            }
            u[i][j][0] /= rho[i][j];
            u[i][j][1] /= rho[i][j];
        }
    }

    // left & right boundary
    for(int j=1;j<NY;j++){
        rho[0][j] = rho[1][j];
        rho[NX][j] = rho[NX-1][j];
        for(int k=0;k<Q;k++){
            f[0][j][k] = feq(k,rho[0][j],u[0][j]) + (1-1/tau) * (f[1][j][k] - feq(k,rho[1][j],u[1][j]));
            f[NX][j][k] = feq(k,rho[NX][j],u[NX][j]) + (1-1/tau) * (f[NX-1][j][k] - feq(k,rho[NX-1][j],u[NX-1][j]));
        }
    }

    // top & bottom boundary
    for(int i=0;i<=NX;i++){
        rho[i][0] = rho[i][1];
        rho[i][NY] = rho[i][NY-1];
        u[i][NY][0] = U;
        for(int k=0;k<Q;k++){
            f[i][0][k] = feq(k,rho[i][0],u[i][0]) + (1-1/tau) * (f[i][1][k] - feq(k,rho[i][1],u[i][1]));
            f[i][NY][k] = feq(k,rho[i][NY],u[i][NY]) + (1-1/tau) * (f[i][NY-1][k] - feq(k,rho[i][NY-1],u[i][NY-1]));
        }
    }
}

inline double ERR(){
    double tp1=0,tp2=0;
    for(int i=1;i<NX;i++){
        for(int j=1;j<NY;j++){
            tp1 += ((u[i][j][0]-u0[i][j][0])*(u[i][j][0]-u0[i][j][0])+(u[i][j][1]-u0[i][j][1])*(u[i][j][1]-u0[i][j][1]));
            tp2 += (u[i][j][0]*u[i][j][0]+u[i][j][1]*u[i][j][1]);
        }
    }
    tp1 = sqrt(tp1);
    tp2 = sqrt(tp2);
    return tp1/(tp2+1e-30);
}

inline void output(int num){
    ostringstream name;
    name<<"couette_flow"<<num<<".dat";
    ofstream out(name.str().c_str());
    out<<"Title= \"LBM Couette Flow\"\n"<<"VARIABLES=\"X\",\"Y\",\"U\",\"V\"\n"<<"ZONE T=\"BOX\",I="<<NX+1<<",J="<<NY+1<<",F=POINT"<<endl;
    for(int i=0;i<=NX;i++){
        for(int j=0;j<=NY;j++){
            out<<double(i)/Lx<<" "<<double(j)/Ly<<" "<<u[i][j][0]<<" "<<u[i][j][1]<<endl;
        }
    }
}

int main(){
    init();
    for(int i=1;;i++){
        evolution();
        if(i%100==0){
            double error = ERR();
            printf("%dth error = %e\n",i,error);
            if(i>=1000){
                if(i%1000==0) output(i);
                if(error<1.0e-5) break;
            }
        }
    }
}