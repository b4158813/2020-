#include<bits/stdc++.h>
using namespace std;

const int Q = 9;
const int NX = 100;
const int NY = 100;
int e[Q][2] = {{0,0},{1,0},{0,1},{-1,0},{0,-1},{1,1},{-1,1},{-1,-1},{1,-1}};
double w[Q] = {4.0/9,1.0/9,1.0/9,1.0/9,1.0/9,1.0/36,1.0/36,1.0/36,1.0/36};
double gamma[Q] = {0,1.0/3,1.0/3,1.0/3,1.0/3,1.0/12,1.0/12,1.0/12,1.0/12};
double rho[NX+1][NY+1], u[NX+1][NY+1][2], u0[NX+1][NY+1][2], f[NX+1][NY+1][Q], F[NX+1][NY+1][Q];
double U,Re,dx,dy,dt,Lx,Ly,rho0,tau,niu,c,cs,f_x,f_y;

struct rectangle
{
    int x1,y1,x2,y2;
}rec[1];

inline double feq(int k,double rho,double u[2]){
    double res;
    double eu = e[k][0]*u[0]+e[k][1]*u[1];
    double uu = u[0]*u[0]+u[1]*u[1];
    res = w[k] * rho * (1 + eu/(cs*cs) + (eu*eu)/(2*cs*cs*cs*cs) - uu/(2*cs*cs));
    return res;
}

inline bool inrec(int x,int y,int id){
    return x>=rec[id].x1&&x<=rec[id].x2&&y>=rec[id].y1&&y<=rec[id].y2;
}

inline void set_single_rec_geo(int id,int x1,int y1,int x2,int y2){
    rec[id] = rectangle{x1,y1,x2,y2};
}

inline void init(){
    dx = 1.0;
    dy = 1.0;
    Lx = 1.0*dx*NX;
    Ly = 1.0*dy*NY;
    dt = 1.0;
    rho0 = 1.0;
    f_x = 1e-5;
    f_y = 0;
    Re = 3000;
    niu = 0.01;
    c = dx/dt;
    cs = c/sqrt(3.0);
    tau = niu/(cs*cs) + 0.5*dt;

    set_single_rec_geo(0,30,30,60,60);

    memset(u,0,sizeof(u));
    for(int i=0;i<=NX;i++){
        for(int j=0;j<=NY;j++){
            rho[i][j] = rho0;
            for(int k=0;k<Q;k++){
                f[i][j][k] = feq(k,rho[i][j],u[i][j]);
            }
        }
    }
}

inline void evolution(){
    // collision & streaming
    for(int i=0;i<=NX;i++){
        for(int j=1;j<NY;j++){
            if(inrec(i,j,0)) continue;
            for(int k=0;k<Q;k++){
                int ip = i - e[k][0];
                int jp = j - e[k][1];
                // periodic boundary
                if(ip<0) ip=NX;
                if(ip>NX) ip=0;
                F[i][j][k] = f[ip][jp][k] + (feq(k,rho[ip][jp],u[ip][jp]) - f[ip][jp][k])/tau;
            }
        }
    }

    memcpy(u0,u,sizeof(u));
    memset(rho,0,sizeof(rho));
    memset(u,0,sizeof(u));

    // get_macro_value
    for(int i=0;i<=NX;i++){
        for(int j=1;j<NY;j++){
            if(inrec(i,j,0)) rho[i][j] = 1.0;
            else{
                for(int k=0;k<Q;k++){
                	// pressure boundary
                    f[i][j][k] = F[i][j][k] + gamma[k] * (f_x * e[k][0] + f_y * e[k][1]);
                    rho[i][j] += f[i][j][k];
                    u[i][j][0] += e[k][0] * f[i][j][k];
                    u[i][j][1] += e[k][1] * f[i][j][k];
                }
            }
            u[i][j][0] /= rho[i][j];
            u[i][j][1] /= rho[i][j];
        }
    }

    // top & bottom boundary
    for(int i=0;i<=NX;i++){
        rho[i][0] = rho[i][1];
        rho[i][NY] = rho[i][NY-1];
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
    name<<"poiseuille_flow_with_rec"<<num<<".dat";
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