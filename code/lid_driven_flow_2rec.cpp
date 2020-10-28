#include<bits/stdc++.h>
using namespace std;

const int Q = 9;
const int NX = 100;
const int NY = 100;
int e[Q][2] = {{0,0},{1,0},{0,1},{-1,0},{0,-1},{1,1},{-1,1},{-1,-1},{1,-1}};
double w[Q] = {4.0/9,1.0/9,1.0/9,1.0/9,1.0/9,1.0/36,1.0/36,1.0/36,1.0/36};
double rho[NX+1][NY+1], u[NX+1][NY+1][2], u0[NX+1][NY+1][2], f[NX+1][NY+1][Q], F[NX+1][NY+1][Q];
double U,Re,dx,dy,dt,Lx,Ly,rho0,tau,niu,c,cs;
struct rectangle
{
    int x1,y1,x2,y2;
}rec[2];


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

inline void set_rec_boundary(int id){
    // rectangle boundary: left & right side
    for(int j=rec[id].y1;j<=rec[id].y2;j++){
        rho[rec[id].x1][j] = rho[rec[id].x1-1][j];
        rho[rec[id].x2][j] = rho[rec[id].x2+1][j];
        for(int k=0;k<Q;k++){
            f[rec[id].x1][j][k] = feq(k,rho[rec[id].x1][j],u[rec[id].x1][j]) + (1-1/tau) * (f[rec[id].x1-1][j][k] - feq(k,rho[rec[id].x1-1][j],u[rec[id].x1-1][j]));
            f[rec[id].x2][j][k] = feq(k,rho[rec[id].x2][j],u[rec[id].x2][j]) + (1-1/tau) * (f[rec[id].x2+1][j][k] - feq(k,rho[rec[id].x2+1][j],u[rec[id].x2+1][j]));
        }
    }
    // rectangle boundary: bottom & top side
    for(int i=rec[id].x1;i<=rec[id].x2;i++){
        rho[i][rec[id].y1] = rho[i][rec[id].y1-1];
        rho[i][rec[id].y2] = rho[i][rec[id].y2+1];
        for(int k=0;k<Q;k++){
            f[i][rec[id].y1][k] = feq(k,rho[i][rec[id].y1],u[i][rec[id].y1]) + (1-1/tau) * (f[i][rec[id].y1-1][k] - feq(k,rho[i][rec[id].y1-1],u[i][rec[id].y1-1]));
            f[i][rec[id].y2][k] = feq(k,rho[i][rec[id].y2],u[i][rec[id].y2]) + (1-1/tau) * (f[i][rec[id].y2+1][k] - feq(k,rho[i][rec[id].y2+1],u[i][rec[id].y2+1]));
        }
    }
}

inline void init(){
    dx = 1.0;
    dy = 1.0;
    Lx = 1.0*dx*NX;
    Ly = 1.0*dy*NY;
    dt = 1.0;
    U = 0.1;
    rho0 = 1.0;
    Re = 3000;
    niu = U*Lx/Re;
    c = dx/dt;
    cs = c/sqrt(3.0);
    tau = niu/(cs*cs) + 0.5*dt;
    // 1st rectangle
    rec[0].x1 = 50;
    rec[0].y1 = 60;
    rec[0].x2 = 60;
    rec[0].y2 = 80;

    // 2nd rectangle
    rec[1].x1 = 60;
    rec[1].y1 = 20;
    rec[1].x2 = 80;
    rec[1].y2 = 40;

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
            if(inrec(i,j,0)||inrec(i,j,1)){
                continue;
            }else{
                for(int k=0;k<Q;k++){
                    int ip = i - e[k][0];
                    int jp = j - e[k][1];
                    F[i][j][k] = f[ip][jp][k] + (feq(k,rho[ip][jp],u[ip][jp]) - f[ip][jp][k])/tau;
                }
            }
        }
    }

    memcpy(u0,u,sizeof(u));
    memset(rho,0,sizeof(rho));
    memset(u,0,sizeof(u));

    // streaming
    for(int i=1;i<NX;i++){
        for(int j=1;j<NY;j++){
            if(inrec(i,j,0)||inrec(i,j,1)){
                rho[i][j] = 1.0;
            }else{
                for(int k=0;k<Q;k++){
                    f[i][j][k] = F[i][j][k];
                    rho[i][j] += f[i][j][k];
                    u[i][j][0] += e[k][0] * f[i][j][k];
                    u[i][j][1] += e[k][1] * f[i][j][k];
                }
            }
            u[i][j][0] /= rho[i][j];
            u[i][j][1] /= rho[i][j];
        }
    }

    // left boundary
    for(int j=1;j<NY;j++){
        rho[0][j] = rho[1][j];
        for(int k=0;k<Q;k++){
            f[0][j][k] = feq(k,rho[0][j],u[0][j]) + (1-1/tau) * (f[1][j][k] - feq(k,rho[1][j],u[1][j]));
        }
    }

    // right boundary
    for(int j=1;j<NY;j++){
        rho[NX][j] = rho[NX-1][j];
        // if(j>c1y){
            for(int k=0;k<Q;k++){
                f[NX][j][k] = feq(k,rho[NX][j],u[NX][j]) + (1-1/tau) * (f[NX-1][j][k] - feq(k,rho[NX-1][j],u[NX-1][j]));
            } 
        // }
        
    }

    // top boundary
    for(int i=0;i<=NX;i++){
        rho[i][NY] = rho[i][NY-1];
        u[i][NY][0] = U;
        for(int k=0;k<Q;k++){
            f[i][NY][k] = feq(k,rho[i][NY],u[i][NY]) + (1-1/tau) * (f[i][NY-1][k] - feq(k,rho[i][NY-1],u[i][NY-1]));
        }
    }

    // bottom boundary
    for(int i=0;i<=NX;i++){
        rho[i][0] = rho[i][1];
        // if(i<c1x){
            for(int k=0;k<Q;k++){
                f[i][0][k] = feq(k,rho[i][0],u[i][0]) + (1-1/tau) * (f[i][1][k] - feq(k,rho[i][1],u[i][1]));
            }
        // }
    }

    set_rec_boundary(0);
    set_rec_boundary(1);

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
    name<<"lid_driven_flow"<<num<<".dat";
    ofstream out(name.str().c_str());
    out<<"Title= \"LBM Lid Driven Flow\"\n"<<"VARIABLES=\"X\",\"Y\",\"U\",\"V\"\n"<<"ZONE T=\"BOX\",I="<<NX+1<<",J="<<NY+1<<",F=POINT"<<endl;
    for(int i=0;i<=NX;i++){
        for(int j=0;j<=NY;j++){
            out<<double(i)<<" "<<double(j)<<" "<<u[i][j][0]<<" "<<u[i][j][1]<<endl;
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