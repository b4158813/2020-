#include<bits/stdc++.h>
using namespace std;

const int Q = 9;
const int NX = 500;
const int NY = 500;
int e[Q][2] = {{0,0},{1,0},{0,1},{-1,0},{0,-1},{1,1},{-1,1},{-1,-1},{1,-1}};
double w[Q] = {4.0/9,1.0/9,1.0/9,1.0/9,1.0/9,1.0/36,1.0/36,1.0/36,1.0/36};
double gamma[Q] = {0,1.0/3,1.0/3,1.0/3,1.0/3,1.0/12,1.0/12,1.0/12,1.0/12};
double rho[NX+1][NY+1], u[NX+1][NY+1][2], u0[NX+1][NY+1][2], f[NX+1][NY+1][Q], F[NX+1][NY+1][Q];
double U,Re,dx,dy,dt,Lx,Ly,rho0,tau,niu,c,cs,f_x,f_y;
const int rec_num = 53;
double zoom = 5;
struct rectangle
{
    int x1,y1,x2,y2;
}rec[rec_num];


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

inline bool in_any_rec(int x,int y){
    for(int i=0;i<rec_num;i++){
        if(inrec(x,y,i)) return true;
    }
    return false;
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

inline void set_single_rec_geo(int id,int x1,int y1,int x2,int y2){
    rec[id] = rectangle{x1,y1,x2,y2};
}

inline void set_rec_geo(){
    for(int i=0;i<4;i++) set_single_rec_geo(i,10+20*i,5,20+20*i,9);
    for(int i=4;i<6;i++) set_single_rec_geo(i,45,21+(i-4)*20,55,31+(i-4)*20);
    for(int i=6;i<11;i++) set_single_rec_geo(i,23+(i-6)*10,74,30+(i-6)*10,90);
    set_single_rec_geo(11,10,80,20,90);
    set_single_rec_geo(12,75,80,85,90);
    set_single_rec_geo(13,10,69,18,74);
    set_single_rec_geo(14,47,59,52,69);
    for(int i=15;i<19;i++) set_single_rec_geo(i,5,15+(i-15)*13,10,25+(i-15)*13);
    for(int i=19;i<23;i++) set_single_rec_geo(i,13,15+(i-19)*13,18,25+(i-19)*13);
    for(int i=23;i<27;i++) set_single_rec_geo(i,21,20+(i-23)*13,26,30+(i-23)*13);
    for(int i=27;i<31;i++) set_single_rec_geo(i,29,20+(i-27)*13,34,30+(i-27)*13);
    for(int i=31;i<35;i++) set_single_rec_geo(i,37,20+(i-31)*13,42,30+(i-31)*13);
    for(int i=35;i<39;i++) set_single_rec_geo(i,58,20+(i-35)*13,63,30+(i-35)*13);
    for(int i=39;i<43;i++) set_single_rec_geo(i,66,20+(i-39)*13,71,30+(i-39)*13);
    for(int i=43;i<48;i++) set_single_rec_geo(i,74,15+(i-43)*13,79,25+(i-43)*13);
    for(int i=48;i<53;i++) set_single_rec_geo(i,82,15+(i-48)*13,87,25+(i-48)*13);

    for(int i=0;i<rec_num;i++){
        rec[i].x1 *= zoom;
        rec[i].y1 *= zoom;
        rec[i].x2 *= zoom;
        rec[i].y2 *= zoom;
    }
}

inline void init(){
    dx = 1.0;
    dy = 1.0;
    Lx = 1.0*dx*NX;
    Ly = 1.0*dy*NY;
    dt = 1.0;
    // U = 0.1;
    rho0 = 1.0;
    f_x = 1e-5;
    f_y = 0;
    // Re = 1000;
    niu = 0.01;
    c = dx/dt;
    cs = c/sqrt(3.0);
    tau = niu/(cs*cs) + 0.5*dt;

    set_rec_geo();

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
            if(in_any_rec(i,j)) continue;
            for(int k=0;k<Q;k++){
                int ip = i - e[k][0];
                int jp = j - e[k][1];
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
            if(in_any_rec(i,j)){
                rho[i][j] = 1.0;
            }else{
                for(int k=0;k<Q;k++){
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

    // // left boundary
    // for(int j=1;j<NY;j++){
    //     rho[0][j] = rho[1][j];
    //     for(int k=0;k<Q;k++){
    //         f[0][j][k] = feq(k,rho[0][j],u[0][j]) + (1-1/tau) * (f[1][j][k] - feq(k,rho[1][j],u[1][j]));
    //     }
    // }

    // // right boundary
    // for(int j=1;j<NY;j++){
    //     rho[NX][j] = rho[NX-1][j];
    //     for(int k=0;k<Q;k++){
    //         f[NX][j][k] = feq(k,rho[NX][j],u[NX][j]) + (1-1/tau) * (f[NX-1][j][k] - feq(k,rho[NX-1][j],u[NX-1][j]));
    //     }
    // }

    // top boundary
    for(int i=0;i<=NX;i++){
        rho[i][NY] = rho[i][NY-1];
        for(int k=0;k<Q;k++){
            f[i][NY][k] = feq(k,rho[i][NY],u[i][NY]) + (1-1/tau) * (f[i][NY-1][k] - feq(k,rho[i][NY-1],u[i][NY-1]));
        }
    }

    // bottom boundary
    for(int i=0;i<=NX;i++){
        rho[i][0] = rho[i][1];
        for(int k=0;k<Q;k++){
            f[i][0][k] = feq(k,rho[i][0],u[i][0]) + (1-1/tau) * (f[i][1][k] - feq(k,rho[i][1],u[i][1]));
        }
    }

    for(int i=0;i<rec_num;i++)
        set_rec_boundary(i);

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
    name<<"cpu_simulation"<<num<<".dat";
    ofstream out(name.str().c_str());
    out<<"Title= \"LBM CPU SIMULATION\"\n"<<"VARIABLES=\"X\",\"Y\",\"U\",\"V\"\n"<<"ZONE T=\"BOX\",I="<<NX+1<<",J="<<NY+1<<",F=POINT"<<endl;
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