#include "/home/xdengae/xdrfile-1.1.4/include/xdrfile_xtc.h"
#include "/home/xdengae/xdrfile-1.1.4/include/xdrfile_trr.h"
#include "/home/xdengae/xdrfile-1.1.4/include/xdrfile.h"
#include <stdlib.h>
#include <stdio.h>                      
#include <iostream>
#include "device_launch_parameters.h"
#include "cuda_runtime.h"
#include <math.h>

using namespace std;

class Complex{
    float a,b;
    // (a+bi)*(c+di)=(ac-bd)+(ad+bc)i       (a+bi)/(c+di)={(ac+bd)+(bc-ad)i}/(c^2+d^2)
public:
    __device__ __host__ Complex();
    __device__ __host__ Complex(float real,float imag){
        a=real;
        b=imag;
    }
    __device__ __host__ Complex operator * (Complex c){
        Complex ctmp;
        ctmp.a=(a*c.a-b*c.b);
        ctmp.b=(a*c.b+b*c.a);
        return ctmp;
    }
    __device__ __host__ Complex operator / (Complex c){
        Complex ctmp;
        ctmp.a=(a*c.a+b*c.b)/(c.a*c.a+c.b*c.b);
        ctmp.b=(b*c.a-a*c.b)/(c.a*c.a+c.b*c.b);
        return ctmp;
    }
    __device__ __host__ Complex operator + (Complex c){
        Complex ctmp;
        ctmp.a=a+c.a;
        ctmp.b=b+c.b;
        return ctmp;
    }
    __device__ __host__ Complex operator - (Complex c){
        Complex ctmp;
        ctmp.a=a-c.a;
        ctmp.b=b-c.b;
        return ctmp;
    }
    __device__ __host__ float real(){
        return a;
    }
    __device__ __host__ float imag(){
        return b;
    }
};

void getNumOfFrameAtom(char* path,int* numOfAtom,int* numOfFrame){
    int numOfAtoms=0;
    int step=0;
    float mdtime,lambda;
    matrix box;
    XDRFILE *trr=xdrfile_open(path,"r");
    int read_return=read_trr_natoms(path,&numOfAtoms);
    if(read_return!=exdrOK){
        cout<<"fail to read trr file"<<endl;
        exit(0);
    }
    rvec* x=(rvec * )calloc(numOfAtoms,sizeof(x[0]));
    rvec* v=(rvec * )calloc(numOfAtoms,sizeof(x[0]));
    rvec* f=(rvec * )calloc(numOfAtoms,sizeof(x[0]));
    int ac=0;
    while((read_return=read_trr(trr,numOfAtoms,&step,&mdtime,&lambda,box,x,v,f))==exdrOK){
        ac++;
    }
    *numOfFrame=ac;
    *numOfAtom=numOfAtoms;
    free(x);
    free(v);
    free(f);
    xdrfile_close(trr);
}
// counting base start from 0
void extract_trr_file(char* path,float* vx,float* vy,float* vz,float* coor_x,float* coor_z,int startFrame,int endFrame,float* H,float* L){
    int numOfAtoms=0;
    int step=0;
    float mdtime,lambda;
    matrix box;
    XDRFILE *trr=xdrfile_open(path,"r");
    int read_return=read_trr_natoms(path,&numOfAtoms);
    if(read_return!=exdrOK){
        cout<<"fail to read trr file"<<endl;
        exit(0);
    }
    rvec* x=(rvec * )calloc(numOfAtoms,sizeof(x[0]));
    rvec* v=(rvec * )calloc(numOfAtoms,sizeof(x[0]));
    rvec* f=(rvec * )calloc(numOfAtoms,sizeof(x[0]));
    int frameIdx=0;
    int frame_tag=0; 
    while((read_return=read_trr(trr,numOfAtoms,&step,&mdtime,&lambda,box,x,v,f))==exdrOK){
        if(step==0){
            cout<<"Basic properties of current MD system: "<<endl;
            printf("\tatom num: %d\t",numOfAtoms);
            printf("\tbox size is : %f * %f * %f\n",box[0][0],box[1][1],box[2][2]);
            *H=box[2][2]/2.;
            *L=box[0][0]/2.;
            printf("\tThe fluid atom's number density is %f\n",numOfAtoms/(box[0][0]*box[1][1]*box[2][2]));
        }
        if(frame_tag>=startFrame&&frame_tag<=endFrame){
            for(int i=0;i<numOfAtoms;i++){
                vx[frameIdx*numOfAtoms+i]=v[i][0];
                vy[frameIdx*numOfAtoms+i]=v[i][1];
                vz[frameIdx*numOfAtoms+i]=v[i][2];
                coor_x[frameIdx*numOfAtoms+i]=x[i][0];
                coor_z[frameIdx*numOfAtoms+i]=x[i][2];
            }
            frameIdx++;
        }
        if(frame_tag>endFrame){
            break;
        }
        frame_tag++;
    }
    printf("%d frames have been extracted from trr file to calculate the correlation function\n",frameIdx);
    xdrfile_close(trr);
}
/*
Complex omega(float kx,float kz,float R){
    Complex kx_(kx,0);
    Complex kz_(kz,0);
    Complex i_(0,1);
    Complex R_(R,0);
    Complex ctmp=(kx_*kx_+kz_*kz_)/i_/R_;
    return ctmp;
}
*/
__device__ __host__ float velocityX(float kx,float kz,float coor_x,float coor_z,float H,float alpha){
    // calculate lambda
    float lambda=-kx/(kx*kx+kz*kz)*cos(kz*H)*sin(alpha)/cosh(kx*H);
    float velocity=(lambda*sinh(kx*coor_z)+kz/(kx*kx+kz*kz)*cos(kz*coor_z+alpha))*cos(kx*coor_x);
    return velocity;
}
__device__ __host__ float velocityZ(float kx,float kz,float coor_x,float coor_z,float H,float alpha){
    // calculate lambda
    
    return 1.;
}
void generate_kz(float* kz,int N,float dk){
    for(int i=0;i<N;i++){
        kz[i]=0.001+i*dk;
    }
}
bool calCorrelationXKernel_TEST(float kx,float kz,int numOfFrame,int numOfAtom, float* vx,float* vz,float* coor_x,float* coor_z,int deltatIdx,float result_gpu_x,float result_gpu_z,float H,float alpha){
    float corX_up=0.;
    float corX_down=0.;
    float corZ_up=0.;
    float corZ_down=0.;
    for(int frameIdx=0;frameIdx<numOfFrame-deltatIdx;frameIdx++){
        float vx_sum_t0=0.;
        float vx_sum_deltat=0.;
        float vz_sum_t0=0.;
        float vz_sum_deltat=0.;
        for(int atomIdx=0;atomIdx<numOfAtom;atomIdx++){
            int idx_t0=frameIdx*numOfAtom+atomIdx;
            int idx_deltat=(frameIdx+deltatIdx)*numOfAtom+atomIdx;
            vx_sum_t0+=vx[idx_t0]*velocityX(kx,kz,coor_x[idx_t0],coor_z[idx_t0],H,alpha);
            vx_sum_deltat+=vx[idx_deltat]*velocityX(kx,kz,coor_x[idx_deltat],coor_z[idx_deltat],H,alpha);
            vz_sum_t0+=vz[idx_t0]*velocityZ(kx,kz,coor_x[idx_t0],coor_z[idx_t0],H,alpha);
            vz_sum_deltat+=vz[idx_deltat]*velocityZ(kx,kz,coor_x[idx_deltat],coor_z[idx_deltat],H,alpha);
        }
        corX_up+=vx_sum_t0*vx_sum_deltat;
        corX_down+=vx_sum_t0*vx_sum_t0;
        corZ_up+=vz_sum_t0*vz_sum_deltat;
        corZ_down+=vz_sum_t0*vz_sum_t0;
    }
    float corX=corX_up/corX_down;
    float corZ=corZ_up/corZ_down;
    float relative_error_x=fabs(corX-result_gpu_x)/corX;
    float relative_error_z=fabs(corZ-result_gpu_z)/corZ;
    if(relative_error_x<0.0001&&relative_error_z<0.0001){
        return true;
    }
    return false;
}
float calCorrelationXKernel_CPU(float kx,float kz,int numOfFrame,int numOfAtom, float* vx,float* vz,float* coor_x,float* coor_z,int deltatIdx,float H,float alpha){
    float corX_up=0.;
    float corX_down=0.;
    for(int frameIdx=0;frameIdx<numOfFrame-deltatIdx;frameIdx++){
        float vx_sum_t0=0.;
        float vx_sum_deltat=0.;
        for(int atomIdx=0;atomIdx<numOfAtom;atomIdx++){
            int idx_t0=frameIdx*numOfAtom+atomIdx;
            int idx_deltat=(frameIdx+deltatIdx)*numOfAtom+atomIdx;
            vx_sum_t0+=vx[idx_t0]*velocityX(kx,kz,coor_x[idx_t0],coor_z[idx_t0],H,alpha);
            vx_sum_deltat+=vx[idx_deltat]*velocityX(kx,kz,coor_x[idx_deltat],coor_z[idx_deltat],H,alpha);
        }
        corX_up+=vx_sum_t0*vx_sum_deltat;
        corX_down+=vx_sum_t0*vx_sum_t0;
    }
    float corX=corX_up/corX_down;
    return corX;
}
void write_temperature_per_frame(float* table,int num){ // frame number
    FILE* fp=fopen("/home/xdengae/LJ_Fluid/DATA/temperature_table","w+");
    if(fp==NULL){
        cout<<"Fail to write temperature! file io error "<<endl;
        exit(0);
    }
    for(int i=0;i<num;i++){
        fprintf(fp,"%d\t%f\n",i,table[i]);
    }
    fclose(fp);
    cout<<"Temperature per frame has been written to file [temperature_table]"<<endl;
}