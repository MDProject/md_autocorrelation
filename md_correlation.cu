#include "/home/xdengae/xdrfile-1.1.4/include/xdrfile_xtc.h"
#include "/home/xdengae/xdrfile-1.1.4/include/xdrfile_trr.h"
#include "/home/xdengae/xdrfile-1.1.4/include/xdrfile.h"
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "correlation_lib.h"

#define threadsPerBlock 8
#define blocksPerGrid 8

using namespace std;
#define PI 3.14159265359

void handleError(cudaError_t cu) {
	if (cu != cudaSuccess) {
		printf("%s\n", cudaGetErrorString(cu));
		system("pause");
		exit(0);
	}
}

__global__ void calCorrelationXKernel(float kx,float kz,int numOfFrame,int numOfAtom, float* vx,float* vz,float* coor_x,float* coor_z,int deltatIdx,float* corXup,float* corXdown,float* corZup,float* corZdown,float H,float alpha){
    // corX : store the partial summation of correlation function for kz at different frame      deltatIdx:    time interval's index dt=deltatIdx*delta_t
    // using shared memory to calculate the sum(corXFrame)
    // blockIdx --> frameIdx    threadIdx --> particle idx
    // corXPerFrame1 --> accumulate the summation of vx(z)[t0] * vec_x(z)[t0] of all atoms; corXPerFrame2 --> accumulate the summation of vx(z)[t0+dt] * vec_x(z)[t0+dt]
	__shared__ float corXPerFrame1[threadsPerBlock]; // cacheIdx12
	__shared__ float corXPerFrame2[threadsPerBlock]; // cacheIdx12
    __shared__ float corZPerFrame1[threadsPerBlock]; // cacheIdx12
	__shared__ float corZPerFrame2[threadsPerBlock];
	float corXFrameTmp_up = 0.;
	float corXFrameTmp_down = 0.;
    float corZFrameTmp_up = 0.;
	float corZFrameTmp_down = 0.;
	int frameIdx = blockIdx.x;
	while (frameIdx<numOfFrame - deltatIdx) {
        int cacheIdx12 = threadIdx.x;
        float corXPerFrame1Tmp = 0.;
	    float corXPerFrame2Tmp = 0.;
        float corZPerFrame1Tmp = 0.;
        float corZPerFrame2Tmp = 0.;
		while (cacheIdx12<numOfAtom) {
			int atomIdx_t0 = frameIdx*numOfAtom + cacheIdx12;
			int atomIdx_deltat = (frameIdx + deltatIdx)*numOfAtom + cacheIdx12;
			corXPerFrame1Tmp += vx[atomIdx_t0] * velocityX(kx, kz, coor_x[atomIdx_t0], coor_z[atomIdx_t0],H,alpha);
			corXPerFrame2Tmp += vx[atomIdx_deltat] * velocityX(kx, kz, coor_x[atomIdx_deltat], coor_z[atomIdx_deltat],H,alpha);
            corZPerFrame1Tmp += vz[atomIdx_t0] * velocityZ(kx, kz, coor_x[atomIdx_t0], coor_z[atomIdx_t0],H,alpha);
			corZPerFrame2Tmp += vz[atomIdx_deltat] * velocityZ(kx, kz, coor_x[atomIdx_deltat], coor_z[atomIdx_deltat],H,alpha);
			cacheIdx12 += blockDim.x;
		}
		corXPerFrame1[threadIdx.x] = corXPerFrame1Tmp;
		corXPerFrame2[threadIdx.x] = corXPerFrame2Tmp;
        corZPerFrame1[threadIdx.x] = corZPerFrame1Tmp;
		corZPerFrame2[threadIdx.x] = corZPerFrame2Tmp;
		__syncthreads();
		if (threadIdx.x == 0) {
			float sum_vx_t0 = 0.;
			float sum_vx_delta = 0.;
            float sum_vz_t0 = 0.;
            float sum_vz_delta = 0.;
			for (int i = 0; i < threadsPerBlock; i++) {
				sum_vx_t0 += corXPerFrame1[i];
				sum_vx_delta += corXPerFrame2[i];
                sum_vz_t0 += corZPerFrame1[i];
                sum_vz_delta += corZPerFrame2[i];
			}
			corXFrameTmp_up += sum_vx_t0*sum_vx_delta;
            corXFrameTmp_down += sum_vx_t0*sum_vx_t0;
            corZFrameTmp_up += sum_vz_t0*sum_vz_delta;
            corZFrameTmp_down += sum_vz_t0*sum_vz_t0;
		}
		__syncthreads();
		frameIdx += gridDim.x;
	}
	if (threadIdx.x == 0) {
        corXup[blockIdx.x] = corXFrameTmp_up;
        corXdown[blockIdx.x]=corXFrameTmp_down;
        corZup[blockIdx.x]=corZFrameTmp_up;
        corZdown[blockIdx.x]=corZFrameTmp_down;
        //printf("%f\t%f\n",corXFrameTmp_up,corXFrameTmp_down);
	}
}

int main(int argc,char** argv){ // first argument: startFrame fraction   2nd argument: endFrame
    char path[]="/home/xdengae/LJ_Fluid/DATA/traj.trr"; // *.trr file
    // prepare host memory to store the particle's info
    int numOfAtoms=0;
    int numOfFrames=0;
    getNumOfFrameAtom(path,&numOfAtoms,&numOfFrames);
    int startFrame=numOfFrames*atof(argv[1]);
    int endFrame=numOfFrames*atof(argv[2]);
    numOfFrames=endFrame-startFrame+1;
    float* vx=(float*)calloc(numOfAtoms*numOfFrames,sizeof(float));
    float* vy=(float*)calloc(numOfAtoms*numOfFrames,sizeof(float));
    float* vz=(float*)calloc(numOfAtoms*numOfFrames,sizeof(float));
    float* coor_x=(float*)calloc(numOfAtoms*numOfFrames,sizeof(float));
    float* coor_z=(float*)calloc(numOfAtoms*numOfFrames,sizeof(float));
    float* corXFrameUp=(float*)calloc(blocksPerGrid,sizeof(float));
    float* corXFrameDown=(float*)calloc(blocksPerGrid,sizeof(float));
    float* corZFrameUp=(float*)calloc(blocksPerGrid,sizeof(float));
    float* corZFrameDown=(float*)calloc(blocksPerGrid,sizeof(float));
    float H,L;
    extract_trr_file(path,vx,vy,vz,coor_x,coor_z,startFrame,endFrame,&H,&L);
    unsigned int memsize_cpu=(5*numOfAtoms*numOfFrames*sizeof(float)+4*blocksPerGrid*sizeof(float))/1000000.;
    unsigned int memsize_gpu=(4*numOfAtoms*numOfFrames*sizeof(float)+4*blocksPerGrid*sizeof(float))/1000000.;
    printf("Atom's info takes up around %d MB memory space on host RAM and %d MB memory on GPU\n",memsize_cpu,memsize_gpu);
    // copy particle info from host to device
    float* dev_vx,* dev_vz,* dev_coor_x,* dev_coor_z,* dev_corXFrameUp,* dev_corXFrameDown,* dev_corZFrameUp,* dev_corZFrameDown;
    handleError(cudaMalloc(&dev_vx,numOfAtoms*numOfFrames*sizeof(float)));
    handleError(cudaMalloc(&dev_vz,numOfAtoms*numOfFrames*sizeof(float)));
    handleError(cudaMalloc(&dev_coor_x,numOfAtoms*numOfFrames*sizeof(float)));
    handleError(cudaMalloc(&dev_coor_z,numOfAtoms*numOfFrames*sizeof(float)));
    handleError(cudaMalloc(&dev_corXFrameUp,blocksPerGrid*sizeof(float)));
    handleError(cudaMalloc(&dev_corXFrameDown,blocksPerGrid*sizeof(float)));
    handleError(cudaMalloc(&dev_corZFrameUp,blocksPerGrid*sizeof(float)));
    handleError(cudaMalloc(&dev_corZFrameDown,blocksPerGrid*sizeof(float)));
    handleError(cudaMemcpy(dev_coor_x,coor_x,numOfAtoms*numOfFrames*sizeof(float),cudaMemcpyHostToDevice));
    handleError(cudaMemcpy(dev_coor_z,coor_z,numOfAtoms*numOfFrames*sizeof(float),cudaMemcpyHostToDevice));
    handleError(cudaMemcpy(dev_vx,vx,numOfAtoms*numOfFrames*sizeof(float),cudaMemcpyHostToDevice));
    handleError(cudaMemcpy(dev_vz,vz,numOfAtoms*numOfFrames*sizeof(float),cudaMemcpyHostToDevice));
    // launch GPU test kernel
    H-=3.6;
    dim3 grid(blocksPerGrid,1);
    dim3 block(threadsPerBlock,1);
    float kx_test=PI/L;
    float kz_test=0.01;
    int deltatIdx_test=4;
    calCorrelationXKernel<<<grid,block>>>(kx_test,kz_test,numOfFrames,numOfAtoms,dev_vx,dev_vz,dev_coor_x,dev_coor_z,deltatIdx_test,dev_corXFrameUp,dev_corXFrameDown,dev_corZFrameUp,dev_corZFrameDown,H,PI/2);
    handleError(cudaMemcpy(corXFrameUp,dev_corXFrameUp,blocksPerGrid*sizeof(float),cudaMemcpyDeviceToHost));
    handleError(cudaMemcpy(corXFrameDown,dev_corXFrameDown,blocksPerGrid*sizeof(float),cudaMemcpyDeviceToHost));
    handleError(cudaMemcpy(corZFrameUp,dev_corZFrameUp,blocksPerGrid*sizeof(float),cudaMemcpyDeviceToHost));
    handleError(cudaMemcpy(corZFrameDown,dev_corZFrameDown,blocksPerGrid*sizeof(float),cudaMemcpyDeviceToHost));
    float corX_up=0.;
    float corX_down=0.;
    float corZ_up=0.;
    float corZ_down=0.;
    for(int i=0;i<blocksPerGrid;i++){
        corX_up+=corXFrameUp[i];
        corX_down+=corXFrameDown[i];
        corZ_up+=corZFrameUp[i];
        corZ_down+=corZFrameDown[i];
    }
    float corX=corX_up/corX_down;
    float corZ=corZ_up/corZ_down;
    if(calCorrelationXKernel_TEST(kx_test,kz_test,numOfFrames,numOfAtoms,vx,vz,coor_x,coor_z,deltatIdx_test,corX,corZ,H,PI/2)){
        cout<<"GPU kernel test PASS"<<endl;
    }else{
        cout<<"GPU kernel test FAIL"<<endl;
        exit(0);
    }
    // temperature analysis
    float sum_kinetic=0.;
    float* temperature_table=(float*)malloc(numOfFrames*sizeof(float));
    for(int i=0;i<numOfFrames;i++){
        for(int j=0;j<numOfAtoms;j++){
            int idx=i*numOfAtoms+j;
            sum_kinetic+=vx[idx]*vx[idx]+vy[idx]*vy[idx]+vz[idx]*vz[idx];
        }
        float tmp=sum_kinetic/3./numOfAtoms/(i+1);
        temperature_table[i]=tmp;
    }
    sum_kinetic=sum_kinetic/3./numOfAtoms/numOfFrames;
    cout<<"frames average temperature: "<<sum_kinetic<<endl;
    write_temperature_per_frame(temperature_table,numOfFrames);
    free(temperature_table);
    // start to calculate the correaltion function versus time interval deltat
    cout<<"velocity X correaltion function versus time interval : "<<endl;
    for(int i=0;i<20;i++){
        calCorrelationXKernel<<<grid,block>>>(kx_test,kz_test,numOfFrames,numOfAtoms,dev_vx,dev_vz,dev_coor_x,dev_coor_z,i,dev_corXFrameUp,dev_corXFrameDown,dev_corZFrameUp,dev_corZFrameDown,H,PI/2);
        handleError(cudaMemcpy(corXFrameUp,dev_corXFrameUp,blocksPerGrid*sizeof(float),cudaMemcpyDeviceToHost));
        handleError(cudaMemcpy(corXFrameDown,dev_corXFrameDown,blocksPerGrid*sizeof(float),cudaMemcpyDeviceToHost));
        float corX_up=0.;
        float corX_down=0.;
        for(int i=0;i<blocksPerGrid;i++){
            corX_up+=corXFrameUp[i];
            corX_down+=corXFrameDown[i];
        }
        float corX_gpu=corX_up/corX_down;
        //float corX_cpu=calCorrelationXKernel_CPU(kx_test,kz_test,numOfFrames,numOfAtoms,vx,vz,coor_x,coor_z,i,H,PI/2);
        cout<<corX_gpu<<endl;
    }
}
