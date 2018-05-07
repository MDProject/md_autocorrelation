//#include "cuda_runtime.h"
//#include "device_launch_parameters.h"
#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <fstream>

using namespace std;

#define DELTA_W 1.8//1.30254
#define DELTA_F 1.273

class Atom {
public:
	int id;
	float xc, yc, zc;
	float vx, vy, vz;
};
float max(float a,float b){
	
}
int main() {
	char path[] = "input.gro";
	float H, L, T;
	cout << "enter the H L and T:" << endl;
	cin >> H;
	cin >> L;
	cin >> T;
	H -= 1. / 2.*DELTA_F;
	int num_f_x = 2 * (int)(L / DELTA_F) + 1;
	int num_f_y = 2 * (int)(T / DELTA_F) + 1;
	int num_f_z = 2 * (int)(H / DELTA_F) + 1;
	int N_F = num_f_x*num_f_y*num_f_z;
	Atom* atom_f = (Atom*)malloc(N_F * sizeof(Atom));
	// x,y,z >=0
	int idx = 0;
	for (int nz = 0; nz <= (int)((num_f_z - 1) / 2); nz++) {
		for (int nx = 0; nx <= (int)((num_f_x - 1) / 2); nx++) {
			for (int ny = 0; ny <= (int)((num_f_y - 1) / 2); ny++) {
				atom_f[idx].id = 2;
				atom_f[idx].vx = 0; atom_f[idx].vy = 0; atom_f[idx].vz = 0;
				atom_f[idx].xc = nx*DELTA_F; atom_f[idx].yc = ny*DELTA_F; atom_f[idx].zc = nz*DELTA_F;
				idx++;
			}
		}
	}
	// duplicate the 1/8 partition to whole system
	// x>0,y<0,z>0
	for (int nz = 0; nz <= (int)((num_f_z - 1) / 2); nz++) {
		for (int nx = 0; nx <= (int)((num_f_x - 1) / 2); nx++) {
			for (int ny = 1; ny <= (int)((num_f_y - 1) / 2); ny++) {
				atom_f[idx].id = 2;
				atom_f[idx].vx = 0; atom_f[idx].vy = 0; atom_f[idx].vz = 0;
				atom_f[idx].xc = nx*DELTA_F; atom_f[idx].yc = -ny*DELTA_F; atom_f[idx].zc = nz*DELTA_F;
				idx++;
			}
		}
	}
	// x>0,y>0,z<0
	for (int nz = 1; nz <= (int)((num_f_z - 1) / 2); nz++) {
		for (int nx = 0; nx <= (int)((num_f_x - 1) / 2); nx++) {
			for (int ny = 0; ny <= (int)((num_f_y - 1) / 2); ny++) {
				atom_f[idx].id = 2;
				atom_f[idx].vx = 0; atom_f[idx].vy = 0; atom_f[idx].vz = 0;
				atom_f[idx].xc = nx*DELTA_F; atom_f[idx].yc = ny*DELTA_F; atom_f[idx].zc = -nz*DELTA_F;
				idx++;
			}
		}
	}
	// x>0,y<0,z<0
	for (int nz = 1; nz <= (int)((num_f_z - 1) / 2); nz++) {
		for (int nx = 0; nx <= (int)((num_f_x - 1) / 2); nx++) {
			for (int ny = 1; ny <= (int)((num_f_y - 1) / 2); ny++) {
				atom_f[idx].id = 2;
				atom_f[idx].vx = 0; atom_f[idx].vy = 0; atom_f[idx].vz = 0;
				atom_f[idx].xc = nx*DELTA_F; atom_f[idx].yc = -ny*DELTA_F; atom_f[idx].zc = -nz*DELTA_F;
				idx++;
			}
		}
	}
	// x<0,y>0,z>0
	for (int nz = 0; nz <= (int)((num_f_z - 1) / 2); nz++) {
		for (int nx = 1; nx <= (int)((num_f_x - 1) / 2); nx++) {
			for (int ny = 0; ny <= (int)((num_f_y - 1) / 2); ny++) {
				atom_f[idx].id = 2;
				atom_f[idx].vx = 0; atom_f[idx].vy = 0; atom_f[idx].vz = 0;
				atom_f[idx].xc = -nx*DELTA_F; atom_f[idx].yc = ny*DELTA_F; atom_f[idx].zc = nz*DELTA_F;
				idx++;
			}
		}
	}
	// x<0,y<0,z>0
	for (int nz = 0; nz <= (int)((num_f_z - 1) / 2); nz++) {
		for (int nx = 1; nx <= (int)((num_f_x - 1) / 2); nx++) {
			for (int ny = 1; ny <= (int)((num_f_y - 1) / 2); ny++) {
				atom_f[idx].id = 2;
				atom_f[idx].vx = 0; atom_f[idx].vy = 0; atom_f[idx].vz = 0;
				atom_f[idx].xc = -nx*DELTA_F; atom_f[idx].yc = -ny*DELTA_F; atom_f[idx].zc = nz*DELTA_F;
				idx++;
			}
		}
	}
	// x<0,y>0,z<0
	for (int nz = 1; nz <= (int)((num_f_z - 1) / 2); nz++) {
		for (int nx = 1; nx <= (int)((num_f_x - 1) / 2); nx++) {
			for (int ny = 0; ny <= (int)((num_f_y - 1) / 2); ny++) {
				atom_f[idx].id = 2;
				atom_f[idx].vx = 0; atom_f[idx].vy = 0; atom_f[idx].vz = 0;
				atom_f[idx].xc = -nx*DELTA_F; atom_f[idx].yc = ny*DELTA_F; atom_f[idx].zc = -nz*DELTA_F;
				idx++;
			}
		}
	}
	// x<0,y<0,z<0
	for (int nz = 1; nz <= (int)((num_f_z - 1) / 2); nz++) {
		for (int nx = 1; nx <= (int)((num_f_x - 1) / 2); nx++) {
			for (int ny = 1; ny <= (int)((num_f_y - 1) / 2); ny++) {
				atom_f[idx].id = 2;
				atom_f[idx].vx = 0; atom_f[idx].vy = 0; atom_f[idx].vz = 0;
				atom_f[idx].xc = -nx*DELTA_F; atom_f[idx].yc = -ny*DELTA_F; atom_f[idx].zc = -nz*DELTA_F;
				idx++;
			}
		}
	}
	if ((idx - 1) != N_F - 1) {
		cout << "fluid particle number is wrong" << endl;
		system("pause");
	}
	H += 1. / 2.*DELTA_F;
	int wall_layer = 3;
	int num_w_x = 2 * L / DELTA_W + 1;
	int num_w_y = 2 * T / DELTA_W + 1;
	// BOX SIZE
	float BOX_X=2*L;
	float BOX_Y=2*T;
	float BOX_Z=2*H;//+(wall_layer-1)*2*(float)DELTA_W;
	int N_W = num_w_x*num_w_y*wall_layer; // up or low
	Atom* atom_w_up = (Atom*)malloc(N_W * sizeof(Atom));
	Atom* atom_w_low = (Atom*)malloc(N_W * sizeof(Atom));
	idx = 0;
	for (int nz = 0; nz < wall_layer; nz++) {
		for (int nx = 0; nx < num_w_x; nx++) {
			for (int ny = 0; ny < num_w_y; ny++) {
				atom_w_up[idx].id = 1;
				atom_w_up[idx].vx = 0; atom_w_up[idx].vy = 0; atom_w_up[idx].vz = 0;
				atom_w_up[idx].xc = -L + nx*DELTA_W; atom_w_up[idx].yc = -T + ny*DELTA_W; atom_w_up[idx].zc = H + nz*DELTA_W;
				idx++;
			}
		}
	}
	idx = 0;
	for (int nz = 0; nz < wall_layer; nz++) {
		for (int nx = 0; nx < num_w_x; nx++) {
			for (int ny = 0; ny < num_w_y; ny++) {
				atom_w_low[idx].id = 1;
				atom_w_low[idx].vx = 0; atom_w_low[idx].vy = 0; atom_w_low[idx].vz = 0;
				atom_w_low[idx].xc = -L + nx*DELTA_W; atom_w_low[idx].yc = -T + ny*DELTA_W; atom_w_low[idx].zc = -H - nz*DELTA_W;
				idx++;
			}
		}
	}
	// face centered dot
	int num_centered_x = num_w_x - 1;
	int num_centered_y = num_w_y - 1;
	int num_centered_z = wall_layer - 1;
	Atom* atom_w_centered_up_1 = (Atom*)malloc(num_centered_x*num_centered_y*wall_layer * sizeof(Atom)); // (0,0,1) plane
	Atom* atom_w_centered_up_2 = (Atom*)malloc(num_w_y*num_centered_x*num_centered_z * sizeof(Atom)); // (0,1,0) plane
	Atom* atom_w_centered_up_3 = (Atom*)malloc(num_w_x*num_centered_y*num_centered_z * sizeof(Atom)); // (1,0,0) plane
	idx = 0;
	for (int nz = 0; nz < wall_layer; nz++) {
		for (int nx = 0; nx < num_centered_x; nx++) {
			for (int ny = 0; ny < num_centered_y; ny++) {
				atom_w_centered_up_1[idx].id = 1;
				atom_w_centered_up_1[idx].xc = -L + DELTA_W / 2. + nx*DELTA_W;
				atom_w_centered_up_1[idx].yc = -T + DELTA_W / 2. + ny*DELTA_W;
				atom_w_centered_up_1[idx].zc = H + nz*DELTA_W;
				atom_w_centered_up_1[idx].vx = 0.; atom_w_centered_up_1[idx].vy = 0.; atom_w_centered_up_1[idx].vz = 0.;
				idx++;
			}
		}
	}
	idx = 0;
	for (int nz = 0; nz < num_centered_z; nz++) {
		for (int nx = 0; nx < num_centered_x; nx++) {
			for (int ny = 0; ny < num_w_y; ny++) {
				atom_w_centered_up_2[idx].id = 1;
				atom_w_centered_up_2[idx].xc = -L + DELTA_W / 2. + nx*DELTA_W;
				atom_w_centered_up_2[idx].yc = -T + ny*DELTA_W;
				atom_w_centered_up_2[idx].zc = H + DELTA_W / 2. + nz*DELTA_W;
				atom_w_centered_up_2[idx].vx = 0.; atom_w_centered_up_2[idx].vy = 0.; atom_w_centered_up_2[idx].vz = 0.;
				idx++;
			}
		}
	}
	idx = 0;
	for (int nz = 0; nz < num_centered_z; nz++) {
		for (int nx = 0; nx < num_w_x; nx++) {
			for (int ny = 0; ny < num_centered_y; ny++) {
				atom_w_centered_up_3[idx].id = 1;
				atom_w_centered_up_3[idx].xc = -L + nx*DELTA_W;
				atom_w_centered_up_3[idx].yc = -T + DELTA_W / 2. + ny*DELTA_W;
				atom_w_centered_up_3[idx].zc = H + DELTA_W / 2. + nz*DELTA_W;
				atom_w_centered_up_3[idx].vx = 0.; atom_w_centered_up_3[idx].vy = 0.; atom_w_centered_up_3[idx].vz = 0.;
				idx++;
			}
		}
	}
	Atom* atom_w_centered_low_1 = (Atom*)malloc(num_centered_x*num_centered_y*wall_layer * sizeof(Atom)); // (0,0,1) plane
	Atom* atom_w_centered_low_2 = (Atom*)malloc(num_w_y*num_centered_x*num_centered_z * sizeof(Atom)); // (0,1,0) plane
	Atom* atom_w_centered_low_3 = (Atom*)malloc(num_w_x*num_centered_y*num_centered_z * sizeof(Atom)); // (1,0,0) plane
	idx = 0;
	for (int nz = 0; nz < wall_layer; nz++) {
		for (int nx = 0; nx < num_centered_x; nx++) {
			for (int ny = 0; ny < num_centered_y; ny++) {
				atom_w_centered_low_1[idx].id = 1;
				atom_w_centered_low_1[idx].xc = -L + DELTA_W / 2. + nx*DELTA_W;
				atom_w_centered_low_1[idx].yc = -T + DELTA_W / 2. + ny*DELTA_W;
				atom_w_centered_low_1[idx].zc = -H - nz*DELTA_W;
				atom_w_centered_low_1[idx].vx = 0.; atom_w_centered_low_1[idx].vy = 0.; atom_w_centered_low_1[idx].vz = 0.;
				idx++;
			}
		}
	}
	idx = 0;
	for (int nz = 0; nz < num_centered_z; nz++) {
		for (int nx = 0; nx < num_centered_x; nx++) {
			for (int ny = 0; ny < num_w_y; ny++) {
				atom_w_centered_low_2[idx].id = 1;
				atom_w_centered_low_2[idx].xc = -L + DELTA_W / 2. + nx*DELTA_W;
				atom_w_centered_low_2[idx].yc = -T + ny*DELTA_W;
				atom_w_centered_low_2[idx].zc = -H - DELTA_W / 2. - nz*DELTA_W;
				atom_w_centered_low_2[idx].vx = 0.; atom_w_centered_low_2[idx].vy = 0.; atom_w_centered_low_2[idx].vz = 0.;
				idx++;
			}
		}
	}
	idx = 0;
	for (int nz = 0; nz < num_centered_z; nz++) {
		for (int nx = 0; nx < num_w_x; nx++) {
			for (int ny = 0; ny < num_centered_y; ny++) {
				atom_w_centered_low_3[idx].id = 1;
				atom_w_centered_low_3[idx].xc = -L + nx*DELTA_W;
				atom_w_centered_low_3[idx].yc = -T + DELTA_W / 2. + ny*DELTA_W;
				atom_w_centered_low_3[idx].zc = -H - DELTA_W / 2. - nz*DELTA_W;
				atom_w_centered_low_3[idx].vx = 0.; atom_w_centered_low_3[idx].vy = 0.; atom_w_centered_low_3[idx].vz = 0.;
				idx++;
			}
		}
	}
	// write *.gro file
	FILE* fp = fopen(path, "w+");
	fprintf(fp, "%s\n", "MD at t=0 , input.gro");
	fprintf(fp, "%d\n", N_F + N_W * 2 + 2 * (num_centered_x*num_centered_y*wall_layer + num_w_y*num_centered_x*num_centered_z + num_w_x*num_centered_y*num_centered_z));
	// tp01 tp02 and tp03 in sequence; tp01--low wall tp02--up
	int Ntmp0 = num_centered_x*num_centered_y*wall_layer;
	for (int i = 1; i <= Ntmp0; i++) {
		fprintf(fp, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n", i, "tp01", "tp01", i, atom_w_centered_low_1[i - 1].xc, atom_w_centered_low_1[i - 1].yc, atom_w_centered_low_1[i - 1].zc, atom_w_centered_low_1[i - 1].vx, atom_w_centered_low_1[i - 1].vy, atom_w_centered_low_1[i - 1].vz);
	}
	for (int i = Ntmp0 + 1; i <= Ntmp0 + num_w_y*num_centered_x*num_centered_z; i++) {
		fprintf(fp, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n", i, "tp01", "tp01", i, atom_w_centered_low_2[i - Ntmp0 - 1].xc, atom_w_centered_low_2[i - Ntmp0 - 1].yc, atom_w_centered_low_2[i - Ntmp0 - 1].zc, atom_w_centered_low_2[i - Ntmp0 - 1].vx, atom_w_centered_low_2[i - Ntmp0 - 1].vy, atom_w_centered_low_2[i - Ntmp0 - 1].vz);
	}
	int Ntmp1 = Ntmp0 + num_w_y*num_centered_x*num_centered_z;
	for (int i = Ntmp1 + 1; i <= Ntmp1 + num_w_x*num_centered_y*num_centered_z; i++) {
		fprintf(fp, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n", i, "tp01", "tp01", i, atom_w_centered_low_3[i - Ntmp1 - 1].xc, atom_w_centered_low_3[i - Ntmp1 - 1].yc, atom_w_centered_low_3[i - Ntmp1 - 1].zc, atom_w_centered_low_3[i - Ntmp1 - 1].vx, atom_w_centered_low_3[i - Ntmp1 - 1].vy, atom_w_centered_low_3[i - Ntmp1 - 1].vz);
	}
	int Ntmp2 = Ntmp1 + num_w_x*num_centered_y*num_centered_z;
	for (int i = Ntmp2 + 1; i <= Ntmp2 + N_W; i++) {
		fprintf(fp, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n", i, "tp01", "tp01", i, atom_w_low[i - Ntmp2 - 1].xc, atom_w_low[i - Ntmp2 - 1].yc, atom_w_low[i - Ntmp2 - 1].zc, atom_w_low[i - Ntmp2 - 1].vx, atom_w_low[i - Ntmp2 - 1].vy, atom_w_low[i - Ntmp2 - 1].vz);
	}
	int Ntmp3 = Ntmp2 + N_W;
	for (int i = Ntmp3 + 1; i <= num_centered_x*num_centered_y*wall_layer + Ntmp3; i++) {
		fprintf(fp, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n", i, "tp02", "tp02", i, atom_w_centered_up_1[i - Ntmp3 - 1].xc, atom_w_centered_up_1[i - Ntmp3 - 1].yc, atom_w_centered_up_1[i - Ntmp3 - 1].zc, atom_w_centered_up_1[i - Ntmp3 - 1].vx, atom_w_centered_up_1[i - Ntmp3 - 1].vy, atom_w_centered_up_1[i - Ntmp3 - 1].vz);
	}
	int Ntmp4 = num_centered_x*num_centered_y*wall_layer + Ntmp3;
	for (int i = Ntmp4 + 1; i <= Ntmp4 + num_w_y*num_centered_x*num_centered_z; i++) {
		fprintf(fp, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n", i, "tp02", "tp02", i, atom_w_centered_up_2[i - Ntmp4 - 1].xc, atom_w_centered_up_2[i - Ntmp4 - 1].yc, atom_w_centered_up_2[i - Ntmp4 - 1].zc, atom_w_centered_up_2[i - Ntmp4 - 1].vx, atom_w_centered_up_2[i - Ntmp4 - 1].vy, atom_w_centered_up_2[i - Ntmp4 - 1].vz);
	}
	int Ntmp5 = Ntmp4 + num_w_y*num_centered_x*num_centered_z;
	for (int i = Ntmp5 + 1; i <= Ntmp5 + num_w_x*num_centered_y*num_centered_z; i++) {
		fprintf(fp, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n", i, "tp02", "tp02", i, atom_w_centered_up_3[i - Ntmp5 - 1].xc, atom_w_centered_up_3[i - Ntmp5 - 1].yc, atom_w_centered_up_3[i - Ntmp5 - 1].zc, atom_w_centered_up_3[i - Ntmp5 - 1].vx, atom_w_centered_up_3[i - Ntmp5 - 1].vy, atom_w_centered_up_3[i - Ntmp5 - 1].vz);
	}
	int Ntmp6 = Ntmp5 + num_w_x*num_centered_y*num_centered_z;
	for (int i = Ntmp6 + 1; i <= Ntmp6 + N_W; i++) {
		fprintf(fp, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n", i, "tp02", "tp02", i, atom_w_up[i - Ntmp6 - 1].xc, atom_w_up[i - Ntmp6 - 1].yc, atom_w_up[i - Ntmp6 - 1].zc, atom_w_up[i - Ntmp6 - 1].vx, atom_w_up[i - Ntmp6 - 1].vy, atom_w_up[i - Ntmp6 - 1].vz);
	}
	int Ntmp7 = Ntmp6 + N_W;
	for (int i = Ntmp7 + 1; i <= Ntmp7 + N_F; i++) {
		fprintf(fp, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n", i, "tp03", "tp03", i, atom_f[i - Ntmp7 - 1].xc, atom_f[i - Ntmp7 - 1].yc, atom_f[i - Ntmp7 - 1].zc, atom_f[i - Ntmp7 - 1].vx, atom_f[i - Ntmp7 - 1].vy, atom_f[i - Ntmp7 - 1].vz);
	}
	fprintf(fp, "%f\t%f\t%f", BOX_X, BOX_Y, BOX_Z);
	fclose(fp);
	free(atom_f);
	free(atom_w_low);
	free(atom_w_up);
	free(atom_w_centered_low_3); free(atom_w_centered_low_2); free(atom_w_centered_low_1);
	free(atom_w_centered_up_3); free(atom_w_centered_up_2); free(atom_w_centered_up_1);
	//system("pause");
}




