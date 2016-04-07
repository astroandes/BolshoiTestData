#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define HALO_POS_X 0
#define HALO_POS_Y 1
#define HALO_POS_Z 2
#define HALO_RADIUS 3
#define HALO_COL_NUM 4

#define PART_POS_X 0
#define PART_POS_Y 1
#define PART_POS_Z 2
#define PART_COL_NUM 3

long long count_lines(char *filename);
float *load_halo_data(char *filename, long long *n_points);
float *load_particle_data(char *filename, long long *n_points);
void dump_halo_particles(float *particle_data, float *halo_data, long long n_particles, long long n_halos, FILE **out);
FILE **open_halo_files(long long n_halos);
void close_halo_files(FILE **out, long long n_halos);

int main(int argc, char **argv){
  long long n_halos;
  long long n_particles;
  long long n_part_files = 10;
  float *halo_data;
  float *part_data;
  char particle_filename[512];
  long long i;
  FILE **halo_file_list;
  float *full_part_data;

  full_part_data = load_particle_data("/lustre/home/ciencias/fisica/je.forero/BolshoiParticleData/fullPM.0416", &n_particles);

  return 0;
}
  /*
  halo_data = load_halo_data(argv[1], &n_halos);

  n_halos = 1000;
  halo_file_list = open_halo_files(n_halos);

  for(i=0; i<n_part_files; i++){
    sprintf(particle_filename, "%s%04lld", 
	    "/lustre/home/ciencias/fisica/je.forero/BolshoiParticleData/subsetPM.0416.", i);
    fprintf(stdout, "%s\n", particle_filename);
    part_data = load_particle_data(particle_filename, &n_particles);
    dump_halo_particles(part_data, halo_data, n_particles, n_halos, halo_file_list);
    free(part_data);    
  }

  close_halo_files(halo_file_list, n_halos);

  return 0;
}
  */

FILE **open_halo_files(long long n_halos){
  char halo_filename[512];
  FILE **out;
  long long i;
  if(!(out = malloc(sizeof(FILE *) * n_halos))){
    fprintf(stderr, "problem with file allocation\n");
    exit(1);
  }

  for(i=0; i<n_halos;i++){
    sprintf(halo_filename, 
	    "/lustre/home/ciencias/fisica/je.forero/BolshoiTestData/data/halo_%06lld.dat", i);
    fprintf(stdout, "opening file %s\n", halo_filename);
    if(!(out[i] = fopen(halo_filename, "w"))){
      fprintf(stderr, "problem opening file %s\n", halo_filename);
      exit(1);
    }
  }
  return out;
}


void close_halo_files(FILE **out, long long n_halos){
  long long i;
  for(i=0; i<n_halos;i++){
    fclose(out[i]);
  }
}


void dump_halo_particles(float *particle_data, float *halo_data, long long n_particles, long long n_halos, FILE **out){
  long long i, j, l;
  long long found;
  float dist;
  float x_part, y_part, z_part, x_halo, y_halo, z_halo, r_halo;
  l = 0;
  for(i=0;i<n_particles;i++){
     x_part = 	   particle_data[i*PART_COL_NUM + PART_POS_X], 
     y_part = 	   particle_data[i*PART_COL_NUM + PART_POS_Y], 
     z_part = 	   particle_data[i*PART_COL_NUM + PART_POS_Z];


     if(!(i%(n_particles/10))){
       l++;
       fprintf(stdout, "finished %lld percent\n", l*10);
     }

     j = 0;
     found = 0;
     do{
       x_halo = halo_data[j*HALO_COL_NUM + HALO_POS_X];
       y_halo = halo_data[j*HALO_COL_NUM + HALO_POS_Y];
       z_halo = halo_data[j*HALO_COL_NUM + HALO_POS_Z];
       r_halo = halo_data[j*HALO_COL_NUM + HALO_RADIUS];
       dist = (x_halo - x_part) * (x_halo - x_part) +
	 (y_halo - y_part) * (y_halo - y_part) +
	 (z_halo - z_part) * (z_halo - z_part);
       dist = sqrt(dist);
       if(dist <  1.3 * r_halo){
	 found = 1;
	 fprintf(out[j], "%f %f %f\n", x_part, y_part, z_part);
       }
       j++;
     }while((j<n_halos) && (!found));
  }
}

float *load_halo_data(char *filename, long long *n_points){
  long long n;
  float *data;
  FILE *in;
  long long i;

  
  n = count_lines(filename);
  fprintf(stdout, "there are %lld halos\n", n);

  if(!(data=malloc(sizeof(float) * n * HALO_COL_NUM))){
    fprintf(stderr, "problem with data allocation\n");
    exit(1);
  }

  if(!(in=fopen(filename, "r"))){
    fprintf(stdout, "problem reading %s\n", filename);
    exit(1);
  }
  for(i=0;i<n;i++){
    fscanf(in, "%f %f %f %f\n",
	   &data[i*HALO_COL_NUM + HALO_POS_X], 
	   &data[i*HALO_COL_NUM + HALO_POS_Y], 
	   &data[i*HALO_COL_NUM + HALO_POS_Z], 
	   &data[i*HALO_COL_NUM + HALO_RADIUS]
	   );
  }
  
  *n_points = n;
  fclose(in);
  return data;
}

float *load_particle_data(char *filename, long long *n_points){
  long long n;
  float *data;
  FILE *in;
  long long i;
  long long dumb_ll;
  n = count_lines(filename);
  fprintf(stdout, "there are %lld points\n", n);

  if(!(data=malloc(sizeof(float) * n * PART_COL_NUM))){
    fprintf(stderr, "problem with data allocation\n");
    exit(1);
  }

  if(!(in=fopen(filename, "r"))){
    fprintf(stdout, "problem reading %s\n", filename);
    exit(1);
  }
  for(i=0;i<n;i++){
    fscanf(in, "%lld %f %f %f\n",
	   &dumb_ll,
	   &data[i*PART_COL_NUM + PART_POS_X], 
	   &data[i*PART_COL_NUM + PART_POS_Y], 
	   &data[i*PART_COL_NUM + PART_POS_Z]
	   );
  }
  
  *n_points = n;
  fclose(in);
  fprintf(stdout, "finished loading %lld points\n", n);
  return data;
}

long long count_lines(char *filename){
  FILE *in;
  long long n, c;
  if(!(in=fopen(filename, "r"))){
    fprintf(stdout, "problem reading %s\n", filename);
    exit(1);
  }
  n = 0;
  do{
    if(c=='\n')
      n++;
    c = fgetc(in); 
  }while(c!=EOF);
  fclose(in);
  return n;
}
