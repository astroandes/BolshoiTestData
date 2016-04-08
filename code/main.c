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
void dump_halo_particles(float *particle_data, float *halo_data, long long n_particles, long long n_halos, FILE **out, int halo_min, int halo_max);
FILE **open_halo_files(int halo_min, int halo_max);
void close_halo_files(FILE **out, long long n_halos);

int main(int argc, char **argv){
  long long n_halos;
  long long n_particles;
  long long n_part_files = 10; //581
  float *halo_data;
  float *part_data;
  char particle_filename[512];
  long long i;
  FILE **halo_file_list;
  float *full_part_data;
  int halo_min, halo_max;
  //  full_part_data = load_particle_data("/lustre/home/ciencias/fisica/je.forero/BolshoiParticleData/fullPM.0416", &n_particles);

  if(argc!=4){
    fprintf(stderr, "USAGE: ./a.out halo_file halo_id_min halo_id_max\n");
    exit(1);
  }

  halo_data = load_halo_data(argv[1], &n_halos);
  halo_min = atoi(argv[2]);
  halo_max = atoi(argv[3]);

  if(halo_max < 0){
    halo_max = n_halos;
  }

  halo_file_list = open_halo_files(halo_min, halo_max);

  for(i=0; i<n_part_files; i++){
    sprintf(particle_filename, "%s%04lld", 
	    "/lustre/home/ciencias/fisica/je.forero/BolshoiParticleData/subsetPM.0416.", i);
    fprintf(stdout, "%s\n", particle_filename);
    part_data = load_particle_data(particle_filename, &n_particles);
    dump_halo_particles(part_data, halo_data, n_particles, n_halos, halo_file_list, halo_min, halo_max);
    free(part_data);            
  }
  close_halo_files(halo_file_list, (halo_max - halo_min));    
  return 0;
}

FILE **open_halo_files(int halo_min, int halo_max){
  char halo_filename[512];
  FILE **out;
  long long i;
  if(!(out = malloc(sizeof(FILE *) * (halo_max - halo_min)))){
    fprintf(stderr, "problem with file allocation\n");
    exit(1);
  }

  for(i=halo_min; i<halo_max;i++){
    sprintf(halo_filename, 
	    "/lustre/home/ciencias/fisica/je.forero/BolshoiTestData/data/halo_%06lld.dat", i);
    //    fprintf(stdout, "opening file %s\n", halo_filename);
    if(!(out[i-halo_min] = fopen(halo_filename, "a"))){
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


void dump_halo_particles(float *particle_data, float *halo_data, long long n_particles, long long n_halos, FILE **out, int halo_min, int halo_max){
  long long i, j, l, id;
  long long *halo_in;
  long long n_in;
  float dist;
  float f=1.3;
  float x_part, y_part, z_part, x_halo, y_halo, z_halo, r_halo;
  float min_x, max_x, min_y, max_y, min_z, max_z;
  l = 0;
 
  min_x = min_y = min_z = 1E10;
  max_x = max_y = max_z = -1E10;

  for(i=0;i<n_particles;i++){
     x_part = 	   particle_data[i*PART_COL_NUM + PART_POS_X], 
     y_part = 	   particle_data[i*PART_COL_NUM + PART_POS_Y], 
     z_part = 	   particle_data[i*PART_COL_NUM + PART_POS_Z];
     if(x_part < min_x) min_x = x_part;
     if(y_part < min_y) min_y = y_part;
     if(z_part < min_z) min_z = z_part;

     if(x_part > max_x) max_x = x_part;
     if(y_part > max_y) max_y = y_part;
     if(z_part > max_z) max_z = z_part;
  }
  fprintf(stdout, "x[min max]: %f %f\n", min_x, max_x);
  fprintf(stdout, "y[min max]: %f %f\n", min_y, max_y);
  fprintf(stdout, "z[min max]: %f %f\n", min_z, max_z);

  /*make list of halos to be considered*/
  n_in = 0;
  if(!(halo_in = malloc(sizeof(long long) * n_halos))){
    printf("problem doing halo halo list allocation\n");
    exit(1);
  }
  for(j=0;j<n_halos;j++){
       x_halo = halo_data[j*HALO_COL_NUM + HALO_POS_X];
       y_halo = halo_data[j*HALO_COL_NUM + HALO_POS_Y];
       z_halo = halo_data[j*HALO_COL_NUM + HALO_POS_Z];
       r_halo = halo_data[j*HALO_COL_NUM + HALO_RADIUS];
       
       if( (x_halo>(min_x-2*f*r_halo))&&( x_halo< (max_x+2*f*r_halo)) 
	   && (y_halo>(min_y - 2*f*r_halo))&&( y_halo< (max_y+2*f*r_halo)) 
	   && (z_halo>(min_z -2*f*r_halo))&&( z_halo<  (max_z+2*f*r_halo))
	   && (j>=halo_min && j<halo_max)){
	 halo_in[n_in] =  j;
	 n_in++;
       }      
  }
  fprintf(stdout, "%lld halos to be considered\n", n_in);

  if(n_in>0){
    for(i=0;i<n_particles;i++){
      x_part = 	   particle_data[i*PART_COL_NUM + PART_POS_X];
      y_part = 	   particle_data[i*PART_COL_NUM + PART_POS_Y]; 
      z_part = 	   particle_data[i*PART_COL_NUM + PART_POS_Z];
      
      if(!(i%(n_particles/10))){
	l++;
	fprintf(stdout, "finished %lld percent\n", l*10);
      }
      
      for(j=0;j<n_in;j++){
	id = halo_in[j];
	x_halo = halo_data[id * HALO_COL_NUM + HALO_POS_X];
	y_halo = halo_data[id * HALO_COL_NUM + HALO_POS_Y];
	z_halo = halo_data[id * HALO_COL_NUM + HALO_POS_Z];
	r_halo = halo_data[id * HALO_COL_NUM + HALO_RADIUS];
	
	dist = (x_halo - x_part) * (x_halo - x_part) +
	  (y_halo - y_part) * (y_halo - y_part) +
	  (z_halo - z_part) * (z_halo - z_part);
	dist = sqrt(dist);
	if(dist <  1.3 * r_halo){
	  //	 fprintf(stdout, "writing to item %d\n", id-halo_min);
	  fprintf(out[id - halo_min], "%f %f %f\n", x_part, y_part, z_part);
	}
      }    
    }
  }
  free(halo_in);
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
