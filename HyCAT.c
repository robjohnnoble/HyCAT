// A model of expansion of two cell subpopulations dependent on a diffused resource
// Stops after 'maxloops' updates or after a cell reaches the edge of the grid

// Requires: gnuplot and libconfig
// Requires: properly formatted configuration file
// Output: data files and png images

// Author: Rob Noble

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stdbool.h>
#include <libconfig.h>

#define MAXITS 9999 // maximum iterations of diffusion equation solver before aborting
#define EPS 1.0e-4 // precision parameter of diffusion equation solver; lower value = more precise
#define GNUPLOT "gnuplot -persist"
#define INVROOT2 0.70710678118
#define ROOT2 1.41421356237
#define PI 3.14159265359

// parameters of random number generator:
#define MASK 123459876
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836

#define CELL_CLASS 0 // these index numbers are used in the "occupied" array
#define CELL_DEATHRATE 1
#define CELL_MAX_PROLIFRATE 2
#define CELL_SURVLIMIT 3
#define CELL_COMPRATE 4
#define CELL_ER 5
#define CELL_RESISTANCE 6
#define CELL_O2COEFF 7
#define CELL_DIVISIONS 8
#define CELL_STEM 9 // this is the last of the cell properties that are the same for all members of a class
#define CELL_NOW_PROLIFRATE 10 // this is the first of the cell properties that vary individually
#define CELL_NOW_PROLIFLIMIT 11
#define CELL_NOW_EXCESSDEATHRATE 12
#define CELL_NOW_SUMRATES 13
#define CELL_XCOORD 14
#define CELL_YCOORD 15
#define NUM_CELL_PROPERTIES 16 // first dimension of the "occupied" array
#define NUM_FIXED_CELL_PROPERTIES 10

// possible values of the "current" array:
#define EMPTY -5 // previously occupied by a cell that died due to the inhibitor
#define NORMAL -4 // either normal tissue or medium (which may produce oxygen etc.)
#define DEAD -3 // previously occupied by a cell that died due to hypoxia
#define MEMBRANE -2 // membrane that tumour cells can't cross
#define VESSEL -1 // vessel that produces oxygen etc.
#define TUMOUR 0 // live tumour cell

#define CONS_NORMAL 0 // for the consumption rate arrays (oxygen, estrogen, tamoxifen)
#define CONS_DEAD 1
#define CONS_MEMBRANE 2
#define CONS_VESSEL 3
#define CONS_PROLIF 4
#define CONS_QUI 5

void update(int **, int, long *, int, double **, float *, float *, float *, float *, float, int *, float **, float *, int *, int, int, int *, int **,
	int **, float **, int *, float *, float, double **, float, double **, double **, float *, int *, float *, float *, float *, float *, float *, int, float, int *, float, float, int,
	float, float, float, int, float, int, float, int, float);
void celldeath_resource(int **, double **, float **, int *, int, int *, int *, float *, int, float, int, double **, float, double **, double **, int, float *, float,
	int);
void celldeath(int **, float **, int *, int, int, int, int *, int *, float *, int, double **, double **, double **, double **, float *, int, float, int);
void cellreplacement(int **, int *, int, int, int *, int, int, int, long *, int **, int **, float **, int *, float **, float *, double **,
	double **, double **, double **, float *, int, int *, float);
void cellproliferation(int **, double **, float *, int *, int, int, int *, int, int, int, long *, int **, int **, float **, int *, int, int *, float **, int *,
	float, int, float, float *, int *, float *, float *, float *, float *, float *, float *, float *, float *, double **, double **, double **, int, float, int *, float, float, int,
	float, float, float);
void add_mutated_cells(float, int, long *, float **, int *, float *, float *, float *, float *, float *, float *, float *, float *, float *);
void mutate_cell(float **, int, float *, float *, float *, float *, float *, float *, float *, float *, int);
void set_proliferation_rate(double **, double **, double **, double **, float **, int, float *, int, int **, float, int);
void sor(double **, double **, double **, double **, double **, double **, double **, double **, int, double, int **, long *, int *, bool, float, int, float **);
void write_matrix(int **, double **, int, int, int, char *, int, FILE *, int, float *, float **, double **, double **, float, float);
void write_matrix_occupied(int **, double **, int, int, int, char *, int, FILE *, int, int, float **, int, int);
void write_matrix_stem(int **, double **, int, int, int, char *, int, FILE *, int, int, float **);
void write_matrix_maxprolif(int **, double **, int, int, int, char *, int, FILE *, int, int, float **);
void write_matrix_to_file(int **, double **, int, int, int, FILE *, float *, int);
void assign_coefficients(double **, double **, double **, double **, double **, double **, double **, int **, int, float, float, float, float,
	float *, double **, float **, float **, float *, float, bool);
void define_neighbourhood(int **, int **, float **, int *, int);
float consumption(int, float *, float);
void create_vessels(float, int **, int, float, long *, int, float **);
void create_membrane(float, int, int **);
void create_tumour(float, int, int, int **, int, int *, float **, int *, long *, float, float *, float *, float *, float *, float *, float *, float *, float *, float, int, int, FILE *);
float find_radius(int, int **);
char *preamble(int, int);
char *preamble2(int);
char *preamble3(int, int);

// start of generic functions
void generate_poisson(int, int, int *, int *, long *, int *);
void generateRandomPointAround(int *, int *, int, long *);
int inRectangle(int *, int);
int inNeighbourhood(int **, int *, int, float, int, int *, int *);
int image_to_grid(int, float);
char *zeroes(int, int);
char *concat(char *, char *);
float ran0(long *);
float maxi(float,float);
float mini(float,float);
void shuffle(int *,int,long *);
int plusmod(int, int);
int sign(float);
void shell_one(int *, int, bool);
void shell_two(float *, int *, int, bool);
int weighted_random (float *, int, long *);
float expdev(long *);
float gasdev(long *);
float truncated_log_gasdev(long *, float);
void SampleWithoutReplacement(int, int, int *, long *);
 
int main(int argc, char **argv)
{
	char* config_file = "parameters.cfg";

	long t1 = (long)time(NULL); // for timing how long the program takes to run in seconds
	int prog_timer = (int)clock(); // for timing how long the program takes to run in processor time
	int loop_t1; // for timing how parts of the program take to run in processor time
	int write_timer, update_timer, sor_timer, all_timer; // for timing how long writing to disc and the update and sor functions take to run in processor time
	int i,j,k,loop_count; // counters
	int **current; // current state of the grid
	float **occupied; // coordiates of sites currently occupied by tumour cells
	int num_occupied=0; // number of sites currently occupied by tumour cells
	int pop[2]; // number of live tumour cells in each subpopulation
	double **u; // oxygen level
	double **est; // estrogen level
	double **tam; // tamoxifen level
	double **inhibitor; // CDK inhibitor level
	double **a,**b,**c,**d,**e1,**e2,**f; // coefficients of the diffusion equation
	FILE *gp; // gnuplot pipe
	FILE *output; // subpopulation sizes
	FILE *matrix, *current_file; // record grid values at one time point for debugging
	FILE *doubling_file; // record times taken for tumour to double in size
	FILE *birthday_file;
	FILE *oxygen_gradient;
	FILE *mean_prolif_rate;
	FILE *outermost_resistant_cell_location;
	char *preamble_text = malloc(999+1); // gnuplot code needed to set up each plot
	char *preamble_text2 = malloc(999+1); // gnuplot code needed to set up each plot
	char *preamble_text3 = malloc(999+1); // gnuplot code needed to set up each plot
	int sum1,sum2; // subpopulation sizes
	float sum_gens_elapsed=0; // elapsed time in terms of cell generations
	float gens_elapsed=0;
	int stop_condition=0; // 0 = tumour has not reached the edge of the grid; 1 = tumour has reached the edge of the grid; 2 = grid appears to be static
	int n_size; // number of sites in neighbourhood
	int **n_x, **n_y; // x and y coordinates of each neigbourhood site
	float **n_p; // probability of moving to each neigbourhood site
	int *n_n; // number of cells per neighbourhood (i.e. second dimensions of n_x, n_y, n_p)
	float doubling_times[40]; // times taken for tumour to double in size
	float *birthdays;
	int num_birthdays;
	float sum_prolif_rates, max_prolif_rate;
	bool include_vessels;
	config_t cfg;
	const config_setting_t *buff_array;
	double buff; // used when reading parameter values from configuration file
	long long buff2;
	char buffer[99];
	int type_index, rate_index;
	float radius;

	//////////// cell parameters:

	float initial_propn; // initial proportion of type 0 (blue) cells
	float min_resource_survival[2]; // oxygen threshold for survival of cell types 0 and 1, relative to serum oxygen level
	float prolif_rate_limit_ratio[2]; // proliferation rate thresholds for cell types 0 and 1, relative to max proliferation rates
	float growth_rate[2]; // reproduction rates of cell types 0 and 1
	float comp_rate[2]; // cell-cell competition rates of cell types 0 and 1
	float death_rate[2]; // death rates of cell types 0 and 1 regardless of resource level
	float er_values[2]; // ER expression levels of cell types 0 and 1
	float resistance_values[2]; // inhibitor resistance of cell types 0 and 1
	float o2coeff_values[2]; // coefficents determining how proliferation rates depend on oxygen level for cell types 0 and 1
	int prolif_radius; // radius of neighbourhood for proliferation
	int replace_radius; // radius of neighbourhood for cell-cell competition
	int dim; // linear dimension of the square grid
	float membrane_radius; //dim/2-5; // radius of circular membrane that will limit tumour growth; 0 = no membrane
	int prolif_hypothesis; // 0: proliferation rates unaffected by chemical concs; 1: proliferation rate limited by estradiol and oxygen levels; 
							// 2: proliferation rate limited by estradiol but not oxygen level; 3: proliferation rate determined by inhibitor and oxygen levels
	float qui_pro_rate_ratio; // if cell's prolif rate > this fraction of max prolif rate then it is proliferative
	float stem_prob; // probability that a stem cell division will give rise to a new stem cell
	int max_divisions; // maximum number of divisions for non-stem cells
	int stem_cells_immortal; // // 0:stem cells are immortal (good for CSC models); 1:stem cells die same as others (good for non-CSC models)

	//////////// mutation parameters:

	float mutation_rate; // rate of mutation from one cell type to the other, per cell division
	float pro_mutation_rate; // probability of a mutation that changes proliferation rate
	float pro_mutation_mean; // mean of distribution of changes in proliferation rate due to mutation
	float pro_mutation_sd; // std dev of distribution of changes in proliferation rate due to mutation

	//////////// diffusion parameters:

	int solve_for_oxygen; // 0 = don't solve diffusion eqns; 1 = do solve
	int solve_for_estrogen; // 0 = don't solve diffusion eqns; 1 = do solve
	int solve_for_tamoxifen; // 0 = don't solve diffusion eqns; 1 = do solve
	int solve_for_inhibitor; // 0 = don't solve diffusion eqns; 1 = do solve

	float permeability; // vessel permeability for oxygen
	float diffusion; // oxygen diffusion rate in cm^2 s^-1
	float serum; // oxygen level in vessels or in the environment (if no vessels)
		// relative to 2.2E-7 mol cm^-3
	float cons_rate[6]; // maximum oxygen consumption rate,
		// relative to max oxygen consumption rate of proliferating cells
	
	float permeability_est; // vessel permeability for estrogen
	float diffusion_est; // estrogen diffusion rate in cm^2 s^-1
	float serum_est; // estrogen level in vessels or in the environment (if no vessels)
		// relative to 1E-15 mol cm^-3
	float cons_rate_est[6]; // maximum estrogen consumption rate,
		// relative to max oxygen consumption rate of proliferating cells

	float permeability_tam; // vessel permeability for tamoxifen
	float diffusion_tam; // tamoxifen diffusion rate in cm^2 s^-1
	float serum_tam; // tamoxifen level in vessels or in the environment (if no vessels)
		// relative to 1E-15 mol cm^-3
	float cons_rate_tam[6]; // maximum tamoxifen consumption rate,
		// relative to max oxygen consumption rate of proliferating cells

	float permeability_inhibitor; // vessel permeability for inhibitor
	float diffusion_inhibitor; // inhibitor diffusion rate in cm^2 s^-1
	float serum_inhibitor; // inhibitor level in vessels or in the environment (if no vessels)
		// relative to 1E-15 mol cm^-3
	float cons_rate_inhibitor[6]; // maximum inhibitor consumption rate,
		// relative to max oxygen consumption rate of proliferating cells

	float tam_start; // default 20
	int tam_duration;

	float inhibitor_start; // default 20
	int inhibitor_duration;

	float estrogen_curve[30];
	
	float delta; // cell diameter, in cm
	float vessel_density; // proportion of sites occupied by vessels
	int pattern; // pattern of vessel placement: -1=none (uniform resource level), 0=random Poisson disk sampling, 1=random uniform, any other number=regular
	float **vessel_size; // relative blood supply of each vessel

	//////////// other parameters:
	
	int write[6]; // number of generations before writing to disc (in order): write_matrix, write_matrix_occupied, write_matrix_stem, write_matrix_maxprolif, write_matrix_to_file, max_prolif[]
	long idum; // seed for random number generator
	int maxloops; // maximum number of updates before program terminates
	int maxage; // maximum number of generations before program terminates
	int max_colour_scale = 25; // affects colour scales
	int initial_diameter; // width of the intitial square population
	float rjac; // estimated spectral radius of the Jacobi iteration, for the diffusion equation solver; must be >0 and <1; 0.999 seems to be optimal
	float update_limit; // upper limit on number of changes to the grid before re-solving diffusion equations, relative to live tumour cell population
	int display_type; // 0 = tumour cells shaded according to age; other value = tumour cells shaded according to resource level
	int mutated; // for artificially induced mutation

	//////////// read parameters from configuration file:

	config_init(&cfg);

	if (!config_read_file(&cfg, config_file)) {
		fprintf(stderr, "%s:%d - %s\n",
			config_error_file(&cfg),
			config_error_line(&cfg),
			config_error_text(&cfg));
		config_destroy(&cfg);
		return(EXIT_FAILURE);
	}

	config_lookup_float(&cfg, "cell.initial_propn", &buff);
	initial_propn = (float)buff;
	buff_array = config_lookup(&cfg, "cell.min_resource_survival");
	for(i=0; i<2; i++) min_resource_survival[i] = config_setting_get_float_elem(buff_array, i);
	buff_array = config_lookup(&cfg, "cell.prolif_rate_limit_ratio");
	for(i=0; i<2; i++) prolif_rate_limit_ratio[i] = config_setting_get_float_elem(buff_array, i);
	buff_array = config_lookup(&cfg, "cell.growth_rate");
	for(i=0; i<2; i++) growth_rate[i] = config_setting_get_float_elem(buff_array, i);
	buff_array = config_lookup(&cfg, "cell.comp_rate");
	for(i=0; i<2; i++) comp_rate[i] = config_setting_get_float_elem(buff_array, i);
	buff_array = config_lookup(&cfg, "cell.death_rate");
	for(i=0; i<2; i++) death_rate[i] = config_setting_get_float_elem(buff_array, i);
	buff_array = config_lookup(&cfg, "cell.er_values");
	for(i=0; i<2; i++) er_values[i] = config_setting_get_float_elem(buff_array, i);
	buff_array = config_lookup(&cfg, "cell.resistance_values");
	for(i=0; i<2; i++) resistance_values[i] = config_setting_get_float_elem(buff_array, i);
	buff_array = config_lookup(&cfg, "cell.o2coeff_values");
	for(i=0; i<2; i++) o2coeff_values[i] = config_setting_get_float_elem(buff_array, i);
	config_lookup_int(&cfg, "cell.prolif_radius", &prolif_radius);
	config_lookup_int(&cfg, "cell.replace_radius", &replace_radius);
	config_lookup_int(&cfg, "cell.dim", &dim);
	config_lookup_float(&cfg, "cell.membrane_radius", &buff);
	membrane_radius = (float)buff;
	config_lookup_int(&cfg, "cell.prolif_hypothesis", &prolif_hypothesis);
	config_lookup_float(&cfg, "cell.qui_pro_rate_ratio", &buff);
	qui_pro_rate_ratio = (float)buff;
	config_lookup_float(&cfg, "cell.stem_prob", &buff);
	stem_prob = (float)buff;
	config_lookup_int(&cfg, "cell.max_divisions", &max_divisions);
	config_lookup_int(&cfg, "cell.stem_cells_immortal", &stem_cells_immortal);

	config_lookup_float(&cfg, "mutation.mutation_rate", &buff);
	mutation_rate = (float)buff;
	config_lookup_float(&cfg, "mutation.pro_mutation_rate", &buff);
	pro_mutation_rate = (float)buff;
	config_lookup_float(&cfg, "mutation.pro_mutation_mean", &buff);
	pro_mutation_mean = (float)buff;
	config_lookup_float(&cfg, "mutation.pro_mutation_sd", &buff);
	pro_mutation_sd = (float)buff;

	config_lookup_int(&cfg, "diffusion.solve_for_oxygen", &solve_for_oxygen);
	config_lookup_int(&cfg, "diffusion.solve_for_estrogen", &solve_for_estrogen);
	config_lookup_int(&cfg, "diffusion.solve_for_tamoxifen", &solve_for_tamoxifen);
	config_lookup_int(&cfg, "diffusion.solve_for_inhibitor", &solve_for_inhibitor);

	config_lookup_float(&cfg, "diffusion.permeability", &buff);
	permeability = (float)buff;
	config_lookup_float(&cfg, "diffusion.diffusion", &buff);
	diffusion = (float)buff;
	config_lookup_float(&cfg, "diffusion.serum", &buff);
	serum = (float)buff;
	buff_array = config_lookup(&cfg, "diffusion.cons_rate");
	for(i=0; i<6; i++) cons_rate[i] = config_setting_get_float_elem(buff_array, i);

	config_lookup_float(&cfg, "diffusion.permeability_est", &buff);
	permeability_est = (float)buff;
	config_lookup_float(&cfg, "diffusion.diffusion_est", &buff);
	diffusion_est = (float)buff;
	config_lookup_float(&cfg, "diffusion.serum_est", &buff);
	serum_est = (float)buff;
	buff_array = config_lookup(&cfg, "diffusion.cons_rate_est");
	for(i=0; i<6; i++) cons_rate_est[i] = config_setting_get_float_elem(buff_array, i);

	config_lookup_float(&cfg, "diffusion.permeability_tam", &buff);
	permeability_tam = (float)buff;
	config_lookup_float(&cfg, "diffusion.diffusion_tam", &buff);
	diffusion_tam = (float)buff;
	config_lookup_float(&cfg, "diffusion.serum_tam", &buff);
	serum_tam = (float)buff;
	buff_array = config_lookup(&cfg, "diffusion.cons_rate_tam");
	for(i=0; i<6; i++) cons_rate_tam[i] = config_setting_get_float_elem(buff_array, i);

	config_lookup_float(&cfg, "diffusion.permeability_inhibitor", &buff);
	permeability_inhibitor = (float)buff;
	config_lookup_float(&cfg, "diffusion.diffusion_inhibitor", &buff);
	diffusion_inhibitor = (float)buff;
	config_lookup_float(&cfg, "diffusion.serum_inhibitor", &buff);
	serum_inhibitor = (float)buff;
	buff_array = config_lookup(&cfg, "diffusion.cons_rate_inhibitor");
	for(i=0; i<6; i++) cons_rate_inhibitor[i] = config_setting_get_float_elem(buff_array, i);

	config_lookup_float(&cfg, "diffusion.tam_start", &buff);
	tam_start = (float)buff;
	config_lookup_int(&cfg, "diffusion.tam_duration", &tam_duration);
	config_lookup_float(&cfg, "diffusion.inhibitor_start", &buff);
	inhibitor_start = (float)buff;
	config_lookup_int(&cfg, "diffusion.inhibitor_duration", &inhibitor_duration);

	buff_array = config_lookup(&cfg, "diffusion.estrogen_curve");
	for(i=0; i<30; i++) estrogen_curve[i] = config_setting_get_float_elem(buff_array, i);

	config_lookup_float(&cfg, "diffusion.delta", &buff);
	delta = (float)buff;
	config_lookup_float(&cfg, "diffusion.vessel_density", &buff);
	vessel_density = (float)buff;
	config_lookup_int(&cfg, "diffusion.pattern", &pattern);

	buff_array = config_lookup(&cfg, "other.write");
	for(i=0; i<6; i++) write[i] = config_setting_get_int_elem(buff_array, i);
	config_lookup_int64(&cfg, "other.idum", &buff2);
	idum = (long)buff2;
	config_lookup_int(&cfg, "other.maxloops", &maxloops);
	config_lookup_int(&cfg, "other.maxage", &maxage);
	config_lookup_int(&cfg, "other.initial_diameter", &initial_diameter);
	config_lookup_float(&cfg, "other.rjac", &buff);
	rjac = (float)buff;
	config_lookup_float(&cfg, "other.update_limit", &buff);
	update_limit = (float)buff;
	config_lookup_int(&cfg, "other.display_type", &display_type);
	config_lookup_int(&cfg, "other.mutated", &mutated);

	config_destroy(&cfg);

	//////////// 

	int max_radius = maxi(prolif_radius,replace_radius); // max radius of neighbourhood for movement and cell-cell competition

	// max number of cells in neighbourhood for movement and cell-cell competition:
	int max_area = maxi(2*(prolif_radius+1.0/2),2*(replace_radius+1.0/2))*maxi(2*(prolif_radius+1.0/2),2*(replace_radius+1.0/2));

	preamble_text = preamble(dim, max_colour_scale); // gnuplot code needed to set up each plot
	preamble_text2 = preamble2(dim); // gnuplot code needed to set up each plot
	preamble_text3 = preamble3(dim, max_divisions); // gnuplot code needed to set up each plot
	
	gp = popen(GNUPLOT,"w"); // 'gp' is the pipe descriptor
	if (gp==NULL) {printf("Error opening pipe to GNU plot.\n"); exit(0);}

	output=fopen("output.dat","w+"); // subpopulation sizes
	fprintf(output, "Generations\tType0_pop\tType1_pop\tTotal_pop\tRadius\n");

	mean_prolif_rate=fopen("mean_prolif_rate.dat","w+");
	fprintf(mean_prolif_rate, "Generations\tMeanProlifRate\tMaxProlifRate\n");

	doubling_file=fopen("doubling_file.dat","w+");
	birthday_file=fopen("birthdays.dat","w+");
	oxygen_gradient=fopen("oxygen_gradient.dat","w+");
	outermost_resistant_cell_location = fopen("outermost_resistant_cell_location.dat","w+");
	
	//////////// memory allocation:
	
	current = malloc(dim * sizeof *current);
	if(current == NULL) {printf("Memory problem!\n"); exit(0);}
	for (i=0; i<dim; i++) {
		current[i] = malloc(dim * sizeof *current[i]);
		if(current[i] == NULL) {printf("Memory problem!\n"); exit(0);}
		}

	occupied = malloc(NUM_CELL_PROPERTIES * sizeof *occupied);
	if(occupied == NULL) {printf("Memory problem!\n"); exit(0);}
	for (i=0; i<NUM_CELL_PROPERTIES; i++) {
		occupied[i] = malloc(dim*dim * sizeof *occupied[i]);
		if(occupied[i] == NULL) {printf("Memory problem!\n"); exit(0);}
		}

	u = malloc(dim * sizeof *u);
	if(u == NULL) {printf("Memory problem!\n"); exit(0);}
	for (i=0; i<dim; i++) {
		u[i] = malloc(dim * sizeof *u[i]);
		if(u[i] == NULL) {printf("Memory problem!\n"); exit(0);}
		}

	est = malloc(dim * sizeof *est);
	if(est == NULL) {printf("Memory problem!\n"); exit(0);}
	for (i=0; i<dim; i++) {
		est[i] = malloc(dim * sizeof *est[i]);
		if(est[i] == NULL) {printf("Memory problem!\n"); exit(0);}
		}

	tam = malloc(dim * sizeof *tam);
	if(tam == NULL) {printf("Memory problem!\n"); exit(0);}
	for (i=0; i<dim; i++) {
		tam[i] = malloc(dim * sizeof *tam[i]);
		if(tam[i] == NULL) {printf("Memory problem!\n"); exit(0);}
		}

	inhibitor = malloc(dim * sizeof *inhibitor);
	if(inhibitor == NULL) {printf("Memory problem!\n"); exit(0);}
	for (i=0; i<dim; i++) {
		inhibitor[i] = malloc(dim * sizeof *inhibitor[i]);
		if(inhibitor[i] == NULL) {printf("Memory problem!\n"); exit(0);}
		}

	a = malloc(dim * sizeof *a);
	if(a == NULL) {printf("Memory problem!\n"); exit(0);}
	for (i=0; i<dim; i++) {
		a[i] = malloc(dim * sizeof *a[i]);
		if(a[i] == NULL) {printf("Memory problem!\n"); exit(0);}
		}

	b = malloc(dim * sizeof *b);
	if(b == NULL) {printf("Memory problem!\n"); exit(0);}
	for (i=0; i<dim; i++) {
		b[i] = malloc(dim * sizeof *b[i]);
		if(b[i] == NULL) {printf("Memory problem!\n"); exit(0);}
		}

	c = malloc(dim * sizeof *c);
	if(c == NULL) {printf("Memory problem!\n"); exit(0);}
	for (i=0; i<dim; i++) {
		c[i] = malloc(dim * sizeof *c[i]);
		if(c[i] == NULL) {printf("Memory problem!\n"); exit(0);}
		}

	d = malloc(dim * sizeof *d);
	if(d == NULL) {printf("Memory problem!\n"); exit(0);}
	for (i=0; i<dim; i++) {
		d[i] = malloc(dim * sizeof *d[i]);
		if(d[i] == NULL) {printf("Memory problem!\n"); exit(0);}
		}

	e1 = malloc(dim * sizeof *e1);
	if(e1 == NULL) {printf("Memory problem!\n"); exit(0);}
	for (i=0; i<dim; i++) {
		e1[i] = malloc(dim * sizeof *e1[i]);
		if(e1[i] == NULL) {printf("Memory problem!\n"); exit(0);}
		}

	e2 = malloc(dim * sizeof *e2);
	if(e2 == NULL) {printf("Memory problem!\n"); exit(0);}
	for (i=0; i<dim; i++) {
		e2[i] = malloc(dim * sizeof *e2[i]);
		if(e2[i] == NULL) {printf("Memory problem!\n"); exit(0);}
		}

	f = malloc(dim * sizeof *f);
	if(f == NULL) {printf("Memory problem!\n"); exit(0);}
	for (i=0; i<dim; i++) {
		f[i] = malloc(dim * sizeof *f[i]);
		if(f[i] == NULL) {printf("Memory problem!\n"); exit(0);}
		}

	n_x = malloc(max_radius * sizeof *n_x);
	if(n_x == NULL) {printf("Memory problem!\n"); exit(0);}
	for (i=0; i<max_radius; i++) {
		n_x[i] = malloc(max_area * sizeof *n_x[i]);
		if(n_x[i] == NULL) {printf("Memory problem!\n"); exit(0);}
		}

	n_y = malloc(max_radius * sizeof *n_y);
	if(n_y == NULL) {printf("Memory problem!\n"); exit(0);}
	for (i=0; i<max_radius; i++) {
		n_y[i] = malloc(max_area * sizeof *n_y[i]);
		if(n_y[i] == NULL) {printf("Memory problem!\n"); exit(0);}
		}

	n_p = malloc(max_radius * sizeof *n_p);
	if(n_p == NULL) {printf("Memory problem!\n"); exit(0);}
	for (i=0; i<max_radius; i++) {
		n_p[i] = malloc(max_area * sizeof *n_p[i]);
		if(n_p[i] == NULL) {printf("Memory problem!\n"); exit(0);}
		}

	n_n = malloc(max_radius * sizeof *n_n);
	if(n_n == NULL) {printf("Memory problem!\n"); exit(0);}

	birthdays = malloc(dim*dim * sizeof *birthdays);
	if(birthdays == NULL) {printf("Memory problem!\n"); exit(0);}

	vessel_size = malloc(dim * sizeof *vessel_size);
	if(vessel_size == NULL) {printf("Memory problem!\n"); exit(0);}
	for (i=0; i<dim; i++) {
		vessel_size[i] = malloc(dim * sizeof *vessel_size[i]);
		if(vessel_size[i] == NULL) {printf("Memory problem!\n"); exit(0);}
		}
	
	printf("\nMemory allocated OK\n");
	
	//////////// initial conditions:

	define_neighbourhood(n_x, n_y, n_p, n_n, max_radius);

	for(i=0;i<dim;i++) for(j=0;j<dim;j++) current[i][j]=NORMAL; // empty grid
	num_occupied=0;
	pop[0]=0;
	pop[1]=0;
	for(i=0;i<40;i++) doubling_times[i]=-1; // -1 translates to "NA"

	if(membrane_radius>0) create_membrane(membrane_radius, dim, current); // assign membrance cells to sites

	if(pattern>=0) {
		include_vessels=true;
		create_vessels(membrane_radius,current,dim,vessel_density,&idum,pattern,vessel_size); // assign vessels to sites
		}
	else include_vessels=false;

	create_tumour(initial_propn, dim, initial_diameter, current, max_colour_scale, pop, occupied, &num_occupied, &idum, membrane_radius, death_rate,
		 growth_rate, min_resource_survival, prolif_rate_limit_ratio, comp_rate, er_values, resistance_values, o2coeff_values, stem_prob, max_divisions, mutated, outermost_resistant_cell_location); // assign tumour cells to sites

	for(i=0;i<dim;i++) for(j=0;j<dim;j++) u[i][j]=1; // intitial oxygen levels; important if solve_for_oxygen = 0

	// find equilibrium oxygen concentrations and proliferation and death rates:
	if(solve_for_oxygen != 0) for(i=0; i<5; i++) {
		// assign coefficients for oxygen:
		assign_coefficients(a,b,c,d,e1,e2,f,current,dim,permeability,diffusion,serum,delta,cons_rate,u,occupied,vessel_size,estrogen_curve,sum_gens_elapsed,false);
		// solve diffusion equations for oxygen:
		sor(a,b,c,d,e1,e2,f,u,dim,rjac,current,&idum,&stop_condition,include_vessels,serum,num_occupied,occupied);
		for(j=0;j<num_occupied;j++) set_proliferation_rate(u, est, tam, inhibitor, occupied, j, prolif_rate_limit_ratio, prolif_hypothesis, current, qui_pro_rate_ratio, 0);
		}

	printf("Initialised\n\n");

	//////////// main loop:

	for(loop_count=0;loop_count<maxloops;loop_count++) {

		all_timer=(int)clock();

		loop_t1 = (int)clock()-prog_timer;
		printf("Update %d, Cell time %.1f, Elapsed %lds, %.1f megaticks, ",loop_count, sum_gens_elapsed, (long)time(NULL)-t1, (float)loop_t1/1E6);

		sor_timer=(int)clock();

		if(solve_for_oxygen != 0) {
			// assign coefficients for oxygen:
			assign_coefficients(a,b,c,d,e1,e2,f,current,dim,permeability,diffusion,serum,delta,cons_rate,u,occupied,vessel_size,estrogen_curve,sum_gens_elapsed,false);
			for(i=0;i<dim;i++) for(j=0;j<dim;j++) u[i][j]=0; // intitial diffusion levels
			// solve diffusion equations for oxygen:
			sor(a,b,c,d,e1,e2,f,u,dim,rjac,current,&idum,&stop_condition,include_vessels,serum,num_occupied,occupied);
			}

		// exit if there is a problem with solving the diffusion equations:
		if(stop_condition==3) {
			printf("\nMAXITS exceeded (problem solving oxygen diffusion equations)\n");
			break;
			}

		if(solve_for_estrogen != 0) {
			// assign coefficients for estrogen:
			assign_coefficients(a,b,c,d,e1,e2,f,current,dim,permeability_est,diffusion_est,serum_est,delta,cons_rate_est,u,occupied,vessel_size,estrogen_curve,sum_gens_elapsed,true);
			for(i=0;i<dim;i++) for(j=0;j<dim;j++) est[i][j]=0; // intitial diffusion levels
			// solve diffusion equations for estrogen:
			if(solve_for_estrogen != 0) sor(a,b,c,d,e1,e2,f,est,dim,rjac,current,&idum,&stop_condition,include_vessels,serum_est,num_occupied,occupied);
			}

		// exit if there is a problem with solving the diffusion equations:
		if(stop_condition==3) {
			printf("\nMAXITS exceeded (problem solving estrogen diffusion equations)\n");
			break;
			}

		if(solve_for_tamoxifen != 0) {
			// assign coefficients for tamoxifen:
			assign_coefficients(a,b,c,d,e1,e2,f,current,dim,permeability_tam,diffusion_tam,serum_tam,delta,cons_rate_tam,u,occupied,vessel_size,estrogen_curve,sum_gens_elapsed,false);
			for(i=0;i<dim;i++) for(j=0;j<dim;j++) tam[i][j]=0; // intitial diffusion levels
			// solve diffusion equations for tamoxifen:
			if(sum_gens_elapsed>=tam_start && sum_gens_elapsed<=tam_start+tam_duration) sor(a,b,c,d,e1,e2,f,tam,dim,rjac,current,&idum,&stop_condition,
				include_vessels,serum_tam,num_occupied,occupied);
			}

		// exit if there is a problem with solving the diffusion equations:
		if(stop_condition==3) {
			printf("\nMAXITS exceeded (problem solving tamoxifen diffusion equations)\n");
			break;
			}

		if(solve_for_inhibitor != 0 && serum_inhibitor > 0) {
			// assign coefficients for inhibitor:
			assign_coefficients(a,b,c,d,e1,e2,f,current,dim,permeability_inhibitor,diffusion_inhibitor,serum_inhibitor,delta,cons_rate_inhibitor,u,occupied,vessel_size,estrogen_curve,sum_gens_elapsed,false);
			for(i=0;i<dim;i++) for(j=0;j<dim;j++) inhibitor[i][j]=0; // intitial diffusion levels
			// solve diffusion equations for inhibitor:
			if(sum_gens_elapsed>=inhibitor_start && sum_gens_elapsed<=inhibitor_start+inhibitor_duration) sor(a,b,c,d,e1,e2,f,inhibitor,dim,rjac,current,&idum,&stop_condition,
				include_vessels,serum_inhibitor,num_occupied,occupied);
			}

		// exit if there is a problem with solving the diffusion equations:
		if(stop_condition==3) {
			printf("\nMAXITS exceeded (problem solving inhibitor diffusion equations)\n");
			break;
			}

		sor_timer=(int)clock()-sor_timer;

		// set proliferation rates and proliferation rate limits:
		if(num_occupied < 10) for(i=0;i<num_occupied;i++) set_proliferation_rate(u, est, tam, inhibitor, occupied, i, prolif_rate_limit_ratio, prolif_hypothesis, current, qui_pro_rate_ratio, 1);
		else for(i=0;i<num_occupied;i++) set_proliferation_rate(u, est, tam, inhibitor, occupied, i, prolif_rate_limit_ratio, prolif_hypothesis, current, qui_pro_rate_ratio, 0);

		// write current[][] matrix to gnuplot:
		write_timer=(int)clock();
		if(!(loop_count%write[0])) write_matrix(current,u,dim,loop_count,maxloops,preamble_text,max_colour_scale,gp,display_type,min_resource_survival,occupied,tam,inhibitor,gens_elapsed,sum_gens_elapsed);
		if(!(loop_count%write[1])) for(rate_index = 0; rate_index <= 1; rate_index++) for(type_index = 0; type_index <= 2; type_index++) write_matrix_occupied(current,u,dim,loop_count,maxloops,preamble_text2,max_colour_scale,gp,display_type,num_occupied,occupied,rate_index,type_index);
		if(!(loop_count%write[2])) write_matrix_stem(current,u,dim,loop_count,maxloops,preamble_text3,max_colour_scale,gp,display_type,num_occupied,occupied);
		if(!(loop_count%write[3])) write_matrix_maxprolif(current,u,dim,loop_count,maxloops,preamble_text2,max_colour_scale,gp,display_type,num_occupied,occupied);
		if(!(loop_count%write[4])) write_matrix_to_file(current, u, dim, loop_count, max_colour_scale, matrix, min_resource_survival, maxloops);
		write_timer=(int)clock()-write_timer;

		// write subpopulation sizes to file:
		radius = find_radius(dim, current);
		fprintf(output,"%f\t%d\t%d\t%d\t%f\n",sum_gens_elapsed,pop[0],pop[1],num_occupied,radius);

		// write mean and max proliferation rates to file:
		sum_prolif_rates = 0;
		max_prolif_rate = -1;
		for(i=0; i<num_occupied; i++) {
			sum_prolif_rates += occupied[CELL_MAX_PROLIFRATE][i];
			if(occupied[CELL_MAX_PROLIFRATE][i] > max_prolif_rate) max_prolif_rate = occupied[CELL_MAX_PROLIFRATE][i];
			}
		fprintf(mean_prolif_rate, "%f\t%f\t%f\n", sum_gens_elapsed, sum_prolif_rates / num_occupied, max_prolif_rate);

		if(!(loop_count%write[5])) { // write all max prolif rates to file
			sprintf(buffer, "max_prolif_rates/max_prolif_%s%d.dat", zeroes(loop_count,maxloops), loop_count);
			current_file = fopen(buffer, "w+");
			fprintf(current_file, "Generations\tx\ty\tMaxProlifRate\n");
			for(i=0;i<num_occupied;i++) fprintf(current_file, "%f\t%d\t%d\t%.4f\n", sum_gens_elapsed, (int)occupied[CELL_XCOORD][i], (int)occupied[CELL_YCOORD][i], occupied[CELL_MAX_PROLIFRATE][i]);
			}

		num_birthdays = 0;

		// update values of current[][]:
		update_timer=(int)clock();
		update(current, dim, &idum, max_colour_scale, u, min_resource_survival, prolif_rate_limit_ratio, growth_rate, death_rate, update_limit, &num_occupied,
			occupied, &sum_gens_elapsed, pop, prolif_radius, replace_radius, &stop_condition, n_x, n_y, n_p, n_n, doubling_times, membrane_radius,
			est, serum, tam, inhibitor, birthdays, &num_birthdays, &gens_elapsed, er_values, resistance_values, o2coeff_values, comp_rate, prolif_hypothesis, mutation_rate, &mutated, qui_pro_rate_ratio,
			stem_prob, max_divisions, pro_mutation_rate, pro_mutation_mean, pro_mutation_sd, stem_cells_immortal, tam_start, tam_duration, inhibitor_start, inhibitor_duration, initial_propn);
		update_timer=(int)clock()-update_timer;

		// write number of birthdays to file:
		printf("num_birthdays = %d\n",num_birthdays);
		for(i=0; i<num_birthdays; i++) fprintf(birthday_file, "%f\t", birthdays[i]);

		// exit if any cell has reached the edge:
		if(stop_condition==1) {
			printf("\nHit edge!\n\n");
			break;
			}

		// exit if grid appears to be static:
		if(stop_condition==2) {
			printf("\nReached loops limit (grid appears to be static).\n\n");
			break;
			}

		if(num_occupied==0) { // exit if no live cells remain
			printf("\nNo live cells.\n\n");
			break;
			}

		if(sum_gens_elapsed >= maxage) { // exit if number of generations has reached the limit
			printf("\nReached generations limit.\n\n");
			break;
			}

		all_timer=(int)clock()-all_timer;
		printf("sor took %.1f%%, writing took %.1f%%, update took %.1f%%\n",100*(float)sor_timer/all_timer,100*(float)write_timer/all_timer,100*(float)update_timer/all_timer);

		}
	
	if(loop_count==maxloops) printf("\nReached loops limit.\n\n");

	
	// write to disc when finished:
	radius = find_radius(dim, current);
	fprintf(output,"%f\t%d\t%d\t%d\t%f\n",sum_gens_elapsed,pop[0],pop[1],num_occupied,radius);
	write_matrix(current,u,dim,loop_count,maxloops,preamble_text,max_colour_scale,gp,display_type,min_resource_survival,occupied,tam,inhibitor,gens_elapsed,sum_gens_elapsed);
	for(rate_index = 0; rate_index <= 1; rate_index++) for(type_index = 0; type_index <= 2; type_index++) write_matrix_occupied(current,u,dim,loop_count,maxloops,preamble_text2,max_colour_scale,gp,display_type,num_occupied,occupied,rate_index,type_index);
	write_matrix_stem(current,u,dim,loop_count,maxloops,preamble_text3,max_colour_scale,gp,display_type,num_occupied,occupied);
	write_matrix_maxprolif(current,u,dim,loop_count,maxloops,preamble_text2,max_colour_scale,gp,display_type,num_occupied,occupied);
	write_matrix_to_file(current, u, dim, loop_count, max_colour_scale, matrix, min_resource_survival, maxloops);
	for(i=0;i<39;i++) {
		if(doubling_times[i+1]==-1) fprintf(doubling_file, "NA\n");
		else  fprintf(doubling_file, "%f\n", doubling_times[i+1]-doubling_times[i]);
		}
	for(i=0;i<dim;i++) fprintf(oxygen_gradient, "%f\t", u[i][(int)(dim/2)]);
	
	fclose(gp);
	fclose(matrix);
	fclose(current_file);
	fclose(output);
	fclose(birthday_file);
	
	printf("Completed in %ld seconds.\n",(long)time(NULL)-t1);
	
	return 0;
}

void update(int **current, int dim, long *idum, int max_colour_scale, double **u, float *min_resource_survival, float *prolif_rate_limit_ratio, float *growth_rate,
 float *death_rate, float update_limit, int *num_occupied, float **occupied, float *sum_gens_elapsed, int *pop, int prolif_radius, int replace_radius,
  int *stop_condition, int **n_x, int **n_y, float **n_p, int *n_n, float *doubling_times, float membrane_radius, double **est,
   float serum, double **tam, double **inhibitor, float *birthdays, int *num_birthdays, float *gens_elapsed, float *er_values, float *resistance_values, float *o2coeff_values, float *comp_rate,
   int prolif_hypothesis, float mutation_rate, int *mutated, float qui_pro_rate_ratio, float stem_prob, int max_divisions, float pro_mutation_rate,
   float pro_mutation_mean, float pro_mutation_sd, int stem_cells_immortal, float tam_start, int tam_duration, float inhibitor_start, int inhibitor_duration, float initial_propn)
{
	// updates values of current[][]

	int n_start;
	float event_type;
	int i,j,k,it,loops,site,subpop;
	int neighbour;
	int loops_limit=100*dim*dim; // maximum number of update loops (just in case the grid ends up static)
	int cell_num;
	int occ=*num_occupied;
	int num_bir = *num_birthdays;
	int event_count=0;
	int xnew,ynew;
	int subtype;
	int initial_pop=pop[0]+pop[1]; // needed for tracking how many sites have changed value as a proportion of initial population
	float bar1,bar2,bar3,barsum; // thresholds for the event types
	int x, y;
	float sumsums;
	int treatment_on = 0, treatment_on_init = 0;
	
	loops = 0;
	*gens_elapsed = 0;

	// if tamixofen model:
	if(prolif_hypothesis == 1 || prolif_hypothesis == 2) {
		if(*sum_gens_elapsed == tam_start) {*sum_gens_elapsed = 1.000001 * *sum_gens_elapsed; return;}
		if(*sum_gens_elapsed == tam_start + tam_duration) {*sum_gens_elapsed = 1.000001 * *sum_gens_elapsed; return;}
		if(*sum_gens_elapsed >= tam_start && *sum_gens_elapsed < tam_start + tam_duration) treatment_on_init = 1;
	}
	// if CDK inhibitor model:
	if(prolif_hypothesis == 3) {
		if(*sum_gens_elapsed == inhibitor_start) {*sum_gens_elapsed = 1.000001 * *sum_gens_elapsed; return;}
		if(*sum_gens_elapsed == inhibitor_start + inhibitor_duration) {*sum_gens_elapsed = 1.000001 * *sum_gens_elapsed; return;}
		if(*sum_gens_elapsed >= inhibitor_start && *sum_gens_elapsed < inhibitor_start + inhibitor_duration) treatment_on_init = 3;

	// death due to insufficient resources:
	celldeath_resource(current, u, occupied, pop, subtype, &event_count, &occ, min_resource_survival, initial_pop, update_limit, max_colour_scale, est,
		serum, tam, inhibitor, prolif_hypothesis, prolif_rate_limit_ratio, qui_pro_rate_ratio, stem_cells_immortal);
	}

	if (event_count<update_limit*initial_pop) do{

		sumsums = 0; // sumsums is used to calculate time between events
		for(i=0; i<occ; i++) sumsums += occupied[CELL_NOW_SUMRATES][i];

		// time passing: add number of generations elapsed in this loop
		// (i.e. generations of a cell line with proliferation rate 1, death rate 0, and competition rate 0)
		*gens_elapsed += expdev(idum)/sumsums; // Gillespie algorithm

		// if tamixofen model:
		if(prolif_hypothesis == 1 || prolif_hypothesis == 2) {
			if(*sum_gens_elapsed + *gens_elapsed >= tam_start && *sum_gens_elapsed + *gens_elapsed < tam_start + tam_duration) treatment_on = 1;
			else treatment_on = 0;
		}
		// if CDK inhibitor model:
		if(prolif_hypothesis == 3) {
			if(*sum_gens_elapsed + *gens_elapsed >= inhibitor_start && *sum_gens_elapsed + *gens_elapsed < inhibitor_start + inhibitor_duration) treatment_on = 3;
			else treatment_on = 0;
		}

		if(treatment_on != treatment_on_init) {
			//printf("treatment_on switched from %d to %d\n", treatment_on_init, treatment_on);
			if(treatment_on - treatment_on_init == 1) *sum_gens_elapsed = tam_start;
			else if(treatment_on - treatment_on_init == -1) *sum_gens_elapsed = tam_start + tam_duration;
			if(treatment_on - treatment_on_init == 3) *sum_gens_elapsed = inhibitor_start;
			else if(treatment_on - treatment_on_init == -3) *sum_gens_elapsed = inhibitor_start + inhibitor_duration;
			
			*num_birthdays = num_bir;
			*num_occupied=occ;

			// mutate randomly selected cells from susceptible to resistant:
			if(*mutated == 1 && treatment_on > 0) add_mutated_cells(initial_propn, *num_occupied, idum, occupied, pop, &sumsums, er_values, resistance_values, o2coeff_values, 
				min_resource_survival, prolif_rate_limit_ratio, growth_rate, death_rate, comp_rate);

			return; // exit this function and resolve diffusion equations when drug is about to switch on or off
		}

		cell_num = weighted_random(occupied[CELL_NOW_SUMRATES], occ, idum); // select a cell at random
		
		x = occupied[CELL_XCOORD][cell_num]; // x coordinate of selected cell
		y = occupied[CELL_YCOORD][cell_num]; // y coordinate of selected cell
		bar1 = occupied[CELL_DEATHRATE][cell_num] + occupied[CELL_NOW_EXCESSDEATHRATE][cell_num]; // relative risk of dying
		bar2 = occupied[CELL_COMPRATE][cell_num]; // relative risk of competing
		bar3 = occupied[CELL_NOW_PROLIFRATE][cell_num]; // relative risk of proliferating
		barsum = bar1 + bar2 + bar3;

		event_type=ran0(idum)*barsum; // determines whether event is death, birth, predation, or nothing

		subtype=occupied[CELL_CLASS][cell_num];

		// death regardless of resource level:
		if(event_type < bar1) celldeath(current, occupied, pop, x, y, cell_num, &event_count, &occ, &sumsums, 0, u, est, tam, inhibitor, prolif_rate_limit_ratio,
			prolif_hypothesis, qui_pro_rate_ratio, stem_cells_immortal);

		// replacement (cell-cell competition):
		else if(event_type < bar1+bar2) cellreplacement(current, pop, x, y, &event_count, max_colour_scale, dim,
			replace_radius, idum, n_x, n_y, n_p, n_n, occupied, &sumsums, u, est, tam, inhibitor, prolif_rate_limit_ratio, prolif_hypothesis, &occ, qui_pro_rate_ratio);

		// reproduction:
		else cellproliferation(current, u, prolif_rate_limit_ratio, pop, x, y, &event_count, max_colour_scale, dim, prolif_radius, idum, n_x, n_y, n_p, n_n, subtype,
			 stop_condition, occupied, &occ, membrane_radius, cell_num, *sum_gens_elapsed, birthdays, &num_bir, &sumsums, er_values, resistance_values, o2coeff_values, min_resource_survival,
			 growth_rate, death_rate, comp_rate, est, tam, inhibitor, prolif_hypothesis, mutation_rate, mutated, qui_pro_rate_ratio, stem_prob, max_divisions,
			 pro_mutation_rate, pro_mutation_mean, pro_mutation_sd);

		if(sumsums == 0) { // if all event rates are zero - i.e. if tumour is entirely dormant - then just progress the timer
			*sum_gens_elapsed += 1; // the size of this increment is arbritrary
			*num_birthdays = num_bir;
			*num_occupied=occ;
			return;
			}

		// record doubling times:
		if(log2(pop[1])==(int)log2(pop[1]) && log2(pop[1])<40 && doubling_times[(int)log2(pop[1])-1]==-1) doubling_times[(int)log2(pop[1])-1]=*sum_gens_elapsed;

		if((*stop_condition)==1) {
			*num_occupied=occ;
			return; // tumour has reached edge of grid
			}

		loops++;
		if(loops>=loops_limit) {
			*stop_condition=2; // grid appears to be static
			*num_occupied=occ;
			return;
			}

		}while(event_count<update_limit*initial_pop); // diffusion equations must be re-solved after the grid has substantially changed

	*num_birthdays = num_bir;
	*num_occupied=occ;
	*sum_gens_elapsed += *gens_elapsed;
}

void celldeath_resource(int **current, double **u, float **occupied, int *pop, int subtype, int *event_count, int *occ, float *min_resource_survival,
 int initial_pop, float update_limit, int max_colour_scale, double **est, float serum, double **tam, double **inhibitor, int prolif_hypothesis, float *prolif_rate_limit_ratio,
 float qui_pro_rate_ratio, int stem_cells_immortal)
{
	int deathlist[*occ];
	float deathlist_val[*occ];
	int deathlist_now[*occ];
	float deathlist_now_val[*occ];
	int deathlist_length = 0;
	int deathevents;
	int x, y;
	int cell_num, cell_to_kill;
	float prolif_energy, balance;
	float energy[2];
	int kill_it;

	for(cell_num=0;cell_num<(*occ);cell_num++) {
		x = occupied[CELL_XCOORD][cell_num]; // x coordinate of selected cell
		y = occupied[CELL_YCOORD][cell_num]; // y coordinate of selected cell
		subtype = occupied[CELL_CLASS][cell_num];
		
		deathlist_val[cell_num] = u[x][y];
		deathlist[cell_num] = cell_num;

		kill_it = 0;
		
		// if cell is trying to proliferate faster than the oxygen supply will allow:
		if(prolif_hypothesis == 2 && occupied[CELL_NOW_PROLIFRATE][cell_num] > occupied[CELL_NOW_PROLIFLIMIT][cell_num]) {
			kill_it = 1;
			deathlist_val[cell_num] = 0; // these cells will be killed first
			printf("\nToo fast (%.2f,\t%.2f)\n",occupied[CELL_NOW_PROLIFRATE][cell_num], occupied[CELL_NOW_PROLIFLIMIT][cell_num]);
			}
		
		// if cell has insufficient oxygen to survive:
		if(deathlist_val[cell_num] < occupied[CELL_SURVLIMIT][cell_num]) kill_it = 1;

		if(kill_it == 1) deathlist_length++;
		
		}
	shell_two(deathlist_val, deathlist, *occ, true); // sorts deathlist_val and deathlist in ascending order of deathlist_val

	deathevents = mini(deathlist_length, maxi(update_limit*initial_pop, 1)); // number of cells to be killed
	(*event_count) = deathevents;

	shell_one(deathlist, deathevents, false); // sorts deathlist (indices < deathevents) in descending order
	// this ensures that when slots are removed from occupied[][], they don't affect the pointers stored in deathlist

	printf("deathevents = %d (%d, %f)\n", deathevents, deathlist_length, update_limit*initial_pop);
	
	for(cell_num = 0; cell_num < deathevents; cell_num++) {
		cell_to_kill = deathlist[cell_num];
		x = occupied[CELL_XCOORD][cell_to_kill]; // x coordinate of selected cell
		y = occupied[CELL_YCOORD][cell_to_kill]; // y coordinate of selected cell
		celldeath(current, occupied, pop, x, y, cell_to_kill, event_count, occ, 0, 1, u, est, tam, inhibitor, prolif_rate_limit_ratio, prolif_hypothesis, qui_pro_rate_ratio,
			stem_cells_immortal);
		}
}

void celldeath(int **current, float **occupied, int *pop, int x, int y, int cell_to_kill, int *event_count, int *occ, float *sumsums, int death_due_to_resource,
	double **u, double **est, double **tam, double **inhibitor, float *prolif_rate_limit_ratio, int prolif_hypothesis, float qui_pro_rate_ratio, int stem_cells_immortal)
{
	int i,j;
	int x_endoflist = occupied[CELL_XCOORD][(*occ)-1];
	int y_endoflist = occupied[CELL_YCOORD][(*occ)-1];
	int subtype = occupied[CELL_CLASS][cell_to_kill];
	int xneigh, yneigh;
	int neigh_cell_num;

	if(death_due_to_resource == 0 && occupied[CELL_STEM][cell_to_kill] == 1 && stem_cells_immortal == 1) return; // stem cells are immortal

	if(death_due_to_resource == 0) (*sumsums) -= occupied[CELL_NOW_SUMRATES][cell_to_kill];
	// subtract from the sum of all sums of rates, unless this function is being called from celldeath_resource

	for(i = 0; i < NUM_CELL_PROPERTIES; i++) occupied[i][cell_to_kill]=occupied[i][(*occ)-1];
	current[x_endoflist][y_endoflist] = current[x][y];
	(*occ)--;

	if(death_due_to_resource == 0) current[x][y] = DEAD; // necrotic cells remain on the grid
	else current[x][y] = EMPTY; // cells killed by inhibitor are removed

	(*event_count)++;
	pop[subtype]--;

	// set proliferation rate for 8 nearest neighbours:
	for(xneigh = x-1; xneigh <= x+1; xneigh++) for(yneigh = y-1; yneigh <= y+1; yneigh++) {
		if((xneigh != x || yneigh != y) && current[xneigh][yneigh] >= TUMOUR) {
			neigh_cell_num = current[xneigh][yneigh];
			if(death_due_to_resource == 0) (*sumsums) -= occupied[CELL_NOW_SUMRATES][neigh_cell_num];
			set_proliferation_rate(u, est, tam, inhibitor, occupied, neigh_cell_num, prolif_rate_limit_ratio, prolif_hypothesis, current, qui_pro_rate_ratio, 0);
			if(death_due_to_resource == 0) (*sumsums) += occupied[CELL_NOW_SUMRATES][neigh_cell_num];
			}
		}
}

void cellreplacement(int **current, int *pop, int x, int y, int *event_count, int max_colour_scale, int dim, int replace_radius, long *idum, int **n_x,
 int **n_y, float **n_p, int *n_n, float **occupied, float *sumsums, double **u, double **est, double **tam, double **inhibitor,
 float *prolif_rate_limit_ratio, int prolif_hypothesis, int *occ, float qui_pro_rate_ratio)
{
	int i,j,ring;
	int site;
	int xnew, ynew;
	int diameter = replace_radius*2+1;
	int neigh[diameter*diameter]; // sequence of neighbouring cell indices, to be shuffled before use
	int predator = current[x][y];
	int prey;
	int subtype_predator = occupied[CELL_CLASS][predator];
	int subtype_prey;

	if(occupied[CELL_NOW_PROLIFRATE][predator] == 0) return; // quiescent cells don't replace others

	for(ring=0;ring<replace_radius;ring++) {
		for(i=0;i<n_n[ring];i++) neigh[i]=i; // sequence of neighbouring cell indices, to be shuffled before use
		shuffle(neigh,n_n[ring],idum);
		for(i=0;i<n_n[ring];i++) {
			site=neigh[i]; // cycle through all neighbourhood sites
			xnew=x+n_x[ring][site];
			ynew=y+n_y[ring][site];
			prey = current[xnew][ynew];
			subtype_prey = occupied[CELL_CLASS][prey];
			
			// if neighbour belongs to other subpopulation then neighbour is replaced:
			if(subtype_prey != subtype_predator && ran0(idum) < n_p[ring][site]) {
				
				(*sumsums) -= occupied[CELL_NOW_SUMRATES][prey];
				
				// copy all fixed properties:
				for(j = 0; j < NUM_FIXED_CELL_PROPERTIES; j++) occupied[j][prey]=occupied[j][predator];
				// set prolif rate and prolif rate limit:
				set_proliferation_rate(u, est, tam, inhibitor, occupied, prey, prolif_rate_limit_ratio, prolif_hypothesis, current, qui_pro_rate_ratio, 0);
				
				(*sumsums) += occupied[CELL_NOW_SUMRATES][prey];
				
				(*event_count)++;
				pop[subtype_prey]--;
				pop[subtype_predator]++;

				return;
				}
			}
		}
}

void cellproliferation(int **current, double **u, float *prolif_rate_limit_ratio, int *pop, int x, int y, int *event_count, int max_colour_scale, int dim,
	int prolif_radius, long *idum, int **n_x, int **n_y, float **n_p, int *n_n, int subtype, int *stop_condition, float **occupied, int *occ,
	float membrane_radius, int cell_num, float generations, float *birthdays, int *num_bir, float *sumsums, float *er_values, float *resistance_values, float *o2coeff_values,
	float *min_resource_survival, float *growth_rate, float *death_rate, float *comp_rate, double **est, double **tam, double **inhibitor, int prolif_hypothesis,
	float mutation_rate, int *mutated, float qui_pro_rate_ratio, float stem_prob, int max_divisions, float pro_mutation_rate,
	float pro_mutation_mean, float pro_mutation_sd)
{
	int i,j,ring;
	int site;
	int xnew, ynew;
	int x_to_fill, y_to_fill;
	int xnext, ynext;
	int x_base, y_base;
	float xdiff, ydiff;
	int xneigh, yneigh;
	int diameter = prolif_radius*2+1;
	int neigh[diameter*diameter]; // sequence of neighbouring cell indices, to be shuffled before use
	float hypoten;
	int outside_membrane;
	int loops;
	int new_cell_num;
	int neigh_cell_num;
	int do_mutation;

	if(occupied[CELL_STEM][cell_num] == 0 && occupied[CELL_DIVISIONS][cell_num] == 0) return;
	// if cell is a non-stem cell and has reached its divisions limit then it cannot proliferate

	for(ring=0;ring<prolif_radius;ring++) {
		for(i=0;i<n_n[ring];i++) neigh[i]=i; // sequence of neighbouring cell indices, to be shuffled before use
		shuffle(neigh,n_n[ring],idum);

		for(i=0;i<n_n[ring];i++) {
			site=neigh[i]; // cycle through all neighbourhood sites
			
			xnew=x+n_x[ring][site]; // coordinates of candidate site
			ynew=y+n_y[ring][site];
			
			outside_membrane = 0;
			if(membrane_radius > 0 && sqrt((xnew-dim/2)*(xnew-dim/2) + (ynew-dim/2)*(ynew-dim/2)) >= membrane_radius) outside_membrane = 1;
			
			if((current[xnew][ynew] == DEAD || current[xnew][ynew] == EMPTY || current[xnew][ynew] == NORMAL) && ran0(idum)<n_p[ring][site] && outside_membrane==0) {
				// if site is not occupied by a tumour cell
				// n_p[][] = proportion of the site that is within the neighbourhood

				x_to_fill = xnew; // coordinates of site that will receive a cell through budging or proliferation
				y_to_fill = ynew;
				
				loops = 0;
				do {
					x_base = x_to_fill;
					y_base = y_to_fill;
					do {
						xdiff = x - x_base;
						ydiff = y - y_base;
						hypoten = sqrt(xdiff*xdiff + ydiff*ydiff); // distance between sites
						if(fabs(xdiff/hypoten) >= ROOT2-1) xnext = x_base + sign(xdiff);
						else xnext = x_base;
						if(fabs(ydiff/hypoten) >= ROOT2-1) ynext = y_base + sign(ydiff);
						else ynext = y_base;
						// coordinates of the cell that is to be budged into the site (x_to_fill,y_to_fill)
						// (xnext,ynext) is between (x,y) and (x_to_fill,y_to_fill)
						x_base = xnext;
						y_base = ynext;
						} while(current[xnext][ynext] == VESSEL);

					if(xnext == x && ynext == y) break; // no more budging to be done
					else { // budge the cell from (xnext,ynext) to (x_to_fill,y_to_fill)
						if(current[xnext][ynext] >= TUMOUR) {
							// record the budging of a cell from (xnext,ynext) to (x_to_fill,y_to_fill):
							occupied[CELL_XCOORD][current[xnext][ynext]] = x_to_fill;
							occupied[CELL_YCOORD][current[xnext][ynext]] = y_to_fill;
							current[x_to_fill][y_to_fill] = current[xnext][ynext];
							current[xnext][ynext] = NORMAL;
							}
						x_to_fill = xnext;
						y_to_fill = ynext;
						}
					loops++;
					if(loops>999) {
						printf("Problem with proliferation\n");
						exit(0);
						}
					} while(1); // end of do loop
				
				// place a new daughter cell:
				new_cell_num = *occ;
				occupied[CELL_XCOORD][new_cell_num]=x_to_fill; 
				occupied[CELL_YCOORD][new_cell_num]=y_to_fill;
				
				do_mutation = 0;

				if(ran0(idum) < mutation_rate || do_mutation == 1) { // mutation of cell type
					mutate_cell(occupied, new_cell_num, er_values, resistance_values, o2coeff_values, min_resource_survival, prolif_rate_limit_ratio, growth_rate, death_rate, comp_rate, subtype);
					subtype = occupied[CELL_CLASS][new_cell_num];
					}
				// copy all fixed properties:
				else for(j = 0; j < NUM_FIXED_CELL_PROPERTIES; j++) occupied[j][new_cell_num] = occupied[j][current[x][y]];

				if(ran0(idum) < pro_mutation_rate) { // mutation of proliferation rate
					occupied[CELL_MAX_PROLIFRATE][new_cell_num] += pro_mutation_sd * gasdev(idum) + pro_mutation_mean;
					if(occupied[CELL_MAX_PROLIFRATE][new_cell_num] < 0) occupied[CELL_MAX_PROLIFRATE][new_cell_num] = 0;
					}
				
				if(occupied[CELL_STEM][current[x][y]] == 1) { // if the parent is a stem cell
					if(ran0(idum) < stem_prob) {
						occupied[CELL_STEM][new_cell_num] = 1; // daughter is a stem cell
						occupied[CELL_DIVISIONS][new_cell_num] = 99; // this is arbitrary
						}
					else {
						occupied[CELL_STEM][new_cell_num] = 0; // daughter is a non-stem cell
						occupied[CELL_DIVISIONS][new_cell_num] = max_divisions;
						}
					}
				else { // if the parent is a non-stem cell
					occupied[CELL_STEM][new_cell_num] = 0; // daughter is a non-stem cell
					occupied[CELL_DIVISIONS][current[x][y]]--; // parent has one less division remaining
					occupied[CELL_DIVISIONS][new_cell_num] = occupied[CELL_DIVISIONS][current[x][y]]; // daughter inherits parent's new division number
					}

				set_proliferation_rate(u, est, tam, inhibitor, occupied, new_cell_num, prolif_rate_limit_ratio, prolif_hypothesis, current, qui_pro_rate_ratio, 0);

				(*sumsums) += occupied[CELL_NOW_SUMRATES][new_cell_num];
				current[x_to_fill][y_to_fill] = new_cell_num;
				(*occ)++;
				(*event_count)++;
				pop[subtype]++;
				birthdays[*num_bir] = generations;
				(*num_bir)++;

				// set proliferation rate for 8 nearest neighbours:
				for(xneigh = xnew-1; xneigh <= xnew+1; xneigh++) for(yneigh = ynew-1; yneigh <= ynew+1; yneigh++) {
					if((xneigh != xnew || yneigh != ynew) && current[xneigh][yneigh] >= TUMOUR) {
						neigh_cell_num = current[xneigh][yneigh];
						(*sumsums) -= occupied[CELL_NOW_SUMRATES][neigh_cell_num];
						set_proliferation_rate(u, est, tam, inhibitor, occupied, neigh_cell_num, prolif_rate_limit_ratio, prolif_hypothesis, current, qui_pro_rate_ratio, 0);
						(*sumsums) += occupied[CELL_NOW_SUMRATES][neigh_cell_num];
						}
					}

				if(xnew == 1 || xnew == dim-2 || ynew == 1 || ynew == dim-2) *stop_condition=1; // tumour has reached edge of grid

				return;
				} // end of if block
			} // end of loop through neighbourhood sites
		} // end of loop through neighbourhoods (rings of increasing size)
}

void add_mutated_cells(float initial_propn, int num_occupied, long *idum, float **occupied, int *pop, float *sumsums, float *er_values, float *resistance_values, float *o2coeff_values, 
	float *min_resource_survival, float *prolif_rate_limit_ratio, float *growth_rate, float *death_rate, float *comp_rate)
{
	int i, subtype;
	int num_to_mutate = (int)((1 - initial_propn) * num_occupied + 0.5);
	int samples[num_to_mutate];

	SampleWithoutReplacement(num_occupied, num_to_mutate, samples, idum);
	for(i = 0; i < num_to_mutate; i++) {
		subtype = occupied[CELL_CLASS][samples[i]];
		(*sumsums) -= occupied[CELL_NOW_SUMRATES][samples[i]];
		pop[subtype]--;
		mutate_cell(occupied, samples[i], er_values, resistance_values, o2coeff_values, min_resource_survival, prolif_rate_limit_ratio, growth_rate, death_rate, comp_rate, subtype);
		subtype = occupied[CELL_CLASS][samples[i]];
		(*sumsums) += occupied[CELL_NOW_SUMRATES][samples[i]];
		pop[subtype]++;
	}
}

void mutate_cell(float **occupied, int cell_num, float *er_values, float *resistance_values, float *o2coeff_values, float *min_resource_survival, float *prolif_rate_limit_ratio,
	float *growth_rate, float *death_rate, float *comp_rate, int subtype)
{
	int class = abs(1 - subtype); // opposite of parent's subtype

	printf("\nMUTATION!!\n");

	occupied[CELL_CLASS][cell_num] = class;
	occupied[CELL_DEATHRATE][cell_num]=death_rate[class];
	occupied[CELL_MAX_PROLIFRATE][cell_num]=growth_rate[class];
	occupied[CELL_SURVLIMIT][cell_num]=min_resource_survival[class];
	occupied[CELL_COMPRATE][cell_num]=comp_rate[class];
	occupied[CELL_ER][cell_num] = er_values[class];
	occupied[CELL_RESISTANCE][cell_num] = resistance_values[class];
	occupied[CELL_O2COEFF][cell_num] = o2coeff_values[class];
}

void set_proliferation_rate(double **u, double **est, double **tam, double **inhibitor, float **occupied, int cell_num, float *prolif_rate_limit_ratio,
	int prolif_hypothesis, int **current, float qui_pro_rate_ratio, int print_output)
{
	int x,y;
	int subtype;
	float e,p,g,t,h,r,Pmax,G;
	int xnew,ynew;
	float pro_rate, excess_death_rate;
	int qui;

	x = occupied[CELL_XCOORD][cell_num];
	y = occupied[CELL_YCOORD][cell_num];

	subtype = occupied[CELL_CLASS][cell_num];

	e = est[x][y];
	p = occupied[CELL_ER][cell_num];
	g = u[x][y];
	t = tam[x][y];
	h = inhibitor[x][y];
	r = occupied[CELL_RESISTANCE][cell_num];
	Pmax = occupied[CELL_MAX_PROLIFRATE][cell_num];
	G = occupied[CELL_O2COEFF][cell_num];

	excess_death_rate = 0;
	
	// Hypothesis 0:
	if(prolif_hypothesis==0) pro_rate = Pmax;

	// Hypothesis 1:
	if(prolif_hypothesis==1) pro_rate = Pmax * G*g / (G*g+1) * 2000*p*e / (2000*p*e + t + 2000);

	// Hypothesis 2:
	if(prolif_hypothesis==2) pro_rate = Pmax * G*g / (G*g+1) * 2000*p*e / (2000*e + t + 2000);

	// Hypothesis 3:
	if(prolif_hypothesis==3) {
		pro_rate = Pmax * G*g / (G*g+1) * maxi(1 - 4.99E-8 * h / r, 0);
		if(print_output > 0) {
			printf("Pmax = %f, G = %f, g = %f, h = %f, r = %f\n", Pmax, G, g, h, r);
			printf("pro_rate = %f * %f\n", Pmax * G*g / (G*g+1), maxi(1 - 4.99E-8 * h / r, 0));
		}
		excess_death_rate = -mini((0.361 - 1.82E-8 * h) / r, 0);
	}

	occupied[CELL_NOW_PROLIFLIMIT][cell_num] = prolif_rate_limit_ratio[subtype] * Pmax * G*g / (G*g+1);

	if(pro_rate >= qui_pro_rate_ratio * Pmax) qui = 0;
	else {
		qui = 1; // examine 8 nearest neighbours: if cell is surrounded on all sides then it is quiescent
		for(xnew = x-1; xnew <= x+1; xnew++) for(ynew = y-1; ynew <= y+1; ynew++) {
			if((xnew != x || ynew != y) && (current[xnew][ynew] == NORMAL || current[xnew][ynew] == EMPTY)) {
				qui = 0;
				break;
				}
			}
		}

	if(qui == 1) printf("Found a quiescent cell\n");

	if(print_output > 0) printf("qui = %d\n", qui);

	if(qui == 1) occupied[CELL_NOW_PROLIFRATE][cell_num] = 0;
	else occupied[CELL_NOW_PROLIFRATE][cell_num] = pro_rate;
	occupied[CELL_NOW_EXCESSDEATHRATE][cell_num] = excess_death_rate;

	occupied[CELL_NOW_SUMRATES][cell_num] = occupied[CELL_NOW_PROLIFRATE][cell_num] + occupied[CELL_DEATHRATE][cell_num] + occupied[CELL_NOW_EXCESSDEATHRATE][cell_num] + occupied[CELL_COMPRATE][cell_num];
}

void sor(double **a, double **b, double **c, double **d, double **e1, double **e2, double **f, double **u, int dim, double rjac, int **current, long *idum,
 int *stop_condition, bool include_vessels, float serum, int num_occupied, float **occupied)
	// Successive overrelaxation solution with Chebyshev acceleration
	// a, b, c, d, e, and f are input as the coeﬃcients of the equation, each dimensioned to the grid size
	// u is input as the initial guess to the solution, usually zero
	// u returns with the ﬁnal value
	// rjac is input as the spectral radius of the Jacobi iteration, or an estimate of it
	{
	int n,x,y,i,j;
	double anorm,anormf=0.0,omega=1.0,resid;
	double ua,ub,uc,ud; // Double precision is a good idea for dim bigger than about 25
	int minx, maxx, miny, maxy;
	float e;

	if(include_vessels == false) { // find the limits of the tumour
		minx = dim-1;
		miny = dim-1;
		maxx = 0;
		maxy = 0;
		for(i=0; i<dim; i++) for(j=0; j<dim; j++) if(current[i][j] >= TUMOUR || current[i][j] == DEAD || current[i][j] == EMPTY) {
			if(i<minx) minx = i;
			if(i>maxx) maxx = i;
			if(j<miny) miny = j;
			if(j>maxy) maxy = j;
			}
		minx = maxi(minx-1, 0);
		maxx = mini(maxx+1, dim-1);
		miny = maxi(miny-1, 0);
		maxy = mini(maxy+1, dim-1);
		}
	else {
		minx = 0;
		maxx = dim-1;
		miny = 0;
		maxy = dim-1;
	}

	if(include_vessels == false) for(x = 0; x < dim; x++) for(y = 0; y < dim; y++) if(current[x][y] == NORMAL) u[x][y] = serum; // uniform resource level

	for (n=1;n<=MAXITS;n++) {

		anorm=0.0;

		for(x = minx; x <= maxx; x++) for(y = miny; y <= maxy; y++) {

			if(include_vessels == true && current[x][y] == VESSEL) u[x][y] = 0; // don't solve for vessels

			else if(include_vessels == false && current[x][y] == NORMAL) u[x][y] = serum; // uniform resource level

			else {

				if(x==maxx) ua = u[minx][y]; // at the edge
				else ua = u[x+1][y]; // not at the edge

				if(x==minx) ub = u[maxx][y];
				else ub = u[x-1][y];

				if(y==maxy) uc = u[x][miny];
				else uc = u[x][y+1];

				if(y==miny) ud = u[x][maxy];
				else ud = u[x][y-1];

				e = e1[x][y] + e2[x][y];

				resid = a[x][y]*ua + b[x][y]*ub + c[x][y]*uc + d[x][y]*ud + e*u[x][y] + f[x][y]; // residual at x,y
				if(isnan(resid)) {
					printf("\nProblem in solving diffusion equations!\nresid = %f * %f\n+ %f * %f\n+ %f * %f\n+ %f * %f\n+ %f * %f\n+ %f\nx = %d, y = %d,\ni = %d\nn = %d\n",
							a[x][y],ua, b[x][y],ub, c[x][y],uc, d[x][y],ud, e,u[x][y], f[x][y], x, y, i, n);
					getchar();
					}
				anorm += fabs(resid); // sum of residuals, which is to be minimised
				if(n==1) anormf+=fabs(resid); // initial sum of residuals

				u[x][y] -= omega*resid/e; // update u[][]

				}

			}

		omega=(n == 1 ? 1.0/(1.0-0.5*rjac*rjac) : 1.0/(1.0-0.25*rjac*rjac*omega)); // update omega

		if (anorm < EPS*anormf) {
			for(i=minx;i<=maxx;i++) for(j=miny;j<=maxy;j++) u[i][j]=maxi(u[i][j],0); // because u[][] values can sometimes be very slightly negative
			printf("%d sor iterations\n",n);
			return;
			}
		}
	
	printf("\nresiduals sum = %f (target is %f); latest residual = %f\n",anorm,EPS*anormf,resid);
	*stop_condition=3;
}

void write_matrix(int **current, double **u, int dim, int loop_count, int maxloops, char *preamble_text, int max_colour_scale, FILE *gp, int display_type,
 float *min_resource_survival, float **occupied, double **tam, double **inhibitor, float gens_elapsed, float sum_gens_elapsed)
{
	int i,j;
	float maxu; // maximum value of u[][] at the end of the diffusion solver, should equal the serum level
	float minu; // minimum value of u[][] for cell survival
	float maxtam, maxinhibitor;
	float u_temp[dim][dim];
	char *buffer = malloc(20+1); // string buffer
	if(buffer == NULL) {printf("Memory problem (write_matrix)!\n"); exit(0);}
	int delay;

	delay = (int)(10.0*gens_elapsed);

	maxu=0;
	for(i=0;i<dim;i++) for(j=0;j<dim;j++) if(current[i][j] == NORMAL) if(u[i][j]>maxu) maxu=u[i][j];

	maxtam=0;
	for(i=0;i<dim;i++) for(j=0;j<dim;j++) if(current[i][j] == NORMAL) if(tam[i][j]>maxtam) maxtam=tam[i][j];

	maxinhibitor=0;
	for(i=0;i<dim;i++) for(j=0;j<dim;j++) if(current[i][j] == NORMAL) if(inhibitor[i][j]>maxinhibitor) maxinhibitor=inhibitor[i][j];
	
	minu=mini(min_resource_survival[0],min_resource_survival[1]);

	for(i=0;i<dim;i++) for(j=0;j<dim;j++) {
		if(current[i][j] == VESSEL) { // sets a reasonable u[][] value at vessels, for better plotting
			if(maxtam > 0) u_temp[i][j]=maxtam;
			if(maxinhibitor > 0) u_temp[i][j]=maxinhibitor;
			else u_temp[i][j]=maxu;
			}
		else if(current[i][j] >= TUMOUR && u[i][j]<minu) u_temp[i][j]=minu; // minu is the threshold for cell survival
		else u_temp[i][j] = u[i][j];
		}

	fprintf(gp, "set terminal png\n");
	fprintf(gp, "set output '");
	fprintf(gp, "animation_frames_type/ani%s%d_%d_%d.png'\n",zeroes(loop_count,maxloops),loop_count,(int)sum_gens_elapsed,delay);
	fprintf(gp, preamble_text, 0);
	for(i=0;i<dim;i++) {
		for(j=0;j<dim;j++) {
			if(current[i][j] == VESSEL || current[i][j] == DEAD || current[i][j] == EMPTY || current[i][j] == NORMAL) {
				if(tam[i][j] > 0) sprintf(buffer, "%.4f ", 3.001*max_colour_scale+0.999*max_colour_scale*tam[i][j]/maxtam); // if is not a tumour cell or membrane then write the resource level
				if(inhibitor[i][j] > 0) sprintf(buffer, "%.4f ", 3.001*max_colour_scale+0.999*max_colour_scale*inhibitor[i][j]/maxinhibitor); // if is not a tumour cell or membrane then write the resource level
				else sprintf(buffer, "%.4f ", 2.001*max_colour_scale+0.999*max_colour_scale*u_temp[i][j]/maxu); // if is not a tumour cell or membrane then write the resource level
				}
			else if(current[i][j] == MEMBRANE) sprintf(buffer, "%d ", max_colour_scale*3+1); // if it is a membrane cell
			else { // if it is a tumour cell
				 // write its resource level:
				sprintf(buffer, "%.4f ", 1 + occupied[CELL_CLASS][current[i][j]]*max_colour_scale + 0.999*(max_colour_scale-1)*(1-(u_temp[i][j]-minu)/(maxu-minu))*(1-(u_temp[i][j]-minu)/(maxu-minu)));
				// NB: the proportion is squared to produce a more gradual shading transition
				}
			fprintf(gp, buffer, 0);
			}
		fprintf(gp, "\n");
		}
	fprintf(gp, "e\n");
	fprintf(gp, "e\n");
}

void write_matrix_occupied(int **current, double **u, int dim, int loop_count, int maxloops, char *preamble_text2, int max_colour_scale, FILE *gp, int display_type,
 int num_occupied, float **occupied, int death_instead_of_birth, int type_to_plot)
{
	int i,j;
	int binary_occupied[dim][dim];
	char *buffer = malloc(20+1); // string buffer
	if(buffer == NULL) {printf("Memory problem (write_matrix)!\n"); exit(0);}
	int x,y,type;
	float scale_birth = 2, scale_death = 0.2;
	float val_to_write;
	char *suffix = malloc(50+1);

	for(i=0;i<dim;i++) for(j=0;j<dim;j++) binary_occupied[i][j] = 0;

	if(!death_instead_of_birth) sprintf(suffix, "_birth");
	else sprintf(suffix, "_death");
	if(type_to_plot < 2) sprintf(suffix + strlen(suffix), "_type%d", type_to_plot);

	fprintf(gp, "set terminal png\n");
	fprintf(gp, "set output '");
	fprintf(gp, "animation_frames%s/occupied%s%d.png'\n", suffix, zeroes(loop_count, maxloops), loop_count);
	fprintf(gp, preamble_text2, 0);

	for(i=0;i<num_occupied;i++) {
		x = occupied[CELL_XCOORD][i];
		y = occupied[CELL_YCOORD][i];
		type = occupied[CELL_CLASS][i];
		if(type == type_to_plot || type_to_plot == 2) binary_occupied[x][y] = 1;
		}

	for(i=0;i<dim;i++) {
		for(j=0;j<dim;j++) {
			if(!death_instead_of_birth) {
				if(binary_occupied[i][j] == 0) val_to_write = 1.5;
				else if(occupied[CELL_NOW_PROLIFRATE][current[i][j]] == 0) val_to_write = 2.5;
				else val_to_write = 3.1 + 8.9 * occupied[CELL_NOW_PROLIFRATE][current[i][j]] / scale_birth;
				sprintf(buffer, "%f ", val_to_write); // write its occupancy state
				fprintf(gp, buffer, 0);
				}
			else {
				if(binary_occupied[i][j] == 0) val_to_write = 1.5;
				else if(occupied[CELL_NOW_EXCESSDEATHRATE][current[i][j]] == 0) val_to_write = 2.5;
				else val_to_write = 3.1 + 8.9 * occupied[CELL_NOW_EXCESSDEATHRATE][current[i][j]] / scale_death;
				sprintf(buffer, "%f ", val_to_write); // write its occupancy state
				fprintf(gp, buffer, 0);
				}
			if(val_to_write > 12.0) val_to_write = 12.0;
			}
		fprintf(gp, "\n");
		}
	fprintf(gp, "e\n");
	fprintf(gp, "e\n");
}

void write_matrix_maxprolif(int **current, double **u, int dim, int loop_count, int maxloops, char *preamble_text2, int max_colour_scale, FILE *gp, int display_type,
 int num_occupied, float **occupied)
{
	int i,j;
	int binary_occupied[dim][dim];
	char *buffer = malloc(20+1); // string buffer
	if(buffer == NULL) {printf("Memory problem (write_matrix)!\n"); exit(0);}
	int x,y;
	float val_to_write;
	float max = 2; // value at which the scale saturates

	for(i=0;i<dim;i++) for(j=0;j<dim;j++) binary_occupied[i][j] = 0;

	fprintf(gp, "set terminal png\n");
	fprintf(gp, "set output '");
	fprintf(gp, "animation_frames_maxbirth/maxprolif%s%d.png'\n",zeroes(loop_count,maxloops),loop_count);
	fprintf(gp, preamble_text2, 0);

	for(i=0;i<num_occupied;i++) {
		x = occupied[CELL_XCOORD][i];
		y = occupied[CELL_YCOORD][i];
		binary_occupied[x][y] = 1;
		}

	for(i=0;i<dim;i++) {
		for(j=0;j<dim;j++) {
			if(binary_occupied[i][j] == 0) val_to_write = 1.5;
			else if(occupied[CELL_MAX_PROLIFRATE][current[i][j]] == 0) val_to_write = 2.5;
			else val_to_write = 3.1 + 8.9 * occupied[CELL_MAX_PROLIFRATE][current[i][j]] / max;
			if(val_to_write > 12.0) val_to_write = 12.0;
			sprintf(buffer, "%f ", val_to_write); // write its occupancy state
			fprintf(gp, buffer, 0);
			}
		fprintf(gp, "\n");
		}
	fprintf(gp, "e\n");
	fprintf(gp, "e\n");
}

void write_matrix_stem(int **current, double **u, int dim, int loop_count, int maxloops, char *preamble_text2, int max_colour_scale, FILE *gp, int display_type,
 int num_occupied, float **occupied)
{
	int i,j;
	char *buffer = malloc(20+1); // string buffer
	if(buffer == NULL) {printf("Memory problem (write_matrix)!\n"); exit(0);}
	float val_to_write;

	fprintf(gp, "set terminal png\n");
	fprintf(gp, "set output '");
	fprintf(gp, "animation_frames_stem/stem%s%d.png'\n",zeroes(loop_count,maxloops),loop_count);
	fprintf(gp, preamble_text2, 0);

	for(i=0;i<dim;i++) {
		for(j=0;j<dim;j++) {
			if(current[i][j] < TUMOUR) val_to_write = 0.5; // if not a tumour cell
			else if(occupied[CELL_STEM][current[i][j]] == 1) val_to_write = 1.5;
			else val_to_write = occupied[CELL_DIVISIONS][current[i][j]] + 2.5;
			sprintf(buffer, "%f ", val_to_write); // write its state
			fprintf(gp, buffer, 0);
			}
		fprintf(gp, "\n");
		}
	fprintf(gp, "e\n");
	fprintf(gp, "e\n");
}

void write_matrix_to_file(int **current, double **u, int dim, int loop_count, int max_colour_scale, FILE *matrix, float *min_resource_survival, int maxloops)
{
	int i,j;
	float maxu; // maximum value of u[][] at the end of the diffusion solver, should equal the serum level
	float minu; // minimum value of u[][] at the end of the diffusion solver
	char *buffer = malloc(50+1); // string buffer
	if(buffer == NULL) {printf("Memory problem (write_matrix_to_file)!\n"); exit(0);}
	float u_temp[dim][dim];

	sprintf(buffer, "matrices/matrix_%s%d.dat", zeroes(loop_count,maxloops), loop_count);
	matrix = fopen(buffer, "w+");

	maxu=0;
	for(i=0;i<dim;i++) for(j=0;j<dim;j++) if(current[i][j] == NORMAL) if(u[i][j]>maxu) maxu=u[i][j];
	
	minu=mini(min_resource_survival[0],min_resource_survival[1]);

	for(i=0;i<dim;i++) for(j=0;j<dim;j++) {
		if(current[i][j] == VESSEL) u_temp[i][j]=maxu; // sets a reasonable u[][] value at vessels, for better plotting
		else if(current[i][j] >= TUMOUR && u[i][j]<minu) u_temp[i][j]=minu; // minu is the threshold for cell survival
		else u_temp[i][j] = u[i][j];
		}

	for(i=0;i<dim;i++) {
		for(j=0;j<dim;j++) {
			 // if is a normal cell or vessel then write the diffusion level:
			if(current[i][j] == VESSEL || current[i][j] == DEAD || current[i][j] == EMPTY || current[i][j] == NORMAL) sprintf(buffer, "%.4f ", u_temp[i][j]);
			else sprintf(buffer, "%d ", current[i][j]); // if it is a tumour or membrane cell then write its value
			fprintf(matrix, buffer, 0);
			}
		fprintf(matrix, "\n");
		}
}

void assign_coefficients(double **a, double **b, double **c, double **d, double **e1, double **e2, double **f, int **current, int dim, float permeability, float diffusion,
 float serum, float delta, float *cons_rate, double **u, float **occupied, float **vessel_size, float *estrogen_curve, float sum_gens_elapsed, bool estrogen)
{
	int i,j,k;
	int x,y;
	int next_to_vessel;
	float cons;
	float size;

	for(i=0;i<dim;i++) for(j=0;j<dim;j++) if(current[i][j] != VESSEL) { // doesn't assign values for blood vessels
		
		cons = consumption(current[i][j], cons_rate, occupied[CELL_NOW_PROLIFRATE][current[i][j]]); // look up the consumption rate
		next_to_vessel=0;

		// default coefficient values:
		a[i][j]=1;
		b[i][j]=1;
		c[i][j]=1;
		d[i][j]=1;

		y = j;
		if(i==0) x = dim - 1; // if the cell is at an edge
		else x = i - 1; // if not at an edge
		if(current[x][y] == VESSEL) { // if next to a vessel
			next_to_vessel = 1;
			size = vessel_size[x][y];
			b[i][j] = 0;
			}

		y = j;
		if(i==dim-1) x = 0;
		else x = i + 1;
		if(current[x][y] == VESSEL) {
			next_to_vessel = 1;
			size = vessel_size[x][y];
			a[i][j] = 0;
			}

		x = i;
		if(j==0) y = dim - 1;
		else y = j - 1;
		if(current[x][y] == VESSEL) {
			next_to_vessel = 1;
			size = vessel_size[x][y];
			d[i][j] = 0;
		}

		x = i;
		if(j==dim-1) y = 0;
		else y = j + 1;
		if(current[x][y] == VESSEL) {
			next_to_vessel = 1;
			size = vessel_size[x][y];
			c[i][j] = 0;
		}

		if(next_to_vessel) { // if cell is next to a vessel
			e1[i][j] = -3;
			e2[i][j] = (-delta * delta * 48 * cons / (48 * u[i][j] + 1) - delta * permeability) / diffusion;
			f[i][j] = permeability * delta * serum * size / diffusion;
			}
		else { // if not next to a vessel
			e1[i][j] = -4;
			e2[i][j] = - delta * delta / diffusion * 48 * cons / (48 * u[i][j] + 1);
			if(estrogen) e2[i][j] *= estrogen_curve[((int)sum_gens_elapsed)%30];
			f[i][j] = 0;
			}
	}
}

void define_neighbourhood(int **n_x, int **n_y, float **n_p, int *n_n, int max_radius)
{
	int i,j,k,ring;
	float d;
	float R;
	float r=0.5; // cell radius
	int diameter = 2*max_radius+1;
	int x,y;
	float p;
	int **x_temp,**y_temp;
	float **p_temp;
	int index;
	float part;

	x_temp = malloc(max_radius * sizeof *x_temp);
	if(x_temp == NULL) {printf("Memory problem!\n"); exit(0);}
	for (i=0; i<max_radius; i++) {
		x_temp[i] = malloc(diameter*diameter * sizeof *x_temp[i]);
		if(x_temp[i] == NULL) {printf("Memory problem!\n"); exit(0);}
		}

	y_temp = malloc(max_radius * sizeof *y_temp);
	if(y_temp == NULL) {printf("Memory problem!\n"); exit(0);}
	for (i=0; i<max_radius; i++) {
		y_temp[i] = malloc(diameter*diameter * sizeof *y_temp[i]);
		if(y_temp[i] == NULL) {printf("Memory problem!\n"); exit(0);}
		}

	p_temp = malloc(max_radius * sizeof *p_temp);
	if(p_temp == NULL) {printf("Memory problem!\n"); exit(0);}
	for (i=0; i<max_radius; i++) {
		p_temp[i] = malloc(diameter*diameter * sizeof *p_temp[i]);
		if(p_temp[i] == NULL) {printf("Memory problem!\n"); exit(0);}
		}

	for(ring = 0; ring<max_radius; ring++) {
		R = ring+1+0.5; // +1 because smallest neighbourhood has radius 1; +0.5 because the central site has radius 0.5
		for(i=0;i<diameter;i++) for(j=0;j<diameter;j++) {
			x = i - (diameter-1)/2;
			y = j - (diameter-1)/2;
			d = sqrt(x*x + y*y); // distance to the site
			
			if(x==0 && y==0) p = 0; // don't count the site at the centre of the neighbourhood
			else if(d>R+r) p = 0; // site is wholly outside neighbourhood
			else if(R>=d+r) p = 1; // site is wholly inside neighbourhood
			else p = r*r*acos((d*d + r*r - R*R)/(2*d*r)) + R*R*acos((d*d + R*R -r*r)/(2*d*R)) - 1.0/2*sqrt((-d+r+R)*(d+r-R)*(d-r+R)*(d+r+R));
			// weight assigned to site is proportional to how much of site is within neighbourhood
			// see http://mathworld.wolfram.com/Circle-CircleIntersection.html

			p_temp[ring][i*diameter+j] = p;

			// remove any sites wholly included in smaller neighbourhoods
			for(k=0;k<ring;k++) if(p_temp[ring-k-1][i*diameter+j] == 1) p_temp[ring][i*diameter+j] = 0;
			
			x_temp[ring][i*diameter+j] = x;
			y_temp[ring][i*diameter+j] = y;
			}
		}

	// store only the attributes of sites with n_p > 0:
	for(ring = 0; ring<max_radius; ring++) {
		index=0;
		for(i=0;i<diameter;i++) for(j=0;j<diameter;j++) {
			if(p_temp[ring][i*diameter+j] > 1E-5) {
				n_x[ring][index] = x_temp[ring][i*diameter+j];
				n_y[ring][index] = y_temp[ring][i*diameter+j];
				n_p[ring][index] = p_temp[ring][i*diameter+j];
				index++;
				}
			}
		n_n[ring] = index;
		}
}

float consumption(int state, float *cons_rate, float prolif_rate)
{
	float ans;

	if(state == NORMAL) ans = cons_rate[CONS_NORMAL];
	else if(state == DEAD || state == EMPTY) ans = cons_rate[CONS_DEAD];
	else if(state == MEMBRANE) ans = cons_rate[CONS_MEMBRANE];
	else if(state == VESSEL) ans = cons_rate[CONS_VESSEL];
	else if(state >= TUMOUR) {
		if(prolif_rate > 0) ans = cons_rate[CONS_PROLIF];
		else ans = cons_rate[CONS_QUI];
		}
	else {
		printf("\nInvalid state (%d, %f)\n", state, prolif_rate);
		exit(0);
		}

	return ans;
}

void create_vessels(float membrane_radius, int **current, int dim, float vessel_density, long *idum, int pattern, float **vessel_size)
{
	int x1,x2,y1,y2; // coordinates of corners of a prospective vessel site
	int xx1,xx2,yy1,yy2; // coordinates of neighbours of neighbours of a prospective vessel site
	int rnd,rnd2; // random integer variables
	int i,j,k;
	int num_steps=(int)rint(dim*sqrt(vessel_density));
	float step_size=(float)dim/num_steps;
	float centre = dim/2;
	float a, b;
	int x, y;
	float area = dim*dim - PI * membrane_radius*membrane_radius;
		
	// for Poisson Disk Sampling:
	int mindist = rint(sqrt(0.75/vessel_density));
	float cell_size = (float)mindist/ROOT2;
	int griddim = ceil(dim/cell_size);
	int samplePointsX[griddim*griddim];
	int samplePointsY[griddim*griddim];
	int samplePointsLength = 0;

	// random Poisson Disk Sampling vessel placement:
	if(pattern==0) {
		generate_poisson(dim, mindist, samplePointsX, samplePointsY, idum, &samplePointsLength);
		for(i=1;i<samplePointsLength;i++) {
			x = samplePointsX[i];
			y = samplePointsY[i];
			if(membrane_radius==0 || sqrt(x*x + y*y) < membrane_radius){
				current[x][y]=VESSEL;
				vessel_size[x][y] = truncated_log_gasdev(idum, 2.5); // assign vessel size
				}
			}
		}

	// random uniform vessel placement:
	else if(pattern==1) for(i=0;i<area * vessel_density;i++) {
		rnd = ran0(idum) * dim;
		rnd2 = ran0(idum) * dim;

		a = rnd - centre;
		b = rnd2 - centre;

		if(sqrt(a*a + b*b) > membrane_radius) {

			if(rnd==0) x1=dim-1; // if at the edge
			else x1=rnd-1;
			if(rnd<2) xx1=dim-(2-rnd); // if at the edge or one site away from the edge
			else xx1=rnd-2;

			if(rnd==dim-1) x2=0;
			else x2=rnd+1;
			if(rnd>dim-3) xx2=2-(dim-rnd);
			else xx2=rnd+2;

			if(rnd2==0) y1=dim-1;
			else y1=rnd2-1;
			if(rnd2<2) yy1=dim-(2-rnd2);
			else yy1=rnd2-2;

			if(rnd2==dim-1) y2=0;
			else y2=rnd2+1;
			if(rnd2>dim-3) yy2=2-(dim-rnd2);
			else yy2=rnd2+2;

			// place the vessel only if it is not too close to another vessel
			// too close means a diagonal neighbour is a vessel, or a neighbour of a neighbour (parallel to an axis) is a vessel
			if(current[x1][y1] != VESSEL && current[x1][y2]!= VESSEL && current[x2][y1]!= VESSEL && current[x2][y2]!= VESSEL
				&& current[rnd][yy1]!= VESSEL && current[rnd][yy2]!= VESSEL &&
				current[xx1][rnd2]!= VESSEL && current[xx2][rnd2]!= VESSEL) {
					current[rnd][rnd2] = VESSEL;
					vessel_size[rnd][rnd2] = truncated_log_gasdev(idum, 2.5); // assign vessel size
					printf("%f, ", vessel_size[rnd][rnd2]);
					}
			}
		}

	// regular vessel placement:
	else for(i=0;i<num_steps;i++) for(j=0;j<num_steps;j++) {
		a = i*step_size - centre;
		b = j*step_size - centre;
		//printf("%f, %f; %f, %f",a,b, sqrt(a*a + b*b), membrane_radius);
		//getchar();
		if(sqrt(a*a + b*b) > membrane_radius) {
			x = (int)rint(i*step_size);
			y = (int)rint(j*step_size);
			current[x][y] = VESSEL;
			vessel_size[x][y] = truncated_log_gasdev(idum, 2.5); // assign vessel size
			}
		}
}

void create_membrane(float membrane_radius, int dim, int **current)
{
	int i,j;
	float centre = dim/2;
	float a, b;

	for(i=0;i<dim;i++) for(j=0;j<dim;j++){
		a = i-centre;
		b = j-centre;
		if((int)sqrt(a*a + b*b) == membrane_radius || (int)sqrt(a*a + b*b) == membrane_radius-1) current[i][j] = MEMBRANE;
		}
}

float find_radius(int dim, int **current)
{
	int i, j;
	int direction;
	int x_add, y_add;
	float sum, radius;

	sum = 0;
	for(direction=0;direction<8;direction++) {
		for(i=0;i<dim/2;i++) {
			if(direction > 4) x_add = -1;
			else if(direction == 0 || direction == 4) x_add = 0;
			else x_add = 1;
			if(direction == 7 || direction < 2) y_add = -1;
			else if(direction == 6 || direction == 2) y_add = 0;
			else y_add = 1;
			if(current[dim/2 + x_add * i][dim/2 + y_add * i] >= TUMOUR || current[dim/2 + x_add * i][dim/2 + y_add * i] == DEAD) radius = i + 0.5;
		}
		if(x_add != 0 && y_add != 0) radius = radius * ROOT2;
		sum = sum + radius;
	}

	return(sum / 8);
}

void create_tumour(float initial_propn, int dim, int initial_diameter, int **current, int max_colour_scale, int *pop, float **occupied, int *num_occupied,
 long *idum, float membrane_radius, float *death_rate, float *growth_rate, float *min_resource_survival, float *prolif_rate_limit_ratio, float *comp_rate, float *er_values, float *resistance_values,
 float *o2coeff_values,float stem_prob, int max_divisions, int mutated, FILE *outermost_resistant_cell_location)
{
	int i,j;
	int initial_sequence[initial_diameter*initial_diameter]; // enables intial subpopulations to be randomly mixed in space
	int start_point_x, start_point_y; // x and y coordinates of the top corner of the initial population
	int class;
	float centre_spot;
	int initial_number_of_cells;
	int x_coord[initial_diameter*initial_diameter], y_coord[initial_diameter*initial_diameter];
	float distance, max_distance;

	if(membrane_radius == 0) {
		start_point_x = (dim-initial_diameter)/2; // centre of grid
		start_point_y = (dim-initial_diameter)/2; // centre of grid
		}
	else {
		start_point_x = (dim-initial_diameter)/2; // centre of grid
		start_point_y = (int)(dim/2 - membrane_radius + 1); // on inner edge of membrane
		}

	initial_number_of_cells = 0;
	centre_spot = (float)(initial_diameter-1)/2;
	for(i=0;i<initial_diameter;i++) for(j=0;j<initial_diameter;j++) {
		distance = sqrt((i - centre_spot)*(i - centre_spot) + (j - centre_spot)*(j - centre_spot));
		if(distance <= (float)initial_diameter/2) {
			x_coord[initial_number_of_cells] = start_point_x + i;
			y_coord[initial_number_of_cells] = start_point_y + j;
			initial_number_of_cells++;
		}
	}

	// create shuffled vector containing subpopulation types of the initial population:
	for(i=0;i<initial_number_of_cells;i++) {
		if(mutated == 1 || i <= (int)(initial_propn * initial_number_of_cells + 0.5) - 1) initial_sequence[i]=0;
		else initial_sequence[i]=1;
		}
	shuffle(initial_sequence, initial_number_of_cells, idum);

	// place cells in a disc:
	max_distance = -1;
	for(i=0;i<initial_number_of_cells;i++) {
		
		class = initial_sequence[i];
		
		current[x_coord[i]][y_coord[i]] = *num_occupied;
		pop[class]++;

		distance = sqrt((x_coord[i] - start_point_x - centre_spot)*(x_coord[i] - start_point_x - centre_spot) + (y_coord[i] - start_point_y - centre_spot)*(y_coord[i] - start_point_y - centre_spot));
		
		if(class == 1 && distance > max_distance) max_distance = distance;
		
		occupied[CELL_CLASS][*num_occupied]=class;
		occupied[CELL_XCOORD][*num_occupied]=x_coord[i];
		occupied[CELL_YCOORD][*num_occupied]=y_coord[i];
		occupied[CELL_DEATHRATE][*num_occupied]=death_rate[class];
		occupied[CELL_MAX_PROLIFRATE][*num_occupied]=growth_rate[class];
		occupied[CELL_SURVLIMIT][*num_occupied]=min_resource_survival[class];
		occupied[CELL_COMPRATE][*num_occupied]=comp_rate[class];
		occupied[CELL_ER][*num_occupied]=er_values[class];
		occupied[CELL_RESISTANCE][*num_occupied]=resistance_values[class];
		occupied[CELL_O2COEFF][*num_occupied]=o2coeff_values[class];
		occupied[CELL_STEM][*num_occupied]=1; // all initial cells are stem cells
		occupied[CELL_DIVISIONS][*num_occupied]=99; // this is arbitray, as stem cells don't have a max number of divisions
		(*num_occupied)++;
		}

	fprintf(outermost_resistant_cell_location, "%f\n", max_distance);
}

char *preamble(int dim, int max_colour_scale)
{
	// returns a string containing gnuplot code needed to set up each plot

	char *text = malloc(999+1);
	char limit[15];
	char max1[15];
	char max2[15];
	char max2a[15];
	char max3[15];
	char max4[15];
	char ratio[15];
	
	if(text == NULL) {printf("Memory problem!\n"); exit(0);}
	
	sprintf(limit, "%d", dim-1);
	
	sprintf(max1, "%d", max_colour_scale); // tumour cells of type 0 have colour values between 1 and max_colour_scale
	sprintf(max2, "%d", 2*max_colour_scale); // tumour cells of type 1 have colour values between max_colour_scale+1 and 2*max_colour_scale
	sprintf(max2a, "%f", 2.5*max_colour_scale); 
	sprintf(max3, "%d", 3*max_colour_scale); // limit for normal cells
	sprintf(max4, "%d", 4*max_colour_scale); // limit for membrane cells
	
	strcpy(text, "");
	
	strcpy(text, concat(text,"set size square\n")); // sets the canvas to be square

	strcpy(text, concat(text,"unset key\n")); // disables a key (or legend) describing plots on a plot
	strcpy(text, concat(text,"set tic scale 0\n")); // control of the major (labelled) tics on all axes at once
	
	// set the range of values to be coloured using the current palette:
	strcpy(text, concat(text,"set cbrange [1:"));
	strcpy(text, concat(text,max4));
	strcpy(text, concat(text,"]\n"));
	
	// define the divisions in the palette
	strcpy(text, concat(text,"set palette model RGB defined ("));

	strcpy(text, concat(text,"1 \"skyblue\", ")); // tumour cells of type 0 are blue
	strcpy(text, concat(text,max1));
	strcpy(text, concat(text," \"dark-blue\", "));

	strcpy(text, concat(text,max1));
	strcpy(text, concat(text," \"pink\", ")); // tumour cells of type 1 are red
	strcpy(text, concat(text,max2));
	strcpy(text, concat(text," \"dark-red\", "));

	strcpy(text, concat(text,max2));
	strcpy(text, concat(text," \"dark-grey\", ")); // normal cells and vessels are green
	strcpy(text, concat(text,max2a));
	strcpy(text, concat(text," \"dark-green\", "));
	strcpy(text, concat(text,max3));
	strcpy(text, concat(text," \"green\", "));

	strcpy(text, concat(text,max3));
	strcpy(text, concat(text," \"dark-grey\", ")); // membrane cells are yellow
	strcpy(text, concat(text,max4));
	strcpy(text, concat(text," \"orange\""));

	strcpy(text, concat(text,")\n"));
	
	strcpy(text, concat(text,"unset cbtics\n")); // removes major (labelled) tics on the color box axis
	
	// x range of the plot:
	strcpy(text, concat(text,"set xrange [-0.5:"));
	strcpy(text, concat(text,limit));
	strcpy(text, concat(text,".5]\n"));
	
	// y range of the plot:
	strcpy(text, concat(text,"set yrange [-0.5:"));
	strcpy(text, concat(text,limit));
	strcpy(text, concat(text,".5]\n"));
	
	strcpy(text, concat(text,"set view map\n")); // converts surface to a two-dimensional 'map' style view
	strcpy(text, concat(text,"splot '-' matrix with image\n")); // makes the plot (splot plots 3-d surfaces and data)

	return text;
}

char *preamble2(int dim)
{
	// returns a string containing gnuplot code needed to set up each plot

	char *text = malloc(999+1);
	char limit[15];
	
	if(text == NULL) {printf("Memory problem!\n"); exit(0);}

	sprintf(limit, "%d", dim-1);
	
	strcpy(text, "");
	
	strcpy(text, concat(text,"set size square\n")); // sets the canvas to be square

	strcpy(text, concat(text,"unset key\n")); // disables a key (or legend) describing plots on a plot
	strcpy(text, concat(text,"set tic scale 0\n")); // control of the major (labelled) tics on all axes at once
	
	// set the range of values to be coloured using the current palette:
	strcpy(text, concat(text,"set cbrange [1:12]\n"));
	
	// define the divisions in the palette
	strcpy(text, concat(text,"set palette model RGB defined ("));

	strcpy(text, concat(text,"1 \"black\", 2 \"black\", "));

	strcpy(text, concat(text,"2 \"dark-grey\", 3 \"dark-grey\", "));

	strcpy(text, concat(text,"3 \"purple\", 4.5 \"red\", "));

	strcpy(text, concat(text,"4.5 \"red\", 9 \"yellow\", "));

	strcpy(text, concat(text,"9 \"yellow\", 12 \"white\""));

	strcpy(text, concat(text,")\n"));
	
	strcpy(text, concat(text,"unset cbtics\n")); // removes major (labelled) tics on the color box axis
	
	// x range of the plot:
	strcpy(text, concat(text,"set xrange [-0.5:"));
	strcpy(text, concat(text,limit));
	strcpy(text, concat(text,".5]\n"));
	
	// y range of the plot:
	strcpy(text, concat(text,"set yrange [-0.5:"));
	strcpy(text, concat(text,limit));
	strcpy(text, concat(text,".5]\n"));
	
	strcpy(text, concat(text,"set view map\n")); // converts surface to a two-dimensional 'map' style view
	strcpy(text, concat(text,"splot '-' matrix with image\n")); // makes the plot (splot plots 3-d surfaces and data)

	return text;
}

char *preamble3(int dim, int max_divisions)
{
	// returns a string containing gnuplot code needed to set up each plot

	char *text = malloc(999+1);
	char limit[15];
	char max_div[15];
	
	if(text == NULL) {printf("Memory problem!\n"); exit(0);}

	sprintf(limit, "%d", dim-1);
	sprintf(max_div, "%f", max_divisions+2.5);
	
	strcpy(text, "");
	
	strcpy(text, concat(text,"set size square\n")); // sets the canvas to be square

	strcpy(text, concat(text,"unset key\n")); // disables a key (or legend) describing plots on a plot
	strcpy(text, concat(text,"set tic scale 0\n")); // control of the major (labelled) tics on all axes at once
	
	// set the range of values to be coloured using the current palette:
	strcpy(text, concat(text,"set cbrange [0:"));
	strcpy(text, concat(text,max_div));
	strcpy(text, concat(text,"]\n"));
	
	// define the divisions in the palette
	strcpy(text, concat(text,"set palette model RGB defined ("));

	strcpy(text, concat(text,"0 \"white\", 1 \"white\", "));

	strcpy(text, concat(text,"1 \"red\", 2 \"red\", "));

	strcpy(text, concat(text,"2 \"black\", "));
	strcpy(text, concat(text,max_div));
	strcpy(text, concat(text," \"dark-grey\""));

	strcpy(text, concat(text,")\n"));
	
	strcpy(text, concat(text,"unset cbtics\n")); // removes major (labelled) tics on the color box axis
	
	// x range of the plot:
	strcpy(text, concat(text,"set xrange [-0.5:"));
	strcpy(text, concat(text,limit));
	strcpy(text, concat(text,".5]\n"));
	
	// y range of the plot:
	strcpy(text, concat(text,"set yrange [-0.5:"));
	strcpy(text, concat(text,limit));
	strcpy(text, concat(text,".5]\n"));
	
	strcpy(text, concat(text,"set view map\n")); // converts surface to a two-dimensional 'map' style view
	strcpy(text, concat(text,"splot '-' matrix with image\n")); // makes the plot (splot plots 3-d surfaces and data)

	return text;
}

////////////////////////////////
// Start of generic functions //
////////////////////////////////

// Start of Poisson disc sampling functions

void generate_poisson(int dim, int mindist, int *samplePointsX, int *samplePointsY, long *idum, int *samplePointsLength)
{
  // fills the samplePoints array with the coordinates of the chosen points

  //Create the grid
  float cell_size = (float)mindist/ROOT2;
  int griddim = ceil(dim/cell_size);
  int **grid;
  int processList[griddim*griddim];
  int processListLength = 0;
  int rnd;
  int max_tries = 1000; // for faster performance, use 30
  int point[2];
  int newpoint[2];
  int tries_count;
  int success;
  int i, j;

  grid = malloc(griddim * sizeof *grid);
  if(grid == NULL) {printf("Memory problem!\n"); exit(0);}
  for (i=0; i<griddim; i++) {
	grid[i] = malloc(griddim * sizeof *grid[i]);
	if(grid[i] == NULL) {printf("Memory problem!\n"); exit(0);}
	}

  for(i=0;i<griddim;i++) for(j=0;j<griddim;j++) grid[i][j] = -1;

  //generate the first point randomly
  samplePointsX[0] = ran0(idum)*dim;
  samplePointsY[0] = ran0(idum)*dim;
  grid[image_to_grid(samplePointsX[0],cell_size)][image_to_grid(samplePointsY[0],cell_size)] = 0;
  processList[0] = 0;
  processListLength = 1;
  (*samplePointsLength) = 1;

  //generate other points from points in queue.
  while (processListLength > 0) {
	
	rnd = ran0(idum)*processListLength;
	point[0] = samplePointsX[processList[rnd]];
	point[1] = samplePointsY[processList[rnd]];
	success = 0;

	for (tries_count = 0; tries_count < max_tries; tries_count++) {

	  generateRandomPointAround(point, newpoint, mindist, idum);
	  
	  //check that the point is in the image region
	  //and no points exists in the point's neighbourhood
	  if (inRectangle(newpoint, dim) == 1 && inNeighbourhood(grid, newpoint, mindist, cell_size, griddim, samplePointsX, samplePointsY) == 0) {
		
		processList[processListLength] = *samplePointsLength;
		samplePointsX[*samplePointsLength] = newpoint[0];
		samplePointsY[*samplePointsLength] = newpoint[1];
		grid[image_to_grid(newpoint[0],cell_size)][image_to_grid(newpoint[1],cell_size)] = *samplePointsLength;

		processListLength++;
		(*samplePointsLength)++;

		success = 1;
		break;
		}
	  }

	if(success == 0) {
	  processList[rnd] = processList[processListLength-1];
	  processListLength--;
	  }
  }
}

void generateRandomPointAround(int *point, int *newpoint, int mindist, long *idum)
{ //non-uniform, favours points closer to the inner ring, leads to denser packings
  float r1 = ran0(idum);
  float r2 = ran0(idum);
  float radius, angle;

  radius = mindist * (r1 + 1); //random radius between mindist and 2 * mindist
  angle = 2 * PI * r2; //random angle

  newpoint[0] = point[0] + radius * cos(angle);
  newpoint[1] = point[1] + radius * sin(angle);
}

int inRectangle(int *newpoint, int dim)
{
  if(newpoint[0] < dim && newpoint[0] >= 0 && newpoint[1] < dim && newpoint[1] >= 0) return 1;
  else return 0;
}

int inNeighbourhood(int **grid, int *newpoint, int mindist, float cell_size, int griddim, int *samplePointsX, int *samplePointsY)
{
  int i, j;
  int gridPointX = newpoint[0] / cell_size; //get the neighbourhood of the point in the grid
  int gridPointY = newpoint[1] / cell_size;
  int minx, maxx, miny, maxy;
  float dist;
  int otherx, othery;

  if(gridPointX < 2) minx = 0;
  else minx = gridPointX-2;
  if(gridPointX > griddim-2) maxx = griddim;
  else maxx = gridPointX+2;
  if(gridPointY < 2) miny = 0;
  else miny = gridPointY-2;
  if(gridPointY > griddim-2) maxy = griddim;
  else maxy = gridPointY+2;

  for(i=minx;i<maxx;i++) for(j=miny;j<maxy;j++) if(grid[i][j]>=0) {
	otherx = samplePointsX[grid[i][j]];
	othery = samplePointsY[grid[i][j]];
	dist = sqrt((otherx - newpoint[0])*(otherx - newpoint[0]) + (othery - newpoint[1])*(othery - newpoint[1]));
	if(dist < mindist) return 1;
	}

  return 0;
}

int image_to_grid(int x, float cell_size) { // used for Poisson disk sampling
  return (int)((float)x/cell_size);
}

// End of Poisson disc sampling functions

char *concat(char *s1, char *s2)
{
	// returns the result of concatenating two strings

	char *text = malloc(strlen(s1)+strlen(s2)+1); // +1 for the zero-terminator
	
	if(text == NULL) {printf("Memory problem!\n"); exit(0);}
	
	strcpy(text, s1);
	strcat(text, s2);
	
	return text;
}

char *zeroes(int loop_count, int maxloops)
{
	// creates a string of zeroes to prefix the number in the file name,
	// so files are animated in the correct order
	
	char *text = malloc(20+1);
	int i,j,k;
	
	if(text == NULL) {printf("Memory problem!\n"); exit(0);}
	
	strcpy(text, "");
	for(i=(int)maxi((int)(log(loop_count + 0.1)/log(10)),0);i<(int)(log(maxloops + 0.1)/log(10));i++) {
	// the maxi function is invoked to avoid attempted calculation of log(0),
	// which the compiler equates to some negative number of large magnitude;
	// 0.1 is added to avoid rounding errors (e.g. to avoid log10(1000) = 2.99999)
		if(i>5) {printf("Too many zeroes!\n"); exit(0);}
		text=concat(text,"0");
		}
	
	return text;
}

float ran0(long *idum)
{
	// a random number generator

	long k;
	float ans;

	*idum ^= MASK;
	k= (*idum)/IQ;
	*idum= IA*(*idum-k*IQ)-IR*k;
	if(*idum<0) *idum+=IM;
	ans= AM*(*idum);
	*idum ^= MASK;
	
	return ans;
}

float maxi(float a, float b)
{
	// returns the maximum of two numbers

	float ans;
	
	if(a>=b) ans=a;
	else ans=b;
	
	return ans;
}

float mini(float a, float b)
{
	// returns the minimum of two numbers

	float ans;
	
	if(a<=b) ans=a;
	else ans=b;
	
	return ans;
}

void shuffle(int *sequence,int length,long *idum)
{
	// randomly shuffles a sequence of integers

	int i,j,buff;
	 
	for(i=length-1;i>0;i--) {
		do {j=(int)(ran0(idum)*(i+1));} while(j==i+1);
		buff=sequence[i];
		sequence[i]=sequence[j];
		sequence[j]=buff;
		}
}

int plusmod(int num, int max)
{
	// increments num (modulo max)

	int ans = num+1;

	if(ans>max) ans=0;

	return ans;
}

int sign(float x)
{
	if (x > 0) return 1;
	if (x < 0) return -1;
	return 0;
}

void shell_one(int *a, int dim, bool ascending)
{ // sorts the array a
	int i,j,inc;
	int v;
	
	inc=1;
	do {
		inc*=3;
		inc++;
	}while (inc<=dim);

	if(ascending==1) do {
		inc/=3;
		for(i=inc;i<dim;i++) {
			v=a[i];
			j=i;
			while(a[j-inc]>v) {
				a[j]=a[j-inc];
				j-=inc;
				if(j<inc) break;
			}
			a[j]=v;
		}
	} while(inc>1);

	else do {
		inc/=3;
		for(i=inc;i<dim;i++) {
			v=a[i];
			j=i;
			while(a[j-inc]<v) {
				a[j]=a[j-inc];
				j-=inc;
				if(j<inc) break;
			}
			a[j]=v;
		}
	} while(inc>1);
}

void shell_two(float *a, int *b, int dim, bool ascending)
{ // sorts the arrays a and b in with respect to a
	int i,j,inc;
	float v;
	int w;
	
	inc=1;
	do {
		inc*=3;
		inc++;
	}while (inc<=dim);

	if(ascending==1) do {
		inc/=3;
		for(i=inc;i<dim;i++) {
			v=a[i];
			w=b[i];
			j=i;
			while(a[j-inc]>v) {
				a[j]=a[j-inc];
				b[j]=b[j-inc];
				j-=inc;
				if(j<inc) break;
			}
			a[j]=v;
			b[j]=w;
		}
	} while(inc>1);

	else do {
		inc/=3;
		for(i=inc;i<dim;i++) {
			v=a[i];
			w=b[i];
			j=i;
			while(a[j-inc]<v) {
				a[j]=a[j-inc];
				b[j]=b[j-inc];
				j-=inc;
				if(j<inc) break;
			}
			a[j]=v;
			b[j]=w;
		}
	} while(inc>1);
}

int weighted_random (float *weights, int options, long *idum)
{
	float totalWeight = 0; // this stores sum of weights of all elements before current
	int selected; // currently selected element
	float weight;
	float rnd;
	int i, j;
	
	for(i=0; i<options; i++) {
		weight = weights[i]; // weight of current element
		rnd = ran0(idum)*(totalWeight + weight); // random value
		if (rnd >= totalWeight) // probability of this is weight/(totalWeight+weight)
			selected = i; // it is the probability of discarding last selected element and selecting current one instead
		totalWeight += weight; // increase weight sum
		}

	return selected; // when iterations end, selected is some element of sequence. 
}

float expdev(long *idum)
{ // Returns an exponentially distributed, positive, random deviate of unit mean, using ran1(idum) as the source of uniform deviates.
	float dum;

	do {dum=ran0(idum);} while (dum == 0.0);

	return -log(dum);
}

float gasdev(long *idum)
// Returns a normally distributed deviate with zero mean and unit variance
{
	static int iset=0;
	static float gset;
	float fac,rsq,v1,v2;
	
	if (iset == 0) {
		do {
			v1=2.0*ran0(idum)-1.0;
			v2=2.0*ran0(idum)-1.0;
			rsq=v1*v1+v2*v2;
			} while (rsq >= 1.0 || rsq == 0.0);
		fac=sqrt(-2.0*log(rsq)/rsq);
		gset=v1*fac;
		iset=1;
		return v2*fac;
		}
	else {
		iset=0;
		return gset;
		}
}

float truncated_log_gasdev(long *idum, float max)
// Returns a deviate from a log-normal distribution with unit mean, truncated at "max"
{
	float result;

	do {
		result = exp(gasdev(idum) - 0.5);
		} while(result > max);

	return result;
}

void SampleWithoutReplacement(int N, int n, int *samples, long *idum)
// returns a vector of n numbers between 0 and N-1 (inclusive)
{
	int t = 0; // total input records dealt with
	int m = 0; // number of items selected so far
	double u;

	while (m < n) {
		u = ran0(idum);

		if ((N - t)*u >= n - m) t++;
		else {
			samples[m] = t;
			t++;
			m++;
			}
		}
}