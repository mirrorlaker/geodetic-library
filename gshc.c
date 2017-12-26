/* Generate Gravity Gradient from Spherical Harmonic Coefficients
*	Junjun Yang
*	May 24, 2017	*/
#include<stdlib.h>
#include<argp.h>
#include<math.h>
#include<string.h>
#include<errno.h>
#include<omp.h>

#define TRUE 1
#define FALSE 0

int countlines(char *filename);
void gshc(int FLG, double *CNM, double *SNM, int NMIN, int NMAX,
	double LON, double PHI, double HT, double *TTT,
	double *P_bar_nm, double *P_bar_nm_d1, double *P_bar_nm_d2, double *SQT,
	double *cos_m_lambda, double *sin_m_lambda, double GM, double AE);

/* START PARSING PROGRAM ARGUMENTS */
const char *argp_program_version =
	"gshc 0.5.1";
const char *argp_program_bug_address =
	"<osujjy@gmail.com>";

/* A description of the arguments we accept. */
static const char args_doc[] =
	"FILE";
	
/* Program documentation. */
static const char doc[] =
	"gshc -- Compute the value, the 3 first-order derivatives, and the 6 second-order derivatives of "
	"disturbing geopotential model, given by truncated spherical harmonic model, at points stored in FILE."
	"\n\nInput: Geodetic longitude, latitude in degree, and geodetic height in meter."
	"\n  Format: \"%lf %lf %lf\". Example: 123.45 -56.78 8901.234"
	"\n\nOutput: Disturbing gravitational potential V and its derivatives "
	"in spherical North(n)-East(e)-Down(d) frame."
	"\n  V [m^2/s^s]"
	"\n  Vn, Ve, Vd [mGal]"
	"\n  Vnn, Vne, Vnd, Vee, Ved, Vdd [Eotvos]"
	"\n\nOptions:"
	"\vWritten by Junjun Yang.\nMay 24, 2017 at IGG.";

/* The options we understand. */
static const struct argp_option options[] = {
	{"verbose",	'v',	0,		0,	"Produce verbose output",1},
	{"gravitation",	'g',0,		0,	"Compute gravitational potential rather than disturbing potential",1},
	{"check",	'c',	0,		0,	"Check spherical harmonic coefficient file",1},
	{"output",	'o',	"FILE",	0,	"Output to FILE instead of standard output",1},
	{"list",	'l',  "LIST",	0,	"List output, default is -l vgd. Options for LIST:"
									"\nc: long lat height"
									"\nv: V"
									"\ng: Vn Ve Vd"
									"\nd: Vnn Vne Vnd Vee Ved Vdd",1},
	{"format",	'f',  "FORMAT",	0,	"Define output format",1},
	{"human-readable",	'h', 0,	0,	"without -f, print human readable format \"%11.3lg\"",1},
	{"model",	'm',	"FILE",	0,	"Load model, default is EGM2008 spherical harmonic coefficients",1},
	{0,0,0,0, "Getting help",2},
	{0}
};

/* Used by main to communicate with parse_opt. */
struct arguments {
	char *args[1];
	int verbose, check, gravitation, human;
	char *output_file;
	char *output_list;
	char *output_format;
	char *gravitational_model;
};

/* Argp parser function */
static error_t parser(int key, char *arg, struct argp_state *state) {
	struct arguments *arguments = state->input;

	switch (key) {
		case 'v':
			arguments->verbose = 1;
			break;
		case 'c':
			arguments->check = 1;
			break;
		case 'g':
			arguments->gravitation = 1;
			break;
		case 'o':
			arguments->output_file = arg;
			break;
		case 'l':
			arguments->output_list = arg;
			break;
		case 'f':
			arguments->output_format = arg;
			break;
		case 'h':
			arguments->human = 1;
			break;
		case 'm':
			arguments->gravitational_model = arg;
			break;
		case ARGP_KEY_ARG:
			// If too many arguments
			if (state->arg_num >= 1)
				argp_usage(state);
			
			arguments->args[state->arg_num] = arg;
			break;
		case ARGP_KEY_END:
			// If no enough arguments
			if (state->arg_num < 1)
				argp_usage(state);
			break;
		default:
			return ARGP_ERR_UNKNOWN;
	}
	return 0;
}

/* Argp parser */
static struct argp argp={options,parser,args_doc,doc};
/* END PARSING PROGRAM ARGUMENTS */

int main(int argc, char **argv)
{
	/* START PARSING PROGRAM ARGUMENTS */
	struct arguments arguments;
	/* default options */
	arguments.verbose = 0;
	arguments.check = 0;
	arguments.gravitation = 0;
	arguments.output_file = NULL;
	arguments.output_list = NULL;
	arguments.output_format = NULL;
	arguments.human = 0;
	arguments.gravitational_model = "EGM2008_to2190_TideFree";
	
	argp_parse(&argp,argc,argv,0,0,&arguments);
	/* END PARSING PROGRAM ARGUMENTS */
	
	/* Check arguments for the option -l */
	int j, length_output_list;
	if (arguments.output_list != NULL){
		length_output_list = strlen(arguments.output_list);
			for (j=0; j < length_output_list; j++) {
				switch(arguments.output_list[j]) {
					case 'c':
						break;
					case 'v':
						break;
					case 'g':
						break;
					case 'd':
						break;
					default:
						printf("Error! ****** output list unknown!\n");
						exit(202);
				}
			}
	}
	
	/* Check whether files exist */
	FILE *ptr_shc;
	if((ptr_shc = fopen(arguments.gravitational_model,"r")) == NULL) {
		printf("Error! ****** The gravitational model does not exist!\n");
		return ENOENT;
	}
	FILE *ptr_geodetic;
	if((ptr_geodetic = fopen(arguments.args[0], "r")) == NULL) {
		printf("Error! ****** The file for computation points does not exist!\n");
		fclose(ptr_shc);
		return ENOENT;
	}
	FILE *ptr_g;
	if (arguments.output_file != NULL) {
		if((ptr_g = fopen(arguments.output_file, "w")) == NULL) {
			printf("Error! ****** Open output file failed!\n");
			fclose(ptr_shc);
			fclose(ptr_geodetic);
			return ENOENT;
		}
	}
	else{
		ptr_g = stdout;
	}
	
	/* Check the degree and order of the shc file */
	int shc_nmin,shc_mmin,shc_nmax,shc_mmax;
	fscanf(ptr_shc,"%d%d",&shc_nmin,&shc_mmin);
	fseek(ptr_shc,-101,SEEK_END);
	fscanf(ptr_shc,"%d%d",&shc_nmax,&shc_mmax);
	fseek(ptr_shc,0,SEEK_SET);
	if (arguments.check) {
		if (arguments.verbose == 1) {
			printf("Checking spherical harmonic coefficient files: %s ...\n",arguments.gravitational_model);
		}
		if((shc_nmin != 2) | (shc_mmin != 0) | (shc_nmax != shc_mmax) | (shc_nmax < 2) ) {
			printf("Error! ****** [N][M] does not start from [2][0], or Nmax != Mmax\n");
			return ENOENT; // ?? other code
		}
		if (countlines(arguments.gravitational_model) != ((shc_nmax+1)*(shc_nmax+2)/2 - (shc_nmin)*(shc_nmin+1)/2)) {
			printf("Error! ****** Spherical harmonic coefficients incomplete!\n");
			return ENOENT; // ?? other code
		}
		if(shc_nmax > 2190) {
			printf("Attention! ****** NMAX > 2190. If it is not the case, please check your code.\n");
		}
		if (arguments.verbose == 1) {
			printf("--Done!\n");
		}
	}
	int NMIN = 0;					// the min degree that the computation will use
	int NMAX = shc_nmax;
	int ND = (NMAX+1)*(NMAX+2)/2;	// dimension of vectors that store spherical harmonic coefficients
	
	/* Set parameters for normal field */    // ?? option to change normal field
	double GM=3.986004415e14;
	double C2=108262.982131e-8;			// J2 for WGS84 computed from J2=-C20=-C20bar*sqrt(5), need to verify its value
	double AE=6378136.3;				// a in the EGM2008 model
	double A84=6378137.0;				// Semi-major axis for WGS84
	double RF=298.257223563;			// reciprocal of flattening for WGS84
	/* Useful constants */
	double DTR=M_PI/180.0; // use M_PIl if long double pi is preferred. It is only available if _GNU_SOURCE is defined
	double ESQ=2/RF-1.0/pow(RF,2.0);			// square of the first eccentricity
	int i;
	/* process bar */
	int process_rate = 0;
	int process_rate_last = 0;
	char process_buff[102] = {0};
	const char* process_arrow = "-\\|/";
	
	/* Read spherical harmonic coefficients */
	if (arguments.verbose == 1) {
		printf("Reading spherical harmonic coefficients from %s ...\n",arguments.gravitational_model); 
	}
	/* int *N, *M; */
	double *CNM, *SNM;	// {2i5, 2d25.15, 2d20.10}
/* 	N = (int *) calloc (ND, sizeof(int));
	M = (int *) calloc (ND, sizeof(int)); */
	CNM = (double *) calloc (ND, sizeof(double));
	SNM = (double *) calloc (ND, sizeof(double));
	
	char shc_buff[110]; // work buff for converting exponent from D to E, because C language doesnot provide this format
	/* check the format of exponent */
	int mark_exp_format = 0; 
	fscanf(ptr_shc,"%*d%*d");
	fscanf(ptr_shc,"%s",shc_buff);
	fseek(ptr_shc,0,SEEK_SET);
	for(j=0;j<110;j++) {
		if(shc_buff[j] == 'e' | shc_buff[j] == 'E') {
			mark_exp_format = 1;
			break;
		}
		else if(shc_buff[j] == 'D' | shc_buff[j] == 'd'){
			mark_exp_format = 2;
			break;
		}
	}
	if(mark_exp_format == 0){
		printf("Error! ****** shc exponent format is not one of following: D,e,E,d\n");
		exit(201);
	}
	else if(mark_exp_format == 1)
	{
		for(i=3; i<ND; i++) {
			fscanf(ptr_shc,"%*d%*d%lf%lf%*f%*f",&(CNM[i]),&(SNM[i]));
			// process bar
			if(arguments.verbose == 1) {
				process_rate = (i+1)*100/ND;
				if (process_rate > process_rate_last) {
					for(j=process_rate_last;j<process_rate;j++)
						process_buff[j] = '=';
					process_buff[process_rate] = '-';
					process_buff[process_rate + 1] = '\0';
					process_rate_last = process_rate;
					printf("[%-100.100s],%d%%,[%c]\r",process_buff,process_rate,process_arrow[process_rate%4]);
					fflush(stdout);
				}
			}
		}
	}
	else if(mark_exp_format == 2){
		for(i=3; i<ND; i++) {
			fscanf(ptr_shc,"%*d%*d");
			
			fscanf(ptr_shc,"%s",shc_buff);
			for(j=0;j<110;j++) {
				if(shc_buff[j] == 'D')
					break;
			}
			shc_buff[j]='E';
			sscanf(shc_buff,"%lf",&(CNM[i]));
			
			fscanf(ptr_shc,"%s",shc_buff);
			for(j=0;j<110;j++) {
				if(shc_buff[j] == 'D')
					break;
			}
			shc_buff[j]='E';
			sscanf(shc_buff,"%lf",&(SNM[i]));
			
			fscanf(ptr_shc,"%*s%*s");
			
			// process bar
			if(arguments.verbose == 1) {
				process_rate = (i+1)*100/ND;
				if (process_rate > process_rate_last) {
					for(j=process_rate_last;j<process_rate;j++)
						process_buff[j] = '=';
					process_buff[process_rate] = '-';
					process_buff[process_rate + 1] = '\0';
					process_rate_last = process_rate;
					printf("[%-100.100s],%d%%,[%c]\r",process_buff,process_rate,process_arrow[process_rate%4]);
					fflush(stdout);
				}
			}
		}
	}
	if(arguments.verbose == 1){
		printf("\n");
	}
	fclose(ptr_shc);
		
	/* Define zero- and first-degree coefficients */
	CNM[0] = 1.0;
	SNM[0] = 0.0;
	CNM[1] = 0.0;
	SNM[1] = 0.0;
	CNM[2] = 0.0;
	SNM[2] = 0.0;
/* 	N[0] = 0;
	M[0] = 0;
	N[1] = 1;
	M[1] = 0;
	N[2] = 1;
	M[2] = 1; */
	if (arguments.verbose == 1){
		printf("--degree ranges from %d to %d.\n",shc_nmin,NMAX);
		printf("--Done!\n");
	}
		
	/* remove normal field */
	if (!arguments.gravitation) {		
		if (arguments.verbose == 1){
			printf("Removing normal field...\n");
		}
		double C2N[5];
		double CB2, CB4, CB6, CB8, CB10;
		
		for(i=2; i<=5; i++)
			C2N[i-1]=pow(-1.0,(double)i)*(3.0*pow(ESQ,(double)i)/((2.0*i+1.0)*(2.0*i+3.0)))*(1.0-i+5.0*i*C2/ESQ);
		
		CB2=-1.0*C2/sqrt(5.0);
		CB4=C2N[1]/3.0;
		CB6=C2N[2]/sqrt(13.0);
		CB8=C2N[3]/sqrt(17.0);
		CB10=C2N[4]/sqrt(21.0);
		
		// do not remove the following 6 CNM if compute gravitational potential rather than disturbing potential
		// ?? option: add centrifugal potential if compute gravity potential, -g pre-requsit, yes no
		CNM[0]=0.0;
		CNM[3]  -= CB2;
		CNM[10] -= CB4;
		CNM[21] -= CB6;
		CNM[36] -= CB8;
		CNM[55] -= CB10;
		if (arguments.verbose == 1){
			printf("--Done!\n");
		}
	}
	
	/* Read geodetic coordinates */ // ?? option: geodetic, geocentric
	if (arguments.verbose == 1){
		printf("Loading computation points from %s ...\n",arguments.args[0]);
	}
	int NUM_PT;
	NUM_PT = countlines(arguments.args[0]);
	double *LON, *LAT, *GH;
	LON = (double *) calloc (NUM_PT, sizeof(double));
	LAT = (double *) calloc (NUM_PT, sizeof(double));
	GH  = (double *) calloc (NUM_PT, sizeof(double));
	for(i=0; i<NUM_PT; i++)
		fscanf(ptr_geodetic, "%lf%lf%lf", &(LON[i]),&(LAT[i]),&(GH[i]));
	fclose(ptr_geodetic);
	if (arguments.verbose == 1){
		printf("--Done!\n");
		printf("--Loaded %d points.\n",NUM_PT);
	}
	
	/* Transform geodetic coordinates to spherical coordinates */
	if (arguments.verbose == 1){
		printf("Transforming from geodetic to geocentric coordinates...\n");
	}
	double *PHI, *HT;
	double phir, rn, p, z, r;
	PHI = (double *) calloc (NUM_PT, sizeof(double));
	HT  = (double *) calloc (NUM_PT, sizeof(double));
	for(i = 0; i < NUM_PT; i++) {
		phir=LAT[i]*DTR;
		rn=A84/sqrt(1.0-ESQ*pow(sin(phir),2.0));	// radius of curvature                           
		p=(rn+GH[i])*cos(phir);						// distance from polar axis                                                   
		z=(rn*(1.0-ESQ)+GH[i])*sin(phir);
		PHI[i]=atan(z/p)/DTR;						// geocentric latitude (deg) ?? check atan
		r=sqrt(p*p+z*z);
		HT[i]=r-AE;									// height relative to the AE sphere
	}
	if (arguments.verbose == 1){
		printf("--Done!\n");
	}
	
	// ?? sort lonlatgh
	
	/* Compute */
	if (arguments.verbose == 1){
		printf("Computing...\n");
	}
	int flg;
	double r2, theta, cost, sint, cott, GVT[3], GDT[6], TTT[10];
	double *P_bar_nm, *P_bar_nm_d1, *P_bar_nm_d2, *SQT, *cos_m_lambda, *sin_m_lambda;
	P_bar_nm  = (double *) calloc (ND, sizeof(double));
	P_bar_nm_d1 = (double *) calloc (ND, sizeof(double));
	P_bar_nm_d2 = (double *) calloc (ND, sizeof(double));
	SQT   = (double *) calloc (2*NMAX+2, sizeof(double));
	cos_m_lambda    = (double *) calloc (NMAX+1, sizeof(double));
	sin_m_lambda    = (double *) calloc (NMAX+1, sizeof(double));
	
	// process bar
	process_rate = 0;
	process_rate_last = 0;
	
	/* Option to choose output and provide format control */
	char output_format_c[256];
	char output_format_v[256]; // ?? why 20
	char output_format_g[256];
	char output_format_d[256];
	
	if (arguments.output_format == NULL){ // no format defined. automatic define
		if (arguments.output_list == NULL){ // default output_list == "vgd"
			if (arguments.human == 0)
				arguments.output_format = 
				"%25.15le%25.15le%25.15le%25.15le%25.15le%25.15le%25.15le%25.15le%25.15le%25.15le\n";
			else
				arguments.output_format = 
				"%11.3lg%11.3lg%11.3lg%11.3lg%11.3lg%11.3lg%11.3lg%11.3lg%11.3lg%11.3lg\n";
		}
		else{
			if (arguments.human == 0) {
				sprintf(output_format_c,"%s","%16.10lf%16.10lf%16.10lf");
				sprintf(output_format_v,"%s","%25.15le");
				sprintf(output_format_g,"%s","%25.15le%25.15le%25.15le");
				sprintf(output_format_d,"%s","%25.15le%25.15le%25.15le%25.15le%25.15le%25.15le");
			}
			else{
				sprintf(output_format_c,"%s","%11.3lg%11.3lg%11.3lg");
				sprintf(output_format_v,"%s","%11.3lg");
				sprintf(output_format_g,"%s","%11.3lg%11.3lg%11.3lg");
				sprintf(output_format_d,"%s","%11.3lg%11.3lg%11.3lg%11.3lg%11.3lg%11.3lg");
			}
		}
	}
	else{ // format is defined in command line, correct \n 
		int length_output_format = strlen(arguments.output_format); // length does not include the '\0' at tail
		if((arguments.output_format[length_output_format-2] == '\\') & (arguments.output_format[length_output_format-1] == 'n')){
			arguments.output_format[length_output_format-2] = '\n';
			arguments.output_format[length_output_format-1] = '\0';
			length_output_format--;
		}
		else if ( (arguments.output_format[length_output_format-1] == 'n') &
		((arguments.output_format[length_output_format-2] == 'f') | (arguments.output_format[length_output_format-2] == 'e') | 
		 (arguments.output_format[length_output_format-2] == 'E') | (arguments.output_format[length_output_format-2] == 'g') |
		 (arguments.output_format[length_output_format-2] == 'G') | (arguments.output_format[length_output_format-2] == 'a') |
		 (arguments.output_format[length_output_format-2] == 'A')) ){
			arguments.output_format[length_output_format-1] = '\n'; // The format send to command line is assumed without '' or ""
			// and casue the end of format string to be n rather than \ n. what if they just want a n at the end ??
		}
		if (arguments.output_list != NULL){
			// ?? format is defined, list is defined. extract %f and save them to output_format_*[20]
			// temperarily set it to default
			if (arguments.human == 0) {
				sprintf(output_format_c,"%s","%16.10lf%16.10lf%16.10lf");
				sprintf(output_format_v,"%s","%25.15le");
				sprintf(output_format_g,"%s","%25.15le%25.15le%25.15le");
				sprintf(output_format_d,"%s","%25.15le%25.15le%25.15le%25.15le%25.15le%25.15le");
			}
			else{
				sprintf(output_format_c,"%s","%11.3lg%11.3lg%11.3lg");
				sprintf(output_format_v,"%s","%11.3lg");
				sprintf(output_format_g,"%s","%11.3lg%11.3lg%11.3lg");
				sprintf(output_format_d,"%s","%11.3lg%11.3lg%11.3lg%11.3lg%11.3lg%11.3lg");
			}
		}
	}
	for(i=0; i<NUM_PT; i++) {
		if (i == 0)
			flg = TRUE;
		else if (PHI[i] != PHI[i-1])
			flg = TRUE;
		else
			flg = FALSE;
		gshc(flg,CNM,SNM,NMIN,NMAX,LON[i],PHI[i],HT[i],TTT,P_bar_nm,
			P_bar_nm_d1,P_bar_nm_d2,SQT,cos_m_lambda,sin_m_lambda,GM,AE);
/* 		TTT contains:
		  T			[m2/s2]
		 dT/dLam	[m2/s2],	dT/dthet		[m2/s2],	dT/dr		[m/s2],
		d2T/dlam2	[m2/s2],   d2T/dlam*dthet	[m2/s2],   d2T/dlam*dr	[m/s2],
		d2T/dthet2	[m2/s2],   d2T/dthet*dr		[m/s2],	   d2T/dr2		[1/s2] */
		
		/* Compute local NED derivatives */
		r = HT[i]+AE;
		r2 = r*r;
		theta = (90.0-PHI[i])*DTR;        // polar distance (rad)
		cost = cos(theta);
		sint = sin(theta);
		cott = cost/sint;							// xyz -> ned
		GVT[0] = -TTT[2]/r;		// Tx
		GVT[1] = TTT[1]/sint/r;	// Ty
		GVT[2] = -TTT[3];		// Tz
		GDT[0] = TTT[3]/r+TTT[7]/r2;							// Txx
		GDT[1] = (TTT[1]*cott-TTT[5])/r2/sint;					// Txy
		GDT[2] = -TTT[2]/r2+TTT[8]/r;							// Txz
		GDT[3] = TTT[3]/r+(TTT[2]*cott+TTT[4]/pow(sint,2))/r2;	// Tyy
		GDT[4] = (TTT[1]/r2-TTT[6]/r)/sint;						// Tyz
		GDT[5] = TTT[9];			                            // Tzz
		
		/* Change units to mGal & Eotvos */
		for(j=0;j<3;j++) {
			GVT[j] = GVT[j]*1e5;
		}
		for(j=0;j<6;j++) {
			GDT[j] = GDT[j]*1e9;
			}
			
		if (arguments.output_list == NULL)
			fprintf(ptr_g,arguments.output_format,
				TTT[0],GVT[0],GVT[1],GVT[2],GDT[0],GDT[1],GDT[2],GDT[3],GDT[4],GDT[5]);
		else {
			for (j=0; j < length_output_list; j++) {
				switch(arguments.output_list[j]) {
					case 'c':
						fprintf(ptr_g,output_format_c,LON[i],LAT[i],GH[i]);
						break;
					case 'v':
						fprintf(ptr_g,output_format_v,TTT[0]);
						break;
					case 'g':
						fprintf(ptr_g,output_format_g,GVT[0],GVT[1],GVT[2]);
						break;
					case 'd':
						fprintf(ptr_g,output_format_d,GDT[0],GDT[1],GDT[2],GDT[3],GDT[4],GDT[5]);
						break;
				}
			}
			fprintf(ptr_g,"\n");
		}
		
		// process bar
		if((arguments.verbose == 1) & (arguments.output_file != NULL)){
			process_rate = (i+1)*100/NUM_PT;
			if (process_rate > process_rate_last) {
				for(j=process_rate_last;j<process_rate;j++)
					process_buff[j] = '=';
				process_buff[process_rate] = '-';
				process_buff[process_rate + 1] = '\0';
				process_rate_last = process_rate;
			}
			printf("[%-100.100s],%d%%,[%c]\r",process_buff,process_rate,process_arrow[i%4]);
			fflush(stdout);
		}
	}
	// process bar
	if((arguments.verbose == 1) & (arguments.output_file != NULL)){
		printf("[%-100.100s],%d%%,[%c]\r",process_buff,process_rate,'-');
		printf("\n");
	}
	if (arguments.output_file != NULL) {
		fclose(ptr_g);
	}
	if (arguments.verbose == 1){
		printf("--Finished!\n");
	}
	
	// result is already close to fortran.  Fortran: real 2m31.145s, user 2m30.547s, sys 0m0.406s
										//  C      : real 3m54.154s, user 3m49.547s, sys 0m4.359s
	// check pbar1 pbar2 formula
	// c speed is slower than fortran
	// quality control.
	// check denominator, atan, etc, may need to output PHI rather than LAT
	// sort c1000 -g -k1 -k2|awk '{if($7<9.6e5 || $7>10e5){print}}' >tmp
	// sort c1000 -g -k1 -k2| grep nan >>tmp
	// ?? openmp
	// ?? output in local spheroidal coordinates
	// complete readme.pdf
	
	/* END PROGRAM NORMALLY */
	exit(0);
}

int countlines(char *filename)
{
	FILE *ptr;
	if((ptr = fopen(filename, "r")) == NULL) {
		printf("Error! ****** countlines open file failed!\nexit!\n");
		exit(ENOENT);
	}
	int linec = 0;
	while (EOF != (fscanf(ptr,"%*[^\n]"), fscanf(ptr,"%*c")))
		linec++;
	fclose(ptr);
	return linec;
}

int indx(int n,int m)
{
	/* Index locator function
	input: degree n, order m; output: index in the C language vector Cnm Snm */
	return ( m + n*(n+1)/2 );
}

void lgndr(double t, int NMAX, double *P_bar_nm, double *SQT, int *IX)
{
/* 	t is the cosine of polar distance in radian,
	at which the Fully Normalized Associated Legendre
	Function is evaluated up to degree NMAX using a recursion on the order.
	SQT is a work vector of dimension 2*NMAX+2. SQT[0] is not used.
	From [1] to [2*NMAX+1], SQT[i] contains the square root of positive integers i.
	Set IX to 0 in the first call to this subroutine.
	The values of the function are stored with order increasing before degree:
	(n=0,m=0),(n=1,m=0),(n=1,m=1),(n=2,m=0),(n=2,m=1),(n=2,m=2),(n=3,m=0),(n=3,m=1),... */
	
	int n, m;
	double t2 = t*t;
	if (*IX == 0) {
	  *IX = 1;
	  int NSQRTS = NMAX+NMAX+1;
	  #pragma omp parallel for
	  for(m = 1; m <= NSQRTS; m++)
		  SQT[m] = sqrt(m);
	}
	
	P_bar_nm[0] = 1.0;						// P_bar_nm[indx(0,0)]
	P_bar_nm[1] = SQT[3]*t;					// P_bar_nm[indx(1,0)]
	P_bar_nm[2] = SQT[3]*sqrt(1.0 - t2);	// P_bar_nm[indx(1,1)]
	if(NMAX == 1)
		return;
	
	for(n = 2; n <= NMAX; n++)
	{
		for(m = 0; m < n-1; m++)
		{
			P_bar_nm[indx(n,m)] = 
				( (SQT[2*n-1]/SQT[n-m])*
				  (SQT[2*n+1]/SQT[n+m]) )*t*P_bar_nm[indx(n-1,m)] -
				( (SQT[2*n+1]/SQT[2*n-3])*
				  (SQT[n+m-1]/SQT[n+m])*
				  (SQT[n-m-1]/SQT[n-m]) )*P_bar_nm[indx(n-2,m)];
		}
		P_bar_nm[indx(n,n-1)] = SQT[2*n+1]*t*P_bar_nm[indx(n-1,n-1)];
		P_bar_nm[indx(n,n)] = sqrt((2*n+1.0)*(1-t2)/(2*n))*P_bar_nm[indx(n-1,n-1)];
	}
}

void gshc(int FLG, double *CNM, double *SNM, int NMIN, int NMAX,
	double LON, double phi, double HT, double *TTT,
	double *P_bar_nm, double *P_bar_nm_d1, double *P_bar_nm_d2, double *SQT,
	double *cos_m_lambda, double *sin_m_lambda, double GM, double AE)
{
	/* Compute the value, the 3 first-order derivatives and the 6 second-order
	 derivatives of disturbing geopotential model, given by truncated spherical
	 harmonic model, at a point.  That is:

	 T
	 dT/dlam, dT/dthet, dT/dr
	 d2T/dlam2, d2T/dlam*dthet, d2T/dlam*dr, d2T/dthet2, d2T/dthet*dr, d2T/dr2

	 The Legendre functions are computed for each point, unless a flag indicates
	 that the previously used latitude holds.

	 Input:
	 FLG = should be set to .TRUE. only if a new latitude is being used
	 CNM = vector of cosine spherical harmonic coefficients (C00, C10, C11, C20, etc.)
	 SNM = vector of sine spherical harmonic coefficients (S00, S10, S11, S20, etc.)
	 NMIN= minimum degree limit in summation (=2 for complete model)
	 NMAX= maximum degree limit in summation
	 LON = longitude of point [decimal degrees]
	 phi = geocentric latitude of point [decimal degrees]
	 HT  = radial height of point above sphere of radius AE [meters] 

	 Output:
	 TTT = vector of 10 gradients: T [m2/s2], dT/dLam [m2/s2], dT/dthet [m2/s2],	dT/dr [m/s2],
		   d2T/dlam2 [m2/s2], d2T/dlam*dthet [m2/s2], d2T/dlam*dr [m/s2], d2T/dthet2 [m2/s2],
		   d2T/dthet*dr [m/s2], d2T/dr2 [1/s2]

	 Work vectors: P_bar_nm, P_bar_nm_d1, P_bar_nm_d2, SQT, cos_m_lambda, sin_m_lambda

	 Note: CNM and SNM should have dimension (NMAX+1)*(NMAX+2)/2 and should refer to disturbing field 
		   P_bar_nm, P_bar_nm_d1, P_bar_nm_d2 should also have this dimension and will contain the values of the
		   fully normalized Legendre function and its first and second derivatives w.r.t. theta
		   SRT should have dimension 2*NMAX+1
		   cos_m_lambda and sin_m_lambda should have dimension NMAX+1 */

	// Normal field scale factors: GM, AE

	static int IX = 0; // flag for lgndr subroutine
	
	int i,n,m,n1;
	double DTR = M_PI/180.0;
	double r = AE+HT;
	double SCALE = GM/AE;
	double a_r = AE/r;

	/* Check if it's a new latitude */
	if(FLG) {
		double theta = (90.0 - phi)*DTR; // in radian
		double cost = cos(theta);
		double sint = sin(theta);
		double cott = cost/sint;

		/* Compute Legendre functions at latitude phi */
		lgndr(cost,NMAX,P_bar_nm,SQT,&IX);

		/* Compute the first and second derivatives of Legendre Function w.r.t. theta=90-phi */
		for(n=2;n<=NMAX;n++) {
			n1=n+1;
			#pragma omp parallel for
			for(m=0;m<=n;m++) {
				if(m == 0)
					P_bar_nm_d1[indx(n,0)]=-sqrt(n1*n/2.0)*P_bar_nm[indx(n,1)];
				else if(m == 1)
					P_bar_nm_d1[indx(n,1)]=-cott*P_bar_nm[indx(n,1)]+sqrt(2.0*n*n1)*P_bar_nm[indx(n,0)];
				else
					P_bar_nm_d1[indx(n,m)]=-m*cott*P_bar_nm[indx(n,m)]+sqrt((n-m+1)*(n+m))*P_bar_nm[indx(n,m-1)];
				
				P_bar_nm_d2[indx(n,m)]=-cott*P_bar_nm_d1[indx(n,m)]-(n*n1-m*m/pow(sint,2.0))*P_bar_nm[indx(n,m)];
			}
		}
	}
	
	/* Compute cos_m_lambda[m] = cos(m*lambda), sin_m_lambda[m] = sin(m*lambda) */
	double lambda = LON*DTR; // in radian
	sin_m_lambda[0] = 0.0;
	cos_m_lambda[0] = 1.0;
	sin_m_lambda[1] = sin(lambda);
	cos_m_lambda[1] = cos(lambda);
	for(i = 2; i <= NMAX; i++) {
		sin_m_lambda[i] = 2.0*cos_m_lambda[1]*sin_m_lambda[i-1] - sin_m_lambda[i-2];
		cos_m_lambda[i] = 2.0*cos_m_lambda[1]*cos_m_lambda[i-1] - cos_m_lambda[i-2];
	}
	
	/* Initialize accumulators */
	for(i = 0; i < 10; i++)
		TTT[i] = 0.0;

	/* Perform summation over degree and order to compute derivatives w.r.t. spherical coordinates */
	int k = NMIN*(NMIN+1)/2; //index of Cnm Snm in C vector
	double a_r_n1 = pow(a_r,(double)NMIN); //(a/r)^(n+1)
	double n1_r_a_r_n1,m_a_r_n1;
	double CCpSS,CSmSC,CCpSS_PB,CCpSS_PB1,CSmSC_PB;
	
	for(n = NMIN; n <= NMAX; n++) {
		a_r_n1 = a_r_n1*a_r;
		n1_r_a_r_n1 = (n+1)/r*a_r_n1;
		for(m = 0; m <= n; m++) {
			m_a_r_n1 = m*a_r_n1;
			
			CCpSS = CNM[k]*cos_m_lambda[m] + SNM[k]*sin_m_lambda[m];
			CSmSC = CNM[k]*sin_m_lambda[m] - SNM[k]*cos_m_lambda[m];

			CCpSS_PB  = CCpSS*P_bar_nm[k];
			CCpSS_PB1 = CCpSS*P_bar_nm_d1[k];
			CSmSC_PB  = CSmSC*P_bar_nm[k];

			TTT[0] +=              a_r_n1*CCpSS_PB;				// T
			TTT[1] -=            m_a_r_n1*CSmSC_PB;				// dT/dLamda
			TTT[2] +=              a_r_n1*CCpSS_PB1;			// dT/dTheta
			TTT[3] -=         n1_r_a_r_n1*CCpSS_PB;				// dT/dr
			TTT[4] -=          m*m_a_r_n1*CCpSS_PB;				// d^2T/dLamda^2
			TTT[5] -=            m_a_r_n1*CSmSC*P_bar_nm_d1[k];	// d^2T/dLamda*dTheta
			TTT[6] +=       m*n1_r_a_r_n1*CSmSC_PB;				// d^2T/dLamda*dr
			TTT[7] +=              a_r_n1*CCpSS*P_bar_nm_d2[k];	// d^2T/dTheta^2
			TTT[8] -=         n1_r_a_r_n1*CCpSS_PB1;			// d^2T/dTheta*dr
			TTT[9] += (n+2)/r*n1_r_a_r_n1*CCpSS_PB;				// d^2T/dr^2
			
			k++;
		}
	}
	/* Apply scale factors to get potential units */
	for(i=0;i<10;i++)
		TTT[i]=SCALE*TTT[i];
}
