/* Program Prot_evol.
   Needs target structure(s) and an alignment
   Attempts DNA mutations with given mutational bias.
   alpha, Z, energy are computed for the mutated sequence.
   The mutation is accepted or not according to selection criteria.
   The base composition is evaluated at every time step.
*/

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <stdlib.h>
#include <time.h>

char DIR_TNM[100]="/data/ubastolla/BIN/SSCPE/TNM_DATA";
char FILE_STR[200]="/data/ubastolla/BIN/SSCPE/structures.in";
char tnm[100]="/data/ubastolla/BIN/SSCPE/tnm";
char tnm_mut_para[100]="/data/ubastolla/BIN/SSCPE/Mut_para.in";

int Naa=20; // Amino acid number; you may consider gaps

int OPT_REG=1;  // Optimize the regularization parameter REG?
int CV_DOWN=4;  // Number of CV decreasing steps for declaring maximum
int REMUT=0;    // Repeat the fit of mutation parameters?
int SCORE_KL=1; // Score: KL (1) or likelihood (0)
int SCORE_EVO=0; // Consider in score also lik_flux
int LAMBDA_ANALYTIC=0; // Optimize Lambda by vanishing derivative

float REG_INI=0.02; // Initial value of regularization
float REG_STEP=0.02;
float REG_MIN=0.00;
float EPS_REG=0.01; // For computing Cv=-dscore/dReg
int N_REG=20;
float REG_DISTR=0.0; // Regularization of amino acid distr.
float Reg;

int repeat=1;
int clean_REM=1;

int IT_REMUT=7; // Attempt to change mutation distribution
int OPT_MUT=1;  // Determine P_glob fitting f_reg_aa
int it_max_mut=3; // Number of iterations for imposing sum_i Pi(a)/L=f_reg_aa(a)
float Lambda_start[2]={1.0,0.0}; // Initial value of Lambda
float LAMBDA_MF=2.0; // Lambda parameter for meanfield if OPT_LAMBDA=0
int IWT=0;      // For computing wild-type structural deformation
int NORMALIZE=1;  // Normalize the exponent of Psel?
float Lambda_iwt[2]={1.5,0.0}; // Lamdba for wt fitness of structural models
float Sta_IWT=0.25;
int RECOMPUTE_IWT=1;

// Computation of evolutionary likelihood
int SCORE_EVO_TMP=0;
float t_ave=0;
float ***flux_iab=NULL;


// Global variables for computing P_ia
float **exp_str=NULL;
short *aa_seq=NULL;
float *f_glob=NULL;
int N_SAM_LAM=11;
float LAM_MIN=0.1;
float LAM_STEP=0.15;
float **Lambda_sam; // 0.1 0.3 0.5 0.7 0.9 1.1 1.3 1.5 1.7 1.9 
float *score_sam;
int Sample_Lambda=1;

float PMIN=0.001; // Minimal allowed value of predicted a.a. frequencies
int PMIN_ZERO=0; // Set Pmin to zero when optimizing REG?
int SCORE_CV=1; // Score for optimizing REG: 1=CV 0=-|KL_mod-KL_reg|
int MF_COMP=1;  // Perform Mean-field computations?
int STR_COMP=1; // Perform RMSD and DE computations with TNM

// Thermodynamic parameters
float TEMP=1.0;
float sC1=0.065, sC0=0, sU1=0.130; // Configuration entropy
//float sC1=0.032, sC0=0, sU1=0.065; // Configuration entropy
float SEC_STR=1.0; // Coefficient for local interactions // 1

float N_stab_thr=5;  //2(20 sel) 3(18 sel) 4(11 sel) 5(4 sel) 8(3 sel) 
float f_stab_thr=0.5;
float N_stab_thr_2;
float thr_stab=0.95; // minimum sum of stability scores for any MSA sequence

int UPDATE_LAMBDA=1; // Update initial value of lambda for each REG
//
int PRINT_TN=0; // Print Tajima-Nei divergence?
int ALL_MUTS=0; // Compute the effect of all possible DNA mutations?

int ini_print;
int ini_lik=1;
double lik_const=0;

extern int N_SS;
char SEC_EL_cp[16]="CTS01234HcXABEYZ"; // S=T here

// MSA
char file_ali[100]="", name_ali[100]="";
int N_seq=0;
float N_stab_target=0;
double DG_ave_target=0;
float Seq_id_ave_target=0;
float *f_reg_aa=NULL;
float **lf_ia=NULL;

// Structural domain
int ini_dom, end_dom, L_dom;

#define FILE_CODE_DEF "gen_code_ATGC.in"
// #define FILE_ENE_DEF  "energy.in"
#define VBT 0 // Verbatim

#include "REM.h"
#include "coord.h"
#include "gen_code.h"
#include "allocate.h"
#include "protein3.h"
#include "read_pdb.h"
#include "random3.h"           /* Generating random numbers */
#include "mutation.h"
#include "codes.h"
#include "input.h"
#include "meanfield.h"
#include "get_para_pop_dyn.h"
#include "subst_exhaustive.h"
#include "read_str_mut.h"
#include "jtt.h"
#include "wag.h"
#include "lg.h"
#include "Print_meanfield.h"
#include "fits.h"
#include "read_ali.h"
//#include "Sec_str_16_NR.h";
/*
   Energy parameters and genetic codes are given in header files
   energy_BKV.h gen_code.h
*/

//#define FILE_IN "Prot_evol.in"
#define N_CHAR 300          // Max. length of file names
#define EXT_MSA ".MSA.fasta" // Extension for MSA
#define EXT_SEL ".sel.dat"  // Extension for selection statistics
#define EXT_AVE ".ave.dat"  // Extension for output file with averages
#define EXT_DNA ".dna.dat"  // Extension for output file with dna statistics
#define EXT_OUT ".stab.dat" // Extension for output file with folding
                            // thermodynamics and fitness
#define EXT_FIN ".final.dat" // Extension for final results of simulation
#define EXT_SUM ".summary.dat"

int PRINT_MSA=1;
int MAX_MSA=1000;
int nseq_msa=0;

// Global variables in externals.h
// Mutation parameters
/*char coded_aa[64]="FFLLLLLLIIIMVVVVSSSSPPPPTTTTAAAAYYHHQQNNKKDDEECCWRRRRSSRRGGGG***";
  char *codon[64]={"TTT","TTC","TTA","TTG","CTT","CTC","CTA","CTG","ATT","ATC","ATA","ATG","GTT","GTC","GTA","GTG","TCT","TCC","TCA","TCG","CCT","CCC","CCA","CCG","ACT","ACC","ACA","ACG","GCT","GCC","GCA","GCG","TAT","TAC","CAT","CAC","CAA","CAG","AAT","AAC","AAA","AAG","GAT","GAC","GAA","GAG","TGT","TGC","TGG","CGT","CGC","CGA","CGG","AGT","AGC","AGA","AGG","GGT","GGC","GGA","GGG","TAA","TAG","TGA"};*/ 
#define MUTPAR 10
float mut_par[MUTPAR];
int NPAR=6;
// contact definition
float cont_thr;
// contact statistics
/*float *Cont_L=NULL;
float *CNc1_L=NULL;
float *CNc2_L=NULL;
float *Cnc1_L=NULL;
float **nc_nc_L=NULL;*/

// Contact energy function
char AA_code[21];
//extern char *AA_PAML;
char AA_PAML_1[21]="ARNDCQEGHILKMFPSTWYV-";
float Econt[21][21];
float **Econt_T=NULL, T_ratio=1;
//float **E_loc_over_T=NULL;
float hydro[21];
//char SEC_EL[16];

float Conf_entropy, K2Thr;
// general
int Verbose=0;
// AA distr
float Entr_ave=0, Entr_reg=0;

// Input codes
void Read_ene_par(char *, float **interactions);
int Read_ene_new(char *, float **interactions);
char *Read_sequence(int *len_dna, int *nseq, int **ini_seq, int **len_seq,
		    char *inputseq);
int Read_pdb_files(char ***name_pdb, char ***chain_pdb, char *file_pdb_list);
void Run_tnm(char *file_pdb, char *chain, char *pdbid);
int Copy_pdb_name(char ***name_pdb, char ***chain_pdb,
		  char *file_pdb, char *chain);
//float **Enlarge_matrix(float **C_nat, int L, int *ali_PDB, int len_amm);

unsigned long randomgenerator(void);
float *Get_counts(short *seq, int L, int NAA);
int Match_dna(char *dna_seq, int nseq, int *ini_seq, int *len_seq,
	      struct protein pdb, char **codon, char *coded_aa);

// Output
void Output_name(char *file_name, char *dir_out, char *prot_name,
		 float TEMP, float sU1, float sC1,  int MEANFIELD,
		 char *MODEL, float LAMBDA, int OPT_LAMBDA,
		 int NEUTRAL, int N_pop, float *mut_par);
FILE *open_file(char *, char *, short *, int, char *fit_def);
void Remove_extension(char *name1, char *name2);
void Remove_path(char *name);
int Print_dna(char *, FILE *, int);
void Print_ave(FILE *file_ave,
	       long it_sum, long t_indip,
	       int N_pop,
	       double f_sum, double f_dev,
	       double E_sum, double E_dev,
	       double DG_sum, double DG_dev,
	       long num_syn_subst, long num_aa_subst,
	       long num_syn_mut, long num_aa_mut,
	       float seq_entr,
	       struct load mut_load,
	       struct load trans_load);
void Print_TN_div(short *aa_seq, short *aa_seq0, int L,
		  int num_aa_subst, FILE *file_out);
void Print_seq(FILE *file_msa, short *aa_seq, int len_amm, int *nseq_msa,
	       float DeltaG);


void Print_matrix(struct protein target);
FILE *Open_summary(char *name_file, int L, struct REM E_wt,
		   float *mut_par, short *aa_seq, int npdb);
void Get_mean(double *ave, double *err, float sum, float dev,
	      long it_sum, float t_indep);
void Record_evo(float *fit_evo, double fitness,
		float *DG_evo,  double DG,
		float *Enat_evo, double E_nat, int *nsam);
void Print_final(char *name_file, long it_sum,
		 float TEMP, float sU1, float sC1, float sC0,
		 int MEANFIELD, float LAMBDA, int N_pop,
		 float *mut_par, float tt_ratio,
		 //
		 float *fit_evo, float *DG_evo, 
		 float *Enat_evo, int nsam, int step,
		 double *f_ave, double f_dev, 
		 double *E_ave, double E_dev,
		 double *DG_ave, double DG_dev,
		 float seq_entr, double seq_entr_dev,
		 struct load *mut_load, struct load *trans_load,
		 long num_syn_subst, long num_aa_subst,
		 long num_syn_mut, long num_aa_mut,
		 float *dN_dS, float *accept,
		 long **nuc_evo, int len_dna, int sample);

// Calculations
void Compute_nuc_mut(char *dna_seq, int len_dna,
		     short *aa_seq, int len_amm,
		     char **codon, char *coded_aa,
		     char *SSC_TYPE, float **exp1, float **exp2,
		     float *Lambda, FILE *file_mut);
void Regularize(float **f_reg_ia, float **n_msa_ia, float *f_aa,
		float w_max, float reg, int len_amm, int Naa);
float Normalize_exponent(float **SSC_mod, int L, char *model);
void Str_mut_wt(float **exp_MF, float **P_MF, short *aa_seq, int L, int Naa);
void Initialize_load(struct load *load);
void Sum_loads(struct load *load, struct load *load1);

int WildType_DDG(float **DG_mut, float **C_nat, int *i_sec, char *c_sec,
		 short *aa_seq, struct REM *E_wt);
void  Sum_DDG(float **exp_ia, float *norm, float **exp_pdb, float w,
	      int *ali_PDB, int len_amm, int Naa);
void Normalize_DDG(float **exp_ia, float *norm, int len_amm, int Naa,
		   float Lambda);
void Full_prot(float **C_nat_all, float **i_sec_count, float **aa_seq_count,
	       int L_ali, float *norm_c, float w,
	       float **C_nat, int *i_sec, short *aa_seq, int *ali_PDB);
int Set_domain(int *ini_dom, int *end_dom, int *ali_PDB, int L_ali);
void  Align_MSA_PDB(short **ali_seq, char **msa_seq, int *ali_PDB,
		    int N_seq, int L_ali, int end_dom);
int Align_str_MSA(int *ini_dom, int *end_dom, struct protein *prot,
		  int *ali_PDB, int L_ali);

void Make_prot(float **C_nat_all, struct res_short **res,
	       int **i_sec, char **c_sec, short **aa_seq,
	       float **i_sec_count, float **aa_seq_count,
	       int len_amm);

float Compute_stab_seq(float *DG,
		       float *N_stab_str,
		       float *DG_ave,
		       char *buffer,
		       struct REM E_REM,
		       struct protein prot,
		       float **C_nat, int *i_sec, char *c_sec,
		       short **ali_seq, char **name_seq,
		       float *seq_id, int N_seq);
int Write_DG(float **DG,
	     float **seq_id,
	     struct protein *prot,
	     int N_PDB,
	     float *N_stab_str,
	     char **name_seq,
	     int N_seq,
	     char *buffer,
	     char *name_ali);
float **Compute_t_dist(float *t_ave,
		       char **msa_seq, int *l_stab, int N_seq, int L_ali);
void  Get_P_glob(float *P_glob_a, float **P_exp, float *wi,
		int L, float *f_aa, int Naa);

int Optimize_distr(struct MF_results *opt_res, float **P_WT_ia,
		   float *P_glob_a, float **exp1, float **exp2, 
		   float **f_reg_ia, float *wi, int L, int Naa,
		   struct REM *E_wt, float **C_nat, int *i_sec, char *c_sec,
		   char *name_file, FILE *file_summ, char *label,
		   int print);

int Optimize_reg(struct MF_results *opt_res, float **P_WT_ia,
		 float *P_glob_a, float **exp1, float **exp2, 
		 float **f_reg_ia, float *wi, int L, int Naa);
//struct REM *E_wt, float **C_nat, int *i_sec, char *c_sec,
//		 char *name_file,FILE *file_summary,char *label,
//		 int print);
int Analytic_Lambda_KL(struct MF_results *opt_res, float **P_WT_ia,
		       float *P_glob_a, float **exp1, float **exp2,
		       float **f_reg_ia, float *wi, int L, int Naa);
int Analytic_Lambda_lik(struct MF_results *opt_res, float **P_WT_ia,
			float *P_glob_a, float **exp1, float **exp2,
			float **f_reg_ia, float *wi, int L, int Naa,
			char *label);
//float **n_msa_ia, struct REM *E_wt,
//			float **C_nat,int *i_sec,char *c_sec,
//			char *name_file, FILE *file_summ, char *label,
//			int repeat);
int Maximize_Lambda(struct MF_results *opt_res, float **P_WT_ia,
		    float *P_glob_a, float **exp1, float **exp2,
		    float **f_reg_ia, float *wi, int L, int Naa,
		    int repeat);

float Analytic_Lambda_lik_1(float *Cv, float *Lambda, float **P_WT_ia,
			  float **exp1, float **exp2, float *P_glob_a, int L,
			  float *wi, float *cn, float Reg);
void Compute_score(struct MF_results *res,
		   float **P_MF_ia, float *Lambda,
		   float **f_reg_ia, float *wi, int L, int Naa, int comp_all);
float Compute_lik(float **P_MF_ia, float **freq_ia, float *w, int L, int Naa);

void Inverse_corr_reg(float **Corr_reg_inv, float **Ave,
		      float **P_WT_ia, float ***phi, float *wi,
		      int L, float Reg, int nmod, float *Lambda);
float Compute_Cv_lik(float *Lambda, float *dlik_dL, float **Corr_reg_inv,
		     int nmod);
float Compute_lik_cn(float *Lambda, float *cn, float *Z, float *wi,
		     int nmod, int L);
float Compute_average(float *P_WT_a, float *phi);
float KL_symm(float *P, float *Q, int n);

double Optimize_reg_old(struct MF_results *opt_res, float **P_opt_ia,
		    float w_max, float *f_aa,
		    struct MF_results *res, float **P_ia,
		    float **exp1, float **exp2, float *P_glob_a,
		    float **f_reg_ia, float **n_msa_ia, float *wi,
		    struct REM *E_wt, float **C_nat, int *i_sec, char *c_sec,
		    char *name_file, FILE *file_summ, char *label,
		    int L, int Naa);
double Compute_score_reg(struct MF_results *SSC_res, float reg,
			 float w_max, float *f_aa, float **P_SSC_ia,
			 float **exp1, float **exp2, float *P_glob_a,
			 float **f_reg_ia, float **n_msa_ia, float *wi,
			 struct REM *E_wt, float **C_nat, int *i_sec,
			 char *c_sec, char *name_file, FILE *file_summ,
			 char *label, int L, int Naa);
double Compute_Cv(struct MF_results *SSC_res, float reg, float step,
		  float w_max, float *f_aa, float **P_SSC_ia,
		  float **exp1, float **exp2, float *P_glob_a,
		  float **f_reg_ia, float **n_msa_ia, float *wi,
		  struct REM *E_wt, float **C_nat, int *i_sec, char *c_sec,
		  char *name_file, FILE *file_summ, char *label,
		  int L, int Naa);
int Compute_P(float **P_WT_ia, float *Z, float *P_glob_a,
	      float *Lambda, float **exp1, float **exp2,
	      float *f_reg_aa, float *wi, short *aa_seq, 
	      float Pmin, int L, int Naa);
int Selection(float fitness, float fitness_old, int N_pop);
int Detailed_balance(float *p, int xend, int xini);
float Sequence_entropy(double **aa_distr, int L, int Naa);
float Sequence_entropy_mut(float *mut_par, char **codon, char *coded_aa);
void Compute_freq_codons(float *mut_par, float *freq_aa,
			 char **codon, char *coded_aa);
void Compute_load(double *Tload_sum, double *Tload_dev,
		  double *Mload_sum, double *Mload_dev, int *Nload,
		  struct REM *E_wt,
		  float **C_nat, int *i_sec, char *c_sec, short *aa_seq,
		  float fitness_wt, char *dna_seq, int len_dna, 
		  char **codon, char *coded_aa);
extern float Find_max_quad(float x1, float x2, float x3,
			   float y1, float y2, float y3,
			   float MIN, float MAX);
void Remut(float *P_glob_a, float *P_glob_bk,
	   struct MF_results *k_res, float **P_ia,
	   float **exp1, float **exp2, int is_str,
	   float **f_reg_ia, float **n_msa_ia, float *wi,
	   struct REM *E_wt, float **C_nat, int *i_sec, char *c_sec,
	   short *aa_seq, char *name_file, FILE *file_summ,
	   char *mname, FILE *out, int L, int Naa);

// Input parameters
// Mean-field model
// A: Input files defining the protein
static char dir_out[N_CHAR]; //file[N_CHAR], 
static char seq_name[N_CHAR];
static char *file_str_mut[4]; // File with structural mutation scores
// B: Thermodynamic parameters
int REM=2;   // Use 1st (1), 2nd (2) and 3rd (3) moment of misfolded energy 
float S_C, S_U;

// C: Selection model
int Samples=10; // Number of simulated samples of evolutionary trajectories
static int IT_MAX=0; // Number of iterations for each trajectory
static int NEUTRAL=1; // Neutral (1) versus continous (0) fitness
float DG_THR_COEFF=0.95; // Threshold in DeltaG for neutral selection
static int N_pop=100;  // Population size for continuous selection model
// C2: Mean-field model
int MEANFIELD=1;  // Compute site-specific mean-field distributions
int OPT_LAMBDA=0; // 1=Optimize Lambda by maximum likelihood
char MODEL[N_CHAR]="ALL"; // Type of mean-field model (ALL, NAT, DG)
float DG_OPT=-1;  // DG target of the optimization if MODEL="DG"
// D1: Mutation model, P_glob
int CpG=1;       // Enhance mutation rate at CpG dinucleotides
int GET_FREQ=3;  // 0= Get aa frequencies from input DNA model
                 // 1= Get aa frequencies from fitted DNA model
                 // 2= Get aa frequencies from fitted DNA model + observed
                 // 3= Get aa frequencies from observed ones

float tt_ratio=4, kCpG=4; //mut_par[MUTPAR], 


// D1: Mutation model, exchangeability
char MATRIX[40]="OPT"; //"OPT"; // Empirical exchangeability matrix
char MAT_DEF[40]="WAG";
char EXCHANGE='F'; // exchangeability model. M=MUT F=FLUX Q=RATE E=EXCH
float TWONUC=0; // Rate of two to one nucleotide substitutions if MUT
// E: Output
int FORMAT=1;   // PAML format
int PRINT_E=1; // Compute exchangeability matrix and rate for all sites?
int PRINT_GLOB=0; // Print global exchangeability matrices
int PRINT_ALL_EXCH=0; // Print exchangeability matrix for all sites
int PRINT_MUT=0;   // Print mutation exchangeability matrix

unsigned long iran;
// Data read from input files
static int len_amm=0, len_dna=0;
static int L_str_PDB=0; // length of conformation
static int L_ali=0;
static int L_noali=0;
static int L_nomod=0;
static int N_PDB_used;
static int N_PDB_read;
static int N_PDB;
static int count[4];
static float rate[4];
float DG_thr;

int N_disulf_seq=0;
int Len_disulf_seq=0;
float Log_len_disulf_seq=0;
int N_helix_seq=0;
int N_strand_seq=0;

// time
double nbtops = CLOCKS_PER_SEC;

int main(int argc, char **argv){

  
  /***********************
          INPUT
  ************************/
  // Input files
  // char Input_dir[N_CHAR];
  char single_pdb[100]="", chain[80]="", file_pdb_list[80]="", pdb_dir[100]="";
  char file_dna[80]="";
  
  char name_file[200];
  char FILE_CODE[N_CHAR];
  double t0=clock();

  // Genetic code
  // char *codon[64], coded_aa[64], name_code[200];
  
  /***********************
         SEQUENCES
  ************************/
  char *dna_seq=NULL;
  short *dna_seq_short=NULL;
  
  /***********************
         Wild Type
  ************************/
  double fitness_wt;
  
  /***********************
          DUMMY
  ************************/
  int i, j, a;


  // Default parameters of mutation model
  mut_par[0]=0.25; mut_par[1]=0.25; mut_par[2]=0.25; mut_par[3]=0.25;
  mut_par[4]=1; tt_ratio=1; mut_par[5]=1; kCpG=1; mut_par[6]=0; TWONUC=0;

  /******************** Input operations   ************************/
  sprintf(FILE_CODE, FILE_CODE_DEF);
  int N_str_mut=0; //strcpy(MATRIX, MAT_DEF);
  for(i=0; i<4; i++)file_str_mut[i]=malloc(N_CHAR*sizeof(char));
  Get_para(argc, argv, single_pdb, chain, file_pdb_list, pdb_dir,
	   file_ali, file_dna, tnm, tnm_mut_para,
	   file_str_mut, &N_str_mut, &IWT,
	   &MEANFIELD, &MF_COMP, &STR_COMP, &OPT_REG, &SCORE_CV, &REG_INI, 
	   FILE_STR, &TEMP, &sU1, &sC0, &sC1, &REM, &SEC_STR, 
	   &REMUT, &GET_FREQ, mut_par, &tt_ratio, &kCpG, &TWONUC,
	   &PRINT_E, &EXCHANGE, MATRIX, &FORMAT, &ALL_MUTS,
	   &IT_MAX, &Samples, &NEUTRAL, &N_pop,
	   &OPT_LAMBDA, &LAMBDA_MF, &DG_OPT, MODEL,
	   dir_out);
  if(IWT>=20)IWT=19;

  /*********************************************************
      Read multiple sequence alignment (amino acids) if any
  **********************************************************/
  char **name_seq=NULL, **msa_seq=NULL;
  char tmp_buff[1000];
  if(file_ali[0]){
    N_seq=Read_ali(&msa_seq, &L_ali, &name_seq, file_ali, tmp_buff);
    char *ptr=file_ali, *p1=ptr;
    while(*ptr!='\0'){if(*ptr=='/')p1=ptr+1; ptr++;}
    strcpy(name_ali, p1);
  }

  /*********************************************************
      Read PDB files
  **********************************************************/  
  N_PDB_read=0; N_PDB=0; int read_list=0;
  char **name_pdb=NULL; char **chain_pdb=NULL;
  if(file_pdb_list[0]){
    N_PDB_read=Read_pdb_files(&name_pdb, &chain_pdb, file_pdb_list);
    if(N_PDB_read)read_list=1;
    printf("Reading %d pdbs in file %s, path= %s\n",
	   N_PDB_read, file_pdb_list, pdb_dir);
  }
  if(N_PDB_read==0 && single_pdb[0]){
    N_PDB_read=Copy_pdb_name(&name_pdb, &chain_pdb, single_pdb, chain);
  }
  if(N_PDB_read==0){
    printf("ERROR, PDB file not given\n"); exit(8);
  }

  /*****************************************************
                        Log file
  **************************************************/

  char *file_pdb[N_PDB_read]; for(i=0; i<N_PDB_read; i++)file_pdb[i]=NULL;
  char buffer[N_PDB_read*3000];
  char prog_def[]=
    "Program Prot_evol, author Ugo Bastolla (CSIC-UAM) ubastolla@cbm.csic.es\nIt simulates protein evolution with selection on the thermodynamic stability of the folded state and the conservation of the native structure and it computes site-specific substitution matrices (frequencies and exchangeability) according to different structure- and stability-constrained models of protein evolution. It optionally performs simulations of protein evolution with stability constraints and predicts the stability change of all possible DNA mutations.\n\n";
  printf("%s", prog_def);
  strcpy(buffer, prog_def);
  strcat(buffer, tmp_buff);

  if(L_ali){
    sprintf(tmp_buff, "%d aligned sequences L_ali=%d read in %s\n",
	    N_seq, L_ali, file_ali);
    printf("%s", tmp_buff); strcat(buffer, tmp_buff);
  }
  double tt=t0;
  sprintf(tmp_buff, "Read MSA. Time= %.2lf sec.\n",(clock()-tt)/nbtops);
  tt=clock(); strcat(buffer, tmp_buff);

  sprintf(tmp_buff, "%d PDB files to read: ",N_PDB_read);
  if(read_list){strcat(tmp_buff, file_pdb_list);}
  else if(single_pdb[0]){strcat(tmp_buff, name_pdb[0]);}
  strcat(tmp_buff, "\n");
  printf("%s", tmp_buff); strcat(buffer, tmp_buff);
  printf("PDB_dir: %s\n", pdb_dir);

  /************************************************************
                  Read PDB files
  *************************************************************/
  struct protein prot_str[N_PDB_read];
  int L_PDB[N_PDB_read], npdb[N_PDB_read]; 
  // ali_PDB[i] is column of MSA aligned with PDB residue i, i if no MSA
  int *ali_PDB[N_PDB_read];
  float Seq_id_ave[N_PDB_read];
  float *seq_id[N_PDB_read];
  if(L_ali){
    for(i=0; i<N_PDB_read; i++){
      seq_id[i]=malloc(N_seq*sizeof(float));
      ali_PDB[i]=malloc(L_ali*sizeof(int));
    }
  }
  struct residue *res=NULL;
  struct res_short *res_short=NULL;

  //int *res_index=NULL;
  //len_amm=Get_pdb(&target, &res, &res_short, file_pdb, chain, &res_index, 1);

  int L_max=0, str_opt=-1;
  int n=0, n_not=0, n_noalign=0;
  for(int k=0; k<N_PDB_read; k++){
    if(res){free(res); free(res_short);}
    
    // Name of PDB file
    if(file_pdb[n]==NULL)file_pdb[n]=malloc(120*sizeof(char));
    strcpy(file_pdb[n], "");
    if(pdb_dir[0])sprintf(file_pdb[n], "%s/", pdb_dir);
    strcat(file_pdb[n], name_pdb[k]);
    if(name_pdb[k][4]=='\0')strcat(file_pdb[n], ".pdb");
    printf("Reading %s\n",  file_pdb[n]);
    
    // Read PDB file, compute contact matrix
    L_PDB[n]=Get_pdb(prot_str+n, &res, &res_short, file_pdb[n],
		     chain_pdb[k], &(prot_str[n].res_index), 1);
    if(L_PDB[n] <= 0){
      sprintf(tmp_buff, "WARNING, no protein found in file %s chain %s\n",
	      file_pdb[n], chain_pdb[k]);
      printf("%s", tmp_buff); strcat(buffer, tmp_buff);
      n_not++; continue;
    }
    if(L_ali){
      npdb[n]=Align_ali_PDB(ali_PDB[n], seq_id[n], Seq_id_ave+n,  //
			    msa_seq, name_seq, N_seq, L_ali,
			    prot_str[n].aa_seq, prot_str[n].length,
			    prot_str[n].L_PDB);
      if(npdb[n]<0){
	sprintf(tmp_buff, "WARNING, protein %s chain %s not aligned to MSA\n",
		name_pdb[k], chain_pdb[k]);
	printf("%s", tmp_buff); strcat(buffer, tmp_buff);
	n_noalign++; continue;
      }
      if(L_PDB[n]>L_max){L_max=L_PDB[n];}

      if(n!=k){
	chain_pdb[n]=chain_pdb[k];
	strcpy(name_pdb[n], name_pdb[k]);
      }
      n++;
    }
  }
  N_PDB=n;

  sprintf(tmp_buff, "%d PDB files read MSA and aligned. Time= %.2lf sec.\n",
	  N_PDB, (clock()-tt)/nbtops);
  tt=clock(); strcat(buffer, tmp_buff);

  /****************************************************************
        Set log file 
  ****************************************************************/
  char out_file[200]="", list[150]="", tmp[150];
  if(N_PDB >1 || file_pdb_list[0]){
    strcpy(list, file_pdb_list);
  }else{
    strcpy(list, name_pdb[0]);
    strcat(list, chain_pdb[0]);
  }
  Remove_path(list);
  Remove_extension(tmp, list);
  strcpy(list, tmp);
  sprintf(out_file, "%s.SSCPE.log", list);
  FILE *out=fopen(out_file, "w");
  fprintf(out, "%s", buffer);

  /****************************************************************
        Set length of model
  ****************************************************************/

  if(N_PDB==0){
    sprintf(buffer, "ERROR, none of the %d PDB files is valid: "
	   "%d do not exist and %d cannot be aligned\n",
	   N_PDB_read, n_not, n_noalign);
    printf("%s", buffer); fprintf(out, "%s", buffer);
    exit(8);
  }
  if(L_ali) {
    len_amm=L_ali;
    sprintf(buffer,"Length of model=%d equal to length of MSA, L_max= %d\n",
	    len_amm, L_max);
  }else if(N_PDB>1){
    sprintf(buffer, "ERROR, multiple PDB not allowed without MSA\n");
    printf("%s", buffer); fprintf(out, "%s", buffer);
    exit(8);
  }else{
    len_amm=L_max;
     sprintf(buffer, "Unique PDBs, length of model=%d L_ali= %d\n",
	     L_max, L_ali);
  }
  printf("%s", buffer); fprintf(out, "%s", buffer);
  sprintf(buffer, "\n");

  /****************************************************************
        Set preliminary global amino-acid model 
  ****************************************************************/

  float P_glob_a[Naa], *f_tmp;
  if(strncmp(MAT_DEF, "JTT", 3)==0){f_tmp=fJTT;}
  else if(strncmp(MAT_DEF, "WAG", 3)==0){f_tmp=fWAG;}
  else {f_tmp=f_LG;}
  for(a=0; a<Naa; a++){
    int ia=Code_AA(AA_PAML_1[a]);
    P_glob_a[ia]=f_tmp[a];
  }
  
  /**************************************************************
                Set SSCPE matrices
  **************************************************************/
  // 9 models: mut WT MF RMSD DE RMSDWT DEWT RMSD2WT2 DE2WT2
  if(N_str_mut==0)N_str_mut=2;
  int Nmodel=3*N_str_mut, kstr, kmod;
  if(MF_COMP){kstr=3;}else{kstr=2;} Nmodel+=kstr;
  
  int L_str_max=len_amm; if(L_max>L_str_max)L_str_max=L_max;
  struct MF_results MF_res[Nmodel], *k_res;
  char *mod_name[Nmodel];
  char *name_str[N_str_mut];
  float **P_pdb=Allocate_mat2_f(L_str_max, Naa), **P_ia=NULL;
  float **exp_pdb=Allocate_mat2_f(L_str_max, Naa), **exp_ia=NULL;
  float **P_MF_ia[Nmodel];
  float **exp_MF_ia[Nmodel];
  //float **exp2_MF_ia[N_str_mut];
  float **norm_mod=Allocate_mat2_f(Nmodel, len_amm);
  //float ***DDG1[Nmodel];
  //int allpair[Nmodel];
  int num_str[N_str_mut];
  for(kmod=0; kmod<Nmodel; kmod++){
    mod_name[kmod]=malloc(30*sizeof(char));
    P_MF_ia[kmod]=Allocate_mat2_f(len_amm, Naa);
    exp_MF_ia[kmod]=Allocate_mat2_f(len_amm, Naa);
    if(kmod<N_str_mut){
      name_str[kmod]=malloc(20*sizeof(char));
      num_str[kmod]=0;
    }
    k_res=MF_res+kmod;
    k_res->Lambda[0]=0; k_res->Lambda[1]=0;
    k_res->L=len_amm;
    k_res->KL_reg=0;  k_res->KL_mod=0; k_res->KL_mut=0;
    k_res->lik_reg=0; k_res->lik_flux=0;
    k_res->entropy=0;
    k_res->HB_opt=-1; k_res->flux_opt=-1;
    k_res->rate_opt=-1; 
  }
  score_sam=malloc(N_SAM_LAM*sizeof(float));
  Lambda_sam=malloc(N_SAM_LAM*sizeof(float *));
  for(i=0; i<N_SAM_LAM; i++){
    Lambda_sam[i]=malloc(2*sizeof(float));
  }

  /*******************************************
              Folding computations
  *************************************************/
  float *DG[N_PDB], DG_ave[N_PDB];
  float N_stab_seq[N_PDB], N_stab_str[N_seq];
  for(int i=0; i<N_seq; i++){N_stab_str[i]=0;}
  int N_PDB_unst=0;
  int n_DG=N_seq; if(n_DG==0){n_DG=1;}
  for(i=0; i<N_PDB; i++){
      DG[i]=malloc(n_DG*sizeof(float));
  }
  
  float norm_c=0;
  //short *aa_seq;
  int *i_sec;
  char *c_sec;
  float **C_nat=NULL;
  float **C_nat_all=Allocate_mat2_f(len_amm, len_amm);
  float **i_sec_count=Allocate_mat2_f(len_amm, N_SS);
  float **aa_seq_count=Allocate_mat2_f(len_amm, Naa);
  
  struct REM E_wt; E_wt.c1U1=NULL;

  /******************************************************
          Loop over all PDBs
  *******************************************************/
  N_stab_thr_2=f_stab_thr*N_seq;
  N_PDB_used=N_PDB;
  for(n=0; n<N_PDB; n++){

    // Set domain boundaries
    if(L_ali==0){
      ini_dom=0; end_dom=L_PDB[n]; L_dom=end_dom;
    }else{
      //L_dom=Align_str_MSA(&ini_dom,&end_dom,prot_str+n,ali_PDB[n],L_ali);
      L_dom=Set_domain(&ini_dom, &end_dom, ali_PDB[n], L_ali);
    }

    struct protein *prot=prot_str+n;
    printf("\nPDB %s %s L_seq=%d L_str: %d L_dom: %d\n",
	   name_pdb[n], chain_pdb[n], L_PDB[n], prot->L_PDB, L_dom);
    aa_seq=prot->aa_seq; // Based on PDB numbering
    i_sec=prot->i_sec;
    c_sec=prot->sec_str;
    C_nat=Fill_C_nat(end_dom, prot->contact, prot->num_chain);
    
    // REM computation of misfolded state
    E_wt.L=end_dom;
    S_C=sC0+L_dom*sC1; S_U=L_dom*sU1;
    Initialize_E_REM(&E_wt, end_dom, REM, TEMP, S_C, S_U, FILE_STR);
    E_wt.DeltaG=
      Compute_DG_overT_contfreq(&E_wt, aa_seq, C_nat, i_sec, c_sec, 0);
    sprintf(tmp_buff,
	    "PDB: %s%s DG/T= %.2f E_nat= %.4g "
	    "G_misf= %.4g G_unf= %.4g S_C= %.3g\n",
	    name_pdb[n], chain_pdb[n],
	    E_wt.DeltaG, E_wt.E_nat, E_wt.G_misf, -E_wt.S_U, E_wt.S_C);
    printf("%s", tmp_buff); strcat(buffer, tmp_buff);
    
    if(0){
      Test_contfreq(&E_wt, aa_seq, C_nat, i_sec, c_sec, prot_str[n].name);
      char nameout[300]; sprintf(nameout, "%s.DeltaG.dat", prot_str[i].name);
      double T0=Print_DG_contfreq(&E_wt, nameout);
      sprintf(tmp_buff,"Temperature with minimum DG: T=%.2f printed in %s\n",
	      T0, nameout);
      printf("%s", tmp_buff); strcat(buffer, tmp_buff);
    }
    
    if(E_wt.DeltaG > 0){
      sprintf(tmp_buff,"WARNING, unstable target structure!\n");
      printf("%s", tmp_buff); strcat(buffer, tmp_buff);
      N_PDB_unst++;
      goto end_pdb;
    }
    
    float w;
    if(L_ali){
      // ali_seq[n][i]: amino acid at column i in the aligned protein n
      // ali_PDB[n][i]: PDB residue at column i in the aligned protein n
      short *ali_seq[N_seq];
      Align_MSA_PDB(ali_seq, msa_seq, ali_PDB[n], N_seq, L_ali, end_dom);
      N_stab_seq[n]=Compute_stab_seq(DG[n], N_stab_str, DG_ave+n, buffer,
				     E_wt, prot_str[n], C_nat, i_sec, c_sec,
				     ali_seq, name_seq, seq_id[n], N_seq);
      for(int k=0; k<N_seq; k++){free(ali_seq[k]);}

      w=N_stab_seq[n];
      if(w<N_stab_thr && w<N_stab_thr_2){ //// NEW
	sprintf(tmp_buff,
		"%s will be eliminated from the MSA for further calculations\n",
		prot_str[n].name); strcat(buffer, tmp_buff);
	w=0; N_PDB_used--;
	goto end_pdb;
      }
      //w=w*w; //// NEW
      if(str_opt<0 ||
	 N_stab_seq[n] > N_stab_seq[str_opt] ||
	 (N_stab_seq[n]==N_stab_seq[str_opt] && DG_ave[n]<DG_ave[str_opt])){
	str_opt=n;
      }
      
    }else{
      for(i=0; i<len_amm; i++){ali_PDB[n][i]=i;}
      N_stab_seq[n]=1;
      N_stab_str[n]=1;
      Seq_id_ave[0]=1;
      seq_id[n][0]=1;
      DG_ave[0]=E_wt.DeltaG;
      DG[n][0]=E_wt.DeltaG;
      w=1;
    }
    
    if(N_PDB > 1 || L_PDB[n]!=len_amm){
      // The full aa_seq is made here
      Full_prot(C_nat_all, i_sec_count, aa_seq_count, L_ali, &norm_c, w,
		C_nat, i_sec, aa_seq, ali_PDB[n]); //
    }

    // Compute_DDG
    kmod=1;
    WildType_DDG(exp_pdb, C_nat, i_sec, c_sec, aa_seq, &E_wt);
    Sum_DDG(exp_MF_ia[kmod],norm_mod[kmod],exp_pdb,w,ali_PDB[n],len_amm,Naa);
    
    if(MF_COMP){
      kmod=2; k_res=MF_res+kmod;
      fprintf(out,"Mean-field model computed with fix Lambda= %.3f\n",
	      LAMBDA_MF);
      meanfield(P_pdb, k_res, P_glob_a, NULL, C_nat, i_sec, c_sec,
		LAMBDA_MF, E_wt, end_dom, Naa, 1); // Also misfolded state

      float l_Pglob[Naa];
      for(a=0; a<Naa; a++)l_Pglob[a]=log(P_glob_a[a]);
      for(int k=0; k<end_dom; k++){
	float *P_k=P_pdb[k];
	for(a=0; a<Naa; a++){
	  if(P_k[a]<PMIN){P_k[a]=PMIN;}
	  exp_pdb[k][a]=-(log(P_k[a])-l_Pglob[a]);
	  exp_pdb[k][a]/=LAMBDA_MF;
	}
      }
      Sum_DDG(exp_MF_ia[kmod], norm_mod[kmod], exp_pdb, w,
	      ali_PDB[n], len_amm, Naa);
    }
    kmod++;

    // compute str_mut
    char file_tnm[150]="";
    if(DIR_TNM[0]){sprintf(file_tnm, "%s/", DIR_TNM);}
    strcat(file_tnm, prot_str[n].name);
    sprintf(file_str_mut[0], "%s.mut_RMSD.dat", file_tnm);
    sprintf(file_str_mut[1], "%s.mut_DE.dat", file_tnm);
    FILE *file_in=fopen(file_str_mut[0], "r");
    if(file_in==NULL){
      if(STR_COMP){
	printf("Predict structural deformation with the program TNM (slow)\n");
	Run_tnm(file_pdb[n], chain_pdb[n], prot_str[n].name);
      }else{
	printf("WARNING, structural deformations do not exist "
	       "but I will not compute them because you required -nostr\n");
	goto end_pdb;
      }
    }else{
      printf("TNM results found in %s\n", file_str_mut[0]);
      fclose(file_in);
    }
    for(j=0; j<N_str_mut; j++){
      Read_name_str(name_str[j], file_str_mut[j]);
      int read=Read_str_mut(exp_pdb,file_str_mut[j],prot_str[n],end_dom,IWT);
      if(read==0){
	printf("WARNING, could not read model %s in file %s\n",
	       name_str[j], file_str_mut[j]); kmod++; continue;
      }
      num_str[j]++;
      Normalize_exponent(exp_pdb, end_dom, name_str[j]);
      if(IWT==0){
	Compute_P(P_pdb, NULL, P_glob_a,
		  Lambda_iwt, exp_pdb, NULL, NULL, NULL,
		  aa_seq, PMIN, end_dom, Naa);
	Str_mut_wt(exp_pdb, P_pdb, aa_seq, end_dom, Naa);
      }
      Sum_DDG(exp_MF_ia[kmod], norm_mod[kmod], exp_pdb, w,
	      ali_PDB[n], len_amm, Naa);
      kmod++;
    }
    printf("%d fitness models from PDB %s\n", kmod-1, prot_str[n].name);
    
  end_pdb:
    if(N_PDB > 1 || len_amm != L_max){
      Empty_matrix_f(C_nat, end_dom);
    }
  } // end PDBs
  
  /****************** Write DG for all PDBs and MSA seq *********/
  if(L_ali){
    Write_DG(DG, seq_id, prot_str, N_PDB, N_stab_str, name_seq, N_seq,
	     buffer, name_ali);
  }
  sprintf(tmp_buff, "End of PDB-based fitness predictions\n");
  printf("%s", tmp_buff); strcat(buffer, tmp_buff);

  sprintf(tmp_buff, "%d PDB used, %d eliminated if N_stab < min(%.3g,%.3g)\n",
	  N_PDB_used, N_PDB-N_PDB_used, N_stab_thr, N_stab_thr_2);
  printf("%s", tmp_buff); strcat(buffer, tmp_buff);
  
  sprintf(tmp_buff, "Folding computations. Time= %.2lf sec.\n",
	  (clock()-tt)/nbtops);
  tt=clock(); strcat(buffer, tmp_buff);
  fprintf(out, "%s", buffer);
  sprintf(buffer, "\n");
  
  /****************************************************
                  Set target structure
  ******************************************************/
  
  if(N_PDB_used==0){
    printf("ERROR, none of %d PDB satisfied the condition N_stab >= "
	   "min(%.1f,%.1f)\n", N_PDB, N_stab_thr, N_stab_thr_2);
    exit(8);
  }
  
  if(str_opt<0){
    if(L_ali){
      sprintf(buffer, "ERROR, no stable structure found\n");
      printf("%s", buffer); fprintf(out, "%s", buffer);
      exit(8);
    }else{
      str_opt=0;
    }
  }
  
  if(N_PDB > 1 || L_PDB[0]!=len_amm){
    printf("Making protein with %d a.a.\n",len_amm);
    fprintf(out, "Making protein with %d a.a.\n",len_amm);
    Make_prot(C_nat_all, &res_short, &i_sec, &c_sec, &aa_seq,
	      i_sec_count, aa_seq_count, len_amm);
    C_nat=C_nat_all;
    S_C=sC0+L_dom*sC1;
    S_U=L_dom*sU1;
    Initialize_E_REM(&E_wt, len_amm, REM, TEMP, S_C, S_U, FILE_STR);
  }
  
  /****************** Normalize fitness for MSA columns ************/
  float w_col[len_amm]; // It is 0 for columns without structure
  for(i=0; i<len_amm; i++)w_col[i]=0;
  for(kmod=1; kmod<(kstr+N_str_mut); kmod++){
    float *w=norm_mod[kmod], Lam=1.0;
    Normalize_DDG(exp_MF_ia[kmod], w, len_amm, Naa, Lam);
    for(i=0; i<len_amm; i++){if(w[i]>w_col[i])w_col[i]=w[i];}
  }
  L_nomod=0; for(i=0; i<len_amm; i++){if(w_col[i]==0)L_nomod++;}
  
  /* // Print target
     if(file_pdb_list[0]){
     char out_file[200]="";
     sprintf(out_file, "%s.target", list);
     FILE *file_out=fopen(out_file, "w");
     fprintf(file_out, "%s is most stabilizing PDB structure in %s\n",
     prot_str[str_opt].name, file_pdb_list);
     fclose(file_out);
     }*/
  
  target=prot_str[str_opt];
  L_str_PDB=target.L_PDB;
  N_stab_target=N_stab_seq[str_opt];
  DG_ave_target=DG_ave[str_opt];
  Seq_id_ave_target=Seq_id_ave[str_opt];
  int npdb_target=npdb[str_opt];
  if(target.name[4]=='\0'){
    char c=chain_pdb[str_opt][0];
    if(c && c!=' '){
      target.name[4]=c; target.name[5]='\0';
    }
  }
  
  sprintf(tmp_buff,
	  "# Structure that confers maximum stability to the MSA: %s (%d) "
	  "%.1f stable seq. Mean DG: %.2f\n",
	  target.name, npdb_target, N_stab_target, DG_ave_target);
  fprintf(out, "%s", tmp_buff); printf("%s", tmp_buff);
  
  if(N_PDB>1 || file_pdb_list[0]){
    char tmp[100]; Remove_extension(tmp, file_pdb_list);
    sprintf(name_file, "%s.SSCPE", tmp);
  }else{
    sprintf(name_file, "%s.SSCPE", target.name);
  }
  
  /******************** SSCPE computations **********************/
  sprintf(buffer, "\n\nRunning SSCPE computations\n");
  printf("%s", buffer); fprintf(out, "%s", buffer);

  /********************************************************
 	 Read or extract DNA sequence
  ***************************************************/	

  // Random numbers
  iran=randomgenerator();
  InitRandom( (RANDOMTYPE)iran);

  if(file_dna[0]!='\0'){
    // Read from file
    char inputseq[N_CHAR]; strcpy(inputseq, file_dna);
    if(Check_file(inputseq)==0){
      printf("WARNING, DNA sequence file %s not found\n", file_dna);
    }else if((aa_seq)&&(target.length>0)){
      int ns, *ini_seq, *len_seq;
      dna_seq=Read_sequence(&len_dna, &ns, &ini_seq, &len_seq, inputseq);
      if(dna_seq){
	if(Match_dna(dna_seq,ns,ini_seq,len_seq,target,codon,coded_aa)==0){
	  len_dna=0; free(dna_seq); dna_seq=NULL;
	}
      }
    }
    if(len_dna)
      fprintf(out, "DNA sequence of %d nuc. read in file %s\n",
	      len_dna, file_dna);
  }
  // If DNA sequence not given, get it from AA sequence
  int idna=1;
  float *num_dna=NULL;
  if(len_dna<=0){
    idna=0;
    printf("Drawing DNA sequence from AA sequence\n");
    fprintf(out,"Drawing DNA sequence from AA sequence\n");
    dna_seq=Extract_dna(&len_dna, len_amm, aa_seq, codon, coded_aa);
  }
  dna_seq_short=malloc(len_dna*sizeof(short));
  for(i=0; i<len_dna; i++)dna_seq_short[i]=Code_nuc(dna_seq[i]);
  if(idna)num_dna=Get_counts(dna_seq_short, len_dna, 4);

  /**************************  End INPUT  ****************************/

  /****************************************************************
                  Site-specific frequency and entropy
  *****************************************************************/
  /* n_msa[i][a] is the number of a.a. a in column i of MSA or PDB
     num_aa[a] = sum_i n_msa[i][a] 
     wi[i]= sum_a n_msa[i][a]
     f_reg_ia[i][a]=n_msa_ia[i][a]/max(w[i])+reg*f_aa[a]; */

  int l_stab[N_seq], nseq=0;
  for(int n=0; n<N_seq; n++){
    if(N_stab_str[n]<thr_stab){l_stab[n]=0;}else{l_stab[n]=1; nseq++;}
  }
  printf("%d sequences in the alignment could be stabilized\n", nseq);

  float wi[len_amm], w_max=0;
  float sum_msa=0, num_aa[Naa];
  for(a=0; a<Naa; a++){num_aa[a]=0;}
  float **n_msa_ia=Allocate_mat2_f(len_amm, Naa);
  float **flux_ab=NULL, flux_a[Naa];
  float sum_pair=0;
  if(L_ali){
    flux_iab=malloc(len_amm*sizeof(float **));
    flux_ab=Allocate_mat2_f(Naa, Naa);
    float **t_dist=Compute_t_dist(&t_ave, msa_seq, l_stab, N_seq, L_ali);
    for(i=0; i<len_amm; i++){
      flux_iab[i]=Allocate_mat2_f(Naa, Naa);
      int sum=0;
      for(int n=0; n<N_seq; n++){
	if(l_stab[n]==0)continue;
	if(msa_seq[n][i]=='-')continue;
	a=Code_AA(msa_seq[n][i]);
	if(a>=0 && a <Naa){n_msa_ia[i][a]++; sum++;}
	for(int m=0; m<n; m++){
	  if(l_stab[m]==0)continue;
	  int b=Code_AA(msa_seq[m][i]);
	  if(b>=0 && b <Naa){
	    flux_iab[i][a][b]++; //+=t_dist[n][m];
	  } 
	}
      }
      wi[i]=sum;
      if(wi[i]==0){printf("WARNING, column %d has no stable seq\n",i);}
      if(wi[i]>w_max)w_max=wi[i];
      sum_msa+=sum;
      for(a=0; a<Naa; a++){
	num_aa[a]+=n_msa_ia[i][a];
	for(int b=0; b<=a; b++){
	  if(b<a){
	    flux_iab[i][a][b]+=flux_iab[i][b][a];
	    flux_iab[i][b][a]=flux_iab[i][a][b];
	  }
	  flux_ab[a][b]+=flux_iab[i][a][b];
	}
      }
    } // end i
    Empty_matrix_f(t_dist, N_seq);
    int b; sum_pair=0;
    for(a=0; a<Naa; a++){flux_a[a]=0;}
    for(a=0; a<Naa; a++){
      sum_pair+=flux_ab[a][a];
      for(b=0; b<a; b++){
	flux_a[a]+=flux_ab[a][b];
	flux_a[b]+=flux_ab[a][b];
	sum_pair+=flux_ab[a][b];
      }
    }
    if(sum_pair){
      for(a=0; a<Naa; a++){
	flux_a[a]/=sum_pair;
	for(b=0; b<=a; b++){
	  flux_ab[a][b]/=sum_pair;
	  for(i=0; i<len_amm; i++){flux_iab[i][a][b]/=sum_pair;}
	}
      }
    }
  }

  if(sum_msa==0){
    printf("WARNING, alignment not present or contains no stable protein\n"
	   "Obtaining amino acid frequencies from PDB sequence\n");
    for(i=0; i<len_amm; i++){
      if((aa_seq[i]<0)||(aa_seq[i]>=20)){
	printf("WARNING, wrong aa code %d at site %d\n", aa_seq[i], i);
	if(aa_seq[i]>20)exit(8);
      }else{
	n_msa_ia[i][aa_seq[i]]=1;
	wi[i]=1;
      }
    }
    w_max=1;
  }
 
  sprintf(buffer, "L_ali= %d P_glob: ", L_ali);
  for(a=0; a<Naa; a++){
    sprintf(tmp_buff, " %.3f", P_glob_a[a]); strcat(buffer, tmp_buff);
  }
  strcat(buffer,"\n");
  fprintf(out, "%s", buffer);
  printf("%s", buffer);
  sprintf(buffer, "\n");

  float **f_reg_ia=Allocate_mat2_f(len_amm, Naa);
  // Sites-global amino acid frequencies
  float sum=0, f_aa[Naa];
  for(a=0; a<Naa; a++){
    f_aa[a]=0;
    for(i=0; i<len_amm; i++)f_aa[a]+=n_msa_ia[i][a];
    sum+=f_aa[a];
  }
  float norm_a=1.0;
  for(a=0; a<Naa; a++){
    f_aa[a]/=sum; if(f_aa[a]<PMIN){f_aa[a]=PMIN; norm_a+=PMIN;}
  }
  if(norm_a>1){for(a=0; a<Naa; a++)f_aa[a]/=norm_a;}

  Reg=REG_INI;
  // Regularize and normalize the frequencies
  //float reg=REG_DISTR/(sqrt(w_max)+REG_DISTR);
  char reg_out[100];
  float reg_d=REG_DISTR*(1.-Entr_ave/log(Naa));
  sprintf(reg_out,
	  "Regularization= %.2g(1-Mean(entropy)/ln(Naa))=%.3g, Entr=%.3g\n",
	  REG_DISTR, reg_d, Entr_ave);
  printf("%s", reg_out); fprintf(out,"%s", reg_out);
  Regularize(f_reg_ia, n_msa_ia, f_aa, w_max, reg_d, len_amm, Naa);

  // Normalize the weights wi so that their sum is 1
  float norm_w=0; for(i=0; i<len_amm; i++){norm_w+=wi[i];}
  for(i=0; i<len_amm; i++){wi[i]/=norm_w;}

  // Global observed frequency
  f_reg_aa=malloc(Naa*sizeof(float));
  lf_ia=Allocate_mat2_f(len_amm, Naa);
  sum=0;
  for(a=0; a<Naa; a++){
    float s=0;
    for(i=0; i<len_amm; i++){
      if(f_reg_ia[i][a]<PMIN){f_reg_ia[i][a]=PMIN;}
      lf_ia[i][a]=log(f_reg_ia[i][a]);
      s+=wi[i]*f_reg_ia[i][a];
    }
    f_reg_aa[a]=s; sum+=s;
  }
  if(sum>0){for(a=0; a<Naa; a++){f_reg_aa[a]/=sum;}}
  else{for(a=0; a<Naa; a++){f_reg_aa[a]=P_glob_a[a];}}

  // Compute entropy of each site from MSA (not regularized)
  Entr_ave=0;
  double norm=0; int NCmax=32;
  float Entr_ali[len_amm], h_ave[NCmax], num_nc[NCmax], norm_nc[NCmax];
  for(i=0; i<NCmax; i++){h_ave[i]=0; num_nc[i]=0; norm_nc[i]=0;}

  for(i=0; i<len_amm; i++){
    if(wi[i]==0){
      printf("WARNING, MSA position %d not aligned to stable sequences\n", i);
      Entr_ali[i]=-1; L_noali++; continue;
    }else if(L_ali){
      Entr_ali[i]=Entropy(n_msa_ia[i], Naa);
      Entr_ave+=wi[i]*Entr_ali[i]; norm+=wi[i];
    }else{
      Entr_ali[i]=0; norm++;
    }

    double h=0;
    for(a=0; a<Naa; a++){h+=f_reg_ia[i][a]*hydro[a];}
    int nc=0; for(j=0; j<len_amm; j++)nc+=C_nat[i][j];
    if(nc>=NCmax)nc=NCmax-1;
    h_ave[nc]+=wi[i]*h; norm_nc[nc]+=wi[i]; num_nc[nc]++;
  }

  if(L_noali){
    char o[80];
    sprintf(o,"%d MSA positions have not been aligned to stable sequences",
	    L_noali);
    printf("%s\n", o);
    fprintf(out, "%s: ", o);
    for(i=0; i<len_amm; i++){if(wi[i]==0)fprintf(out, " %d", i);}
    fprintf(out, "\n");
  }

  if(L_ali){
    Entr_ave/=norm;
    char name_ent[250]; sprintf(name_ent, "%s.entropy.dat", name_file);
    fprintf(out, "Average entropy of the MSA: %.2f\n", Entr_ave);
    fprintf(out, "Sequence entropy from alignment printed in %s\n", name_ent);
    FILE *file_ent=fopen(name_ent, "w");
    printf("Writing seq. entropy from alignment in %s\n", name_ent);
    for(i=0; i<len_amm; i++){
      char aa;
      if(aa_seq[i]>=0){aa=AMIN_CODE[aa_seq[i]];}else{aa='-';}
      fprintf(file_ent, "%d%c %.2f\n", i, aa, Entr_ali[i]);
    }
    fclose(file_ent);
  }
  // Entropy of regularized frequencies
  Entr_reg=0; norm=0;
  for(i=0; i<len_amm; i++){
    Entr_reg+=wi[i]*Entropy(f_reg_ia[i], Naa);
  }

  // Print hydrophobicity
  char name_h[200]; sprintf(name_h, "%s.hydrophobicity.dat", target.name);
  fprintf(out, "Mean hydroph. from MSA or PDB printed in %s\n", name_h);
  FILE *file_h=fopen(name_h, "w");
  printf("Writing mean hydro from alignment in %s\n", name_h);
  fprintf(file_h, "# ncont ave_hydro num\n");
  for(i=0; i<NCmax; i++){
    if(num_nc[i]==0)continue;
    fprintf(file_h, "%d %.3f %.0f\n", i, h_ave[i]/norm_nc[i], num_nc[i]);
  }
  fclose(file_h);

  /********************************************
   Mutation model
   ********************************************/
  // Mutation parameters
  float P_mut[Naa];
  if(kCpG<tt_ratio)kCpG=tt_ratio;
  if(TWONUC>1 || TWONUC<0){
    printf("WARNING, twonuc= %.1g not allowed\n",TWONUC); TWONUC=0;
  }
  mut_par[4]=kCpG; mut_par[5]=tt_ratio; mut_par[6]=TWONUC;
  if(TWONUC){NPAR=7;}else{NPAR=6;}
  float **Q_cod=Allocate_mat2_f(64, 64), P_cod[64];
  // Output: global amino acid frequencies

  /* Obtain global amino acid frequencies:
     if(GET_FREQ==0), from input mutation model
     if(GET_FREQ==1), from optimized mutation model
     if(GET_FREQ==2), combine optimized mutation model and a.a. frequencies
     if(GET_FREQ==3), from a.a. frequencies in MSA or protein sequence
     if (idna) Obtain nucleotide frequencies from input DNA sequence
     else get mutation parameters fitting amino acid frequencies
  */
  // Mutation model
  Get_mut_par(mut_par, P_mut, P_cod, Q_cod, GET_FREQ,
	      num_aa, len_amm, num_dna, target.name, 0);
  kCpG=mut_par[4]; tt_ratio=mut_par[5]; TWONUC=mut_par[6];
  printf("Mutation model ready\n");
  if(GET_FREQ!=3){
    for(a=0; a<Naa; a++)P_glob_a[a]=P_mut[a];
  }
  char Freq_mod[400]; int N_fit_par=0;
  sprintf(Freq_mod, "# Global a.a. frequencies");
  if(GET_FREQ==3){
    sprintf(tmp," obtained from %.1f stable MSA seq. (+F)\n",N_stab_target);
    strcat(Freq_mod, tmp);
    N_fit_par=19;
  }else{
    strcat(Freq_mod, " computed from mutation model");
    if(GET_FREQ){
      strcat(Freq_mod," optimized from the data\n");
      N_fit_par=3;
    }else{
      strcat(Freq_mod," given from input\n");
    }
    strcat(Freq_mod,"# Nucleotide frequencies: ");
    for(i=0; i<4; i++){
      sprintf(tmp, " %c %.3f", Nuc_code(i), mut_par[i]);
      strcat(Freq_mod, tmp);
    }
    sprintf(tmp," trasition-transversion ratio: %.3f", tt_ratio);
    strcat(Freq_mod, tmp);
    if(GET_FREQ && tt_ratio!=1)N_fit_par++;
    sprintf(tmp, " CpG ratio: %.3f", kCpG);
    strcat(Freq_mod, tmp);
    if(GET_FREQ && kCpG!=1)N_fit_par++;
    sprintf(tmp, "Double nucleotide mutations: %.3f\n", TWONUC);
    strcat(Freq_mod, tmp);
    if(GET_FREQ && TWONUC)N_fit_par++;
  }

  sprintf(tmp_buff, "Mutation model. Time= %.2lf sec.\n",(clock()-tt)/nbtops);
  tt=clock();
  fprintf(out, "%s", tmp_buff); printf("%s", tmp_buff);

  /**********************************************************
       Output
  ****************************************************************/
  fprintf(out,
	  "Representative protein structure: "
	  "%s %d residues %d structured %d dis.\n"
	  "%d MSA columns, %d PDB residues have no associated column\n",
	  target.name,target.length,L_str_PDB,target.length-L_str_PDB,
	  L_ali,L_noali);
  fprintf(out, "Stability model based on contact interactions, parameters:\n");
  fprintf(out, "Temperature: %.3g\n", TEMP);
  fprintf(out, "Conformational entropy unfolded: %.3f L\n", sU1);
  fprintf(out, "Conformational entropy misfolded: %.3f L + %.3f\n", sC1, sC0);
  fprintf(out, "Misfolding model REM= %d\n", REM);
  fprintf(out, "Coefficient of secondary structure free energy: %.2g\n\n",
	  SEC_STR);

  /****************************************************************
                   Mean-field computations
  *****************************************************************/
  int StaO=-1;
  int StrO=-1;
  int SSCO=-1; // Optimal stability and structure constrained model
  int SSC1=-1, SSC2=-1;
  float P_glob_bk[Naa];

  /*
    float ***DDG2[Nmodel];
    DDG2[kmod]=NULL; //allpair[kmod]=0;
    float **exp2_tmp=Allocate_mat2_f(len_amm, Naa);
  struct MF_results res_tmp;
  float **DfRMSD[len_amm], **DfDE[len_amm];
  for(i=0; i<len_amm; i++){
    DfRMSD[i]=Allocate_mat2_f(Naa, Naa);
    DfDE[i]=Allocate_mat2_f(Naa, Naa);
    }*******/

  /************************* Output files ***************************/
  //Output_name(name_file, dir_out, target.name, TEMP, sU1, sC1, MEANFIELD,
  //	      MODEL, LAMBDA_MF, OPT_LAMBDA, NEUTRAL, N_pop, mut_par);

  FILE *file_summ=Open_summary(name_file, target.length, E_wt, mut_par,
			       aa_seq, npdb_target);
  char dumm[200];
  strcpy(dumm, name_file);
  strcat(dumm, EXT_SUM);
  fprintf(out, "Summary results printed in file %s\n", dumm);

  if(MEANFIELD==0 && (ALL_MUTS==0 || idna==0))goto Simulate;

  /**************************************
         Empirical substitution matrix
   *****************************************/

  // Determine global frequencies
  strcpy(buffer, "");
  //printf("MATRIX: %s %.0f\n", MATRIX, sum_pair); exit(8);
  if(strncmp(MATRIX, "OPT", 3)==0){
    if(sum_pair){
      sprintf(tmp_buff, "Looking for optimal echangeability matrix\n");
      strcat(buffer, tmp_buff);
      char TMP_MAT[10];
      float Max_lik=0;
      printf("t_ave= %.2g\n", t_ave);
      for(int k=0; k<3; k++){
	if(k==0){strcpy(TMP_MAT, "JTT");}
	else if(k==1){strcpy(TMP_MAT, "LG");}
	else if(k==2){strcpy(TMP_MAT, "WAG");}
	float f_tmp[Naa],
	  **exch_tmp=Empirical_exchangeability(f_tmp, TMP_MAT, Naa);
	Normalize_exchange(exch_tmp, f_tmp, Naa, 1);
	// The amino acid index is the same as used by Prot_evol
	float lik_1, lik_log, lik=
	  Flux_exch(&lik_1, &lik_log, flux_ab, flux_a, exch_tmp, f_tmp, 20);
	lik=Lik_flux(lik, lik_1, lik_log);
	sprintf(tmp_buff, "%3s : flux lik= %.4g\n", TMP_MAT, lik);
	strcat(buffer, tmp_buff);
	fprintf(file_summ, "# %s", tmp_buff);
	if(k==0 || lik >Max_lik){
	  Max_lik=lik; sprintf(MATRIX, TMP_MAT);
	}
      }
      sprintf(tmp_buff, "The optimal exchangeability matrix is: %s %.4g\n",
	      MATRIX,Max_lik);
      fprintf(file_summ, "# %s", tmp_buff);
      strcat(buffer, tmp_buff);
    }else{
      if(strncmp(MATRIX, "OPT", 3)==0){
	strcat(buffer,"WARNING, no MSA for optimizing exchangeability matrix");
	sprintf(tmp_buff, ", using default %s\n", MAT_DEF);
	strcat(buffer, tmp_buff);
	fprintf(file_summ, "# %s", buffer);
	strcpy(MATRIX, MAT_DEF);
      }
    }
  }
  fprintf(out, "%s", buffer);
  printf("%s", buffer);
  sprintf(buffer, "\n");

  fprintf(file_summ, "# Adopted substitution matrix: %s\n", MATRIX);     
  if(strncmp(MATRIX, "JTT", 3)==0){f_tmp=fJTT;}
  else if(strncmp(MATRIX, "WAG", 3)==0){f_tmp=fWAG;}
  else {f_tmp=f_LG;}
  f_glob=malloc(Naa*sizeof(float));
  for(a=0; a<Naa; a++){
    int ia=Code_AA(AA_PAML_1[a]);
    P_glob_a[ia]=f_tmp[a];
    f_glob[ia]=f_tmp[a];
  }
  
  /**************************************************************

              Computation of SSCPE substitution matrices

  ****************************************************************/

  /********************* Print parameter *******************************/
  // Print model set-up
  fprintf(out,"Stability and structure constrained (SSC)");
  fprintf(out," site-specific substitution processes\n");

  int nstr=0;
  strcpy(mod_name[1], "WT");
  if(MF_COMP)strcpy(mod_name[2], "MF");
  for(i=0; i<N_str_mut; i++){
    j=kstr+i;
    strcpy(mod_name[j], name_str[i]);
    if(IWT){
      char dum[20]; sprintf(dum,".iwt%d", IWT);
      strcat(mod_name[j], dum);
    }
    if(num_str[i]==0)continue;
    nstr++;
  }
  if(nstr!=N_str_mut){
    printf("WARNING, %d str_mut files given as input but only %d found\n",
	   N_str_mut, nstr);
  }
  fprintf(out, "\nSelection Models:\n");
  fprintf(out, "Stability constrained: wild-type (WT) ");
  if(MF_COMP)fprintf(out, " Mean-field (MF)");
  fprintf(out, "\n");
  if(N_str_mut){
    fprintf(out, "Structure constrained:");
    j=kstr;
    for(i=0; i<N_str_mut; i++){
      if(num_str[i]==0)continue;
      fprintf(out, " %s", mod_name[j]); j++;
    }
    char moddef[200]="Wild-type RMSD and DE ", dum[80];
    if(IWT){
      sprintf(dum, "averaged over IWT=%d smallest ones",IWT);
    }else{
      sprintf(dum, "weighted mean with w=exp(-x)");
    }
    strcat(moddef, dum);
    fprintf(out, "\n%s\n",moddef);
    fprintf(file_summ,"# %s\n",moddef);
    fprintf(out, "Structure and stability constrained: ");
    for(i=0; i<N_str_mut; i++){
      if(num_str[i]==0)continue;
      for(int k=1; k<kstr; k++){
	fprintf(out, " %s%s", name_str[i], mod_name[k]);
	if(IWT)fprintf(out, ".iwt%d", IWT);
      }
    }
  }
  fprintf(out, "\n");
  printf("%d structural effect models read\n", nstr);

  if(REMUT && GET_FREQ)
    strcat(Freq_mod,"# Global frequencies optimized iteratively\n");
  sprintf(tmp, "# Number of fitted parameters: %d+1\n", N_fit_par);
  strcat(Freq_mod, tmp);
  fprintf(out, "%s", Freq_mod);
  fprintf(file_summ, "%s", Freq_mod);
  fprintf(file_summ, "# Empirical substitution model: %s\n", MATRIX);
  fprintf(file_summ,
	  "# Model Lambda0 Lambda1 lik(MSA) KL(mod,obs) "
	  "KL(obs,mod) -score KL(mod,mut) entropy(mod) DG Tf h\n");
  fflush(file_summ);
  fflush(out);

  ///////////////////////////////////////////////////////////////

  // Mutation model
  kmod=0; strcpy(mod_name[kmod], "mut");
  k_res=MF_res+kmod; P_ia=P_MF_ia[kmod];
  for(i=0; i<len_amm; i++)for(a=0; a<Naa; a++)P_ia[i][a]=f_reg_aa[a];
  //P_glob_a[a];
  printf("P_glob: "); for(a=0; a<Naa; a++)printf(" %.3f", f_reg_aa[a]);
  printf("\n");
  Test_distr(k_res, P_ia, f_reg_ia, wi,
	     C_nat, i_sec, c_sec, E_wt, len_amm, Naa);
  Print_results(*k_res, mod_name[kmod], file_summ);
  fflush(file_summ);

  // Wild-type model
  kmod=1; strcpy(mod_name[kmod], "WT");
  k_res=MF_res+kmod; P_ia=P_MF_ia[kmod]; exp_ia=exp_MF_ia[kmod];
  printf("Optimizing Lambda for model %s\n", mod_name[kmod]);
  if(Normalize_exponent(exp_ia, len_amm, "WT")<0)goto end_WT;
  k_res->Lambda[0]=Lambda_start[0]; k_res->Lambda[1]=0;
  Sample_Lambda=1;
  Optimize_distr(k_res, P_ia, P_glob_a, exp_ia, NULL,
		 f_reg_ia,wi,len_amm,Naa,&E_wt, C_nat, i_sec, c_sec,
		 name_file, file_summ, mod_name[kmod], 1);
  // Change mutational frequencies
  if(REMUT && GET_FREQ){
    Remut(P_glob_a, P_glob_bk, k_res, P_ia, exp_ia, NULL, 0, 
	  f_reg_ia, n_msa_ia, wi, &E_wt, C_nat,i_sec, c_sec, aa_seq,
	  name_file, file_summ, mod_name[kmod], out, len_amm, Naa);
  }

  // Print exchangeability matrices, compute lik_flux
  Print_profiles(P_ia, mod_name[kmod], k_res->DG, k_res->Lambda[0],
		 P_glob_a, aa_seq, len_amm, MATRIX, name_file,
		 PRINT_GLOB, w_col, out);
  Print_exchange(P_ia, mod_name[kmod], res_short, len_amm,
		 P_glob_a, P_cod, Q_cod, tt_ratio, TWONUC, name_file, FORMAT,
		 EXCHANGE, MATRIX, PRINT_E, PRINT_GLOB, PRINT_MUT,w_col,out,
		 k_res, flux_iab, t_ave, file_summ);
  //for(a=0; a<Naa; a++)P_glob_a[a]=P_glob_bk[a];

  // Consider evolutionary likelihood
  if(SCORE_EVO){k_res->score_reg += k_res->lik_flux;}
  if(StaO < 0 || k_res->score_reg > MF_res[StaO].score_reg){StaO=kmod;}
    
 end_WT:
  printf("%s model finished\n", mod_name[kmod]);
  sprintf(tmp_buff, "Model %s. Time= %.2lf sec.\n\n", mod_name[kmod],
	  (clock()-tt)/nbtops);
  tt=clock();
  fprintf(out, "%s", tmp_buff); printf("%s", tmp_buff);

  if(MEANFIELD==0 && ALL_MUTS==0 && idna)goto Compute_all_muts;

  // Mean-field model
  if(MF_COMP){
    // Fitness matrix derived from MF model MF_mod
    kmod++; strcpy(mod_name[kmod], "MF");
    k_res=MF_res+kmod; P_ia=P_MF_ia[kmod]; exp_ia=exp_MF_ia[kmod];

    if(Normalize_exponent(exp_ia, len_amm, "MF")<0)goto end_MF;
    k_res->Lambda[0]=Lambda_start[0]; k_res->Lambda[1]=0;
    Sample_Lambda=1;
    Optimize_distr(k_res, P_ia, P_glob_a, exp_ia, NULL,
		   f_reg_ia, wi, len_amm, Naa, &E_wt, C_nat,i_sec,c_sec,
		   name_file, file_summ, "MF ", 1);
    if(REMUT && GET_FREQ){
      Remut(P_glob_a, P_glob_bk, k_res, P_ia, exp_ia, NULL, 0,
	    f_reg_ia, n_msa_ia, wi,&E_wt,C_nat,i_sec,c_sec, aa_seq,
	    name_file, file_summ, mod_name[kmod], out, len_amm, Naa);
    }
    //end_MF:
    // Print exchangeability matrices, complute lik_flux
    Print_profiles(P_ia, mod_name[kmod], k_res->DG, k_res->Lambda[0],
		   P_glob_a, aa_seq, len_amm, MATRIX, name_file,
		   PRINT_GLOB, w_col, out);
    Print_exchange(P_ia, mod_name[kmod], res_short, len_amm,
		   P_glob_a, P_cod, Q_cod, tt_ratio, TWONUC, name_file, FORMAT,
		   EXCHANGE, MATRIX, PRINT_E, PRINT_GLOB, PRINT_MUT,w_col,out,
		   k_res, flux_iab, t_ave, file_summ);
    //for(a=0; a<Naa; a++)P_glob_a[a]=P_glob_bk[a];
    // Consider evolutionary likelihood
    if(SCORE_EVO){k_res->score_reg += k_res->lik_flux;}
    if(StaO < 0 || k_res->score_reg > MF_res[StaO].score_reg){StaO=kmod;}

    sprintf(tmp_buff, "Model %s. Time= %.2lf sec.\n\n", mod_name[kmod],
	    (clock()-tt)/nbtops);
    tt=clock();
    fprintf(out, "%s", tmp_buff); printf("%s", tmp_buff);
    //printf("Opt MF: %.3g %.3g\n", k_res->score_reg, k_res->Lambda[0]);
  } // end MF_comp

  // Structural effects model
  // Read structural effects of mutations
 end_MF:
  printf("MF computations finished\n");
  SSCO=StaO;

  // Structural effect model
  for(i=0; i<N_str_mut; i++){
    if(num_str[i]==0)continue;
    kmod++;
    k_res=MF_res+kmod; P_ia=P_MF_ia[kmod]; exp_ia=exp_MF_ia[kmod];
    printf("Optimizing Lambda for model %s\n", mod_name[kmod]);

    if(Normalize_exponent(exp_ia, len_amm, mod_name[kmod])<0)continue;
    k_res->Lambda[0]=Lambda_start[0]; k_res->Lambda[1]=0;
    Sample_Lambda=1;
    Optimize_distr(k_res, P_ia, P_glob_a, exp_ia, NULL,
		   f_reg_ia, wi, len_amm, Naa, &E_wt, C_nat,i_sec,c_sec,
		   name_file, file_summ, mod_name[kmod], 1);
    // Change mutational frequencies
    if(REMUT && GET_FREQ){
      Remut(P_glob_a, P_glob_bk, k_res, P_ia, exp_ia, NULL, 1,
	    f_reg_ia, n_msa_ia, wi, &E_wt,C_nat,i_sec,c_sec, aa_seq,
	    name_file, file_summ, mod_name[kmod], out, len_amm, Naa);
    }

    // Print exchangeability matrices, compute flux_lik
    Print_profiles(P_ia, mod_name[kmod], k_res->DG, k_res->Lambda[0],
		   P_glob_a, aa_seq, len_amm, MATRIX, name_file,
		   PRINT_GLOB, w_col, out);
    Print_exchange(P_ia, mod_name[kmod], res_short, len_amm,
		   P_glob_a, P_cod, Q_cod, tt_ratio, TWONUC, name_file, FORMAT,
		   EXCHANGE, MATRIX, PRINT_E, PRINT_GLOB, PRINT_MUT,w_col,out,
		   k_res, flux_iab, t_ave, file_summ);
    //for(a=0; a<Naa; a++)P_glob_a[a]=P_glob_bk[a];
    // Consider evolutionary likelihood
    if(SCORE_EVO){k_res->score_reg += k_res->lik_flux;}
    if(StrO < 0 || k_res->score_reg > MF_res[StrO].score_reg){StrO=kmod;}
    if(SSCO < 0 || k_res->score_reg > MF_res[SSCO].score_reg){SSCO=kmod;}

    sprintf(tmp_buff, "Model %s. Time= %.2lf sec.\n", mod_name[kmod],
	    (clock()-tt)/nbtops);
    tt=clock();
    fprintf(out, "%s", tmp_buff); printf("%s", tmp_buff);
    printf("Opt %s: %.3g %.3g\n\n",
	   mod_name[kmod], k_res->score_reg, k_res->Lambda[0]);
  }

  float **exp2_opt=NULL;
  
  // Structure and stability constrained models SSCPE
  SSC1=SSCO;
  Lambda_start[1]=Lambda_start[0];
  for(i=0; i<N_str_mut; i++){
    if(num_str[i]==0)continue;
    int jstr=kstr+i, jsta; //=jstab[i];
    if(isnan(MF_res[jstr].score_reg))continue;
    //exp2_MF_ia[i]=Allocate_mat2_f(len_amm, Naa);
    //float **exp2_ia=exp2_MF_ia[i];
    //Copy_P(exp2_ia, exp_MF_ia[jstr], len_amm, Naa);
    float **exp2_ia=exp_MF_ia[jstr];
    for(jsta=1; jsta<kstr; jsta++){
      kmod++;
      k_res=MF_res+kmod;
      P_ia=P_MF_ia[kmod]; exp_ia=exp_MF_ia[kmod];
      sprintf(mod_name[kmod], "%s%s", name_str[i], mod_name[jsta]);
      if(IWT){
	char tmp[20]; sprintf(tmp, ".iwt%d", IWT);
	strcat(mod_name[kmod],  tmp);
      }
      printf("Optimizing Lambda for %s model\n", mod_name[kmod]);
      k_res->Lambda[0]=MF_res[jsta].Lambda[0]*Sta_IWT;    
      k_res->Lambda[1]=MF_res[jstr].Lambda[0];
      Lambda_start[0]=k_res->Lambda[0];
      Lambda_start[1]=k_res->Lambda[1];
 
      Copy_P(exp_ia, exp_MF_ia[jsta], len_amm, Naa);
      Sample_Lambda=1;
      Optimize_distr(k_res, P_ia, P_glob_a, exp_ia, exp2_ia,
		     f_reg_ia, wi, len_amm, Naa, &E_wt, C_nat, i_sec, c_sec,
		     name_file, file_summ, mod_name[kmod], 1);
      if(REMUT && GET_FREQ){
	Remut(P_glob_a, P_glob_bk, k_res, P_ia, exp_ia, exp2_ia, 1,
	      f_reg_ia, n_msa_ia, wi, &E_wt,C_nat, i_sec, c_sec, aa_seq,
	      name_file, file_summ, mod_name[kmod], out, len_amm, Naa);
      }

      // Print exchangeability matrices, compute lik_flux
      Print_profiles(P_ia, mod_name[kmod], k_res->DG, k_res->Lambda[0],
		     P_glob_a, aa_seq, len_amm, MATRIX, name_file,
		     PRINT_GLOB, w_col, out);
      Print_exchange(P_ia, mod_name[kmod], res_short, len_amm, P_glob_a,
		     P_cod, Q_cod, tt_ratio, TWONUC, name_file,
		     FORMAT,EXCHANGE, MATRIX, PRINT_E, PRINT_GLOB,PRINT_MUT,
		     w_col,out, k_res, flux_iab, t_ave, file_summ);
      //for(a=0; a<Naa; a++)P_glob_a[a]=P_glob_bk[a];
      if(SCORE_EVO){k_res->score_reg += k_res->lik_flux;}
      if(SSCO < 0 || k_res->score_reg > MF_res[SSCO].score_reg){
	SSCO=kmod; SSC1=jsta; SSC2=jstr;
      }

      sprintf(tmp_buff, "Model %s. Time= %.2lf sec.\n\n", mod_name[kmod],
	      (clock()-tt)/nbtops);
      tt=clock();
      fprintf(out, "%s", tmp_buff); printf("%s", tmp_buff);
      printf("Opt: score= %.3g L1= %.3g L2= %.3g\n", k_res->score_reg,
	     k_res->Lambda[0], k_res->Lambda[1]);
      if(isnan(k_res->score_reg))continue;
    }// end stab mod
  } // end struct mod
 
  kmod++; 
  if(kmod!=Nmodel){
    fprintf(out, "WARNING, %d models expected but only %d could be computed\n",
	    Nmodel, kmod); Nmodel=kmod;
  }
  fprintf(out, "\nSite-specific amino-acid frequencies and "
	  "exchangeability matrices have been computed and printed\n");
  
  P_ia=P_MF_ia[SSCO]; exp_ia=exp_MF_ia[SSCO]; k_res=MF_res+SSCO;
  fprintf(file_summ, "# Best model: %s HB= %d flux= %d rate= %d\n",
	  mod_name[SSCO], k_res->HB_opt, k_res->flux_opt, k_res->rate_opt);
  if(0){
    Sample_Lambda=1;
    Optimize_distr(k_res, P_ia, P_glob_a, exp_ia, exp2_opt,
		   f_reg_ia, wi, len_amm, Naa, &E_wt,C_nat,i_sec,c_sec,
		   name_file, file_summ, mod_name[SSCO], 1);
    if(REMUT && GET_FREQ){
      Remut(P_glob_a, P_glob_bk, k_res, P_ia, exp_ia, exp2_opt, 1,
	    f_reg_ia, n_msa_ia, wi, &E_wt,C_nat, i_sec, c_sec,aa_seq,
	    name_file, file_summ, mod_name[SSCO], out, len_amm, Naa);
    }
  }
  
  // Print site specific substitution matrices, frequencies and
  // exchangeabilities
  /*fprintf(file_summ,
    "# The best stability constrained, structure constrained ");
    fprintf(file_summ,"and combined model are:\n");
    if(strncmp(STAB_TYPE, "MF", 2)==0){
    fprintf(out, "Mean-field model (MF)\n");
    }else{
    fprintf(out, "Wild-type model (WT). DDG of all possible mutations ");
    fprintf(out, "of wild-type sequence\n");
  }
  fprintf(file_summ, "%s\n", STAB_TYPE);
  Print_profiles(P_ia, STAB_TYPE,Stab_res->DG,Stab_res->Lambda[0],
		 P_glob_a, aa_seq, len_amm, MATRIX, name_file,
		 PRINT_GLOB, w_col, out);
  Print_exchange(P_ia, STAB_TYPE, res_short, len_amm,
		 P_glob_a, P_cod, Q_cod, tt_ratio, TWONUC, name_file, FORMAT,
		 EXCHANGE, MATRIX, PRINT_E, PRINT_GLOB, PRINT_MUT, w_col, out,
		 k_res, flux_iab, t_ave, file_summ);
		 for(a=0; a<Naa; a++)P_glob_a[a]=P_glob_bk[a];
  if(SCO>=0){
    fprintf(file_summ, "%s\n",STR_MUT_T[SCO]);
    fprintf(file_summ, "%s\n",SSC_TYPE_OPT);
    fprintf(out, "Best structure-constrained model (%s). %s ",
	    STR_MUT_T[SCO],STR_MUT_T[SCO]);
    fprintf(out, "of all possible mutations of wild-type sequence\n");
    Print_profiles(P_Str_ia[SCO], STR_MUT_T[SCO], Str_res[SCO].DG,
                   Str_res[SCO].Lambda[0], P_glob_a,  aa_seq, len_amm,
		   MATRIX, name_file, PRINT_GLOB, w_col, out);
    Print_exchange(P_Str_ia[SCO], STR_MUT_T[SCO], res_short, len_amm,
		   P_glob_a, P_cod, Q_cod, tt_ratio, TWONUC, name_file, FORMAT,
		   EXCHANGE, MATRIX, PRINT_E, PRINT_GLOB, PRINT_MUT,w_col,out,
		   k_res, flux_iab, t_ave, file_summ);
		   for(a=0; a<Naa; a++)P_glob_a[a]=P_glob_bk[a];
		   }*/ 
  // Optimal regularization, maximizing specific heat

  if(0 && OPT_REG && SCORE_KL && SSCO>=0){
    struct MF_results res_opt, *res=MF_res+SSCO;
    float **P_opt=Allocate_mat2_f(len_amm, Naa);
    for(int n=0; n<2; n++){
      if(res->Lambda[n]>0)Lambda_start[n]=res->Lambda[n];
    }
    Optimize_reg_old(&res_opt, P_opt, w_max, f_aa,
		     res, P_MF_ia[SSCO], exp_MF_ia[SSCO], exp2_opt, P_glob_a,
		     f_reg_ia, n_msa_ia, wi, &E_wt, C_nat,i_sec,c_sec,
		     name_file, file_summ, mod_name[SSCO], len_amm, Naa);
    
    fprintf(file_summ, "%s\n", mod_name[SSCO]);
    fprintf(out, "Structure and stability constrained model (%s). %s ",
	    mod_name[SSCO], mod_name[SSCO]);
    fprintf(out, "of all possible mutations of wild-type sequence\n");
    Print_profiles(P_opt, mod_name[SSCO], MF_res[SSCO].DG,
		   MF_res[SSCO].Lambda[0], P_glob_a, aa_seq, len_amm,
		   MATRIX, name_file, PRINT_GLOB, w_col, out);
    Print_exchange(P_opt, mod_name[SSCO], res_short, len_amm,
		   P_glob_a, P_cod, Q_cod, tt_ratio, TWONUC, name_file, FORMAT,
		   EXCHANGE, MATRIX, PRINT_E, PRINT_GLOB, PRINT_MUT,w_col,out,
		   k_res, flux_iab, t_ave, file_summ);
    Copy_P(P_MF_ia[SSCO], P_opt, len_amm, Naa);
    Empty_matrix_f(P_opt, len_amm);
    sprintf(tmp_buff, "Optimize REG. Time= %.2lf sec.\n",
	  (clock()-tt)/nbtops);
    tt=clock();
    fprintf(out, "%s", tmp_buff); printf("%s", tmp_buff);
    
  }
  fprintf(file_summ, "# SSCPE computations, Time= %.2lf sec.\n",
	  (clock()-t0)/nbtops);
  fclose(file_summ);
  printf("Output file name: %s %d pdb %d used\n", name_file,N_PDB,N_PDB_used);

 Compute_all_muts:
  if(ALL_MUTS && idna){
    char name_mut[400];
    sprintf(name_mut, "%s.all_nuc_mutations.dat", file_dna);
    FILE *file_mut=fopen(name_mut, "w");
    Compute_nuc_mut(dna_seq, len_dna, aa_seq, len_amm, codon, coded_aa,
		    mod_name[StaO],exp_MF_ia[StaO],NULL,MF_res[StaO].Lambda,
		    file_mut);
    if(StrO>=0)
      Compute_nuc_mut(dna_seq, len_dna, aa_seq, len_amm, codon, coded_aa,
		      mod_name[StrO],exp_MF_ia[StrO],NULL,MF_res[StrO].Lambda,
		      file_mut);
    if(SSCO>=0 && SSC2>=0)
      Compute_nuc_mut(dna_seq, len_dna, aa_seq, len_amm, codon, coded_aa,
		      mod_name[SSCO],exp_MF_ia[SSC1],exp_MF_ia[SSC2],
		      MF_res[SSCO].Lambda, file_mut);
  }

  /******************** End Mean-field  *********************/

  if(clean_REM){
    char cmd[80];
    sprintf(cmd, "rm -f E_loc*\n"); system(cmd);
    sprintf(cmd, "rm -f Local_interaction*\n"); system(cmd);
    sprintf(cmd, "rm -f REM.txt\n"); system(cmd);
  }

 Simulate:
  if(IT_MAX==0)return(0);
  
  /***********************************************************

                      Simulations of evolution
 
  *************************************************************/
  // S samples are simulated each for IT_MAX iterations
  fprintf(out, "Simulations: %d samples of %d substitutions each\n",
	  Samples, IT_MAX);
  char fit_def[100];
  if(NEUTRAL==0){
    sprintf(fit_def, "Fitness=1./(1+exp(DeltaG)), Population size= %d\n",N_pop);
  }else{
    sprintf(fit_def,"Neutral fitness=Theta(%.2f*DeltaG_PDB-DeltaG)\n",
	    DG_THR_COEFF);
  }
  fprintf(out, "%s\n", fit_def);
  fprintf(out, "Mutations are simulated at the DNA level ");
  fprintf(out, "using the genetic code %s\n", FILE_CODE);

  char name_simul[200]; sprintf(name_simul, "%s.Simul", target.name);
  FILE *file_ave =open_file(name_simul, EXT_AVE, aa_seq, 0, fit_def);
  FILE *file_stab=open_file(name_simul, EXT_OUT, aa_seq, 1, fit_def);

  // Accumulate samples
  double p_accept_all=0, p2_accept_all=0;
  double p_stab_incr_all=0, p2_stab_incr_all=0;

  double fit1_all=0, fit2_all=0, DG1_all=0, DG2_all=0,
    Enat1_all=0, Enat2_all=0, dNdS1_all=0, dNdS2_all=0,
    accept1_all=0, accept2_all=0, seqentr1_all=0, seqentr2_all=0;
  long it_all=0;
  struct load mut_load_all, trans_load_all;
  Initialize_load(&mut_load_all);
  Initialize_load(&trans_load_all);

  // Global variables valid for all independent samples
  // Number of iterations
  int it_print=10; if(it_print > IT_MAX*0.2)it_print= IT_MAX*0.2;
  if(it_print==0)it_print=1;
  int it_trans=it_print;
  int NPRINT= 1; if(it_print)NPRINT+=(IT_MAX/it_print);
  int MUTMAX= len_amm;  // Exhaustive search after MUTMAX unfixed mutations  

  // Amino acid distributions and entropy
  // Reference sequence
  short *aa_seq0=malloc(len_amm*sizeof(short));
  for(i=0; i<len_amm; i++)aa_seq0[i]=aa_seq[i];
  double **aa_distr=malloc(len_amm*sizeof(double *));
  double **aa_distr0=malloc(len_amm*sizeof(double *));
  double **aa_distr_all=malloc(len_amm*sizeof(double *));
  for(i=0; i<len_amm; i++){
     aa_distr[i]=malloc(Naa*sizeof(double));
     for(j=0; j<Naa; j++)aa_distr[i][j]=0;
     aa_distr0[i]=malloc(Naa*sizeof(double));
     for(j=0; j<Naa; j++)aa_distr0[i][j]=0;
     aa_distr_all[i]=malloc(Naa*sizeof(double));
     for(j=0; j<Naa; j++)aa_distr_all[i][j]=0;
  }
  float seq_entr_mut=Sequence_entropy_mut(mut_par, codon, coded_aa);

  // Nucleotide distributions at each codon position
  long **nuc_evo=malloc(len_dna*sizeof(long *));
  for(i=0; i<len_dna; i++){
    nuc_evo[i]=malloc(4*sizeof(long));
    for(j=0; j<4; j++)nuc_evo[i][j]=0;
  }

  // Averages
  float *fit_evo=malloc(NPRINT*sizeof(float));
  float *DG_evo=malloc(NPRINT*sizeof(float));
  float *Enat_evo=malloc(NPRINT*sizeof(float));

  // Fitness
  if(MEANFIELD)Divide_by_Pmut(P_MF_ia[StaO], P_glob_a, len_amm, Naa);
  int EXHAUSTIVE=0;
  if((MEANFIELD==0)&&(N_pop > len_amm*10))EXHAUSTIVE=1;
  int INI_EXH=EXHAUSTIVE;
  float Inverse_pop=1./N_pop;

  /***********************
         MUTANT
  ************************/
  double fitness_mut;
  char dna_new;
  int nuc_mut, nuc_new, res_mut=0, aa_new;
  
  // Thermodynamic computations
  struct REM E_mut; E_mut.c1U1=NULL;
  Initialize_E_REM(&E_mut, E_wt.L, E_wt.REM,
		   E_wt.T, E_wt.S_C, E_wt.S_U, FILE_STR);

  /***********************
          OUTPUT
  ************************/

  // FILE *file_dna=open_file(name_simul, EXT_DNA, aa_seq, 2, fit_def);


  FILE *file_msa=NULL; char name_MSA[300];
  if(PRINT_MSA==0 || IT_MAX==0){
    MAX_MSA=0;
  }else{
    sprintf(name_MSA, "%s%s", name_simul, EXT_MSA);
    file_msa=fopen(name_MSA, "w");
    Print_seq(file_msa, aa_seq, len_amm, &nseq_msa, E_wt.DeltaG);
    fprintf(out,"MSA from simulation printed in %s\n", name_MSA);
  }
  char name_div[300]; FILE *file_div=NULL;
  if(PRINT_TN){
    sprintf(name_div, "%s.TN_div.dat", name_simul);
    file_div=fopen(name_div, "a");
    fprintf(file_div, "# Protein %s L=%d\n", target.name, len_amm);
    fprintf(file_div, "#Tajima-Nei_divergence num_subst\n");
    fprintf(out,"Writing Tajima-Nei divergence versus time in %s\n",name_div);
  }
  strcpy(dumm, name_simul); strcat(dumm, EXT_FIN);
  fprintf(out,"Final results of evolutionary simulations printed in %s\n",dumm);
  strcpy(dumm, name_simul); strcat(dumm, EXT_SEL);
  fprintf(out, "Selection summary printed in %s\n", dumm);
  strcpy(dumm, name_simul); strcat(dumm, EXT_AVE);
  fprintf(out, "Average quantities obtained in simulations printed in %s\n",
	  dumm);
  strcpy(dumm, name_simul); strcat(dumm, EXT_OUT);
  fprintf(out, "Results for every accepted substitution printed in %s\n",dumm);


  /****************************************************
                 S samples of evolution
  *****************************************************/

  for(int is=0; is<Samples; is++){

    // Start the trajectory

    // Number of substitutions
    int aa_subst=0, synonymous=0;
    int it1=0;
    long it_subst=0; //int msa_subst=0;
    // Mutations per one aa substitution
    long tot_mut=0, naa_mut=0, nsyn_mut=0, nsyn_subst=0;

    // Count substitutions
    long num_syn_subst=0, num_aa_subst=0;
    long num_syn_mut=0, num_aa_mut=0;
    int num_stab_incr=0;

    // Average stability and fitness
    long it_sum=0;
    double E_nat_ave=0, f_ave=0, DG_ave=0;
    double E_nat_dev=0, f_dev=0, DG_dev=0;

    // Sequence entropy
    double seq_entr_sum=0, seq_entr_dev=0, entr_dev=0;
    float seq_entr;
    int n_seq_entr=0;

    for(i=0; i<len_amm; i++)aa_seq[i]=aa_seq0[i];
    for(i=0; i<len_amm; i++){
      for(j=0; j<Naa; j++)aa_distr[i][j]=0;
    }

    // Loads
    struct load mut_load, trans_load;
    Initialize_load(&mut_load);
    Initialize_load(&trans_load);

    // Wild-type stability
    E_wt.DeltaG=
      Compute_DG_overT_contfreq(&E_wt, aa_seq, C_nat, i_sec, c_sec, 0);    
    // Wild-type fitness
    if(NEUTRAL==0){
      fitness_wt=1./(1+exp(E_wt.DeltaG));
    }else{
      fitness_wt=1;
      DG_thr=E_wt.DeltaG*DG_THR_COEFF;
    }
    fprintf(file_stab, "#\n#\n# Evolutionary trajectory %d\n#\n#\n",is+1);
    fprintf(file_stab, "%3d  %c%c  %.3f  %.3f %7.4f",
	    0, AMIN_CODE[aa_seq[0]], AMIN_CODE[aa_seq[0]],
	    E_wt.E_nat, E_wt.DeltaG, fitness_wt);

    // Record fitness, energy and DeltaG
    int nsam=0;
    Record_evo(fit_evo, fitness_wt, DG_evo, E_wt.DeltaG,
	       Enat_evo, E_wt.E_nat, &nsam);
    
    // Mutation rates
    Ini_count(dna_seq, len_dna, count);
    Compute_rates(rate, mut_par, tt_ratio);

    // Simulations
    while(1){
    
      if(EXHAUSTIVE){
	// Extract next substitution (requires full exploration of sequence
	// neighbors, only convenient if N is large).
	Substitute_exhaustive(&res_mut, &aa_new, &nuc_mut, &nuc_new,
			      &mut_load, &trans_load, it_subst-it_trans,
			      &fitness_mut, &E_mut, &E_wt, 
			      &naa_mut, &nsyn_mut, &nsyn_subst,
			      NEUTRAL, DG_thr, N_pop, C_nat, i_sec, c_sec,
			      aa_seq, dna_seq, dna_seq_short, len_dna,
			      codon, coded_aa, tt_ratio, mut_par);
	tot_mut=naa_mut+nsyn_subst+nsyn_mut+1;
	dna_new=Nuc_code(nuc_new);
	if(it_subst > it_trans){
	  num_aa_mut += naa_mut;
	  num_syn_mut+= nsyn_mut;
	  num_syn_subst += nsyn_subst;
	}
	goto Update_aa;
      }

      tot_mut++;
      if(tot_mut>=MUTMAX){
	printf("WARNING, too many unfixed mutations > %d\n", MUTMAX);
	printf("T= %.2f MEANFIELD= %d NEUTRAL= %d", TEMP,MEANFIELD,NEUTRAL); 
	if(NEUTRAL==0)printf(" N= %d",N_pop);
	printf(" subst= %ld transient=%d\n", it_subst, it_trans);
	printf("DeltaG = %.3g ",E_wt.DeltaG);
	E_wt.DeltaG=
	  Compute_DG_overT_contfreq(&E_wt, aa_seq, C_nat, i_sec, c_sec, 0);
	printf("DeltaG recomputed = %.3g\n",E_wt.DeltaG);
	if(MEANFIELD==0){
	  printf("Performing exhaustive search of mutants\n");
	  EXHAUSTIVE=1; INI_EXH=1;
	  aa_subst=1; synonymous=0;
	  continue;
	}else{
	  exit(8);
	}
      }

      // Individual mutation
      synonymous=Mutate_seq(dna_seq, len_dna, codon, coded_aa,
			    aa_seq, len_amm, mut_par, tt_ratio, count, rate,
			    &nuc_mut, &dna_new, &res_mut, &aa_new);

      if(synonymous > 0){
	nsyn_mut++;
	if(it_subst >= it_trans)num_syn_mut++; 
	//if((MEANFIELD==0)&&(RandomFloating() > Inverse_pop))continue;
	if((RandomFloating() > Inverse_pop))continue;
	nsyn_subst++;
	if(it_subst >= it_trans)num_syn_subst++;
	nuc_new=Code_nuc(dna_new);
	goto Update_dna;

      }else if(synonymous==0){

	if(it_subst >= it_trans)num_aa_mut++;

	// Compute stability and fitness
	Copy_E_REM(&E_mut, &E_wt);
	E_mut.DeltaG=
	  Mutate_DG_overT_contfreq(&E_mut,aa_seq,C_nat,i_sec,c_sec,
				   res_mut,aa_new);
	if(0){
	  int aa_old=aa_seq[res_mut]; aa_seq[res_mut]=aa_new;
	  struct REM E_mut2; E_mut2.c1U1=NULL;
	  Initialize_E_REM(&E_mut2,E_wt.L, E_wt.REM,
			   E_wt.T, E_wt.S_C, E_wt.S_U, FILE_STR);
	  Copy_E_REM(&E_mut2, &E_wt);
	  float DeltaG=
	    Compute_DG_overT_contfreq(&E_mut2, aa_seq, C_nat, i_sec, c_sec, 0);
	  aa_seq[res_mut]=aa_old;
	  if(fabs(DeltaG-E_mut.DeltaG)>1){
	    printf("WARNING, DG= %.4g (full) %.4g (increment)\n",
		   DeltaG, E_mut.DeltaG);
	    printf("DEnat= %.3g DE1= %.3g DE2= %.3g "
		   "DE21= %.3g DE22= %.3g DE2s1= %.3g DE2s2= %.3g DE3= %.3g ",
		   E_mut2.E_nat-E_mut.E_nat, E_mut2.E1-E_mut.E1,
		   E_mut2.E2-E_mut.E2,
		   E_mut2.E2cont1-E_mut.E2cont1,
		   E_mut2.E2cont2-E_mut.E2cont2,
		   E_mut2.E2site1-E_mut.E2site1,
		   E_mut2.E2site2-E_mut.E2site2,
		   E_mut2.E3-E_mut.E3);
	    double sum=0; int i;
	    for(i=0; i<E_wt.L; i++)sum+=fabs(E_mut2.c1U1[i]-E_mut.c1U1[i]);
	    printf(" sum|DE1_i|= %.3g\n", sum);
	    E_mut.DeltaG=DeltaG;
	  }
	}
	if(NEUTRAL){
	  if(E_mut.DeltaG<DG_thr){fitness_mut=1;}else{fitness_mut=0;}
	  aa_subst=fitness_mut;
	}else{
	  fitness_mut=1./(1+exp(E_mut.DeltaG));
	  aa_subst=Selection(fitness_mut, fitness_wt, N_pop);
	}
	if(aa_subst<=0)continue;

      }else{ 
	if(it_subst >= it_trans)num_aa_mut++;
	continue;  // Stop codon
      }

    Update_aa:
      // A non-synonymous substitution has occurred
      
      // Print thermodynamic properties (also in transient)
      it_subst++;
      fprintf(file_stab, " %2ld %3ld\n", nsyn_subst, tot_mut);
      if(INI_EXH){
	fprintf(file_stab, "# Exhaustive search\n"); INI_EXH=0;
      }
      fprintf(file_stab, "%3d  %c%c  %.3f  %.3f %7.4f",
	      res_mut, AMIN_CODE[aa_seq[res_mut]], AMIN_CODE[aa_new],
	      E_mut.E_nat, E_mut.DeltaG, fitness_mut);
      fflush(file_stab);

      if((it_subst/it_print)*it_print == it_subst){
	Record_evo(fit_evo, fitness_wt, DG_evo, E_mut.DeltaG,
		   Enat_evo, E_mut.E_nat, &nsam);
      }
    
      // Write MSA 
      if(nseq_msa<MAX_MSA ){ //&& it_subst-msa_subst>5){
	//msa_subst=it_subst;
	Print_seq(file_msa, aa_seq, len_amm, &nseq_msa, E_mut.DeltaG);
      }
      
      // Statistics after the transient
      if(it_subst >= it_trans){

	// Average old amino acid sequence, weight=tot_mut
      
	// Update amino acid distribution
	for(i=0; i<len_amm; i++)aa_distr0[i][aa_seq[i]]+=tot_mut;
      
	// Averages
	E_nat_ave +=tot_mut*E_wt.E_nat;
	E_nat_dev+=tot_mut*E_wt.E_nat*E_wt.E_nat;
	DG_ave+=tot_mut*E_wt.DeltaG;
	DG_dev+=tot_mut*E_wt.DeltaG*E_wt.DeltaG;
	f_ave+=tot_mut*fitness_wt;
	f_dev+=tot_mut*fitness_wt*fitness_wt;
	it_sum+=tot_mut;
     
	num_aa_subst++; it1++;
	if(E_mut.DeltaG<E_wt.DeltaG)num_stab_incr++;

	// Print averages
	if(it1==it_print){
	  it1=0;
	  seq_entr=Sequence_entropy(aa_distr0, len_amm, Naa)-seq_entr_mut;
	  seq_entr_sum+=seq_entr; seq_entr_dev+=seq_entr*seq_entr;
	  n_seq_entr++;
	  for(i=0; i<len_amm; i++){
	    double *aa0=aa_distr0[i];
	    for(j=0; j<Naa; j++){aa_distr[i][j]+=aa0[j]; aa0[j]=0;}
	  }
	  Print_ave(file_ave, it_sum, num_aa_subst, N_pop,
		    f_ave, f_dev, E_nat_ave, E_nat_dev, DG_ave, DG_dev,
		    num_syn_subst, num_aa_subst, num_syn_mut, num_aa_mut,
		    seq_entr, mut_load, trans_load);
	  if(PRINT_TN){
	    aa_seq[res_mut]=aa_new;
	    Print_TN_div(aa_seq, aa_seq0, len_amm, it_subst, file_div);
	  }
	    
	} // end print
      }

      // Update_AA:
      aa_seq[res_mut]=aa_new;
      Copy_E_REM(&E_wt, &E_mut);
      fitness_wt=fitness_mut;
      tot_mut=0; nsyn_subst=0;

    Update_dna:
      if(aa_subst ||(synonymous > 0)){
	count[nuc_new]++;
	count[dna_seq_short[nuc_mut]]--;
	for(i=0; i<len_dna; i++)nuc_evo[i][dna_seq_short[i]]+=nsyn_mut;
	nsyn_mut=0;
	dna_seq_short[nuc_mut]=nuc_new;
	dna_seq[nuc_mut]=dna_new;
      }

      if(num_aa_subst == IT_MAX)break;
      
    }
    fprintf(file_stab, " %2ld %3ld\n", nsyn_subst, tot_mut);

    /*********************  End of simulation  **************************/

    // Loads
    /*Compute_load(&Tload_sum, &Tload_dev, &Mload_sum, &Mload_dev, &Nload,
      E_wt, C_nat, i_sec, c_sec, aa_seq, fitness_wt, dna_seq, len_dna,
      codon, coded_aa);*/

    // Entropy
    for(i=0; i<len_amm; i++){
      double *aa0=aa_distr0[i];
      for(j=0; j<Naa; j++){aa_distr[i][j]+=aa0[j]; aa0[j]=0;}
    }
    seq_entr=Sequence_entropy(aa_distr, len_amm, Naa)-seq_entr_mut;
    for(i=0; i<len_amm; i++){
      double *aa1=aa_distr[i];
      for(j=0; j<Naa; j++){aa_distr_all[i][j]+=aa1[j]; aa1[j]=0;}
    }

    if(n_seq_entr > 1){
      entr_dev=seq_entr_dev-seq_entr_sum*seq_entr_sum/n_seq_entr;
      entr_dev=sqrt(entr_dev)/n_seq_entr;
    }
    
    it_sum+=tot_mut;
    float dN_dS, accept;
    Print_final(name_simul, it_sum, TEMP, sU1, sC1, sC0, MEANFIELD,LAMBDA_MF,
		N_pop, mut_par, tt_ratio,
		fit_evo, DG_evo, Enat_evo, nsam, it_print,
		&f_ave, f_dev, &E_nat_ave, E_nat_dev, &DG_ave, DG_dev,
		seq_entr, entr_dev,
		&mut_load, &trans_load,
		num_syn_subst, num_aa_subst, num_syn_mut, num_aa_mut,
		&dN_dS, &accept, nuc_evo, len_dna, is);

    // Average samples
    it_all+=it_sum;
    fit1_all+=f_ave; fit2_all+=f_ave*f_ave;
    DG1_all+=DG_ave; DG2_all+=DG_ave*DG_ave;
    Enat1_all+=E_nat_ave; Enat2_all+=E_nat_ave*E_nat_ave;
    Sum_loads(&trans_load_all, &trans_load);
    Sum_loads(&mut_load_all, &mut_load);
    dNdS1_all+=dN_dS; dNdS2_all+=dN_dS*dN_dS;
    accept1_all+=accept; accept2_all+=accept*accept;
    seqentr1_all+=seq_entr; seqentr2_all+=seq_entr*seq_entr;

    {
      float p=num_aa_subst/(float)num_aa_mut;
      p_accept_all+=p; p2_accept_all+=p*p;
      p=num_stab_incr/(float)num_aa_subst;
      p_stab_incr_all+=p; p2_stab_incr_all+=p*p;
    }
  }

  /***************************************************
                     End of samples
  ***************************************************/

  if(MAX_MSA){
    printf("MSA printed in %s\n", name_MSA); fclose(file_msa);
  }
  
  {// Print summary on selection
    char name_sel[300];
    sprintf(name_sel, "%s%s", name_simul, EXT_SEL);
    FILE *file_sel=fopen(name_sel, "w");
    fprintf(out, "\nSimulating %d samples of evolutionary trajectories with %d substitution each and printing accepted mutations and fraction of substitutions with increased stability in file %s\n", Samples, IT_MAX, name_sel);
    float p=p_accept_all/Samples, var;
    fprintf(file_sel,
	    "# averages based on %d samples of %d substitutions each\n",
	    Samples, IT_MAX);
    if(Samples>1){var=(p2_accept_all-Samples*p*p)/(Samples-1);}
    else{var=p*p;}
    fprintf(file_sel, "Accepted mutations: %.3g s.e.= %.2g\n",
	    p, sqrt(var));
    p=p_stab_incr_all/Samples;
    if(Samples>1){var=(p2_stab_incr_all-Samples*p*p)/(Samples-1);}
    else{var=p*p;}
    fprintf(file_sel,
	    "Selected substitutions with increased stability: "
	    "%.3g s.e.= %.2g\n", p, sqrt(var));
    fclose(file_sel);
  }
  
  FILE *file_out=open_file(name_simul,"_samples.dat",aa_seq,0, fit_def);
  char name_sam[240]; sprintf(name_sam, "%s.samples.dat", name_simul);
  fprintf(out, "Printing averages for each sample (DeltaG, E_nat, fitness, "
	  "mutation load, translation load, dN/dS, acceptance ratio, "
	  "sequence entropy) in file %s\n", name_sam);
  
  fprintf(file_out, "# %d samples\n", Samples);
  fprintf(file_out, "# what mean s.e.\n");
  double x, dx; int S=Samples;
  Get_mean(&x, &dx, DG1_all, DG2_all, S, S); DG1_all=x;
  fprintf(file_out, "DeltaG/T\t%.5g\t%.2g\n", x, dx);
  fprintf(file_out, "DeltaG\t%.5g\t%.2g\n", x*TEMP, dx*TEMP);
  Get_mean(&x, &dx, Enat1_all, Enat2_all, S, S);
  fprintf(file_out, "E_nat\t%.5g\t%.2g\n", x, dx);
  Get_mean(&x, &dx, fit1_all, fit2_all, S, S);
  fprintf(file_out, "Fitness\t%.5g\t%.2g\n", x, dx);
  
  Get_mean(&x, &dx,mut_load_all.df_ave, mut_load_all.df_dev, S, S);
  fprintf(file_out, "Mut.load.fitness\t%.5g\t%.2g\n", x, dx);
  Get_mean(&x, &dx, mut_load_all.dG_ave, mut_load_all.dG_dev, S, S);
  fprintf(file_out, "Mut.load.DeltaG\t%.5g\t%.2g\n", x, dx);
  
  Get_mean(&x, &dx, trans_load_all.df_ave, trans_load_all.df_dev, S, S);
  fprintf(file_out, "Trans.load.fitness\t%.5g\t%.2g\n", x, dx);
  Get_mean(&x, &dx, trans_load_all.dG_ave, trans_load_all.dG_dev, S, S);
  fprintf(file_out, "Trans.load.DeltaG\t%.5g\t%.2g\n", x, dx);
  
  Get_mean(&x, &dx, dNdS1_all, dNdS2_all, Samples, Samples);
  fprintf(file_out, "dN/dS\t%.5g\t%.2g\n", x, dx);
  Get_mean(&x, &dx, accept1_all, accept2_all, Samples, Samples);
  fprintf(file_out, "acceptance\t%.5g\t%.2g\n", x, dx);
  Get_mean(&x, &dx, seqentr1_all, seqentr2_all, Samples, Samples);
  fprintf(file_out, "Seq.entropy\t%.5g\t%.2g\n", x, dx);
  fclose(file_out);
  Print_profile_evo(name_simul,aa_distr_all,aa_seq0,len_amm,DG1_all,it_all);

  strcpy(dumm, name_simul);
  strcat(dumm, ".AA_profiles_evo.txt");
  fprintf(out, "Printing evolutionary profiles in file %s\n", dumm);
  
  printf("Writing output files and explanations in %s\n", out_file);
  fclose(out);
  
  // Free memory
  if(aa_seq)free(aa_seq);
  if(aa_seq0)free(aa_seq0);
  if(aa_distr){for(i=0; i<len_amm; i++)free(aa_distr[i]); free(aa_distr);}
  if(aa_distr0){for(i=0;i<len_amm; i++)free(aa_distr0[i]);free(aa_distr0);}
  if(dna_seq)free(dna_seq);
  if(dna_seq_short)free(dna_seq_short);
  if(nuc_evo)free(nuc_evo);
  Empty_E_REM(&E_wt);
  Empty_E_REM(&E_mut);
  return(0);
}

void Print_ave(FILE *file_ave,
	       long it_sum, long n_subst,
	       int N_pop,
	       double f_sum, double f_dev,
	       double E_sum, double E_dev,
	       double DG_sum, double DG_dev,
	       long num_syn_subst, long num_aa_subst,
	       long num_syn_mut, long num_aa_mut,
	       float seq_entr,
	       struct load mut_load, struct load trans_load)
{
//   float t_indip=n_subst/Naa;
   float t_indep=1.0;

  if(ini_print==0){
    fprintf(file_ave, "# Fitness (sd)     Enat  (sd)    DG  (sd)");
    fprintf(file_ave, " seq.entropy(mut)-seq_entropy(sel)     ");
    //fprintf(file_ave, " Trans.load   (sd)  Mut.load  (sd)");
    fprintf(file_ave, "  nonsyn/syn accept  N_subst\n");
    ini_print=1;
  }

  double x, dx;
  Get_mean(&x, &dx, f_sum, f_dev, it_sum, t_indep);
  fprintf(file_ave, " %.4g\t%.2g\t", x, dx);
  Get_mean(&x, &dx, E_sum, E_dev, it_sum, t_indep);
  fprintf(file_ave, " %.4g\t%.2g\t", x, dx);
  Get_mean(&x, &dx, DG_sum, DG_dev, it_sum, t_indep);
  fprintf(file_ave, " %.4g\t%.2g\t", x, dx);
  fprintf(file_ave, " %.4f\t", -seq_entr);

  fprintf(file_ave, " %.3f\t",
	  (float)(num_aa_subst*num_syn_mut)/(float)(num_syn_subst*num_aa_mut));
  fprintf(file_ave, " %.4g\t", (float)(num_syn_subst+num_aa_subst)/it_sum);
  fprintf(file_ave, " %d\n", (int)num_aa_subst);
  fflush(file_ave);
}



void Get_mean(double *x, double *dx, float sum, float dev,
	      long it_sum, float t_indep)
{
  *x=sum/it_sum;
  double x2=dev-it_sum*(*x)*(*x);
  if(dev<=0){*dx=0; return;}
  *dx=sqrt(x2/(it_sum*t_indep));
}
    
void Print_final(char *name_file, long it_sum,
		 float TEMP, float sU1, float sC1, float sC0,
		 int MEANFIELD, float Lambda, int N_pop,
		 float *mut_par, float tt_ratio,
		 //
		 float *fit_evo, float *DG_evo, 
		 float *Enat_evo, int nsam, int step,
		 double *f_ave, double f_dev, 
		 double *E_ave, double E_dev,
		 double *DG_ave, double DG_dev,
		 float seq_entr, double seq_entr_dev,
		 struct load *mut_load, struct load *trans_load,
		 long num_syn_subst, long num_aa_subst,
		 long num_syn_mut, long num_aa_mut,
		 float *dN_dS, float *accept,
		 long **nuc_evo, int len_dna, int sample)
{
  char name_out[200];
  sprintf(name_out, "%s%s", name_file, EXT_FIN);
  printf("Writing %s\n", name_out);
  FILE *file_out=NULL;

  if(sample==0){
    file_out=fopen(name_out, "w");
    int L=len_dna/3;
    // Headers   
    fprintf(file_out, "# Input thermodynamic parameters:\n");
    fprintf(file_out, "L= %d residues\n", L);
    fprintf(file_out, "TEMP=\t%.2f\n", TEMP);
    fprintf(file_out, "Entropy unfolded=\t%.1f (%.2f *L)\n", sU1*L, sU1);
    fprintf(file_out, "Entropy compact=\t%.1f (%.2f + %.2f *L)\n",
	    sC0+sC1*L, sC0, sC1);
    if(SEC_STR)
      fprintf(file_out, "#Local interactions used, factor= %.3f\n", SEC_STR);
    fprintf(file_out, "# Input population parameters:\n");
    if(MEANFIELD){
      fprintf(file_out, "Meanfield model, Lambda=\t%.2f\n", LAMBDA_MF);
    }else{
      fprintf(file_out, "Population model, Npop=\t%d\n", N_pop);
    }
    fprintf(file_out, "# Input mutation parameters:\n");
    for(int i=0; i<4; i++)
      fprintf(file_out,"%c %.3f ",NUC_CODE[i],mut_par[i]);
    fprintf(file_out, " TT_ratio= %.2f\n", tt_ratio);
  }else{
    file_out=fopen(name_out, "a");
  }
  fprintf(file_out, "#\n#\n# sample %d ", sample+1);
  fprintf(file_out, "Attempted mutations= %ld\n", it_sum);
  fprintf(file_out, "# Results:\n");

  double x, dx, tau, r;  float t_indep=num_aa_subst/10;
  Get_mean(&x, &dx, *f_ave, f_dev, it_sum, t_indep);
  fprintf(file_out, "Fitness ave=\t%.4g %.4g\n", x, dx);
  *f_ave=x;

  x=Logistic_fit(&tau, &r, fit_evo, step, nsam);
  fprintf(file_out, "Fitness, logistic fit= %.4g tau= %.4g r=%.2g\n",x,tau,r);

  Get_mean(&x, &dx, *E_ave, E_dev, it_sum, t_indep);
  fprintf(file_out, "E_nat ave=\t%.4g %.2g\n", x, dx);
  *E_ave=x;

  x=Logistic_fit(&tau, &r, Enat_evo, step, nsam);
  fprintf(file_out, "E_nat, logistic fit= %.4g tau= %.4g r=%.2g\n",x,tau,r);

  Get_mean(&x, &dx, *DG_ave, DG_dev, it_sum, t_indep);
  fprintf(file_out, "DeltaG ave=\t%.4g %.2g\n", x, dx);
  *DG_ave=x;

  x=Logistic_fit(&tau, &r, DG_evo, step, nsam);
  fprintf(file_out, "DeltaG, logistic fit= %.4g tau= %.4g r=%.2g\n",x,tau,r);

  // Synonymous, mutation load
  *dN_dS=(float)(num_aa_subst*num_syn_mut)/(float)(num_syn_subst*num_aa_mut);
  fprintf(file_out, "dN/dS=\t%.3f\n", *dN_dS);

  *accept=(float)(num_syn_subst+num_aa_subst)/it_sum;
  fprintf(file_out, "Acceptance rate=\t%.4g\n", *accept);

  fprintf(file_out, "Seq_entropy(sel)-Seq_entropy(mut)=\t%7.4g %7.4g\n",
	  seq_entr, seq_entr_dev);

  int num=mut_load->num, indep=num/30;
  Get_mean(&x, &dx, mut_load->df_ave, mut_load->df_dev, num, indep);
  fprintf(file_out, "Mut.load.fitness=\t%.5g\t%.2g\n", x, dx);
  mut_load->df_ave=x;
  Get_mean(&x, &dx, mut_load->dG_ave, mut_load->dG_dev, num, indep);
  fprintf(file_out, "Mut.load.DeltaG=\t%.5g\t%.2g\n", x, dx);
  mut_load->dG_ave=x;

  Get_mean(&x, &dx, trans_load->df_ave, trans_load->df_dev, num, indep);
  fprintf(file_out, "Trans.load.fitness=\t%.5g\t%.2g\n", x, dx);
  trans_load->df_ave=x;
  Get_mean(&x, &dx, trans_load->dG_ave, trans_load->dG_dev, num, indep);
  fprintf(file_out, "Trans.load.DeltaG=\t%.5g\t%.2g\n", x, dx);
  trans_load->dG_ave=x;

  // Nucleotide content
  fprintf(file_out, "#\n# Final nucleotide frequencies (1, 2, 3)\n");
  int i, j, n; long nuc_count[4][3];
  for(n=0; n<4; n++){
    for(j=0; j<3; j++)nuc_count[n][j]=0;
  }
  j=0;
  for(i=0; i<len_dna; i++){
    for(n=0; n<4; n++)nuc_count[n][j]+=nuc_evo[i][n];
    j++; if(j==3)j=0;
  }
  double sum=0;
  for(n=0; n<4; n++)sum+=nuc_count[n][0];
  for(n=0; n<4; n++){
    fprintf(file_out, "%c ", NUC_CODE[n]);
    for(j=0; j<3; j++){
      fprintf(file_out," %.3f", nuc_count[n][j]/sum);
    }
    fprintf(file_out, "\n");
  }
  fclose(file_out);
}


int Read_ene_new(char *file_ene, float **interactions)
{
  FILE *file_in=fopen(file_ene,"r");
  int i, j, n=0; float ene;
  char string[200], aa1[8], aa2[8];

  if(file_in==NULL){
    printf("ERROR, energy parameter file %s not found\n", file_ene);
    exit(8);
  }
  while(fgets(string, sizeof(string), file_in)!=NULL){
    sscanf(string, "%s%s%f", aa1, aa2, &ene); n++;
    i=Code_AA(aa1[0]); j=Code_AA(aa2[0]);
    interactions[i][j]=ene; interactions[j][i]=ene;
  }
  return(n);
}


void Read_ene_par(char *file_ene, float **interactions)
{
  FILE *file_in=fopen(file_ene,"r");
  int i; char string[200];
  if(file_in==NULL){
    printf("ERROR, energy parameter file %s not found\n", file_ene);
    exit(8);
  }
  
  for(i=0; i<20; i ++){
    float *MAT=interactions[i];
    fgets(string, sizeof(string), file_in);
    sscanf(string, 
	   "%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f", 
	   &MAT[0],&MAT[1],&MAT[2],&MAT[3],&MAT[4],
	   &MAT[5],&MAT[6],&MAT[7],&MAT[8],&MAT[9],&MAT[10],
	   &MAT[11],&MAT[12],&MAT[13],&MAT[14],&MAT[15],
	   &MAT[16],&MAT[17],&MAT[18],&MAT[19]);
  }
  /*fgets(string, sizeof(string), matrix); // Disulfide bonds
    sscanf(string, "%f", &interactions[210]);
    interactions[210]-=interactions[label[17][17]];*/
  (void)fclose(file_in);
  return;
}

char *Read_sequence(int *len_dna, int *nseq, int **ini_seq, int **len_seq,
		    char *inputseq)
{  
  char *sequence, string[1000], *ptr, tmp[1000];
  FILE *file_in=fopen(inputseq, "r");
  int i=0;

  if(file_in==NULL){
    printf("WARNING, sequence file %s does not exist\n", inputseq);
    return(NULL);
  }

  *nseq=0;
  printf("Reading %s\n",inputseq);
  strcpy(seq_name,"");
  while(fgets(string, sizeof(string), file_in)!=NULL){
    if(string[0]=='>'){
      for(i=0; i<sizeof(string); i++)
	if(string[i]=='\n' || string[i]=='\r'|| string[i]=='\b')
	  string[i]='\0';
      sprintf(tmp, "%s", string+1); strcat(seq_name, tmp);
      (*nseq)++; continue;
    }
    ptr=string;
    while((ptr!=NULL)&&(*ptr!='\n')&&(*ptr!='\r')&&(*ptr!='\0')&&(*ptr!='\b')){
      if((*ptr=='a')||(*ptr=='A')||(*ptr=='t')||(*ptr=='T')||
	 (*ptr=='g')||(*ptr=='G')||(*ptr=='c')||(*ptr=='C')){
	(*len_dna)++;
      }else if((*ptr!=' ')){
	printf("Wrong character %d in DNA sequence %s: %c\n",
	       *len_dna, inputseq, *ptr);
	printf("%s\n", string); 
	exit(8);
      }
      ptr++;
    }
  }
  if(*nseq==0)*nseq=1;
  *ini_seq=malloc(*nseq*sizeof(int));
  *len_seq=malloc(*nseq*sizeof(int));

  if(seq_name[0]=='\0'){
    strcpy(seq_name, inputseq);
    for(i=0; i<sizeof(seq_name); i++){
      if((seq_name[i]=='_')||(seq_name[i]=='.')){
	seq_name[i]='\0'; break;
      }
    }
  }
  printf("DNA sequence %s\n", seq_name);
  printf("%d sequences of total length %d (%d a.a.)\n",
	 *nseq, *len_dna, *len_dna/3);
  fclose(file_in);

  if(*len_dna==0)return(NULL);

  // Reading
  sequence=(char *)malloc(*len_dna *sizeof(char));
  file_in=fopen(inputseq, "r"); i=0; *nseq=0;
  while(fgets(string, sizeof(string), file_in)!=NULL){
    if(string[0]=='>'){
      (*ini_seq)[*nseq]=i;
      if(*nseq)(*len_seq)[*nseq-1]=i;
      (*nseq)++;
      continue;
    }
    ptr=string;
    while((ptr!=NULL)&&(*ptr!='\n')&&(*ptr!='\r')&&(*ptr!='\0')&&(*ptr!='\b')){
      if((*ptr!=' ')){
	*ptr=Maiuscule(*ptr); sequence[i]=*ptr;
	//int i_nuc=Code_nuc(*ptr);
	i++;
      }
      ptr++;
    }
  }
  fclose(file_in);
  if(*nseq)(*len_seq)[*nseq-1]=i;
  for(i=0; i<*nseq; i++){
    (*len_seq)[i]-=(*ini_seq)[i];
    printf("Seq %d ini= %d len= %d\n", i, (*ini_seq)[i], (*len_seq)[i]);
  }
  
  return(sequence);
}

FILE *open_file(char *name_file, char *ext, short *seq, int lf, char *fit_def)
{
  FILE *file_out; char name[N_CHAR];
  int i, j; float gc=0;
  for(i=0; i<4; i++)if((NUC_CODE[i]=='G')||(NUC_CODE[i]=='C'))gc+=mut_par[i];
  sprintf(name, "%s%s", name_file, ext);
  file_out=fopen(name, "w");
  printf("Writing %s\n", name);
  fprintf(file_out, "# File %s, sequence %s,  PDB %s,  length=%d\n",
	  name_file, seq_name, target.name, len_amm);
  fprintf(file_out, "# fitness: %s\n", fit_def);
  fprintf(file_out, "# %d iterations, random seed: %ld\n", IT_MAX, iran);
  fprintf(file_out, "# Stationary frequencies: ");
  for(i=0; i<4; i++)fprintf(file_out, "%c %.3f ", NUC_CODE[i], mut_par[i]);
  fprintf(file_out, "\n# Transition-transversion ratio= %.2f\n", tt_ratio);
  fprintf(file_out, "# GC bias: %.3f\n", gc);
  fprintf(file_out, "# T= %.3f\n", TEMP);
  fprintf(file_out, "# Config. entropy per residue, unfolded: %.3f", sU1);
  fprintf(file_out, " misfolded: %.3f\n", sC1);
  if(SEC_STR)
    fprintf(file_out, "#Local interactions used, factor= %.3f\n", SEC_STR);
  fprintf(file_out, "# N_pop= %d\n#", N_pop);

  for(i=0; i<len_amm; i++){
    if(i==(i/60)*60)fprintf(file_out,"\n# ");
    fprintf(file_out,"%c", AMIN_CODE[seq[i]]);
  }
  fprintf(file_out,"\n");
  
  if(lf==1){
    fprintf(file_out,
	    "# site aa_old-aa_new  Enat DG fitness syn_subst attempts\n");
    //D(0,n)
  }else if(lf==2){
    fprintf(file_out, "# iter [GC]_3");
    for(i=0; i<3; i++){
      fprintf(file_out, " ");
      for(j=0; j<3; j++)fprintf(file_out, "%c_%d ", NUC_CODE[i], j+1);
    }
    fprintf(file_out, "\n\n");
    /*fprintf(file_out, "  Energy alpha\n");*/
  }
  
  fflush(file_out);
  return(file_out);
}

int Print_dna(char *seq, FILE *file_out, int iter)
{
  short i, j, i_nuc, nuc_count[4][3];
  float length=(float)len_amm, gc3=0;

  for(i=0; i<4; i++)
    for(j=0; j<3; j++)
      nuc_count[i][j]=0;

  j=0;
  for(i=0; i<len_dna; i++){
    i_nuc=Code_nuc(seq[i]);
    nuc_count[i_nuc][j]++; j++;
    if(j==3)j=0;
  }

  /* GC3 */
  for(i=0; i<4; i++)
    if((NUC_CODE[i]=='G')||(NUC_CODE[i]=='C'))gc3+=nuc_count[i][2];
  fprintf(file_out,"%3d  %.3f ", iter, (gc3)/length);
  for(i=0; i<3; i++){
    fprintf(file_out," ");
    for(j=0; j<3; j++){
      fprintf(file_out," %.3f", nuc_count[i][j]/length);
    }
  }
  fprintf(file_out, "\n");
  fflush(file_out);

  return(0);
}

int Get_name(char *name, char *name_seq, int N){
  int i;
  for(i=0; i<N; i++){
    if(strncmp(name_seq+i, ".dna", 4)==0){
       name[i]='\0'; 
       return(0);
    }
    name[i]=name_seq[i];
  }
  printf("ERROR in file name %s, extension .dna not found\n", name);
  exit(8);
}

void Output_name(char *file_name, char *dir_out, char *prot_name,
		 float TEMP, float sU1, float sC1, int MEANFIELD,
		 char *MODEL, float LAMBDA, int OPT_LAMBDA, 
		 int NEUTRAL, int N_pop, float *mut_par)
{
  char name[400];
  if(dir_out[0]!='\0'){sprintf(name, "%s/%s", dir_out, prot_name);}
  else{sprintf(name, "%s", prot_name);}
  if(MEANFIELD){
    if(OPT_LAMBDA){
      sprintf(file_name, "%s.MF.LAMBDA_OPT.%s", name, MODEL);
    }else{
      sprintf(file_name, "%s.MF.LAMBDA%.2f.%s", name, LAMBDA, MODEL);
    }
    //sprintf(file_name, "%s.REM%d.T%.2f.SU%.3f", file_name, REM, TEMP, sU1);
  }else if (NEUTRAL){
    sprintf(file_name, "%s.NEUTRAL.T%.2f.SU1%.2f.SC1%.2f",
	    name, TEMP, sU1, sC1);
  }else{
    sprintf(file_name, "%s.N%d.T%.2f.SU1%.2f.SC1%.2f",
	    name, N_pop, TEMP, sU1, sC1);
  }
  //sprintf(file_name, "%s.GC%.2f", file_name,
  //	  mut_par[Code_nuc('G')]+mut_par[Code_nuc('C')]);
  //if(CpG)sprintf(file_name, "%s.CpG%.0f", file_name, mut_par[4]);
}
 
float Sequence_entropy(double **aa_distr, int L, int Naa){
  int i, j; float S_sum=0, S, p; double norm=0;

  for(j=0; j<Naa; j++)norm+=aa_distr[0][j];
  for(i=0; i<L; i++){
    S=0;
    for(j=0; j<Naa; j++){
      p=aa_distr[i][j]/norm; if(p)S-=p*log(p);
    }
    S_sum+=S;
  }
  return(S_sum/L);
}

void Compute_freq_codons(float *mut_par, float *freq_aa,
			 char **codon, char *coded_aa)
{
  int i, i_aa; float w;
  for(i_aa=0; i_aa<Naa; i_aa++)freq_aa[i_aa]=0;
  for(i=0; i<64; i++){
    w=mut_par[Code_nuc(codon[i][0])];
    w*=mut_par[Code_nuc(codon[i][1])];
    w*=mut_par[Code_nuc(codon[i][2])];
    if(coded_aa[i]=='*')continue;
    i_aa=Code_AA(coded_aa[i]);
    if(i_aa<0)continue;
    freq_aa[i_aa]+=w;
  }
}

float Sequence_entropy_mut(float *mut_par, char **codon, char *coded_aa)
{
  float freq_aa[Naa]; int i; float norm=0, S=0, p;

  // Calculating amino acid distribution under mutation alone
  Compute_freq_codons(mut_par, freq_aa, codon, coded_aa);

  // Compute entropy
  for(i=0; i<Naa; i++)norm+=freq_aa[i];
  for(i=0; i<Naa; i++){
    if(freq_aa[i]){p=freq_aa[i]/norm; S-=p*log(p);}
  }
  return(S);
}

void Compute_load(double *Tload_sum, double *Tload_dev,
		  double *Mload_sum, double *Mload_dev, int *Nload,
		  struct REM *E_wt,
		  float **C_nat, int *i_sec, char *c_sec, short *aa_seq,
		  float fitness_wt, char *dna_seq, int len_dna, 
		  char **codon, char *coded_aa)
{
  double translation_load=0, mutation_load=0;

  // Folding thermodynamics
  float fitness, DeltaG;
  struct REM E_mut; E_mut.c1U1=NULL;
  Initialize_E_REM(&E_mut, E_wt->L,E_wt->REM,
		   E_wt->T,E_wt->S_C,E_wt->S_U, FILE_STR);

  /*
    Load=Sum_j P(nat->Seq_j)[finess(nat)-fitness(Seq_j)]
    P(nat->Seq_j) is non zero only if sequence j is one base mutation from
    the native sequence.
    Translation load: P(nat->Seq_j)=1
    Mutation load: P(nat->Seq_j)=Mutation probability,
    i.e. f(new base) times 1 if transversion, times tt_ratio if transition    
   */

  // Mutations
  int pos=-1, res_mut=0, aa_old=-1;
  char *codon_nat=dna_seq, codon_mut[3];
  int sum_trans=0, sum_mut=0;

  for(int nuc_mut=0; nuc_mut<len_dna; nuc_mut++){
    
    pos++;
    if(pos==3){
      pos=0; res_mut++; codon_nat=dna_seq+res_mut*3;
    }
    for(int j=0; j<3; j++)codon_mut[j]=codon_nat[j];
    int nuc_wt=Code_nuc(dna_seq[nuc_mut]);

    for(int base=0; base<4; base++){

      // Mutate amino acid
      if(base==nuc_wt)continue;
      codon_mut[pos]=NUC_CODE[base];
      int aa_new=Coded_aa(codon_mut, codon, coded_aa);
      if(aa_new==aa_old)continue;         // Synonymous
      if(aa_new<0){fitness=0; goto load;} // Stop codon

      // Folding thermodynamics
      Copy_E_REM(&E_mut, E_wt);
      DeltaG=
	Mutate_DG_overT_contfreq(&E_mut, aa_seq, C_nat, i_sec, c_sec,
				 res_mut, aa_new);
      fitness=1./(1+exp(DeltaG));
      
    load:
      // Translation load
      sum_trans++;
      translation_load+=fitness;

      // Mutation load
      float p=mut_par[base];
      if(base==Transition(NUC_CODE[nuc_wt]))p*=tt_ratio;
      sum_mut+=p;
      mutation_load+=fitness*p;
    }
  }
  Empty_E_REM(&E_mut);

  translation_load = fitness_wt-(translation_load)/sum_trans;
  mutation_load = fitness_wt-(mutation_load)/sum_mut;

  *Tload_sum += translation_load;
  *Tload_dev += translation_load*translation_load;
  *Mload_sum += mutation_load;
  *Mload_dev += mutation_load*mutation_load;
  (*Nload)++;
}

unsigned long randomgenerator(void){
     
     unsigned long tm;
     time_t seconds;
     
     time(&seconds);
     srand((unsigned)(seconds % 65536));
     do   /* waiting time equal 1 second */
       tm= clock();
     while (tm/CLOCKS_PER_SEC < (unsigned)(1));
     
     return((unsigned long) (rand()));

}

int Selection(float fitness, float fitness_old, int N_pop)
{
  // Moran's process:
  /* P_fix = (1-exp(-Df))/(1-exp(-N*Df)) */
  if(fitness<=0)return(0);
  //if((int)N_pop==1)return(1);
  double f_ratio= fitness_old / fitness;
  double x= (1.-f_ratio)/(1.-pow(f_ratio, N_pop));
  if(VBT){
    float rand=RandomFloating();
    printf("es= %.2f x= %.2g rand= %.2g\n", f_ratio, x, rand);
    if(rand < x)return(1);
  }else{
    if(RandomFloating() < x)return(1);
  }
  return(0);
}

int Detailed_balance(float *p, int xend, int xini)
{
  if(p[xend]>=p[xini])return(1);
  if(RandomFloating() < p[xend]/p[xini])return(1);
  return(0);
}

float *Get_counts(short *seq, int L, int Naa)
{
  float *num=malloc(Naa*sizeof(float));
  int i; for(i=0; i<Naa; i++)num[i]=0;
  for(i=0; i<L; i++){
    if(seq[i]>=Naa){
      printf("ERROR, wrong symbol %d at site %d\n", seq[i], i);
      exit(8);
    }else if(seq[i]>=0){
      num[seq[i]]++;
    }
  }
  //for(i=0; i<Naa; i++)num[i]/=L;
  return(num);
 }

void Record_evo(float *fit_evo, double fitness,
		float *DG_evo,  double DG,
		float *Enat_evo, double E_nat, int *nsam)
{
  fit_evo[*nsam]=fitness;
  DG_evo[*nsam]=DG;
  Enat_evo[*nsam]=E_nat;
  (*nsam)++;
}

int WildType_DDG(float **DG_mut, float **C_nat, int *i_sec, char *c_sec,
		 short *aa_seq, struct REM *E_wt)
{
  float DG;
  struct REM E_mut; E_mut.c1U1=NULL;
  Initialize_E_REM(&E_mut, E_wt->L, E_wt->REM,
		   E_wt->T, E_wt->S_C, E_wt->S_U, FILE_STR);
  E_wt->DeltaG=
    Compute_DG_overT_contfreq(E_wt, aa_seq, C_nat, i_sec, c_sec, 0);
  printf("Wild type mutants. DG_wt=%.2f /L= %.4f REM=%d\n",
	 E_wt->DeltaG, E_wt->DeltaG/E_wt->L, E_wt->REM);

  int L=0;
  for(int res_mut=ini_dom; res_mut<end_dom; res_mut++){
    float *DDG=DG_mut[res_mut];
    if(aa_seq[res_mut]<0){
      for(int aa=0; aa<Naa; aa++)DDG[aa]=0;
      continue;
    }

    //printf("%d  ", res_mut);
    L++;
    for(int aa=0; aa<Naa; aa++){
      if(aa==aa_seq[res_mut]){
	DDG[aa]=0;
      }else{
	Copy_E_REM(&E_mut, E_wt);
	DG=Mutate_DG_overT_contfreq(&E_mut, aa_seq, C_nat,
				   i_sec, c_sec, res_mut, aa);
	DDG[aa]=DG-E_wt->DeltaG;
      }
    }
  }
  printf("%d mutants computed\n", L*19);
  Empty_E_REM(&E_mut);
  return(0);
}

int Optimize_distr(struct MF_results *opt_res, float **P_WT_ia,
		   float *P_glob_a, float **exp1, float **exp2, 
		   float **f_reg_ia, float *wi, int L, int Naa,
		   struct REM *E_wt, float **C_nat, int *i_sec, char *c_sec,
		   char *name_file,FILE *file_summary,char *label, int print)
{
  opt_res->score_reg=-1000000;

  printf("## 1: Model %s L= %d\n", label, L);
  printf("#Setting Lambda by optimizing ");
  if(SCORE_KL){printf("KL");}else{printf("lik");}
  if(LAMBDA_ANALYTIC){printf(" by vanishing derivative\n");}
  else{printf(" numerically\n");}
  printf("#score Lambda\n");

  // New:
  exp_str=NULL;
  if(RECOMPUTE_IWT && IWT==0){
    if(exp2){exp_str=exp2;}
    else if(strcmp(label,"RMSD")==0 || strcmp(label,"DE")==0){exp_str=exp1;}
  }

  if(OPT_REG){
    SCORE_EVO_TMP=0;
    Optimize_reg(opt_res, P_WT_ia, P_glob_a, exp1, exp2, f_reg_ia, wi, L, Naa);
    fprintf(file_summary, "# Model %s Reg= %.3g\n", label, Reg);
  }else if(LAMBDA_ANALYTIC==0){
    SCORE_EVO_TMP=SCORE_EVO;
    Maximize_Lambda(opt_res,P_WT_ia,P_glob_a,exp1, exp2, f_reg_ia, wi, L,
		    Naa, repeat);
  }else{
    Compute_P(P_WT_ia, NULL, P_glob_a, opt_res->Lambda, exp1, exp2,
	      f_reg_aa, wi, aa_seq, PMIN, L,Naa);
    if(SCORE_KL){
      Analytic_Lambda_KL(opt_res, P_WT_ia, P_glob_a,
			 exp1,exp2,f_reg_ia,wi,L,Naa);
    }else{
      Analytic_Lambda_lik(opt_res, P_WT_ia, P_glob_a,
			  exp1, exp2, f_reg_ia, wi, L, Naa, label);
    }
  }

  if(isnan(opt_res->Lambda[0])|| (exp2 && isnan(opt_res->Lambda[1]))){
    printf("ERROR in Optimize_distr, Lambda is nan\n");
  }

  int zero=Compute_P(P_WT_ia,NULL,P_glob_a,opt_res->Lambda,exp1,exp2,
		     f_reg_aa,wi,aa_seq,PMIN,L,Naa);

  // f_reg_aa
  opt_res->L=L;
  if(print){
    Test_distr(opt_res, P_WT_ia, f_reg_ia, wi,
	       C_nat, i_sec, c_sec, *E_wt, L, Naa);
    double KLD=0;
    for(int i=0; i<L; i++)KLD+=wi[i]*KL(P_WT_ia[i], P_glob_a, Naa);
    opt_res->KL_mut=KLD;
    Print_results(*opt_res, label, file_summary);
    if(OPT_REG==0)
      printf("## Model %s Score= %.3f Lam1=%.3f Lam2=%.3f "
	     "Entropy= %.3f DG= %.3g L= %d\n",
	     label, opt_res->score_reg, opt_res->Lambda[0],
	     opt_res->Lambda[1], opt_res->entropy, opt_res->DG, opt_res->L);
  }else{
    Compute_score(opt_res, P_WT_ia, opt_res->Lambda,f_reg_ia, wi, L, Naa, 1);
  }

  return(zero);
}

void  Get_P_glob(float *P_glob_a, float **P_exp, float *wi,
		int L, float *f_aa, int Naa)
{
  int i, a, wrong=0; 
  float P_sel[Naa]; for(a=0; a<Naa; a++){P_sel[a]=0;}
  float u=1./Naa;

  int ns=L; 
  for(i=0; i<L; i++){
    float *P_i=P_exp[i], Z=0;
    for(a=0; a<Naa; a++){Z+=P_i[a];}
    if(wi[i]==0){ns--; continue;}
    if(isnan(Z) || Z<=0){wrong++; continue;}
    Z/=wi[i];
    for(a=0; a<Naa; a++){P_i[a]/=Z; P_sel[a]+=P_i[a];}
  }
  if(wrong){printf("WARNING, %d sites have Z nan or <=0\n", wrong);}
  float sum=0; for(a=0; a<Naa; a++){sum+=P_sel[a];}
  if(sum){
    for(a=0; a<Naa; a++){P_sel[a]/=sum; if(P_sel[a]<PMIN)P_sel[a]=PMIN;}
  }else{
    for(a=0; a<Naa; a++){P_sel[a]=u;}
  }
  sum=0;
  for(a=0; a<Naa; a++){
    P_glob_a[a]=f_aa[a]/P_sel[a];
    sum+=P_glob_a[a];
  }
  for(a=0; a<Naa; a++){
    P_glob_a[a]/=sum;
    if(P_glob_a[a]<PMIN)P_glob_a[a]=PMIN;
  }

  int iter=0;
  float P_tmp[Naa], Pmin=0.001, d_min=100;
  double P_glob_tmp[Naa];
  for(a=0; a<Naa; a++){P_glob_tmp[a]=P_glob_a[a];}

  for(iter=0; iter<it_max_mut; iter++){
    for(a=0; a<Naa; a++){P_sel[a]=0;} wrong=0;
    for(i=0; i<L; i++){
      if(wi[i]==0)continue;
      double Z=0; float *P_e=P_exp[i];
      for(a=0; a<Naa; a++){
	P_tmp[a]=P_glob_tmp[a]*P_e[a];
	Z+=P_tmp[a];
      }
      if(isnan(Z) || Z<=0){wrong++; continue;}
      Z/=wi[i];
      for(a=0; a<Naa; a++){P_tmp[a]/=Z; P_sel[a]+=P_tmp[a];}
    }
    sum=0; for(a=0; a<Naa; a++){sum+=P_sel[a];}
    if(sum){
      for(a=0; a<Naa; a++){P_sel[a]/=sum; if(P_sel[a]<PMIN)P_sel[a]=PMIN;}
    }else{
      break;
    }
    sum=0;
    double d=0;
    for(a=0; a<Naa; a++){
      d+=fabs(f_reg_aa[a]-P_sel[a]);
    }
    if(iter==0 || d<d_min){
      d_min=d;
      for(a=0; a<Naa; a++){P_glob_a[a]=P_glob_tmp[a];}
    }else{
      break;
    }
    // Update
    sum=0;
    for(a=0; a<Naa; a++){
      P_glob_tmp[a]*=(f_reg_aa[a]/P_sel[a]);
      sum+=P_glob_tmp[a];
    }
    for(a=0; a<Naa; a++){
      P_glob_tmp[a]/=sum; if(P_glob_tmp[a]<Pmin){P_glob_tmp[a]=Pmin;}
    }
  }
  
}

int Analytic_Lambda_lik(struct MF_results *opt_res, float **P_WT_ia,
			float *P_glob_a, float **exp1, float **exp2, 
			float **f_reg_ia, float *wi, int L, int Naa,
			char *label)
//float **n_msa_ia, 
//			struct REM *E_wt, float **C_nat, int *i_sec,char *c_sec,
//			char *name_file, FILE *file_summary, char *label,
//			int repeat)
{
  /* New routine 17/06/2022
     The selection parameter L is obtained by maximizing log-likelihood lik
     plus a regularization potential with parameter R.
     We consider two selection parameters L1 L2 that multiply the exponents
     exp1 exp2: P_ia=P_glob(a)*exp(L1*exp1_ia+L2*exp2_ia)/Z_i
     We set c1=sum_i fmsa_ia*exp1_ia independent of L1, L2
     d(lik(L)-RL^2)/dL=0 => 
     L1=sum_i wi*sum_a<exp1_ia>(L1,L2)/(2R+c1) same for L2
     The starting value of R is set to 0.05*max(c1,c2)
     The regularization parameter R is set by maximiing the specific heat
     Cv=-d(lik)/dR=-sum_k d(lik)/dL_k dL_k/dR
     d(lik)/dR= -2 M^(-1)L (vector notation)
     M_kl= sum_ia n_i(<(expk)(expl)>-<expk><expl>)+2R*delta_kl
     <y_ia>=sum_ia P_ia(L1,L2)*y_ia
   */

  float Lambda[2];
  int nmod=1; if(exp2){nmod=2; Lambda[1]=1;}else{Lambda[1]=0;}
  printf("# Optimizing Lambda analytically for %s model\n", label);

  // Initialize likelihood constant
  if(ini_lik){
    ini_lik=0; lik_const=0;
    for(int a=0; a<Naa; a++){
      double w_msa_a=0; //for(int i=0; i<L; i++)w_msa_a+=n_msa_ia[i][a];
      for(int i=0; i<L; i++)w_msa_a+=wi[i]*f_reg_ia[i][a];
      lik_const+=w_msa_a*log(P_glob_a[a]);
    }
    //norm_w=0; for(int i=0; i<L; i++)norm_w+=wi[i];
    //lik_const/=norm_w;
  }
  //float wn[L]; for(int i=0; i<L; i++)wn[i]=wi[i]/norm_w;

  // Compute correlations
  float cn[2]; cn[0]=0; cn[1]=0;
  for(int k=0; k<nmod; k++){
    float **phi;
    if(k==0){phi=exp1;}else{phi=exp2;}
    double sum=0;
    for(int i=0; i<L; i++){
      for(int a=0; a<Naa; a++){
	//sum+=n_msa_ia[i][a]*phi[i][a];
	sum+=wi[i]*f_reg_ia[i][a]*phi[i][a];
      }
    }
    cn[k]=sum; // /norm_w;
  }

  //Reg=cn[0];  if(Reg<cn[1])Reg=cn[1]; Reg*=REG_INI;
  //Reg=REG_INI;

  // Initialize likelihood
  float Cv, lik=
    Analytic_Lambda_lik_1(&Cv, Lambda,P_WT_ia,exp1,exp2,P_glob_a,L,wi,cn,Reg);
  printf("Reg=%.4g L0=%.4g L1=%.4g lik=%.4g dlik_dR=%.4g\n",
	 Reg, Lambda[0], Lambda[1], lik, Cv);

  if(OPT_REG){
    printf("# Determining Reg by maximizing Cv for model %s\n",label);
    int IT_MAX=100, iter; float x[3], y[3];
    x[1]=Reg; y[1]=Cv;
    opt_res->Lambda[0]=Lambda[0];
    opt_res->Lambda[1]=Lambda[1];
    float ymax=Cv, Rmax=Reg, x_max=100*x[1], x_min=0.01*x[1], yy;
    for(iter=0; iter<IT_MAX; iter++){
      for(int j=0; j<3; j++){
	if(j==0){Reg=x[1]*0.85;}
	else if(j==2){Reg=x[1]*1.15;}
	else{continue;}
	Analytic_Lambda_lik_1(&yy,Lambda,P_WT_ia,exp1,exp2,P_glob_a,
			      L,wi,cn,Reg);
	if(yy>ymax){
	  ymax=yy; Rmax=Reg;
	  opt_res->Lambda[0]=Lambda[0];
	  opt_res->Lambda[1]=Lambda[1];
	}
	x[j]=Reg; y[j]=yy;
      }
      Reg=Find_max_quad(x[0], x[1], x[2], y[0], y[1], y[2], x_min, x_max);
      Analytic_Lambda_lik_1(&yy,Lambda,P_WT_ia,exp1,exp2,P_glob_a,L,wi,cn,Reg);
      if(yy>ymax){
	ymax=yy; Rmax=Reg;
	opt_res->Lambda[0]=Lambda[0]; opt_res->Lambda[1]=Lambda[1];
      }
      if(ymax<=y[1])break;
      x[1]=Rmax; y[1]=ymax;
    }
    if(iter==IT_MAX)
      printf("WARNING, optimal Cv not reached in %d iterations\n", iter);
    Reg=Rmax; //Reg/=1000;
    lik=Analytic_Lambda_lik_1(&Cv,opt_res->Lambda,P_WT_ia,exp1,exp2,P_glob_a,
			      L,wi,cn,Reg);
    printf("Reg=%.4g L0=%.4g L1=%.4g lik=%.4g Cv=%.4g\n",
	   Reg, opt_res->Lambda[0], opt_res->Lambda[1], lik, Cv);
  }

  return(0);
}


float Analytic_Lambda_lik_1(float *Cv_opt, float *Lambda, float **P_WT_ia,
			    float **exp1, float **exp2, float *P_glob_a, int L,
			    float *wi, float *cn, float Reg)
{
  float conv_thr=0.000001, Reg2=2*Reg;
  int IT_MAX=200, n_fail=0, iter;
  int n, nmod; if(exp2){nmod=2;}else{nmod=1;}
  float *f_aa=NULL; //f_reg_aa; // NULL;

  float **phi[nmod]; phi[0]=exp1; if(nmod>1)phi[1]=exp2;
  float *Ave[nmod], *Corr_reg_inv[nmod], dlik_dL[nmod], dPhi_dL[nmod];
  for(n=0; n<nmod; n++){
    Corr_reg_inv[n]=malloc(nmod*sizeof(float));
    Ave[n]=malloc(L*sizeof(float));
  }

  float Zi[L];
  Compute_P(P_WT_ia,Zi,P_glob_a, Lambda,exp1,exp2,f_aa,wi,aa_seq,PMIN,L,Naa);
  float lik=Compute_lik_cn(Lambda, cn, Zi, wi, nmod, L);

  float score_opt=lik-Reg*Lambda[0]*Lambda[0];
  if(nmod>1)score_opt-=Reg*Lambda[1]*Lambda[1];
  float score_ini=score_opt; //lik_ini=lik;
  float Lambda_opt[2], Lambda_new[2];
  for(n=0; n<2; n++)Lambda_opt[n]=Lambda[n];

  float lik_opt=lik, Cv, score=score_ini;
  printf("#"); for(n=0; n<nmod; n++)printf(" Lambda%d", n);
  //for(n=0; n<nmod; n++)printf(" dPhi_dL%d", n);
  printf(" lik score Reg= %.4g\n", Reg);
  printf("%.4g", Lambda[0]);
  if(nmod>1)printf("\t%.4g ", Lambda[n]);
  printf("\t%.5g\t%.5g\n", lik, score_opt);

  for(iter=0; iter<IT_MAX; iter++){
    Inverse_corr_reg(Corr_reg_inv,Ave,P_WT_ia,phi,wi,L,Reg2,nmod,Lambda);
    for(n=0; n<nmod; n++){
      double phi_ave=0; for(int i=0; i<L; i++)phi_ave+=wi[i]*Ave[n][i];
      dlik_dL[n]=(phi_ave-cn[n]);
      dPhi_dL[n]=dlik_dL[n]-Reg2*Lambda[n];
    }
    Cv=Compute_Cv_lik(Lambda, dlik_dL, Corr_reg_inv, nmod); 
    if(score>=score_opt)*Cv_opt=Cv;

    float conv=0;
    for(n=0; n<nmod; n++){
      Lambda_new[n]=Lambda[n]+Corr_reg_inv[n][0]*dPhi_dL[0];
      if(nmod>1)Lambda_new[n]+=Corr_reg_inv[n][1]*dPhi_dL[1];
      float d=Lambda_new[n]-Lambda[n]; conv+=d*d;
    }
    Compute_P(P_WT_ia,Zi,P_glob_a,Lambda_new,exp1,exp2,
	      f_aa,wi,aa_seq,PMIN,L,Naa);
    lik=Compute_lik_cn(Lambda_new, cn, Zi, wi, nmod, L);
    score=lik-Reg*Lambda_new[0]*Lambda_new[0];
    if(nmod>1)score-=Reg*Lambda_new[1]*Lambda_new[1];

    if(nmod==1){
      printf("%.4g", Lambda_new[0]);
      //printf("\t%.5g", dPhi_dL[0]);
    }else{
      printf("%.4g\t%.4g",Lambda_new[0], Lambda_new[1]);
      //printf("\t%.5g\t%.5g",dPhi_dL[0],dPhi_dL[1]);
    }
    printf("\t%.5g\t%.5g\n",lik,score);

    //if(lik>lik_opt){
    if(score>score_opt){
      score_opt=score;
      lik_opt=lik; n_fail=0;
      for(n=0; n<nmod; n++){
	Lambda[n]=Lambda_new[n]; Lambda_opt[n]=Lambda[n];
      }
    }else if(n_fail>0){
      break;
    }else if(conv<conv_thr){
      break;
    }else{ // lik increases but no convergence
      n_fail++;
      for(n=0; n<nmod; n++)Lambda[n]=0.5*Lambda_new[n]+0.5*Lambda[n];
      Compute_P(P_WT_ia,Zi,P_glob_a, Lambda,exp1,exp2,
		f_aa,wi,aa_seq,PMIN,L,Naa);
    }
    
  }
  //printf("Reg= %.4g Optimal lambda = %.4g  %.4g lik= %.4g Cv= %.4g\n",
  //	 Reg, Lambda_opt[0], Lambda_opt[1], lik_opt, *Cv_opt);

  if(iter==IT_MAX){
    printf("WARNING, Lambda could not be optimized in %d iterations\n",iter);
  }
  if(score_opt<=score_ini)printf("WARNING, score could not be optimized\n");
  //if(lik_opt<=lik_ini)printf("WARNING, likelihood could not be optimized\n");


  for(n=0; n<nmod; n++){
    free(Corr_reg_inv[n]); free(Ave[n]);
    Lambda[n]=Lambda_opt[n];
  }
  // Recompute P with PMIN>0
  // Compute_P(P_WT_ia,Zi,P_glob_a, Lambda,exp1,exp2,f_aa,wi,aa_seq,PMIN,L,Naa);
  // lik_opt=Compute_lik_cn(Lambda, cn, Zi, wi, nmod, L);
  return(lik_opt);
}

int Compute_P(float **P_ia, float *Zi, float *P_glob_a,
	      float *Lambda, float **exp1, float **exp2,
	      float *f_aa, float *wi, short *aa_seq,
	      float Pmin, int L, int Naa)
{
  int nm; if(exp2){nm=2;}else{nm=1;}
  int i, a, zero=0, error=0;
  if(isnan(Lambda[0]) || (exp2 && isnan(Lambda[1]))){
    printf("WARNING in Compute_P, Lambda is nan: %.3g %.3g\n",
	   Lambda[0],Lambda[1]);
    for(i=0; i<2; i++)if(isnan(Lambda[i]))Lambda[i]=Lambda_start[i];
    printf("Setting Lambda to  %.3g %.3g\n", Lambda[0],Lambda[1]);
  }

  // Compute exponential
  for(i=0; i<L; i++){
    float *P_a=P_ia[i], *e1=exp1[i];
    int a_wt=aa_seq[i];
    for(a=0; a<Naa; a++){
      if(exp_str && a==a_wt)continue;
      double ee=Lambda[0]*e1[a];
      if(exp2)ee+=Lambda[1]*exp2[i][a];
      P_a[a]=exp(-ee);
    }
  }
  if(exp_str){ // Compute structural deformation of wild type a.a.
    for(i=0; i<L; i++){
      int a_wt=aa_seq[i];
      if(a_wt<0 || a_wt>=Naa){continue;}
      float *ei=exp_str[i], *Pi=P_ia[i];
      double e_ave=0, norm=0;
      for(a=0; a<Naa; a++){
	if(a!=a_wt){
	  float P=Pi[a]; if(f_glob)P*=f_glob[a]; // Weight with global freq
	  e_ave+=P*ei[a]; norm+=P;
	}
      }
      if(norm)e_ave/=norm;
      ei[a_wt]=e_ave;
      double ee=Lambda[0]*exp1[i][a_wt];
      if(exp2)ee+=Lambda[1]*exp2[i][a_wt];
      Pi[a_wt]=exp(-ee);
    }
  }

  if(OPT_MUT && f_aa){
    // Determine P_glob based on Lambda
    Get_P_glob(P_glob_a, P_ia, wi, L, f_aa, Naa);
  }
  for(i=0; i<L; i++){
    float *P_a=P_ia[i]; double Z=0;
    for(a=0; a<Naa; a++){
      P_a[a]=P_a[a]*P_glob_a[a];
      Z+=P_a[a];
    }
    if(isnan(Z) || Z<=0){
      error++;
      printf("ERROR Z, i=%d wi= %.0f Z= %.3g Lambda= %.3g",
	     i,wi[i],Z,Lambda[0]);
      if(nm==2){printf(" %.3g", Lambda[1]);} printf("\n");
      for(int n=0; n<nm; n++){
	printf("exp= ");
	float *e1; if(n==0){e1=exp1[i];}else{e1=exp2[i];}
	for(a=0; a<Naa; a++)printf(" %.2g", e1[a]);
	printf("\n");
      }
    }else{
      float sum=0;
      for(a=0; a<Naa; a++){
	P_a[a]/=Z; if(P_a[a]<Pmin){P_a[a]=Pmin; zero++; sum+=Pmin-P_a[a];}
      }
      if(sum){sum+=1; for(a=0; a<Naa; a++){P_a[a]/=sum;}}
    }
    if(Zi)Zi[i]=Z;
  }
  if(error){
    for(int n=0; n<nm; n++)printf("Lam= %.3g ",Lambda[n]);
    printf("\nP_glob:");
    for(a=0; a<Naa; a++)printf(" %.2g", P_glob_a[a]);
    printf("\n");
    exit(8);
  }
  return(zero);
}

void Inverse_corr_reg(float **Corr_reg_inv, float **Ave,
		      float **P_WT_ia, float ***phi, float *wi,
		      int L, float Reg2, int nmod, float *Lambda)
{
  float Corr_reg[nmod][nmod]; int i, k, a;
  float y[Naa];

  for(k=0; k<nmod; k++){
    for(i=0; i<L; i++){
      Ave[k][i]=Compute_average(P_WT_ia[i], phi[k][i]);
    }
    for(int l=0; l<=k; l++){
      // Correlation matrix
      double sum=0;
      for(i=0; i<L; i++){
	float *xk=phi[k][i], *xl=phi[l][i];
	for(a=0; a<Naa; a++)y[a]=xk[a]*xl[a];
	sum+=wi[i]*(Compute_average(P_WT_ia[i], y)-Ave[k][i]*Ave[l][i]);
      }
      Corr_reg[k][l]=sum;
      if(l!=k){Corr_reg[l][k]=Corr_reg[k][l];}
      else{Corr_reg[k][k]+=Reg2;}
    }
  }

  if(nmod==1){
    Corr_reg_inv[0][0]=1./Corr_reg[0][0];
  }else{
    double det=Corr_reg[0][0]*Corr_reg[1][1]-Corr_reg[0][1]*Corr_reg[1][0];
    if(det<=0){
      printf("ERROR, determinant of regularized matrix= %.3g\n", det);
      printf("2*Reg= %.3g", Reg2);
      for(k=0; k<nmod; k++)printf(" Lambda%d = %.3g",k+1,Lambda[k]);
      printf("\n");
      for(k=0; k<nmod; k++){
	for(int l=0; l<nmod; l++){
	  printf("C%d%d = %.3g ", k,l,Corr_reg[k][l]);
	}
      }
      printf("\n");
      exit(8);
    }
    Corr_reg_inv[0][0]=Corr_reg[1][1]/det;
    Corr_reg_inv[1][1]=Corr_reg[0][0]/det;
    Corr_reg_inv[0][1]=-Corr_reg[1][0]/det;
    Corr_reg_inv[1][0]=-Corr_reg[0][1]/det;
  }
}


float Compute_Cv_lik(float *Lambda, float *dlik_dL, float **Corr_reg_inv,
		     int nmod)
{
  // Cv=-(dlik/dL1 dL1/dR + dlik/dL2 dL2/dR)
  // dlik/dL=-L*cn + sum_i wi <exp_ia>_{L}
  // dL/dR =2*M^-1 L
  // -M_kl = <exp^k_ia exp^l_ia>-<exp^k_ia><exp^l_ia> + 2R delta_kl
  // Cv= dlik/dL 2M^{-1} L

  double Cv=0;
  for(int k=0; k<nmod; k++){
    float dL_dR=Corr_reg_inv[k][0]*Lambda[0];
    if(nmod>1)dL_dR+=Corr_reg_inv[k][1]*Lambda[1];
    Cv+=dlik_dL[k]*dL_dR;
  }
  Cv*=(2);
  return(Cv);

}

float Compute_average(float *P_WT_a, float *phi){
  double sum=0;
  for(int a=0; a<Naa; a++)sum+=P_WT_a[a]*phi[a];
  return(sum);
}

float Compute_lik_cn(float *Lambda, float *cn, float *Z, float *wi,
		     int nmod, int L)
{
  double lik=-Lambda[0]*cn[0];
  if(nmod==2)lik-=Lambda[1]*cn[1];
  for(int i=0; i<L; i++)lik-=wi[i]*log(Z[i]);
  return(lik+lik_const);
}


float KL_symm(float *P, float *Q, int n){
  double KL=0;
  for(int i=0; i<n; i++){
    if(P[i]<=0 || Q[i]<=0)continue;
    KL+=(P[i]-Q[i])*log(P[i]/Q[i]);
  }
  return(KL);
}

void Print_matrix(struct protein target){
  char name_out[200];
  sprintf(name_out, "%s.Contact_matrix.cm", target.name);
  FILE *file_out=fopen(name_out, "w"); int i;
  printf("Writing contact matrix in file %s\n", name_out);
  fprintf(file_out, "# %d %d %s\n",
	  target.length, target.n_cont, target.name);
  for(i=0; i<target.length; i++){
    short *Ci=target.contact[i];
    while(*Ci>=0){
      fprintf(file_out, "%d %d\n", i, *Ci); Ci++;
    }
  }
  fclose(file_out);
}

void Initialize_load(struct load *load){
  load->df_ave=0; load->df_dev=0;
  load->dG_ave=0; load->dG_dev=0;
  load->num=0;
}

void Sum_loads(struct load *load, struct load *load1){
  load->df_ave+=load1->df_ave;
  load->df_dev+=(load1->df_ave)*(load1->df_ave);
  load->dG_ave+=load1->dG_ave;
  load->dG_dev+=(load1->dG_ave)*(load1->dG_ave);
  load->num++;
}

FILE *Open_summary(char *name_file, int L, struct REM E_wt,
		   float *mut_par, short *aa_seq, int npdb)
{
  // Print summary of results
  char name[100];
  sprintf(name, "%s%s", name_file, EXT_SUM);
  printf("Writing %s\n", name);
  FILE *file_out=fopen(name, "w");
  name[5]='\0';
  // Wild type sequence
  double h_wt=0; int n=0, i;
  for(i=0; i<L; i++){if(aa_seq[i]>=0){h_wt+=hydro[aa_seq[i]]; n++;}}
  h_wt/=n;
  fprintf(file_out, "# Loc.interactions: %.3f ", SEC_STR);
  fprintf(file_out, "T= %.2f SC=%.3f SU= %.3f REM= %d T_freezing= %.2f\n",
	  E_wt.T, E_wt.S_C, E_wt.S_U, E_wt.REM, E_wt.Tf);
  if(L_ali){fprintf(file_out,"# Most stable PDB=");}
  else{fprintf(file_out,"# PDB=");}
  fprintf(file_out, " %s L_seq=%3d DG_wt= %.1f /L= %.3f hydro= %.3f\n",
	  target.name, L, E_wt.DeltaG, E_wt.DeltaG/L, h_wt);
  fprintf(file_out,
	  "# %d residues in PDB sequence, %d structured %d disordered\n",
	  L, L_str_PDB, L-L_str_PDB);

  /*fprintf(file_out, "#L=%3d   nuc_freq= %.3f %.3f %.3f %.3f", L,
	  mut_par[Code_nuc('G')], mut_par[Code_nuc('C')],
	  mut_par[Code_nuc('A')], mut_par[Code_nuc('T')]);
  fprintf(file_out, "   tt= %.1f CpG= %.1f twonuc= %.1g\n",
  mut_par[4], mut_par[5], mut_par[6]);*/

  // MSA
  if(L_ali){
    fprintf(file_out, "# MSA %s: %d sequences, %d columns\n",
	    file_ali, N_seq, L_ali);
    fprintf(file_out,
	    "# %d structures read, %d aligned with MSA, "
	    "%d eliminated because N_stab <  (%.3g,%.3g), %d used\n",
	    N_PDB_read, N_PDB, N_PDB-N_PDB_used,
	    N_stab_thr, N_stab_thr_2, N_PDB_used);
    fprintf(file_out,
	    "# Most stable structure has %.1f stable aligned sequences "
	    "in file %s Mean DG: %.4g\n",
	    N_stab_target, file_ali, DG_ave_target);

    fprintf(file_out, "# Order of PDB seq in MSA: %d\n", npdb+1);
    //fprintf(file_out, "# Mapping written in %s\n", file_map);
    fprintf(file_out, "# Mean seq. identity between PDB seq and MSA: %.3f\n",
	    Seq_id_ave_target);
    fprintf(file_out,
	    "# %d columns in MSA, %d columns not modelled, "
	    "%d not associated with stable PDB\n", L_ali,L_nomod,L_noali);
  }else{
    fprintf(file_out,"# No stable aligned sequences were found\n");
  }
  fprintf(file_out,
	  "# MSA profiles regularized as f_ia=(n(MSA_ia)+REG*f_a)/N");
  fprintf(file_out, "  REG=%.3f\n",REG_DISTR);
  fprintf(file_out,"# Entropy per site of the alignment: %.3f\n",Entr_ave);
  fprintf(file_out,"# Regularized entropy per site: %.3f\n", Entr_reg);

  // Selection parameter
  fprintf(file_out, "# Lambda is obtained ");
  if(SCORE_KL){
    fprintf(file_out,"minimizing KL(mod,obs)+KL(obs,mod)");
  }else{
    fprintf(file_out,"maximizing log(likelihood)-Reg*Lambda^2");
  }
  if(SCORE_EVO){
    fprintf(file_out," plus approximate evolutionary likelihood");
  }
  fprintf(file_out,"\n");
  if(OPT_REG){
    fprintf(file_out,"# Regularization parameter optimized by ");
    //if(SCORE_KL==0 ||SCORE_CV){
    if(1){
      fprintf(file_out,"mazimizing dlik/dREG " );
      fprintf(file_out,"Max when %d consecutive decreasing steps\n", CV_DOWN);
    }else{
      fprintf(file_out,"symmetrizing KL divergence\n");
    }

    //fprintf(file_out," initial value= %.3g\n", REG_INI);
  }else{
    fprintf(file_out,"# Fix regularization parameter Reg = %.3g\n", REG_INI);
  }
  fprintf(file_out, "# Lambda determined by ");
  if(SCORE_KL==0 || LAMBDA_ANALYTIC){
    fprintf(file_out, "solving for zero derivative\n");
  }else{
    fprintf(file_out, "numerically minimizing KL_symm\n");
  }
  fprintf(file_out, "# KL(obs,mod): -likelihood(obs|mod)-Entropy(obs)\n");
  fprintf(file_out, "# KL(mod,obs): -likelihood(mod|obs)-Entropy(mod)\n");
  fflush(file_out);
  return(file_out);
}

void Print_TN_div(short *aa_seq, short *aa_seq0, int L,
		  int num_aa_subst, FILE *file_out)
{
  float SI0=0.06;
  float SI=0;
  for(int i=0; i<L; i++)if(aa_seq[i]==aa_seq0[i])SI++;
  SI/=L; if(SI<=SI0)return;
  float TN=-log((SI-SI0)/(1-SI0));
  fprintf(file_out, "%.3f\t%d\n", TN, num_aa_subst);
}

void Print_seq(FILE *file_msa, short *aa_seq, int len_amm, int *nseq_msa,
	       float DeltaG)
{
  if(file_msa==NULL)return;

  fprintf(file_msa,">Seq%d_DG=%.4g", *nseq_msa, DeltaG);
  int j=0;
  for(int i=0; i<len_amm; i++){
    if(j==0)fprintf(file_msa,"\n");
    j++; if(j==60)j=0;
    fprintf(file_msa,"%c",AMIN_CODE[aa_seq[i]]);
  }
  fprintf(file_msa,"\n");
  (*nseq_msa)++;
}

void Regularize(float **f_reg_ia, float **n_msa_ia, float *f_aa,
		float w_max, float reg, int len_amm, int Naa)
{
  int nchange=0, a;
  for(int i=0; i<len_amm; i++){
    double sum=0; //float reg=w_max-wi[i]+w_max*REG_DISTR;
    for(a=0; a<Naa; a++){
      f_reg_ia[i][a]=n_msa_ia[i][a]+reg*f_aa[a];
      //f_reg_ia[i][a]=n_msa_ia[i][a]/wi[i]+reg*f_aa[a];
      sum+=f_reg_ia[i][a];
    }
    if(sum){
      int change=0;
      for(a=0; a<Naa; a++){
	f_reg_ia[i][a]/=sum;
	if(f_reg_ia[i][a]<PMIN){f_reg_ia[i][a]=PMIN; change++;}
      }
      if(change){
	sum=0; for(a=0; a<Naa; a++){sum+=f_reg_ia[i][a];}
	for(a=0; a<Naa; a++){f_reg_ia[i][a]/=sum;}
	nchange+=change;
      }
    }else{
      for(a=0; a<Naa; a++){f_reg_ia[i][a]=f_aa[a];}
      nchange+=20;
    }
  }

  printf("%.1f amino acids absent on average at each site\n",
	 (float)nchange/len_amm);
}

void Compute_nuc_mut(char *dna_seq, int len_dna,
		     short *aa_seq, int len_amm,
		     char **codon, char *coded_aa,
		     char *SSC_TYPE, float **exp1, float **exp2,
		     float *Lambda, FILE *file_out)
{
  int n, nmod=1; if(exp2)nmod=2;

  int a, b, i, j;
  double DDG_sum[4][4]; int m_sum[4][4];
  for(a=0; a<4; a++){
    for(b=0; b<4; b++){DDG_sum[a][b]=0; m_sum[a][b]=0;}
  }

  int in=0; char New_cod[3];
  for(i=0; i<len_amm; i++){
    float DDG_max[2]; DDG_max[1]=0;
    for(n=0; n<nmod; n++){
      float *ee; if(n==0){ee=exp1[i];}else{ee=exp2[i];}
      DDG_max[n]=ee[0];
      for(j=1; j<Naa; j++)if(ee[j]>DDG_max[n])DDG_max[n]=ee[j];
    }
    for(j=0; j<3; j++)New_cod[j]=dna_seq[in+j];
    for(j=0; j<3; j++){
      char old_nuc=dna_seq[in];
      a=Code_nuc(old_nuc);
      for(b=0; b<4; b++){
	if(b==a)continue;
	New_cod[j]=NUC_CODE[b];
	int aa_new=Coded_aa(New_cod, codon, coded_aa);
	if(aa_new<0){
	  DDG_sum[a][b]+=Lambda[0]*(len_amm-i)*DDG_max[0];
	  if(exp2)DDG_sum[a][b]+=Lambda[1]*(len_amm-i)*DDG_max[1];
	}else{
	  DDG_sum[a][b]+=Lambda[0]*exp1[i][aa_new];
	  if(exp2)DDG_sum[a][b]+=Lambda[1]*exp2[i][aa_new];
	}
	m_sum[a][b]++;
      }
      New_cod[j]=old_nuc; in++;
    }
  }

  fprintf(file_out, "# Effect of mutations of type %s\n", SSC_TYPE);
  fprintf(file_out,"#a\tb\tDDG(a,b)/Ldna\tnmut\tDDG(a,b)/nmut\n");
  fprintf(file_out, "#DDG represents the total estimated fitness cost of ");
  fprintf(file_out, " mutations a->b, whose number is nmut\n");
  fprintf(file_out, "# Ldna= %d Lambda= %.3g", len_dna, Lambda[0]);
  if(exp2)fprintf(file_out, " %.3g", Lambda[1]);
  fprintf(file_out, "\n");

  for(a=0; a<4; a++){
    for(b=0; b<4; b++){
      if(b==a)continue;
      float df=DDG_sum[a][b]/len_dna;
      fprintf(file_out, "%c\t%c", NUC_CODE[a], NUC_CODE[b]);
      fprintf(file_out, "\t%.5f", df);
      fprintf(file_out, "\t%d", m_sum[a][b]);
      fprintf(file_out, "\t%.5f\n", DDG_sum[a][b]/m_sum[a][b]);
    }
  }

}

int Match_dna(char *dna_seq, int nseq, int *ini_seq, int *len_seq,
	      struct protein pdb, char **codon, char *coded_aa)
{
  int nchain=pdb.nchain, ic, i;
  if(nchain!=nseq){
    printf("ERROR, %d different chains but %d dna sequences,", nchain, nseq);
    printf(" discarding the DNA sequences\n"); return(0);
  }
  int match[nchain]; for(int i=0; i<nchain; i++)match[i]=-1;
  char *aa_trans[nseq]; int len_trans[nseq];
  for(i=0; i<nseq; i++){
    len_trans[i]=Translate_aa(&(aa_trans[i]),dna_seq+ini_seq[i],len_seq[i],
			      codon, coded_aa);
    printf("Translated sequence %d: ",i);
    for(int j=0; j<len_trans[i]; j++)printf("%c", aa_trans[i][j]);
    printf("\n");
  }

  for(ic=0; ic<nchain; ic++){
    short *aa_seq=pdb.aa_seq+pdb.ini_chain[ic];
    int L=pdb.len_chain[ic];
    char aa_char[L];
    for(i=0; i<L; i++)aa_char[i]=AMIN_CODE[aa_seq[i]];
    for(i=0; i<nseq; i++){
      if(Compare_amm(aa_char, L, aa_trans[i], len_trans[i])){
	match[ic]=i; break;
      }
    }
    if(match[ic]<0){
      printf("ERROR, no match found for PDB seq %d %d a.a.\n", ic, len_amm);
      printf("Candidate matches (a.a.): ");
      for(i=0; i<nseq; i++)printf(" %d", len_seq[i]/3);
      printf("\n");
      return(0);
    }
  }
  int equal=1;
  for(ic=0; ic<nchain; ic++)if(match[ic]!=ic){equal=0; break;}
  if(equal)return(1); // Order must not be changed

  // Change order of sequences
  int len_dna=0; for(i=0; i<nseq; i++)len_dna+=len_seq[i];
  char seq_tmp[len_dna]; int len_tmp[nseq];
  for(i=0; i<len_dna; i++)seq_tmp[i]=dna_seq[i];
  for(i=0; i<nseq; i++){len_tmp[i]=len_seq[i];}
  int ini=0; 
  for(ic=0; ic<nchain; ic++){
    int idna=match[ic], j=ini;
    ini_seq[ic]=ini; len_seq[ic]=len_tmp[idna];
    for(i=0; i<len_seq[ic]; i++){dna_seq[j]=seq_tmp[i]; j++;}
    ini+=len_seq[ic];
  }
  return(1);
}

float Normalize_exponent(float **SSC_exp, int L, char *model){
  int i, a; double sum=0;
  for(i=0; i<L; i++){
    float *SSC=SSC_exp[i], min_SSC=SSC[0];
    for(a=1; a<Naa; a++)if(SSC[a]<min_SSC)min_SSC=SSC[a];
    for(a=0; a<Naa; a++){SSC[a]-=min_SSC; sum+=SSC[a];}
  }
  sum/=(Naa*L);
  printf("Average minus fitness %s model: %.3g L=%d", model, sum, L);
  if(NORMALIZE)printf(" After normalization: 1.");
  printf("\n");

  if(isnan(sum) || sum<0){
    printf("WARNING, error in Normalize_exponent\n");
    for(i=0; i<L; i++){
      for(a=0; a<Naa; a++)if(SSC_exp[i][a]<0 || isnan(SSC_exp[i][a]))break;
      if(a<Naa){
	printf("i= %d Q: ", i);
	for(a=0; a<Naa; a++)printf(" %.2g", SSC_exp[i][a]);
	printf("\n");
      }
    }
    return(-1);
  }
  if(NORMALIZE==0)return(1);
  if(sum==0)return(0);

  float max_exp=0;
  for(i=0; i<L; i++){
    float *SSC=SSC_exp[i];
    for(a=0; a<Naa; a++){
      SSC[a]/=sum;
      if(SSC[a]>max_exp)max_exp=SSC[a];
    }
  }
  return(max_exp);
}

///// Routines needed if REGULARIZE==0
int Maximize_Lambda(struct MF_results *opt_res, float **P_WT_ia,
		    float *P_glob_a, float **exp1, float **exp2,
		    float **f_reg_ia, float *wi, int L, int Naa,
		    int repeat)
{
  int REDO=0, PRINT=0;
  float STEP=0.1, step_rate=0.75, step, eps=0.000001;
  // Step of lambda and its reduction at each iteration
  int IT_MAX=200, iter, k, mid;

  int redo; if(REDO && repeat){redo=1;}else{redo=0;}

  float x[4], y[4], Lam[2];
  int n, nmod; if(exp2){nmod=2;}else{nmod=1;}
  struct MF_results res, tmp_res;
  opt_res->score_reg=-100000;

  if(Sample_Lambda){
    Sample_Lambda=0;
    printf("# Lambda score\n");
    float Lambda=LAM_MIN;
    for(int i=0; i<N_SAM_LAM; i++){
      Lambda_sam[i][0]=Lambda; 
      if(nmod==2){Lambda_sam[i][1]=Lambda;}else{Lambda_sam[i][1]=0;}
      Lambda+=LAM_STEP;
      Compute_P(P_WT_ia,NULL,P_glob_a,Lambda_sam[i],
		exp1,exp2,f_reg_aa,wi,aa_seq,PMIN,L,Naa);
      Compute_score(&res,P_WT_ia,Lambda_sam[i],f_reg_ia,wi,L,Naa,0);
      score_sam[i]=res.score;
      printf("%.3f %.3f\n", Lambda_sam[i][0], score_sam[i]);
    }
  }

  float score_max=-1000000; Lam[0]=0; Lam[1]=0;
  for(int i=0; i<N_SAM_LAM; i++){
    res.score_reg=score_sam[i]-Reg*Lambda_sam[i][0]*Lambda_sam[i][0];
    if(nmod==2)res.score_reg -=Reg*Lambda_sam[i][1]*Lambda_sam[i][1];
    if(i==0 || res.score_reg > score_max){
      Lam[0]=Lambda_sam[i][0];
      Lam[1]=Lambda_sam[i][1];
      score_max=res.score_reg;
    }
  }
  tmp_res.Lambda[0]=Lam[0];
  tmp_res.Lambda[1]=Lam[1];
  tmp_res.score_reg=score_max;

 restart:
  step=STEP; mid=0; 
  // Initialize

  int n_down=0;
  if(PRINT){
    printf("Reg= %.3f\nscore=%.4f %.3f %.3f %d\n",Reg, tmp_res.score_reg,
	   tmp_res.Lambda[0],tmp_res.Lambda[1],n_down);
  }

  for(iter=0; iter<IT_MAX; iter++){
    if(n_down>=2)break;
    for(n=0; n<nmod; n++){
      x[1]=tmp_res.Lambda[n];
      y[1]=tmp_res.score_reg;
      for(k=0; k<4; k++){
	if(k==1){
	  continue;
	}else if(k==0){
	  Lam[n]=x[1]*(1-step);
	}else if(k==2){
	  Lam[n]=x[1]*(1+step);
	}else{ // k==3
	  if(y[1] > y[0] && y[1] > y[2]){mid++; if(mid==nmod)break;}
	  else{mid=0;}
	  Lam[n]=Find_max_quad(x[0],x[1],x[2],y[0],y[1],y[2],0,20);
	  if(isnan(Lam[n])){printf("ERROR in find_max\n"); Lam[n]=x[1]*1.1;}
	}
	Compute_P(P_WT_ia,NULL,P_glob_a, Lam,exp1,exp2,
		  f_reg_aa,wi,aa_seq,PMIN,L,Naa);
	Compute_score(&res, P_WT_ia, Lam, f_reg_ia, wi, L, Naa,0);
	//if(res.entropy<=Entropy_min)break;
	if(isnan(res.score_reg))break;
	x[k]=Lam[n];
	y[k]=res.score_reg;
	if(k==3){ // k=3
	  if(res.score_reg > tmp_res.score_reg){n_down=0;}
	  else{n_down++;}
	  //if(OPT_REG==0)
	  //printf("%.4f %.3f %.3f %d\n",res.score_reg,Lam[0],Lam[1],n_down);
	  //if(fabs(res.score_reg-tmp_res.score_reg)<eps){mid++;}else{mid=0;}
	}
	if(res.score_reg > tmp_res.score_reg){
	  tmp_res=res; tmp_res.Lambda[0]=Lam[0]; tmp_res.Lambda[1]=Lam[1]; 
	}
      } // end k (4 samples for quadratic optimization)
      if(mid==nmod)break;
      Lam[n]=tmp_res.Lambda[n];
    } // end n (model components)
    if(step>0.01)step*=step_rate;
  } // end iter
  if(iter==IT_MAX)
    printf("WARNING, optimization of Lambda did not converge\n");
  if(PRINT){
    printf("%.4f %.3f %.3f %d\n",tmp_res.score_reg,
	   tmp_res.Lambda[0],tmp_res.Lambda[1],n_down);
  }

  // Optimize by bisection
  float L0[2], L1[2], y0[2], y1[2], *yy, *LL;
  for(n=0; n<nmod; n++)Lam[n]=tmp_res.Lambda[n];
  for(n=0; n<nmod; n++){
    if(n==0){LL=L0; yy=y0;}else{LL=L1; yy=y1;}
    Lam[n]=tmp_res.Lambda[n]; LL[1]=Lam[n]; yy[1]=tmp_res.score_reg;
    for(k=0; k<=2; k+=2){
      if(k==0){LL[0]=LL[1]*(1-step); Lam[n]=LL[0];}
      else{LL[2]=LL[1]*(1+step); Lam[n]=LL[2];}
      Compute_P(P_WT_ia,NULL,P_glob_a, Lam,exp1,exp2,
		f_reg_aa,wi,aa_seq,PMIN,L,Naa);
      Compute_score(&res, P_WT_ia, Lam, f_reg_ia, wi, L, Naa,0);
      yy[k]=res.score_reg;
    }
    if(yy[1] < yy[0] || yy[1] < yy[2]){
      printf("WARNING, optimization of Lambda by bisection not possible\n");
      printf("n=%d x: %.6g %.6g %.6g y: %.6g %.6g %.6g it= %d d= %d mid= %d\n",
	     n, LL[0],LL[1],LL[2],yy[0],yy[1],yy[2],iter,n_down,mid);
      Lam[n]=LL[1];
      goto last_lambda;
    }
  }
  
  // Look for optimum dividing intermediate interval
  mid=0; int wrong=0;
  for(iter=0; iter<30; iter++){
    for(n=0; n<nmod; n++){
      if(n==0){LL=L0; yy=y0;}else{LL=L1; yy=y1;}
      if(yy[1]<yy[0] || yy[1]<yy[2] || LL[0]>LL[1] || LL[1]>LL[2]){
	wrong=1; break;
      }
      float xa=0.5*(LL[0]+LL[1]); Lam[n]=xa;
      Compute_P(P_WT_ia,NULL,P_glob_a, Lam,exp1,exp2,
		f_reg_aa,wi,aa_seq,PMIN,L,Naa);
      Compute_score(&res, P_WT_ia, Lam, f_reg_ia, wi, L, Naa,0);
      if(res.score_reg > tmp_res.score_reg){tmp_res=res; tmp_res.Lambda[n]=xa;}
      float ya=res.score_reg;
      float xb=0.5*(LL[1]+LL[2]); Lam[n]=xb;
      Compute_P(P_WT_ia,NULL,P_glob_a, Lam,exp1,exp2,
		f_reg_aa,wi,aa_seq,PMIN,L,Naa);
      Compute_score(&res, P_WT_ia, Lam, f_reg_ia, wi, L, Naa,0);
      if(res.score_reg > tmp_res.score_reg){tmp_res=res; tmp_res.Lambda[n]=xb;}
      float yb=res.score_reg;
      if(ya > yb && ya > yy[1]){
	LL[2]=LL[1]; yy[2]=yy[1]; LL[1]=xa; yy[1]=ya; // max at xa
      }else if(yb > ya && yb > yy[1]){
	LL[0]=LL[1]; yy[0]=yy[1]; LL[1]=xb; yy[1]=yb; // max at xb
      }else{
	LL[0]=xa; yy[0]=ya; LL[2]=xb; yy[2]=yb; // max at x1
      }
      Lam[n]=LL[1];
      if(fabs(yy[1]-yy[2])< eps || fabs(yy[1]-yy[0])< eps){
	mid++; if(mid==nmod)break;
      }else{
	mid=0;
      }
    }
    if(wrong)break;
  }

 last_lambda:
  if(tmp_res.score_reg > opt_res->score_reg){
    *opt_res=tmp_res;
    for(n=0; n<nmod; n++){
      opt_res->Lambda[n]=tmp_res.Lambda[n];
      if(isnan(opt_res->Lambda[n]))
	printf("ERROR in numeric, Lambda[%d]= %.3g\n",n,opt_res->Lambda[n]);
    }
  }
  if(PRINT){
    printf("%.4f %.3f %.3f %d\n",tmp_res.score_reg,
	   tmp_res.Lambda[0],tmp_res.Lambda[1],n_down);
  }
  if(redo){
    redo=0; 
    Lam[0]=tmp_res.Lambda[0]*1.5;
    if(nmod==2){Lam[1]=tmp_res.Lambda[1]*1.5;}
    Compute_P(P_WT_ia,NULL,P_glob_a,Lam,
	      exp1,exp2,f_reg_aa,wi,aa_seq,PMIN,L,Naa);
    Compute_score(&tmp_res, P_WT_ia,Lam, f_reg_ia, wi, L, Naa,0);
    goto restart;
  }

  Compute_P(P_WT_ia,NULL,P_glob_a,opt_res->Lambda,
	    exp1,exp2,f_reg_aa,wi,aa_seq,PMIN,L,Naa);
  Compute_score(opt_res, P_WT_ia, opt_res->Lambda,f_reg_ia, wi, L, Naa,0);

  return(0);
}

int Analytic_Lambda_KL(struct MF_results *opt_res, float **P_WT_ia,
		       float *P_glob_a, float **exp1, float **exp2,
		       float **f_reg_ia, float *wi, int L, int Naa)
{
  int IT_MAX=200, iter; float eps=0.000001;
  float Lambda_opt[2], KL_opt;
  float *f_aa=NULL; //f_reg_aa; // =NULL;
  
  int n, nmod; if(exp2){nmod=2;}else{nmod=1;}
  double wQFreg[2], **lFreg[2]; int i,a;
  for(n=0; n<nmod; n++){
    wQFreg[n]=0; lFreg[n]=Allocate_mat2_d(L, Naa);
  }
  for(i=0; i<L; i++){
    for(n=0; n<nmod; n++){
      double QF=0, *lF=lFreg[n][i]; float *freg_i=f_reg_ia[i];
      float *qi; if(n==0){qi=exp1[i];}else{qi=exp2[i];}
      for(a=0; a<Naa; a++){
	QF+=freg_i[a]*qi[a];
	lF[a]=log(freg_i[a]/P_glob_a[a]);
      }
      wQFreg[n]+=wi[i]*QF;
    }
  }

  float Lambda_new[2], Lambda[2];
  for(n=0; n<2; n++){
    Lambda[n]=Lambda_start[n];
    Lambda_new[n]=Lambda[n];
    Lambda_opt[n]=Lambda[n];
  }
  double KL_old=-1;
  for(iter=0; iter<IT_MAX; iter++){
    double wQ1[2], wQ2[2], wPF[2], wQQ=0, wQ1Q2=0, KL=0;
    for(n=0; n<nmod; n++){wQ1[n]=0; wQ2[n]=0; wPF[n]=0;}
    Compute_P(P_WT_ia,NULL,P_glob_a, Lambda,exp1,exp2,
	      f_aa,wi,aa_seq,0,L,Naa);
    for(i=0; i<L; i++){
      float *P_i=P_WT_ia[i]; double q11=0;
      for(n=0; n<nmod; n++){
	double q1=0, q2=0, PF=0, PFq=0, q1q2=0;
	float *qi, *qqi=NULL; double *lF=lFreg[n][i];
	if(n==0){qi=exp1[i];}else{qi=exp2[i]; qqi=exp1[i];}
	for(a=0; a<Naa; a++){
	  double Pq=P_i[a]*qi[a]; q1+=Pq; q2+=Pq*qi[a];
	  double Pl=P_i[a]*lF[a]; PF+=Pl; PFq+=Pl*qi[a];
	  if(n)q1q2+=Pq*qqi[a];
	}
	wQ1[n]+=wi[i]*q1;
	wQ2[n]+=wi[i]*(q2-q1*q1);
	wPF[n]+=wi[i]*(q1*PF-PFq);
	if(n){wQQ+=wi[i]*q1q2; wQ1Q2+=wi[i]*q1*q11;}
	else{q11=q1;}
      }
      double KLi=KL_symm(P_i, f_reg_ia[i], Naa); KL+=wi[i]*KLi;
    }
    double det, wQ12, B[2];
    if(nmod==1){
      det=wQ2[0];
    }else{
      wQ12=wQQ-wQ1Q2;
      det=wQ2[0]*wQ2[1]-wQ12*wQ12;
    }
    for(n=0; n<nmod; n++){
      B[n]=(wPF[n]+wQ1[n]-wQFreg[n]);
      if(wQ2[n]<=0){
	printf("ERROR in Lambda computation, n=%d Lambda=%.3g variance %.3g\n",
	       n, Lambda[n], wQ2[n]); break;
      }
    }
    for(n=0; n<nmod; n++){
      if(nmod==1){Lambda_new[0]=B[0];}
      else if(n==0){Lambda_new[0]=B[0]*wQ2[1]-B[1]*wQ12;}
      else if(n==1){Lambda_new[1]=B[1]*wQ2[0]-B[0]*wQ12;}
      Lambda_new[n]/=det;
      if(Lambda_new[n]<0){
	/*printf("ERROR, n=%d det=%.3g Lambda_new= %.3g, Lambda= %.3g",
	       n, det, Lambda_new[n], Lambda[n]);
	Lambda_new[n]=Lambda[n]/2;
	printf(" -> %.3g\n", Lambda_new[n]); */
	printf("WARNING, n=%d det=%.3g Lambda_new= %.3g, Lambda= %.3g",
	       n, det, Lambda_new[n], Lambda[n]);
	// break; //exit(8);
      }
    }

    if(OPT_REG==0){
      printf("%.6g %.6g", KL, Lambda[0]);
      if(exp2)printf(" %.6g", Lambda[1]);
      printf("\n");
    }
    if(iter==0 || KL<KL_opt){
      KL_opt=KL; for(n=0; n<nmod; n++)Lambda_opt[n]=Lambda[n];
    }
    if(fabs(KL-KL_old)<eps*KL_old)break;
    if(KL >= KL_old){ // KL increases
      for(n=0; n<nmod; n++)Lambda_new[n]=0.3*Lambda_new[n]+0.7*Lambda[n];
    }
    KL_old=KL;
    for(n=0; n<nmod; n++)Lambda[n]=Lambda_new[n];
  }
  printf("Optimizing Lambda, KL=%.4g L1=%.4g", KL_opt, Lambda_opt[0]);
  if(exp2)printf(" L2=%.4g", Lambda_opt[1]);
  printf("\n");
  if(iter==IT_MAX){
    printf("WARNING, Lambda could not be optimized in %d iterations\n",iter);
  }
  for(n=0; n<nmod; n++)Empty_matrix_d(lFreg[n], L);
  for(n=0; n<nmod; n++)opt_res->Lambda[n]=Lambda_opt[n];
  return(0);
}

int Optimize_reg(struct MF_results *opt_res, float **P_WT_ia,
		 float *P_glob_a, float **exp1, float **exp2, 
		 float **f_reg_ia, float *wi, int L, int Naa)
//float **n_msa_ia, 
//		 struct REM *E_wt, float **C_nat, int *i_sec, char *c_sec,
//		 char *name_file,FILE *file_summary,char *label,
//		 int print)
{
  // Determine Reg that maximizes Cv=-dscore/dReg
  printf("Determining Reg by maximizing Cv=-dscore/dReg\n");
  float EPS2=EPS_REG/2, Reg1;
  Sample_Lambda=1; Reg=REG_MIN;
  float Cv_reg[N_REG], Reg_sam[N_REG];  
  struct MF_results res; int imax=-1;
  Reg=REG_MIN;
  int n_down=0;
  for(int i=0; i<N_REG; i++){
    Maximize_Lambda(&res,P_WT_ia,P_glob_a,exp1,exp2,f_reg_ia,wi,L,Naa,1);
    float score1=res.score;
    //printf("Reg= %.4f score= %.4f\n", Reg, score1);
    Reg1=Reg; Reg+=EPS_REG;
    Maximize_Lambda(&res,P_WT_ia,P_glob_a,exp1,exp2,f_reg_ia,wi,L,Naa,1);
    float score2=res.score;
    //printf("Reg= %.4f score= %.4f\n", Reg, score2);
    Reg_sam[i]=Reg1+EPS2; // Mean of Reg1 and Reg2
    Cv_reg[i]=(score1-score2)/EPS_REG; // -dscore/dReg
    printf("Reg= %.4f score= %.4f Cv= %.3g\n",
	   Reg_sam[i], (score1+score2)/2, Cv_reg[i]);
    if(i==0 || Cv_reg[i] > Cv_reg[imax]){imax=i; n_down=0;}
    else{n_down++; if(n_down>=CV_DOWN && i>=2)break;}
    Reg=Reg1+REG_STEP;
  }
  float x[3], y[3]; int i1;
  if(imax==0){i1=1;}
  else if(imax==N_REG-1){i1=N_REG-2;}
  else{i1=imax;}
  for(int i=0; i<3; i++){
    int k; if(i==0){k=i1-1;}else if(i==1){k=i1;}else{k=i1+1;}
    x[i]=Reg_sam[k]; y[i]=Cv_reg[k];
  }
  Reg1=Find_max_quad(x[0], x[1], x[2], y[0], y[1], y[2], 0, 10);
  Reg=Reg1+EPS_REG;
  Maximize_Lambda(&res,P_WT_ia,P_glob_a,exp1,exp2,f_reg_ia,wi,L,Naa,1);
  float score2=res.score;
  Reg=Reg1;
  Maximize_Lambda(opt_res,P_WT_ia,P_glob_a,exp1,exp2,f_reg_ia,wi,L,Naa,1);
  float score1=opt_res->score;
  if(SCORE_EVO){
    SCORE_EVO_TMP=1;
    Sample_Lambda=1;
    printf("Adding the evolutionary score\n"); 
    Maximize_Lambda(opt_res,P_WT_ia,P_glob_a,exp1,exp2,f_reg_ia,wi,L,Naa,1);
  }
  printf("Optimal Reg: %.4f score= %.4f Cv= %.3f Lambda= %.3f %.3f\n",
	 Reg1, opt_res->score, (score1-score2)/EPS_REG,
	 opt_res->Lambda[0], opt_res->Lambda[1]);
  //exit(8);
  return(0);
}

double Optimize_reg_old(struct MF_results *opt_res, float **P_opt_ia,
			float w_max, float *f_aa,
			struct MF_results *res, float **P_ia,
			float **exp1,float **exp2, float *P_glob_a,
			float **f_reg_ia, float **n_msa_ia,
			float *wi, struct REM *E_wt,
			float **C_nat, int *i_sec, char *c_sec,
			char *name_file, FILE *file_summ, char *label,
			int L, int Naa)
{
  float step=0.01, reg_max=0.4; int itmax=1;
  int N_reg=reg_max/step, k;
  float reg=-step/2, reg_step=step/REG_INI;
  float Lambda[N_reg][2], L_opt[2], reg_opt=-1, score_opt=-10000;
  float reg_k[N_reg], score_k[N_reg], lik_k[N_reg], symm_k[N_reg], KL[N_reg];

  char name_out[200]; sprintf(name_out, "%s.Cv.dat", name_file);
  FILE *file_out=fopen(name_out, "w");
  printf("Writing %s\n", name_out);

  int nmod=1, n; if(exp2){nmod=2;}else{L_opt[1]=0; res->Lambda[1]=0;}
  for(n=0; n<nmod; n++)res->Lambda[n]=opt_res->Lambda[n];

  char out[200], tmp[80];
  strcpy(out, "# Optimizing the parameter reg by ");
  if(SCORE_CV){
    strcat(out," maximizing d(lik_MSA)/d(reg)\n");
  }else{
    strcat(out, "minimizing |KL_mod-KL_reg|\n");
  }
  fprintf(file_summ, "%s", out);
  fprintf(file_out, "%s", out);
  fprintf(file_out, "#1=reg 2=Cv 3=lik 4=KL 5=|KL_mod-KL_reg| 6=Lambda");
  if(nmod==2)fprintf(file_out, " 7=Lambda2");
  fprintf(file_out, "\n");

  for(k=0; k<N_reg; k++){
    reg+=reg_step; reg_k[k]=reg;
    Compute_score_reg(res, reg, w_max, f_aa, P_ia, exp1, exp2, P_glob_a,
		      f_reg_ia,n_msa_ia,wi, E_wt, C_nat, i_sec, c_sec,
		      name_file, file_summ, label, L, Naa);
    lik_k[k]=res->lik_MSA;
    symm_k[k]=fabs(res->KL_mod-res->KL_reg);
    KL[k]=res->KL_mod+res->KL_reg;
    for(n=0; n<nmod; n++)Lambda[k][n]=res->Lambda[n];
  }
  for(k=0; k<(N_reg-1); k++){
    float Cv;
    if(k){
      Cv=(lik_k[k-1]-lik_k[k+1])/(2*reg_step);
    }else{
      Cv=(lik_k[k]-lik_k[k+1])/reg_step;
    }
    if(SCORE_CV){score_k[k]=Cv;}else{score_k[k]=-symm_k[k];}
    fprintf(file_out, "%.4g\t%.4g\t%.4g\t%.4g\t%.4g\t%.4g",
	    reg_k[k],Cv,lik_k[k],KL[k],symm_k[k],Lambda[k][0]);
    if(exp2)fprintf(file_out, "\t%.4g",Lambda[k][1]);
    fprintf(file_out, "\n");
  }

  // Find maximum by bisection
  //Lambda_start=Lambda[0]; res->Lambda=Lambda_start;
  for(k=0; k<(N_reg-1); k++){	    

    float Lambda_thr[2];
    for(n=0; n<nmod; n++)
      if(k){Lambda_thr[n]=0.5*Lambda[k-1][n];}else{Lambda_thr[n]=0;}
    if(k && (Lambda[k][0]<=Lambda_thr[0]||Lambda[k][1]<=Lambda_thr[1])
       && (lik_k[k-1]>lik_k[k+1])){ // Cv>0
      // phase transition! Small Lambda, mutational distribution
      fprintf(file_out,"# Sudden drop in Lambda at reg= %.2g",
	      reg_k[k]);
      for(n=0; n<nmod; n++)
	fprintf(file_out," %.4g -> %.4g",Lambda[k-1][n], Lambda[k][n]);
      fprintf(file_out," exiting\n"); break; 
    }

    if(score_k[k]>score_opt){
      score_opt=score_k[k]; reg_opt=reg_k[k];
      for(n=0; n<nmod; n++)L_opt[n]=Lambda[k][n];
    }

    if(k && (score_k[k]>0 || SCORE_CV==0) &&
       score_k[k]>score_k[k-1] && score_k[k]>score_k[k+1]){
      int i=k;float x2=reg_k[i], y2=score_k[i];
      i=k-1;  float x1=reg_k[i];
      i=k+1;  float x3=reg_k[i];
      for(int iter=0; iter<itmax; iter++){
	float xa=0.5*(x1+x2), ya, Cv, symm;
	Cv=Compute_Cv(res, xa, reg_step, w_max,f_aa,
		      P_ia,exp1,exp2,P_glob_a,f_reg_ia,n_msa_ia,wi,
		      E_wt,C_nat,i_sec, c_sec, name_file, file_summ,
		      label, L, Naa);
	symm=fabs(res->KL_mod-res->KL_reg);
	if(SCORE_CV==0){ya=-symm;}else{ya=Cv;}
	float La[2]; for(n=0; n<nmod; n++)La[n]=res->Lambda[n]; 

	fprintf(file_out, "%.4g\t%.4g\t%.4g\t%.4g\t%.4g\t%.4g",
		xa,Cv,res->lik_MSA,-res->score_reg,symm,res->Lambda[0]);
	if(exp2)fprintf(file_out, "\t%.4g",res->Lambda[1]);
	fprintf(file_out, "\n");
	sprintf(out, "# reg_par= %.5g score= %.3g Lambda_end= %.3g",
		xa, ya, res->Lambda[0]);
	if(exp2){
	  sprintf(tmp, "\t%.4g",res->Lambda[1]); strcat(out, tmp);
	}
	fprintf(file_summ, "%s\n", out); printf("%s\n", out);

	float xb=0.5*(x2+x3), yb;
	Cv=Compute_Cv(res, xb, reg_step, w_max,f_aa,
		      P_ia,exp1,exp2,P_glob_a, f_reg_ia,n_msa_ia,wi,
		      E_wt,C_nat,i_sec, c_sec, name_file, file_summ,
		      label, L, Naa);
	symm=fabs(res->KL_mod-res->KL_reg);
	if(SCORE_CV==0){yb=-symm;}else{yb=Cv;}
	float Lb[2]; for(n=0; n<nmod; n++)Lb[n]=res->Lambda[n]; 

	fprintf(file_out, "%.4g\t%.4g\t%.4g\t%.4g\t%.4g\t%.4g",
		xb,Cv,res->lik_MSA,-res->score_reg,symm,res->Lambda[0]);
	if(exp2)fprintf(file_out, "\t%.4g",res->Lambda[1]);
	fprintf(file_out, "\n");
	sprintf(out, "# reg_par= %.5g score= %.3g Lambda_end= %.3g",
		xb, yb, res->Lambda[0]);
	if(exp2){
	  sprintf(tmp, "\t%.4g",res->Lambda[1]); strcat(out, tmp);
	}
	fprintf(file_summ, "%s\n", out); printf("%s\n", out);

	if(ya > yb && ya > y2){
	  x3=x2; x2=xa; y2=ya; // max at xa
	  for(n=0; n<nmod; n++)L_opt[n]=La[n];
	}else if(yb > ya && yb > y2){
	  x1=x2; x2=xb; y2=yb; // max at xb
	  for(n=0; n<nmod; n++)L_opt[n]=Lb[n];
	}else{
	  x1=xa; x3=xb; // max at x2
	}
	if(y2>score_opt &&
	   L_opt[0]>Lambda_thr[0] && L_opt[1]>=Lambda_thr[1]){
	  score_opt=y2; reg_opt=x2;
	}
      }
    } // end bisection
  }

  if(SCORE_CV==0 || reg_opt>=0){
    sprintf(out,"# optimal regularization parameter: %.3g", reg_opt);
    for(n=0; n<nmod; n++)if(L_opt[n]>0)Lambda_start[n]=L_opt[n];
  }else{
    sprintf(out,"# WARNING, optimal regularization parameter");
    sprintf(tmp," could not be determined, using default %.3g\n#",
	    REG_INI);
    strcat(out, tmp);
    reg_opt=REG_INI;
  }

  sprintf(tmp," initial value: %.4g\n", REG_INI); strcat(out, tmp);
  fprintf(file_summ,"%s",out); printf("%s",out); fprintf(file_out,"%s",out);
  Compute_score_reg(res, reg_opt, w_max, f_aa, P_ia, exp1, exp2, P_glob_a,
		    f_reg_ia,n_msa_ia,wi, E_wt, C_nat, i_sec, c_sec,
		    name_file, file_summ, label, L, Naa);
  float symm=fabs(res->KL_mod-res->KL_reg), Cv=-1;
  if(SCORE_CV==0){score_opt=symm;}else{Cv=score_opt;}
  *opt_res=*res;
  for(n=0; n<nmod; n++)L_opt[n]=res->Lambda[n];

  fprintf(file_out, "%.4g\t%.4g\t%.4g\t%.4g\t%.4g\t%.4g",
	  reg_opt, Cv, res->lik_reg,
	  res->KL_mod+res->KL_reg, symm, res->Lambda[0]);
  if(exp2)fprintf(file_out, "\t%.4g",res->Lambda[1]);
  fprintf(file_out, "\n");

  sprintf(out,"# score_max= %.4g\n", score_opt);
  sprintf(tmp,"# optimal Lambda= %.4g %.4g starting: %.4g %.4g\n",
	  L_opt[0], L_opt[1], Lambda_start[0], Lambda_start[1]);
  strcat(out, tmp);
  fprintf(file_summ,"%s",out); printf("%s",out);
  fprintf(file_out,"%s",out); fclose(file_out);

  int zero=Compute_P(P_ia,NULL,P_glob_a, L_opt,exp1,exp2,
		     f_reg_aa,wi,aa_seq,PMIN,L,Naa);
  printf("%d probability values over %d were < %.2g\n",
	 zero, E_wt->L*Naa, PMIN);
  Copy_P(P_opt_ia, P_ia, L, Naa);
  return(reg_opt);
}

double Compute_score_reg(struct MF_results *res, float reg,
			 float w_max, float *f_aa, float **P_ia,
			 float **exp1, float **exp2, float *P_glob_a,
			 float **f_reg_ia, float **n_msa_ia, float *wi,
			 struct REM *E_wt,
			 float **C_nat, int *i_sec, char *c_sec,
			 char *name_file, FILE *file_summ, char *label,
			 int L, int Naa)
{
  int zero; char out[200], tmp[80];
  int n, nmod=1; if(exp2)nmod=2;
  if(UPDATE_LAMBDA){
    for(n=0; n<nmod; n++){
      if(res->Lambda[n]>0)Lambda_start[n]=res->Lambda[n];
    }
  }
  //Regularize(f_reg_ia, n_msa_ia, f_aa, w_max, reg, L, Naa);
  zero=Optimize_distr(res, P_ia, P_glob_a, exp1, exp2,
		      f_reg_ia, wi, L, Naa, E_wt, C_nat, i_sec, c_sec,
		      name_file,file_summ,label,1);
  sprintf(out,
	  "# reg= %.4g zero=%d lik= %.3g KL_mod+KL_reg= %.3g Lambda= %.3g",
	  reg, zero, res->lik_reg, res->score_reg, res->Lambda[0]);
  if(nmod==2){
    sprintf(tmp, " %.3g", res->Lambda[1]); strcat(out, tmp);
  }
  fprintf(file_summ, "%s\n", out); printf("%s\n", out);
  return(res->score_reg);
}

double Compute_Cv(struct MF_results *SSC_res, float reg, float step, 
		  float w_max, float *f_aa, float **P_SSC_ia,
		  float **exp1, float **exp2, float *P_glob_a,
		  float **f_reg_ia, float **n_msa_ia, float *wi,
		  struct REM *E_wt,
		  float **C_nat, int *i_sec, char *c_sec,
		  char *name_file, FILE *file_summ, char *label,
		  int L, int Naa)
{
  double reg1=reg-step, reg2=reg+step;
  int n, nmod=1; if(exp2)nmod=2;
  if(UPDATE_LAMBDA)
    for(n=0; n<nmod; n++){
      if(SSC_res->Lambda[n]>0)Lambda_start[n]=SSC_res->Lambda[n];
    }
  Regularize(f_reg_ia, n_msa_ia, f_aa, w_max, reg1, L, Naa);
  Optimize_distr(SSC_res,P_SSC_ia,P_glob_a,exp1, exp2,
		 f_reg_ia,wi, L, Naa, E_wt, C_nat, i_sec, c_sec,
		 name_file, file_summ, label, 0);
  double E1=SSC_res->lik_MSA, KL_mod=SSC_res->KL_mod, KL_reg=SSC_res->KL_reg;
  float Lambda[2]; for(n=0; n<nmod; n++)Lambda[n]=SSC_res->Lambda[n];

  Regularize(f_reg_ia, n_msa_ia, f_aa, w_max, reg2, L, Naa);
  Optimize_distr(SSC_res,P_SSC_ia,P_glob_a, exp1,exp2,
		 f_reg_ia,wi, L, Naa, E_wt, C_nat, i_sec, c_sec,
		 name_file, file_summ, label, 0);
  double Cv=(E1-SSC_res->lik_MSA)/(step*2);
  SSC_res->lik_MSA=(SSC_res->lik_MSA+E1)/2;
  SSC_res->KL_mod=(SSC_res->KL_mod+KL_mod)/2;
  SSC_res->KL_reg=(SSC_res->KL_reg+KL_reg)/2;
  for(n=0; n<nmod; n++)SSC_res->Lambda[n]=(SSC_res->Lambda[n]+Lambda[n])/2;

  //char output[200];
  //sprintf(output,
  //	  "# reg_par= %.5g dlik/dreg= %.3g Lambda_start= %.4g zero=%d\n",
  //	  reg, Cv, Lambda_start[0], zero);
  //	  printf("%s", output); //fprintf(file_summ, "%s", output);

  return(Cv);
}

void Remut(float *P_glob_a, float *P_glob_bk, 
	   struct MF_results *k_res, float **P_ia,
	   float **exp1, float **exp2, int is_str,
	   float **f_reg_ia, float **n_msa_ia, float *wi,
	   struct REM *E_wt, float **C_nat, int *i_sec, char *c_sec,
	   short *aa_seq, char *name_file, FILE *file_summ,char *mname,
	   FILE *out, int L, int Naa)
{
  Reg=REG_INI;
  int fail_max=1, nfail=0;
  int a, i;
  for(a=0; a<Naa; a++)P_glob_bk[a]=P_glob_a[a];
  float P_glob_opt[Naa];

  float **P_tmp=Allocate_mat2_f(L, Naa);
  float **exp_str=NULL, **exp_store=NULL;
  if(exp2){exp_str=exp2;}else if(is_str){exp_str=exp1;}
  if(exp_str){
    exp_store=Allocate_mat2_f(L, Naa);
    Copy_P(exp_store, exp_str, L, Naa);
  }

  // Change mutation model
  ini_lik=0;
  int success=0;
  fprintf(file_summ, "# Changing the mutation model, %d trials\n", IT_REMUT);
  for(int iter=0; iter<IT_REMUT; iter++){

    for(i=0; i<2; i++){
      if(k_res->Lambda[i]>0)Lambda_start[i]=k_res->Lambda[i];
    }
    if(k_res->Lambda[0]>0 && (exp2==NULL || k_res->Lambda[1]>0)){
      Reg-=REG_STEP; // Decrease Reg if all Lambda>0
      if(Reg<REG_MIN)Reg=REG_MIN; 
      printf("Reg decreased to %.3g\n", Reg);
    }

    struct MF_results remut_res=*k_res; remut_res.score_reg=-1000000;
    Optimize_distr(&remut_res, P_tmp, P_glob_a, exp1, exp2,
		   f_reg_ia, wi, L, Naa, E_wt, C_nat, i_sec, c_sec,
		   name_file, file_summ, mname, 1);
    printf("Opt: score= %.3g L1= %.3g L2= %.3g\n", remut_res.score_reg,
	   remut_res.Lambda[0], remut_res.Lambda[1]);
    fprintf(out, "Regularization= %.2g score= %.4g Lambda= %.3g %.3g\n",
	    Reg, remut_res.score_reg, remut_res.Lambda[0], remut_res.Lambda[1]);

    if(remut_res.score_reg > k_res->score_reg){
      success=1;
      *k_res=remut_res;
      Copy_P(P_ia, P_tmp, L, Naa);
      if(exp_str)Copy_P(exp_store, exp_str, L, Naa);
      for(a=0; a<Naa; a++)P_glob_opt[a]=P_glob_a[a];
    }else{
      nfail++;
      if(nfail>=fail_max)break;
    }
    //if(Reg==REG_MIN)break;
  }
  if(success){
    for(a=0; a<Naa; a++){P_glob_a[a]=P_glob_opt[a];}
  }else{
    for(a=0; a<Naa; a++){P_glob_a[a]=P_glob_bk[a];}
  }
  Reg=REG_INI;

  fprintf(file_summ, "# New a.a. frequencies ");
  if(success==0)fprintf(file_summ, "do not ");
  fprintf(file_summ, "improve the score\n");

  Empty_matrix_f(P_tmp, len_amm);
  if(exp_str)Copy_P(exp_str, exp_store, L, Naa);
  if(exp_store)Empty_matrix_f(exp_store, len_amm);
}

void Compute_score(struct MF_results *res,
		   float **P_MF_ia, float *Lambda,
		   float **f_reg_ia, float *wi, int L, int Naa, int comp_all)
{
  //Entr_reg=0;
  double log_lik=0, KL_reg=0, entropy=0;
  double Pmin=0.001;
  int Comp_KL=0; if(SCORE_KL || comp_all){Comp_KL=1;}
  for(int i=0; i<L; i++){
    if(wi[i]==0)continue;
    //Entr_reg+=wi[i]*Entropy(f_reg_ia[i], Naa);
    float ll_i=0, ll_mod_i=0, ent_i=0;
    for(int a=0; a<Naa; a++){
      float P=P_MF_ia[i][a]; if(P<Pmin)P=Pmin;
      float l_mod=log(P);
      ll_i+=f_reg_ia[i][a]*l_mod;  // Likelihood of data
      if(Comp_KL){
	ent_i+=P*l_mod;  // Model entropy
	ll_mod_i+=P*lf_ia[i][a]; //log(f_reg_ia[i][a]); // Likelihood of model
      }
    }
    log_lik += wi[i]*ll_i;
    if(Comp_KL){
      entropy-=wi[i]*ent_i;
      KL_reg+=wi[i]*(ent_i - ll_mod_i);
    }
  }
  //res->lik_MSA=Compute_lik(P_MF_ia, f_msa_ia, wi, L, Naa);
  res->lik_reg=log_lik;
  if(Comp_KL){
    res->entropy=entropy;
    res->KL_mod=-Entr_reg-res->lik_reg;
    res->KL_reg=KL_reg;
  }
  if(SCORE_KL){
    res->score=-res->KL_mod-res->KL_reg;
  }else{
    res->score=res->lik_reg;
  }
  if(SCORE_EVO_TMP){
    Print_exchange(P_MF_ia, NULL, NULL, L,
		   NULL, NULL, NULL, 0, 0, NULL, FORMAT,
		   EXCHANGE, MATRIX, 0, 0, 0, NULL, NULL,
		   res, flux_iab, t_ave, NULL);
    res->score+=res->lik_flux;
  }

  res->score_reg=res->score-Reg*(Lambda[0]*Lambda[0]+Lambda[1]*Lambda[1]);
  res->Lambda[0]=Lambda[0];
  res->Lambda[1]=Lambda[1];
}

float Compute_lik(float **P_MF_ia, float **freq_ia, float *wi, int L, int Naa)
{

  double log_lik=0; int i; float eps=0.00000001;
  for(i=0; i<L; i++){
    if(wi[i]==0)continue;
    float ll=0, *f=freq_ia[i];
    for(int a=0; a<Naa; a++){
      float p=P_MF_ia[i][a];
      if(p<eps)p=eps;
      ll+=f[a]*log(p);
    }
    log_lik+=wi[i]*ll;
  }
  return(log_lik);
}

void Remove_extension(char *name1, char *name2){
  strcpy(name1, name2);
  char *s=name1; while(*s!='\0'){s++;} s--;
  char *s_last=s;
  while(s!=name1){if(*s=='.'){*s='\0'; break;} s--;}
  if(s!=name1){while(s!=s_last){*s='\0'; s++;}}
}

void Remove_path(char *name){
  char *s=name, *s_ini=NULL;
  while(*s!='\0'){if(*s=='/'){s_ini=s;} s++;}
  if(s_ini)strcpy(name, s_ini+1);
}

 int Copy_pdb_name(char ***name_pdb, char ***chain_pdb,
		   char *file_pdb, char *chain)
 {
    *name_pdb=malloc(1*sizeof(char *));
    *chain_pdb=malloc(1*sizeof(char *));
    (*name_pdb)[0]=malloc(80*sizeof(char));
    (*chain_pdb)[0]=malloc(80*sizeof(char));
    strcpy((*name_pdb)[0], file_pdb);
    strcpy((*chain_pdb)[0], chain);
    return(1);
 }

int Read_pdb_files(char ***name_pdb, char ***chain_pdb, char *file_pdb_list)
{
  FILE *file_in=fopen(file_pdb_list, "r");
  if(file_in==NULL){
    printf("WARNING, file %s does not exist\n", file_pdb_list);
    return(0);
  }
  char string[200]; int n=0;
  while(fgets(string, sizeof(string), file_in)!=NULL){
    if(string[0]=='#'){continue;} n++;
  }
  fclose(file_in);

  *name_pdb=malloc(n*sizeof(char *));
  *chain_pdb=malloc(n*sizeof(char *));
  file_in=fopen(file_pdb_list, "r");
  n=0;
  while(fgets(string, sizeof(string), file_in)!=NULL){
    if(string[0]=='#'){continue;} 
    (*name_pdb)[n]=malloc(80*sizeof(char));
    (*chain_pdb)[n]=malloc(80*sizeof(char));
    sscanf(string, "%s%s", (*name_pdb)[n], (*chain_pdb)[n]);
    n++;
  }
  fclose(file_in);
  return(n);
}

int Write_DG(float **DG,
	     float **seq_id,
	     struct protein *prot,
	     int N_PDB,
	     float *N_stab_str,
	     char **name_seq,
	     int N_seq,
	     char *buffer,
	     char *name_ali)
{
  char nameout[200], tmp_buff[300]; int i, k;
  sprintf(nameout, "%s.DeltaG", name_ali);
  FILE *file_out=fopen(nameout,"w");
  printf("Writing DeltaG of multiple sequence alignment in %s\n", nameout);
  sprintf(tmp_buff,"Writing DeltaG of multiple sequence alignment %s in %s\n",
	  name_ali, nameout);
  strcat(buffer, tmp_buff);
  
  fprintf(file_out, "# REM, T= %.2g sU1= %.2g sC1= %.2g sC0= %.2g\n",
	  TEMP, sU1, sC1, sC0);
  fprintf(file_out, "#seq ");
  k=2;

  float DG_ave[N_PDB], seq_id_ave[N_PDB], stab_ave=0;
  for(i=0; i<N_PDB; i++){
    DG_ave[i]=0; seq_id_ave[i]=0;
    fprintf(file_out, " %d=DG(%s) %d=seq.id", k, prot[i].name, k+1);
    k+=2;
  }
  fprintf(file_out, "\n");
  
  for(k=0; k<N_seq; k++){
    fprintf(file_out, "%s", name_seq[k]);
    for(i=0; i<N_PDB; i++){
      DG_ave[i]+=DG[i][k]; seq_id_ave[i]+=seq_id[i][k];
      fprintf(file_out, "\t%.3g\t%.3f", DG[i][k], seq_id[i][k]);
    }
    fprintf(file_out, "\t%.2g\n", N_stab_str[k]);
    stab_ave+=N_stab_str[k];
  }
  fprintf(file_out, "#Mean:");
  for(i=0; i<N_PDB; i++){
    fprintf(file_out, "\t%.4g\t%.3f", DG_ave[i]/N_seq, seq_id_ave[i]/N_seq);
  }
  fprintf(file_out, "\t%.2g\n", stab_ave/N_seq);
  fclose(file_out);
  return(0);
}

float Compute_stab_seq(float *DG_seq,
		       float *N_stab_str,
		       float *DG_ave,
		       char *buffer,
		       struct REM E_mut,
		       struct protein prot,
		       float **C_nat, int *i_sec, char *c_sec,
		       short **ali_seq, char **name_seq,
		       float *seq_id, int N_seq)
{
  printf("PDB=%s L=%d T= %.2g sU1= %.2g sC1= %.2g sC0= %.2g\n",
	 prot.name, E_mut.L, TEMP, sU1, sC1, sC0);
  float N_stab=0, DG, x; int nstab=0;
  *DG_ave=0;
  for(int k=0; k<N_seq; k++){
    DG=Compute_DG_overT_contfreq(&E_mut, ali_seq[k], C_nat, i_sec, c_sec, 0);
    if(isnan(DG)){printf("ERROR DG= %.3g\n", DG); exit(8);}
    if(DG>20){x=-1;}
    else if(DG < -20){x=1;}
    else{x=2./(1+exp(DG))-1;}
    if(x > 0){N_stab_str[k]+=x; N_stab+=x; nstab++;}
    printf("seq %s DG= %.2g x= %.1g\n", name_seq[k], DG, x);
    DG_seq[k]=DG;
   (*DG_ave)+=DG;
  }
  (*DG_ave)/=N_seq;
  //if(E_mut.c1U1)free(E_mut.c1U1);

  char tmp_buff[300];
  sprintf(tmp_buff, "%d stable and %d unstable proteins with structure %s "
	  "Mean DG: %.4g Ns= %.3g\n",
	  nstab, N_seq-nstab, prot.name, *DG_ave, N_stab);
  strcat(buffer, tmp_buff);
  printf("End DeltaG computations for %s %d stable seq (%.2g) Mean DG: %.4g\n",
	 prot.name, nstab, N_stab, *DG_ave);
  return(N_stab);
}

void Str_mut_wt(float **exp_ia, float **P_ia, short *aa_seq, int L, int Naa)
{
  for(int i=0; i<L; i++){
    int a_wt=aa_seq[i];
    if(a_wt>=0 && a_wt<Naa){
      float *ei=exp_ia[i], *Pi=P_ia[i];
      double e_ave=0, norm=0;
      for(int a=0; a<Naa; a++){
	if(a!=a_wt){e_ave+=Pi[a]*ei[a]; norm+=Pi[a];}
      }
      if(norm)e_ave/=norm;
      ei[a_wt]=e_ave;
    }
  }
}

void  Sum_DDG(float **exp_ia, float *norm_mod, float **exp_pdb, float w,
	      int *ali_PDB, int len_amm, int Naa)
{
  for(int i=0; i<len_amm; i++){
    int k=ali_PDB[i]; if(k<0)continue;
    norm_mod[i]+=w;
    float *e_a=exp_ia[i], *ep_a=exp_pdb[k];
    for(int a=0; a<Naa; a++){
      e_a[a]+=w*ep_a[a]; ep_a[a]=0;
    }
  }
}

void Normalize_DDG(float **exp_ia, float *norm_mod, int len_amm, int Naa,
		   float Lambda)
{
  for(int k=0; k<len_amm; k++){
    if(norm_mod[k]==0)continue;
    float *e_a=exp_ia[k], fact=1/(norm_mod[k]*Lambda);
    for(int a=0; a<Naa; a++){e_a[a]*=fact;}
  }
}

void Run_tnm(char *file_pdb, char *chain, char *pdbid)
{
  char cmd[400];
  sprintf(cmd, "%s -p1 %s -c1 %s -pred_mut -mut_para %s\n",
	  tnm, file_pdb, chain, tnm_mut_para);
  printf("Running %s\n", cmd);
  system(cmd);
  sprintf(cmd, "mv *.mut_DE.dat %s\n", DIR_TNM); system(cmd);
  sprintf(cmd, "mv *.mut_RMSD.dat %s\n", DIR_TNM); system(cmd);
  sprintf(cmd, "mv *prof_DE.dat %s\n", DIR_TNM); system(cmd);
  sprintf(cmd, "mv *prof_RMSD.dat %s\n", DIR_TNM); system(cmd);
  sprintf(cmd, "rm -f *Inter_Force*\n"); system(cmd);
  sprintf(cmd, "rm -f *Modes.dat\n"); system(cmd);
  sprintf(cmd, "rm -f *Bfact_fit.dat\n"); system(cmd);
  sprintf(cmd, "rm -f *MSF.dat\n"); system(cmd);
  sprintf(cmd, "rm -f *MSF_tors.dat\n"); system(cmd);
}

void Full_prot(float **C_nat_all, float **i_sec_count, float **aa_seq_count,
	       int L_ali, float *norm_c, float w,
	       float **C_nat,
	       int *i_sec, short *aa_seq, int *ali_PDB)
  
{
  // ali_PDB[i] is SEQRES residue aligned with column i of MSA
  (*norm_c)+=w;
  for(int i=0; i<L_ali; i++){
    int k=ali_PDB[i]; if(k<0)continue;
    int a=aa_seq[k]; if(a<0 || a>=Naa)continue;
    aa_seq_count[i][a]+=w;
    a=i_sec[k]; if(a>=0 && a<N_SS){i_sec_count[i][a]+=w;}
    for(int j=0; j<i; j++){
      int l=ali_PDB[j]; if(l<0)continue;
      if(C_nat[k][l]){
	  C_nat_all[i][j]=1;
	  //C_nat_all[i][j]+=w*C_nat[i][j];
      }
    }
  } 
}

void Make_prot(float **C_nat_all, struct res_short **res,
	       int **i_sec, char **c_sec, short **aa_seq,
	       float **i_sec_count, float **aa_seq_count,
	       int len_amm)

 {
   *c_sec=malloc(len_amm*sizeof(char));
   *i_sec=malloc(len_amm*sizeof(int));
   *aa_seq=malloc(len_amm*sizeof(short));
   if(*res)free(*res);
   *res=malloc(len_amm*sizeof(struct res_short));
   int nc[len_amm], i, a; for(i=0; i<len_amm; i++){nc[i]=0;}
   for(i=0; i<len_amm; i++){
     int amax=-1; float cmax=0, *count=i_sec_count[i];
     for(a=0; a<N_SS; a++)if(count[a]>cmax){amax=a; cmax=count[a];}
     (*i_sec)[i]=amax;
     if(amax>=0 && amax<N_SS){(*c_sec)[i]=SEC_EL_cp[amax];}
     else{(*c_sec)[i]='-';}
     (*res)[i].sec_str=(*c_sec)[i];
     amax=-1; cmax=0; count=aa_seq_count[i];
     for(a=0; a<Naa; a++)if(count[a]>cmax){amax=a; cmax=count[a];}
     (*aa_seq)[i]=amax;
     (*res)[i].aa=amax;
     if(amax<0){(*res)[i].seq='-';}
     else{(*res)[i].seq=AMIN_CODE[amax];}
     for(int j=0; j<i; j++){
       if(C_nat_all[i][j]){
	 C_nat_all[j][i]=C_nat_all[i][j];
	 nc[i]++; nc[j]++;
       }
     }
   }
   for(i=0; i<len_amm; i++){(*res)[i].n_cont=nc[i];}
}

int Set_domain(int *ini_dom, int *end_dom, int *ali_PDB, int L_ali){
  // ali_PDB[i] is PDB residue (SEQRES) aligned with MSA column i
  *ini_dom=-1;
  for(int i=0; i<L_ali; i++){
    int k=ali_PDB[i]; if(k<0)continue;
    if(*ini_dom<0){*ini_dom=k;}
    *end_dom=k;
  }
  (*end_dom)++;
  return(*end_dom-*ini_dom);
}

int Align_str_MSA(int *ini_dom, int *end_dom, struct protein *prot,
		  int *ali_PDB, int L_ali)
  // Align PDB residues to ali
{
  // ali_PDB[i] is PDB residue (SEQRES) aligned with MSA column i
  printf("Aligning protein structure with multiple alignmenmt\n");
  short *aa_seq=malloc(L_ali*sizeof(short));
  int *i_sec=malloc(L_ali*sizeof(int));
  char *c_sec=malloc(L_ali*sizeof(char));
  short **contact=malloc(L_ali*sizeof(short *)), cmax=30;
  int L_seq_PDB=prot->length, i;
  int ali_col[L_seq_PDB];
  *ini_dom=-1; *end_dom=-1;
  for(i=0; i<L_seq_PDB; i++){ali_col[i]=-1;}
  for(i=0; i<L_ali; i++){
    contact[i]=malloc(cmax*sizeof(short));
    int k=ali_PDB[i];
    if(k<0){
      aa_seq[i]=-1;
      i_sec[i]=-1;
      c_sec[i]='-';
      contact[i][0]=-1;
    }else{
      if(*ini_dom<0)*ini_dom=i;
      *end_dom=i;
      ali_col[k]=i;
      aa_seq[i]=prot->aa_seq[k];
      i_sec[i]=prot->i_sec[k];
      c_sec[i]=prot->sec_str[k];
      short *Ck=prot->contact[k], *Ci=contact[i], j=0;
      while(*Ck>=0){*Ci=*Ck; Ck++; Ci++; j++;}
      if(j>=cmax){
	printf("ERROR, too many contacts: %d\n", j); exit(8);
      }
      *Ci=-1;
    }
  }
  int *res_index=malloc(prot->L_PDB*sizeof(int));
  for(int i=0; i<prot->L_PDB; i++){
    int k=prot->res_index[i];
    if(k>=0 && k<L_seq_PDB){
      int c=ali_col[k];
      res_index[i]=c;
    }else{
      res_index[i]=-1;
    }
  }

  free(prot->res_index); prot->res_index=res_index;
  free(prot->aa_seq); prot->aa_seq=aa_seq;
  free(prot->i_sec); prot->i_sec=i_sec;
  free(prot->sec_str); prot->sec_str=c_sec;
  for(i=0; i<L_seq_PDB; i++)free(prot->contact[i]); 
  free(prot->contact);
  prot->contact=contact;

  prot->ini_dom=*ini_dom;
  (*end_dom)++;
  prot->L_dom=*end_dom-*ini_dom;
  return(prot->L_dom);
}

void  Align_MSA_PDB(short **ali_seq, char **msa_seq, int *ali_PDB,
		    int N_seq, int L_ali, int end_dom)
{
  for(int n=0; n<N_seq; n++){
    ali_seq[n]=malloc(end_dom*sizeof(short));
    short *as=ali_seq[n], k;
    for(k=0; k<end_dom; k++)as[k]=-1;
    char *aa=msa_seq[n];
    for(int i=0; i<L_ali; i++){
      k=ali_PDB[i];  if(k<0){aa++; continue;}
      if(*aa=='-'||*aa=='.'){aa++; continue;}
      int a=Code_AA(*aa);
      if(a>=0 && a <20){as[k]=a;}
      aa++;
    }
  }
}

float **Compute_t_dist(float *t_ave,
		       char **msa_seq, int *l_stab, int N_seq, int L_ali)
{
  (*t_ave)=0; int ns=0;
  float **t_dist=Allocate_mat2_f(N_seq, N_seq);
  float S0=0.055, S_thr=S0+0.015;
  float t_thr=-log((S_thr-S0)/(1-S0)), l0=log(1-S0);
  for(int n=0; n<N_seq; n++){
    if(l_stab[n]==0)continue;
    for(int m=0; m<n; m++){
      int s=0, l=0;
      for(int i=0; i<L_ali; i++){
	char a=msa_seq[n][i], b=msa_seq[m][i];
	if(a=='-' || b=='-')continue;
	if(a==b){s++;} l++;
      }
      float Seq_id=0; if(l){Seq_id=s/(float)l;}
      float t;
      if(Seq_id <= S_thr){t=t_thr;}
      else{t=-log(Seq_id-S0)+l0; (*t_ave)+=t; ns++;}
      t_dist[n][m]=t;
    }
  }
  if(ns){(*t_ave)/=ns;}
  printf("t= %.3g ns= %d\n", *t_ave, ns);
  return(t_dist);
}
