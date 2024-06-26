#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "coord.h"
#include "REM.h"
#include "meanfield.h"
#include "Print_meanfield.h"
#include "codes.h"
#include "output.h"
#include "gen_code.h"
#include "allocate.h"
#include "optimization.h"
//#include "wag.h"
//#include "jtt.h"
//#include "lg.h"
#include "externals.h"
//

extern int Naa;
extern int PRINT_ALL_EXCH;
extern float LG_matrix[20][20], f_LG[20];
extern float JTT[20][20], fJTT[20];
extern float WAG[20][20], fWAG[20];
extern char AA_WAG[20];

int AVE_FLUX=0; // 1: Flux with average frequencies equal to average flux
                // 0: Rate with average frequencies equal to average rate

char AA_string[]="ACDEFGHIKLMNPQRSTVWY";
char AA_PAML[]="ARNDCQEGHILKMFPSTWYV";
//#define AMIN_CODE "AEQDNLGKSVRTPIMFYCWHX"

float hscale[]={0.1366,-0.0484,0.0325,-0.1233,-0.0345,0.4251,-0.0464,-0.0101,-0.0432,0.4083,0.0363,0.0589,0.0019,0.4172,0.1747,0.4076,0.3167,0.2745,0.2362,0.0549};

static void Change_AA_order(int *iaa, char *AA);
//static void Compute_Q_sel(float **Q_sel, float *P_MF);
static void Print_matrix(float **exchange, int *iaa, float norm,
			 FILE *file_mat, char *AA_string,
			 float *p, int FORMAT);
void Compute_rate_matrix(float **rate_matrix, float *P_aa,
			 float *P_cod, float **Q_cod);
void Compute_exchange_mut_sym(float **S_mut, float **rate_matrix, float *P_aa);
float **Exch_flux(float **P_MF_ia, int L, float **exch_emp, float *f_emp);
// Output: site-specific rate and exchangeability matrix (may be NULL)
// Either compute exch_HB_flux: pointer to it, exchange_glob, f_glob
// Or exch_HB_flux=NULL, exchange_glob=exch_HB_flux, f_glob=NULL


/*static void Compute_exchange_mut_old(float **S_mut,
				     float *mut_par,
				     float tt_ratio,float TWONUC); */
/* static void Update_F(float *F_mut, float *S, char *cod2, float *w,
		     int a, int c, int j, int n, int nt, float tt_ratio); */

double Mean_hydro(float *P, float *hydro, int n);
void Sum_matrix(float *ncont, int **Cnat, int L);
float Corr_coeff(float *xx, float *yy, int n);
float Compute_rate_abs(float *p, float **exchange, int n);
float Compute_rate_Q(float *p, float **rate, int n);
float Compute_rate_sel(float *p, float *p_mut, float **rate, int n);
float Corr_vM(float *p, float **exch, int n);

/******************* Codes ****************************/
int Print_profiles(float **P_MF_ia, char *tag, double DG_ave, float Lambda,
		   float *P_mut_a, short *aa_seq, int L, char *MATRIX,
		   char *nameout, int PRINT_GLOB, float *wi, FILE *out)
{
  int iaa[Naa], i, a;
  Change_AA_order(iaa, AA_PAML); //, AA_string

  // Normalize
  for(i=0; i<L; i++){
    float *P_i=P_MF_ia[i], sum=0;
    for(a=0; a<Naa; a++){sum+=P_i[a];}
    if(fabs(1.0-sum) > 0.001){
      for(a=0; a<Naa; a++){P_i[a]/=sum;}
    }
  }

  char name_file[100]; sprintf(name_file, "%s.%s", nameout, tag);
  char name[200]; sprintf(name, "%s.AA_profiles.txt", name_file);
  FILE *file_out=Output_file(name_file, "AA_profiles", "txt");
  fprintf(out, "Printing site-specific amino acid frequencies in %s\n", name); 
 
  fprintf(file_out, "# Predicted site specific amino acid frequencies ");
  fprintf(file_out, "from model %s, Lambda= %.3f\n", tag, Lambda);
  fprintf(file_out, "# ave(DeltaG/T)= %.2f\n", DG_ave);
  fprintf(file_out, "# Flux model: %s\n", MATRIX);
  // corresponding to Npop= %.3g, Lambda*exp(-DG_ave)

  // Print amino acid names
  fprintf(file_out, "#pos AA_PDB");
  fprintf(file_out, "    %c", AA_PAML[0]); //AA_string[0]
  for(a=1; a<20; a++)fprintf(file_out, "     %c", AA_PAML[a]); //AA_string[a]
  fprintf(file_out, " Entropy\n");

  // Print mutational distribution
  fprintf(file_out, "#P_MUT P=");
  double S=0;
  for(a=0; a<Naa; a++){
    int ia=iaa[a]; float p=P_mut_a[ia]; if(p)S+=p*log(p);
    fprintf(file_out, " %.3f", p);
  }
  fprintf(file_out, "   %.2f\n", -S);

  // Print mean-field distributions
  for(i=0; i<L; i++){
    if(wi[i]==0)continue;
    fprintf(file_out, "%3d ", i);
    if(aa_seq[i]>=0){fprintf(file_out, "%c ", Amin_code(aa_seq[i]));}
    else{fprintf(file_out, "- ");}
    S=0;
    for(a=0; a<Naa; a++){
      int ia=iaa[a]; float p=P_MF_ia[i][ia]; if(p)S+=p*log(p);
      fprintf(file_out, " %.3f", p);
    }
    fprintf(file_out, "   %.2f\n", -S);
  }
  fclose(file_out);

  // Print global profile
  if(PRINT_GLOB){
    sprintf(name, "%s.AA_profile_global.txt", name_file);
    fprintf(out, "Printing global amino acid frequencies in %s\n", name);
    file_out=Output_file(name_file, "AA_profile_global", "txt");
    fprintf(file_out, "#AA P_obs P_MF P_mut diff_MF diff_mut\n");
    double Psum[Naa], Pobs[Naa]; int ll=L;
    for(a=0; a<Naa; a++){Psum[a]=0; Pobs[a]=0;}
    for(i=0; i<L; i++){
      if(aa_seq[i]>=0 && aa_seq[i]<Naa){Pobs[aa_seq[i]]++;}
      else{ll--;}
      for(a=0; a<Naa; a++)Psum[a]+=P_MF_ia[i][a];
    }
    for(a=0; a<Naa; a++){Psum[a]/=L; Pobs[a]/=ll;}
    for(a=0; a<Naa; a++){
      int ia=iaa[a];
      fprintf(file_out, "%c %.2f %.2f %.2f  %5.2f %5.2f\n", 
	      AA_string[a], Pobs[ia]*100, Psum[ia]*100, P_mut_a[ia]*100,
	      (Pobs[ia]-Psum[ia])*100, (Pobs[ia]-P_mut_a[ia])*100);
    }
  }

  return(0);
}

int Print_profile_evo(char *name, double **P_ia, short *aa_seq, int L,
		      double DG_ave, long it_sum)
{
  FILE *file_out=Output_file(name, "AA_profiles_evo", "txt");
  int *iaa=malloc(Naa*sizeof(int)), i, a;
  Change_AA_order(iaa, AA_string);

  fprintf(file_out, "# Evolutionary distribution, %ld substitutions\n",
	  it_sum);
  fprintf(file_out, "# ave(DeltaG/T)= %.2f\n", DG_ave/it_sum);
  fprintf(file_out, "#pos AA");
  for(a=0; a<Naa; a++)fprintf(file_out, "    %c", AA_string[a]);
  fprintf(file_out, " Entropy\n");

  // Print distribution
  for(i=0; i<L; i++){
    fprintf(file_out, "%3d ", i+1);
    if(aa_seq[i]>=0){fprintf(file_out, "%c ", Amin_code(aa_seq[i]));}
    else{fprintf(file_out, "- ");}
    double S=0, norm=0;
    for(a=0; a<Naa; a++)norm+=P_ia[i][a];
    for(a=0; a<Naa; a++){
      int ia=iaa[a]; float p=P_ia[i][ia]/norm; if(p)S+=p*log(p);
      fprintf(file_out, " %.3f", p);
    }
    fprintf(file_out, "   %.2f\n", -S);
  }
  fclose(file_out);
  free(iaa);
  return(0);
}

/* int Print_subst_rate(float **P_MF_ia, float *P_mut_a,
		     float *P_cod, float **Q_cod,
		     float tt_ratio, short *aa_seq,
		     int **Cnat, char *sec_str, int L,
		     char *name_file, float TWONUC)
{
  int iaa[20], i;
  Change_AA_order(iaa, AA_string);
  float **exchange=Allocate_mat2_f(20, 20);
  float ncont[L];
  Sum_matrix(ncont, Cnat, L);
  float entropy[L], Rate[L], hydro_i[L], hydro_ave[L];

  FILE *file_out=Output_file(name_file, "rate_profile", "dat");
  fprintf(file_out, "#SITE AA sec.str. ncont rate entropy hydro ave_hydro\n");
  float P_mut[20], **rate_mut=Allocate_mat2_f(20, 20);
  Compute_rate_matrix(rate_mut, P_mut, P_cod, Q_cod);
  Compute_exchange_mut_sym(exchange, rate_mut, P_mut);
  for(i=0; i<L; i++){
    float *p=P_MF_ia[i];
    Rate[i]=Compute_rate_sel(p, P_mut, exchange, 20);
    entropy[i]=Entropy(P_MF_ia[i], 20);
    hydro_i[i]=hscale[aa_seq[i]];
    hydro_ave[i]=Mean_hydro(P_MF_ia[i], hscale, 20);
    fprintf(file_out, "%3d %c %c %2.0f   %5.3f  %5.2f %5.2f %5.2f\n",
	    i+1, Amin_code(aa_seq[i]), sec_str[i],  ncont[i], Rate[i],
	    entropy[i], hydro_i[i], hydro_ave[i]);
  }
  float r=Corr_coeff(ncont, hydro_i, L);
  fprintf(file_out, "# Corr(ncont, hydro)=     %6.3f %d sites\n", r, L);
  r=Corr_coeff(ncont, hydro_ave, L);
  fprintf(file_out, "# Corr(ncont, ave_hydro)= %6.3f %d sites\n", r, L);
  r=Corr_coeff(ncont, entropy, L);
  fprintf(file_out, "# Corr(ncont, entropy)=   %6.3f %d sites\n", r, L);
  r=Corr_coeff(ncont, Rate, L);
  fprintf(file_out, "# Corr(ncont, rate)=      %6.3f %d sites\n", r, L);
  r=Corr_coeff(entropy, Rate, L);
  fprintf(file_out, "# Corr(rate, entropy)=    %6.3f %d sites\n", r, L);

  fclose(file_out);
  Empty_matrix_f(exchange, 20);
  Empty_matrix_f(rate_mut, 20);
  return(0);
}
*/

int Print_exchange(float **P_MF_ia, char *TAG, struct res_short *res, int L,
		   float *P_glob_a, float *P_cod, float **Q_cod,
		   float tt_ratio, float TWONUC,
		   char *nameout, int FORMAT, char EXCHANGE, char *MATRIX,
		   int PRINT_EXCH, int PRINT_GLOB, int PRINT_MUT,
		   float *wi, FILE *out,
		   struct MF_results *k_res, float ***flux_iab, float t_ave,
		   FILE *file_summ)
{
  int Comp_all; if(file_summ){Comp_all=1;}else{Comp_all=0;}
  //int HB1=1, FL1=1, rate1=1; // if(Comp_all==0) Compute only HB=1 FL=1 rate=1

  // Amino acid order
  int iaa[Naa], i, a, b;
  if(FORMAT==0){Change_AA_order(iaa, AA_string);}
  else{Change_AA_order(iaa, AA_PAML);}
  char name_file[300];
  if(out){sprintf(name_file, "%s.%s", nameout, TAG);}

  // Type of exchangeability model
  char EXCH=EXCHANGE, exc_type[200], namexc[20], type[100];
 //float **exchange_global;
  if(out){
    if(EXCH=='M'){
      strcpy(namexc, "MUT");
      sprintf(type, "exchangeability_HB_MUT");
      sprintf(exc_type,
	      "# Exchangeability: genetic code, tt_ratio=%.2f TWONUC=%.3f\n",
	      tt_ratio, TWONUC); 
      //exchange_global=exchange_mut;
    }else if(EXCH=='E'){
      strcpy(namexc, "EXCHANGE");
      sprintf(type, "exchangeability_HB.%s_EMP", MATRIX);
      sprintf(exc_type, "# Exchangeability: %s\n", MATRIX);
      //exchange_global=exchange_emp;
    }else{
      if(EXCH!='F'){
	printf("WARNING, exchangeability model %c not defined\n", EXCH);
	EXCH='F';printf("Using default %c\n", EXCH);
      }
      strcpy(namexc, "FLUX");
      sprintf(exc_type,"# Exchangeability such that Mean FLUX as %s model\n",
	      MATRIX);
      //exchange_global=exchange_HB_flux;
    }
    fprintf(out, "%s", exc_type);
  }

  FILE *file_mat=NULL; char name[500];
  
  // Empirical model
  float f_emp[Naa],
    **exchange_emp=Empirical_exchangeability(f_emp, MATRIX, Naa);
  Normalize_exchange(exchange_emp, f_emp, Naa, 1);
  float Rate_HB_noflux[L], Rate_noHB_flux[L], Rate_noHB_noflux[L];

  // Average_frequencies sum_i P_i[a]
  float F_all[Naa], sum=0;
  for(a=0; a<20; a++){
    double P=0; for(i=0; i<L; i++)P+=P_MF_ia[i][a]; //wi[i]*
    F_all[a]=P; sum+=P;
  }
  for(a=0; a<20; a++)F_all[a]/=sum;

  float **flux_emp=Allocate_mat2_f(Naa, Naa);
  for(a=1; a<Naa; a++){
    for(b=0; b<a; b++){
      flux_emp[a][b]=f_emp[a]*f_emp[b]*exchange_emp[a][b];
    }
  }
  // Selection matrix sum_i P_i[a]*P_i[b] for computing flux
  float **exchange_flux=NULL;
  if(Comp_all){
    exchange_flux=Exch_flux(P_MF_ia, L, exchange_emp, f_emp);
    Normalize_exchange(exchange_emp, F_all, Naa, 1);
  }

  float norm_rate=1;
  // FILE <>_rate_mut.dat
  // Print exchangeability matrices of mutation model
  // Mutation model
  float R_mut[Naa], Rate_HB_mut[L]; 
  if(out && PRINT_MUT){
    // Compute rates of mutation model
    sprintf(name, "%s.rate_mut.dat", nameout);
    fprintf(out, "Printing Mutational exchangeability ");
    fprintf(out, "matrices in %s\n", name);
    file_mat=fopen(name, "w"); printf("Writing %s\n", name);
    fprintf(file_mat,"#ab E^mut(a,b) ba  E^mut(b,a) E_flux(%s)\n", MATRIX);

    for(i=0; i<L; i++)Rate_HB_mut[i]=0;
    float **rate_mut=Allocate_mat2_f(Naa,Naa), P_mut[Naa];
    Compute_rate_matrix(rate_mut, P_mut, P_cod, Q_cod);

    float **exchange_mut=Allocate_mat2_f(Naa, Naa);
    Compute_exchange_mut_sym(exchange_mut, rate_mut, P_mut);
    Normalize_exchange(exchange_mut, P_mut, Naa, 1);

    // Halpern-Bruno rate obtained with the _mut exchangeability matrix
    Exch_Halpern(Rate_HB_mut, NULL, NULL,
		 P_MF_ia, P_mut, L, exchange_mut, P_mut);
    // Print matrix
    for(a=0; a<Naa; a++){
      R_mut[a]=-rate_mut[a][a];
      for(b=a+1; b<Naa; b++){
	fprintf(file_mat, "%c%c %.2f %c%c %.2f %.2f\n", 
		Amin_code(a), Amin_code(b), rate_mut[a][b]/P_mut[b],
		Amin_code(b), Amin_code(a), rate_mut[b][a]/P_mut[a],
		exchange_flux[b][a]);
      }
    }
    fclose(file_mat);

    Empty_matrix_f(rate_mut, Naa);
    Empty_matrix_f(exchange_mut, Naa);
  }

  //// SSCPE models
  // Print exchangeability matrix in <>_exchangeability_sites_LG_FLUX.txt
  // Print rates in rate_profile.dat


  float Rate_HB_flux[L], **P_fix_HB[L], **exchange_HB_flux=NULL;
  if(PRINT_EXCH || out==NULL){
    exchange_HB_flux=Allocate_mat2_f(Naa,Naa);

    //Exch_Halpern(Rate_HB_noflux, NULL, exchange_HB_flux,
    //		 P_MF_ia, P_mut_a, L, exchange_emp, f_emp);

    for(i=0; i<L; i++){
      float *P_i=P_MF_ia[i];
      float **P_fix_HB_i=Allocate_mat2_f(Naa,Naa);
      P_fix_HB[i]=P_fix_HB_i;
      Get_P_fix(P_fix_HB_i, P_i, F_all); //f_emp P_mut_a
      /* It is important to compute P_fix with HB formula using f_emp instead
	 of P_mut, because now P_mut is biased towards polar amino acids */
      // Sum matrix
      for(a=1; a<Naa; a++){
	for(b=0; b<a; b++){
	  exchange_HB_flux[a][b]+=P_fix_HB_i[a][b]*(P_i[a]*P_i[b]);
	}
      }
    }
    // Flux matrix
    for(a=1; a<Naa; a++){
      for(b=0; b<a; b++){
	exchange_HB_flux[a][b]=flux_emp[a][b]*L/exchange_HB_flux[a][b];
	exchange_HB_flux[b][a]=exchange_HB_flux[a][b];
      }
    }

    // Print global exchangebility matrix
    if(PRINT_GLOB && out){
      sprintf(type, "exchangeability_global.%s_%s", MATRIX, "FLUX");
      file_mat=Output_file(name_file, type, "txt");
      sprintf(name, "%s.%s.txt", name_file, type);
      fprintf(out, "Printing global exchangeability matrix of type ");
      fprintf(out, "%s %c in %s\n", MATRIX, EXCH, name);
      fprintf(file_mat, exc_type);
      if(FORMAT==0){
	fprintf(file_mat, "AA ");
	for(a=0; a<Naa; a++)fprintf(file_mat, "%c\t", AA_string[a]);
      }else{
	for(a=0; a<Naa; a++)fprintf(file_mat, "%c\t", AA_PAML[a]);
      }
      fprintf(file_mat, "\n");
      Print_matrix(exchange_HB_flux, iaa, norm_rate, file_mat,
		   AA_string, P_glob_a, FORMAT);
      if(FORMAT){
	fprintf(file_mat, "\n// end of data. The rest are notes\n");
	for(a=0; a<Naa; a++)fprintf(file_mat, "%c\t", AA_PAML[a]);
	fprintf(file_mat, "\n");
      }
      fclose(file_mat);
    }

    float **exchange_HB_i=Allocate_mat2_f(Naa,Naa), **exchange_HB_noflux=NULL;
    if(Comp_all){exchange_HB_noflux=Allocate_mat2_f(Naa,Naa);}

    // Print site-specific exchange matrices
    if(PRINT_ALL_EXCH && out){
      sprintf(type, "exchangeability_sites.%s_%s", MATRIX, "FLUX");
      file_mat=Output_file(name_file, type, "txt");
      sprintf(name, "%s.%s.txt", name_file, type);
      fprintf(out, "Printing site-specific exchangeability matrices of type ");
      fprintf(out, "%s %c in %s\n", MATRIX, EXCH, name);
      fprintf(file_mat, exc_type);
      if(FORMAT==0){
	fprintf(file_mat, "AA ");
	for(a=0; a<Naa; a++)fprintf(file_mat, "%c\t", AA_string[a]);
      }else{
	for(a=0; a<Naa; a++)fprintf(file_mat, "%c\t", AA_PAML[a]);
      }
      fprintf(file_mat, "\n");
    }

    float **exch_tmp=Allocate_mat2_f(Naa, Naa);
    float lik_flux[8], lik_flux_1[8], lik_flux_log[8];
    int HB_FL[8], FL_FL[8], rate_FL[8], k;
    for(k=0; k<8; k++){
      HB_FL[k]=-1;
      lik_flux[k]=0; lik_flux_1[k]=0; lik_flux_log[k]=0;
    }
    int kk=-1;

    for(i=0; i<L; i++){
      //if(wi[i]==0){continue;}
      /*Exch_Halpern(Rate_HB_flux+i, exchange_HB_i, NULL,
	P_MF_ia+i, P_mut_a, 1, exchange_HB_flux, NULL);
      */
      float **P_fix_HB_i=P_fix_HB[i], *P_i=P_MF_ia[i];
      for(a=1; a<Naa; a++){
	for(b=0; b<a; b++){
	  exchange_HB_i[a][b]=exchange_HB_flux[a][b]*P_fix_HB_i[a][b];
	}
	if(Comp_all){
	  for(b=0; b<a; b++){
	    exchange_HB_noflux[a][b]=exchange_emp[a][b]*P_fix_HB_i[a][b];
	  }
	}
      }
      Empty_matrix_f(P_fix_HB_i, Naa);
      
      if(PRINT_ALL_EXCH && out){
	fprintf(file_mat, "SITE %3d %c\n", i, res[i].seq);
	Print_matrix(exchange_HB_i, iaa, norm_rate, file_mat,
		     AA_string, P_i, FORMAT);
      }
      Rate_HB_flux[i]=Compute_rate_abs(P_i, exchange_HB_i, Naa);
      if(Comp_all){
	Rate_HB_noflux[i]=Compute_rate_abs(P_i, exchange_HB_noflux, Naa);
	Rate_noHB_flux[i]=Compute_rate_abs(P_i, exchange_flux, Naa);
	Rate_noHB_noflux[i]=Compute_rate_abs(P_i, exchange_emp, Naa);
      }

      // Compute flux
      float **flux_ab=flux_iab[i];
      float flux_a[Naa]; for(a=0; a<Naa; a++){flux_a[a]=0;}
      for(a=0; a<Naa; a++){
	for(b=0; b<a; b++){flux_a[a]+=flux_ab[a][b]; flux_a[b]+=flux_ab[a][b];}
      } 
      for(int k=0; k<4; k++){
	int HB, flux, k2=2*k; float **exch;
	if(Comp_all==0 && k!=3)continue;
	if(k==0)     {HB=0; flux=0; exch=exchange_emp;}
	else if(k==1){HB=0; flux=1; exch=exchange_flux;}
	else if(k==2){HB=1; flux=0; exch=exchange_HB_noflux;}
	else if(k==3){HB=1; flux=1; exch=exchange_HB_i;}
	for(a=0; a<Naa; a++)for(b=0; b<a; b++)exch_tmp[a][b]=exch[a][b];
	for(int irate=1; irate>=0; irate--){
	  if(Comp_all==0 && irate==0)continue;
	  if(irate==1){Normalize_exchange(exch_tmp, P_i, Naa, 0);}
	  else{Normalize_exchange(exch_tmp, P_i, Naa, 1);}
	  // irate==0: rates are normalized individually so that rate(i)=1
	  float l1, l_log;
	  lik_flux[k2]+=
	    Flux_exch(&l1, &l_log, flux_ab, flux_a, exch_tmp, P_i, Naa);
	  lik_flux_1[k2]+=l1;
	  lik_flux_log[k2]+=l_log;
	  if(HB_FL[k2]<0){HB_FL[k2]=HB; FL_FL[k2]=flux; rate_FL[k2]=irate;}
	  if(kk<0)kk=k2;
	  k2++;
	}
      }
    }
    Empty_matrix_f(exchange_HB_i, Naa);
    Empty_matrix_f(exch_tmp, Naa);



    if(file_summ){
      for(k=0; k<8; k++){
	lik_flux[k]=Lik_flux(lik_flux[k],lik_flux_1[k],lik_flux_log[k]);
	fprintf(file_summ, "# lik_flux %.4f HB= %d flux= %d rate= %d\n",
		lik_flux[k], HB_FL[k], FL_FL[k], rate_FL[k]);
	if(k==0 || lik_flux[k] > k_res->lik_flux){
	  k_res->HB_opt=HB_FL[k]; k_res->flux_opt=FL_FL[k];
	  k_res->rate_opt=rate_FL[k];
	  k_res->lik_flux=lik_flux[k];
	}
      }
      fprintf(file_summ, "# Highest_lik_flux %.3f HB= %d flux= %d rate= %d\n",
	      k_res->lik_flux, k_res->HB_opt, k_res->flux_opt, k_res->rate_opt);
    }else{
      k=kk;
      lik_flux[k]=Lik_flux(lik_flux[k],lik_flux_1[k],lik_flux_log[k]);
      k_res->HB_opt=HB_FL[k]; k_res->flux_opt=FL_FL[k];
      k_res->rate_opt=rate_FL[k];
      k_res->lik_flux=lik_flux[k];
      goto end_Print;
    }

    if(out==NULL)goto end_Print;

    if(PRINT_ALL_EXCH){
      if(FORMAT){
	fprintf(file_mat, "\n// end of data. The rest are notes\n");
	for(a=0; a<Naa; a++)fprintf(file_mat, "%c\t", AA_PAML[a]);
	fprintf(file_mat, "\n");
      }
      fclose(file_mat);
    }

    sprintf(name, "%s.rate_profile.dat", name_file);
    fprintf(out, "Printing rates of all positions in %s\n", name);
    FILE *file_out=Output_file(name_file, "rate_profile", "dat");
    fprintf(file_out, "# Protein %s L= %d\n", nameout, L);
    fprintf(file_out, "# Local interactions coefficient=%.3f\n", SEC_STR);
    fprintf(file_out, "# Frequencies rate entropy corr(f_a,sum_b E_ab)\n");
		float p[Naa], pp=1./Naa; for(a=0; a<Naa; a++)p[a]=pp;
    fprintf(file_out, "# Equal       %.3f %.3f %6.3f\n",
	    Compute_rate_abs(p, exchange_flux, Naa),
	    log(Naa),
	    Corr_vM(p, exchange_flux, Naa));
    fprintf(file_out, "# Global  %.3f %.3f %6.3f\n",
	    Compute_rate_abs(F_all, exchange_flux, Naa),
	    Entropy(F_all, Naa),
	    Corr_vM(F_all, exchange_flux, Naa));
    fprintf(file_out, "# %s_%s     %.3f %.3f %6.3f\n", MATRIX, namexc,
	    Compute_rate_abs(f_emp, exchange_flux, Naa),
	    Entropy(f_emp, 20),
	    Corr_vM(f_emp, exchange_flux, Naa));
    fprintf(file_out, "#SITE AA sec.str. ncont  entropy ave_hydro");
    fprintf(file_out," rate_noHB_flux rate_noHB_noF rate_HB_flux rate_HB_noflux");
    //fprintf(file_out, " rate_abs_mut rate_abs_emp rate_abs_flux");
    fprintf(file_out, "\n");

    struct res_short *rs=res;
    float entropy[L], hydro_i[L], hydro_ave[L];
    double hydro_tot=0;
    float h_min=1; for(i=0; i<Naa; i++)if(hscale[i]<h_min)h_min=hscale[i];
    for(i=0; i<L; i++){
      float *P_ia=P_MF_ia[i];
      entropy[i]=Entropy(P_ia, Naa);
      hydro_ave[i]=Mean_hydro(P_ia, hscale, Naa);
      hydro_tot+=hydro_ave[i];
      rs=res+i;
      if(rs->aa >=0 && rs->aa < Naa){hydro_i[i]=hscale[rs->aa];}
      else{hydro_i[i]=h_min;}
      if(wi[i]==0){continue;}
      fprintf(file_out, "%d %c %c %d ",i, rs->seq, rs->sec_str, rs->n_cont);
      fprintf(file_out, "%.3f %5.2f ", entropy[i], hydro_ave[i]);
      fprintf(file_out, "%.4f %.4f %.4f %.4f\n",
	      Rate_noHB_flux[i], Rate_noHB_noflux[i],
	      Rate_HB_flux[i], Rate_HB_noflux[i]);
      //Rate_abs_mut[i]=Compute_rate_Q(P_ia, rate_mut, 20);
      //Rate_abs_mut[i]=Compute_rate_sel(P_ia, P_mut, exchange_mut, 20);
      /*Rate_abs_mut[i]=Compute_rate_abs(P_ia, exchange_mut, 20);
	Rate_abs_emp[i]=Compute_rate_abs(P_ia, exchange_emp, 20);
	Rate_abs_flux[i]=Compute_rate_abs(P_ia, exchange_flux, 20);
	fprintf(file_out, "  %.4f %.4f %.4f",
	Rate_abs_mut[i], Rate_abs_emp[i], Rate_abs_flux[i]);*/
    }
    fprintf(file_out, "# %d sites\n", L);
    fprintf(file_out, "# Mean averaged hydrophobicity= %.4f\n", hydro_tot/L);
    float r, ncont[L]; for(i=0; i<L; i++)ncont[i]=res[i].n_cont;
    r=Corr_coeff(ncont, hydro_i, L);
    fprintf(file_out, "# Corr(ncont, hydro)=     %6.3f\n", r);
    r=Corr_coeff(ncont, hydro_ave, L);
    fprintf(file_out, "# Corr(ncont, ave_hydro)= %6.3f\n", r);
    r=Corr_coeff(ncont, entropy, L);
    fprintf(file_out, "# Corr(ncont, entropy)=   %6.3f\n", r);
    r=Corr_coeff(ncont, Rate_HB_flux, L);
    fprintf(file_out, "# Corr(ncont, rate_HB_flux)=  %6.3f\n", r);
    r=Corr_coeff(ncont, Rate_HB_noflux, L);
    fprintf(file_out, "# Corr(ncont, rate_HB_noflux)=  %6.3f\n", r);
    r=Corr_coeff(ncont, Rate_noHB_flux, L);
    fprintf(file_out, "# Corr(ncont, rate_noHB_flux)=  %6.3f\n", r);
    r=Corr_coeff(ncont, Rate_noHB_noflux, L);
    fprintf(file_out, "# Corr(ncont, rate_noHB_noflux)=  %6.3f\n", r);
    //r=Corr_coeff(ncont, Rate_HB_mut, L);
    //fprintf(file_out, "# Corr(ncont, rate_HB_mut)=  %6.3f\n", r);
    r=Corr_coeff(entropy, Rate_HB_flux, L);
    fprintf(file_out, "# Corr(entropy, rate_HB_flux)=%6.3f\n", r);
    r=Corr_coeff(entropy, Rate_HB_noflux, L);
    fprintf(file_out, "# Corr(entropy, rate_HB_noflux)=%6.3f\n", r);
    r=Corr_coeff(entropy, Rate_noHB_flux, L);
    fprintf(file_out, "# Corr(entropy, rate_noHB_flux)=%6.3f\n", r);
    r=Corr_coeff(entropy, Rate_noHB_noflux, L);
    fprintf(file_out, "# Corr(entropy, rate_noHB_noflux)=%6.3f\n", r);
    //r=Corr_coeff(entropy, Rate_HB_mut, L);
    //fprintf(file_out, "# Corr(entropy, rate_HB_mut)=%6.3f\n", r);
    fclose(file_out);
    Empty_matrix_f(exchange_HB_noflux, Naa);
  }

  // FILE exchangeability_glob_%s_FLUX
  // Global exchangeability matrix with average frequencies F_all
  if(PRINT_GLOB){
    sprintf(type, "exchangeability_glob.%s_%s", MATRIX, "FLUX");
    sprintf(name, "%s.%s.txt", name_file, type);
    fprintf(out,
	    "Printing global exchangeability matrix of type %s FLUX in %s\n",
	    MATRIX, name);
    file_mat=Output_file(name_file, type, "txt");
    fprintf(file_mat,"# %s: %s matrix\n", "FLUX", MATRIX);

    float **exchange_ave=Allocate_mat2_f(20,20);
    for(a=0; a<20; a++){
      for(b=0; b<a; b++){
	exchange_ave[a][b]= exchange_emp[a][b]*f_emp[a]*f_emp[b]
	  /(F_all[a]*F_all[b]);
	exchange_ave[b][a]=exchange_ave[a][b];
      }
    }
    Print_matrix(exchange_ave,iaa,norm_rate,file_mat,AA_string,F_all,FORMAT);
    if(exchange_ave)Empty_matrix_f(exchange_ave, 20);
    if(FORMAT){
      fprintf(file_mat,"// end of data. The rest are notes\n");
      for(a=0; a<20; a++)fprintf(file_mat, "%c\t", AA_PAML[a]);
      fprintf(file_mat, "\n");
      for(a=0; a<20; a++)fprintf(file_mat, "%d\t", iaa[a]);
      fprintf(file_mat, "\nL= %d\n", L);
      fprintf(file_mat, "P_mut:\n");
      for(a=0; a<20; a++)fprintf(file_mat, "%.4f\t", P_glob_a[iaa[a]]);
      fprintf(file_mat, "\n");
    }
    fclose(file_mat);
  }

  // FILE <>_AArates.dat
  // Print rates per amino acid in matrix <>_AArates.dat
  if(PRINT_GLOB){
    sprintf(name, "%s.AArates.dat", name_file);
    file_mat=fopen(name, "w"); printf("Writing %s\n", name);
    fprintf(out, "Printing rates for all amino acids in %s\n", name);
    fprintf(file_mat, "#(s)=rescaled to reduce selection\n");
    fprintf(file_mat,"#AA Q^mut(a,a) sum_b(E^emp(s)_ab/19) ");
    fprintf(file_mat,"sum_b(E^emp(s)_ab*f_b) sum_b(E^emp_ab*f_b) hydro\n");
    float R_emp_s[20], Rf_emp_s[20], Rf_emp[20];
    for(a=0; a<20; a++){
      double R=0; for(b=0; b<20; b++)if(b!=a)R+=exchange_flux[a][b];
      R_emp_s[a]=R/19;
      R=0; for(b=0; b<20; b++)if(b!=a)R+=exchange_flux[a][b]*f_emp[b];
      Rf_emp_s[a]=R;
      R=0; for(b=0; b<20; b++)if(b!=a)R+=exchange_emp[a][b]*f_emp[b];
      Rf_emp[a]=R;
      fprintf(file_mat, "%c %.3f %.3f %.3f  %.3f %.2f\n",
	      Amin_code(a), R_mut[a], R_emp_s[a], Rf_emp_s[a],
	      Rf_emp[a], hscale[a]);
    }
    fprintf(file_mat, "# Corr(R_mut,R_emp(s))= %.3f\n",
	    Corr_coeff(R_mut, R_emp_s, 20));
    fprintf(file_mat, "# Corr(R_mut,Rf_emp(s))= %.3f\n",
	    Corr_coeff(R_mut, Rf_emp_s, 20));
    fprintf(file_mat, "# Corr(R_mut,Rf_emp)= %.3f\n",
	    Corr_coeff(R_mut, Rf_emp, 20));
    fprintf(file_mat, "# Corr(R_mut, h)=    %.3f\n",
	    Corr_coeff(R_mut, hscale, 20));
    fprintf(file_mat, "# Corr(R_emp(s), h)=    %.3f\n",
	    Corr_coeff(R_emp_s, hscale, 20));
    fprintf(file_mat, "# Corr(Rf_emp(s), h)=   %.3f\n",
	    Corr_coeff(Rf_emp_s, hscale, 20));
    fprintf(file_mat, "# Corr(Rf_emp, h)=   %.3f\n",
	    Corr_coeff(Rf_emp, hscale, 20));
    fclose(file_mat);
  }
  
 end_Print:
  if(flux_emp)Empty_matrix_f(flux_emp, Naa);
  if(exchange_emp)Empty_matrix_f(exchange_emp, Naa);
  if(exchange_flux)Empty_matrix_f(exchange_flux, Naa);
  if(exchange_HB_flux)Empty_matrix_f(exchange_HB_flux, Naa);
  return(0);
}

float Compute_rate_abs(float *p, float **exchange, int n){
  // F= sum_{a!=b} P_a P_b E_ab
  double F=0; int a, b;
  for(a=1; a<n; a++){
    double Fa=0; for(b=0; b<a; b++)Fa += exchange[a][b]*p[b];
    F+=p[a]*Fa;
  }
  return(2*F);
}

float Compute_rate_Q(float *p, float **rate, int n){
  // R= sum_{a!=b} P_a Q_ab = -sum_a P_a Q_aa
  double R=0;
  for(int a=0; a<n; a++)R+=p[a]*rate[a][a];
  return(-R);
}

float Compute_rate_sel(float *p, float *p_mut, float **exchange, int n)
{
  double R=0; int a, b; float P_mix[n];
  for(a=0; a<n; a++)P_mix[a]=sqrt(p[a]*p_mut[a]);
  for(a=0; a<n; a++){
    double Ra=0;
    for(b=0; b<n; b++){
      if(b!=a)Ra+=exchange[a][b]*P_mix[b];
    }
    R+=P_mix[a]*Ra;
  }
  return(R);
}


float Corr_vM(float *p, float **exch, int n){
  float *e=malloc(n*sizeof(float)); int a, b;
  for(a=0; a<n; a++){
    e[a]=0; for(b=0; b<n; b++)if(b!=a)e[a]+=exch[a][b];
  }
  return(Corr_coeff(p, e, n));
}

void Change_AA_order(int *iaa, char *AA){
  int a; for(a=0; a<20; a++)iaa[a]=Code_AA(AA[a]);
}



/* void Compute_exchange_mut_old(float **S_mut,
			      float *mut_par, float tt_ratio, float TWONUC)
{
  // S(a,b)=sum_{c\in C(a) c'\in C(b)} w_c*Q(c,c')/(W_a*W_b)
   int a, b, c;
  for(a=0; a<20; a++)for(b=0; b<20; b++)S_mut[a][b]=0;
  double *norm=malloc(20*sizeof(double));
  for(a=0; a<20; a++)norm[a]=0;
  float w[64];
  for(c=0; c<64; c++)w[c]=Weight_codon(codon[c], mut_par);
  float r=1; if(TWONUC>0)r=TWONUC; // rate

  for(c=0; c<64; c++){
    if(coded_aa[c]=='*')continue; // Stop codon
    a=Code_AA(coded_aa[c]); norm[a]+=w[c];
    //printf("%2d a=%c b= ", c, Amin_code(a));
    float *Fa=S_mut[a];
    int j1, n1; char cod2[3], *cod=codon[c];
    for(j1=0; j1<3; j1++)cod2[j1]=cod[j1];
    for(j1=0; j1<3; j1++){
      char nj=cod[j1]; int nt=Transition(nj);
      for(n1=0; n1<4; n1++){
	if(Nuc_code(n1)==nj)continue;
	// cod2[j]=n
	float S=r;
	Update_F(Fa,&S,cod2,w,a,c,j1,n1,nt,tt_ratio,codon,coded_aa);
	if(TWONUC<=0)continue;
	int j2, n2, j3, n3;
	for(j2=j1+1; j2<3; j2++){
	  int nt2=Transition(cod[j2]);
	  for(n2=0; n2<4; n2++){
	    if(Nuc_code(n2)==cod[j2])continue;
	    float S2=S*r;
	    Update_F(Fa,&S2,cod2,w,a,c,j2,n2,nt2,tt_ratio,codon,coded_aa);
	    for(j3=j2+1; j3<3; j3++){
	      int nt3=Transition(cod[j3]);
	      for(n3=0; n3<4; n3++){
		if(Nuc_code(n3)==cod[j3])continue;
		float S3=S2*r;
		Update_F(Fa,&S3,cod2,w,a,c,
			 j3,n3,nt3,tt_ratio,codon,coded_aa);
	      } // End n3
	      cod2[j3]=cod[j3];
	    } // End j3
	  } // End n2
	  cod2[j2]=cod[j2];
	} // End j2
      } // End n1
      cod2[j1]=cod[j1];
    } // End j
  } // End c loop
  
  // Impose symmetry
  double sum=0;
  for(a=0; a<20; a++)sum+=norm[a];
  for(a=0; a<20; a++)norm[a]/=sum;
  for(a=0; a<20; a++){
    for(b=0; b<a; b++){
      S_mut[a][b]=(S_mut[a][b]+S_mut[b][a])/(2*norm[a]*norm[b]);
      S_mut[b][a]=S_mut[a][b];
    }
  }

  // Impose normalization Q_aa=-sum_b Q_ab
  //~ for(a=0; a<20; a++){
  //~   double sum=0;
  //~   for(b=0; b<20; b++){
  //~     if(a!=b)sum+=S_mut[a][b]*norm[b];
  //~   }
  //~   S_mut[a][a]=-sum/norm[a];
  //~   }
}*/

/* void Update_F(float *F_mut, float *S, char *cod2, float *w,
	      int a, int c, int j, int n, int nt, float tt_ratio)
{
  cod2[j]=Nuc_code(n);
  int c2=Code_codon(cod2, codon);
  if(coded_aa[c2]=='*')return; // Stop codon
  int b=Code_AA(coded_aa[c2]); if(b==a)return;
  if(n==nt)(*S)*=tt_ratio;
  F_mut[b]+=w[c]*w[c2]*(*S);
  //printf("%c ", Amin_code(b));
} */

/*  
void Compute_Q_sel(float **Q_sel, float *P_MF)
{
  int a, b;
  for(a=0; a<20; a++){
    double sum=0;
    for(b=0; b<20; b++){
      if(b==a)continue;
      if(P_MF[b]>=P_MF[a]){Q_sel[a][b]=1;}
      else{Q_sel[a][b]=P_MF[b]/P_MF[a];}
      sum+=Q_sel[a][b];
    }
    Q_sel[a][a]=-sum;
  }
}
*/

void Sum_matrix(float *ncont, int **Cnat, int L){
  int i, j; float sum;
  for(i=0; i<L; i++){
    sum=0; for(j=0; j<L; j++)sum+=Cnat[i][j];
    ncont[i]=sum;
  }
}

float Average_entropy(float **P_MF_ia, int L, int Naa){
  double entropy=0;
  for(int i=0; i<L; i++)entropy+=Entropy(P_MF_ia[i], Naa);
  return(entropy/=L);
}

float Entropy(float *PP, int n){
  double S=0, norm=0; float *p=PP; int i; 
  for(i=0; i<n; i++){
    if(*p){S+=(*p)*log(*p); norm+=(*p);} p++;
  }
  if(norm==0)return(-1); //Columns without observations
  S = S/norm -log(norm);
  if(S==0){return(S);}
  else{return(-S);}
}

float Corr_coeff(float *xx, float *yy, int n){
  double x1=0, x2=0, y1=0, y2=0, xy=0;
  int i; float *x=xx, *y=yy;
  for(i=0; i<n; i++){
    x1 += *x; x2+= (*x)*(*x);
    y1 += *y; y2+= (*y)*(*y);
    xy += (*x)*(*y); x++; y++;
  }
  float r=(n*xy-y1*x1);
  x2=(n*x2-x1*x1);
  y2=(n*y2-y1*y1);
  if(x2 && y2)r/=sqrt(x2*y2);
  return(r);
}

double Mean_hydro(float *P, float *hydro, int n){
  double ave=0; int a;
  for(a=0; a<n; a++)ave+=hydro[a]*P[a];
  return(ave);
}

void Print_matrix(float **exchange, int *iaa, float norm,
		  FILE *file_mat, char *AA_string, float *p, int FORMAT)
{
  int a,b;
  if(FORMAT==0){
    for(a=0; a<20; a++){
      fprintf(file_mat, "%c", AA_string[a]);
      float *e=exchange[iaa[a]];
      for(b=0; b<20; b++)fprintf(file_mat, "\t%.4f", e[iaa[b]]/norm);
      fprintf(file_mat, "\n");
    }
  }else{
    for(a=1; a<20; a++){
      float *e=exchange[iaa[a]];
      for(b=0; b<a; b++)fprintf(file_mat, "%.4f\t", e[iaa[b]]/norm);
      fprintf(file_mat, "\n");
    }
    fprintf(file_mat, "\n");
    for(a=0; a<20; a++)fprintf(file_mat, "%.4f\t", p[iaa[a]]);
    fprintf(file_mat, "\n");
  }
}

void Compute_exchange_mut_sym(float **S_mut, float **rate_matrix, float *P_aa)
{
  // S(a,b)=sum_{c\in C(a) c'\in C(b)} w_c*Q(c,c')/(W_a*W_b)
  // From flux to exchangeability; symmetrize
  for(int a=0; a<20; a++){
    for(int b=0; b<a; b++){
      S_mut[a][b]=(P_aa[a]*rate_matrix[a][b]+P_aa[b]*rate_matrix[b][a])
	/(2*P_aa[a]*P_aa[b]);
      S_mut[b][a]=S_mut[a][b];
    }
    S_mut[a][a]=rate_matrix[a][a]/P_aa[a];
  }
}
  
void Compute_rate_matrix(float **rate_matrix, float *P_aa,
			 float *P_cod, float **Q_cod)
{
  int a, b, c, d;

  /*float P_cod[64], **Q_cod=Allocate_mat2_f(64, 64); int CpG=1; 
  Compute_P_mut(P_aa, P_cod, Q_cod, mut_par, codon, coded_aa, NULL, NULL);*/
  
  for(a=0; a<20; a++){
    P_aa[a]=0;
    for(b=0; b<20; b++)rate_matrix[a][b]=0;
  }
  for(c=0; c<64; c++){
    if(coded_aa[c]=='*')continue;
    a=Code_AA(coded_aa[c]);
    float p=P_cod[c], *Q=Q_cod[c];
    P_aa[a]+=p;
    for(d=0; d<64; d++){
      if(coded_aa[d]=='*')continue;
      b=Code_AA(coded_aa[d]);
      rate_matrix[a][b]+=p*Q[d];
    }
  } 
  double sum=0; // sum_a!=b P(a)Q(a,b) 
  for(a=0; a<20; a++){
    for(b=0; b<20; b++){
      if(a==b)continue;
      sum+=rate_matrix[a][b];
      rate_matrix[a][b]/=P_aa[a];
    }
  }
  for(a=0; a<20; a++){
    double R=0;
    for(b=0; b<20; b++){
      if(a==b)continue;
      rate_matrix[a][b]/=sum;
      R+=rate_matrix[a][b];
    }
    rate_matrix[a][a]=-R;
  }
}


float **Empirical_exchangeability(float *f_emp, char *MATRIX, int Naa)
{
  float **exch_emp=Allocate_mat2_f(Naa, Naa);
  // Choose empirical exchangeability matrix
  int a, b; float *f, *exch[20];
  if(strncmp(MATRIX, "WAG", 3)==0){
    f=fWAG; for(a=0; a<20; a++)exch[a]=WAG[a];
  }else if(strncmp(MATRIX, "JTT", 3)==0){
    strcpy(MATRIX, "JTT");
    f=fJTT; for(a=0; a<20; a++)exch[a]=JTT[a];
  }else{
    strcpy(MATRIX, "LG");
    f=f_LG; for(a=0; a<20; a++)exch[a]=LG_matrix[a];
  }
  int iwag[20]; Change_AA_order(iwag, AA_WAG);
  for(a=0; a<20; a++){
    int ia=iwag[a]; f_emp[ia]=f[a];
    for(b=0; b<a; b++){
      int ib=iwag[b];
      exch_emp[ia][ib]=exch[a][b];
      exch_emp[ib][ia]=exch[a][b];
    }
  }
  return(exch_emp);
}

float Flux_exch(float *Lik1, float *Lik_log,
		float **flux_ab, float *flux_a,
		float **exch, float *f, int Naa)
{
  /*
    Approximate computation of log-likelihood
    L_i(a,b,t)=f_i(a)exp(t*Q_i(a,b)) (dropping site i for simplicity)
    ~ f(a)(1+tQ(a,a))*delta(a,b) + tf(a)Q(a,b) (1-delta(a,b))
    log(L) ~ delta(a,b) * tf(a)Q(a,a)  +
             (1-delta(a,b))*[log(t)+log(f_a)+log(f_b)+log(E_ab)]
    = t*Lik1 +log(t)*Lik_log *log(t) + Lik
    Lik1    = delta(a,b) * f(a)*Q(a,a)
    Lik     =  (1-delta(a,b))*[log(f_a)+log(f_b)+log(E_ab)]
    Lik_log =  (1-delta(a,b))
    The exchangebility matrix is assumed to be symmetric
    and normalized such that sum_ab E_ab*f_b=1, E_aa*f_a=-sum_b E_ab*f_b
    The amino acid indexes are the same as used for the observed flux
    int iwag[20]; Change_AA_order(iwag, AA_WAG);
    Choose the value of t that maximizes the likelihood
  */
  double Lik=0; *Lik1=0; *Lik_log=0;
  for(int a=0; a<Naa; a++){
    for(int b=0; b<a; b++){
      Lik+=flux_ab[a][b]*log(exch[a][b]);
    }
    Lik+=flux_a[a]*log(f[a]);
    (*Lik_log)+=(flux_a[a]-flux_ab[a][a]);
    (*Lik1)+=flux_ab[a][a]*exch[a][a]*f[a]; //*f[a]
  }
  return(Lik);
}

float Lik_flux(float lik, float lik_1, float lik_log){
  /*
    Approximate computation of log-likelihood
    log(L) ~  t*Lik1 +log(t)*Lik_log *log(t) + Lik
    Choose the value of t that maximizes the likelihood
  */
  float tmax=-lik_log/lik_1;
  float ymax=lik_1*tmax+lik_log*log(tmax);
  ymax+=lik;
  printf("Approximate likelihood: %.4g t= %.2g l=%.2g l1=%.2g l_log= %.2g\n",
	 ymax, tmax, lik, lik_1, lik_log);
  //exit(8);
  return(ymax);
}

float Flux_exch_MSA(float **f_msa, int L, float sum_msa,
		    float exch_emp[20][20], float fe[20])
{
  // Symmetrize
  int a, b, i;
  float exch[20][20];
  for(a=0; a<20; a++){
    for(b=0; b<a; b++){
      exch[a][b]=exch_emp[a][b];
      exch[b][a]=exch_emp[a][b];
    }
  }

  // Normalize exchangebility matrix
  float rate=0;
  for(a=0; a<20; a++){
    float Q_aa=0;
    for(b=0; b<20; b++){if(b!=a)Q_aa+=exch[a][b]*fe[b];}
    exch[a][a]=-Q_aa/fe[a];
    rate+=fe[a]*Q_aa;
  }
  for(a=0; a<20; a++){
    for(b=0; b<20; b++){exch[a][b]/=rate;}
  }

  // Choose empirical exchangeability matrix
  double Lik=0;
  int iwag[20];
  Change_AA_order(iwag, AA_WAG);
  for(i=0; i<L; i++){
    float *f=f_msa[i];
    for(a=0; a<20; a++){
      int ia=iwag[a]; if(f[ia]==0)continue;
      for(b=0; b<a; b++){
	int ib=iwag[b]; if(f[ib]==0)continue;
	Lik+=f[ia]*f[ib]*exch[a][b];
      }
    }
  }
  return(Lik/(sum_msa*sum_msa));
}

float Normalize_exchange(float **exch, float *f, int n, int norm_rate)
{
  double sum=0; int a, b;
  float Ra[n]; for(a=0; a<n; a++){Ra[a]=0;}
  for(a=1; a<n; a++){
    for(b=0; b<a; b++){
      Ra[a]+=exch[a][b]*f[b]; Ra[b]+=exch[a][b]*f[a];
    }
    sum+=f[a]*Ra[a];
  }
  sum*=2;
  for(a=0; a<20; a++){
    exch[a][a]=-Ra[a]/f[a];
    if(norm_rate){exch[a][a]/=sum;}
    for(b=0; b<a; b++){
      if(norm_rate){exch[a][b]/=sum;}
      exch[b][a]=exch[a][b];
    }
  }
  return(sum);
}


float **Exch_flux(float **P_MF_ia, int L, float **exch_emp, float *f_emp)
{
  // Impose that average flux (F) is as in empirical model
  float **exchange_flux=Allocate_mat2_f(20,20);
  int a, b, i;
  for(a=0; a<20; a++){
    for(b=0; b<a; b++){
      double P=0; for(i=0; i<L; i++)P+=P_MF_ia[i][a]*P_MF_ia[i][b];
      exchange_flux[a][b]= exch_emp[a][b]*f_emp[a]*f_emp[b]*L/P;
      exchange_flux[b][a]=exchange_flux[a][b];
    }
  }
  return(exchange_flux);
}

void Get_P_fix(float **P_fix_HB_i, float *Pi, float *P_mut)
{
  float P_MIN=0.001, log_MIN=log(P_MIN);
  int a,b;

  float P_sel[Naa], log_P_sel[Naa];
  for(a=0; a<Naa; a++){
    P_sel[a]=Pi[a]/P_mut[a];
    if(P_sel[a]>P_MIN){log_P_sel[a]=log(P_sel[a]);}
    else{log_P_sel[a]=log_MIN;}
  }
  for(a=1; a<20; a++){
    // Fixation probability from a to b
    for(b=0; b<a; b++){
      if(P_sel[a]!=P_sel[b]){
	P_fix_HB_i[a][b]=(log_P_sel[a]-log_P_sel[b])/(P_sel[a]-P_sel[b]);
      }else{
	P_fix_HB_i[a][b]=1;
      }
    }
    //P_fix_HB_i[b][a]=P_fix_HB_i[a][b];
  }
}


void Exch_Halpern(float *rate_HB, float **exchange_HB,
		  float **exchange_HB_flux, // output
		  float **P_MF_ia, float *P_mut, int L,  // Input
		  float **exchange_emp, float *f_emp)
{
  float P_MIN=0.001, log_MIN=log(P_MIN);
  int i,a,b;
  if(exchange_HB_flux){
    for(a=0; a<20; a++)for(b=0; b<20; b++)exchange_HB_flux[a][b]=0;
  }
  for(i=0; i<L; i++){
    if(rate_HB)rate_HB[i]=0;
    float *Pi=P_MF_ia[i];
    float P_sel[20], log_P_sel[20];
    for(a=0; a<20; a++){
      P_sel[a]=Pi[a]/P_mut[a];
      if(P_sel[a]>P_MIN){log_P_sel[a]=log(P_sel[a]);}
      else{log_P_sel[a]=log_MIN;}
    }
    for(a=1; a<20; a++){
      double R_a=0, Fix;
      // Fixation probability from a to b
      for(b=0; b<a; b++){
	if(P_sel[a]!=P_sel[b]){
	  Fix=(log_P_sel[a]-log_P_sel[b])/(P_sel[a]-P_sel[b]);
	}else{
	  Fix=1;
	}
	double e_fix=exchange_emp[a][b]*Fix;
	R_a+=Pi[b]*e_fix;
	if(exchange_HB){
	  //float **exchange_HB_i=exchange_HB[i];
	  exchange_HB[a][b]=e_fix;
	  exchange_HB[b][a]=e_fix;
	}
	if(exchange_HB_flux)
	  exchange_HB_flux[a][b]+=Pi[a]*Pi[b]*Fix;
      }
      if(rate_HB){rate_HB[i]+=Pi[a]*R_a;}
    }
    if(rate_HB)rate_HB[i]*=2;
  }
  if(exchange_HB_flux){
    for(a=0; a<20; a++){
      for(b=0; b<a; b++){
	exchange_HB_flux[a][b]=
	  L*f_emp[a]*f_emp[b]*exchange_emp[a][b]/exchange_HB_flux[a][b];
      }
    }
  }
  return;
}
