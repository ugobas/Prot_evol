int Print_profiles(float **P_MF_ia, char *tag, double DG_ave, float Lambda,
		   float *P_mut_a, short *aa_seq, int L, char *MATRIX,
		   char *name_file, int PRINT, float *wi, FILE *out);
float Flux_exch(float *Lik1, float *Lik_log,
		float **flux_ab, float *flux_a,
		float **exch, float *fe, int Naa);
float Lik_flux(float lik, float lik_1, float lik_log);
float Flux_exch_MSA(float **f_msa, int L, float sum_msa,
		    float exch[20][20],float fe[20]);

int Print_exchange(float **P_MF_ia, char *TAG,
		   struct res_short *res, int L,
		   float *P_mut_a, float *P_cod, float **Q_cod,
		   float tt_ratio, float TWONUC,
		   char *nameout, int FORMAT, char EXCHANGE, char *MATRIX,
		   int PRINT_EXCH, int PRINT_GLOB, int PRINT_MUT,
		   float *wi, FILE *out,
		   struct MF_results *k_res, float ***flux_iab, float t_ave,
		   FILE *file_summ);
// When called not for printing, out should be set to NULL

/*int Print_subst_rate(float **P_MF_ia, float *P_mut_a,
		     float *P_cod, float **Q_cod,
		     float tt_ratio, short *aa_seq,
		     int **Cnat, char *sec_str, int L,
		     char *name_file, float TWONUC);*/
int Print_profile_evo(char *name, double **P_ia, short *aa_seq, int L,
		      double DG_ave, long it_sum);
float Entropy(float *P, int n);
float **Empirical_exchangeability(float *f_emp, char *MATRIX, int Naa);
float Normalize_exchange(float **exch, float *f, int n, int norm_rate);
float Flux_lik(float **f_msa, int L, float sum_msa,
	       float exch_emp[20][20], float fe[20]);
void Get_P_fix(float **P_fix_HB_i, float *Pi, float *P_mut);
void Exch_Halpern(float *rate_HB, float **exchange_HB_i,
		  float **exchange_HB_flux, // output
		  float **P_MF_ia, float *P_mut, int L, // Input
		  float **exchange_glob, float *f_glob);
