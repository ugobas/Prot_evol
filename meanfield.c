#include "REM.h"
//#include "energy_BKV.h"  // Contact energy parameters
#include "protein3.h"
#include "allocate.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "meanfield.h"
#include "externals.h"
extern float Entropy(float *P, int n);
extern void Compute_score(struct MF_results *res,
			  float **P_MF_ia, float *Lambda,
			  float **f_reg_ia, float *wi, int L, int Naa,
			  int comp_all);
static float Average(float *u, float *p, int n);
void Average_P(float **P1, float **P2, int L, int Naa);
extern int N_SS;

#define IT_MAX 50
#define EPS 0.0000001
//float **Econt_T=NULL, T_ratio=1; // Not used anymore

int L_old=0;
float *Cont2_L=NULL;
float **P_new_ia=NULL, **P_opt_ia=NULL,
  **u1_ja=NULL, **u2_ja=NULL, **u3_ja=NULL, **c1U1_ja=NULL, **c1U2_ja=NULL;
// Not used anymore
float **u1U_ja=NULL, **Uc1U1_ja=NULL, **Uc1U2_ja=NULL, **U2c1U1_ja=NULL;

int INI_MF=0;
extern int SCORE_KL;
extern int DTAIL_MAX;

void Average_U(float **u1_ja, float **u2_ja, float **u3_ja,
	       float **P_MF_ia, int L, int Naa);
void  Sum_ene(float **c1U1_ja, float **c1U2_ja,
	      float **Uc1U1_ja, float **Uc1U2_ja,
	      float **U2c1U1_ja, float **u1U_ja,
	      float **u1_ja, float **u2_ja,
	      float **P_MF_ia, float *Cont_L, 
	      int L, int Naa);

double Compute_DG_overT_REM(struct REM *E, float *T_eff, float *factor,
			    float **u1_ja, float **u2_ja, float **u3_ja,
			    float **c1U1_ja,float **c1U2_ja,
			    float **P_MF_ia,
			    float **C_nat, int *i_sec, char *c_sec,
			    float *Cont_L, int L, int Naa);
static void Update_P(float **P_new_ia, float **P_old_ia, float *P_mut_a,
		     float **u1_ja, float **u2_ja, float **u3_ja,
		     float **c1U1_ja, float **c1U2_ja,
		     float **Uc1U1_ja, float **Uc1U2_ja,
		     float **U2c1U1_ja, float **u1U_ja,
		     struct REM *E, float **C_nat, int *i_sec, char *c_sec,
		     float *Cont_L, float *T_eff, float *factor,
		     float Lambda, int L, int Naa);
static void Update_P_nat(float **P_new_ia, float **P_old_ia, float *P_mut_a,
			 float **u1_ja, struct REM *E,
			 float **C_nat, int *i_sec, char *c_sec,
			 float Lambda, int L, int Naa);

static float Test_convergence(float **P1_ia, float **P2_ia, int L, int Naa);
void Boltzmann_distr(float *P, float beta, float *hydro, int N, float equal);
void Initialize_P(float **P_MF_ia, float **C_nat, int *i_sec, char *c_sec,
		  float *Cont_L, float *hydro, float Lambda, int L, int Naa);
static void Allocate_mean_field(float ***P_new_ia, float ***P_opt_ia,
				float ***u1_ja, float ***u2_ja, float ***u3_ja,
				float ***c1U1_ja, float ***c1U2_ja,
				float ***Uc1U1_ja, float ***Uc1U2_ja,
				float ***U2c1U1_ja, float ***u1U_ja,
				int REM, int L, int Naa);

static void Empty_mean_field(float ***P_new_ia, float ***P_opt_ia,
			     float ***u1_ja, float ***u2_ja, float ***u3_ja,
			     float ***c1U1_ja, float ***c1U2_ja, 
			     float ***Uc1U1_ja, float ***Uc1U2_ja,
			     float ***U2c1U1_ja, float ***u1U_ja,
			     int L);


/*****************************************************************************/

int meanfield(float **P_MF_ia, struct MF_results *res,
	      //double *DG_ave, float *Tfreeze, float *Kul_Leib
	      float *P_mut_a, float **P_ini_ia,
	      float **C_nat, int *i_sec, char *c_sec,
	      float Lambda, struct REM E_wt, int L, int Naa, int ALL) //Npop
{
  printf("# Mean field, Lambda=%.3f:\n", Lambda);
  int SCORE_DEF=1; // 1=Minimize KLD+Lambda*DG 0=Minimize convergence

  // Allocate
  Allocate_mean_field(&P_new_ia, &P_opt_ia,
		      &u1_ja, &u2_ja, &u3_ja, &c1U1_ja,&c1U2_ja,
		      &Uc1U1_ja, &Uc1U2_ja,&U2c1U1_ja,&u1U_ja,
		      E_wt.REM, L, Naa);

  struct REM E; E.c1U1=NULL;
  Initialize_E_REM(&E, E_wt.L, E_wt.REM, E_wt.T, E_wt.S_C, E_wt.S_U,
		   FILE_STR);

  // Initialize distribution
  if(P_ini_ia){
    Copy_P(P_MF_ia, P_ini_ia, L, Naa);
  }else{
    Initialize_P(P_MF_ia, C_nat, i_sec, c_sec, Cont_L, hydro, Lambda, L, Naa);
  }
  Copy_P(P_opt_ia, P_MF_ia, L, Naa);

  // Iteration
  int convergence=0, iter, it_opt=0, i, nfail=1, nfail_max=4;
  float conv=0, conv_opt=0, score, score_min=10000000, free_e, KLD=0, KLD_opt=0,
    f_opt=0, DG_prev=0, DG_opt=0, Tf_opt=0, factor=1, T_eff;

  printf("#iteration Kul-Leib DeltaG/L (KL+Lam*DG)/L factor convergence\n");
  //double T_eff, Tf_opt, free_e; //T_old //Tfreeze
  //if(ALL){printf("KL+Lambda*DG");}else{printf("KL+Lambda*E_nat");}
  //printf(" (Lambda=%.2f T=%.2f SC=%.2f SU=%.2f)",
  //	 Lambda, E_wt.T, E_wt.S_C, E_wt.S_U);

  for(iter=0; iter<IT_MAX; iter++){

    //if(iter==0)T_old=T_eff;
    //T_new=(T_old+T_eff)/2; T_old=T_eff;
    //float T_new=T_eff;

    KLD=0; for(i=0; i<L; i++)KLD+=KL(P_MF_ia[i], P_mut_a, Naa);
    if(ALL){
      Update_P(P_new_ia, P_MF_ia, P_mut_a,
	       u1_ja, u2_ja, u3_ja, c1U1_ja, c1U2_ja,
	       Uc1U1_ja, Uc1U2_ja, U2c1U1_ja, u1U_ja,
	       &E, C_nat, i_sec, c_sec, Cont_L,
	       &T_eff, &factor, Lambda, L, Naa);
    }else{
      Update_P_nat(P_new_ia, P_MF_ia, P_mut_a, u1_ja, &E,
		   C_nat, i_sec, c_sec, Lambda, L, Naa);
      E.DeltaG=E.E_nat+E.S_U;
    }
    conv=Test_convergence(P_new_ia, P_MF_ia, L, Naa);
    free_e = KLD+Lambda*E.DeltaG;
    printf("%2d\t%.3g\t%.3g\t%.4g\t%.2f\t%.2g\n", //\t%.2f\t%.2f
	   iter, KLD/L, E.DeltaG/L, free_e/L, factor, sqrt(conv));

    if(SCORE_DEF){score= free_e;}
    else{score=conv;}
    if(score < score_min){
      nfail=0; it_opt=iter;
      score_min = score;
      DG_opt=E.DeltaG; KLD_opt=KLD; Tf_opt=E.Tf;
      f_opt=factor; conv_opt=conv;
      Copy_P(P_opt_ia, P_MF_ia, L, Naa);
    }else{
      nfail++; if(nfail>=nfail_max)break;
    }
    if((conv<EPS)||(fabs(E.DeltaG-DG_prev)<0.0001)){convergence=1; break;}

    Copy_P(P_MF_ia, P_new_ia, L, Naa);
    DG_prev=E.DeltaG;
  }
  if(convergence==0){
    printf("WARNING, mean field computation did not converge in %d", iter+1);
    printf(" steps.\nUsing distribution with minimum ");
    if(SCORE_DEF){
      printf("Kullback-Leiber div. + Lambda*DG= %.3g (/L) ",score_min/L);
    }else{
      printf("convergence error %.3g ", sqrt(score_min));
    }
    printf("DG= %.3f\n", DG_opt);
  }

  //Average_P(P_opt_ia, P_opt_ia, L, Naa);
   //KLD=KLD_opt; DG=DG_opt; E.Tf=Tf_opt;
  Copy_P(P_MF_ia, P_opt_ia, L, Naa);
  printf("%2d\t%.3g\t%.3g\t%.4g\t%.2f\t%.2g\n", //\t%.2f\t%.2f
	 it_opt, KLD_opt/L, DG_opt/L, (KLD_opt+Lambda*DG_opt)/L,
	 f_opt, sqrt(conv_opt));

  res->KL_mut=KLD_opt/L;
  res->DG=DG_opt;
  res->Tf=Tf_opt;
  res->Lambda[0]=Lambda;
  res->h=Hydro_ave(P_MF_ia, hydro, NULL, L);

  Empty_mean_field(&P_new_ia,&P_opt_ia,&u1_ja,&u2_ja,&u3_ja,&c1U1_ja,&c1U2_ja,
		   &Uc1U1_ja, &Uc1U2_ja, &U2c1U1_ja, &u1U_ja, L);
  Empty_E_REM(&E);
  //exit(8);
  return(convergence);
}


void Average_U(float **u1_ja, float **u2_ja, float **u3_ja,
	       float **P_MF_ia, int L, int Naa)
{
  int i, a, b;
  for(i=0; i<L; i++){
    float *p=P_MF_ia[i];
    float *u1j=u1_ja[i];
    for(a=0; a<Naa; a++){
      float *U=Econt_T[a];
      u1j[a]=0; double U2=0, U3=0; 
      for(b=0; b<Naa; b++){
	if(isnan(p[b])){
	  printf("ERROR in Average_U i= %d b= %d p[b]= %.2g\n",i,b,p[b]);
	  exit(8);
	}
	float x=U[b]*p[b]; u1j[a]+=x;
	if(u2_ja){x*=U[b]; U2+=x;}
	if(u3_ja){x*=U[b]; U3+=x;}
      }
      if(u2_ja){u2_ja[i][a]=U2;}
      if(u3_ja){u3_ja[i][a]=U3;}
    }
  }
}

void  Sum_ene(float **c1U1_ja, float **c1U2_ja,
	      float **Uc1U1_ja, float **Uc1U2_ja,
	      float **U2c1U1_ja, float **u1U_ja, 
	      float **u1_ja, float **u2_ja,
	      float **P_MF_ia, float *Cont_L,
	      int L, int Naa)
{

  int i, j, a;
  float *c1U2=NULL;
  for(i=0; i<L; i++){
    float *c1U1=c1U1_ja[i];
    for(a=0; a<Naa; a++)c1U1[a]=0;
    if(c1U2_ja){
      c1U2=c1U2_ja[i]; for(a=0; a<Naa; a++)c1U2[a]=0;
    }
    for(j=0; j<L; j++){
      int l=abs(i-j); if(l<IJ_MIN)continue;
      float c=Cont_L[l];
      for(a=0; a<Naa; a++){
	c1U1[a]+= c*u1_ja[j][a];
	if(c1U2)c1U2[a]+= c*u2_ja[j][a];
      }
    }
  }

  float *Uc1U1=NULL, *Uc1U2=NULL, *U2c1U1=NULL, *u1U=NULL;
  for(i=0; i<L; i++){
    float *c1U1=c1U1_ja[i]; 
    if(Uc1U1_ja){
    Uc1U2=Uc1U2_ja[i];U2c1U1=U2c1U1_ja[i];Uc1U1=Uc1U1_ja[i],u1U=u1U_ja[i];
    }
    float *p=P_MF_ia[i];
    for(a=0; a<Naa; a++){
      if(Uc1U1_ja){Uc1U1[a]=0; u1U[a]=0;Uc1U2[a]=0; U2c1U1[a]=0;}
      float *U=Econt_T[a];
      for(int b=0; b<Naa; b++){
	float pU=p[b]*U[b];
	if(Uc1U1_ja){
	  Uc1U1[a] += pU*c1U1[b];
	  u1U[a]   += pU*u1_ja[i][b];
	  Uc1U2[a] += pU*c1U2[b];
	  U2c1U1[a]+= pU*c1U1[b]*U[b];
	}
      }
    }
  }
}

double Compute_DG_overT_REM(struct REM *E, float *T_eff, float *factor,
			    float **u1_ja, float **u2_ja, float **u3_ja,
			    float **c1U1_ja, float **c1U2_ja,
			    float **P_MF_ia,
			    float **C_nat, int *i_sec, char *c_sec,
			    float *Cont_L, int L, int Naa)
{
  Set_to_zero(E);
  //double E12_sum=0; 
  int i, j, a;
  for(i=0; i<L; i++)E->c1U1[i]=0;
  for(i=0; i<L; i++){
    float *p=P_MF_ia[i];
    if(SEC_STR && i_sec[i]>=0 && i_sec[i]<N_SS){
      E->E_loc+=Average(E_loc_over_T[i_sec[i]], p, Naa);
    }
    //double u12_c2[Naa];
    //for(a=0; a<Naa; a++)u12_c2[a]=0;
    for(j=0; j<=(i-IJ_MIN); j++){
      int l=abs(i-j);
      float U1ij=Average(u1_ja[j], p, Naa);
      if(C_nat[i][j])E->E_nat+=C_nat[i][j]*U1ij;
      if(E->REM==0)continue;
      float c=Cont_L[l], cU=c*U1ij;
      E->c1U1[i]+=cU;
      E->c1U1[j]+=cU;
      if(E->REM>=2){
	E->E2cont2+=C2_L[l]*cU;
	float U2ij=Average(u2_ja[j], p, Naa);
	E->E2cont1+=U2ij*Cont2_L[l];
	if(E->REM==3){
	  E->E3cont3 += C3_L[l]*cU;
	  E->E3cont2 += c*U2ij*Cont_2L[l];
	  E->E3cont1 += c*Average(u3_ja[j], p, Naa)*Cont_3L[l];
	}
      }
      /*float c2=c*c;
	for(a=0; a<Naa; a++)u12_c2[a]+=c2*u1j[a]*u1j[a];
	if(j>i){    
	double U2ij=0; float *u2j=u2_ja[j];
	for(a=0; a<Naa; a++)U2ij+=p[a]*u2j[a];
	  E->E22 += C221_ij[l]*U2ij;*/
	  /*if(E->REM==3){
	    double U3ij=0; float *u3j=u3_ja[j];
	    for(a=0; a<Naa; a++)U3ij+=p[a]*u3j[a];
	    E->E321+= C321_ij[l]*U3ij;
	    E->e342+= C342_ij[l]*U2ij;
	    E->c332U2[i]+= C332_ij[i][j]*U2ij;
	    E->c332U2[j]+= C332_ij[j][i]*U2ij;
	    }*/
    }
    if(E->REM==0)continue;


    /*if(E->REM > 1){
      E12_sum += E1ii*E1ii;
      double E11=0; float *c1U1=c1U1_ja[i];
      for(a=0; a<Naa; a++)E11+=p[a]*(c1U1[a]*c1U1[a]-u12_c2[a]);
      E->E23  += C232_i[i]*E11;
      if(E->REM==3){
	double E12=0; float *c1U2=c1U2_ja[i];
	for(a=0; a<Naa; a++)E12+=p[a]*c1U1[a]*c1U2[a];
	E332 += M332_i[i]*E12;
	}
	}*/    
  } //end i
  for(i=0; i<L; i++)E->E1+=E->c1U1[i];
  E->E1/=2;

  /*if(E->REM > 1){E->E24=(C242)*E->E1*E->E1-E12_sum;}
    else{E->E24=0;}*/
  E->G_misf=G_misfold(E);
  if(isnan(E->G_misf)){
    printf("WARNING, G_misf= %.2g Tf= %.2g E1= %.2g E2= %.2g ",
	   E->G_misf, E->Tf, E->E1, E->E2);
    printf("E3= %.2g E3cont3= %.2g E3site2= %.2g\n",
	   E->E3, E->E3cont3, E->E3site2);
    for(i=0; i<10; i++){
      printf("i= %d P= ", i);
      for(a=0; a<Naa; a++)printf(" %.2f", P_MF_ia[i][a]);
      printf("\n");
    }
  }
  float DeltaG_overT=DeltaG((float)E->E_nat, E->G_misf, E->S_U);
  if(E->REM==1)E->Tf=0;
  if(E->Tf>1){*T_eff=E->Tf;}else{*T_eff=1;}
  double DG=E->G_misf+E->S_U;
  if(DG < -30){*factor=1;}
  else if(DG > 30){*factor=0;}
  else{*factor=1./(1+exp(DG));}
  return(DeltaG_overT);
}

void Update_P_nat(float **P_new_ia, float **P_old_ia, float *P_mut_a,
		  float **u1_ja, struct REM *E,
		  float **C_nat, int *i_sec, char *c_sec,
		  float Lambda, int L, int Naa)
{
  int i, j, a; double E_nat=0;
  //float x=exp(DG), Lambda=(x/(1+x))*Npop;
  //printf("lambda= %.2g\n", Lambda);
  for(i=0; i<L; i++){
    for(a=0; a<Naa; a++)
      u1_ja[i][a]=Average(Econt_T[a], P_old_ia[i], Naa);
  }
  for(i=0; i<L; i++){
    float *p=P_new_ia[i]; double norm=0;
    for(a=0; a<Naa; a++){
      double h=0;
      if(SEC_STR && i_sec[i]>=0 && i_sec[i]<N_SS){
	h=E_loc_over_T[i_sec[i]][a];
      }
      for(j=0; j<L; j++){
	if(C_nat[i][j])h+=C_nat[i][j]*u1_ja[j][a];
      }
      p[a]=P_mut_a[a]*exp(-Lambda*h);
      norm+=p[a]; E_nat+=h;
    }
    for(a=0; a<Naa; a++)p[a]/=norm;
  }
  E->E_nat=E_nat/2;
}

void Update_P(float **P_new_ia, float **P_old_ia, float *P_mut_a,
	      float **u1_ja, float **u2_ja, float **u3_ja,
	      float **c1U1_ja, float **c1U2_ja,
	      float **Uc1U1_ja, float **Uc1U2_ja,
	      float **U2c1U1_ja, float **u1U_ja,
	      struct REM *E, 
	      float **C_nat, int *i_sec, char *c_sec,
	      float *Cont_L, float *T_ratio, float *factor,
	      float Lambda, int L, int Naa)
{
  float h_a[Naa];
  int i, j, a, error=0;
  //float x=exp(DG), Lambda=(x/(1+x))*Npop;
  //printf("lambda= %.2g\n", Lambda);
  //float T2=T_ratio*T_ratio;
 
  Set_to_zero(E);
  Average_U(u1_ja, u2_ja, u3_ja, P_old_ia, L, Naa);
  //Sum_ene(c1U1_ja, c1U2_ja, Uc1U1_ja, Uc1U2_ja, U2c1U1_ja, u1U_ja,
  //	  u1_ja, u2_ja, P_MF_ia, Cont_L, L, Naa);

  //float *c1U2=NULL;
  double E2cont2=0;
  if(E->REM){
    for(i=0; i<L; i++){
      float E2cont2_i[Naa], *c1U1=c1U1_ja[i]; 
      for(a=0; a<Naa; a++){E2cont2_i[a]=0; c1U1[a]=0;}
      for(j=0; j<L; j++){
	int l=abs(i-j); if(l<IJ_MIN)continue;
	float c=Cont_L[l];
	for(a=0; a<Naa; a++){
	  float cu1=c*u1_ja[j][a];
	  c1U1[a] += cu1;
	  if(j<i)E2cont2_i[a]+=C2_L[l]*cu1;
	}
      }
      E2cont2+=Average(E2cont2_i, P_old_ia[i], Naa);
      E->c1U1[i]=Average(c1U1_ja[i], P_old_ia[i], Naa);
      E->E1+=E->c1U1[i];
    }
  }
  E->E1/=2;

  float hnat_a[Naa];
  float h2cont2_a[Naa],h2cont1_a[Naa],h2site1_a[Naa],h2site2_a[Naa];
  float h3cont1_a[Naa],h3cont2_a[Naa],h3cont3_a[Naa],h3site1_a[Naa],
    h3site2_a[Naa];

  double h2cont2=0, h2cont1=0, h2site1=0, h2site2=0;
  double h3cont1=0, h3cont2=0, h3cont3=0, h3site1=0, h3site2=0;
  int di=0;
  double hnat, h1, h2;
  for(i=0; i<L; i++){
    double h_min=1000;
    if(E->REM>=2){
      di=L-1-i; if(di>i)di=i;
      if(di>DTAIL_MAX)di=DTAIL_MAX;
    }
    for(a=0; a<Naa; a++){
      hnat=0; h1=c1U1_ja[i][a]; h2=0;
      if(SEC_STR && i_sec[i]>=0 && i_sec[i]<N_SS){
	if(isnan(hnat)){error=11; goto Print_error;}
	hnat=E_loc_over_T[i_sec[i]][a];
      }
      if(E->REM>1){
	h2site1=nc_nc_L[di][0]*c1U1_ja[i][a]*c1U1_ja[i][a];
	h2cont1=0; h2cont2=0; h2site2=0;
	if(E->REM>2){
	  h3site1=nc_nc_Nc_L[di][0]*c1U1_ja[i][a]*c1U1_ja[i][a];
	  h3cont1=0; h3cont2=0; h3cont3=0; h3site2=0;
	}
      } 
      for(j=0; j<L; j++){
	if(C_nat[i][j]||C_nat[j][i])
	  hnat+=C_nat[i][j]*u1_ja[j][a];
	if(isnan(hnat)){
	  printf("j= %d C= %.2g u= %.3g\n", j, C_nat[i][j], u1_ja[j][a]);
	  error=12; goto Print_error;
	}
	if(E->REM<2)continue;
	int d, dj=L-1-j; if(dj>j)dj=j;
	if(di<=dj){d=di;}else{d=dj;}
	if(d>DTAIL_MAX)d=DTAIL_MAX;
	int l=abs(i-j);
	h2site2+=nc_nc_L[d][l]*c1U1_ja[j][a];
	if(E->REM==3)h3site2+=nc_nc_Nc_L[di][l]*c1U1_ja[j][a];
	if(l<IJ_MIN)continue;
	float c=Cont_L[l], cU=c*u1_ja[j][a];
	//h1+=cU;
	if(E->REM > 1){
	  if(isnan(Cont2_L[l])==0){
	    h2cont1+=Cont2_L[l]*u2_ja[j][a];
	  }
	  h2cont2+=C2_L[l]*cU;
	  if(E->REM==3){
	    h3cont1+=c*Cont_3L[l]*u3_ja[j][a];
	    h3cont2+=c*Cont_2L[l]*u2_ja[j][a];
	    h3cont3+=C3_L[l]*cU;
	  }
	}
      }

      if(isnan(hnat)){error=0; goto Print_error;}
      hnat_a[a]=hnat;
      h_a[a]=hnat;
      if(E->REM){
	h2cont1_a[a]=h2cont1; h2cont2_a[a]=h2cont2;
	h2site1_a[a]=h2site1; h2site2_a[a]=h2site2;
	h3cont1_a[a]=h3cont1; h3cont2_a[a]=h3cont2;
	h3cont3_a[a]=h3cont3; h3site1_a[a]=h3site1; h3site2_a[a]=h3site2;

	if(isnan(h1)){error=1; goto Print_error;}
	h_a[a]-=h1;
	if(E->REM > 1){
	  h2=(h2cont1+h2cont2*E->E1+E2cont2*h1+2*h2site2*h1+h2site1);
	  if(isnan(h2)){error=2; goto Print_error;}
	  h_a[a] += h2/2;
	  //hnat += (*factor)*h_misf; // Not good idea
	}
	if((a==0)||(h_a[a]< h_min))h_min=h_a[a];
	if(isnan(h_a[a])){error=3; goto Print_error;}
      }
    }
    double norm=0;
    float *p=P_new_ia[i];
    for(a=0; a<Naa; a++){
      double h=h_a[a]-h_min;
      p[a]=P_mut_a[a]*exp(-Lambda*h);
      norm+=p[a];
    }
    if(norm <= 0){
      printf("ERROR, P_MF vanishes!\nLambda= %.2g norm= %.2g f=%.2g i=%d\n",
	     Lambda, norm, *factor, i);
      printf("exponent of P= ");
      for(a=0; a<Naa; a++){printf("%.2g ", h_a[a]-h_min);}; printf("\n");
      printf("P_mut= ");
      for(a=0; a<Naa; a++){printf("%.2g ", P_mut_a[a]);}; printf("\n");
      exit(8);
    }
    for(a=0; a<Naa; a++)p[a]/=norm;
    float *P_old=P_old_ia[i];
    E->E_nat += Average(hnat_a, P_old, Naa);
    if(E->REM>=2){
      E->E2cont1  += Average(h2cont1_a, P_old, Naa);
      E->E2cont2  += Average(h2cont2_a, P_old, Naa);
      E->E2site1  += Average(h2site1_a, P_old, Naa);
      E->E2site2  += Average(h2site2_a, P_old, Naa);
      if(E->REM==3){
	E->E3cont1  += Average(h3cont1_a, P_old, Naa);
	E->E3cont2  += Average(h3cont2_a, P_old, Naa);
	E->E3cont3  += Average(h3cont3_a, P_old, Naa);
	E->E3site1  += Average(h3site1_a, P_old, Naa);
	E->E3site2  += Average(h3site2_a, P_old, Naa);
      }
    }
  }
  // Compute DG
  E->G_misf=G_misfold(E);
  E->DeltaG=DeltaG((float)E->E_nat, E->G_misf, E->S_U);
  if(isnan(E->G_misf)){
    printf("WARNING, G_misf= %.2g Tf= %.2g E1= %.2g E2= %.2g ",
	   E->G_misf, E->Tf, E->E1, E->E2); 
  }
  if(E->Tf>1){*T_ratio=E->Tf;}else{*T_ratio=1;}
  double DG=E->G_misf+E->S_U;
  if(DG < -30){*factor=1;}
  else if(DG > 30){*factor=0;}
  else{*factor=1./(1+exp(DG));}
  return;

 Print_error:
  printf("ERROR (type %d), P_MF is nan!\n"
	 "SITE %d out of %d amm %d h=%.2g hnat=%.2g h1=%.2g h2=%.2g"
	 " f=%.2g\n", error, i,L,a,h_a[a],hnat,h1,h2,*factor);
  printf("h2cont1=%.2g h2cont2=%.2g E1=%.2g E2cont2=%.2g "
	 "h2site2=%.2g h2site1=%.2g\n",
	 h2cont1, h2cont2, E->E1, E2cont2, h2site2, h2site1);
  for(int l=IJ_MIN; l<L; l++){
    if(isnan(Cont2_L[l]))
      printf("l=%d Cont2=%.2g Cont=%.2g\n",l,Cont2_L[l],Cont_L[l]);
  }
  exit(8);

}


float Test_convergence(float **P1_ia, float **P2_ia, int L, int Naa)
{
  double delta=0; int i, a;
  for(i=0; i<L; i++){
    float *p1=P1_ia[i], *p2=P2_ia[i], d;
    for(a=0; a<Naa; a++){
      d=p1[a]-p2[a]; delta+=d*d;
    }
  }
  return(delta/(Naa*L));
}

void Copy_P(float **P1_ia, float **P2_ia, int L, int Naa)
{
  int i, a;
  for(i=0; i<L; i++){
    float *p1=P1_ia[i], *p2=P2_ia[i];
    for(a=0; a<Naa; a++)p1[a]=p2[a];
  }
}
 
void Divide_by_Pmut(float **P_MF_ia, float *P_mut_a, int L, int Naa){
  int i, a;
  for(i=0; i<L; i++){
    for(a=0; a<Naa; a++)P_MF_ia[i][a]/=P_mut_a[a];
  }
}

void Initialize_P(float **P_MF_ia, float **C_nat, int *i_sec, char *c_sec,
		  float *Cont_L, float *hydro, float Lambda, int L, int Naa)
{
  printf("Initializing amino acid distribution as exp(L*(c_i-cave)/s*h[a])\n");
  // Compute average native contact matrix
  float csum[L];
  double ave=0; int i, j;
  for(i=0; i<L; i++){
    double c=0;
    for(j=0; j<L; j++){
      int l=abs(i-j); if(l<IJ_MIN)continue;
      if(C_nat[i][j])c+=C_nat[i][j];
      c-=Cont_L[l];
    }
    csum[i]=c; ave+=c;
  }
  ave/=L;

  // Compute BETA
  float equal=1./Naa;
  for(i=0; i<L; i++){
    float beta=Lambda*(csum[i]-ave);
    Boltzmann_distr(P_MF_ia[i], beta, hydro, Naa, equal);
    if(SEC_STR && i_sec[i]>=0 && i_sec[i]<N_SS){
      float *El=E_loc_over_T[i_sec[i]], *p=P_MF_ia[i];
      double sum=0; int a;
      for(a=0; a<Naa; a++){
	p[a]*=exp(-Lambda*El[a]); sum+=p[a];
      }
      if(sum){for(a=0; a<Naa; a++)p[a]/=sum;}
      else{for(a=0; a<Naa; a++)p[a]=equal;}
    }
  }

}

void Boltzmann_distr(float *P, float beta, float *hydro, int N, float equal)
{
  double sum=0; int a;
  for(a=0; a<N; a++){
    P[a]=exp(beta*hydro[a]);
    sum+=P[a];
  }
  if(sum){for(a=0; a<N; a++)P[a]/=sum;}
  else{for(a=0; a<N; a++)P[a]=equal;}
}

void Test_distr(//float *DG, float *lik, float *Tf,
		struct MF_results *MF_res, float **P_MF_ia,
		float **f_reg_ia, float *wi,
		float **C_nat, int *i_sec, char *c_sec,
		struct REM E_wt, int L, int Naa)
{
  Allocate_mean_field(&P_new_ia, &P_opt_ia,
		      &u1_ja, &u2_ja, &u3_ja, &c1U1_ja, &c1U2_ja,
		      &Uc1U1_ja, &Uc1U2_ja, &U2c1U1_ja, &u1U_ja,
		      E_wt.REM, L, Naa);
  struct REM E; E.c1U1=NULL;
  Initialize_E_REM(&E, E_wt.L, E_wt.REM, E_wt.T, E_wt.S_C, E_wt.S_U,
		   FILE_STR);

  float Lambda=1, P_mut_a[Naa], p=1./(float)Naa;
  for(int a=0; a<Naa; a++)P_mut_a[a]=p;
  float T_eff=1, factor=1;
  printf("Assessing site-specific distribution\n");
  Update_P(P_new_ia, P_MF_ia, P_mut_a, u1_ja, u2_ja, u3_ja, c1U1_ja, c1U2_ja,
	   Uc1U1_ja, Uc1U2_ja, U2c1U1_ja, u1U_ja,
	   &E, C_nat, i_sec, c_sec, Cont_L, &T_eff, &factor, Lambda,
	   L, Naa);
  /*MF_res->DG=Compute_DG_overT(&E, &T_ratio, &factor, u1_ja, u2_ja, u3_ja,
			     c1U1_ja, c1U2_ja, P_MF_ia, L, C_nat, i_sec, c_sec,
			     Cont_L);*/
  MF_res->DG=E.DeltaG;
  MF_res->Tf=E.Tf;
  MF_res->h=Hydro_ave(P_MF_ia, hydro, wi, L);
  //MF_res->entropy=Average_entropy(P_MF_ia, L);
  //MF_res->lik=Compute_lik(P_MF_ia, f_reg_ia, L, Naa);
  Compute_score(MF_res,P_MF_ia,MF_res->Lambda,f_reg_ia,wi,L,Naa,1);

  Empty_mean_field(&P_new_ia,&P_opt_ia,&u1_ja,&u2_ja,&u3_ja,&c1U1_ja,&c1U2_ja,
		   &Uc1U1_ja, &Uc1U2_ja, &U2c1U1_ja, &u1U_ja, L);
  Empty_E_REM(&E);
}


void Allocate_mean_field(float ***P_new_ia, float ***P_opt_ia,
			 float ***u1_ja, float ***u2_ja, float ***u3_ja,
			 float ***c1U1_ja, float ***c1U2_ja,
			 float ***Uc1U1_ja, float ***Uc1U2_ja,
			 float ***U2c1U1_ja, float ***u1U_ja,
			 int REM, int L, int Naa)
{
  if(Cont2_L)free(Cont2_L);
  if(*P_new_ia)Empty_matrix_f(*P_new_ia, L_old);
  if(*P_opt_ia)Empty_matrix_f(*P_opt_ia, L_old);
  if(*u1_ja)Empty_matrix_f(*u1_ja, L_old);
  if(*u2_ja)Empty_matrix_f(*u2_ja, L_old);
  if(*u3_ja)Empty_matrix_f(*u3_ja, L_old);
  if(*c1U1_ja)Empty_matrix_f(*c1U1_ja, L_old);
  if(*c1U2_ja)Empty_matrix_f(*c1U2_ja, L_old);

  Cont2_L=malloc(L*sizeof(float));
  for(int i=0; i<L; i++)Cont2_L[i]=Cont_L[i]*(1-Cont_L[i]);
  *P_new_ia=Allocate_mat2_f(L, Naa);
  *P_opt_ia=Allocate_mat2_f(L, Naa);
  *u1_ja=Allocate_mat2_f(L, Naa);
  *c1U1_ja=Allocate_mat2_f(L, Naa);
  if(REM>=2){
    *u2_ja=Allocate_mat2_f(L, Naa);
    *c1U2_ja=Allocate_mat2_f(L, Naa);
    //*u1U_ja=Allocate_mat2_f(L, Naa);
    //*Uc1U1_ja=Allocate_mat2_f(L, Naa);
  }
  if(REM==3){
    *u3_ja=Allocate_mat2_f(L, Naa);
    //*Uc1U2_ja=Allocate_mat2_f(L, Naa);
    //*U2c1U1_ja=Allocate_mat2_f(L, Naa);
  }
  L_old=L;
}

void Empty_mean_field(float ***P_new_ia, float ***P_opt_ia,
		      float ***u1_ja, float ***u2_ja, float ***u3_ja,
		      float ***c1U1_ja, float ***c1U2_ja,
		      float ***Uc1U1_ja, float ***Uc1U2_ja,
		      float ***U2c1U1_ja, float ***u1U_ja,
		      int L)

{
  if(*P_opt_ia){Empty_matrix_f(*P_opt_ia, L);} *P_opt_ia=NULL;
  if(*P_new_ia){Empty_matrix_f(*P_new_ia, L);} *P_new_ia=NULL;
  if(*u1_ja){Empty_matrix_f(*u1_ja, L);} *u1_ja=NULL;
  if(*u2_ja){Empty_matrix_f(*u2_ja, L);} *u2_ja=NULL;
  if(*u3_ja){Empty_matrix_f(*u3_ja, L);} *u3_ja=NULL;
  if(*u1U_ja){Empty_matrix_f(*u1U_ja, L);} *u1U_ja=NULL;
  if(*c1U1_ja){Empty_matrix_f(*c1U1_ja, L);} *c1U1_ja=NULL;
  if(*c1U2_ja){Empty_matrix_f(*c1U2_ja, L);} *c1U2_ja=NULL;
  if(*Uc1U1_ja){Empty_matrix_f(*Uc1U1_ja, L);} *Uc1U1_ja=NULL;
  if(*Uc1U2_ja){Empty_matrix_f(*Uc1U2_ja, L);} *Uc1U2_ja=NULL;
  if(U2c1U1_ja){Empty_matrix_f(*U2c1U1_ja, L);} *U2c1U1_ja=NULL;
}

float KL(float *P, float *Q, int n){ // Kullback-Leibler divergence
  int i; double KL=0;
  for(i=0; i<n; i++)if(P[i])KL+=P[i]*log(P[i]/Q[i]);
  return(KL);
}

float Average(float *u, float *p, int n){
  double u_ave=0;
  for(int a=0; a<n; a++)u_ave+=p[a]*u[a];
  return(u_ave);
}

void Average_P(float **P1, float **P2, int L, int Naa)
{
  for(int i=0; i<L; i++){
    float *P1i=P1[i], *P2i=P2[i];
    for(int a=0; a<Naa; a++){
      *P1i=0.5*(*P1i+*P2i); P1i++; P2i++;
    }
  }
}
