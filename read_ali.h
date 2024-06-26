int Read_ali(char ***msa_seq, int *L_ali, char ***name_seq, char *file_ali,
	     char *pdb_buff);
int Align_ali_PDB(//short **ali_seq,
		  int *ali_PDB,
		  float *seq_id, float *Seq_id_ave, 
		  char **msa_seq, char **name_seq, int N_seq, int L_ali,
		  short *i_seq, int L_seq_PDB, int L_str_PDB
		  //char *pdbname, int *ali_PDB,struct protein *prot
		  );
