/* 
 * GRBEEFCONTAM - ground beef batch contamination probability
 *
 * Copyright 2015 Allan R. Willms
 *
 * -----------
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * -----------
 *
 * GRBEEFCONTAM computes the probability that each batch of ground beef is 
 * contaminated given that one batch is contaminated.
 *
 * INPUT is the name of a file with the follow information in it:
 *   The number of raw sources, S, (sources are numbered 0 to S-1), 
 *   the fat percentage, f, the susceptibility factor, g, the number of pieces per 
 *   carcass, p, and the average piece size, a, for each raw source.
 *   The carcass spread parameters for each source: K, N_1,...,N_K, Hm_1,...,Hm_K,
 *   and Hp_0,...,Hp_{K-1}
 *   The overlap probabilities V_{s1,s2} between sources s1 and s2, this matrix must 
 *   be symmetric.
 *   The total mass in each source, M, and the mass used from each source, m0, in
 *   batches prior to those being considered here.
 *   The number of batches, B, being considered (batches are numbered 0 to B-1), and 
 *   the mass, m_{b,s}, for each batch b input from source s.
 *   The number of hot batches to do the calculation for, and the list of this
 *   many batch numbers identifying the hot batches in turn.
 *   The layout of the file is:
 *  S
 *  f_1  f_2  ...  f_S
 *  g_1  g_2  ...  g_S
 *  p_1  p_2  ...  p_S
 *  a_1  a_2  ...  a_S
 *  K_1
 *  N_{1,1}  N_{1,2}  ...  N_{1,K_1}
 *  Hm_{1,1}  Hm_{1,2}  ...  Hm_{1,K_1}
 *  Hp_{1,0}  Hp_{1,1}  ...  Hp_{1,K_1-1}
 *  K_2
 *  N_{2,1}  N_{2,2}  ...  N_{2,K}
 *  ...
 *  K_S
 *  N_{S,1}  N_{S,2}  ...  N_{S,K_S}
 *  Hm_{S,1}  Hm_{S,2} ...  Hm_{S,K_S}
 *  Hp_{S,0}  Hp_{S,1} ...  Hp_{S,K_S-1}
 *  V_11
 *  V_21  V_22
 *  V_31  V_32  V_33
 *  ...
 *  V_S1  V_S2  V_S3  ...  V_SS
 *  M_1  M_2  ...  M_S
 *  m0_1  m0_2  ...  m0_S
 *  B
 *  m_11  m_12  ...  m_1S
 *  m_21  m_22  ...  m_2S
 *  ...
 *  m_B1  m_B2  ...  m_BS
 *  nhot
 *  h_1  ...  h_nhot
 *
 *  S : number of raw sources
 *  f : fraction of fat in each raw source
 *  g : contamination susceptibility factor for each raw source
 *  p : number of pieces per carcass in each raw source
 *  a : average piece size in each raw source
 *  K_s : number of sections in carcass distribution function for source s
 *  N_{s,k} : right boundary of section k for source s
 *  Hm_{s,k} : height of function in limit from left at N_{s,k}*p*a
 *  Hp_{s,k} : height of function in limit from right at N_{s,k}*p*a
 *  V_{s1,s2} : source carcass overlap fraction for sources s1 and s2
 *  M : Mass of each raw source
 *  m0 : mass from each raw source input to batches prior to those considered here
 *  B : number of batches being considered
 *  m_{b,s} : mass input to batch b from source s
 *  nhot : number of hot batches to consider
 *  h : hot batch numbers (nhot of them). P will be computed for each in turn.
 *
 *
 * OUTPUT:
 *  P : (B x nhot) matrix, giving the probability (as a percentage) that batch (row) 
 *      is contaminted given that batch (column) is contaminated.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAX(a,b) ((a) >= (b) ? (a) : (b))
#define MIN(a,b) ((a) <= (b) ? (a) : (b))
#define ACALC(a,b) (pow(1.0 - (a), (double) b))

float Qcalc(int n, float *L, float *Hp, float *Hm, float mu, float M, float *set);

int main(int argc, char *argv[]) {
	char *USEAGE_TEXT[] = {
	"",
    "GRBEEFCONTAM",
    "   Usage:  GRBEEFCONTAM input_filename [output_filename]",
    "      input_filename: name of input file specifying parameters for the problem.",
    "      output_filename: name of optional output file for the contamination probabilities.",
	"                       If not given, output is written to 'input_filename_out'",
    " ",
    " GRBEEFCONTAM computes the probability that various batches in ground beef",
	" production are contaminated.",
    " ",
    " Copyright 2015, Allan Willms,",
    " Dept. of Mathematics and Statistics,",
    " University of Guelph, Guelph, ON N1G 2W1, Canada"
    };
    #define NUM_USEAGE_TEXT_LINES (sizeof USEAGE_TEXT / sizeof(char *))

	FILE *file;
    char infname[FILENAME_MAX],outfname[FILENAME_MAX];
	float ***Q,***A,***Pcfromsinb,***Pcfromsinhandb;
	float **mu,**Pc1froms1equivc2froms2,**doubleproduct,**P,**Pcishot,**m,**cumMsb;
	float **L,**Hp,**Hm,**V;
	float *a,*pa,*Psishot,*f,*g,*M,*m0;
	int *p,*K,*N,*C,*h;
    int s,b,c,i;
	int S,B,nhot;
	float T0;

    if (argc < 2 || strncmp(argv[1],"-h",2) == 0) {
		for (i=0; i<NUM_USEAGE_TEXT_LINES; i++)
			fprintf(stdout,"%s\n",USEAGE_TEXT[i]);
		exit(0);
    }
	strncpy(infname,argv[1],FILENAME_MAX);
	if (argc < 3) 
		sprintf(outfname,"%s_out",argv[1]);
	else
		strncpy(outfname,argv[2],FILENAME_MAX);

	/* read in input file */
	if ((file = fopen(infname,"r")) == NULL) {
		fprintf(stderr,"Error opening input file %s.\n",infname);
		exit(1);
	}
	i = 1;
	if (fscanf(file,"%d",&S) != 1) {
		fprintf(stderr,"Error reading S on line %d of %s.\n",i,infname);
		exit(1);
	}
	if (S <= 0) {
		fprintf(stderr,"Error, S must be positive, S value read was %d.\n",S);
		exit(1);
	}
	/* allocate space for input parameters */
	f = (float *) malloc(S*sizeof(float));
	g = (float *) malloc(S*sizeof(float));
	p = (int *) malloc(S*sizeof(int));
	a = (float *) malloc(S*sizeof(float));
	pa = (float *) malloc(S*sizeof(float));
	K = (int *) malloc(S*sizeof(int));
	L = (float **) malloc(S*sizeof(float *));
	Hm = (float **) malloc(S*sizeof(float *));
	Hp = (float **) malloc(S*sizeof(float *));
	M = (float *) malloc(S*sizeof(float));
	m0 = (float *) malloc(S*sizeof(float));
	V = (float **) malloc(S*sizeof(float *));
	i++;
	for (s=0; s<S; s++) {
		if (fscanf(file,"%g",f+s) != 1) {
			fprintf(stderr,"Error reading f[%d] on line %d of %s.\n",s,i,infname);
			exit(1);
		}
		if (f[s] < 0.0 || f[s] > 1.0) {
			fprintf(stderr,"Error, f[%d] must be between 0.0 and 1.0, value read was %f.\n",s,f[s]);
			exit(1);
		}
	}
	i++;
	for (s=0; s<S; s++) {
		if (fscanf(file,"%g",g+s) != 1) {
			fprintf(stderr,"Error reading g[%d] on line %d of %s.\n",s,i,infname);
			exit(1);
		}
		if (g[s] < 0.0) {
			fprintf(stderr,"Error, g[%d] must be >= 0.0, value read was %f.\n",s,g[s]);
			exit(1);
		}
	}
	i++;
	for (s=0; s<S; s++) {
		if (fscanf(file,"%d",p+s) != 1) {
			fprintf(stderr,"Error reading p[%d] on line %d of %s.\n",s,i,infname);
			exit(1);
		}
		if (p[s] <= 0) {
			fprintf(stderr,"Error, p[%d] must be a positive integer, value read was %d.\n",s,p[s]);
			exit(1);
		}
	}
	i++;
	for (s=0; s<S; s++) {
		if (fscanf(file,"%g",a+s) != 1) {
			fprintf(stderr,"Error reading a[%d] on line %d of %s.\n",s,i,infname);
			exit(1);
		}
		if (a[s] <= 0.0) {
			fprintf(stderr,"Error, a[%d] must be positive, value read was %f.\n",s,a[s]);
			exit(1);
		}
	}
	for (s=0; s<S; s++) {
		pa[s] = p[s]*a[s];
		i++;
		if (fscanf(file,"%d",K+s) != 1) {
			fprintf(stderr,"Error reading K[%d] on line %d of %s.\n",s,i,infname);
			exit(1);
		}
		if (K[s] <= 0) {
			fprintf(stderr,"Error, K[%d] must be a positive integer, value read was %d.\n",s,K[s]);
			exit(1);
		}
		N = (int *) malloc(K[s]*sizeof(int));
		L[s] = (float *) malloc((2*K[s]+1)*sizeof(float));
		Hm[s] = (float *) malloc((2*K[s]+1)*sizeof(float));
		Hp[s] = (float *) malloc((2*K[s]+1)*sizeof(float));
		i++;
		for (c=0; c<K[s]; c++) {
			if (fscanf(file,"%d",N+c) != 1) {
				fprintf(stderr,"Error reading N_{%d,%d} on line %d of %s.\n",s,c,i,infname);
				exit(1);
			}
			if (N[c] <= 0) {
				fprintf(stderr,"Error, N[%d][%d] must be a positive integer, value read was %d.\n",s,c,N[c]);
				exit(1);
			}
		}
		L[s][K[s]] = 0.0;
		for (b=0; b<K[s]; b++) {
			L[s][K[s]+b+1] = pa[s]*N[b];
			L[s][K[s]-b-1] = -L[s][K[s]+b+1];
		}
		free(N);
		i++;
		for (c=0; c<K[s]; c++) {
			if (fscanf(file,"%g",Hm[s]+K[s]+1+c) != 1) {
				fprintf(stderr,"Error reading Hm_{%d,%d} on line %d of %s.\n",s,c,i,infname);
				exit(1);
			}
			if (Hm[s][K[s]+c+1] < 0.0) {
				fprintf(stderr,"Error, Hm[%d][%d] must be nonnegative, value read was %f.\n",
						s,c+1,Hm[s][K[s]+c+1]);
				exit(1);
			}
		}
		i++;
		for (c=0; c<K[s]; c++) {
			if (fscanf(file,"%g",Hp[s]+K[s]+c) != 1) {
				fprintf(stderr,"Error reading Hp_{%d,%d} on line %d of %s.\n",s,c,i,infname);
				exit(1);
			}
			if (Hp[s][K[s]+c] < 0.0) {
				fprintf(stderr,"Error, Hp[%d][%d] must be nonnegative, value read was %f.\n",
						s,c,Hp[s][K[s]+c]);
				exit(1);
			}
		}
		Hm[s][K[s]] = Hp[s][K[s]];
		Hp[s][2*K[s]] = 0.0;
		for (b=0; b<K[s]; b++) {
			Hm[s][K[s]-b-1] = Hp[s][K[s]+b+1];
			Hp[s][K[s]-b-1] = Hm[s][K[s]+b+1];
		}
		T0 = 0.0;
		for (b=K[s]+1; b<2*K[s]+1; b++)
			T0 += (L[s][b] - L[s][b-1])*(Hm[s][b] + Hp[s][b-1]);
		if (fabs(T0 - 1.0) > 1e-5) {
			fprintf(stderr,"Error, invalid carcass distribution function for source %d.\n%s\n",
					s,"Sum of (N_i-N_{i-1})*p*a*(Hm_i + Hp_{i-1}) must be 1.0.");
			exit(1);
		}
	}
	for (s=0; s<S; s++) {
		V[s] = (float *) malloc(S*sizeof(float));
		i++;
		for (c=0; c<S; c++) {
			if (fscanf(file,"%g",V[s]+c) != 1) {
				fprintf(stderr,"Error reading V_{%d,%d} on line %d of %s.\n",s,c,i,infname);
				exit(1);
			}
			if (V[s][c] < 0.0 || V[s][c] > 1.0) {
				fprintf(stderr,"Error, V[%d][%d] must be between 0.0 and 1.0, value read was %f.\n",
						s,c,V[s][c]);
				exit(1);
			}
		}
		for (b=0; b<s; b++) 
			if (V[b][s] != V[s][b])
				fprintf(stderr,"WARNING: V[%d][%d] not equal to V[%d][%d]; using first for both.\n",
						b,s,s,b);
	}
	i++;
	for (s=0; s<S; s++) {
		if (fscanf(file,"%g",M+s) != 1) {
			fprintf(stderr,"Error reading M[%d] on line %d of %s.\n",s,i,infname);
			exit(1);
		}
		if (M[s] <= 0.0) {
			fprintf(stderr,"Error, M[%d] must be positive, value read was %f.\n",s,M[s]);
			exit(1);
		}
	}
	i++;
	for (s=0; s<S; s++) {
		if (fscanf(file,"%g",m0+s) != 1) {
			fprintf(stderr,"Error reading m0[%d] on line %d of %s.\n",s,i,infname);
			exit(1);
		}
		if (m0[s] < 0.0) {
			fprintf(stderr,"Error, M[%d] must be nonnegative, value read was %f.\n",s,m0[s]);
			exit(1);
		}
	}
	i++;
	if (fscanf(file,"%d",&B) != 1) {
		fprintf(stderr,"Error reading B on line %d of %s.\n",i,infname);
		exit(1);
	}
	if (B <= 0) {
		fprintf(stderr,"Error, B must be a positive integer, value read was %d.\n",B);
		exit(1);
	}
	m = (float **) malloc(B*sizeof(float *));
	for (b=0; b<B; b++) {
		m[b] = (float *) malloc(S*sizeof(float));
		i++;
		for (s=0; s<S; s++) {
			if (fscanf(file,"%g",m[b]+s) != 1) {
				fprintf(stderr,"Error reading m_{%d,%d} on line %d of %s.\n",b,s,i,infname);
				exit(1);
			}
			if (m[b][s] < 0.0) {
				fprintf(stderr,"Error, m[%d][%d] must be nonnegative, value read was %f.\n",b,s,m[b][s]);
				exit(1);
			}
		}
	}
	for (s=0; s<S; s++) {
		T0 = 0.0;
		for (b=0; b<B; b++)
			T0 += m[b][s];
		if (T0 > M[s]) {
			fprintf(stderr,"WARNING, sum of mass input to all batches from source %d %s%d%s%f.\n%s%d%s%f.\n",
					s,"exceeds specified value of M[",s,"] = ",M[s],"M[",s,"] will be altered to ",T0);
			M[s] = T0;
		}
	}
	i++;
	if (fscanf(file,"%d",&nhot) != 1) {
		fprintf(stderr,"Error reading nhot on line %d of %s.\n",i,infname);
		exit(1);
	}
	if (nhot <= 0) {
		fprintf(stderr,"Error, nhot must be a positive integer, value read was %d.\n",nhot);
		exit(1);
	}
	h = (int *) malloc(nhot*sizeof(int));
	i++;
	for (c=0; c<nhot; c++) {
		if (fscanf(file,"%d",h+c) != 1) {
			fprintf(stderr,"Error reading h[%d] on line %d of %s.\n",c,i,infname);
			exit(1);
		}
		if (h[c] < 0 || h[c] >= B) {
			fprintf(stderr,"Error, h[%d] must be between 0 and %d, value read was %d.\n",c,B-1,h[c]);
			exit(1);
		}
	}
	fclose(file);

	/* allocate space */
	Q = (float ***) malloc(S*sizeof(float **));
	A = (float ***) malloc(S*sizeof(float **));
	Pcfromsinb = (float ***) malloc(S*sizeof(float **));
	mu = (float **) malloc(S*sizeof(float *));
	cumMsb = (float **) malloc(S*sizeof(float *));
	for (s=0; s<S; s++) {
		Q[s] = (float **) malloc(B*sizeof(float *));
		A[s] = (float **) malloc(B*sizeof(float *));
		Pcfromsinb[s] = (float **) malloc(B*sizeof(float *));
		cumMsb[s] = (float *) malloc((B+1)*sizeof(float));
	}
	C = (int *) malloc(S*sizeof(int));

	for (s=0; s<S; s++) {
		/* compute the cumulative mass in source s prior to batch b */
		cumMsb[s][0] = m0[s];
		for (b=0; b<B; b++) {
			cumMsb[s][b+1] = cumMsb[s][b] + m[b][s];
		}
		/* compute the number of carcasses in each source */
		pa[s] = p[s]*a[s];
		C[s] = (int) ceil(M[s]/pa[s]);
		/* Modify M so that there is an integer number of carcasses filling the source */
		M[s] = C[s]*pa[s];
		/* compute centre of carcass distributions */
		mu[s] = (float *) malloc(C[s]*sizeof(float));
		for (c=0; c<C[s]; c++)
			mu[s][c] = (c + 0.5)*pa[s];
		/* compute Q_sc for all possible mass ranges, compute P(c from s in b) */
		for (b=0; b<B; b++) {
			Q[s][b] = (float *) malloc(C[s]*sizeof(float));
			A[s][b] = (float *) malloc(C[s]*sizeof(float));
			Pcfromsinb[s][b] = (float *) malloc(C[s]*sizeof(float));
			for (c=0; c<C[s]; c++) {
				Q[s][b][c] = Qcalc(2*K[s]+1,L[s],Hp[s],Hm[s],mu[s][c],M[s],&(cumMsb[s][b]));
				A[s][b][c] = ACALC(Q[s][b][c], p[s]);
				Pcfromsinb[s][b][c] = 1.0 - A[s][b][c];
			}
		}
		free(mu[s]);
		free(cumMsb[s]);
		free(L[s]);
		free(Hp[s]);
		free(Hm[s]);
	}
	free(mu);
	free(cumMsb);
	free(L);
	free(Hp);
	free(Hm);
	free(K);
	free(pa);

	Pc1froms1equivc2froms2 = (float **) malloc(S*sizeof(float *));
	for (s=0; s<S; s++) {
		Pc1froms1equivc2froms2[s] = (float *) malloc(S*sizeof(float));
		/* compute P(c1 from s1 equiv c2 from s2) */
		for (i=0; i<S; i++)
			Pc1froms1equivc2froms2[s][i] = (V[s][i]*(C[s] + C[i]))/((1+V[s][i])*(C[s]*C[i]));
		free(V[s]);
	}
	free(V);
	/* the variable doubleproduct is the product on s2 and c2 of  (s2 != s1)
	 * (1-Prob(c1 from s1 equiv c2 from s2)*Prob(c2 from s2 in b)
	 * the columns are the different batches b, and the rows are for different s1
	 */
	doubleproduct = (float **) malloc(S*sizeof(float *));
	for (s=0; s<S; s++) {
		doubleproduct[s] = (float *) malloc(B*sizeof(float));
		for (b=0; b<B; b++) {
			doubleproduct[s][b] = 1.0;
			for (i=0; i<S; i++)
				if (i != s)
					for (c=0; c<C[i]; c++)
						doubleproduct[s][b] *= 1.0 - Pc1froms1equivc2froms2[s][i]*Pcfromsinb[i][b][c];
		}
		free(Pc1froms1equivc2froms2[s]);
	}
	free(Pc1froms1equivc2froms2);
	/* compute main probability for each value of h in the vector */
	/* allocate space */
	Pcfromsinhandb = (float ***) malloc(S*sizeof(float **));
	Pcishot = (float **) malloc(S*sizeof(float *));
	for (s=0; s<S; s++) {
		Pcfromsinhandb[s] = (float **) malloc(B*sizeof(float *));
		for (b=0; b<B; b++) {
			Pcfromsinhandb[s][b] = (float *) malloc(C[s]*sizeof(float));
		}
		Pcishot[s] = (float *) malloc(C[s]*sizeof(float));
	}
	P = (float **) malloc(B*sizeof(float *));
	for (b=0; b<B; b++) {
		P[b] = (float *) malloc(nhot*sizeof(float));
	}
	Psishot = (float *) malloc(S*sizeof(float));
	for (i=0; i<nhot; i++) {
		/* compute P(c is hot) */
		for (s=0; s<S; s++) {
			if (m[h[i]][s] == 0.0) {
				for (c=0; c<C[s]; c++) 
					Pcishot[s][c] = 0.0;
			}
			else {
				T0 = 0.0;
				for (c=0; c<C[s]; c++) {
					Pcishot[s][c] = Pcfromsinb[s][h[i]][c];
					T0 += Pcishot[s][c];
				}
				for (c=0; c<C[s]; c++)
					Pcishot[s][c] /= T0;
			}
		}
		/* compute probability that source s is hot: P(s is hot) */
		T0 = 0.0;
		for (s=0; s<S; s++) {
			Psishot[s] = g[s]*f[s]*m[h[i]][s];
			T0 += Psishot[s];
		}
		for (s=0; s<S; s++) 
			Psishot[s] /= T0;
		/* compute P(c from s in h and j) */
		for (s=0; s<S; s++)  {
			for (b=0; b<B; b++)  {
				if (b == h[i]) {
					for (c=0; c<C[s]; c++)
						Pcfromsinhandb[s][h[i]][c] = 1 - A[s][h[i]][c];
				}
				else {
					for (c=0; c<C[s]; c++)
						Pcfromsinhandb[s][b][c] = 1 - (A[s][h[i]][c] + A[s][b][c] -
							ACALC(Q[s][h[i]][c] + Q[s][b][c], p[s]));
				}
			}
		}
		/* compute contamination probability */
		for (b=0; b<B; b++)
			P[b][i] = 0.0;
		for (s=0; s<S; s++) {
			for (b=0; b<B; b++) {
				for (c=0; c<C[s]; c++) {
					if (Pcfromsinb[s][h[i]][c] != 0.0) { /* If 0 then Pcishot will be zero */
						T0 = 1.0 - Pcfromsinhandb[s][b][c]/Pcfromsinb[s][h[i]][c];
						P[b][i] += Psishot[s]*Pcishot[s][c]*(1.0 - doubleproduct[s][b]*T0);
					} 
				}
			}
		}
	}
	/* write output file */
	if ((file = fopen(outfname,"w")) == NULL) {
		fprintf(stderr,"Error opening output file %s.\n",outfname);
		exit(1);
	}
	for (b=0; b<B; b++) {
		for (i=0; i<nhot; i++)
			fprintf(file,"%5.1f ",100.0*P[b][i]);
		fprintf(file,"\n");
	}
	/* free space */
	free(f);
	free(g);
	free(Psishot);
	for (s=0; s<S; s++) {
		free(Pcishot[s]);
		free(doubleproduct[s]);
		for (b=0; b<B; b++) {
			free(A[s][b]);
			free(Q[s][b]);
			free(Pcfromsinb[s][b]);
			free(Pcfromsinhandb[s][b]);
		}
		free(A[s]);
		free(Q[s]);
		free(Pcfromsinb[s]);
		free(Pcfromsinhandb[s]);
	}
	free(Pcishot);
	free(doubleproduct);
	free(A);
	free(Q);
	free(Pcfromsinb);
	free(Pcfromsinhandb);
	for (b=0; b<B; b++) {
		free(m[b]);
		free(P[b]);
	}
	free(m);
	free(P);
		
	
}

float Qcalc(int n, float *L, float *Hp, float *Hm, float mu, float M, float *set) {
	/* The expected fraction of a carcass centred at mu present in the given set of positions
	 * The set should be a column vector with two entries defining the start and end of the 
	 * interval.  The density function is G(x) = F(x-mu) + F(-x-mu) + F(-x-mu+2M) where 
	 * F(x) = ((|x|/p*a - N_{i-1})*Hm_i + (N_i - |x|/pa)*Hp_{i-1}) / (N_i - N_{i-1}),  
	 *       N_{i-1} <= |x|/pa < N_i 
	 * and F is even around zero.  L is N*p*a;  L, Hp, and Hm are 2K+1 vectors giving the 
	 * values of N_{i}*p*a, H^{+}_{i-1}, and H^{-}_{i}.
	 */
	float Q;
	float limits[3][2],lim[2],bnd_dist,height[2];
	int i,startK,endK,k;

	/* Set up limits for the centred piece and the reflections at either end. */
	for (k=0; k<2; k++) {
		limits[0][k] = set[k] + (mu-2*M);
		limits[1][k] = set[k] - mu;
		limits[2][k] = set[k] + mu;
	}
	Q = 0.0;
	for (i=0; i<3; i++) {
		for (k=0; k<2; k++)
			lim[k] = limits[i][k];
		if (lim[0] >= L[n-1] || lim[1] <= L[0])
			continue;
		else {
			lim[0] = MAX(lim[0],L[0]);
			lim[1] = MIN(lim[1],L[n-1]);
		}
		startK = 0;
		while (startK < n-1 && lim[0] >= L[startK + 1])  startK++;
		endK = startK + 1;
		while (endK < n-1 && lim[1] > L[endK])  endK++;
		/* compute any middle full pieces */
		for (k=startK+1; k<endK-2; k++) 
			Q += 0.5*(Hp[k] + Hm[k+1])*(L[k+1] - L[k]);
		/* add in the initial piece */
		height[0] = ((lim[0] - L[startK])*Hm[startK+1] + (L[startK+1] - lim[0])*Hp[startK])/
			(L[startK+1] - L[startK]);
		if (L[startK + 1] < lim[1]) {
			height[1] = Hm[startK + 1];
			bnd_dist = L[startK + 1] - lim[0];
		}
		else {
			height[1] = ((lim[1] - L[startK])*Hm[startK+1] + (L[startK+1] - lim[1])*Hp[startK])/
				(L[startK+1] - L[startK]);
			bnd_dist = lim[1] - lim[0];
		}
		Q += 0.5*(height[0] + height[1])*bnd_dist;
		/* add in a final piece if necessary */
		if (endK > startK + 1) {
			height[0] = Hp[endK-1];
			height[1] = ((lim[1] - L[endK-1])*Hm[endK] + (L[endK] - lim[1])*Hp[endK-1])/
					(L[endK] - L[endK-1]);
			Q +=  0.5*(height[0] + height[1])*(lim[1] - L[endK-1]);
		}
	}
	return Q;
}
