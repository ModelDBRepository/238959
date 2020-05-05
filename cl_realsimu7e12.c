// cl_realsimu7e12.c:
// Based on cl_realsimu7e10.c. Here I modify the presynaptic inhibition mechanism so that it affects the vesicle release
//
// Notes:
// [Oct 15, 14] ver 7e12: modifying the mechanism of presynaptic inhibition. Now the presynaptic inhibition also affects
//              short-term plasticity by reducing the amount of vesicle release.
// [Oct 12, 13] ver 7e11: adding new parameters "RefractT" for neuron. The parameter indicates the refractory period
//              for each neuron in the unit of time step.
// [Aug 26, 10] ver 7e10: [Ca2+]^n, n=3.5
// [Jun 22, 10] ver 7e9: add pre-synaptic inhibition
// [May 7, 10] ver 7e8: remove voltage dependence in NMDA synaptic efficacy. Change recpetor naming:
//              NMDA->GABAB, GABA->GABAA, AMPA->Ach
// [Mar 13, 10] ver 7e7: add FreqExtDecay as a changable event to allow decaying external frequency
// [Mar 11, 10] ver 7e6: increase maximum number of target population and decrease the
//                      axonal delay for spikes to 1ms
// [Mar 11, 10] ver 7e6: change the connectivity to one-to-one
// [Apr 15, 08] make output file names changable
// [Apr 5, 05]  change the algorithm for the noise generator.
// [Mar 3, 05]  Add FreqExtSD into the structure PopDescr to describe
//              the standard deviation of the time variation of FreqExt.
//              Note that due to the design of the program structure,
//              FreqExtSD is incompatible with ExtEffSD. If FreqExtSD!=0
//              ExtEffSD will not function. All external efficacies will
//              be set to MeanExtEff.
// [Mar 2, 05]  change the outward rectifying K channel to be time-
//              and voltage-dependent, just as described in Compte-A-2003a
// [Feb 28, 05] make the synaptic efficacy changeable throug .pro file:
//              modified: structure EventDescr, function ParseProtocol()
//              and function HandleEvent().
// [Feb 03, 05] add distribution width to the synaptic efficacy
// [Oct 07, 04] Change max event number from 100 to 200 (line 255)
// [Oct 06, 04] Change the short-term facilitation and short-term
//              depression from cell-wide properties to synapses properties.
// [Oct 06, 04] Change synaptic delay from 5ms to 2ms
// [Oct 05, 04] Add short-term facilitation
// [Sep 22, 04] change NumberofTraces from 2 to 5 (line 1197)
// [Sep 22, 04] Add potassium outward rectifying conductances.
// [Sep 21, 04] fix a bug in DescribeNetwork() relating to the initialization
//              of the population parameters. This bug makes the program initialize
//              the short-term-depression related parameters only for the first
//              population and neglect the rest populations.
// [Sep 21, 04] add potassium inward rectifying conductances (ADR).
//
// [Jun 08, 04] modyfying SaveSpikes(). Use the sliding window instead of the
//              exponential filter.
// [Jun 07, 04] adding code to initialize parameters for short-term depression
//              in case that they are not givien in the configuration file.
// [Jun 02, 04] adding short-term depression
// [Mar 04, 04] adding spike-frequency adaptation from the paper Wang-X-1999a
//
// by Chung-Chuan Lo
//
//###########################################
// Original notes from realsimu8.c
//
// realsimu8 - Ver 0.8
//
// 23-9-03 Number of synapses in GenerateNetwork: tested
// 28-9-03 Parser for the description introduced
// 29-9-03 Protocol parser added
// 30-9-03 Mg block added. Vaux added
// 1-10-03 NMDA saturation
// 2-10-03 bug (use of Tableofspikes) fixed. Sourceneuron introduced
// 8-10-03 bug fixed: two TimeLastemitted spikes are necessary, otherwise the decay for NMDA s variable is always the delay
// 13-10-03 bug fixed: DeltaLS should be \sum_k \delta_k (1-s_k)alpha w_k and not the sum of all s_k
// 15-10-03 successful test for XJ model
// 15-10-03 multitrial version (ver 0.8)

/* How to add a new variable: 1) add it to the proper structure and make a copy in the
description structure 2) DescribeNetwork should initialize the description structure
3) GenerateNetwork should initialize the variable for each neuron/synapse


CAUTION: no saturation on external NMDA!!!
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdarg.h>
#include "rando2.h"


/*
Some stuff from Numerical Recipes
*/

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)
//#define PI 3.141592654

long *idum;

void sran1(long *iseed)
{
  idum = iseed;
}

float ran1()
{
	int j;
	long k;
	static long iy=0;
	static long iv[NTAB];
	float temp;

	if (*idum <= 0 || !iy) {
		if (-(*idum) < 1) *idum=1;
		else *idum = -(*idum);
		for (j=NTAB+7;j>=0;j--) {
			k=(*idum)/IQ;
			*idum=IA*(*idum-k*IQ)-IR*k;
			if (*idum < 0) *idum += IM;
			if (j < NTAB) iv[j] = *idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ;
	*idum=IA*(*idum-k*IQ)-IR*k;
	if (*idum < 0) *idum += IM;
	j=iy/NDIV;
	iy=iv[j];
	iv[j] = *idum;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
}


float gasdev()
{
	float ran1();
	static int iset=0;
	static float gset;
	float fac,rsq,v1,v2;

	if  (iset == 0) {
		do {
			v1=2.0*ran1()-1.0;
			v2=2.0*ran1()-1.0;
			rsq=v1*v1+v2*v2;
		} while (rsq >= 1.0 || rsq == 0.0);
		fac=sqrt(-2.0*log(rsq)/rsq);
		gset=v1*fac;
		iset=1;
		return v2*fac;
	} else {
		iset=0;
		return gset;
	}
}


//Below is the real stuff


#define MAXTN 10000 // max number of target neurons

// types of receptors

#define MAXRECEPTORS 3

#define ACh 0
#define GABAA 1
#define GABAB 2

// =============================================================
// Structures for the simulation

// --------------------- SYNAPSE structure ---------------------

typedef struct {
  int TargetPop;
  int TargetNeuron;
  int TargetAxon;
  float Efficacy; // change in the conductance due to a single spike (nS) (hence it is always positive, IPSP come out from rev pots)
  float LastConductance; // synaptic conductance when the last spike was emitted (useful for NMDA saturation), -1 if not NMDA (and hence no saturation)
  int TargetReceptor; // 0=Ach, 1=GABAA, 2=GABAB
  //Parameters and variables for the short term depression
  float D;                 // The fraction of available vesicle.
  float pv;                // Each spike reduce D by the factor pv.
  float tauD;             // The time constant for the speed of D recovery.

  //Parameters and variables for the short term facilitation
  float F;                 // The fraction of available vesicle.
  float Fp;                // Each spike increase F by the factor Fp.
  float tauF;             // The time constant for the speed of F decrease.

  // V. 7e9 for pre-synaptic inhibition
  float Ca_pre;   //Relative Ca level (max=1) at the presynaptic terminal
  float last_pre_update[MAXRECEPTORS];  // Time of the last update of Ca_pre
  float St_pre[MAXRECEPTORS];  // Gating variable on the pre-synaptic terminal
  float Ca_pre_tau[MAXRECEPTORS];        // Recovery time constant of Ca_pre (copied from the presynaptic time constant)

} Synapse;

// --------------------- NEURON structure -----------------------


typedef struct {

  // outputs

  int Naxonals;    // axonal tree (total number of synapses on the axon)
  Synapse *Axonals;

  // inputs

  int Nreceptors;
  float Tau[MAXRECEPTORS]; // tau for each conductance
  float LS[MAXRECEPTORS]; // total conductance for internal inputs (spikes generated by the network)
  float RevPots[MAXRECEPTORS]; // reversal potentials
  int   MgFlag[MAXRECEPTORS]; // 1 is magnesium block is active, 0 otherwise

  // external gaussian inputs for each receptor (to be added to S at each time step)
  float ExtS[MAXRECEPTORS]; // nS

  float ExtMuS[MAXRECEPTORS]; // nS/s
  float ExtSigmaS[MAXRECEPTORS]; // nS/s^2

  // here I should consider also the external inputs!!!

  // dynamic variables

  float V;
  int RefrState; // Refractory state counter
  float C; // capacitance in nF
  float Taum; // membrane time constant
  float RestPot; // Resting potential
  float ResetPot; // reset potential
  float Threshold; // threhsold
  float RefractT;    // Refractory period (ms)


  // Spike-frequency adaptation.
  // The properties of calcium-activated potassium current I

  float Ca;  //Concentration of Ca ions
  float Vk;  // resting potential
  float alpha_ca; // Amount of increment of [Ca] with each spike discharge. (muM)
  float tau_ca; // time constant
  float g_ahp; // efficacy

  // the times of the last two spikes
  float TimeLastSpike; // time of the last emitted spike (for NMDA saturation)
  float PTimeLastSpike;



  //Anomalous delayed rectifier (ADR). Simulating the up and down state in striatal neurons.
  float g_adr_max;  //Maximun value of the g
  float Vadr_h; //Potential for g_adr=0.5g_adr_max
  float Vadr_s; //Slop of g_adr at Vadr_h, defining how sharp the shape of g_ard is.
  float ADRRevPot; //reverse potential for ADR.
  float g_k_max;  // Maximun outward rectifying current
  float Vk_h; //Potential for g_k=0.5g_k_max
  float Vk_s; //Slop of g_k at Vk_h, defining how sharp the shape of g_k is.
  float tau_k_max; //maximum value of tau_k
  float n_k;  // gating variable for outward rectifying K channel;


} Neuron;



// ------------------- Population structure ---------------------
#define MAXDELAYINDT 50
#define MAXSPIKESINDT 2000

typedef struct {
  char Label[100];

  int Ncells;
  Neuron *Cell;

  // Table of spikes emitted by the neurons of this population
  int CTableofSpikes; // pointer where we are (time t)
  int DTableofSpikes; // pointer to the t-delay
  int NTableofSpikes[MAXDELAYINDT];
  int TableofSpikes[MAXDELAYINDT][MAXSPIKESINDT];

} Population;

#define MAXP 300
#define MAXEXTMEM 20

int Npop;
Population Pop[MAXP];

// global variables for the simulation

float dt=0.1;      // in ms
float Time;    // absolute time from the beginning of the trial
int delayindt; // transmission delay in dt
int flagverbose=1; // flag to activate verbose messages
int FlagSaveAllSpikes=1;
float ext_freq; // external noise
// =================================================================
// Descriptors to generate the network structure

typedef struct {
  float Connectivity; // mean fraction of randomly connected target neurons
  float TargetReceptor; // 0=Ach, ...
  float MeanEfficacy; // mean efficacy (for initialization)
  float EfficacySD; // Standard deviation of the efficacy distribution.
  float pv;  // Each spike reduce the fraction D of available vesicle by the factor pv.
  float tauD; // The time constant for the speed of D recovery.
  float Fp;   // Each spike increase F by the factor Fp.
  float tauF; // The time constant for the speed of F decrease.
//  float Ca_pre_tau;  // Time time constant for the pre-synaptic inhibition
} SynPopDescr;

typedef struct {
  char Label[100];
  int Ncells;

  int NTargetPops;
  int TargetPops[MAXP];
  int TargetAxons[MAXP];

  SynPopDescr SynP[MAXP];

  int Nreceptors;
  char  ReceptorLabel[MAXRECEPTORS][100]; // label for the receptor (needed for the compiler)
  float Tau[MAXRECEPTORS]; // tau for each conductance
  float RevPots[MAXRECEPTORS]; // reversal potentials
  int   MgFlag[MAXRECEPTORS]; // magnesium block flag

  // external input (externally defined)
  float MeanExtCon[MAXRECEPTORS]; // mean total number of external connections
  float MeanExtEff[MAXRECEPTORS]; // external mean efficacy
  float ExtEffSD[MAXRECEPTORS]; // Deviation of distribution of external efficacy
  float FreqExt[MAXRECEPTORS]; // external frequency in Hz
 float FreqExtTau[MAXRECEPTORS]; // time constant of the decaly from old FreqExt to new FreqExt
  float FreqExtSD[MAXRECEPTORS]; // external frequency in Hz
  float FreqExtMem[MAXRECEPTORS]; // memory in the external noise
  float FreqExtDecay[MAXRECEPTORS]; //decay factor of the external noise
  float FreqExtNorm[MAXRECEPTORS]; //Normalization factor for the memory of the external noise

  // external input (statistics internally calculated)
  float MeanExtMuS[MAXRECEPTORS]; // statistics of the external input nS/s
  float MeanExtSigmaS[MAXRECEPTORS];
  // dynamic variables

  float C; // capacitance
  float Taum; // membrane time constant
  float RestPot; // Resting potential
  float ResetPot; // reset potential
  float Threshold; // threhsold
  float RefractT; // Refractory period (ms)

  float Vk;  // resting potential
  float alpha_ca; // Amount of increment of [Ca] with each spike discharge. (muM)
  float tau_ca; // time constant
  float g_ahp; // efficacy

  //Anomalous delayed rectifier (ADR)
  float g_adr_max;  //Maximun value of the g
  float Vadr_h; //Potential for g_adr=0.5g_adr_max
  float Vadr_s; //Slop of g_adr at Vadr_h, defining how sharp the shape of g_ard is.
  float ADRRevPot; //Reverse potential for ADR
  float g_k_max;  // Maximun outward rectifying current
  float Vk_h; //Potential for g_k=0.5g_k_max
  float Vk_s; //Slop of g_k at Vk_h, defining how sharp the shape of g_k is.
  float tau_k_max; //maximum tau for outward rectifying k current


} PopDescr;

PopDescr PopD[MAXP];

// ===============================================================
// Protocol descriptors

#define MAXEVENTS 200

#define ENDOFTRIAL 1
#define CHANGEEXTFREQ 2
#define CHANGEMEANEFF 3
#define CHANGEEXTFREQSD 4

typedef struct {
  int Type;
  float ETime;
  int PopNumber;
  int TargetPopNumber;
  int ReceptorNumber;
  float FreqExt;
  float FreqExtDecay;
  float FreqExtSD;
  float FreqExtsd;
  float MeanEff;
  char Label[100];
} EventDescr;

int NEvents=0;
int CEvent; // current event
float NextEventTime=0.;
float TrialDuration=1000.; // in ms
EventDescr Events[MAXEVENTS];

// Multitrial version:

int CurrentTrial;
int NumberofTrials;


// ===============================================================================
// AUXILIARY SUBROUTINES
// ===============================================================================

/* --------------------------------------------------------------------------------
   Report: prints on screen and log file (dev_oss filename device)
   -------------------------------------------------------------------------------- */

void report(char *fmt,...)
{
  static FILE *devlog;
  static int initflag=1;
  va_list ap;
  char *s,fmtemp[20];
  int ival,fmtemp_i;
  float fval;
  char *cval;

  if(initflag) {
    devlog=fopen("simu.log","w");
    if(devlog==NULL) return;
    initflag=0;
  }

  va_start(ap,fmt);
  for(s=fmt; *s; s++) {

    if(*s!='%') { if(flagverbose) printf("%c",*s); fprintf(devlog,"%c",*s); continue; }

    fmtemp_i=0; while (*s && !(*s=='s' || *s=='d' || *s=='f')) { fmtemp[fmtemp_i]=*s; s++; fmtemp_i++; }
    if(*s=='d') {  ival=va_arg(ap,int); fmtemp[fmtemp_i]='d'; fmtemp[fmtemp_i+1]=0;
                   if(flagverbose) printf(fmtemp,ival); fprintf(devlog,fmtemp,ival);  }
    if(*s=='f') {  fval=va_arg(ap,double); fmtemp[fmtemp_i]='f'; fmtemp[fmtemp_i+1]=0;
                   if(flagverbose) printf(fmtemp,fval); fprintf(devlog,fmtemp,fval);  }
    if(*s=='s') {  cval=va_arg(ap,char *); fmtemp[fmtemp_i]='s'; fmtemp[fmtemp_i+1]=0;
                   if(flagverbose) printf(fmtemp,cval); fprintf(devlog,fmtemp,cval);  }
  } /* end for */
  va_end(ap);
}



// ===============================================================
//
// INPUT
// ==============================================================
// auxiliary subroutines for parsing the descriptor file
// returns the code of the population with name s (-1 in case of error)

int PopulationCode(char *s)
{
  int p;


  for(p=0;p<Npop;p++)
    {

      if(strcmp(PopD[p].Label,s)==0) {
	return p;
      }
    }
  return -1;
}


// returns the code of the receptor with name s for pop p (-1 in case of error)

int ReceptorCode(char *s,int p)
{
  int r;

  for(r=0;r<PopD[p].Nreceptors;r++)
    {
      if(strcmp(PopD[p].ReceptorLabel[r],s)==0) {
	return r;
      }
    }
  return -1;
}

// =======================================================================
// DescribeNetwork()
// Initializes the descriptors of the network by parsing cl_network.conf
// =======================================================================

#define EXC 0
#define INH 1


int DescribeNetwork()
{

  FILE *devconf;
  char buf[1000],*s,*es;
  int currentpopflag=0;
  int currentreceptorflag=0;
  int line,auxi;
  int currentpop,currentreceptor,currenttarget,currenttargetaxon;
  float aux;
  float std_p, std_tau;
  int flag_d;

  // FIRST PASSAGE
  // -------------------------------------------------------------------
  // the parser has to go over the file twice: the first time reads all
  // the labels and initializes the number of populations and the number
  // of receptors
  report("Parsing network configuration... first passage\n");

  /*  strncpy=(network_conf,prefix,strlen(prefix));
  strcpy=(network_conf+strlen(prefix),"network.conf");
  devconf=fopen(network_conf,"r");*/
  devconf=fopen("network.conf","r");
  if(devconf==NULL) { printf("ERROR:  Unable to read configuration file\n"); return 0;}

  Npop=0; line=-1;

  while(fgets(buf,1000,devconf)!=NULL)
    {
      line++;
      s=buf;
      // trim \n at the end
      es=buf; while(*es && *es!='\n') es++; *es=0;

      while(*s==' ' || *s=='\t') s++; // skip blanks
      if(*s==0) continue; // empty line
      if(*s=='%') continue; // remark

      // commands for defining a new population
      if(strncmp(s,"NeuralPopulation:",17)==0)
	{
	  currentpopflag=1;
	  s+=17; while(*s==' ') s++;
	  strcpy(PopD[Npop].Label,s);
	  report("Population: %s\n",PopD[Npop].Label);
	  PopD[Npop].Nreceptors=0;
	  continue;
	}

      if(strncmp(s,"EndNeuralPopulation",19)==0)
	{
	  currentpopflag=0;
	  Npop++;
	  continue;
	}

      // command for defining a receptor
      if(strncmp(s,"Receptor:",9)==0) {
	currentreceptorflag=1;
	s+=9; while(*s==' ') s++;
	strcpy(PopD[Npop].ReceptorLabel[PopD[Npop].Nreceptors],s);
	report("Receptor %d: %s\n",PopD[Npop].Nreceptors,PopD[Npop].ReceptorLabel[PopD[Npop].Nreceptors]);
      }

      if(strncmp(s,"EndReceptor",11)==0) {
	currentreceptorflag=0;
	PopD[Npop].Nreceptors++;
      }
    }

  fclose(devconf);

  // Second passage: now all the parameters and the target populations are parsed
  // ----------------------------------------------------------------------------

  report("Parsing network configuration... second passage\n");

  //  devconf=fopen(network_conf,"r");
  devconf=fopen("network.conf","r");
  if(devconf==NULL) { printf("ERROR:  Unable to read configuration file\n"); return 0;}

  line=-1;

  while(fgets(buf,1000,devconf)!=NULL)
    {
      line++;
      s=buf;

      // trim \n at the end
      es=buf; while(*es && *es!='\n') es++; *es=0;

      while(*s==' ' || *s=='\t') s++; // skip blanks
      if(*s==0) continue; // empty line
      if(*s=='%') continue; // remark

      // population commands
      if(strncmp(s,"NeuralPopulation:",17)==0)
	{
	  currentpopflag=1;
	  s+=17; while(*s==' ') s++;
	  currentpop=PopulationCode(s);
	  if(currentpop==-1) { printf("Unknown population [%s]: line %d\n",s,line); return -1;}
	  PopD[currentpop].NTargetPops=0;
	  report("-----------------------------------\n    Population: %s\n-----------------------------------\n",PopD[currentpop].Label);

	  //Initialize some population parameters
         PopD[currentpop].RefractT=2.0;
          PopD[currentpop].g_adr_max=0;
          PopD[currentpop].Vadr_h=-100;
          PopD[currentpop].Vadr_s=10;
          PopD[currentpop].ADRRevPot=-90;
          PopD[currentpop].g_k_max=0;
          PopD[currentpop].Vk_h=-34;
	  PopD[currentpop].Vk_s=6.5;
          PopD[currentpop].tau_k_max=8;
	  continue;
	}

      if(strncmp(s,"EndNeuralPopulation",19)==0)
	{
	  report("EndPopulation\n");
	  currentpopflag=0;
	  continue;
	}

      // paramters for the population


      if(strncmp(s,"N=",2)==0 && currentpopflag) {
	PopD[currentpop].Ncells=atoi(s+2); report("  N=%d\n",PopD[currentpop].Ncells); continue;
      }
      if(strncmp(s,"C=",2)==0 && currentpopflag) {
	PopD[currentpop].C=atof(s+2);report("  C=%f nF\n",(double)PopD[currentpop].C); continue;
      }
      if(strncmp(s,"Taum=",5)==0 && currentpopflag) {
	PopD[currentpop].Taum=atof(s+5);report("  Membrane time constant=%f ms\n",(double)PopD[currentpop].Taum); continue;
      }
      if(strncmp(s,"RestPot=",8)==0 && currentpopflag) {
	PopD[currentpop].RestPot=atof(s+8);report("  Resting potential=%f mV\n",(double)PopD[currentpop].RestPot); continue;
      }
      if(strncmp(s,"ResetPot=",9)==0 && currentpopflag) {
	PopD[currentpop].ResetPot=atof(s+9);report("  Reset potential=%f mV\n",(double)PopD[currentpop].ResetPot); continue;
      }
      if(strncmp(s,"Threshold=",10)==0 && currentpopflag) {
	PopD[currentpop].Threshold=atof(s+10);report("  Threshold =%f mV\n",(double)PopD[currentpop].Threshold); continue;
      }

      if(strncmp(s,"RestPot_ca=",11)==0 && currentpopflag) {
	PopD[currentpop].Vk=atof(s+11);report("  RestPot_ca =%f mV\n",(double)PopD[currentpop].Vk); continue;
      }

      if(strncmp(s,"RefractT=",9)==0 && currentpopflag) {
	PopD[currentpop].RefractT=atof(s+9);report("  Refractory period =%f mV\n",(double)PopD[currentpop].RefractT); continue;
      }

      if(strncmp(s,"Alpha_ca=",9)==0 && currentpopflag) {
	PopD[currentpop].alpha_ca=atof(s+9);report("  Alpha_ca =%f mV\n",(double)PopD[currentpop].alpha_ca); continue;
      }

      if(strncmp(s,"Tau_ca=",7)==0 && currentpopflag) {
	PopD[currentpop].tau_ca=atof(s+7);report("  Tau_ca =%f mV\n",(double)PopD[currentpop].tau_ca); continue;
      }

      if(strncmp(s,"Eff_ca=",7)==0 && currentpopflag) {
	PopD[currentpop].g_ahp=atof(s+7);report("  Eff_ca =%f mV\n",(double)PopD[currentpop].g_ahp); continue;
      }

      if(strncmp(s,"g_ADR_Max=",10)==0 && currentpopflag) {
	PopD[currentpop].g_adr_max=atof(s+10);report("  g_adr_max=%f mV\n",(double)PopD[currentpop].g_adr_max); continue;
      }

      if(strncmp(s,"V_ADR_h=",8)==0 && currentpopflag) {
	PopD[currentpop].Vadr_h=atof(s+8);report("  V_ADR_h=%f mV\n",(double)PopD[currentpop].Vadr_h); continue;
      }

      if(strncmp(s,"V_ADR_s=",8)==0 && currentpopflag) {
	PopD[currentpop].Vadr_s=atof(s+8);report("  V_ADR_s=%f mV\n",(double)PopD[currentpop].Vadr_h); continue;
      }

      if(strncmp(s,"ADRRevPot=",10)==0 && currentpopflag) {
	PopD[currentpop].ADRRevPot=atof(s+10);report("  ADRRevPot=%f mV\n",(double)PopD[currentpop].ADRRevPot); continue;
      }

      if(strncmp(s,"g_k_Max=",8)==0 && currentpopflag) {
	PopD[currentpop].g_k_max=atof(s+8);report("  g_k_max=%f mV\n",(double)PopD[currentpop].g_k_max); continue;
      }

      if(strncmp(s,"V_k_h=",6)==0 && currentpopflag) {
	PopD[currentpop].Vk_h=atof(s+6);report("  V_k_h=%f mV\n",(double)PopD[currentpop].Vk_h); continue;
      }

      if(strncmp(s,"V_k_s=",6)==0 && currentpopflag) {
      	PopD[currentpop].Vk_s=atof(s+6);report("  V_k_s=%f mV\n",(double)PopD[currentpop].Vk_h); continue;
         }
      if(strncmp(s,"tau_k_max=",10)==0 && currentpopflag) {
	PopD[currentpop].tau_k_max=atof(s+10);report("  tau_k_max=%f mV\n",(double)PopD[currentpop].tau_k_max); continue;
      }

      // receptor paramters
      if(strncmp(s,"Receptor:",9)==0 && currentpopflag) {
	currentreceptorflag=1;
	s+=9; while(*s==' ') s++;
	currentreceptor=ReceptorCode(s,currentpop);
	if(currentreceptor==-1) { printf("ERROR: Unknown receptor type\n"); return -1; }
	if(strncmp(s,"GABAB",4)==0) { // deactivate magnesium block for GABAB
	  PopD[currentpop].MgFlag[currentreceptor]=0;
	}
	else PopD[currentpop].MgFlag[currentreceptor]=0;
	report("Receptor %d: %s (Mg block: %d)\n",currentreceptor,PopD[currentpop].ReceptorLabel[currentreceptor],PopD[currentpop].MgFlag[currentreceptor]);

        PopD[currentpop].ExtEffSD[currentreceptor]=0;
        PopD[currentpop].FreqExtSD[currentreceptor]=0;
	continue;
      }

      if(strncmp(s,"EndReceptor",11)==0) {
	currentreceptorflag=0;
	continue;
      }

      if(strncmp(s,"Tau=",4)==0 && currentpopflag && currentreceptorflag) {
	PopD[currentpop].Tau[currentreceptor]=atof(s+4);report("  Tau=%f ms\n",(double)PopD[currentpop].Tau[currentreceptor]); continue;
      }
      if(strncmp(s,"RevPot=",7)==0 && currentpopflag && currentreceptorflag) {
	PopD[currentpop].RevPots[currentreceptor]=atof(s+7);report("  Reversal potential=%f mV\n",(double)PopD[currentpop].RevPots[currentreceptor]); continue;
      }
      if(strncmp(s,"FreqExt=",8)==0 && currentpopflag && currentreceptorflag) {
	PopD[currentpop].FreqExt[currentreceptor]=atof(s+8);report("  Ext frequency=%f Hz\n",(double)PopD[currentpop].FreqExt[currentreceptor]); continue;
      }
      if(strncmp(s,"FreqExtSD=",10)==0 && currentpopflag && currentreceptorflag) {
	PopD[currentpop].FreqExtSD[currentreceptor]=atof(s+10);report("  Ext frequency SD=%f Hz\n",(double)PopD[currentpop].FreqExtSD[currentreceptor]); continue;
      }
      if(strncmp(s,"MeanExtEff=",11)==0 && currentpopflag && currentreceptorflag) {
	PopD[currentpop].MeanExtEff[currentreceptor]=atof(s+11);report("  Ext efficacy=%f nS\n",(double)PopD[currentpop].MeanExtEff[currentreceptor]); continue;}
      if(strncmp(s,"ExtEffSD=",9)==0 && currentpopflag && currentreceptorflag) {
	PopD[currentpop].ExtEffSD[currentreceptor]=atof(s+9);report("  Ext efficacy SD=%f nS\n",(double)PopD[currentpop].ExtEffSD[currentreceptor]); continue;
      }
      if(strncmp(s,"MeanExtCon=",11)==0 && currentpopflag && currentreceptorflag) {
	PopD[currentpop].MeanExtCon[currentreceptor]=atof(s+11);report("  Ext connections=%f\n",(double)PopD[currentpop].MeanExtCon[currentreceptor]); continue;
      }

      // target populations

      if(strncmp(s,"TargetPopulation:",17)==0 && currentpopflag)
	{
	  s+=17; while(*s==' ') s++;
	  currenttarget=PopulationCode(s);
	  if(currenttarget==-1) { printf("Unknown target population: line %d\n",line); return -1;}
	  PopD[currentpop].TargetPops[PopD[currentpop].NTargetPops]=currenttarget;

	  report("Synapses targeting population: %s (%d)\n ",PopD[currenttarget].Label,currenttarget);

          PopD[currentpop].TargetAxons[PopD[currentpop].NTargetPops]=-1;
          PopD[currentpop].SynP[PopD[currentpop].NTargetPops].pv=0;
          PopD[currentpop].SynP[PopD[currentpop].NTargetPops].tauD=1000;
          PopD[currentpop].SynP[PopD[currentpop].NTargetPops].Fp=0;
          PopD[currentpop].SynP[PopD[currentpop].NTargetPops].tauF=5000;
          PopD[currentpop].SynP[PopD[currentpop].NTargetPops].EfficacySD=0;
//          PopD[currentpop].SynP[PopD[currentpop].NTargetPops].Ca_pre_tau=dt;
	  continue;
	}

      if(strncmp(s,"EndTargetPopulation",20)==0)
	{
	  PopD[currentpop].NTargetPops++;
	  continue;
	}

      if(strncmp(s,"TargetAxon=",11)==0 && currentpopflag)
	{
	  s+=11; while(*s==' ') s++;
	  currenttargetaxon=PopulationCode(s);
	  if(currenttarget==-1) { printf("Unknown target population: line %d\n",line); return -1;}
	  PopD[currentpop].TargetAxons[PopD[currentpop].NTargetPops]=currenttargetaxon;
          report("  TargetAxon=%s\n",s);
	}
      if(strncmp(s,"Connectivity=",13)==0 && currentpopflag)
	{
	  aux=atof(s+13);
	  PopD[currentpop].SynP[PopD[currentpop].NTargetPops].Connectivity=aux;report("  Connectivity=%f\n",(double)aux);
	}
      if(strncmp(s,"TargetReceptor=",15)==0 && currentpopflag)
	{
	  s+=15; while(*s==' ' || *s=='\t') s++;
	  auxi=ReceptorCode(s,currenttarget);
	  if(auxi==-1) { printf("Unknown target receptor, line %d\n",line); return -1; }
	  PopD[currentpop].SynP[PopD[currentpop].NTargetPops].TargetReceptor=auxi;
	  report("  Target receptor code=%d\n",auxi);
	}

      if(strncmp(s,"MeanEff=",8)==0 && currentpopflag)
	{
	  aux=atof(s+8);
	  PopD[currentpop].SynP[PopD[currentpop].NTargetPops].MeanEfficacy=aux;report("  Mean efficacy=%f\n",(double)aux);
	}
      if(strncmp(s,"EffSD=",6)==0 && currentpopflag)
	{
	  aux=atof(s+6);
	  PopD[currentpop].SynP[PopD[currentpop].NTargetPops].EfficacySD=aux;report("  Efficacy SD=%f\n",(double)aux);
	}

      if(strncmp(s,"STDepressionP=",14)==0 && currentpopflag) {
        aux=atof(s+14);
	PopD[currentpop].SynP[PopD[currentpop].NTargetPops].pv=aux;
        report("  STDepressionP=%f mV\n",(double)aux);
      }

      if(strncmp(s,"STDepressionTau=",16)==0 && currentpopflag) {
        aux=atof(s+16);
	PopD[currentpop].SynP[PopD[currentpop].NTargetPops].tauD=aux;
        report("  STDepressionTau=%f mV\n",(double)aux);
      }
      if(strncmp(s,"STFacilitationP=",16)==0 && currentpopflag) {
        aux=atof(s+16);
	PopD[currentpop].SynP[PopD[currentpop].NTargetPops].Fp=aux;
        report("  STFacilitationP=%f mV\n",(double)aux);
      }

      if(strncmp(s,"STFacilitationTau=",18)==0 && currentpopflag) {
        aux=atof(s+18);
	PopD[currentpop].SynP[PopD[currentpop].NTargetPops].tauF=aux;
        report("  STFacilitationTau=%f mV\n",(double)aux);
      }
    } // end while


  fclose(devconf);

 	return 1;

}


// -------------------------------------------------------------------------------------
// Parse protocol
// -------------------------------------------------------------------------------------

int ParseProtocol()
{
  FILE *devprot;
  char buf[1000],*s,*es;
  int eventflag=0;
  int line,auxi,currentpop;
  float aux;

  report("-------------------------------------------------\nParsing protocol...\n");

  /*  strncpy=(network_pro,prefix,strlen(prefix));
  strcpy=(network_pro+strlen(prefix),"network.pro");
  devconf=fopen(network_pro,"r");*/
  devprot=fopen("network.pro","r");
  if(devprot==NULL) { printf("ERROR:  Unable to read protocol file\n"); return 0;}

  line=-1; NEvents=0;

  while(fgets(buf,1000,devprot)!=NULL)
    {
      line++;
      s=buf;
      // trim \n at the end
      es=buf; while(*es && *es!='\n') es++; *es=0;

      while(*s==' ' || *s=='\t') s++; // skip blanks
      if(*s==0) continue; // empty line
      if(*s=='%') continue; // remark



      // commands for defining a new event
      if(strncmp(s,"EventTime",9)==0)
	{
	  eventflag=1;
	  s+=9;
	  Events[NEvents].ETime=atof(s);
	  if(Events[NEvents].ETime<0.)
	    {
	      printf("ERROR: Invalid event time, line %d\n",line);
	      return -1;
	    }
          Events[NEvents].FreqExtsd=0;    // initialize FreqExtsd
	  report("Event %d: time %f ms\n",NEvents,Events[NEvents].ETime);
	  continue;
	}

      if(strncmp(s,"EndEvent",8)==0)
	{
	  eventflag=0;
	  NEvents++;
	  continue;
	}

      if(strncmp(s,"Type=",5)==0) {
	s+=5; while(*s==' ' || *s=='\t') s++;
	if(strncmp(s,"ChangeExtFreq",13)==0) {
	  Events[NEvents].Type=CHANGEEXTFREQ;
         Events[NEvents].FreqExtDecay=dt;
	}
	if(strncmp(s,"ChangeExtFreqSD",15)==0) {
	  Events[NEvents].Type=CHANGEEXTFREQSD;
	}
	if(strncmp(s,"ChangeMeanEff",13)==0) {
	  Events[NEvents].Type=CHANGEMEANEFF;
	}
	if(strncmp(s,"EndTrial",8)==0) {
	  Events[NEvents].Type=ENDOFTRIAL;
	  TrialDuration=Events[NEvents].ETime;
	}
	continue;
      }

      if(strncmp(s,"FreqExt=",8)==0) {
	s+=8;
	Events[NEvents].FreqExt=atof(s);
	report("  New external frequency: %f Hz\n",Events[NEvents].FreqExt);
	continue;
      }
      if(strncmp(s,"FreqExtDecay=",13)==0) {
	s+=13;
	Events[NEvents].FreqExtDecay=atof(s);
	report("  Time constant for reaching the new external frequency: %f Hz\n",Events[NEvents].FreqExtDecay);
	continue;
      }
      if(strncmp(s,"FreqExtSD=",10)==0) {
	s+=10;
	Events[NEvents].FreqExtSD=atof(s);
	report("  New external frequency SD: %f Hz\n",Events[NEvents].FreqExtSD);
	continue;
      }
      if(strncmp(s,"FreqExtsd=",10)==0) {
	s+=10;
	Events[NEvents].FreqExtsd=atof(s);
	report("  New external frequency offset: %f Hz\n",Events[NEvents].FreqExtsd);
	continue;
      }
      if(strncmp(s,"MeanEff=",8)==0) {
	s+=8;
	Events[NEvents].MeanEff=atof(s);
	report("  New Mean Eff: %f Hz\n",Events[NEvents].MeanEff);
	continue;
      }

      if(strncmp(s,"Label=",6)==0) {
	s+=6;
	strcpy(Events[NEvents].Label,s);
	report("   Label: [%s]\n",Events[NEvents].Label);
	continue;

      }

      if(strncmp(s,"Population:",11)==0 && eventflag) {
	s+=11;	while(*s==' ') s++;
	currentpop=PopulationCode(s);
	if(currentpop==-1) { printf("ERROR: Unknown population: line %d\n",line); return -1;}
	Events[NEvents].PopNumber=currentpop;
	report("  Population code: %d\n",currentpop);
	continue;
      }

      if(strncmp(s,"TargetPopulation:",17)==0 && eventflag) {
	s+=17;	while(*s==' ') s++;
	currentpop=PopulationCode(s);
	if(currentpop==-1) { printf("ERROR: Unknown population: line %d\n",line); return -1;}
	Events[NEvents].TargetPopNumber=currentpop;
	report("  TargetPopulation code: %d\n",currentpop);
	continue;
      }

      if(strncmp(s,"Receptor:",9)==0 && eventflag)
	{
	  s+=9; while(*s==' ' || *s=='\t') s++;
	  auxi=ReceptorCode(s,currentpop);
	  if(auxi==-1) { printf("ERROR: Unknown receptor, line %d\n",line); return -1; }
	  Events[NEvents].ReceptorNumber=auxi;
	  report("  Receptor code:%d\n",auxi);

	}

    } // end while

  // sort events!... (for now I rely on the fact that events are sorted)

  return 1;

}


// ======================================================================================
// GenerateNetwork()
//
// Generates all the structures which are needed to run the simulation
// on the basis of PopD structures. It also allocates the memory for the
// network. PopD is initialized in DescribeNetwork.
// ======================================================================================

int GenerateNetwork()
{

  int p,i,j,r,tp,tpi,tn,k,rd,na,tpa,flag_a;
  int ni[MAXTN],tpni[MAXTN],trni[MAXTN], tpax[MAXTN];
  float eff[MAXTN],ts[MAXTN],timets[MAXTN],D[MAXTN],pv[MAXTN],tauD[MAXTN],F[MAXTN],Fp[MAXTN],tauF[MAXTN];
  float efficacy;

  if(flagverbose) {
    report("\nGenerating network...\n");
  }


  delayindt=1;

  // loop over all the populations
  for(p=0;p<Npop;p++)
    {
      strcpy(Pop[p].Label,PopD[p].Label);
      if(flagverbose) {
	report("Population %2d: %s\n",p,Pop[p].Label);
      }

      for(rd=0;rd<PopD[p].Nreceptors;rd++)
	{
	  // external input: compute first the asymptoic mu and sigma (m,s in nphys notation) of the conductances
	  //          do{efficacy=(gasdev()*PopD[p].ExtEffSD[rd])+PopD[p].MeanExtEff[rd];}while(efficacy<0);
          efficacy=PopD[p].MeanExtEff[rd];
	  PopD[p].MeanExtMuS[rd]=PopD[p].FreqExt[rd]*.001*efficacy*PopD[p].MeanExtCon[rd]*PopD[p].Tau[rd];
	  PopD[p].MeanExtSigmaS[rd]=sqrt(PopD[p].Tau[rd]*.5*PopD[p].FreqExt[rd]*.001*efficacy*efficacy*PopD[p].MeanExtCon[rd]);
          PopD[p].FreqExtTau[rd]=dt;
          PopD[p].FreqExtDecay[rd]=0.999;
          PopD[p].FreqExtNorm[rd]=PopD[p].FreqExtDecay[rd]/(1-PopD[p].FreqExtDecay[rd]);
          PopD[p].FreqExtMem[rd]=PopD[p].MeanExtEff[rd]*PopD[p].FreqExtNorm[rd];
	  report("Pop %d Conductance: %d m=%f nS s=%f nS\n",p,rd,PopD[p].MeanExtMuS[rd],PopD[p].MeanExtSigmaS[rd]);
	}


      // allocate memory for all neurons
      Pop[p].Ncells=PopD[p].Ncells;
      if(flagverbose) { report("  - Number of cells: %d\n",Pop[p].Ncells); }
      Pop[p].Cell=(Neuron *) calloc(Pop[p].Ncells,sizeof(Neuron));

      // init auxiliary variables like the tables of spikes
      if(delayindt>MAXDELAYINDT) { printf("ERROR: delay too long. Increase MAXDELAYINDT and recompile\n"); return -1;}
      Pop[p].CTableofSpikes=delayindt;
      Pop[p].DTableofSpikes=0;
      for(k=0;k<MAXDELAYINDT;k++)
	{
	  Pop[p].NTableofSpikes[k]=0;
	}

      // generate neurons and their axonal trees
      for(i=0;i<Pop[p].Ncells;i++)
	{
	  // single neuron parameters
	  Pop[p].Cell[i].Taum=PopD[p].Taum;
	  Pop[p].Cell[i].ResetPot=PopD[p].ResetPot;
	  Pop[p].Cell[i].C=PopD[p].C;
	  Pop[p].Cell[i].RestPot=PopD[p].RestPot;
	  Pop[p].Cell[i].Threshold=PopD[p].Threshold;
	  Pop[p].Cell[i].RefractT=PopD[p].RefractT;
	  Pop[p].Cell[i].Vk=PopD[p].Vk;
	  Pop[p].Cell[i].alpha_ca=PopD[p].alpha_ca;
	  Pop[p].Cell[i].tau_ca=PopD[p].tau_ca;
	  Pop[p].Cell[i].g_ahp=PopD[p].g_ahp;
          Pop[p].Cell[i].g_adr_max=PopD[p].g_adr_max;
          Pop[p].Cell[i].Vadr_h=PopD[p].Vadr_h;
          Pop[p].Cell[i].Vadr_s=PopD[p].Vadr_s;
          Pop[p].Cell[i].ADRRevPot=PopD[p].ADRRevPot;
          Pop[p].Cell[i].g_k_max=PopD[p].g_k_max;
          Pop[p].Cell[i].Vk_h=PopD[p].Vk_h;
          Pop[p].Cell[i].Vk_s=PopD[p].Vk_s;
          Pop[p].Cell[i].tau_k_max=PopD[p].tau_k_max;
          Pop[p].Cell[i].n_k=0;


	  // receptors
	  Pop[p].Cell[i].Nreceptors=PopD[p].Nreceptors;
	  for(r=0;r<Pop[p].Cell[i].Nreceptors;r++)
	    {
	      // parameters (simply copy if the network is homoegeneous)
	      Pop[p].Cell[i].Tau[r]=PopD[p].Tau[r];
	      Pop[p].Cell[i].RevPots[r]=PopD[p].RevPots[r];
	      Pop[p].Cell[i].MgFlag[r]=PopD[p].MgFlag[r]; // magnesium block
	      // init
	      Pop[p].Cell[i].LS[r]=0.;

	      // external input
	      do{efficacy=(gasdev()*PopD[p].ExtEffSD[r])+PopD[p].MeanExtEff[r];}while(efficacy<0);
  	      Pop[p].Cell[i].ExtMuS[r]=PopD[p].FreqExt[r]*.001*efficacy*PopD[p].MeanExtCon[r]*PopD[p].Tau[r];
	      Pop[p].Cell[i].ExtSigmaS[r]=sqrt(PopD[p].Tau[r]*.5*PopD[p].FreqExt[r]*.001*efficacy*efficacy*PopD[p].MeanExtCon[r]);
	      //if(r==0)printf("%f\n",Pop[p].Cell[i].ExtMuS[r]);
	      Pop[p].Cell[i].ExtS[r]=0.;

	    } // end for r

	  // init
	  Pop[p].Cell[i].V=Pop[p].Cell[i].RestPot;
	  Pop[p].Cell[i].RefrState=0;
	  Pop[p].Cell[i].TimeLastSpike=-10000.; // time of the last emitted spike (for NMDA saturation)
	  Pop[p].Cell[i].PTimeLastSpike=-10000.; // time of spike emitted previous the last emitted spiek (for NMDA saturation)

	  // Axonal tree -------------------------------
	  // loop on target populations
	  Pop[p].Cell[i].Naxonals=0;

     	  for(tpi=0;tpi<PopD[p].NTargetPops;tpi++)
	    {
	      tp=PopD[p].TargetPops[tpi]; // target population (tpi=target pop index)
	      // loop on the neurons of the target population

	      // choose the target neurons/populations and put them in a temporary vector ni/tpni
	      na=Pop[p].Cell[i].Naxonals; // auxiliary variable to speed up

              // targeted axon (for pre-synaptic inhibition)
              tpa=PopD[p].TargetAxons[tpi];


              // Source and target populations have to have the same number of neurons to be
              // connected one-to-one
              if(PopD[tp].Ncells!=PopD[p].Ncells){printf("Error: Source and target populations do not have the same number of neurons. Cannot make one-to-one connection. Exit..."); exit(0);}
              //Neuron cannot connect to itself
              if(tp==p){printf("Error: Neuron cannot connect to itself. Exit..."); exit(0);}

              // To make one-to-one connection, connectivity should be equal to -1
              if(PopD[p].SynP[tpi].Connectivity!=-1)
                {printf("Error: To make one-to-one connection, the value of connectivity should be set to -1. Exit..."); exit(0);}


//	      for(tn=0;tn<PopD[tp].Ncells;tn++)
//		{


//		  if(drand49()<PopD[p].SynP[tpi].Connectivity) { // the connection exists
		    ni[na]=i;   // target neuron
		    tpni[na]=tp; // target population
                    tpax[na]=tpa;
		    trni[na]=PopD[p].SynP[tpi].TargetReceptor; // target receptor
                    do{efficacy=gasdev()*PopD[p].SynP[tpi].EfficacySD+PopD[p].SynP[tpi].MeanEfficacy;}while(efficacy<0);
		    eff[na]= efficacy;// efficacy
		    // printf("%f\n",efficacy);
                    D[na]=1;
                    pv[na]=PopD[p].SynP[tpi].pv;
                    tauD[na]=PopD[p].SynP[tpi].tauD;
                    Fp[na]=PopD[p].SynP[tpi].Fp;
                    if(Fp[na]==0)
                       F[na]=1;
                    else
                       F[na]=0;
                    tauF[na]=PopD[p].SynP[tpi].tauF;

		    if(strncmp(PopD[tp].ReceptorLabel[trni[na]],"NMDA",4)==0) {
		      ts[na]=0.; // initial LastConductance (for NMDA saturation)
		    }
		    else ts[na]=-1.; // not an NMDA type (so, no saturation)
		    na++;
		    Pop[p].Cell[i].Naxonals++;
//		  }
//		} // end for tn
	    } // enf for tpi
	  //
	  Pop[p].Cell[i].Axonals=(Synapse *) calloc(Pop[p].Cell[i].Naxonals,sizeof(Synapse));

	  for(j=0;j<Pop[p].Cell[i].Naxonals;j++)
	    {
	      Pop[p].Cell[i].Axonals[j].TargetPop=tpni[j];
	      Pop[p].Cell[i].Axonals[j].TargetNeuron=ni[j];
              Pop[p].Cell[i].Axonals[j].TargetAxon=tpax[j];
              // check if the targeted axon exists
              if(tpax[j]!=-1){

                // check efficacy value
                if(eff[j]>1 || eff[j]<0){printf("For presynaptic inhibition, the efficacy should be between 0 and 1. Stop! \n"); exit(0);}

                flag_a=0;
                for(k=0;k<Pop[tpni[j]].Cell[i].Naxonals;k++){
                   if(PopD[tpni[j]].TargetPops[k]==tpax[j]){  // need to use PopD becuase the neuron may have not yet created.
                        Pop[p].Cell[i].Axonals[j].TargetAxon=k;
                        flag_a=1;
                     }
                 if(flag_a==1) break;
                  }

              }
		  Pop[p].Cell[i].Axonals[j].Ca_pre=1.0;
	      Pop[p].Cell[i].Axonals[j].Efficacy=eff[j];
	      Pop[p].Cell[i].Axonals[j].TargetReceptor=trni[j];
	      Pop[p].Cell[i].Axonals[j].D=D[j];
	      Pop[p].Cell[i].Axonals[j].pv=pv[j];
              Pop[p].Cell[i].Axonals[j].tauD=tauD[j];
	      Pop[p].Cell[i].Axonals[j].F=F[j];
	      Pop[p].Cell[i].Axonals[j].Fp=Fp[j];
              Pop[p].Cell[i].Axonals[j].tauF=tauF[j];
	      Pop[p].Cell[i].Axonals[j].LastConductance=ts[j];

              for(k=0;k<PopD[p].Nreceptors;k++){
                Pop[p].Cell[i].Axonals[j].St_pre[k]=0;
                Pop[p].Cell[i].Axonals[j].last_pre_update[k]=0;
                Pop[p].Cell[i].Axonals[j].Ca_pre_tau[k]=Pop[p].Cell[i].Tau[k];
              }
	    }
	} // end for i
    } // end for p


  if(flagverbose) { report("Network generated\n"); fflush(stdout); }
}


// same as generate network, but it does not allocate memory. Used in the multitrial mode
// (to change the network from trial to trial, the only way is to start again with a different seed)


int InitializeNetwork()
{

  int p,i,j,r,tp,tpi,tn,k,rd,na;
  int ni[MAXTN],tpni[MAXTN],trni[MAXTN];
  float eff[MAXTN],ts[MAXTN],timets[MAXTN];

  if(flagverbose) {
    report("\nGenerating network...\n");
  }

  dt=.1;
  delayindt=5;

  // loop over all the populations
  for(p=0;p<Npop;p++)
    {
      strcpy(Pop[p].Label,PopD[p].Label);
      if(flagverbose) {
	report("Population %2d: %s\n",p,Pop[p].Label);
      }

      for(rd=0;rd<PopD[p].Nreceptors;rd++)
	{
	  // external input: compute first the asymptoic mu and sigma (m,s in nphys notation) of the conductances
	  // CAREFUL: PROTOCOL FILE SHOULD BE REREAD TO INITIALIZE CORRECTLY EXT FREQS (THEY MIGHT HAVE BEEN CHANGED IN THE PREVIOUS TRIAL

	  PopD[p].MeanExtMuS[rd]=PopD[p].FreqExt[rd]*.001*PopD[p].MeanExtEff[rd]*PopD[p].MeanExtCon[rd]*PopD[p].Tau[rd];
	  PopD[p].MeanExtSigmaS[rd]=sqrt(PopD[p].Tau[rd]*.5*PopD[p].FreqExt[rd]*.001*PopD[p].MeanExtEff[rd]*PopD[p].MeanExtEff[rd]*PopD[p].MeanExtCon[rd]);
	  report("Pop %d Conductance: %d m=%f nS s=%f nS\n",p,rd,PopD[p].MeanExtMuS[rd],PopD[p].MeanExtSigmaS[rd]);
	}

      // init auxiliary variables like the tables of spikes
      if(delayindt>MAXDELAYINDT) { printf("ERROR: delay too long. Increase MAXDELAYINDT and recompile\n"); return -1;}
      Pop[p].CTableofSpikes=delayindt;
      Pop[p].DTableofSpikes=0;
      for(k=0;k<MAXDELAYINDT;k++)
	{
	  Pop[p].NTableofSpikes[k]=0;
	}

      // init     neurons and their axonal trees
      for(i=0;i<Pop[p].Ncells;i++)
	{
	  for(r=0;r<Pop[p].Cell[i].Nreceptors;r++)
	    {
	      Pop[p].Cell[i].LS[r]=0.;

	      // external input
	      Pop[p].Cell[i].ExtMuS[r]=PopD[p].MeanExtMuS[r];
	      Pop[p].Cell[i].ExtSigmaS[r]=PopD[p].MeanExtSigmaS[r];
	      Pop[p].Cell[i].ExtS[r]=0.;
	    } // end for r

	  // init
	  Pop[p].Cell[i].V=Pop[p].Cell[i].RestPot;
          Pop[p].Cell[i].Ca=0;
	  Pop[p].Cell[i].RefrState=0;
	  Pop[p].Cell[i].TimeLastSpike=-10000.; // time of the last emitted spike (for NMDA saturation)
	  Pop[p].Cell[i].PTimeLastSpike=-10000.; // time of spike emitted previous the last emitted spiek (for NMDA saturation)

	} // end for i
    } // end for p

  if(flagverbose) { report("Network initialized\n"); fflush(stdout); }
}



// ---------------------------------------------------------------------------------------------------

/*

Trick for NMDA:

$$ {ds_k \over dt} = -{s_k \over \tau} +\alpha(1-s_k)\sum_j \delta(t-t^k_j)$$

Multiply by $w_k$ and sum:

$$ {dS \over dt} = -S \over \tau+\sum_k \alpha (1-s_k)w_k \delta(t-t^k_j)$$

*/

// float current_freq;

int SimulateOneTimeStep()
{
  int aux;
  int p,i,j,r,sourceneuron,k;
  int tn,tp,tr,ta,tpi, tni, tri;
  float s,saturationfactor,total_St_pre;
  float Vaux; // auxiliary V: during the emission of the spike V is set artificially to 0. This is bad for the reversal potential
  float D,F,Ca_pre,DT;
  float g_adr, g_k, tau_max, alpha, beta, dv, n, tau_n, n_inif, efficacy, ExtMuS, ExtSigmaS;
  static float freq[MAXP][MAXRECEPTORS], mean_freq[MAXP][MAXRECEPTORS];
  static float last_freq[MAXP][MAXRECEPTORS];
  static int flag=0;
  static

  // DEBUG
  float nvalues=0,value=0;

  // END DEBUG

  if(flag==0){
    for(j=0;j<MAXP;j++){
    for(i=0;i<MAXRECEPTORS;i++)
      freq[j][i]=PopD[j].FreqExt[i];
      mean_freq[j][i]=PopD[j].FreqExt[i];
      last_freq[j][i]=PopD[j].FreqExt[i];
    }
   flag=1;
  }


  // Compute the decay of the total conductances and add external input
  // ------------------------------------------------------------------
  for(p=0;p<Npop;p++)
    {

      for(r=0;r<MAXRECEPTORS;r++){
	if(PopD[p].FreqExtSD[r]!=0){
	    // do{freq[p][r]=PopD[p].FreqExt[r]+gasdev()*PopD[p].FreqExtSD[r];}while(freq[p][r]<0);
	  do{freq[p][r]+=-(1-PopD[p].FreqExtDecay[r])*(last_freq[p][r]-PopD[p].FreqExt[r])
             +gasdev()*PopD[p].FreqExtSD[r];}while(freq[p][r]<0);
           last_freq[p][r]=freq[p][r]; ext_freq=freq[1][0];

	}
        else{
           if(PopD[p].FreqExtTau[r]==dt)
             {
//printf("cp1 PopD[%i].FreqExt[%i]=%f, FreqExtTau=%f \n",p,r,PopD[p].FreqExt[r],PopD[p].FreqExtTau[r]);
mean_freq[p][r]=PopD[p].FreqExt[r];}
           else{
//printf("cp2 PopD[%i].FreqExt[%i]=%f, FreqExtTau=%f, FreqExtDecay=%f \n",p,r,PopD[p].FreqExt[r],PopD[p].FreqExtTau[r],PopD[p].FreqExtDecay[r]);
             mean_freq[p][r]=(mean_freq[p][r]-PopD[p].FreqExt[r])*(1-dt/PopD[p].FreqExtTau[r])+PopD[p].FreqExt[r];}
             freq[p][r]=mean_freq[p][r]; ext_freq=freq[1][0];}
      }
      //      current_freq=freq[0][0];
      for(i=0;i<Pop[p].Ncells;i++)
	{
	  for(r=0;r<Pop[p].Cell[i].Nreceptors;r++)
	    {
              efficacy=PopD[p].MeanExtEff[r];
              Pop[p].Cell[i].ExtMuS[r]=freq[p][r]*.001*efficacy*PopD[p].MeanExtCon[r]*Pop[p].Cell[i].Tau[r];
              Pop[p].Cell[i].ExtSigmaS[r]=sqrt(Pop[p].Cell[i].Tau[r]*.5*freq[p][r]*.001*efficacy*efficacy*PopD[p].MeanExtCon[r]);
   	      s=Pop[p].Cell[i].ExtSigmaS[r];
	      if(s!=0.) // to optimize
		{
		  Pop[p].Cell[i].ExtS[r]+=dt/Pop[p].Cell[i].Tau[r]*(-Pop[p].Cell[i].ExtS[r]+Pop[p].Cell[i].ExtMuS[r])
		+s*sqrt(dt*2./Pop[p].Cell[i].Tau[r])*gauss();
		}
	      else {
		Pop[p].Cell[i].ExtS[r]+=dt/Pop[p].Cell[i].Tau[r]*(-Pop[p].Cell[i].ExtS[r]+Pop[p].Cell[i].ExtMuS[r]);
	      }

	      Pop[p].Cell[i].LS[r]*=exp(-dt/Pop[p].Cell[i].Tau[r]); // decay
	    }
	}
    }

  // Update the total conductances (changes provoked by the spikes)
  // --------------------------------------------------------------
  for(p=0;p<Npop;p++)
    {
      // loop over all the spikes emitted at time t-delay (they are received now)
      for(i=0;i<Pop[p].NTableofSpikes[Pop[p].DTableofSpikes];i++)
	{
	  sourceneuron=Pop[p].TableofSpikes[Pop[p].DTableofSpikes][i];

	  // for each spike, loop over the target conductances
	  for(j=0;j<Pop[p].Cell[sourceneuron].Naxonals;j++)
	    {
            if((ta=Pop[p].Cell[sourceneuron].Axonals[j].TargetAxon)!=-1){ //if this is a presynaptic inhibition
                tpi=Pop[p].Cell[sourceneuron].Axonals[j].TargetPop;
                tni=Pop[p].Cell[sourceneuron].Axonals[j].TargetNeuron;
                tri=Pop[p].Cell[sourceneuron].Axonals[j].TargetReceptor;
                if(DT=(Time-Pop[tpi].Cell[tni].Axonals[ta].last_pre_update[tri])>0){  // need to update presynaptic gating variable first
                   Pop[tpi].Cell[tni].Axonals[ta].St_pre[tri]=Pop[tpi].Cell[tni].Axonals[ta].St_pre[tri]*exp(-DT/Pop[tpi].Cell[tni].Axonals[ta].Ca_pre_tau[tri]);
                   Pop[tpi].Cell[tni].Axonals[ta].last_pre_update[tri]=Time;
                   }
                 // Now add the effect of presynaptic spikes
                 Pop[tpi].Cell[tni].Axonals[ta].St_pre[tri]+=Pop[p].Cell[sourceneuron].Axonals[j].Efficacy;
                 if(Pop[tpi].Cell[tni].Axonals[ta].St_pre[tri]>1) Pop[tpi].Cell[tni].Axonals[ta].St_pre[tri]=1.0;

              }
           else{ // Not a presynaptic inhibition, also need to update S_pre
               for(k=0;k<PopD[p].Nreceptors;k++){
                 if((strcmp(PopD[p].ReceptorLabel[k],"GABAA")==0)||(strcmp(PopD[p].ReceptorLabel[k],"GABAB")==0)){ // Only update GABAA and GABAB
                 DT=(Time-Pop[p].Cell[sourceneuron].Axonals[j].last_pre_update[k]);
                 if(DT>0&&(Pop[p].Cell[sourceneuron].Axonals[j].St_pre[k]>0)){
                    Pop[p].Cell[sourceneuron].Axonals[j].St_pre[k]=Pop[p].Cell[sourceneuron].Axonals[j].St_pre[k]*exp(-DT/Pop[p].Cell[sourceneuron].Axonals[j].Ca_pre_tau[k]);
                    Pop[p].Cell[sourceneuron].Axonals[j].last_pre_update[k]=Time;
                      }
                  }
                 }

             // Short-term depression: Calculate the recovery of D first.

               Pop[p].Cell[sourceneuron].Axonals[j].D = 1-(1-Pop[p].Cell[sourceneuron].Axonals[j].D)*exp(-(Time-Pop[p].Cell[sourceneuron].PTimeLastSpike)/Pop[p].Cell[sourceneuron].Axonals[j].tauD);

	        D=Pop[p].Cell[sourceneuron].Axonals[j].D;

               total_St_pre=0;
               for(k=1;k<PopD[p].Nreceptors;k++){
                  total_St_pre+=Pop[p].Cell[sourceneuron].Axonals[j].St_pre[k];
                 }
              Ca_pre=1-total_St_pre;
             //D=1.0;
	     // if(sourceneuron==20) printf(" %.1f %.1f .1%f .3%f \n",Time, Pop[p].Cell[sourceneuron].PTimeLastSpike, Pop[p].Cell[sourceneuron].tauD, Pop[p].Cell[sourceneuron].D);

             // Short-term facilitation: Calculate the F value
              if(Pop[p].Cell[sourceneuron].Axonals[j].Fp==0)
                F=1;
              else{
               Pop[p].Cell[sourceneuron].Axonals[j].F = Pop[p].Cell[sourceneuron].Axonals[j].F*exp(-(Time-Pop[p].Cell[sourceneuron].PTimeLastSpike)/Pop[p].Cell[sourceneuron].Axonals[j].tauF);
               F=Pop[p].Cell[sourceneuron].Axonals[j].F;
	        }
	       tn=Pop[p].Cell[sourceneuron].Axonals[j].TargetNeuron;
	       tp=Pop[p].Cell[sourceneuron].Axonals[j].TargetPop;
	       tr=Pop[p].Cell[sourceneuron].Axonals[j].TargetReceptor;

	       if(Pop[p].Cell[sourceneuron].Axonals[j].LastConductance<0.) { // NO NMDA (no saturation)
		  Pop[tp].Cell[tn].LS[tr]+=powf(Ca_pre,3.5)*D*F*Pop[p].Cell[sourceneuron].Axonals[j].Efficacy; // no saturation
	        }
	       else {

		// Now it should be correct. ALPHA factor to be determined (jump for every spike): now it is the maximum of a single spike
		// TEMPORARY
#define ALPHA 0.6332
		// 0.6332 is the best value for the area at 55 Hz. Best value depends on the frequency!!!
		 Pop[p].Cell[sourceneuron].Axonals[j].LastConductance*=exp(-(Time-Pop[p].Cell[sourceneuron].PTimeLastSpike)/Pop[tp].Cell[tn].Tau[tr]);

		 Pop[tp].Cell[tn].LS[tr]+=D*F*Pop[p].Cell[sourceneuron].Axonals[j].Efficacy*(1.-Pop[p].Cell[sourceneuron].Axonals[j].LastConductance)*ALPHA;

		 Pop[p].Cell[sourceneuron].Axonals[j].LastConductance+=ALPHA*(1.-Pop[p].Cell[sourceneuron].Axonals[j].LastConductance);

	       }
 // Calculate the change of D and F due to the spike.
                Pop[p].Cell[sourceneuron].Axonals[j].F += Pop[p].Cell[sourceneuron].Axonals[j].Fp*(1-Pop[p].Cell[sourceneuron].Axonals[j].F);

               Pop[p].Cell[sourceneuron].Axonals[j].D -= Pop[p].Cell[sourceneuron].Axonals[j].pv*powf(Ca_pre,3.5)*F*Pop[p].Cell[sourceneuron].Axonals[j].D;
                if(Pop[p].Cell[sourceneuron].Axonals[j].D<0)  Pop[p].Cell[sourceneuron].Axonals[j].D=0;


              }
	    } // for j
		//                if(sourceneuron==20) printf("   %f\n",Pop[p].Cell[sourceneuron].D);
                //if(Pop[p].Cell[sourceneuron].D<=0) Pop[p].Cell[sourceneuron].D=0;
	}
    }

  // Update the neuronal variables
  // -----------------------------

  for(p=0;p<Npop;p++)
    {
      Pop[p].NTableofSpikes[Pop[p].CTableofSpikes]=0; // reset the number of emitted spikes for pop p
      for(i=0;i<Pop[p].Ncells;i++)
	{

	  if(Pop[p].Cell[i].V>Pop[p].Cell[i].Threshold)
	    {
	      Pop[p].Cell[i].V=Pop[p].Cell[i].ResetPot; // special state after emission
	      Pop[p].Cell[i].RefrState--;
	      continue;
	    }

	  // refractory period

	  if(Pop[p].Cell[i].RefrState) {
	    Pop[p].Cell[i].RefrState--;
	    continue;
	  }

          //Anomalous delayed rectiflier (ADR)
          if(Pop[p].Cell[i].g_adr_max==0) // No ADR
             g_adr=0;
          else //with ADR
	     g_adr = Pop[p].Cell[i].g_adr_max/(1+exp((Pop[p].Cell[i].V-Pop[p].Cell[i].Vadr_h)/Pop[p].Cell[i].Vadr_s));

          //potassium outward rectifying current
          if(Pop[p].Cell[i].g_k_max==0) // No outward current
             g_k=0;
          else{ // with outward current
	    //  g_k = Pop[p].Cell[i].g_k_max/(1+exp((-Pop[p].Cell[i].V+Pop[p].Cell[i].Vk_h)/Pop[p].Cell[i].Vk_s));
            tau_max=Pop[p].Cell[i].tau_k_max;
            dv=(Pop[p].Cell[i].V+55.0);
	    //  alpha=(-10.0/tau_max)*(dv-49)/(1-exp(-(dv-49)/100));
	    //            beta=(0.17/tau_max)*exp(-(dv-11)/10);
            tau_n=tau_max/(exp(-1*dv/30)+exp(dv/30));
            n_inif=1/(1+exp(-(Pop[p].Cell[i].V-Pop[p].Cell[i].Vk_h)/Pop[p].Cell[i].Vk_s));

            Pop[p].Cell[i].n_k+=-1/tau_n*(Pop[p].Cell[i].n_k-n_inif);
            g_k=Pop[p].Cell[i].g_k_max*Pop[p].Cell[i].n_k;
	    //	    printf("v=%f dv=%f alpha=%f beta=%f g_k=%f \n",Pop[p].Cell[i].V,dv,alpha,beta,g_k);
	  }


          // decay
          if(Pop[p].Cell[i].g_ahp==0) //No spike-frequency adaptation
    	    Pop[p].Cell[i].V+= -dt*(1/Pop[p].Cell[i].Taum*(Pop[p].Cell[i].V-Pop[p].Cell[i].RestPot)
                               + g_adr/Pop[p].Cell[i].C*(Pop[p].Cell[i].V-Pop[p].Cell[i].ADRRevPot)
                               + g_k/Pop[p].Cell[i].C*(Pop[p].Cell[i].V-Pop[p].Cell[i].ADRRevPot));
          else // with spike-frequency adaptation (the factor 1/1000 is needed to convert nS/nF to 1/ms)
     	    Pop[p].Cell[i].V+= -dt*(1/Pop[p].Cell[i].Taum *(Pop[p].Cell[i].V-Pop[p].Cell[i].RestPot)
                + Pop[p].Cell[i].Ca*Pop[p].Cell[i].g_ahp/Pop[p].Cell[i].C*0.001*(Pop[p].Cell[i].V-Pop[p].Cell[i].Vk)
                + g_adr/Pop[p].Cell[i].C*(Pop[p].Cell[i].V-Pop[p].Cell[i].ADRRevPot)
                + g_k/Pop[p].Cell[i].C*(Pop[p].Cell[i].V-Pop[p].Cell[i].ADRRevPot));


          // [Ca] decay -- 1st order approximation
          Pop[p].Cell[i].Ca += -Pop[p].Cell[i].Ca*dt/Pop[p].Cell[i].tau_ca;

	  Vaux=Pop[p].Cell[i].V;
      	  if(Vaux>Pop[p].Cell[i].Threshold) Vaux=Pop[p].Cell[i].Threshold;

	  // now add the synaptic currents (the factor 1/1000 is needed to convert nS/nF to 1/ms)
	  for(r=0;r<Pop[p].Cell[i].Nreceptors;r++)
	    {
	      if(Pop[p].Cell[i].MgFlag[r]) { // magnesium block
		Pop[p].Cell[i].V+=dt*(Pop[p].Cell[i].RevPots[r]-Vaux)*.001*(Pop[p].Cell[i].LS[r]+Pop[p].Cell[i].ExtS[r])
		  /Pop[p].Cell[i].C/(1.+exp(-0.062*Vaux)/3.57);

		/*		// DEBUG DEBUG DEBUG
		if(p==1 && i==0) {
		  printf("[%f]",Pop[p].Cell[i].LS[r]);
		}
		// end DEBUG */

	      } else {
	      Pop[p].Cell[i].V+=dt*(Pop[p].Cell[i].RevPots[r]-Vaux)*.001*(Pop[p].Cell[i].LS[r]+Pop[p].Cell[i].ExtS[r])
              /Pop[p].Cell[i].C;
	      }
	    }
	  // spike condition
	  if(Pop[p].Cell[i].V>Pop[p].Cell[i].Threshold) {

	    // a spike is emitted
	    Pop[p].TableofSpikes[Pop[p].CTableofSpikes][Pop[p].NTableofSpikes[Pop[p].CTableofSpikes]]=i;
	    if(Pop[p].NTableofSpikes[Pop[p].CTableofSpikes]<MAXSPIKESINDT-1) Pop[p].NTableofSpikes[Pop[p].CTableofSpikes]++;
	    else { printf("\nERROR: too many spikes in a dt (change MAXSPIKESINDT and recompile)\n"); fflush(stdout); }

	    Pop[p].Cell[i].V=Pop[p].Cell[i].ResetPot;
	    Pop[p].Cell[i].PTimeLastSpike=Pop[p].Cell[i].TimeLastSpike;
	    Pop[p].Cell[i].TimeLastSpike=Time;

	    Pop[p].Cell[i].V=0.; // spike! (temporary)
	    Pop[p].Cell[i].RefrState=(int)floor(Pop[p].Cell[i].RefractT/dt); // refractory period!!! (temporary)

            // [Ca] increases;
            Pop[p].Cell[i].Ca += Pop[p].Cell[i].alpha_ca;

	  }
	}
    }

  // Update the pointers of the table of spikes
  // ---------------------------------------------------------------

  for(p=0;p<Npop;p++)
    {
      Pop[p].CTableofSpikes++; if(Pop[p].CTableofSpikes>MAXDELAYINDT-1) Pop[p].CTableofSpikes=0;
      Pop[p].DTableofSpikes++; if(Pop[p].DTableofSpikes>MAXDELAYINDT-1) Pop[p].DTableofSpikes=0;
    }
}

// DATA ANALYSIS
// --------------------------------------------------------------------------------------------

// A matlab script is generated at the beginning, just before starting the simulation.

#define TIMEWINDOWFORFREQ 20.  // time window on which the mean pop frequency is estimated (in ms)
#define RASTERPLOTNEURONS 50   // number of neurons in the raster plot
#define NUMBEROFTRACES 2  // number of visualized traces of V
#define STEPSFORPRINTIGFREQS 500 // every ... step, the frequencies are printed
#define STEPSFORSAVINGFREQS 10 // every ... step the frequencies are saved (not always, to save space)
#define STEPSFORFLUSHING 80
#define SPIKEBUFFER 300 //szie of the array SpikeBuffer[] in SaveSpikes(). This value should be no
                        // smaller than TIMEWINDOWFORFREQ/dt

int NumberofTraces=5;


int GenMatLabFile()
{
  int i,j,nsp,r;
  FILE *devmatlab;
  static char color[5]="rgby";

  devmatlab=fopen("mrealsimu.m","w");
  if(devmatlab==NULL) { printf("ERROR: Unable to open mathlab file\n"); return 0;}

  fprintf(devmatlab,"clear all\n");
  fprintf(devmatlab,"load popfreqs.dat\n");
  if(NumberofTraces>0) {
    fprintf(devmatlab,"load poptraces.dat\n");
  }

  nsp=2+2*NumberofTraces; // total number of subplots

  for(i=0;i<Npop;i++)
    {
      // plot the raster
      fprintf(devmatlab,"figure(%d)\nclf\nsubplot(%d,1,1);\n",i+1,nsp);
      fprintf(devmatlab,"load pop%d.dat\n",i+1);
      fprintf(devmatlab,"plotRast(pop%d,0.,%f,%f,%d);\n",i+1,TrialDuration,dt,RASTERPLOTNEURONS);
      fprintf(devmatlab,"title('%s - Raster plot');\n",Pop[i].Label);
      fprintf(devmatlab,"xlabel('Time (ms)');\n");

      // plot the mean frequencies
      fprintf(devmatlab,"subplot(%d,1,2);\nplot(popfreqs(:,1),popfreqs(:,%d));\n",nsp,i+2);
      fprintf(devmatlab,"ylabel('Mean spike freq (Hz)');\nxlabel('Time (ms)');\n");
      fprintf(devmatlab,"set(gca,'XLim',[0 %f]);\n",TrialDuration);

      // plot the traces
      for(j=0;j<NumberofTraces;j++)
	{
	  fprintf(devmatlab,"subplot(%d,1,%d);\n",nsp,3+j*2);
	  fprintf(devmatlab,"plot(poptraces(:,1),poptraces(:,%d));\n",2+j*(Pop[i].Cell[j].Nreceptors+1)+i*(NumberofTraces*(Pop[i].Cell[j].Nreceptors+1)));

	  fprintf(devmatlab,"ylabel('Membrane pot (mV) - Neuron %d')\n",j);
	  fprintf(devmatlab,"set(gca,'YLim',[%f 5]);\nset(gca,'XLim',[0 %f]);\n",PopD[i].RestPot,TrialDuration);
	  fprintf(devmatlab,"subplot(%d,1,%d);\n",nsp,3+j*2+1);

	  for(r=0;r<Pop[i].Cell[j].Nreceptors;r++)
	    {
	      fprintf(devmatlab,"plot(poptraces(:,1),poptraces(:,%d),'%c');\nhold on;\n",r+3+j*(Pop[i].Cell[j].Nreceptors+1)+i*(NumberofTraces*(Pop[i].Cell[j].Nreceptors+1)),color[r]);
	    }
	  fprintf(devmatlab,"set(gca,'XLim',[0 %f]);\n",TrialDuration);
	  fprintf(devmatlab,"title('Tot conductances (nS) - Neuron %d')\n",j);
	}
    }

  fclose(devmatlab);

  return 1;
}


int GenMatLabMultiTrial()
{
  FILE *devmatlab;
  int p,trial;

  devmatlab=fopen("mrealmulti.m","w");
  if(devmatlab==NULL) { printf("ERROR: Unable to open multi trial mathlab file\n"); return 0;}

  fprintf(devmatlab,"clear all\nfigure(1);\nclf\n");

  // load all frequency files
  for(trial=0;trial<CurrentTrial+1;trial++)
    {
      fprintf(devmatlab,"load popfreqs%d.dat;\n",trial);
    }

  for(p=0;p<Npop;p++)
    {
      fprintf(devmatlab,"subplot(%d,1,%d);\n",Npop,p+1);
      for(trial=0;trial<CurrentTrial+1;trial++)
	{
	     fprintf(devmatlab,"plot(popfreqs%d(:,1),popfreqs%d(:,%d));\nhold on\n",trial,trial,p+2);
	}
      fprintf(devmatlab,"ylabel('Mean spike freq (Hz)');\nxlabel('Time (ms)');\n");
      fprintf(devmatlab,"set(gca,'XLim',[0 %f]);\n",TrialDuration);
      fprintf(devmatlab,"title('Population: %s');\n",Pop[p].Label);
    }
  fclose(devmatlab);
}

int SaveSpikes(int eventflag, long sno)
{

  static int InitFlag=1;
  static FILE *devspikes[MAXP],*devfreqs;
  static float meanfreqs[MAXP];
  static float timelastevent;
  static float meanfreqsbetweenevents[MAXP];
  static float SpikeBuffer[MAXP][SPIKEBUFFER];
  static int counter;
  static int lasttrial=0;
  static int currentpt, buffersize;  // for SpikeBuffer[][]
  int i,p,TempPointer;
  float  TempBuffer;
  char TempName[100];

  // initialize if it is the first call
  if(InitFlag || (lasttrial!=CurrentTrial)) {

    // if there is a new trial, close all the files of the previous trial
    if(lasttrial!=CurrentTrial) {
      if(FlagSaveAllSpikes)
	{
	  for(p=0;p<Npop;p++)
	    {
	      fclose(devspikes[p]);
	    }
	}
      lasttrial=CurrentTrial;
      // and then open th new ones

    }

    buffersize=(int)(TIMEWINDOWFORFREQ/dt);
    if(buffersize>=SPIKEBUFFER) {printf("Error: SPIKEBUFFER is too small. Edit the code and recomplie the program. \n"); exit(1);};
    currentpt=0;

    // open all files
    printf("\n Time (ms) ");

    for(p=0;p<Npop;p++)
      {
	sprintf(TempName,"pop%d_%d.%ld.dat",p+1,CurrentTrial,sno);
	if(FlagSaveAllSpikes)
	  {
	    devspikes[p]=fopen(TempName,"w");
	    if(devspikes[p]==NULL) return 0;
	  }
	meanfreqs[p]=0.;
	meanfreqsbetweenevents[p]=0.;

	for(i=0;i<SPIKEBUFFER;i++)
	  SpikeBuffer[p][i]=0;

	printf("%12s",Pop[p].Label);
      }
    timelastevent=0.;

    sprintf(TempName,"popfreqs%d.%ld.dat",CurrentTrial,sno);
    devfreqs=fopen(TempName,"w");
    if(devfreqs==NULL) return 0;
    InitFlag=0;
    counter=0;
    printf("\n---------------------------------------------------------------\n");
  }

  if((counter % STEPSFORSAVINGFREQS)==0) fprintf(devfreqs,"%f\t",Time);

  if((counter % STEPSFORPRINTIGFREQS)==0) printf("%7.1f ms",Time);

  for(p=0;p<Npop;p++)
    {
      TempPointer=Pop[p].CTableofSpikes-1;
      if(TempPointer<0) TempPointer=MAXDELAYINDT-1;

      meanfreqsbetweenevents[p]+=(float)Pop[p].NTableofSpikes[TempPointer]/(float)Pop[p].Ncells/dt*1000.;


      meanfreqs[p] += -SpikeBuffer[p][currentpt] + (TempBuffer=(float)Pop[p].NTableofSpikes[TempPointer]/(float)Pop[p].Ncells*1000/TIMEWINDOWFORFREQ);
      SpikeBuffer[p][currentpt] = TempBuffer;

      if(currentpt<buffersize-1)currentpt++; else currentpt=0;
      //      meanfreqs[p]+=dt/TIMEWINDOWFORFREQ*(-meanfreqs[p]+(float)Pop[p].NTableofSpikes[TempPointer]/(float)Pop[p].Ncells/dt*1000.); // compute mean freq in Hz on a time window of 10 ms

      if((counter % STEPSFORPRINTIGFREQS)==0) printf("%12.1f",meanfreqs[p]);
      //      if((counter % STEPSFORPRINTIGFREQS)==0) printf("%s FreqExtSD=%f \n",PopD[0].ReceptorLabel[0],PopD[0].FreqExtSD[0]);
      if((counter % STEPSFORSAVINGFREQS)==0) fprintf(devfreqs,"%f\t",meanfreqs[p]); // mean frequency in Hz per neuron

      if(FlagSaveAllSpikes)
	{
	  for(i=0;i<Pop[p].NTableofSpikes[TempPointer];i++) {
	    fprintf(devspikes[p],"%d\t%f\n",Pop[p].TableofSpikes[TempPointer][i],Time);
	  }
	}
    }

    if((counter % STEPSFORSAVINGFREQS)==0) fprintf(devfreqs,"\n");
    if((counter % STEPSFORPRINTIGFREQS)==0) printf("\n");

  //  if((counter % STEPSFORSAVINGFREQS)==0) fprintf(devfreqs," %f \n",current_freq);
  //  if((counter % STEPSFORPRINTIGFREQS)==0) printf(" %f \n",current_freq);

  if((counter % STEPSFORFLUSHING)==0) {
    if(FlagSaveAllSpikes)
      {
	for(p=0;p<Npop;p++)
	  {
	    fflush(devspikes[p]);
	  }
      }
    fflush(devfreqs);
  }

  if(eventflag) {

    printf("Average:  ");
    for(p=0;p<Npop;p++)
      {
	printf("%12.1f",meanfreqsbetweenevents[p]/(Time-timelastevent)*dt);
	meanfreqsbetweenevents[p]=0.;
      }
    printf("\n");
    fflush(stdout);
    timelastevent=Time;
  }

  counter++;

  return 1;
}

int SaveTraces(long sno)
{
  int i,p,r;
  static int FlagInit=1;
  static FILE *devtraces;
  static int lasttrial=0;
  char TempName[100];
  float Vaux;

  if(NumberofTraces==0) return 1;

  if(FlagInit || (lasttrial!=CurrentTrial)) {
    if(lasttrial!=CurrentTrial) {
      fclose(devtraces);
      lasttrial=CurrentTrial;
    }

    sprintf(TempName,"poptraces%d.%ld.dat",CurrentTrial,sno);
    devtraces=fopen(TempName,"w");
    if(devtraces==NULL) return 0;
    FlagInit=0;
  }

  fprintf(devtraces,"%f\t",Time);
  for(p=0;p<Npop;p++)
    {
      for(i=0;i<NumberofTraces;i++)
	{
	  Vaux=Pop[p].Cell[i].V;
	  if(Vaux>Pop[p].Cell[i].Threshold) Vaux=Pop[p].Cell[i].Threshold;

	  fprintf(devtraces,"%f\t",Pop[p].Cell[i].V);
	  for(r=0;r<Pop[p].Cell[i].Nreceptors;r++)
	    {
	      if(Pop[p].Cell[i].MgFlag[r]) {
		fprintf(devtraces,"%f\t",Pop[p].Cell[i].LS[r]);  // /(1.+exp(-0.062*Vaux)/3.57)); mg block now not included
	      }
	      else {
		fprintf(devtraces,"%f\t",Pop[p].Cell[i].LS[r]);
	      }
	    }
	}
    }
  fprintf(devtraces,"\n");
}


// Handle event
// -------------------------------------------------------------------------------------------------------------------

int HandleEvent(void)
{
  int i,p,r,q,j;
  float efficacy,MeanEff;

  if(Events[CEvent].Type==ENDOFTRIAL) {
    return 0;
  }

  if(Events[CEvent].Type==CHANGEEXTFREQ) {
    p=Events[CEvent].PopNumber;
    r=Events[CEvent].ReceptorNumber;
    PopD[p].FreqExt[r]=Events[CEvent].FreqExt;
    if(Events[CEvent].FreqExtsd!=0)
       PopD[p].FreqExt[r]=Events[CEvent].FreqExt+gasdev()*Events[CEvent].FreqExtsd;
    else
       PopD[p].FreqExt[r]=Events[CEvent].FreqExt;

 //   printf("Events[CEvent].FreqExtDecay=%f \n",Events[CEvent].FreqExtDecay);
    PopD[p].FreqExtTau[r]=Events[CEvent].FreqExtDecay;
    printf("%7.1f ------------------ Event: %s ----------------\n",Time,Events[CEvent].Label);

    efficacy=PopD[p].MeanExtEff[r];
    PopD[p].MeanExtMuS[r]=PopD[p].FreqExt[r]*.001*efficacy*PopD[p].MeanExtCon[r]*PopD[p].Tau[r];
    PopD[p].MeanExtSigmaS[r]=sqrt(PopD[p].Tau[r]*.5*PopD[p].FreqExt[r]*.001*efficacy*efficacy*PopD[p].MeanExtCon[r]);

    for(i=0;i<Pop[p].Ncells;i++)
      {
        do{efficacy=(gasdev()*PopD[p].ExtEffSD[r])+PopD[p].MeanExtEff[r];}while(efficacy<0);
	Pop[p].Cell[i].ExtMuS[r]=PopD[p].FreqExt[r]*.001*efficacy*PopD[p].MeanExtCon[r]*PopD[p].Tau[r];
	Pop[p].Cell[i].ExtSigmaS[r]=sqrt(PopD[p].Tau[r]*.5*PopD[p].FreqExt[r]*.001*efficacy*efficacy*PopD[p].MeanExtCon[r]);
      }
  }

  if(Events[CEvent].Type==CHANGEEXTFREQSD) {
    p=Events[CEvent].PopNumber;
    r=Events[CEvent].ReceptorNumber;
    PopD[p].FreqExtSD[r]=Events[CEvent].FreqExtSD;
    printf("%7.1f ------------------ Event: %s ----------------\n",Time,Events[CEvent].Label);


  }

 if(Events[CEvent].Type==CHANGEMEANEFF) {
    p=Events[CEvent].PopNumber;
    q=Events[CEvent].TargetPopNumber;
    r=Events[CEvent].ReceptorNumber;
    printf("%7.1f ------------------ Event: %s ----------------\n",Time,Events[CEvent].Label);
    MeanEff=Events[CEvent].MeanEff;
    for(i=0;i<PopD[p].NTargetPops;i++)
       if(PopD[p].TargetPops[i]==q){
	 if(PopD[p].SynP[i].TargetReceptor=r)
           PopD[p].SynP[i].MeanEfficacy=MeanEff;
       }
    for(i=0;i<Pop[p].Ncells;i++){
      for(j=0;j<Pop[p].Cell[i].Naxonals;j++)
	if((Pop[p].Cell[i].Axonals[j].TargetPop==q)&&(Pop[p].Cell[i].Axonals[j].TargetReceptor==r))
          Pop[p].Cell[i].Axonals[j].Efficacy=MeanEff;
    }
 }

  CEvent++;
  NextEventTime=Events[CEvent].ETime;

  return 1;
}


main(int argc, char *argv[])
{
  int ti,runflag,tistepforsaving,rseed;
  long iseed;
  long sno=0;

  iseed=1;
  sran1(&iseed);

  //  char prefix[50];

  // read line commands

  flagverbose=0;
  NumberofTrials=1;
  FlagSaveAllSpikes=1;
  //  strcpy(prefix,"cl_");

  if(argc>1) {
    do {
      if(strncmp(argv[argc-1],"-v",2)==0) { flagverbose=1; argc--; continue; }
      if(strncmp(argv[argc-1],"-h",2)==0) {
	printf("realsimu - Ver. 0.8\n Usage:\n-h  : this help\n-v  : verbose mode\n-t# : number of saved traces per population\n");
	printf("-T# : number of trials (the network is the same for each trial, the realization of the ext noise changes)\n");
	printf("-ns : spikes and traces are not saved. Only the mean frequencies are saved for each trial\n");
        printf("-s# : seed for the random number generator.\n");
	return -1;
      }
      if(strncmp(argv[argc-1],"-t",2)==0) {
	NumberofTraces=atoi(&argv[argc-1][2]);
	printf("Number of saved traces: %d\n",NumberofTraces);
	argc--; continue;
      }
      if(strncmp(argv[argc-1],"-s",2)==0) {
	rseed=atoi(&argv[argc-1][2]);
	srand49(rseed);
        iseed=-1*(long)rseed;
        sran1(&iseed);
	printf("Seed for random generator: %d\n",rseed);
	argc--; continue;
      }
      if(strncmp(argv[argc-1],"-ns",3)==0) {
	FlagSaveAllSpikes=0;
	NumberofTraces=0;
	printf("Spikes are not saved\n");
	argc--; continue;
      }
      if(strncmp(argv[argc-1],"-no",3)==0) {
	sno=atol(&argv[argc-1][3]);
	printf("Output file sequence: %ld\n",sno);
	argc--; continue;
      }
      if(strncmp(argv[argc-1],"-T",2)==0) {
	NumberofTrials=atoi(&argv[argc-1][2]);
	printf("Number of trials: %d\n",NumberofTrials);
	argc--; continue;
      }
      /*
      if(strncmp(argv[argc-1],"-p",2)==0) {
	strcpy=(prefix,&argv[argc-1][2]);
	printf("Prefix to all input/output files: %s\n",prefix);
	argc--; continue;
      }
      */
      printf("ERROR: unrecognized option: %s\n",argv[argc]);
      argc--;
    } while(argc>1);
  }

  if(DescribeNetwork()==-1) {printf("Unrecoverable error in parsing the conf file, exiting...\n"); return -1;}
  if(ParseProtocol()==-1) { printf("Unrecoverable error in parsing the protocol file, exiting...\n"); return -1;}

  GenerateNetwork();
  GenMatLabFile();

  for(CurrentTrial=0;CurrentTrial<NumberofTrials;CurrentTrial++)
    {
      report("Trial #\%d\n=====================================================================\n",CurrentTrial);
      CEvent=0;
      NextEventTime=Events[CEvent].ETime;
      runflag=1;

      GenMatLabMultiTrial();

      if(CurrentTrial) InitializeNetwork();

      for(ti=0,Time=0.;runflag;Time+=dt,ti++)
	{
	  SimulateOneTimeStep();
	  // handle all the events
	  if(Time>=NextEventTime) {
	    SaveSpikes(1,sno);
	    while(Time>=NextEventTime)
	      {
		runflag=HandleEvent();
		if(runflag==0) break;
	      }
	  }
	  else SaveSpikes(0,sno);

	  SaveTraces(sno);
	}
      report("\rEnd of the trial\n");
    }
}
