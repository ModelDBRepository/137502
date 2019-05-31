// ZOH solver version of the striatal network simulator
// NOTE: some of the pre-multipliers contain dopamine-dependent terms - these
// will have to go back in the main update loop if time-dependent dopamine is required


#include <mex.h>
#include <matrix.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

// ===============================================================
// RANDOM NUMBER GENERATOR 
#define SHR3 (jz=jsr, jsr^=(jsr<<13), jsr^=(jsr>>17), jsr^=(jsr<<5),jz+jsr)
#define UNI (.5 + (signed) SHR3 * .2328306e-9)
#define RNOR (hz=SHR3, iz=hz&127, (abs(hz)<kn[iz])? hz*wn[iz] : nfix())
#define REXP (jz=SHR3, iz=jz&255, ( jz <ke[iz])? jz*we[iz] : efix())

static unsigned int iz,jz,jsr=123456789,kn[128],ke[256];
static int hz; static float wn[128],fn[128], we[256],fe[256];

float nfix(void) { /*provides RNOR if #define cannot */
    const float r = 3.442620f; static float x, y;
    for(;;){ x=hz*wn[iz];
    if(iz==0){ do{x=-log(UNI)*0.2904764; y=-log(UNI);} while(y+y<x*x);
    return (hz>0)? r+x : -r-x;
    }
    if( fn[iz]+UNI*(fn[iz-1]-fn[iz]) < exp(-.5*x*x) ) return x;
    hz=SHR3; iz=hz&127;if(abs(hz)<kn[iz]) return (hz*wn[iz]);
    }
}
float efix(void) { /*provides REXP if #define cannot */
    float x;
    for(;;){
        if(iz==0) return (7.69711-log(UNI));
        x=jz*we[iz];
        if( fe[iz]+UNI*(fe[iz-1]-fe[iz]) < exp(-x) ) return (x);
        jz=SHR3; iz=(jz&255);
        if(jz<ke[iz]) return (jz*we[iz]);
    }
}

// == This procedure sets the seed and creates the tables ==
void zigset(unsigned int jsrseed) {
    const double m1 = 2147483648.0, m2 = 4294967296.;
    double dn=3.442619855899,tn=dn,vn=9.91256303526217e-3, q;
    double de=7.697117470131487, te=de, ve=3.949659822581572e-3;
    int i; jsr=jsrseed;
    
   /* Tables for RNOR: */
    q=vn/exp(-.5*dn*dn);
    kn[0]=(dn/q)*m1; kn[1]=0;
    wn[0]=q/m1; wn[127]=dn/m1;
    fn[0]=1.; fn[127]=exp(-.5*dn*dn);
    for(i=126;i>=1;i--) {
        dn=sqrt(-2.*log(vn/dn+exp(-.5*dn*dn)));
        kn[i+1]=(dn/tn)*m1; tn=dn;
        fn[i]=exp(-.5*dn*dn); wn[i]=dn/m1;
    }
     /* Tables for REXP */
    q = ve/exp(-de);
    ke[0]=(de/q)*m2; ke[1]=0;
    we[0]=q/m2; we[255]=de/m2;
    fe[0]=1.; fe[255]=exp(-de);
    for(i=254;i>=1;i--) {
        de=-log(ve/de+exp(-de));
        ke[i+1]= (de/te)*m2; te=de;
        fe[i]=exp(-de); we[i]=de/m2;
    }
}

int mxGetScalarInt32(const mxArray* a, int defaultValue = -2147483648)
{
    if (mxIsEmpty(a))
    {
        if (defaultValue == -2147483648) throw "missing input";
        return defaultValue;
    }
    if (!mxIsInt32(a)) throw "not int32";
    if (mxGetNumberOfDimensions(a) != 2) throw "expected scalar";
    if (mxGetM(a) != 1 || mxGetN(a) != 1) throw "expected scalar";
    return mxGetScalar(a);
    //if (fabs(val) > pow(2,30)) throw "value out of range";
    //if (floor(val) != val) throw "value not an integer";
    //return val;
}

double mxGetScalarDouble(const mxArray* a, double defaultValue = -2147483648.1)
{
    if (mxIsEmpty(a))
    {
        if (defaultValue == -2147483648.1) throw "missing input";
        return defaultValue;
    }
    if (!mxIsDouble(a)) throw "not Double";
    if (mxGetNumberOfDimensions(a) != 2) throw "expected scalar";
    if (mxGetM(a) != 1 || mxGetN(a) != 1) throw "expected scalar";
    return mxGetScalar(a);
}

struct MatrixInt32
{
    int* data;
    int M, N;
};

MatrixInt32 mxGetMatrixInt32(const mxArray* a, int M = -1, int N = -1)
{
    MatrixInt32 ret;
    if (!mxIsInt32(a)) throw "not int32";
    if (mxGetNumberOfDimensions(a) != 2) throw "expected 2D matrix";
    if (mxIsComplex(a)) throw "expected real matrix";
    if (mxIsEmpty(a)) ret.data = NULL;
    else ret.data = (int*)mxGetData(a);
    ret.M = mxGetM(a);
    ret.N = mxGetN(a);
    if (M != -1 && M != ret.M) throw "interger matrix has wrong dimension (M)";
    if (N != -1 && N != ret.N) throw "interger matrix has wrong dimension (N)";
    return ret;
}

struct MatrixDouble
{
    double* data;
    int M, N;
};

MatrixDouble mxGetMatrixDouble(const mxArray* a, int M = -1, int N = -1)
{
    MatrixDouble ret;
    if (!mxIsDouble(a)) throw "not double";
    if (mxGetNumberOfDimensions(a) != 2) throw "expected 2D matrix";
    if (mxIsComplex(a)) throw "expected real matrix";
    if (mxIsEmpty(a)) ret.data = NULL;
    else ret.data = (double*)mxGetData(a);
    ret.M = mxGetM(a);
    ret.N = mxGetN(a);
    if (M != -1 && M != ret.M) throw "double matrix has wrong dimension (M)";
    if (N != -1 && N != ret.N) throw "double matrix has wrong dimension (N)";
    return ret;
}

// ===============================================================
// Main simulation function
void execute(int nlhs, mxArray *plhs[],
int nrhs, const mxArray *prhs[])
{   
    // ===============================================================
    //  The Inputs
    // ---------------------------------------------------------------
    if (nrhs != 72) throw "not enough inputs";
    
    // ---------------------------------------------------------------
    // Get the inputs
    double tstart = mxGetScalarDouble(prhs[0]);
    double tfinal = mxGetScalarDouble(prhs[1]);
    double dt = mxGetScalarDouble(prhs[2]);
//    double dt_fs = mxGetScalarDouble(prhs[3]);       
    int tend = (int)(tfinal / dt);
    
    MatrixDouble MSparams = mxGetMatrixDouble(prhs[3]);
    MatrixDouble FSparams = mxGetMatrixDouble(prhs[4]);

    double Eglu = mxGetScalarDouble(prhs[5]);
    double Egaba = mxGetScalarDouble(prhs[6]);    
    double ts_glu_AMPA = mxGetScalarDouble(prhs[7]);
    double ts_glu_NMDA = mxGetScalarDouble(prhs[8]);
    double ts_gaba = mxGetScalarDouble(prhs[9]);
    double tau_fsgap = mxGetScalarDouble(prhs[10]);
    
    int MSspikebuffer = mxGetScalarInt32(prhs[11]);
    int FSspikebuffer = mxGetScalarInt32(prhs[12]);
    
    MatrixDouble initVms = mxGetMatrixDouble(prhs[13]);
    MatrixDouble initUms = mxGetMatrixDouble(prhs[14], initVms.M, initVms.N);
    MatrixDouble initVfs = mxGetMatrixDouble(prhs[15]);
    MatrixDouble initUfs = mxGetMatrixDouble(prhs[16], initVfs.M, initVfs.N);
    MatrixDouble initVfsgap = mxGetMatrixDouble(prhs[17]);
    
    // get the dimensions of the MS and FS networks
    const int *dims_ms = mxGetDimensions(prhs[13]);
    int ndim_ms = mxGetNumberOfDimensions(prhs[13]);
    int N_MS = initVms.M;
    
    const int *dims_fs = mxGetDimensions(prhs[15]);
    int ndim_fs = mxGetNumberOfDimensions(prhs[15]);
    int N_FS = initVfs.M;
    
    const int *dims_fsgap = mxGetDimensions(prhs[17]);
    int ndim_fsgap = mxGetNumberOfDimensions(prhs[17]);
    int N_fsgap = initVfsgap.M*initVfsgap.N;
    
    MatrixDouble initSEQ_MSglu = mxGetMatrixDouble(prhs[18], N_MS);
    MatrixDouble initSEQ_FSglu = mxGetMatrixDouble(prhs[19], N_FS);
    MatrixDouble initSEQ_MSGABA = mxGetMatrixDouble(prhs[20], N_MS);
    MatrixDouble initSEQ_FSGABA = mxGetMatrixDouble(prhs[21], N_FS);
    MatrixDouble initCTX = mxGetMatrixDouble(prhs[22]);

    // *** NEED TO REMOVE THIS AND FIND ANOTHER WAY TO SET UP MATRICES FURTHER DOWN ***
    const int *initSEQ_MSGABA_dims = mxGetDimensions(prhs[20]);
    int initSEQ_MSGABA_ndim = mxGetNumberOfDimensions(prhs[20]);
    
    const int *initSEQ_FSGABA_dims = mxGetDimensions(prhs[21]);
    int initSEQ_FSGABA_ndim = mxGetNumberOfDimensions(prhs[21]);
    
    // *** Not used currently, but will need bounds check put on dimensions ***
    MatrixDouble Ims = mxGetMatrixDouble(prhs[23]);
    MatrixDouble Ifs = mxGetMatrixDouble(prhs[24]);
    
    MatrixInt32 Cctms = mxGetMatrixInt32(prhs[25], N_MS, 1); // need to set bound on to N_MS, no + 1
    MatrixInt32 Cctms_b = mxGetMatrixInt32(prhs[26], N_MS+1, 1);
    MatrixInt32 Cctms_d = mxGetMatrixInt32(prhs[27], N_MS, 1);
    MatrixDouble Cctms_w = mxGetMatrixDouble(prhs[28], N_MS, 1);
    
//    double a_ms = mxGetScalarDouble(prhs[29]);
    
    MatrixInt32 Cmsms = mxGetMatrixInt32(prhs[29]);
    MatrixInt32 Cmsms_b = mxGetMatrixInt32(prhs[30], N_MS+1, 1);
    MatrixInt32 Cmsms_d = mxGetMatrixInt32(prhs[31], Cmsms.M, Cmsms.N);
    MatrixDouble Cmsms_w = mxGetMatrixDouble(prhs[32], Cmsms.M, Cmsms.N);

    MatrixInt32 Cfsms = mxGetMatrixInt32(prhs[33]);
    MatrixInt32 Cfsms_b = mxGetMatrixInt32(prhs[34], N_FS+1, 1);
    MatrixInt32 Cfsms_d = mxGetMatrixInt32(prhs[35], Cfsms.M, Cfsms.N);
    MatrixDouble Cfsms_w = mxGetMatrixDouble(prhs[36], Cfsms.M, Cfsms.N);
    const int *dims_Cfsms_b = mxGetDimensions(prhs[34]);
    
    MatrixInt32 Cctfs = mxGetMatrixInt32(prhs[37], N_FS, 1);
    MatrixInt32 Cctfs_b = mxGetMatrixInt32(prhs[38], N_FS+1, 1);
    MatrixInt32 Cctfs_d = mxGetMatrixInt32(prhs[39], N_FS, 1);
    MatrixDouble Cctfs_w = mxGetMatrixDouble(prhs[40], N_FS, 1);
    const int *dims_Cctfs_b = mxGetDimensions(prhs[38]);
//    double a_fs = mxGetScalarDouble(prhs[42]);

    MatrixInt32 Cfsfs = mxGetMatrixInt32(prhs[41]);
    MatrixInt32 Cfsfs_b = mxGetMatrixInt32(prhs[42], N_FS+1, 1);
    MatrixInt32 Cfsfs_d = mxGetMatrixInt32(prhs[43], Cfsfs.M, Cfsfs.N);
    MatrixDouble Cfsfs_w = mxGetMatrixDouble(prhs[44], Cfsfs.M, Cfsfs.N);
    const int *dims_Cfsfs_b = mxGetDimensions(prhs[42]);
    
    MatrixInt32 Cgapfs = mxGetMatrixInt32(prhs[45]);
    MatrixInt32 Cgapfs_b = mxGetMatrixInt32(prhs[46], N_FS+1, 1);
    MatrixDouble Cgapfs_w = mxGetMatrixDouble(prhs[47], Cgapfs.M, Cgapfs.N);
    MatrixInt32 Pgapfs = mxGetMatrixInt32(prhs[48]);
    
    MatrixDouble CTX_state = mxGetMatrixDouble(prhs[49]);
    MatrixInt32 CHAN1_MS = mxGetMatrixInt32(prhs[50]);
    MatrixInt32 CHAN1_FS = mxGetMatrixInt32(prhs[51]);
    MatrixInt32 CHAN2_MS = mxGetMatrixInt32(prhs[52]);
    MatrixInt32 CHAN2_FS = mxGetMatrixInt32(prhs[53]);
    const int *dims_CHAN1_MS = mxGetDimensions(prhs[50]);
    int N_CHAN1_MS = dims_CHAN1_MS[0];
    const int *dims_CHAN1_FS = mxGetDimensions(prhs[51]);
    int N_CHAN1_FS = dims_CHAN1_FS[0];
    
    MatrixInt32 N_MSSEG = mxGetMatrixInt32(prhs[54], N_MS, 1);
    MatrixDouble r_MSSEG = mxGetMatrixDouble(prhs[55], N_MS, 1);
    MatrixDouble alpha_MSSEG = mxGetMatrixDouble(prhs[56], N_MS, 1);
    
    MatrixInt32 N_FSSEG = mxGetMatrixInt32(prhs[57], N_FS, 1);
    MatrixDouble r_FSSEG = mxGetMatrixDouble(prhs[58], N_FS, 1);
    MatrixDouble alpha_FSSEG = mxGetMatrixDouble(prhs[59], N_FS, 1);

    double glu_ratio = mxGetScalarDouble(prhs[60]); 
    double DA = mxGetScalarDouble(prhs[61]);
    
    MatrixInt32 RecordChan_MS = mxGetMatrixInt32(prhs[62]);
    const int dims_RecordChan_MS_out [2] = {tend, RecordChan_MS.M};
    
    MatrixDouble PULSE = mxGetMatrixDouble(prhs[63]);
    double Nctx_ms = mxGetScalarDouble(prhs[64]); 
    double Nctx_fs = mxGetScalarDouble(prhs[65]); 
    double ts_spks = mxGetScalarDouble(prhs[66]); 
    
    MatrixInt32 Pt = mxGetMatrixInt32(prhs[67]);
    MatrixInt32 Pch = mxGetMatrixInt32(prhs[68], Pt.M, Pt.N);
    MatrixDouble Phz = mxGetMatrixDouble(prhs[69], Pt.M, Pt.N);

    int RANDSEED = mxGetScalarInt32(prhs[70]);
    
    const char* filename = mxArrayToString(prhs[71]);
    
    // ===============================================================
    // set the random numer seeds
    zigset(RANDSEED);
        
    // ===============================================================
    // Set up easy access to the MS and FS parameters
    #define _MS(__m,__n) (MSparams.data[__n * MSparams.M + __m])    
    #define _FS(__m,__n) (FSparams.data[__n * FSparams.M + __m])
    
    // ===============================================================
    // get some dimensions     
    if (!filename) throw "no filename supplied for log file";
    //  open log file
    FILE* fid = fopen(filename, "w");
    if (!fid) throw "failed open log file";
    
    fprintf(fid, "*************************************************** \n");
    fprintf(fid, "*************** STARTING SIMULATION *************** \n");
    fprintf(fid, "*************************************************** \n");
    fprintf(fid, "tstart = %f msec \n", tstart);
    fprintf(fid, "tfinal = %f msec \n", tfinal);
    fprintf(fid, "dt = %f msec \n", dt);
    fprintf(fid, "tend = %i iterations \n", tend);
    fprintf(fid, " \n");    
    fprintf(fid, "N_MS = %i \n", N_MS);    
    fprintf(fid, "Example MS parameters are C = %f, vr = %f, vt = %f, a = %f, b = %f, c = %f, d = %f, vp = %f, k = %f, EDA = %f, alpha = %f, beta1 = %f, beta2 = %f, gDA = %f \n", _MS(0,0), _MS(0,1), _MS(0,2), _MS(0,3), _MS(0,4), _MS(0,5), _MS(0,6), _MS(0,7), _MS(0,8), _MS(0,9), _MS(0,10), _MS(0,11), _MS(0,12), _MS(0,13), _MS(0,14));    
    fprintf(fid, " \n");    
    fprintf(fid, "N_FS = %i \n", N_FS);    
    fprintf(fid, "Example FS parameters are C = %f, vr = %f, vt = %f, k = %f, a = %f, b = %f, c = %f, d = %f, vpeak = %f, vb = %f, eta = %f, epsilon = %f \n", _FS(0,0), _MS(0,1), _FS(0,2), _FS(0,3), _FS(0,4), _FS(0,5), _FS(0,6), _FS(0,7), _FS(0,8), _FS(0,9), _FS(0,10), _FS(0,11));    
    fprintf(fid, " \n");    
    fprintf(fid, "Eglu = %f \n", Eglu);    
    fprintf(fid, "Egaba = %f \n", Egaba);    
    fprintf(fid, "ts_glu_AMPA = %f \n", ts_glu_AMPA);    
    fprintf(fid, "ts_glu_NMDA = %f \n", ts_glu_NMDA);    
    fprintf(fid, "ts_gaba = %f \n", ts_gaba);    
    fprintf(fid, "tau_fsgap = %f \n", tau_fsgap);    
    fprintf(fid, " \n"); 
    fprintf(fid, "gAMPA_MS = %f \n", Cctms_w.data[0]);    
    fprintf(fid, "gNMDA_MS = %f \n", glu_ratio * Cctms_w.data[0]);    
    fprintf(fid, "gGLUTAMATE_FS = %f \n", Cctfs_w.data[0]);    
    fprintf(fid, " \n"); 
    fprintf(fid, "MSspikebuffer = %i \n", MSspikebuffer);    
    fprintf(fid, "FSspikebuffer = %i \n", FSspikebuffer);    
    fprintf(fid, " \n"); 
    fprintf(fid, "*************************************************** \n");
    fclose(fid);
       
    // ===============================================================
    // Set and init some parameters 
    double T;                               // Simulation Time
//    int extasteps = (int)(dt_ms / dt_fs);   // Number of extra iterations needed to solve the FS neurons
    float VMS,UMS,VFS,UFS,U;                // temp variables for the ms and fs neurons
    
    // parameters for the MS neurons 
    int MSspikecount = 0;
    int SEQind_ms = 0;
    int SpikeEventCycle_MSGABA = 0;
    int maxSEQdelay_MSGABA = initSEQ_MSGABA.N - 1;
    int SpikeQueSize_MS = N_MS * (maxSEQdelay_MSGABA+1);
    double Mg_ms[N_MS];     // Mg block for the MS NMDA channels
    double lambda1_ms;
    double lambda2_ms;
    
    // parameters for the FS neurons 
    int FSspikecount = 0;
    int SEQind_fs;
    int SpikeEventCycle_FSGABA = 0;    
    int maxSEQdelay_FSGABA = initSEQ_FSGABA_dims[1] - 1;
    int SpikeQueSize_FS = N_FS * (maxSEQdelay_FSGABA+1);
    double lambda1_fs;
    double lambda2_fs;
    
    // parameters for the MS neuron spike-event generator
    float p_MSSEG; // adjust for time in msec, not seconds!
    int spk_ms[N_MS];
    float U_MSSEG;
    int R_MSSEG;
    double gluExp_ms_AMPA = exp(-dt / ts_glu_AMPA);
    double gluExp_ms_NMDA = exp(-dt / ts_glu_NMDA);
    double gabaExp_ms = exp(-dt / ts_gaba);
    float EMg = 1; // magnesium ion concentration in mM
    double SpksExp = exp(-dt / ts_spks);
    
    // parameters for the FS neuron spike-event generator
    float p_FSSEG; // adjust for time in msec, not seconds!
    int spk_fs[N_FS];
    float U_FSSEG;
    int R_FSSEG;
    double gluExp_fs = exp(-dt / ts_glu_AMPA);
    double gabaExp_fs = exp(-dt / ts_gaba);
    
    // ZOH temporary variables
    float qsum, A, roots, exp_tc, C0, z;
    float B, Y, X, tmp, cx, sx;
    
    // ===============================================================
    // The Outputs
    double *tout;       // Time stamps
    double *Vmsout;     // example MS membrane potential
    double *Vfsout;     // Example FS membran potential
    double *STms;       // MS neuron spike times
    double *STfs;       // FS neuron spike times
    
    // output states so we can continue simulations
    double *Vms;        // MS membrane potentials
    double *Ums;        // MS U values
    double *Vfs;        // FS membrane potentials
    double *Ufs;        // FS U values
    double *Vfsgap;     // Membrane potentials of the FS gap junctions
    
    double *Gglu_ms_AMPA;     // *** Place holder for current glutamate input conductance to the MS neurons ***
    double *Gglu_ms_NMDA;     // *** Place holder for current glutamate input conductance to the MS neurons ***
    double *Ggaba_ms;    // Current input GABA conductance to the MS neurons
    double *SEQ_MSGABA;  // Incomming GABAergic spike buffer for the FS neurons
    
    double *Gglu_fs;     // *** Place holder for current glutamate input conductance to the MS neurons ***
    double *Ggaba_fs;    // Current input GABA conductance to the MS neurons
    double *SEQ_FSGABA;  // Incomming GABAergic spike buffer for the FS neurons
    
    double *RecordChan_MS_out; // output from the recording electrode
       
    //create the output arrays
    plhs[0] = mxCreateDoubleMatrix(tend, 1, mxREAL);
    tout = mxGetPr(plhs[0]);
    plhs[1] = mxCreateDoubleMatrix(tend, 1, mxREAL);
    Vmsout = mxGetPr(plhs[1]);
    plhs[2] = mxCreateDoubleMatrix(tend, 1, mxREAL);
    Vfsout = mxGetPr(plhs[2]);
    plhs[3] = mxCreateDoubleMatrix(MSspikebuffer, 2, mxREAL);
    STms = mxGetPr(plhs[3]);
    plhs[4] = mxCreateDoubleMatrix(FSspikebuffer, 2, mxREAL);
    STfs = mxGetPr(plhs[4]);
    plhs[5] = mxCreateNumericArray(ndim_ms,dims_ms,mxDOUBLE_CLASS,mxREAL);
    Vms = mxGetPr(plhs[5]);
    plhs[6] = mxCreateNumericArray(ndim_ms,dims_ms,mxDOUBLE_CLASS,mxREAL);
    Ums = mxGetPr(plhs[6]);    
    plhs[7] = mxCreateNumericArray(ndim_fs,dims_fs,mxDOUBLE_CLASS,mxREAL);
    Vfs = mxGetPr(plhs[7]);
    plhs[8] = mxCreateNumericArray(ndim_fs,dims_fs,mxDOUBLE_CLASS,mxREAL);
    Ufs = mxGetPr(plhs[8]);
    plhs[9] = mxCreateNumericArray(ndim_fsgap,dims_fsgap,mxDOUBLE_CLASS,mxREAL);
    Vfsgap = mxGetPr(plhs[9]);
    
    plhs[10] =  mxCreateNumericArray(ndim_ms,dims_ms,mxDOUBLE_CLASS,mxREAL);
    Gglu_ms_AMPA = mxGetPr(plhs[10]);
    plhs[11] =  mxCreateNumericArray(ndim_ms,dims_ms,mxDOUBLE_CLASS,mxREAL);
    Gglu_ms_NMDA = mxGetPr(plhs[11]);
    
    plhs[12] =  mxCreateNumericArray(ndim_ms,dims_ms,mxDOUBLE_CLASS,mxREAL);
    Ggaba_ms = mxGetPr(plhs[12]);
    plhs[13] =  mxCreateNumericArray(initSEQ_MSGABA_ndim,initSEQ_MSGABA_dims,mxDOUBLE_CLASS,mxREAL);
    SEQ_MSGABA = mxGetPr(plhs[13]);    
    
    plhs[14] =  mxCreateNumericArray(ndim_fs,dims_fs,mxDOUBLE_CLASS,mxREAL);
    Gglu_fs = mxGetPr(plhs[14]);
    plhs[15] =  mxCreateNumericArray(ndim_fs,dims_fs,mxDOUBLE_CLASS,mxREAL);
    Ggaba_fs = mxGetPr(plhs[15]);
    plhs[16] =  mxCreateNumericArray(initSEQ_FSGABA_ndim,initSEQ_FSGABA_dims,mxDOUBLE_CLASS,mxREAL);
    SEQ_FSGABA = mxGetPr(plhs[16]);
    
    plhs[17] = mxCreateNumericArray(2, dims_RecordChan_MS_out,mxDOUBLE_CLASS,mxREAL);
    RecordChan_MS_out = mxGetPr(plhs[17]);
    
    // ===============================================================
    // Set up easy access to the cortical state and recording channels
    const int *CTX_state_dims = mxGetDimensions(prhs[58]);  
    #define _CTX_state(__m,__n) (CTX_state[__n*CTX_state_dims[0]+__m]) 
    
    const int *RecordChan_MS_out_dims = mxGetDimensions(plhs[17]);      
    #define _RecordChan_MS_out(__m,__n) (RecordChan_MS_out[__n*RecordChan_MS_out_dims[0]+__m]) 
    
    // ===============================================================
    // init the variables
    for (int i = 0; i < N_MS; i++){
        Vms[i] = initVms.data[i];
        Ums[i] = initUms.data[i];
    }
    for (int i = 0; i < N_FS; i++){
        Vfs[i] = initVfs.data[i];
        Ufs[i] = initUfs.data[i];
    }
    
    for (int i = 0; i < N_fsgap; i++){
        Vfsgap[i] = initVfsgap.data[i];
    }
        
    float IAMPA_ms[N_MS];
    float INMDA_ms[N_MS];
    float IDA_ms[N_MS];
    float IGABA_ms[N_MS];
    float Isyn_ms[N_MS];
    float tmpIDA_ms;
    float kD2[N_MS];
    // ZOH intermediate variables
    float xa_ms[N_MS]; float xb_ms[N_MS]; float x4_ms[N_MS]; float x5_ms[N_MS]; float ux1_ms[N_MS]; float vp_ms[N_MS]; float u_exp_tc_ms[N_MS]; 
    for (int i = 0; i < N_MS; i++){
        IAMPA_ms[i] = 0;
        INMDA_ms[i] = 0;
        IDA_ms[i] = 0;
        IGABA_ms[i] = 0;
        Isyn_ms[i] = 0;
        tmpIDA_ms = 0;
        kD2[i] =  _MS(i,8) * (1 - _MS(i,10)*lambda2_ms);
        // ZOH intermediates
        xa_ms[i] = 2*kD2[i];  // 2*k
        xb_ms[i] = -kD2[i] * (_MS(i,1)+_MS(i,2)); // k*(vr+vt)
        x4_ms[i] = 4 * kD2[i];  // 4 * k
        x5_ms[i] = xb_ms[i] * xb_ms[i] - x4_ms[i] * kD2[i] * _MS(i,1) * _MS(i,2);
        ux1_ms[i] = _MS(i,4) * _MS(i,1); // b * vr
        vp_ms[i] = (_MS(i,1)+_MS(i,2)) / 2;  // (vr + vt)/2 
        u_exp_tc_ms[i] = exp(-_MS(i,3)*dt); 
    }
   
    float IAMPA_fs[N_FS];
    float IGABA_fs[N_FS];
    float IDA_fs[N_FS];
    float Isyn_fs[N_FS];
    float Igapfs[N_FS];
    float vr_fs[N_FS];
    float xa_fs[N_FS]; float xb_fs[N_FS]; float x4_fs[N_FS]; float x5_fs[N_FS]; float ux1_fs[N_FS]; float vp_fs[N_FS]; float u_exp_tc_fs[N_FS]; 
   
    for (int j = 0; j < N_FS; j++){
        IAMPA_fs[j] = 0.0;
        IGABA_fs[j] = 0.0;
        IDA_fs[j] = 0.0;
        Isyn_fs[j] = 0.0;
        Igapfs[j] = 0.0;
        vr_fs[j] = _FS(j, 1) * (1-_FS(j, 10)*lambda1_fs);
        // ZOH intermediates
        xa_fs[j] = 2*_FS(j,3);  // 2*k
        xb_fs[j] = -_FS(j,3) * (vr_fs[j]+_FS(j,2)); // k*(vr+vt)
        x4_fs[j] = 4 * _FS(j,3);  // 4 * k
        x5_fs[j] = xb_fs[j] * xb_fs[j] - x4_fs[j] * _FS(j,3) * vr_fs[j] * _FS(j,2); 
        ux1_fs[j] = _FS(j,5) * vr_fs[j]; // b * vr
        vp_fs[j] = (vr_fs[j]+_FS(j,2)) / 2;  // (vr + vt)/2 
        u_exp_tc_fs[j] = exp(-_FS(j,4)*dt);        
        }
    
    int Vgap1,Vgap2;
    
    // ===============================================================
    // Active D1 and D2 receptors for the MS and FS neurons
    lambda1_ms = DA; lambda2_ms = DA;
    lambda1_fs = DA; lambda2_fs = DA;    
    
    // ===============================================================
    // run the simulation
    int SelectionCounter = 0;
    for (int t = (int)tstart; t < tend; t++){
        T = t * dt;
           
        // reset the random seed
        if (UNI < (0.01*dt)){
            RANDSEED++;
            zigset(RANDSEED);
        }
        
        // ================================================================
        // Look for changes in the selection pulses
        if (Pt.data[SelectionCounter] == t){
            printf("Selection Pulse Detected %f \n", T);
            
            if (Pch.data[SelectionCounter] == 1){
                printf("Selection Pulse Channel %i \n", Pch.data[SelectionCounter]);
                for (int i = 0; i < CHAN1_MS.M; i++){
                    r_MSSEG.data[CHAN1_MS.data[i]] = Phz.data[SelectionCounter];
                }
                for (int i = 0; i < CHAN1_FS.M; i++){
                    r_FSSEG.data[CHAN1_FS.data[i]] = Phz.data[SelectionCounter];
                }
            }
            if (Pch.data[SelectionCounter] == 2){
                printf("Selection Pulse Channel %i \n", Pch.data[SelectionCounter]);
                for (int i = 0; i < CHAN2_MS.M; i++){
                    r_MSSEG.data[CHAN2_MS.data[i]] = Phz.data[SelectionCounter];
                }
                for (int i = 0; i < CHAN2_FS.M; i++){
                    r_FSSEG.data[CHAN2_FS.data[i]] = Phz.data[SelectionCounter];
                }
            } 
            SelectionCounter++;
        }
        
        // ===============================================================
        // Update the MS neurons                                
        for (int i = 0; i < N_MS; i++){            
            // ---------------------------------------------------------------
            // Update the MS spike-event que
            spk_ms[i] = 0;                          // reset the spike que
            int T_MSSEG = 0;
            p_MSSEG = r_MSSEG.data[i] * dt * 0.001;
            U_MSSEG = UNI;
            R_MSSEG = (U_MSSEG <= p_MSSEG);         // the reference spike train
            for (int ist = 0; ist < N_MSSEG.data[i]; ist++){
                if (UNI <= alpha_MSSEG.data[i]){
                    T_MSSEG = R_MSSEG;
                }
                else {
                    T_MSSEG = (UNI <= p_MSSEG);
                }
                spk_ms[i] = spk_ms[i] + T_MSSEG;
            }
            
            // ---------------------------------------------------------------
            //Update the MS neurons cortical input
             Gglu_ms_AMPA[i] = Gglu_ms_AMPA[i] + ((Cctms_w.data[i]) * spk_ms[i] / ts_glu_AMPA);
             Gglu_ms_AMPA[i] = Gglu_ms_AMPA[i] * gluExp_ms_AMPA;
             Gglu_ms_NMDA[i] = Gglu_ms_NMDA[i] + (glu_ratio * Cctms_w.data[i] * spk_ms[i] / ts_glu_NMDA);
             Gglu_ms_NMDA[i] = Gglu_ms_NMDA[i] * gluExp_ms_NMDA;
             
            // ---------------------------------------------------------------
            // Mg block of the MS NMDA channels
            Mg_ms[i] = 1 / ( 1 + (EMg / 3.57) * exp(-Vms[i] * 0.062) ); 
            
            // ---------------------------------------------------------------
            // get PSPs from the other MS and FS neurons
            Ggaba_ms[i] += SEQ_MSGABA[SpikeEventCycle_MSGABA*N_MS + i];	// Adds current GABAergic PSPs to the conductance bin
            Ggaba_ms[i] = Ggaba_ms[i] * gabaExp_ms;
            SEQ_MSGABA[SpikeEventCycle_MSGABA*N_MS + i] = 0;              // reset the PSP buffer
            
            // ---------------------------------------------------------------
            // update the membrane
            VMS = Vms[i];   // save the previous state
            UMS = Ums[i];   // save the previous state
                        
            IGABA_ms[i] = (Ggaba_ms[i] * (Egaba - VMS));                       
            IAMPA_ms[i] = (Gglu_ms_AMPA[i] * (Eglu - VMS)) * (1 + _MS(i,12) * lambda1_ms);
            INMDA_ms[i] = Mg_ms[i] * (Gglu_ms_NMDA[i] * (Eglu - VMS)) * (1 + _MS(i,11) * lambda2_ms);
            IDA_ms[i] = lambda1_ms * _MS(i,13) * (Vms[i] - _MS(i,9));
            Isyn_ms[i] = IAMPA_ms[i] + INMDA_ms[i] + IGABA_ms[i] + IDA_ms[i];
            
            // compute v
            qsum = x5_ms[i] - x4_ms[i] *(Isyn_ms[i]-UMS);
            if (qsum >=0){
                // then root is real
                A = sqrt(qsum); 
                roots = (A-xb_ms[i]) / xa_ms[i];
                exp_tc = exp(-dt*A/_MS(i,0));
                C0 = exp_tc / (VMS-roots);
                z = -kD2[i] / A * (1-exp_tc) + C0;
                Vms[i] = roots + 1/z;
            }
            else {
                // root is complex
                B = sqrt(fabs(qsum));
                Y = VMS - vp_ms[i];
                X = -dt*B/_MS(i,0);
                tmp = kD2[i]*Y/B;
                cx = cos(X); sx = sin(X);
                Vms[i] = vp_ms[i] + (Y*cx + Y*tmp*sx - (B/(4*kD2[i]))*sx) / 
                            (0.5*(1+cx) + 2*tmp*sx - 2*tmp*tmp*(cx-1));
            }
            
            // Vms[i] = VMS + dt * ( (  _MS(i,8) * (1 - _MS(i,10)*lambda2_ms) * (VMS - _MS(i,1)) * (VMS - _MS(i,2)) - Ums[i]) + Isyn_ms[i] ) / _MS(i,0);                        
            
            // compute u
            //Ums[i] = Ums[i] + dt * _MS(i,3) * (_MS(i,4) * (VMS - _MS(i,1)) - Ums[i]);
            Ums[i] = (_MS(i,4)*VMS-ux1_ms[i]) * (1-u_exp_tc_ms[i]) + UMS * u_exp_tc_ms[i];            
            
            if (Vms[i] >= _MS(i,7)){ // Check for spike events
                VMS = _MS(i,7);
                Vms[i] = _MS(i,5);
                Ums[i] = Ums[i] + _MS(i,6);
                
                // ---------------------------------------------------------------
                // save the cell number and spike time
                if (MSspikecount < MSspikebuffer){
                    STms[MSspikecount] = i;
                    STms[MSspikecount+MSspikebuffer] = T;
                    MSspikecount++;
                }
                else if (MSspikecount >= MSspikebuffer){
                    printf("Warning, exceeded MS spike buffer size \n");
                    printf("%i %i \n", MSspikecount, MSspikebuffer);
                    t = tend;
                }
                
                // ---------------------------------------------------------------
                // Update the Spike event que for the target cells
                for (int targcellind = Cmsms_b.data[i]; targcellind <= Cmsms_b.data[i+1]-1; targcellind++){ // for each of the target cells
                    // get the index in the spike que for the target cell, given its delay
                    SEQind_ms = Cmsms.data[targcellind] + ((SpikeEventCycle_MSGABA-1) + Cmsms_d.data[targcellind])*N_MS; 
                    
                    if (SEQind_ms >= SpikeQueSize_MS){ // if SEind is larger than the spike que, wrap around
                        SEQind_ms = SEQind_ms - SpikeQueSize_MS;
                    }
                                        
                    // add the spike event to the spike que for the target cell
                    SEQ_MSGABA[SEQind_ms] = SEQ_MSGABA[SEQind_ms] + (Cmsms_w.data[targcellind] / ts_gaba); // update the Spike event que
                }
            }
        }
        
        // update the FS gap junctions
        for (int k = 0; k < N_fsgap; k++){
            Vgap1 = Pgapfs.data[k];
            Vgap2 = Pgapfs.data[k + N_fsgap];
            
            Vfsgap[k] = Vfsgap[k] + dt * ((Vfs[Vgap1] - Vfsgap[k]) + (Vfs[Vgap2] - Vfsgap[k])) / tau_fsgap;
        }
        
        // ---------------------------------------------------------------
        // for each FS neuron
        // ---------------------------------------------------------------
        for (int j = 0; j < N_FS; j++){
            float testIgapfs = Igapfs[j];
            
            spk_fs[j] = 0;                          // reset the spike que
            int T_FSSEG = 0;
            p_FSSEG = r_FSSEG.data[j] * dt * 0.001;
            U_FSSEG = UNI;
            R_FSSEG = (U_FSSEG <= p_FSSEG);         // the reference spike train
            for (int ist = 0; ist < N_FSSEG.data[j]; ist++){
                if (UNI <= alpha_FSSEG.data[j]){
                    T_FSSEG = R_FSSEG;
                }
                else {
                    T_FSSEG = (UNI <= p_FSSEG);
                }
                spk_fs[j] = spk_fs[j] + T_FSSEG;
            }

            // ---------------------------------------------------------------
            //Update the FS neurons cortical input
            Gglu_fs[j] = Gglu_fs[j] + (Cctfs_w.data[j] * spk_fs[j] / ts_glu_AMPA);
            Gglu_fs[j] = Gglu_fs[j] * gluExp_fs;
            
            // ---------------------------------------------------------------
            // get GABAergic input from the other FS neurons
            Ggaba_fs[j] += SEQ_FSGABA[SpikeEventCycle_FSGABA*N_FS + j];	// Adds current GABAergic PSPs to the conductance bin
            Ggaba_fs[j] = Ggaba_fs[j] * gabaExp_fs;
            if (SEQ_FSGABA[SpikeEventCycle_FSGABA*N_FS + j] > 1000){
                printf("Too many spikes! %i %i \n", SEQ_FSGABA[SpikeEventCycle_FSGABA*N_FS + j], t);
                throw "Too many spikes";
            }
            
            SEQ_FSGABA[SpikeEventCycle_FSGABA*N_FS + j] = 0;              // reset the PSP buffer
            
            // ---------------------------------------------------------------
            // calculate the current from the gap junctions
            Igapfs[j] = 0.0;
            
            for (int source = Cgapfs_b.data[j]; source < Cgapfs_b.data[j+1]; source++){
                Igapfs[j] = Igapfs[j] + Cgapfs_w.data[source] * (Vfsgap[Cgapfs.data[source]] - Vfs[j]);
                if (isinf(Igapfs[j])){
                    printf("\n \n **** Warning, Igapfs is a inf. Vfsgap %f Cgapfs %i source %i **** \n \n", Vfsgap[Cgapfs.data[source]], Cgapfs.data[source], source);
                    throw "Igapfs is a Inf";
                }
            }
            
            // ---------------------------------------------------------------
            //Update the FS neurons
            VFS = Vfs[j];   // save the previous state
            UFS = Ufs[j];   // save the previous state
            
            // DA modulation of the GABA currents
            IAMPA_fs[j] = (Gglu_fs[j] * (Eglu - Vfs[j]));
            IGABA_fs[j] = (Ggaba_fs[j] * (Egaba - Vfs[j])) * (1 - _FS(j, 11)*lambda2_fs);
            Isyn_fs[j] = IAMPA_fs[j] + IGABA_fs[j];
            
            if (isnan(Isyn_fs[j])){
                printf("Warning, Isyn_fs is a NaN %f, IAMPA = %f, IGABA = %f IDA_fs = %f\n", Isyn_fs[j], IAMPA_fs[j], IGABA_fs[j], IDA_fs[j]);
                printf("Warning, Vfs is a NaN \n");
                throw "Vfs is a NaN";
            }
            
            if (isnan(Igapfs[j])){
                printf("Warning, Igapfs is a NaN \n");
                throw "Igapfs is a NaN";
            }
            if (isinf(Igapfs[j])){
                printf("Warning, Igapfs is a inf \n");
                throw "Igapfs is a Inf";
            }
            
            // compute v
            qsum = x5_fs[j] - x4_fs[j] *((Isyn_fs[j]+Igapfs[j])-UFS);
            if (qsum >=0){
                // then root is real
                A = sqrt(qsum); 
                roots = (A-xb_fs[j]) / xa_fs[j];
                exp_tc = exp(-dt*A/_FS(j,0));
                C0 = exp_tc / (VFS-roots);
                z = -_FS(j,3) / A * (1-exp_tc) + C0;
                Vfs[j] = roots + 1/z;
            }
            else {
                // root is complex
                B = sqrt(fabs(qsum));
                Y = VFS - vp_fs[j];
                X = -dt*B/_FS(j,0);
                tmp = _FS(j,3)*Y/B;
                cx = cos(X); sx = sin(X);
                Vfs[j] = vp_fs[j] + (Y*cx + Y*tmp*sx - (B/(4*_FS(j,3)))*sx) / 
                            (0.5*(1+cx) + 2*tmp*sx - 2*tmp*tmp*(cx-1));
            }
            // Vfs[j] = VFS + dt * (_FS(j, 3) * ((VFS - _FS(j, 1) * (1-_FS(j, 10)*lambda1_fs)) * (VFS - _FS(j, 2))) - Ufs[j] + Igapfs[j] + Isyn_fs[j] ) / _FS(j, 0);
            
            if (isnan(Vfs[j])){
                printf("Warning, Isyn_fs = %f, and Igapfs = %f \n", Isyn_fs[j], Igapfs[j]);
                printf("Warning, Vfs is a NaN \n");
                throw "Vfs is a NaN";
            }
            
            // compute u
            if(VFS < _FS(j, 9)){
                U = 0;
            }
            else{
                tmp = VFS - _FS(j, 9);
                U = _FS(j, 5) * tmp * tmp * tmp;    // b*(v-vb)^3
            }
            
            Ufs[j] = U * (1- u_exp_tc_fs[j]) + UFS * u_exp_tc_fs[j];
            
//            // update Ufs
//            if (VFS < _FS(j, 9)){
//                Ufs[j] = UFS + dt * -_FS(j, 4)*UFS;
//            }
//            else if (VFS >= _FS(j, 9)){
//                U = pow((_FS(j, 5) * (VFS - _FS(j, 9))), 3);
//                Ufs[j] = UFS + dt * _FS(j, 4) * (_FS(j, 5) * pow((VFS - _FS(j, 9)), 3) - UFS);
//            }
            
            // ---------------------------------------------------------------
            // Check for spike events
            if (Vfs[j] >= _FS(j, 8)){
                VFS = _FS(j, 8);
                Vfs[j] = _FS(j, 6);
                Ufs[j] = Ufs[j] + _FS(j, 7);
                
                // ---------------------------------------------------------------
                // save the cell number and spike time
                if (FSspikecount < FSspikebuffer){
                    STfs[FSspikecount] = j;
                    STfs[FSspikecount+FSspikebuffer] = T;
                    FSspikecount++;
                }
                // check for overrun of the spike buffer
                else if (FSspikecount >= FSspikebuffer){
                    printf("Warning, exceeded FS spike buffer size \n");
                    throw "exceeded FS spike buffer size";
                }
                
                // ---------------------------------------------------------------
                // Update the MS spike event que
                for (int targcellind = Cfsms_b.data[j]; targcellind <= Cfsms_b.data[j+1]-1; targcellind++){ // for each of the target cells
                    // get the index in the spike que for the target cell, given its delay
                    SEQind_ms = Cfsms.data[targcellind] + ((SpikeEventCycle_MSGABA-1) + Cfsms_d.data[targcellind])*N_MS;
                    
                    if (SEQind_ms >= SpikeQueSize_MS){ // if SEind is larger than the spike que, wrap around
                        SEQind_ms = SEQind_ms - SpikeQueSize_MS;
                    }
                    
                    // add the spike event to the spike que for the target cell
                    SEQ_MSGABA[SEQind_ms] = SEQ_MSGABA[SEQind_ms] + (Cfsms_w.data[targcellind] / ts_gaba); // update the Spike event que
                }
                
                // ---------------------------------------------------------------
                // Update the FS spike event que
                for (int targcellind = Cfsfs_b.data[j]; targcellind <= Cfsfs_b.data[j+1]-1; targcellind++){ // for each of the target cells
                    // get the index in the spike que for the target cell, given its delay
                    SEQind_fs = Cfsfs.data[targcellind] + ((SpikeEventCycle_FSGABA-1) + Cfsfs_d.data[targcellind])*N_FS;
                    
                    // if SEind is larger than the spike que, wrap around
                    if (SEQind_fs >= SpikeQueSize_FS){
                        SEQind_fs = SEQind_fs - SpikeQueSize_FS;
                    }
                    
                    if ((Cfsfs_w.data[targcellind]) > 100){
                        printf("Warning, FS weight too big! %i %f \n", targcellind, Cfsfs_w.data[targcellind]);
                        throw "FS weight too big!";
                    }
                    
                    // add the spike event to the spike que for the target cell
                    SEQ_FSGABA[SEQind_fs] = SEQ_FSGABA[SEQind_fs] + (Cfsfs_w.data[targcellind] / ts_gaba);
                }
            }
        } // for each FS neurons
        
        // ===============================================================
        // update the MS PSP cycle position
        if (SpikeEventCycle_MSGABA < maxSEQdelay_MSGABA){
            SpikeEventCycle_MSGABA++;
        }
        else {
            SpikeEventCycle_MSGABA = 0;
        }
        
        // update the FS PSP cycle position
        if (SpikeEventCycle_FSGABA < maxSEQdelay_FSGABA){
            SpikeEventCycle_FSGABA++;
        }
        else {
            SpikeEventCycle_FSGABA = 0;
        }
        
        // ===============================================================
        // save some example neuron behaviour
        tout[t] = T;
        Vmsout[t] = Vms[0];
        Vfsout[t] = Vfs[0];
        
        for (int RecInd = 0; RecInd < RecordChan_MS.M; RecInd++){
            _RecordChan_MS_out(t, RecInd) = Isyn_ms[RecordChan_MS.data[RecInd]];
        }
        
        _RecordChan_MS_out(t, 0) = Vms[RecordChan_MS.data[1]];
        _RecordChan_MS_out(t, 1) = Ums[RecordChan_MS.data[1]];
        _RecordChan_MS_out(t, 2) = spk_ms[RecordChan_MS.data[1]];
        _RecordChan_MS_out(t, 3) = Gglu_ms_AMPA[RecordChan_MS.data[1]];
        _RecordChan_MS_out(t, 4) = Gglu_ms_NMDA[RecordChan_MS.data[1]];
        _RecordChan_MS_out(t, 5) = IAMPA_ms[RecordChan_MS.data[1]];
        _RecordChan_MS_out(t, 6) = INMDA_ms[RecordChan_MS.data[1]];
        _RecordChan_MS_out(t, 7) = Ggaba_ms[RecordChan_MS.data[1]];
        _RecordChan_MS_out(t, 8) = IGABA_ms[RecordChan_MS.data[1]];
        _RecordChan_MS_out(t, 9) = IDA_ms[RecordChan_MS.data[1]];
        _RecordChan_MS_out(t, 10) = Isyn_ms[RecordChan_MS.data[1]];
        
        _RecordChan_MS_out(t, 11) = Vms[RecordChan_MS.data[10]];
        _RecordChan_MS_out(t, 12) = Ums[RecordChan_MS.data[10]];
        _RecordChan_MS_out(t, 13) = spk_ms[RecordChan_MS.data[10]];
        _RecordChan_MS_out(t, 14) = Gglu_ms_AMPA[RecordChan_MS.data[10]];
        _RecordChan_MS_out(t, 15) = Gglu_ms_NMDA[RecordChan_MS.data[10]];
        _RecordChan_MS_out(t, 16) = IAMPA_ms[RecordChan_MS.data[10]];
        _RecordChan_MS_out(t, 17) = INMDA_ms[RecordChan_MS.data[10]];
        _RecordChan_MS_out(t, 18) = Ggaba_ms[RecordChan_MS.data[10]];
        _RecordChan_MS_out(t, 19) = IGABA_ms[RecordChan_MS.data[10]];
        _RecordChan_MS_out(t, 20) = IDA_ms[RecordChan_MS.data[10]];
        _RecordChan_MS_out(t, 21) = Isyn_ms[RecordChan_MS.data[10]];
        
        _RecordChan_MS_out(t, 22) = Vms[RecordChan_MS.data[15]];
        _RecordChan_MS_out(t, 23) = Ums[RecordChan_MS.data[15]];
        _RecordChan_MS_out(t, 24) = spk_ms[RecordChan_MS.data[15]];
        _RecordChan_MS_out(t, 25) = Gglu_ms_AMPA[RecordChan_MS.data[15]];
        _RecordChan_MS_out(t, 26) = Gglu_ms_NMDA[RecordChan_MS.data[15]];
        _RecordChan_MS_out(t, 27) = IAMPA_ms[RecordChan_MS.data[15]];
        _RecordChan_MS_out(t, 28) = INMDA_ms[RecordChan_MS.data[15]];
        _RecordChan_MS_out(t, 29) = Ggaba_ms[RecordChan_MS.data[15]];
        _RecordChan_MS_out(t, 30) = IGABA_ms[RecordChan_MS.data[15]];
        _RecordChan_MS_out(t, 31) = IDA_ms[RecordChan_MS.data[15]];
        _RecordChan_MS_out(t, 32) = Isyn_ms[RecordChan_MS.data[15]];
        
        _RecordChan_MS_out(t, 33) = Vms[RecordChan_MS.data[20]];
        _RecordChan_MS_out(t, 34) = Ums[RecordChan_MS.data[20]];
        _RecordChan_MS_out(t, 35) = spk_ms[RecordChan_MS.data[20]];
        _RecordChan_MS_out(t, 36) = Gglu_ms_AMPA[RecordChan_MS.data[20]];
        _RecordChan_MS_out(t, 37) = Gglu_ms_NMDA[RecordChan_MS.data[20]];
        _RecordChan_MS_out(t, 38) = IAMPA_ms[RecordChan_MS.data[20]];
        _RecordChan_MS_out(t, 39) = INMDA_ms[RecordChan_MS.data[20]];
        _RecordChan_MS_out(t, 40) = Ggaba_ms[RecordChan_MS.data[20]];
        _RecordChan_MS_out(t, 41) = IGABA_ms[RecordChan_MS.data[20]];
        _RecordChan_MS_out(t, 42) = IDA_ms[RecordChan_MS.data[20]];
        _RecordChan_MS_out(t, 43) = Isyn_ms[RecordChan_MS.data[20]];
        
    } // end of the simulation loop
} // end of the main simulation function


// ========================================================================
// Code to catch exceptions before they crash DCE!!!
// ========================================================================
#include <exception>
using namespace std;

void mexFunction(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[]) {
    try {
        // dodgy code
        execute(nlhs, plhs, nrhs, prhs);
    }
    catch(std::exception& e) {
        printf(e.what());
        return;
    }
    catch(const char* e) {
        printf(e);
        return;
    }
    catch(...) {
        //  report and close gracefully
        printf("an unexpected exception occurred\n");
        return;
    }
}
