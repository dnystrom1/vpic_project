//========================================================================
// Short Pulse Laser 3D
// Modified Jan 7 2017 for testing VPIC head version for trinity KNL open science II
// user defined diag. turned off
//
// no need to "mkdir log field hydro particle restart rundata"
//========================================================================

// ASCII Logging of IO writes to disk.  Since wall clock elapsed time is written,
// and these require an MPI_Allreduce, don't include this macro inside code
// which only executes on a single processor, else deadlock will ensue.
// Turn off logging by setting "#if 0"
  
//??????????????????????????????????????????????????????????????????????

#include "tracer_interp_one_pass.hxx"
#include <mpi.h>
#include <ctime>

#define NUM_TURNSTILES 36  

#if 0
#  define DIAG_LOG(MSG)                                          \
   {                                                             \
     FILE *fptmp=NULL;                                           \
     char fnametmp[256];                                         \
     sprintf( fnametmp, "log/diag_log.%i", rank() );             \
     if ( !(fptmp=fopen(fnametmp, "a")) ) ERROR(("Cannot open file %s", fnametmp));        \
     fprintf( fptmp, "At time step %i: %s\n", step(), MSG ); \
     fclose( fptmp );                                            \
   }
#else
#  define DIAG_LOG(MSG)
#endif


// DJS
# define RANK_TO_INDEX(rank,ix,iy,iz) BEGIN_PRIMITIVE {   \
  int _ix, _iy, _iz;    \
  _ix = (rank);         \
  _iy = _ix/int(topology_x) ;  \
  _ix -= _iy*int(topology_x) ;	       \
  _iz = _iy/int(topology_y) ;  \
  _iy -= _iz*int(topology_y) ;  \
  (ix) = _ix;    \
  (iy) = _iy;    \
  (iz) = _iz;    \
} END_PRIMITIVE


begin_globals {
  double emax;                   // E0 of the laser
  double omega_0;                // w0/wpe
  double vthe;                   // vthe/c   <- these are needed to make movie files
  double vthi_I2;                // vthi_I2/c
  double vthi_I1;                // vthi_I1/c
  int energies_interval;        // how frequently to dump energies
  int    field_interval;         // how frequently to dump field built-in diagnostic
  int    restart_interval; 	 // how frequently to write restart file. 
  int    quota_check_interval;   // how often to check if quote exceeded
  int    fft_ex_interval;        // how frequently to save ex fft data
  int    fft_ey_interval;        // how frequently to save ey fft data
  int    mobile_ions;	         // flag: 0 if ions are not to be pushed
  int    I1_present;             // flag nonzero when H ions are present. 
  int    I2_present;             // flag nonzero when He ions are present.  
  int    launch_wave;            // whether or not to propagate a laser from y-z boundary

  int    eparticle_interval;
  int    I1particle_interval;
  int    I2particle_interval;
  int    load_particles;         // Flag to turn off particle load for testing wave launch. 

  // Parameters for 2d and 3d Gaussian wave launch
  int    run_mode;               // determines how many dimensions, polarization
  double lambda;
  double waist;                  // how wide the focused beam is
  double width;
  double zcenter;		 // center of beam at boundary in z
  double ycenter;		 // center of beam at boundary in y
  double xfocus;                 // how far from boundary to focus
  double mask;			 // # gaussian widths from beam center where I nonzero
  double mask2;			 // # gaussian widths from beam center where I nonzero

  int    pulse_shape;            // 0 for steady, indefinite pulse, 1 for square pulse, 2 for sin^2
  double pulse_length;           // for pulse_shape=1, length of pulse. for pulse_shape=2, FWHM
  double pulse_FWHM;
  double pulse_start;
  double wpe1fs;

  double quota_sec;              // Run quota in seconds
  int    rtoggle;                // Enables save of last two restart dumps for safety


// particle tracking section
  int tracer_interval;
  int tracer2_interval;
  int particle_tracing;
  int particle_tracing_start;

  // particle tracking section
  species_t *tracers_list;

  int num_tracers_total;        // max number of tracers 
  int nbuf_tracer;              // number of buffered writes between file I/O 
  int nspec_tracer;             // how many tracer species 


  // Dump parameters for standard VPIC output formats
  DumpParameters fdParams;
  DumpParameters hedParams;
  DumpParameters hI1dParams;
  DumpParameters hI2dParams;
  std::vector<DumpParameters *> outputParams;

};

begin_initialization {

  // do not change parameters in this block:
  double elementary_charge  = 4.8032e-10;                                // stat coulomb
  double elementary_charge2 = elementary_charge * elementary_charge;    
  double speed_of_light     = 2.99792458e10;                             // cm/sec
  double m_e                = 9.1094e-28;                                // g
  double k_boltz            = 1.6022e-12;                                // ergs/eV
  double mec2               = m_e*speed_of_light*speed_of_light/k_boltz;
  double mpc2               = mec2*1836.0;
  double eps0               = 1;

  double cfl_req            = 0.98;        // How close to Courant should we try to run
  double damp               = 0.0;         // Level of radiation damping

// particle tracking section
  // for this run, we tracked every 4 particles in the beginning and then select interesting particles for the second run
  // For the first run particle_tracing=1, and particle_tracing=2 for the second run
 int particle_tracing = 0; // 0: notracing, 1: forward tracing 2: tracing from particle files
                            // should be 1 for 1-pass tracers
 int particle_tracing_start = 0; //the time step that particle tracking is triggered
                                  // this should be set to 0 for Pass1 and 2


  double n_e_over_n_crit    = 10.;       // n_e/n_crit in solid slab
  double laser_intensity    = 1e20;      // units of W/cm^2, not used for now
  double vacuum_wavelength  = 1000 * 1e-7; // 1/3 micron laser

//???????????????????????????????????????????????????????????????????????????

  // Box size for a single node.
  double box_size_x         = REPLACE_nx_sn * 30.0 * 1e-4 / 6000.0;  // microns
  double box_size_y         = REPLACE_nx_sn * 30.0 * 1e-4 / 6000.0;  // microns (ignored if 1d or 2d out of plane)
  double box_size_z         = REPLACE_nx_sn * 30.0 * 1e-4 / 6000.0;  // microns (ignored if 1d or 2d out of plane)

  // Scale box size for single node to adjust single node memory footprint.
  box_size_x              *= REPLACE_scale_Lx;
  box_size_y              *= REPLACE_scale_Ly;
  box_size_z              *= REPLACE_scale_Lz;

  // Scale box size for multiple nodes.
  box_size_x              *= REPLACE_scale_topology_x;
  box_size_y              *= REPLACE_scale_topology_y;
  box_size_z              *= REPLACE_scale_topology_z;

  // Grid size for a single node.
  double nx                = REPLACE_nx_sn;
  double ny                = REPLACE_ny_sn;
  double nz                = REPLACE_nz_sn;

  // Scale grid size for single node to adjust single node memory footprint.
  nx                      *= REPLACE_scale_Lx;
  ny                      *= REPLACE_scale_Ly;
  nz                      *= REPLACE_scale_Lz;

  // Scale grid size for multiple nodes.
  nx                      *= REPLACE_scale_topology_x;
  ny                      *= REPLACE_scale_topology_y;
  nz                      *= REPLACE_scale_topology_z;

  // Topology for a single node.
  double topology_x        = REPLACE_topology_x;
  double topology_y        = REPLACE_topology_y;
  double topology_z        = REPLACE_topology_z;

  // Scale topology for multiple nodes.
  topology_x              *= REPLACE_scale_topology_x;
  topology_y              *= REPLACE_scale_topology_y;
  topology_z              *= REPLACE_scale_topology_z;

  double FWHM               = 3.0 * 1e-4;         // from Wilks et al.
  double length_vac_l       = 10.0  * 1e-4;       // x length of vac on the left
  double length_Lp          = 1. * 1e-4;          // x length of slab dense plasma
  double length_ramp_l      = 0.0  * 1e-4;        // x length of leftmost density ramp
  double length_ramp_r      = 0.0 * 1e-4;         // x length of rightmost density ramp

//float dfrac               = 0.1932;             // fraction of charge density for n_Al12/ne
  float dfrac               = 0.0;                // fraction of charge density for n_Al12/ne

  double x0=length_vac_l;
  double x1=x0+length_ramp_l;
  double x2=x1+length_Lp;
  double x3=x2+length_ramp_r;

  int launch_wave           = 1;           // whether or not to launch wave from y-z plane
  int rng_seed              = 1;                        // random number generator seed.
  int run_mode = 3;                      // run mode=0, 1, 2, 3 means, respectively, 1D,
  int load_particles = 1;         // Flag to turn off particle load for testing wave launch. 
  //????????????????????????????????????????????????????????????????????
  double nppc        = REPLACE_nppc   ;          // Average number of macro particles/cell of each species
  int mobile_ions         = REPLACE_mobile_ions;           // whether or not to push ions

  //????????????????????????????????????????????????????????????????????
  double quota = 3.8;             // Run quota in hours.  
  double quota_sec = quota*3600;  // Run quota in seconds. 

  // Treat species 2 as being a heavy substrate, species 1 as a deposition layer.
  double A_I1      = 1;                 // proton
  double A_I2     = 12.0;             // carbon
  double Z_I1      = 1;
  double Z_I2     = 6;
  double mic2_I1   = mpc2*A_I1;
  double mic2_I2  = mpc2*A_I2;
  double mime_I1   = mic2_I1/mec2;
  double mime_I2  = mic2_I2/mec2;

  double t_e = 14.e3;                   // electron temp, eV
  double t_i = 10.;                      // ion temp, eV
  double uthe    = sqrt(t_e/mec2);    // vthe/c
  double uthi_I1 = sqrt(t_i/mic2_I1);  // vthi/c
  double uthi_I2 = sqrt(t_i/mic2_I2); // vthi/c

  double delta = (vacuum_wavelength/(2.0*M_PI))/sqrt(n_e_over_n_crit); // c/wpe
  double n_e = speed_of_light*speed_of_light*m_e/(4.0*M_PI*elementary_charge2*delta*delta);
  // n_e is used for emax only
  double debye = uthe*delta;

  double iv_thick          = 2;            // Thickness of impermeable vacuum (in cells)

  //double nx = 20968;  
  //??????????????????????????????????????????????????????????????????????????
  // double nx = 6000;  // for 7 instead of 16 micron
  // double ny = 1;    // was 1280;
  // double nz = 6000;

  //??????????????????????????????????????????????????????????????????????????
  // double topology_x = 6;
  // double topology_y = 1; // 18 cells per proc. domain; 36x36=1296 cores total
  // double topology_z = 6;

  double cell_size_t = box_size_y/(debye*ny);  // in debye
  double cell_size_tz= box_size_z/(debye*nz);
  double cell_size_l = box_size_x/(debye*nx);
  //double cell_size_l =0.216777;  // in debye

  double cells_within_Lp = length_Lp/(debye*cell_size_l);

  double hy = debye*cell_size_t/delta;
  double hz = debye*cell_size_tz/delta;
  double hx = debye*cell_size_l/delta;

  //double nx = box_size_x/(debye*cell_size_l);
  //nx = trunc_granular(nx, 1);
  //double hx=box_size_x/(delta*nx);
  //cell_size_l = box_size_x/(debye*nx);

  double Lx = nx*hx;          // in c/wpe
  double Ly = ny*hy;
  double Lz = nz*hz;

  // DEBUG: want to make this smaller ppc to get through the writes
  // double nppc      = 1       ;          // Average number of macro particles/cell of each species
  double particles_alloc = nppc*ny*nz*((0.5*length_ramp_l+length_Lp+0.5*length_ramp_r)/(cell_size_l*debye));

#if 0
  if ( run_mode==0 || run_mode==1 ) ny=1, Ly=hx; // y dir ignorable 
  if ( run_mode==0 || run_mode==2 ) nz=1, Lz=hx; // z dir ignorable
#endif

  double dt         = cfl_req*courant_length(Lx, Ly, Lz, nx, ny, nz);
  double dt_courant = dt;
  double t_stop     = REPLACE_nstep*dt + 0.001*dt; // Runtime in 1/wpe

  //int pulse_shape=2;                   // sin pulse
  //double wpe1fs = 1.0e-15*speed_of_light/delta;
  //double pulse_length=50.0*wpe1fs;             // in t*wpe, which is = 500 fs
  //double t_stop = 1000.0*wpe1fs;               // runtime in 1/omega_pe

  int pulse_shape=3;                   // sin pulse
  double wpe1fs = 1.0e-15*speed_of_light/delta;
  double pulse_length=2*150.*wpe1fs;   // in t*wpe, FWHM=650 fs
  double pulse_FWHM = 150e-15*speed_of_light / delta; // distance in code units
  // How far in front of the focus should the peak start?
  double pulse_start = 0; //+ 3*pulse_FWHM; // distance in code units
  // double t_stop = 1.*pulse_length;                // runtime in 1/omega_pe
//double t_stop = 0.0012*pulse_length;                // runtime in 1/omega_pe

  // Diagnostics intervals.  
  int energies_interval   = 200;
  int field_interval      = int(5.*wpe1fs/dt);
//int field_interval = int(300.0/dt); 
  int fft_ex_interval     = int(M_PI/(dt*1.0));       //  Num. steps between dumping Ey to resolve w/wpe=1
  int fft_ey_interval     = int(M_PI/(dt*0.4));       //  Num. steps between dumping Ey to resolve w/wpe=0.4
//?????????????????????????????????????????????????????????????????????????????
  int I1particle_interval = field_interval ;
  int I2particle_interval = field_interval ;
  int eparticle_interval  = field_interval;

//MORE TRACER STUFF
  int tracer_int        = int(0.5*wpe1fs/dt);
  int tracer_interval   = tracer_int;
  int tracer2_int       = int(0.5*wpe1fs/dt);
  int tracer2_interval  = tracer2_int;

  int nbuf_tracer       = 2;   // number of writes between file I/O for tracer data; which is 10 * tracer2_interval 
  int nspec_tracer      = 2;   // number of tracer species 

  // Both these intervals must be multiples of nbuf_tracer and the tracer
  // interval, or your tracers will quietly fail upon restart
  int restart_interval         = 25*nbuf_tracer*tracer_interval;   // num steps between restart dumps.
  int quota_check_interval = nbuf_tracer*tracer_interval;

                                                                   // must be an even multiple of 
                                                                   // nbuf_tracer*tracer_interval


  // BJA - total number of tracers per species should be <= 16,777,216 so the index can be stored in a float

  // Restrict tracers to defined region about the center of box. Assume 1D decomposition
  //
  double dy_tracer             = 15. * 1e-4;                       // in cm, radius of tracer region
  double dz_tracer             = 15. * 1e-4;                       // in cm, radius of tracer region

  int tracer_min_iy            = topology_y * ( 0.5 * box_size_y - dy_tracer ) / box_size_y;
  int tracer_max_iy            = topology_y * ( 0.5 * box_size_y + dy_tracer ) / box_size_y;
  int tracer_min_iz            = topology_z * ( 0.5 * box_size_z - dz_tracer ) / box_size_z;
  int tracer_max_iz            = topology_z * ( 0.5 * box_size_z + dz_tracer ) / box_size_z;

  if ( tracer_min_iy < 0          ) tracer_min_iy = 0;
  if ( tracer_max_iy >= topology_y ) tracer_max_iy = topology_y-1;
  if ( tracer_min_iz < 0          ) tracer_min_iz = 0;
  if ( tracer_max_iz >= topology_z ) tracer_max_iz = topology_z-1;

//??????????????????????????????????????????????????????????????????????????????????????????



  int num_tracers_per_proc     = 200;  // num_tracers_per_proc * nproc should be < 16e6
  //int num_tracers_total        = num_tracers_per_proc * nproc(); 
  int num_tracers_total        = num_tracers_per_proc
                                 * ( tracer_max_iy - tracer_min_iy + 1 )
                                 * ( tracer_max_iz - tracer_min_iz + 1 );
  if ( num_tracers_total > 16777216 )
    ERROR(("Too many tracers: %d > 16777216\n", num_tracers_total));


  double omega_0   = sqrt(1.0/n_e_over_n_crit);  // w0/wpe
  double intensity_cgs = 1e7*laser_intensity;    // [ergs/(s*cm^2)]
  double emax = 
    sqrt(2.0*intensity_cgs/
	 (m_e*speed_of_light*speed_of_light*speed_of_light*n_e)); // at focus, of NPIC

  double delta_0 = delta*sqrt(n_e_over_n_crit); // c/w0
  double Lp       = length_Lp    * Lx/box_size_x;
  double L_ramp_l = length_ramp_l* Lx/box_size_x;
  double L_ramp_r = length_ramp_r* Lx/box_size_x; 
  length_vac_l *= Lx/box_size_x;
  x0*=Lx/box_size_x;
  x1*=Lx/box_size_x;
  x2*=Lx/box_size_x;
  x3*=Lx/box_size_x;


  // laser focusing parameters
  double lambda    = vacuum_wavelength/delta;  // Wavelength in c/wpe
  double xfocus    = 10e-4*Lx/box_size_x; // in c/wpe
  double f_number  = 1.; // Not used!                    // f number of beam
  double waist     = lambda*1.25;  // in c/wpe, width of beam at focus
  double ycenter   = 0;         // spot centered in y on lhs boundary
  double zcenter   = 0;         // spot centered in z on lhs boundary
  double mask      = 2.4;       // Set drive I>0 if r>mask*width at boundary.
  double width = waist*sqrt(1+(lambda*xfocus/(M_PI*waist*waist))*(lambda*xfocus/(M_PI*waist*waist)));
// if plane wave:
  emax = emax*sqrt(waist/width); // at entrance if 2D Gaussian
  //emax = emax*(waist/width); // at entrance if 3D Gaussian

// these are not used
  FWHM = FWHM/delta;  // in c/wpe
  double mask2      = 0.5;       // Set drive I>0 if r>mask*width at boundary.
  double wa = 0.5*FWHM*0.5*FWHM/0.15;
  double wb = (lambda*xfocus/M_PI)*(lambda*xfocus/M_PI);
  double wy1 = 0.5*(wa + sqrt(wa*wa - 4.0*wb));
  double wy2 = 0.5*(wa - sqrt(wa*wa - 4.0*wb));
  double wy =( wy1>wy2 ? wy2 : wy1 );
//double waist = sqrt(wy);

  double Ne    = nppc*nx*ny*nz;             // Number of macro electrons in box
  Ne = trunc_granular(Ne, nproc());         // Make Ne divisible by number of processors       
  double Ni    = Ne;
  double Npe   = Lx*Ly*Lz;                  // Number of physical electrons in box, wpe = 1
  double qe    = -Npe/Ne;                   // Charge per macro electron
  double qi_I1 = -dfrac*qe;                 // Charge per macro ion of type 1. Note that species
  double qi_I2 = -(1.0-dfrac)*qe;           // I2 and I1 are separate from one another in the loading.
  if ( load_particles==0 ) Ne=Ni=98.7654;   // A weird number to signal that load_particles turned off. 

  int I1_present=0;
  int I2_present=1;

  // Print stuff that I need for plotters and such, and with enough sig figs!
  // Be very careful modifying this.  Plotters depend on explicit locations of
  // some of these numbers.  Generally speaking, add lines at the end only.
  if(rank() == 0){
    FILE * out;
    out = fopen("params.txt", "w");
    fprintf(out, "# Parameter file used for plotters.\n");
    fprintf(out, "%.14e   Time step (dt), code units\n", dt);
    fprintf(out, "%.14e   Laser wavelength, SI\n", vacuum_wavelength*1e-2);
    fprintf(out, "%.14e   Ratio of electron to critical density\n", n_e_over_n_crit);
    fprintf(out, "%d   Number of cells in x\n", int(nx));
    fprintf(out, "%d   Number of cells in y\n", int(ny));
    fprintf(out, "%d   Number of cells in z\n", int(nz));
    fprintf(out, "%.14e   Box size x, microns\n", box_size_x*1e4);
    fprintf(out, "%.14e   Box size y, microns\n", box_size_y*1e4);
    fprintf(out, "%.14e   Box size z, microns\n", box_size_z*1e4);
    fprintf(out, "%d   Field Interval\n", field_interval);
    fprintf(out, "%d   Tracer Interval\n", tracer_interval);
    fprintf(out, "%d   Number of tracers per species\n", num_tracers_total);
    fprintf(out, "%d   Number of tracer species\n", nspec_tracer);
    fprintf(out, "%d   Number of variables per tracer (possibly wrong)\n", TRACER_NUM_FIELDS_BASE + TRACER_NUM_ADDED_FIELDS);
    fprintf(out, "%d   Number of steps in the entire simulation\n", int(t_stop/(dt)));
    fprintf(out, "%d   This is my rank\n", rank());
    fclose(out);
  }
  
  // PRINT SIMULATION PARAMETERS 
    
  sim_log("***** Simulation parameters *****");
  sim_log("* Processors:                    "<<nproc());
  sim_log("* dt_courant, wpe1fs, cells_within_Lp =  "<<dt_courant<<" "<<wpe1fs<<" "<<cells_within_Lp);
  sim_log("* delta/dx,delta/dz =  "<<1.0/hx<<" "<<1.0/hz);
  sim_log("* Time step, max time, nsteps =  "<<dt<<" "<<t_stop<<" "<<int(t_stop/(dt)));
  sim_log("* Debye length,cell size_l,cell size_t,delta,delta_0 = "<<debye<<" "<<cell_size_l<<" "<<cell_size_t<<" "<<delta<<" "<<delta_0);
  sim_log("* Lx, Ly, Lz =                   "<<Lx<<" "<<Ly<<" "<<Lz);
  sim_log("* nx, ny, nz =                   "<<nx<<" "<<ny<<" "<<nz);
  sim_log("* Lp, L_ramp_l, L_ramp_r =       "<<Lp<<" "<<L_ramp_l<<" "<<L_ramp_r);
  sim_log("* x0, x1, x2, x3 =               "<<x0<<" "<<x1<<" "<<x2<<" "<<x3);
  sim_log("* Charge/macro electron =        "<<qe);
  sim_log("* Charge/macro I2 =              "<<qi_I2);
  sim_log("* Charge/macro I1 =              "<<qi_I1);
  sim_log("* particles_alloc =              "<<particles_alloc);
  sim_log("* Average particles/processor:   "<<Ne/nproc());
  sim_log("* Average particles/cell:        "<<nppc);
  sim_log("* Do we have mobile ions?        "<<(mobile_ions ? "Yes" : "No"));
  sim_log("* Are we launching a laser?      "<<(launch_wave ? "Yes" : "No"));
  sim_log("* Omega_0, Omega_pe:             "<<(omega_0)<<" "<<1);
  sim_log("* Plasma density, ne/nc:         "<<n_e<<" "<<n_e_over_n_crit);
  sim_log("* Vac wavelength, I_laser:       "<<vacuum_wavelength<<" "<<laser_intensity);
  sim_log("* T_e, T_i, m_e, m_i_I1, m_i_I2:  "<<t_e<<" "<<t_i<<" "<<1<<" "<<mime_I1<<" "<<mime_I2);
  sim_log("* Radiation damping:             "<<damp);
  sim_log("* Fraction of courant limit:     "<<cfl_req);
  sim_log("* vthe/c:                        "<<uthe);
  sim_log("* vthi_I1/c, vth_I2/c:            "<<uthi_I1<<" "<<uthi_I2);
  sim_log("* emax at entrance, emax at waist:  "<<emax<<" "<<emax/sqrt(waist/width));
  sim_log("* energies_interval:             "<<energies_interval);
  sim_log("* field_interval:                "<<field_interval);
  sim_log("* restart interval: "             <<restart_interval);
  sim_log("* ex save interval:              "<<fft_ex_interval);
  sim_log("* ey save interval:              "<<fft_ey_interval);
  sim_log("* random number base seed:       "<<rng_seed);
  sim_log("* waist, width, xfocus:          "<<waist<<" "<<width<<" "<<xfocus);
  sim_log("* ycenter, zcenter, mask:        "<<ycenter<<" "<<zcenter<<" "<<mask);
  sim_log("* wy1, wy2, wy:                  "<<wy1<<" "<<wy2<<" "<<wy);
  sim_log("* tracer_interval:               "<<tracer_interval);
  sim_log("* tracer2_interval:              "<<tracer2_interval); 
  sim_log("*********************************");

  // SETUP HIGH-LEVEL SIMULATION PARMETERS
  // FIXME : proper normalization in these units for: xfocus, ycenter, zcenter, waist
  sim_log("Setting up high-level simulation parameters. "); 
  num_step             = int(t_stop/(dt)); 
  status_interval      = REPLACE_status_interval; 
//????????????????????????????????????????????????????????????????????????????????
//sync_shared_interval = status_interval/10;
//clean_div_e_interval = status_interval/10;
//clean_div_b_interval = status_interval/10;
  sync_shared_interval = REPLACE_sync_shared_interval;
  clean_div_e_interval = REPLACE_clean_div_e_interval;
  clean_div_b_interval = REPLACE_clean_div_b_interval;
  verbose = 0;
  global->energies_interval        = energies_interval;
  global->field_interval           = field_interval; 
  global->restart_interval         = restart_interval;
  global->quota_check_interval = quota_check_interval;
  global->fft_ex_interval     = fft_ex_interval;
  global->fft_ey_interval     = fft_ey_interval;
  global->vthe                     = uthe;     // c=1
  global->vthi_I2                  = uthi_I2;  // c=1
  global->vthi_I1                  = uthi_I1;  // c=1
  global->emax                     = emax; 
  global->omega_0                  = omega_0;
  global->mobile_ions              = mobile_ions; 
  global->I1_present                = I1_present;
  global->I2_present               = I2_present;
  global->launch_wave              = launch_wave; 
  global->lambda                   = lambda;
  global->waist                    = waist;
  global->width                    = width;
  global->mask                     = mask;
  global->mask2                    = mask2;
  global->xfocus                   = xfocus; 
  global->ycenter                  = ycenter; 
  global->zcenter                  = zcenter; 
  global->run_mode                 = run_mode; 
  global->wpe1fs                   = wpe1fs;

  global->quota_sec                = quota_sec;
  global->rtoggle                  = 0;

  global->eparticle_interval          = eparticle_interval; 
  global->I1particle_interval           = I1particle_interval; 
  global->I2particle_interval           = I2particle_interval; 
  global->load_particles           = load_particles; 

  global->pulse_shape              = pulse_shape; 
  global->pulse_length             = pulse_length; 
  global->pulse_FWHM               = pulse_FWHM;
  global->pulse_start              = pulse_start;

  
  // particle tracking
  global->particle_tracing         = particle_tracing;
  global->tracer_interval          = tracer_interval;
  global->tracer2_interval         = tracer2_interval;

  global->num_tracers_total        = num_tracers_total;
  global->nbuf_tracer              = nbuf_tracer;
  global->nspec_tracer             = nspec_tracer;


  // SETUP THE GRID
  sim_log("Setting up computational grid."); 
  grid->dx = hx;
  grid->dy = hy;
  grid->dz = hz;
  grid->dt = dt;
  grid->cvac = 1;
  grid->eps0 = eps0;

  // Partition a periodic box among the processors sliced uniformly in z: 
  define_absorbing_grid( 0,      -0.5*Ly,  -0.5*Lz,    // Low corner
                        Lx,      0.5*Ly,   0.5*Lz,    // High corner
                        nx,      ny,       nz,        // Resolution
                        topology_x, topology_y, topology_z, // Topology
                        reflect_particles );             // Default particle boundary condition

// for plane wave:
#if 0
  define_periodic_grid( 0,      -0.5*Ly,  -0.5*Lz,    // Low corner
                        Lx,      0.5*Ly,   0.5*Lz,    // High corner
                        nx,      ny,       nz,        // Resolution
                        topology_x, topology_y, topology_z); // Topology
    set_domain_field_bc( BOUNDARY(-1,0,0), absorb_fields );
    set_domain_field_bc( BOUNDARY( 1,0,0), absorb_fields );
#endif

// X boundary conditions:
//  set_domain_field_bc( BOUNDARY(-1,0,0), absorb_fields );
//  set_domain_field_bc( BOUNDARY( 1,0,0), absorb_fields );

//  set_domain_particle_bc( BOUNDARY(-1,0,0), absorb_particles );
//  set_domain_particle_bc( BOUNDARY( 1,0,0), absorb_particles );

  // Z boundary conditions: reflecting particles, absorbing fields
  // Parallel decomposition in Z direction assumed. 
#if 0
  if ( rank()==nproc()-1 ) {
    set_domain_field_bc( BOUNDARY(0,0, 1), absorb_fields );
    set_domain_particle_bc( BOUNDARY(0,0, 1), absorb_particles );
  }
  if ( rank()==0 ) {
    set_domain_field_bc( BOUNDARY(0,0,-1), absorb_fields );
    set_domain_particle_bc( BOUNDARY(0,0,-1), absorb_particles );
  }
#endif

  // SETUP THE SPECIES - N.B. OUT OF ORDER WITH GRID SETUP IN CANONICAL ORDERING

  // Allow 10% more local_particles in case of non-uniformity.
  // Sort interval is each 20 steps for electrons, 40 for ions.

  // Also, while we install the species id data in the maxwellian reflux
  // particle boundary condition data.

  sim_log("Setting up electrons. ");

#if 0
  double num_particles_per_proc = 2.0*particles_alloc/nproc();
  species_t * electron = define_species("electron", -1, 1, num_particles_per_proc, -1, 20, 1);
  species_t *ion_I1, *ion_I2;
  if ( mobile_ions ) {
  sim_log("Setting up ions. ");
    if ( I1_present  ) ion_I1 = define_species("I1", Z_I1, mime_I1, num_particles_per_proc, -1, 80, 1);
    if ( I2_present  ) ion_I2 = define_species("I2", Z_I2, mime_I2, num_particles_per_proc, -1, 80, 1);
  }
#endif

//???????????????????????????????????????????????????????????????????????
  double max_local_np_e            = 8.0*particles_alloc/nproc();
  double max_local_np_i1            = max_local_np_e;
  double max_local_np_i2            = max_local_np_e;
  double max_local_nm_e            = max_local_np_e / 10.0;
  double max_local_nm_i1            = max_local_nm_e;
  double max_local_nm_i2            = max_local_nm_e;
  species_t * electron = define_species("electron", -1, 1, max_local_np_e, max_local_nm_e, REPLACE_eon_sort_interval, REPLACE_eon_sort_method);
  species_t *ion_I1, *ion_I2;
  if ( mobile_ions ) {
  sim_log("Setting up ions. ");
    if ( I1_present  ) ion_I1 = define_species("I1", Z_I1, mime_I1, max_local_np_i1, max_local_nm_i1, REPLACE_ion_sort_interval, REPLACE_ion_sort_method);
    if ( I2_present  ) ion_I2 = define_species("I2", Z_I2, mime_I2, max_local_np_i2, max_local_nm_i2, REPLACE_ion_sort_interval, REPLACE_ion_sort_method);
  }



 // Add tracers
 
  int num_tracer_particles_per_proc = 3.0*num_tracers_per_proc;

  //  sim_log( "num_particles_per_proc        = "<<num_particles_per_proc );
  sim_log( "num_tracer_particles_per_proc = "<<num_tracer_particles_per_proc );

  sim_log("Setting up tracer electrons.");
  species_t * e_tracer   = define_species("electron_tracer", -1, 1,  num_tracer_particles_per_proc, -1, -1, 0 );

  sim_log("Setting up tracer ions species 2");
  species_t * I2_tracer  = define_species("I2_tracer",  Z_I2,mime_I2, num_tracer_particles_per_proc, -1, -1, 0 );

  hijack_tracers(nspec_tracer);

  sim_log("defined tracer species.");




  // SETUP THE MATERIALS

  sim_log("Setting up materials. "); 
  define_material( "vacuum", 1 );
  // define_material( "impermeable_vacuum", 1 );
//define_material( "impermeable_vacuum_xr", 1 );
  define_field_array( NULL, damp ); 

  // Paint the simulation volume with materials and boundary conditions
# define iv_region (   x<hx*iv_thick || x>Lx-hx*iv_thick                                   \
		       /* || y<-Ly/2+hy*iv_thick || y>Ly/2-hy*iv_thick*/  \
                    || z<-Lz/2+hz*iv_thick || z>Lz/2-hz*iv_thick ) /* all boundaries are i.v. */


#if 0
// for 2D
# define iv_region (   x<hx*iv_thick || x>Lx-hx*iv_thick                               \
                /*    || y<-Ly/2+hy*iv_thick || y>Ly/2-hy*iv_thick     */              \
                    || z<-Lz/2+hz*iv_thick || z>Lz/2-hz*iv_thick ) /* all boundaries are i.v. */

// for plane wave:
# define iv_region ( x<hx*iv_thick || x>Lx -hx*iv_thick )
#endif

  // set_region_material( iv_region, "impermeable_vacuum",     "impermeable_vacuum" );
//set_region_bc( iv_region, reflect_particles, reflect_particles, reflect_particles);
  set_region_bc( iv_region, absorb_particles, absorb_particles, absorb_particles);

#if 0
# define iv_region_xr (  x>Lx-hx*iv_thick )
//set_region_material( iv_region_xr, "impermeable_vacuum_xr",     "impermeable_vacuum_xr" );
  set_region_bc( iv_region_xr, absorb_particles, absorb_particles, absorb_particles);
#endif

  // LOAD PARTICLES

  // x0, x1 are positions of left and right edge of left ramp.
  // x2, x3 are positions of left and right edge of right ramp.

  if ( x1<x0 ) { x1=x0; sim_log("WARNING: x1<x0 encountered--setting x1=x0."); }
  if ( x2<x1 ) { x2=x1; sim_log("WARNING: x2<x1 encountered--setting x2=x1."); }
  if ( x3<x2 ) { x3=x2; sim_log("WARNING: x3<x2 encountered--setting x3=x2."); }


// particle tracking
  int itp=0;
  // particle tracking
  int ip   = 0;                           // particle load counter

  // Load particles using rejection method (p. 290 Num. Recipes in C 2ed, Press et al.)  

  if ( load_particles!=0 ) {
    sim_log( "Loading particles" );
  
    // Fast load of particles
    double xmin = x0;
    double xmax = x3;
    double ymin = grid->y0;
    double ymax = grid->y0+grid->ny*grid->dy;
    double zmin = grid->z0;
    double zmax = grid->z0+grid->nz*grid->dz;
 

    
    int particle_select = (Ne*(x2-x1)/Lx)/(topology_x*topology_y*topology_z) / num_tracers_per_proc;
    //    particle_select = 4;

    sim_log("particle_select " << particle_select);
    // iy_ptag, iz_ptag are used to build up the particle index
    // ix, iy, iz are used to determine whether to make tracer particles 
    int ix, iy, iz, iy_ptag, iz_ptag;
    RANK_TO_INDEX( rank(), ix, iy, iz );
    iy_ptag = iy - tracer_min_iy;
    iz_ptag = iz - tracer_min_iz;





//int blah =  (Ne*(x3-x0)/Lx)/(topology_x*topology_y*topology_z);
//std::cout <<" blah" << blah << " xmin " << xmin << " xmax " << xmax << std::endl;
//std::cout << "hx iv thick " << hx*iv_thick << " lz .. " << Lz/2+hz*iv_thick << std::endl;
//std::cout << "x lower " << hx*iv_thick << "x upper " << Lx-hx*iv_thick << std::endl;
//std::cout << "y lower " << -Ly/2 + hy*iv_thick << "y upper " << Ly/2 -hy*iv_thick << std::endl;
// std::cout << "z lower " << -Lz/2 + hz*iv_thick << "z upper " << Lz/2 -hz*iv_thick << std::endl;
// std::cout << "xmin " << xmin << " xmax "<<xmax<< " ymin "<<ymin<< " ymax "<<ymax<< " zmin " << zmin << " zmax "<<zmax<<std::endl;
 
    repeat( (Ne*(x3-x0)/Lx)/(topology_x*topology_y*topology_z) ) {
      double x = uniform( rng(0), xmin, xmax );     // Sample only from x=x0 to x=x3.
      double y = uniform( rng(0), ymin, ymax );   
      double z = uniform( rng(0), zmin, zmax );
     // std::cout << "iv_region " << iv_region << std::endl;   
      if ( iv_region ) continue;           // Particle fell in iv_region.  Don't load.
  
      // Flat distribution from x1<x<x2, linear ramps: x0<x<x1 and x2<x<x3
      if ( (x>=x1 && x<x2) ) {                      // substrate
//??????????????????????????????????????????????????????????????????????????????????
      // third to last arg is "weight," a positive number
          //std::cout<< " injecting electron " << std::endl;
          inject_particle( electron, x, y, z,
                         normal( rng(0), 0, uthe ),
                         normal( rng(0), 0, uthe ),
                         normal( rng(0), 0, uthe ), fabs(qe), 0, 0 );
  
//        if ( mobile_ions )
//          inject_particle( ion_I2, x, y, z,
//                           normal( rng(0), 0, uthi_I2 ),
//                           normal( rng(0), 0, uthi_I2 ),
//                           normal( rng(0), 0, uthi_I2 ), qi_I2 );

        if ( mobile_ions ) {
#if 0
            inject_particle( ion_I1, x, y, z,
                             normal( rng(0), 0, uthi_I1 ),
                             normal( rng(0), 0, uthi_I1 ),
                             normal( rng(0), 0, uthi_I1 ), fabs(qi_I1)/Z_I1, 0, 0 );
#endif
  
//??????????????????????????????????????????????????????????????????????????????????
            inject_particle( ion_I2, x, y, z,
                             normal( rng(0), 0, uthi_I2 ),
                             normal( rng(0), 0, uthi_I2 ),
                             normal( rng(0), 0, uthi_I2 ), fabs(qi_I2)/Z_I2, 0, 0 );
        }
 

        // particle tracking
        ++ip;
	//std::cout << "Particle Tracing " << particle_tracing << " iy " << iy << " iy_min " <<  tracer_min_iy << " particle_select " << particle_select << " num_tracers_per_proc " << num_tracers_per_proc << " itp " << itp << " ip " << ip << std::endl;
        if ( particle_tracing == 1
             && iy >= tracer_min_iy
             && iy <= tracer_max_iy
             && iz >= tracer_min_iz
             && iz <= tracer_max_iz ) {
          if ( ip%particle_select == 0 ) {
            if ( itp<num_tracers_per_proc ) {
              printf("itp %d \n",itp);
              float ptag = ( iz_ptag * ( tracer_max_iy - tracer_min_iy + 1 ) + iy_ptag )
                           * num_tracers_per_proc + itp; // exact up to ptag=16777216
              tag_tracer( (electron->p + electron->np-1), e_tracer,  ptag );
              tag_tracer( (ion_I2->p + ion_I2->np-1),     I2_tracer, ptag );
              itp++;

            } // if ( itp<num_tracers_per_proc )
          } // if ( ip%particle_select == 0 )
        } // if ( particle_tracing == 1 )



 
      }
  
  
#if 0
      if (( x>=x0 && x<x1 && uniform( rng(0), 0, 1 )<(x-x0)/(x1-x0) )) { // left ramp
        inject_particle( electron, x, y, z,
                         normal( rng(0), 0, uthe ),
                         normal( rng(0), 0, uthe ),
                         normal( rng(0), 0, uthe ), qe );
        if ( mobile_ions )
          inject_particle( ion_I2, x, y, z,
                           normal( rng(0), 0, uthi_I2 ),
                           normal( rng(0), 0, uthi_I2 ),
                           normal( rng(0), 0, uthi_I2 ), qi_I2 );
      }
 
      if ( ( x>=x2 && x<x3 && uniform( rng(0), 0, 1 )<(x3-x)/(x3-x2) )) { // right ramp
        inject_particle( electron, x, y, z,
                         normal( rng(0), 0, uthe ),
                         normal( rng(0), 0, uthe ),
                         normal( rng(0), 0, uthe ), qe );
        if ( mobile_ions )
          inject_particle( ion_I1, x, y, z,
                           normal( rng(0), 0, uthi_I1 ),
                           normal( rng(0), 0, uthi_I1 ),
                           normal( rng(0), 0, uthi_I1 ), qi_I1 );
      }
#endif

    }
 
    
    if ( particle_tracing==2 ) {
      exit(0);  // Not enabled for 1-pass processing 
    } // if ( particle_tracing==2 ) 


 } // if load_particles

 /*--------------------------------------------------------------------------
  * New dump definition
  *------------------------------------------------------------------------*/

 /*--------------------------------------------------------------------------
  * Set data output format
  * 
  * This option allows the user to specify the data format for an output
  * dump.  Legal settings are 'band' and 'band_interleave'.  Band-interleave
  * format is the native storage format for data in VPIC.  For field data,
  * this looks something like:
  * 
  *   ex0 ey0 ez0 div_e_err0 cbx0 ... ex1 ey1 ez1 div_e_err1 cbx1 ...
  *   
  * Banded data format stores all data of a particular state variable as a
  * contiguous array, and is easier for ParaView to process efficiently. 
  * Banded data looks like:
  * 
  *   ex0 ex1 ex2 ... exN ey0 ey1 ey2 ...
  *   
  *------------------------------------------------------------------------*/
  sim_log("Setting up hydro and field diagnostics.");

  global->fdParams.format = band;
  sim_log ( "Field output format          : band" );

  global->hedParams.format = band;
  sim_log ( "Electron hydro output format : band" );

  global->hI1dParams.format = band;
  sim_log ( "I1 hydro output format : band" );

  global->hI2dParams.format = band;
  sim_log ( "I2 hydro output format   : band" );


 /*--------------------------------------------------------------------------
  * Set stride
  * 
  * This option allows data down-sampling at output.  Data are down-sampled
  * in each dimension by the stride specified for that dimension.  For
  * example, to down-sample the x-dimension of the field data by a factor
  * of 2, i.e., half as many data will be output, select:
  * 
  *   global->fdParams.stride_x = 2;
  *
  * The following 2-D example shows down-sampling of a 7x7 grid (nx = 7,
  * ny = 7.  With ghost-cell padding the actual extents of the grid are 9x9.
  * Setting the strides in x and y to equal 2 results in an output grid of
  * nx = 4, ny = 4, with actual extents 6x6.
  *
  * G G G G G G G G G
  * G X X X X X X X G
  * G X X X X X X X G         G G G G G G
  * G X X X X X X X G         G X X X X G
  * G X X X X X X X G   ==>   G X X X X G
  * G X X X X X X X G         G X X X X G
  * G X X X X X X X G         G X X X X G
  * G X X X X X X X G         G G G G G G
  * G G G G G G G G G
  *
  * Note that grid extents in each dimension must be evenly divisible by
  * the stride for that dimension:
  *
  *   nx = 150;
  *   global->fdParams.stride_x = 10; // legal -> 150/10 = 15
  *
  *   global->fdParams.stride_x = 8; // illegal!!! -> 150/8 = 18.75
  *------------------------------------------------------------------------*/

  // Strides for field and hydro arrays.  Note that here we have defined them 
  // the same for fields and all hydro species; if desired, we could use different
  // strides for each.   Also note that strides must divide evenly into the number 
  // of cells in a given domain. 

  // Define strides and test that they evenly divide into grid->nx, ny, nz
//int stride_x = 2, stride_y = 4, stride_z = 4; 
//??????????????????????????????????????????????????????????????????????
//int stride_x = 4, stride_y = 4, stride_z = 4;
//int stride_x = 4, stride_y = 3, stride_z = 3;
  int stride_x = 1, stride_y = 1, stride_z = 1;
//int stride_x = 1, stride_y = 3, stride_z = 3;
  if ( int(grid->nx)%stride_x ) ERROR(("Stride doesn't evenly divide grid->nx."));
  if ( int(grid->ny)%stride_y ) ERROR(("Stride doesn't evenly divide grid->ny."));
  if ( int(grid->nz)%stride_z ) ERROR(("Stride doesn't evenly divide grid->nz."));

  //----------------------------------------------------------------------
  // Fields

  // relative path to fields data from global header
  sprintf(global->fdParams.baseDir, "field");

  // base file name for fields output
  sprintf(global->fdParams.baseFileName, "fields");

  // set field strides
  global->fdParams.stride_x = stride_x;
  global->fdParams.stride_y = stride_y;
  global->fdParams.stride_z = stride_z;
  sim_log ( "Fields x-stride " << global->fdParams.stride_x );
  sim_log ( "Fields y-stride " << global->fdParams.stride_y );
  sim_log ( "Fields z-stride " << global->fdParams.stride_z );

  // add field parameters to list
  global->outputParams.push_back(&global->fdParams);

  //----------------------------------------------------------------------
  // Electron hydro

  // relative path to electron species data from global header
  sprintf(global->hedParams.baseDir, "ehydro");

  // base file name for fields output
  sprintf(global->hedParams.baseFileName, "e_hydro");

  // set electron hydro strides
  global->hedParams.stride_x = stride_x;
  global->hedParams.stride_y = stride_y;
  global->hedParams.stride_z = stride_z;
  sim_log ( "Electron species x-stride " << global->hedParams.stride_x );
  sim_log ( "Electron species y-stride " << global->hedParams.stride_y );
  sim_log ( "Electron species z-stride " << global->hedParams.stride_z );

  // add electron hydro parameters to list
  global->outputParams.push_back(&global->hedParams);

  //----------------------------------------------------------------------
  // ion I1 hydro

  // relative path to electron species data from global header
  sprintf(global->hI1dParams.baseDir, "I1hydro");

  // base file name for fields output
  sprintf(global->hI1dParams.baseFileName, "I1_hydro");

  // set hydrogen hydro strides
  global->hI1dParams.stride_x = stride_x;
  global->hI1dParams.stride_y = stride_y;
  global->hI1dParams.stride_z = stride_z;
  sim_log ( "Ion species x-stride " << global->hI1dParams.stride_x );
  sim_log ( "Ion species y-stride " << global->hI1dParams.stride_y );
  sim_log ( "Ion species z-stride " << global->hI1dParams.stride_z );

  // add hydrogen hydro parameters to list
  global->outputParams.push_back(&global->hI1dParams);

  //----------------------------------------------------------------------
  // ion I2 hydro

  // relative path to electron species data from global header
  sprintf(global->hI2dParams.baseDir, "I2hydro");

  // base file name for fields output
  sprintf(global->hI2dParams.baseFileName, "I2_hydro");

  // set helium hydro strides
  global->hI2dParams.stride_x = stride_x;
  global->hI2dParams.stride_y = stride_y;
  global->hI2dParams.stride_z = stride_z;
  sim_log ( "Ion species x-stride " << global->hI2dParams.stride_x );
  sim_log ( "Ion species y-stride " << global->hI2dParams.stride_y );
  sim_log ( "Ion species z-stride " << global->hI2dParams.stride_z );

  // add helium hydro parameters to list
  global->outputParams.push_back(&global->hI2dParams);

 /*-----------------------------------------------------------------------
  * Set output fields
  *
  * It is now possible to select which state-variables are output on a
  * per-dump basis.  Variables are selected by passing an or-list of
  * state-variables by name.  For example, to only output the x-component
  * of the electric field and the y-component of the magnetic field, the
  * user would call output_variables like:
  *
  *   global->fdParams.output_variables( ex | cby );
  *
  * NOTE: OUTPUT VARIABLES ARE ONLY USED FOR THE BANDED FORMAT.  IF THE
  * FORMAT IS BAND-INTERLEAVE, ALL VARIABLES ARE OUTPUT AND CALLS TO
  * 'output_variables' WILL HAVE NO EFFECT.
  *
  * ALSO: DEFAULT OUTPUT IS NONE!  THIS IS DUE TO THE WAY THAT VPIC
  * HANDLES GLOBAL VARIABLES IN THE INPUT DECK AND IS UNAVOIDABLE.
  *
  * For convenience, the output variable 'all' is defined:
  *
  *   global->fdParams.output_variables( all );
  *------------------------------------------------------------------------*/
 /* CUT AND PASTE AS A STARTING POINT
  * REMEMBER TO ADD APPROPRIATE GLOBAL DUMPPARAMETERS VARIABLE

   output_variables( all );

   output_variables( electric | div_e_err | magnetic | div_b_err |
                     tca      | rhob      | current  | rhof |
                     emat     | nmat      | fmat     | cmat );

   output_variables( current_density  | charge_density |
                     momentum_density | ke_density     | stress_tensor );
  */

  //global->fdParams.output_variables( all );
  global->fdParams.output_variables( electric | magnetic );

  //global->hedParams.output_variables( all );
  //global->hedParams.output_variables( current_density | momentum_density );
  global->hedParams.output_variables(  current_density  | charge_density |
                                       momentum_density | ke_density |
                                       stress_tensor );
  global->hI1dParams.output_variables(  current_density  | charge_density |
                                       momentum_density | ke_density |
                                       stress_tensor );
  global->hI2dParams.output_variables( current_density  | charge_density |
                                       momentum_density | ke_density |
                                       stress_tensor );

 /*--------------------------------------------------------------------------
  * Convenience functions for simlog output
  *------------------------------------------------------------------------*/
  char varlist[256];

  create_field_list(varlist, global->fdParams);
  sim_log ( "Fields variable list: " << varlist );

  create_hydro_list(varlist, global->hedParams);
  sim_log ( "Electron species variable list: " << varlist );

  create_hydro_list(varlist, global->hI1dParams);
  sim_log ( "I1 species variable list: " << varlist );

  create_hydro_list(varlist, global->hI2dParams);
  sim_log ( "I2 species variable list: " << varlist );

 /*------------------------------------------------------------------------*/


  sim_log("*** Finished with user-specified initialization ***"); 

  // Upon completion of the initialization, the following occurs:
  // - The synchronization error (tang E, norm B) is computed between domains
  //   and tang E / norm B are synchronized by averaging where discrepancies
  //   are encountered.
  // - The initial divergence error of the magnetic field is computed and
  //   one pass of cleaning is done (for good measure)
  // - The bound charge density necessary to give the simulation an initially
  //   clean divergence e is computed.
  // - The particle momentum is uncentered from u_0 to u_{-1/2}
  // - The user diagnostics are called on the initial state
  // - The physics loop is started
  //
  // The physics loop consists of:
  // - Advance particles from x_0,u_{-1/2} to x_1,u_{1/2}
  // - User particle injection at x_{1-age}, u_{1/2} (use inject_particles)
  // - User current injection (adjust field(x,y,z).jfx, jfy, jfz)
  // - Advance B from B_0 to B_{1/2}
  // - Advance E from E_0 to E_1
  // - User field injection to E_1 (adjust field(x,y,z).ex,ey,ez,cbx,cby,cbz)
  // - Advance B from B_{1/2} to B_1
  // - (periodically) Divergence clean electric field
  // - (periodically) Divergence clean magnetic field
  // - (periodically) Synchronize shared tang e and norm b
  // - Increment the time step
  // - Call user diagnostics
  // - (periodically) Print a status message
}


begin_diagnostics {
//int mobile_ions=global->mobile_ions, 
//    I1_present=global->I1_present,
//    I2_present=global->I2_present;

  //if ( step()%10==0 ) sim_log("Time step: "<<step); 

# define should_dump(x) \
  (global->x##_interval>0 && remainder(step(),global->x##_interval)==0)

  if ( step()==0 ) {
    // A grid dump contains all grid parameters, field boundary conditions,
    // particle boundary conditions and domain connectivity information. This
    // is stored in a binary format. Each rank makes a grid dump
    //dump_grid("grid");

    // A materials dump contains all the materials parameters. This is in a
    // text format. Only rank 0 makes the materials dump
    //dump_materials("materials");

    // A species dump contains the physics parameters of a species. This is in
    // a text format. Only rank 0 makes the species dump
    //dump_species("species");

    if ( rank()==0 ) {
    dump_mkdir("rundata");
    dump_mkdir("fft");
    dump_mkdir("field");
    dump_mkdir("ehydro");
    dump_mkdir("I1hydro");
    dump_mkdir("I2hydro");
    dump_mkdir("restart");
    dump_mkdir("particle");

    //Also necessary for tracers
    dump_mkdir("restart/restart0");
    dump_mkdir("restart/restart1");

    // Turn off rundata for now
    // dump_grid("rundata/grid");
    // dump_materials("rundata/materials");
    // dump_species("rundata/species");
    global_header("global", global->outputParams);

    } // if 

  }

  // energy in various fields/particles 
  if( should_dump(energies) ) {
            dump_energies( "rundata/energies", step() ==0 ? 0 : 1 );
  } //if

  if ( should_dump(field) ) {
    field_dump( global->fdParams );

    if ( global->load_particles ) {
      hydro_dump( "electron", global->hedParams );
      if ( global->mobile_ions ) {
        if ( global->I1_present ) hydro_dump( "I1", global->hI1dParams );
        if ( global->I2_present ) hydro_dump( "I2", global->hI2dParams );
      }
    }

  }

  // Particle dump data
#if 0
  if ( should_dump(particle) && global->load_particles ) {
    dump_particles( "electron", "particle/eparticle" );
    if ( global->mobile_ions ) {
      if (global->I1_present) dump_particles( "I1", "particle/I1particle" );
      if (global->I2_present) dump_particles( "I2", "particle/I2particle" );
    }
  }
#endif


  //----------------------------------------------------------------------------
  // particle tracking


  // Set TRACER_ACCUM_HYDRO to 1 if we need to accumulate hydro moments before 
  // writing trajectory output. Since this involves a pass through all the particles
  // in the system as well as synchronization (which hits MPI), don't do this step
  // unless we must. 

//??????????????????????????????????????????????????????????????????????
# undef  TRACER_DO_ACCUM_HYDRO
# define TRACER_DO_ACCUM_HYDRO 1 

  // Setup data needed for hydro output
# ifdef TRACER_DO_ACCUM_HYDRO 
  TRACER_HYDRO_SETUP( e,  "electron" )
  TRACER_HYDRO_SETUP( I2, "I2"       )
# endif 

#undef  TRACER_NUM_ADDED_FIELDS
#define TRACER_NUM_ADDED_FIELDS 10

# if 1
// FIXME: interpolate here? 
#undef TRACER_USER_DEFINED_DATA
#define TRACER_USER_DEFINED_DATA                                              \
         pout[index + TRACER_NUM_FIELDS_BASE + 0]  = e_hydro->h[p->i].rho;       \
         pout[index + TRACER_NUM_FIELDS_BASE + 1]  = e_hydro->h[p->i].jx;        \
         pout[index + TRACER_NUM_FIELDS_BASE + 2]  = e_hydro->h[p->i].px;        \
         pout[index + TRACER_NUM_FIELDS_BASE + 3]  = e_hydro->h[p->i].txx;       \
         pout[index + TRACER_NUM_FIELDS_BASE + 4]  = e_hydro->h[p->i].jy;        \
         pout[index + TRACER_NUM_FIELDS_BASE + 5]  = e_hydro->h[p->i].py;        \
         pout[index + TRACER_NUM_FIELDS_BASE + 6]  = e_hydro->h[p->i].tyy;       \
         pout[index + TRACER_NUM_FIELDS_BASE + 7]  = e_hydro->h[p->i].jz;        \
         pout[index + TRACER_NUM_FIELDS_BASE + 8]  = e_hydro->h[p->i].pz;        \
         pout[index + TRACER_NUM_FIELDS_BASE + 9]  = e_hydro->h[p->i].tzz;         
# endif 

#define CLEAR_BUFFER( BUF, LEN ) for ( int itmp=0; itmp<LEN; ++itmp ) BUF[itmp]=0; 

// For debugging 
#define DUMP_BUFFER( BUF, LEN ) for ( int itmp=0; itmp<LEN; ++itmp )          \
                                  sim_log( "DEBUG - "<<itmp<<" : "<<BUF[itmp] ); 


// The following is not used now:
// Macro for storing tracers on one-pass 
// FIXME: put this in header after it's validated 
#ifndef store_tracers
#define store_tracers(BUF)                                                    \
  BEGIN_PRIMITIVE {                                                           \
    float ex, ey, ez, bx, by, bz;                                             \
    float dx0, dy0, dz0;                                                      \
    float ux, uy, uz;                                                         \
    species_t *s = global->tracers_list ;                                     \
    const particle_t     * ALIGNED(32) p;                                     \
    const interpolator_t * ALIGNED(16) f;                                     \
    const grid_t * g = grid;                                                  \
    int ii, j, nvar;                                                          \
    float *pout;                                                              \
    int32_t local_rank[1];                                                    \
    int isp=0;                                                                \
                                                                              \
    nvar = TRACER_NUM_FIELDS_BASE + TRACER_NUM_ADDED_FIELDS;                  \
    *local_rank = rank();                                                     \
                                                                              \
    while( s ) {                                                              \
      if ( s->np > 0 ) {                                                      \
        pout = new float[nvar];                                               \
        for ( p=s->p, j=0; j<s->np; j++, p++ ) {                              \
          int index = 0;                                                      \
          dx0 = p->dx;                                                        \
          dy0 = p->dy;                                                        \
          dz0 = p->dz;                                                        \
          ii  = p->i;                                                         \
          ux  = p->ux;                                                        \
          uy  = p->uy;                                                        \
          uz  = p->uz;                                                        \
          f   = interpolator + ii;                                            \
          ex  = f->ex + dy0*f->dexdy + dz0*(f->dexdz+dy0*f->d2exdydz);        \
          ey  = f->ey + dz0*f->deydz + dx0*(f->deydx+dz0*f->d2eydzdx);        \
          ez  = f->ez + dx0*f->dezdx + dy0*(f->dezdy+dx0*f->d2ezdxdy);        \
          bx  = f->cbx + dx0*f->dcbxdx;                                       \
          by  = f->cby + dy0*f->dcbydy;                                       \
          bz  = f->cbz + dz0*f->dcbzdz;                                       \
          pout[index + 0]  = (float) p->q;                                    \
          pout[index + 1]  = (float) tracer_x ;                               \
          pout[index + 2]  = (float) tracer_y ;                               \
          pout[index + 3]  = (float) tracer_z ;                               \
          pout[index + 4]  = ux;                                              \
          pout[index + 5]  = uy;                                              \
          pout[index + 6]  = uz;                                              \
          pout[index + 7]  = ex;                                              \
          pout[index + 8]  = ey;                                              \
          pout[index + 9]  = ez;                                              \
          pout[index + 10] = bx;                                              \
          pout[index + 11] = by;                                              \
          pout[index + 12] = bz;                                              \
          pout[index + 13] = p->dx;                                           \
          pout[index + 14] = p->dy;                                           \
          pout[index + 15] = p->dz;                                           \
          pout[index + 16] = *reinterpret_cast<float*>( local_rank );         \
          pout[index + 17] = *reinterpret_cast<float*>( &ii );                \
          TRACER_USER_DEFINED_DATA;                                           \
          /* Now write pout[] into buffer array */                            \
          int iloc = ibuf * istride +                                         \
                     ( isp * global->num_tracers_total + p->q )               \
                     * (TRACER_NUM_FIELDS_BASE + TRACER_NUM_ADDED_FIELDS);    \
          for ( int ifield=0;                                                 \
                ifield < TRACER_NUM_FIELDS_BASE + TRACER_NUM_ADDED_FIELDS;    \
                ++ifield ) {                                                  \
            BUF[iloc+ifield] = pout[ifield];                                  \
          }                                                                   \
        }                                                                     \
        delete [] pout;                                                       \
      }                                                                       \
      isp++;                                                                  \
      s = s->next;                                                            \
    }                                                                         \
  } END_PRIMITIVE
#endif // dump_tracers


  if ( global->particle_tracing==1 ) { 
    static int initted=0;
    static float *tracer_buffer = NULL; 
    int istride = global->num_tracers_total 
                  * (TRACER_NUM_FIELDS_BASE + TRACER_NUM_ADDED_FIELDS)
                  * global->nspec_tracer; 
    int len     = global->nbuf_tracer * istride; 
//  this replaces old traj dir:
    char dname[] = "tracers_1pass";

    if ( initted==0 ) { 
      // allocate tracer_buffer and prepare for storing trajectories
      ALLOCATE( tracer_buffer, len, float ); 
      CLEAR_BUFFER( tracer_buffer, len ); 
         
      if ( rank()==0 && step()==0 ) { 

        // create output directory if necessary 
        dump_mkdir(dname);

        // write an ancillary file contining information about the trajectory files
        FILE *fp=NULL; 
        char fname[256];
        sprintf( fname, "%s/%s", dname, "trajectory_data.out" ); 
        fp = fopen( fname, "w" ); 
        if ( !fp ) ERROR(("Cannot open file: %s\n", fname)); 
        fprintf( fp, "nspec_tracer      %d\n", global->nspec_tracer ); 
        fprintf( fp, "num_tracers_total %d\n", global->num_tracers_total ); 
        fprintf( fp, "nvar              %d\n", TRACER_NUM_FIELDS_BASE + TRACER_NUM_ADDED_FIELDS ); 
        fclose( fp ); 

      } // if 

      mp_barrier(); /* avoid race condition */
      initted=1; 
    } // if 

    if (should_dump(tracer) ) {
      species_t *s = global->tracers_list; 
      int ip; 
      //int ibuf    = remainder( step, global->tracer_interval * global->nbuf_tracer );
      int ibuf    = ( step() / global->tracer_interval ) % global->nbuf_tracer;

      // DEBUG 
#if 0
      sim_log( "DEBUG - ibuf = "<<ibuf ); 
      sim_log( "DEBUG - len =  "<<len ); 
      sim_log( "DEBUG - istride =  "<<istride ); 
#endif

      if ( ibuf==0 ) { 
//      sim_log( "DEBUG - clearing tracer_buffer "); 
        CLEAR_BUFFER( tracer_buffer, len ); 
      } 

      if ( TRACER_DO_ACCUM_HYDRO ) { 
        // accumulate electron hydro
        TRACER_ACCUM_HYDRO( e ); 
        // accumulate I2 hydro 
        TRACER_ACCUM_HYDRO( I2 ); 
      } // if 

      // store tracers in buffer 
      //store_tracers(tracer_buffer);
      //

      // override store_tracers() macro for debugging 
      BEGIN_PRIMITIVE {
        float ex, ey, ez, bx, by, bz;
        float dx0, dy0, dz0;
        float ux, uy, uz;
        species_t *s = global->tracers_list ;
        const particle_t     * ALIGNED(32) p;
        const interpolator_t * ALIGNED(16) f;
        const grid_t * g = grid;
        int ii, j, nvar;
        float *pout;
        int32_t local_rank[1];
        int isp=0;
        
        nvar = TRACER_NUM_FIELDS_BASE + TRACER_NUM_ADDED_FIELDS;
        *local_rank = rank();
        
        while( s ) {
          if ( s->np > 0 ) {
            pout = new float[nvar];
            for ( p=s->p, j=0; j<s->np; j++, p++ ) {
              int index = 0;
              dx0 = p->dx;
              dy0 = p->dy;
              dz0 = p->dz;
              ii  = p->i;
              ux  = p->ux;
              uy  = p->uy;
              uz  = p->uz;
              f   = &interpolator(ii);
              ex  = f->ex + dy0*f->dexdy + dz0*(f->dexdz+dy0*f->d2exdydz);
              ey  = f->ey + dz0*f->deydz + dx0*(f->deydx+dz0*f->d2eydzdx);
              ez  = f->ez + dx0*f->dezdx + dy0*(f->dezdy+dx0*f->d2ezdxdy);
              bx  = f->cbx + dx0*f->dcbxdx;
              by  = f->cby + dy0*f->dcbydy;
              bz  = f->cbz + dz0*f->dcbzdz;
              pout[index + 0]  = step()*grid->dt;
              // pout[index + 0]  = (float) p->q; // DEBUG
              pout[index + 1]  = (float) tracer_x ;
              pout[index + 2]  = (float) tracer_y ;
              pout[index + 3]  = (float) tracer_z ;
              pout[index + 4]  = ux;
              pout[index + 5]  = uy;
              pout[index + 6]  = uz;
              pout[index + 7]  = ex;
              pout[index + 8]  = ey;
              pout[index + 9]  = ez;
              pout[index + 10] = bx;
              pout[index + 11] = by;
              pout[index + 12] = bz;
              pout[index + 13] = p->dx;
              pout[index + 14] = p->dy;
              pout[index + 15] = p->dz;
              pout[index + 16] = *reinterpret_cast<float*>( local_rank );
              pout[index + 17] = *reinterpret_cast<float*>( &ii );
              TRACER_USER_DEFINED_DATA;
              /* Now write pout[] into buffer array */
              int iloc = ibuf * istride + 
                         ( isp * global->num_tracers_total + p->w ) * nvar;    
              for ( int ifield=0; ifield < nvar; ++ifield ) {
                //BUF[iloc+ifield] = pout[ifield];
                tracer_buffer[iloc+ifield] = pout[ifield];
              } // for
            } // for 
            delete [] pout;
          } // if ( s->np > 0 )
          isp++;
          s = s->next;
        } // while
      } END_PRIMITIVE;

      // now write separate electron and ion files 
      if ( ibuf == global->nbuf_tracer - 1 ) { 
        species_t    *s = global->tracers_list; 
        int          isp = 0; 
        FileIO       fh; 
        FileIOStatus status;      
        char fname[256];                     

        // DEBUG
        //sim_log("DEBUG - Starting to write tracers."); 


        // tracer_buffer is initialized to zero; tracers on local rank are written to 
        // each local copy of tracer_buffer. mpi all_sum gives every rank a copy of all 
        // tracers

        // Updated by Scott V. Luedtke, XCP-6, to use much less RAM and perform
        // much better, espicially for large runs.  Now uses an MPI ruduce in
        // single precission instead of unnecessary copies and an MPI allruduce
        // in double precision.
        BEGIN_PRIMITIVE { 
          // Do an in-place reduce to the processor
          // that will actually write to disk
          if (rank()==0){
            MPI_Reduce(MPI_IN_PLACE, tracer_buffer, len, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
          }
          else{
            MPI_Reduce(tracer_buffer, tracer_buffer, len, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
          }

          // DEBUG 
//        sim_log("DEBUG - finished exchanging tracers among processors"); 
          //DUMP_BUFFER( tracer_buffer, len ); 

        } END_PRIMITIVE; 
 
        // On rank 0, write tracer_buffer to file
        if ( rank()==0 ) { 
          while ( s ) { 
 
            // compute file offset: 
            //
            // num writes  = step / ( nbuf_tracer * tracer_interval )
            // each write  = nbuf_tracer * global->num_tracers_total records
            // each record = (TRACER_NUM_FIELDS_BASE + TRACER_NUM_ADDED_FIELDS) floats
            //
            // offset is the number of bytes to advance from start of file 
            //  
            // BJA - reworked to (hopefully) prevent overflowing intermediate values
            //
            // uint64_t offset = (int)( step / ( global->nbuf_tracer * global->tracer_interval ) )
            //                        * global->nbuf_tracer 
            //                        * global->num_tracers_total 
            //                        * (TRACER_NUM_FIELDS_BASE + TRACER_NUM_ADDED_FIELDS)
            //                        * sizeof(float); 
                                    
            uint64_t offset  = step() / ( global->nbuf_tracer * global->tracer_interval );
                     offset *= global->nbuf_tracer; 
                     offset *= global->num_tracers_total;
                     offset *= (TRACER_NUM_FIELDS_BASE + TRACER_NUM_ADDED_FIELDS);
                     offset *= sizeof(float); 
 
            // open file 
            sprintf( fname, "%s/%s.out", dname , s->name );                 
            status = fh.open( fname,                                          
                              ( step() > global->nbuf_tracer*global->tracer_interval 
                                ? io_read_write 
                                : io_write ) );        
            if ( status == fail ) {                                           
              ERROR(("Could not open file %s", fname));                     
            }                                                                 
 
            // now unpack buffer and write buffered tracer data to files
            for ( int ib=0; ib<global->nbuf_tracer; ++ib ) { 

              // index within buffer from which to start the write
              int buffer_offset = ib * istride + 
                                  ( isp * global->num_tracers_total ) 
                                  * (TRACER_NUM_FIELDS_BASE + TRACER_NUM_ADDED_FIELDS);  

              // additional file offset from buffer ib
              //
              // BJA - reworked to (hopefully) prevent overflowing intermediate values
              //
              // uint64_t additional_offset = ib 
              //                              * global->num_tracers_total 
              //                              * (TRACER_NUM_FIELDS_BASE + TRACER_NUM_ADDED_FIELDS) 
              //                              * sizeof(float);
                
              uint64_t additional_offset  = ib;
                       additional_offset *= global->num_tracers_total;
                       additional_offset *= (TRACER_NUM_FIELDS_BASE + TRACER_NUM_ADDED_FIELDS); 
                       additional_offset *= sizeof(float);

              fh.seek( offset + additional_offset, SEEK_SET ); 
              fh.write( tracer_buffer + buffer_offset, 
                        global->num_tracers_total * (TRACER_NUM_FIELDS_BASE + TRACER_NUM_ADDED_FIELDS) ); 
            } 
            fh.close();                                                       
 
            isp++; 
            s = s->next; 

          } // while
        } // if ( rank()==0 )
      } // if ( ibuf == global->tracer_interval - 1 ) 
    } // if ( should_dump(tracer) ) 
  } // if ( global->particle_tracing==1 ) 









#define ALLOCATE(A,LEN,TYPE)                                             \
  if ( !((A)=(TYPE *)malloc((size_t)(LEN)*sizeof(TYPE))) ) ERROR(("Cannot allocate."));

// This is the fft diag along x
#if 0
  // ------------------------------------------------------------------------
  // Custom diagnostic where we write out Ex, Ey, and cBz for each cell in grid. 
  // These data are stored for each time step.
  //
  // Refactored for 2D LPI problem.  
# define FFT_HEADER_SIZE (sizeof(int))

# define WRITE_FFT(SUFF,INTERVAL)                                              \
  BEGIN_PRIMITIVE {                                                            \
    status = fileIO_##SUFF.open( fname_##SUFF, io_read_write);                 \
    if ( status==fail ) ERROR(("Could not open file."));                       \
    fileIO_##SUFF.seek( uint64_t( FFT_HEADER_SIZE                              \
                                  + uint64_t((step()/INTERVAL)*stride*sizeof(float)) ),  \
                        SEEK_SET );                                            \
    fileIO_##SUFF.write( SUFF, stride );                                       \
    fileIO_##SUFF.close();                                                     \
  } END_PRIMITIVE

# define WRITE_FFT_HEADER(SUFF)                                                \
  BEGIN_PRIMITIVE {                                                            \
    status = fileIO_##SUFF.open( fname_##SUFF, io_write);                      \
    if ( status==fail ) ERROR(("Could not open file."));                       \
    fileIO_##SUFF.write( &grid->nx, 1 );                                       \
    fileIO_##SUFF.close();                                                     \
  } END_PRIMITIVE


  // BJA - only write FFT on a single processor to reduce data volume

  if ( global->run_mode==1 && rank()==(int)(nproc()/2) ) {
    static int initted=0;
    static float *ex, *ey, *cbz;
    static char fname_ex[256], fname_ey[256], fname_cbz[256];
    FileIO       fileIO_ex, fileIO_ey, fileIO_cbz;
    FileIOStatus status;
    int stride=grid->nx;

    if ( !initted ) {
      // Allocate space for data arrays
      long array_length=grid->nx;
      ALLOCATE(ex,  array_length, float);
      ALLOCATE(ey,  array_length, float);
      ALLOCATE(cbz, array_length, float);
      // Define filenames
      sprintf( fname_ex,  "fft/fft_ex.%i",  (int)rank() );
      sprintf( fname_ey,  "fft/fft_ey.%i",  (int)rank() );
      sprintf( fname_cbz, "fft/fft_cbz.%i", (int)rank() );
      // On first timestep prepend a header with number of x meshpoints to each file
      if ( step()==0 ) {
        WRITE_FFT_HEADER(ex);
        WRITE_FFT_HEADER(ey);
        WRITE_FFT_HEADER(cbz);
      }
      initted=1;
    }

    // *** Ex ***
    if ( !(step()%global->fft_ex_interval) ) {
      // Store data into array ex
      for ( int i=0; i<grid->nx; ++i ) {
     // BJA - I think the indexing might not be right - try to fix
        // int k=INDEX_FORTRAN_3(i+1,1,grid->nz/2+1,0,grid->nx+1,0,grid->ny+1,0,grid->nz+1);
        // ex[i]=field[k].ex;
        ex[i]=field(i+1,1,grid->nz/2+1).ex;
      }
      // Write array to file
      sim_log("Writing FFT data for Ex fields.");
      DIAG_LOG("Startign to dump FFT Ex data.");
      WRITE_FFT(ex, global->fft_ex_interval);
      DIAG_LOG("Finished dumping FFT Ex data.");
    }

    // *** Ey, cBz ***
    if ( !(step()%global->fft_ey_interval) ) {
      // Store data into array ey, cbz
      for ( int i=0; i<grid->nx; ++i ) {
     // BJA - I think the indexing might not be right - try to fix
        // int k=INDEX_FORTRAN_3(i+1,1,grid->nz/2+1,0,grid->nx+1,0,grid->ny+1,0,grid->nz+1);
        // ey[i]  = field[k].ey /            global->emax;
        // cbz[i] = field[k].cbz/(grid->cvac*global->emax);
        ey[i]  = field(i+1,1,grid->nz/2+1).ey /            global->emax;
        cbz[i] = field(i+1,1,grid->nz/2+1).cbz/(grid->cvac*global->emax);
      }
      // Write array to file
      sim_log("Writing FFT data for Ey, cBz fields.");
      DIAG_LOG("Startign to dump FFT Ey, cBz data.");
      WRITE_FFT(ey,  global->fft_ey_interval);
      WRITE_FFT(cbz, global->fft_ey_interval);
      DIAG_LOG("Finished dumping FFT Ey, cBz data.");
    }
  }


#endif


#if 0
   //=============================================================
   // Graceful shutdown of VPIC when a signal is received from a prepsecified path.
   // Written to work around problems with Roadrunner and fitting jobs  in while
   // Bob Weaver is trying to use the machine.   Restart files are  written and
   // execution is terminated when the first character of the file  given by
   // file_flag_path is an ASCII '1'.

   // Full path for signal file
   char flag_file_path[512]       = "/users/lyin/moab_signal";

   // Set this to a reasonable value for ~5 minutes runtime between  checks, perhaps
   int query_moab_signal_interval = 100;

   // Flag to signal whether VPIC should reset the flag after it receives it.  This
   // is probably not desired unless we are only running one job.   Better would be
   // for the cron job to set it to "no action" and set '1' in the  first byte if
   // shutdown is desired.
   int have_vpic_reset_flag       = 0;

   if ( step()%query_moab_signal_interval==0 ) {
     sim_log( "Checking if signal is received for saving restart file and terminating." );
     char moab_signal[32];   // Kludge - no wrapper for formatted ASCII read :(
     FileIO fileIO;
     FileIOStatus status;

     // Read moab_signal from file
     status = fileIO.open( flag_file_path, io_read );
     if ( status==fail ) ERROR(("Could not open file."));
     fileIO.read( &moab_signal[0], 1 );
     fileIO.close();
     if ( moab_signal[0]=='1' ) {
       dump_restart("restart/restart",0);
       sim_log( "Signal received; Dump restart completed.  Terminating." );

       // Reset file to '0' as necessary
       if ( have_vpic_reset_flag ) {
         const char clear_string[] = "0\n";
         status = fileIO.open( flag_file_path, io_write );
         if ( status==fail ) ERROR(("Could not open file."));
         fileIO.write( clear_string, 2 );
         fileIO.close();
       }

       // Terminate VPIC
       mp_barrier(); // Just to be safe
       mp_finalize();
       exit(0);
     }
   }
   //=================================================================
#endif


  // Restart dump filenames are by default tagged with the current timestep.
  // If you do not want to do this add "0" to the dump_restart arguments. One
  // reason to not tag the dumps is if you want only one restart dump (the most
  // recent one) to be stored at a time to conserve on file space. Note: A
  // restart dump should be the _very_ _last_ dump made in a diagnostics
  // routine. If not, the diagnostics performed after the restart dump but
  // before the next timestep will be missed on a restart. Restart dumps are
  // in a binary format. Each rank makes a restart dump.

  // Restarting from restart files is accomplished by running the executable 
  // with "restart restart" as additional command line arguments.  The executable
  // must be identical to that used to generate the restart dumps or else 
  // the behavior may be unpredictable. 

  // Note:  restart is not currently working with custom boundary conditions
  // (such as the reflux condition) and has not been tested with emission 
  // turned on.  
  
  // if ( step()!=0 && !(step()%global->restart_interval) ) dump_restart("restart",0); 

#if 0
  if ( step() && !(step()%global->restart_interval) ) {
    if ( !global->rtoggle ) {
      DIAG_LOG("Starting to dump restart0.");
      sim_log("Starting restart0 dump");
      global->rtoggle=1;
      dump_restart("restart/restart0",0);
      DIAG_LOG("Finished dumping restart0.");
      sim_log("Restart dump restart0 completed.");
    } else {
      DIAG_LOG("Starting to dump restart1.");
      sim_log("Starting restart1 dump");
      global->rtoggle=0;
      dump_restart("restart/restart1",0);
      DIAG_LOG("Finished dumping restart1.");
      sim_log("Restart dump restart1 completed.");
    }
  }
#endif

  // Shut down simulation if wall clock time exceeds global->quota_sec.
  // Note that the mp_elapsed() is guaranteed to return the same value for all
  // processors (elapsed time on proc #0), and therefore the abort should
  // be synchronized across processors.

#if 1
  if ( global->quota_check_interval && (step()%global->quota_check_interval)==(global->quota_check_interval - global->tracer_interval) ) {
    if ( uptime() > global->quota_sec ) {
//  if ( uptime() > 11.6*3600.0 ) {

    BEGIN_TURNSTILE(NUM_TURNSTILES) {
//    checkpt( "restart/restart", step() );
      checkpt( "restart/restart", global->rtoggle );
      dump_tracer_restart(global->rtoggle);   
    } END_TURNSTILE;

      // SVL - This barrier needs to be before the sim_log, else the sim_log is
      // probably a lie.
      mp_barrier();
      if (world_rank==0) {
        // We are done checkpointing.  Leave a note saying what the newest
        // checkpoint is.
        FILE * newckp = fopen("newest_checkpoint.txt", "w");
        fprintf(newckp, "restart/restart.%i", global->rtoggle);
        fclose(newckp);
      }
      sim_log( "Restart dump invoked by quota completed." );
      sim_log( "Allowed runtime exceeded for this job.  Terminating." );
      mp_barrier(); // Just to be safe
      halt_mp();
      exit(0);
    }
  }
#endif

  if ( global->restart_interval>0 && ((step()%global->restart_interval)==(global->restart_interval - global->tracer_interval)) ) {
    //double dumpstart = uptime();

    // NOTE: If you want an accurate start time, you need an MPI barrier.
    // This can feasibly slow down the simulation, so you might not want it
    // unless you are benchmarking the dump.
    //mp_barrier();
    std::time_t dumpstart = std::time(NULL);

//???????????????????????????????????????????????????????????????????????
    BEGIN_TURNSTILE(NUM_TURNSTILES) {
//    checkpt( restart_fbase[global->rtoggle], step() );
      checkpt("restart/restart", global->rtoggle);
      dump_tracer_restart(global->rtoggle);   
 } END_TURNSTILE;

    // SVL - This is a kind of messy way to calculate this.
    //double dumpelapsed = uptime() - dumpstart;

    // Do put a barrier here so we don't say we're done before we are
    mp_barrier();

    if (world_rank==0) {
      double dumpelapsed = std::difftime(std::time(NULL), dumpstart);
      sim_log("Restart duration "<<dumpelapsed<<" seconds");

      // We are done checkpointing.  Leave a note saying what the newest
      // checkpoint is.
      FILE * newckp = fopen("newest_checkpoint.txt", "w");
      fprintf(newckp, "restart/restart.%i", global->rtoggle);
      fclose(newckp);
    }

    global->rtoggle^=1;
  }


#if 0
    {
    double elapsed_time = mp_elapsed();
    if ( elapsed_time > global->quota_sec ) {
      // elapsed_time equiv. across processors, so no danger of DIAG_LOG deadlock
      /* Writing restart file.  Make sure buffered IO writes don't rely on this */
      dump_restart("restart/restart",0);
      sim_log( "Restart dump restart completed." );
      sim_log( "Allowed runtime exceeded for this job.  Terminating....\n");
      mp_barrier(); // just to be safe
      mp_finalize();
      exit(0);
    }
  }
#endif

}


begin_particle_injection {
  // No particle injection for this simulation
static int initted=0;

advance_tracers();
}


begin_current_injection {
  // No current injection for this simulation
}

begin_field_injection { 
  // Inject a light wave from lhs boundary with E aligned along y
  // For 2d, 3d, use scalar diffraction theory for the Gaussian beam source. 
  // See, e.g., Lin's notes on this or Brian's notes of 12 March 2005. 
  // (Note the sqrt(2) discrepancy in the two regarding definition of waist. 
  // The driver below assumes Brian's definition). 

  // Note, for quiet startup (i.e., so that we don't propagate a delta-function noise pulse
  // at time t=0) we multiply by a constant phase term exp(i phi) where: 
  //   phi = k*global->xfocus+atan(h)/2    (2D)
  //   phi = k*global->xfocus+atan(h)      (3D)

# define loop_over_left_boundary \
    for ( int iz=1; iz<=grid->nz+1; ++iz ) for ( int iy=1; iy<=grid->ny; ++iy )
# define DY    ( grid->y0 + (iy-0.5)*grid->dy - global->ycenter )
# define DZ    ( grid->z0 +  (iz-1) *grid->dz - global->zcenter )
# define R2    ( DY*DY+DZ*DZ )
# define PHASE ( global->omega_0*t + h*R2/(global->width*global->width) )
# define MASK  ( R2<=pow(global->mask*global->width,2) ? 1 : 0 )


  if ( global->launch_wave == 0 ) return;

  if ( grid->x0==0 ) { // Node is on left boundary
    double alpha = grid->cvac*grid->dt/grid->dx;
//  double emax_coeff = ((4/(1+alpha))*global->omega_0*grid->dt*global->emax)/1.7; //1.7 factor is correction from test runs
    double emax_coeff = ((4/(1+alpha))*global->omega_0*grid->dt*global->emax); //no 1.7 factor but sqrt(2/pi) is used in prefactor
    double t=grid->dt*step();

#if 0
    double pulse_shape_factor=1;
    float pulse_length = 50.0*global->wpe1fs; // in units of 1/wpe, = 50 fs
    float sin_t_tau = sin(0.5*t*M_PI/pulse_length);
    pulse_shape_factor=( t<pulse_length ? sin_t_tau : 1 );
#endif

#if 0
    double pulse_shape_factor=1;
    float pulse_length = 70; // in units of 1/wpe
    float sin_t_tau = sin(0.5*t*M_PI/pulse_length);
    pulse_shape_factor=( t<pulse_length ? sin_t_tau : 1 );
#endif 

#if 1
    double pulse_shape_factor=1;
    if ( global->pulse_shape>=1 ) pulse_shape_factor=( t<global->pulse_length ? 1 : 0 );
    if ( global->pulse_shape==2 ) {
    float sin_t_tau = sin(t*M_PI/global->pulse_length);
      pulse_shape_factor =( t<global->pulse_length ? sin_t_tau : 0 );
  }
    // Hyperbolic secant profile
    if (global->pulse_shape==3){
        double fac = log(1.+sqrt(2.))/global->pulse_FWHM;
        // I think code time and code length are related by c=1
        double z = t - global->pulse_start;
        z *= fac;
        pulse_shape_factor = 2./(exp(z)+exp(-z));
    }
#endif 

    double prefactor = emax_coeff*sqrt(2.0/M_PI);               // Wave norm at edge of box
//  double prefactor = emax_coeff;                             // Wave norm at edge of box
    double rl     = M_PI*global->waist*global->waist/global->lambda;// in c/wpe, Rayleigh length
    double h = global->xfocus/rl;                 // distance/Rayleigh length

#if 0 // cir
    for ( int iz=1; iz<=grid->nz+1; ++iz )
      for ( int iy=1; iy<=grid->ny; ++iy )
        field(1,iy,iz).ey+=prefactor*cos(PHASE)*exp(-R2Z/(global->width*global->width))*MASK*pulse_shape_factor/sqrt(2.0);
      for ( int iy=1; iy<=grid->ny+1; ++iy )
        field(1,iy,iz).ez+=prefactor*sin(PHASE)*exp(-R2Z/(global->width*global->width))*MASK*pulse_shape_factor/sqrt(2.0);
#endif

#if 0 // Linear
    for ( int iz=1; iz<=grid->nz+1; ++iz )
      for ( int iy=1; iy<=grid->ny; ++iy )
        field(1,iy,iz).ey+=prefactor*cos(PHASE)*exp(-R2Z/(global->width*global->width))*MASK*pulse_shape_factor;
#endif

#if 1 // Linear
    for ( int iz=1; iz<=grid->nz+1; ++iz )
      for ( int iy=1; iy<=grid->ny; ++iy )
        field(1,iy,iz).ey+=prefactor*cos(PHASE)*exp(-R2/(global->width*global->width))*MASK*pulse_shape_factor;
#endif

#if 0 // Linear plane wave
    for ( int iz=1; iz<=grid->nz+1; ++iz )
      for ( int iy=1; iy<=grid->ny; ++iy )
        field(1,iy,iz).ey+=prefactor*cos(PHASE)*pulse_shape_factor;
# endif

  }

} 

begin_particle_collisions {
  // No particle collisions for this simulation
}

