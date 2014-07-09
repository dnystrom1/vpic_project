//========================================================================
//
// LPI 3D deck - Linearly polarized (in y) plane wave incident from left 
//               boundary 
//
// Adapted from Albright's Lightning 3D LPI deck.  
// B. Albright, X-1-PTA; 28 Jan. 2007
// Lin Yin, X-1-PTA, 23 Feb 2009, for Cerrillos test 
// 
// Executable creates its own directory structure.  Remove the old with: 
//
// rm -rf rundata ehydro Hhydro Hehydro restart poynting velocity particle field
//========================================================================

// Employ turnstiles to partially serialize the high-volume file writes. 
// In this case, the restart dumps.  Set NUM_TURNSTILES to be the desired
// number of simultaneous writes. 
#define NUM_TURNSTILES 256 

begin_globals {
  double e0;                   // peak amplitude of oscillating electric field
  double omega;                // angular freq. of the beam
  int field_interval;          // how often to dump field and hydro
  int particle_interval;       // how often to dump particle data
  int poynting_interval;       // how often to compute poynting flux on boundary
  int restart_interval;        // how often to write restart data
  int quota_check_interval;    // how often to check if quota exceeded
  double quota_sec;            // run quota in sec
  int rtoggle;                 // enable save of last 2 restart files for safety
  int load_particles;          // were particles loaded? 
  double topology_x;           // domain topology needed to normalize Poynting diagnostic 
  double topology_y;
  double topology_z;
  int mobile_ions;
  int H_present; 
  int He_present; 

  // Parameters for 3d Gaussian wave launch
  double lambda;               
  double waist;                // how wide the focused beam is
  double width;                
  double zcenter;              // center of beam at boundary in z 
  double ycenter;              // center of beam at boundary in y
  double xfocus;               // how far from boundary to focus
  double mask;                 // # gaussian widths from beam center where I is nonzero

  // Ponyting diagnostic flags - which output to turn on

  // write_backscatter_only flag:  when this flag is nonzero, it means to only compute 
  // poynting data for the lower-x surface.  This flag affects both the summed poynting 
  // data as well as the surface data. 

  int write_poynting_data;     // Whether to write poynting data to file (or just stdout)

  int write_backscatter_only;  // nonzero means we only write backscatter diagnostic for fields

  // write_poynting_sum is nonzero if you wish to write a file containing the integrated
  // poynting data on one or more surfaces.  If write_backscatter_only is nonzero, then 
  // the output file will be a time series of double precision numbers containing the 
  // integrated poynting flux at that point in time.  If write_backscatter_only is zero,
  // then for each time step, six double precision numbers are written, containing the 
  // poynting flux through the lower x, upper x, lower y, upper y, lower z, and upper z
  // surfaces. 

  int write_poynting_sum;      // nonzero if we wish to write integrated Poynting data

  // write_poynting_faces is probably useless, but I put it here anyway in case it's not. 
  // When this flag is nonzero, it will print out the poynting flux at each of a 2D array 
  // of points on the surface.  When this flag is turned on and write_backscatter_only is 
  // nonzero, then only the 2D array of points on the lower-x boundary surface are written
  // for each time step.  When this flag is turned on and write_bacscatter_only is 
  // zero, then it will write 2D array data for the lower x, upper x, lower y, upper y, 
  // lower z, upper z surfaces for each time step. 

  int write_poynting_faces;    // nonzero if we wish to write Poynting data on sim boundaries 

  // write_eb_faces is nonzero when we wish to get the raw e and b data at the boundary
  // (e.g., to be used with a filter to distinguish between SRS and SBS backscatter).  
  // When this flag is on and write_backscatter_only is nonzero, then only the 2D array
  // of points on the lower-x boundary surface are written for each time step.  When 
  // this flag is turned on and write_backscatteR_only is zero, then it will write 2D
  // array data for the lower x, upper x, lower y, upper y, lower z, upper z surfaces for
  // each time step.  When turned on, four files are produced: e1, e2, cb1, cb2.  The 
  // values of the quantities printed depend on the face one is considering:  for the 
  // x faces, e1 = ey, e2 = ez, cb1 = cby, cb2 = cbz.  Similarly for y and z surfaces, 
  // but where 1 and 2 stand for (z,x) and (x,y) coordinates, respectively.  

  int write_eb_faces;          // nonzero if we wish to write E and B data on sim boundaries

  // Needed for movie files
  float vthe;                   // vthe/c  
  float vthi_He;                // vthi_He/c
  float vthi_H;                 // vthi_H/c
  int   velocity_interval;      // how frequently to dump binned velocity space data

  // Dump parameters for standard VPIC output formats
  DumpParameters fdParams;
  DumpParameters hedParams;
  DumpParameters hHdParams;
  DumpParameters hHedParams;
  std::vector<DumpParameters *> outputParams;
};

begin_initialization {
  
  // System of units
  double ec         = 4.8032e-10;          // stat coulomb
  double c_vac      = 2.99792458e10;       // cm/sec
  double m_e        = 9.1094e-28;          // g
  double k_b        = 1.6022e-12;          // erg/eV
  double mec2       = m_e*c_vac*c_vac/k_b; 
  double mpc2       = mec2*1836.0; 

  double cfl_req    = 0.98;                // How close to Courant should we try to run
  double damp       = 0;                   // How much radiation damping
  double iv_thick   = 2;                   // Thickness of impermeable vacuum (in cells)

  // Experimental parameters

  double t_e               = 2600;         // electron temperature, eV
  double t_i               = 1300;         // ion temperature, eV
  double n_e_over_n_crit   = 0.1415;         // n_e/n_crit
  double vacuum_wavelength = 527 * 1e-7;   // third micron light (cm)
  double laser_intensity   = 6.0e15 * 1e7;   // in ergs/cm^2 (note: 1 W = 1e7 ergs)

  // Simulation parameters 

#if 0
  // 2048 processors
  double Lx                = 4*35.0* 1e-4;   // In cm (note: 1 micron = 1e-4 cm)   
  double Ly                = 2*6.0 * 1e-4;              
  double Lz                = 2*6.0 * 1e-4;                 
  double nx                = 8192;
  double ny                = 512;
  double nz                = 512; 
  double topology_x        = 64;
  double topology_y        = 4;
  double topology_z        = 8;            
  // single-processor mesh = 128 x 128 x 64
#endif

#if 1
  // 
  double Lx                = 16 * 12.0 * 1e-4;   // In cm (note: 1 micron = 1e-4 cm)   
  double Ly                =      12.0 * 1e-4;              
  double Lz                =      12.0 * 1e-4;                 
//???????????????????????????????????????????????????????
  double nx                = 44*16;
  double ny                = 44;
  double nz                = 44; 
  double topology_x        = 16;
  double topology_y        = 1;
  double topology_z        = 1;            
  // single-processor mesh = 16 * 44 x 44 x 44
#endif

  double nppc               = 150;     // Ave. number of particles/cell in ea. species
  int load_particles        = 1;       // Flag to turn on/off particle load 
  int mobile_ions           = 0;       // Whether or not to push ions
  double f_He               = 0;       // Ratio of number density of He to total ion density
  double f_H                = 1-f_He;  // Ratio of number density of H  to total ion density
  int H_present             = ( (f_H !=0) ? 1 : 0 );  
  int He_present            = ( (f_He!=0) ? 1 : 0 );  

  // Precompute some useful variables. 
  double A_H                = 1;
  double A_He               = 4;
  double Z_H                = 1;
  double Z_He               = 2; 
  double mic2_H             = mpc2*A_H;
  double mic2_He            = mpc2*A_He;
  double mime_H             = mic2_H /mec2; 
  double mime_He            = mic2_He/mec2; 
                            
  double uth_e              = sqrt(t_e/mec2);     // vthe/c
  double uthi_H             = sqrt(t_i/mic2_H);   // vthi/c for H
  double uthi_He            = sqrt(t_i/mic2_He);  // vthi/c for He

  // Plasma skin deptth in cm
  double delta = (vacuum_wavelength / (2*M_PI) ) / sqrt( n_e_over_n_crit ); 

  double n_e   = c_vac*c_vac*m_e/(4*M_PI*ec*ec*delta*delta); // electron density in cm^-3
  double debye = uth_e*delta;                     // electron Debye length (cm)
  double omega = sqrt( 1/n_e_over_n_crit );       // laser beam freq. in wpe

  // Peak instantaneous E field in "natural units" 
  double e0    = sqrt( 2*laser_intensity / (m_e*c_vac*c_vac*c_vac*n_e) );  

  // Set up local mesh resolution and time step
  Lx /= delta;                                    // Convert box size to skin depths
  Ly /= delta;   
  Lz /= delta;   

  double hx = Lx/nx; 
  double hy = Ly/ny; 
  double hz = Lz/nz; 

  double cell_size_x       = hx*delta/debye;      // Cell size in Debye lengths
  double cell_size_y       = hy*delta/debye;
  double cell_size_z       = hz*delta/debye;

  double f_number          = 6;                   // f/# of beam
  double lambda            = vacuum_wavelength/delta; // vacuum wavelength in c/wpe
  double waist             = f_number*lambda;     // width of beam at focus in c/wpe
  double xfocus            = Lx/2;                // in c/wpe
  double ycenter           = 0;                   // center of spot in y on lhs boundary
  double zcenter           = 0;                   // center of spot in z on lhs boundary
  double mask              = 1.5;                 // set drive I=0 outside r>mask*width at lhs boundary
  double width = waist*sqrt( 1 + (lambda*xfocus/(M_PI*waist*waist))*(lambda*xfocus/(M_PI*waist*waist))); 
  e0                       = e0*(waist/width);    // at entrance (3D Gaussian) 
 
  double dt                = cfl_req*courant_length(Lx,Ly,Lz,nx,ny,nz); // in 1/wpe; n.b. c=1 in nat. units
  double nsteps_cycle      = trunc_granular(2*M_PI/(dt*omega),1)+1; 
  dt                       = 2*M_PI/omega/nsteps_cycle; // nsteps_cycle time steps in one laser cycle

  //double t_stop            = 11;               // Runtime in 1/wpe
  double t_stop            = 101;               // Runtime in 1/wpe
  //double t_stop            = 600001;               // Runtime in 1/wpe
  int particle_interval    = 0; 
  //int poynting_interval    = int(M_PI/((omega+1.5)*dt));     // Num. steps between dumping poynting flux
  int poynting_interval    = 0;                              // Num. steps between dumping poynting flux
//????????????????????????????????????????????????????????
  //int field_interval       = int(100.0/dt);                  // Num. steps between saving field, hydro data
  int field_interval       = 0;                              // Num. steps between saving field, hydro data
  int velocity_interval    = int(100.0/dt);                  // How frequently to dump binned velocity space data
  int restart_interval     = int(150.0/dt);                  // Num. steps between restart dumps
  int quota_check_interval = 20;

  // Ben:  This quota thing gracefully terminates after writing a final restart after 
  // 11.5 hours; if you want a longer time before shutdown, set this value larger.  If 
  // you want the code to just run all weekend, then set it to very long (2400.*3500, e.g.) 

//????????????????????????????????????????????????????????
  double quota_sec         = 11.6*3600;           // Run quota in sec. 

  // Turn on integrated backscatter poynting diagnostic - right now there is a bug in this, so we 
  // only write the integrated backscatter time history on the left face. 

  int write_poynting_data = 1;                    // Whether to write poynting data to file (or just stdout)

  int write_backscatter_only = 0;                 // Nonzero means only write lower x face
  int write_poynting_sum   = 0;                   // Whether to write integrated Poynting data 
  int write_poynting_faces = 0;                   // Whether to write poynting data on sim boundary faces
  int write_eb_faces       = 0;                   // Whether to write e and b field data on sim boundary faces

  double N_e               = nppc*nx*ny*nz;       // Number of macro electrons in box
  double Np_e              = Lx*Ly*Lz;            // "Number" of "physical" electrons in box (nat. units)
  double q_e               = -Np_e/N_e;           // Charge per macro electron
  double N_i               = N_e;                 // Number of macro ions of each species in box 
  double Np_i              = Np_e/(Z_H*f_H+Z_He*f_He); // "Number" of "physical" ions of each sp. in box
  double qi_H              = Z_H *f_H *Np_i/N_i;  // Charge per H  macro ion
  double qi_He             = Z_He*f_He*Np_i/N_i;  // Charge per He macro ion 
  
  // Print simulation parameters

  sim_log("***** Simulation parameters *****");
  sim_log("* Processors:                     "<<nproc());
  sim_log("* Topology:                       "<<topology_x<<" "<<topology_y<<" "<<topology_z); 
  sim_log("* nsteps_cycle =                 "<<nsteps_cycle);
  sim_log("* Time step, max time, nsteps:    "<<dt<<" "<<t_stop<<" "<<int(t_stop/(dt))); 
  sim_log("* Debye length, XYZ cell sizes:   "<<debye<<" "<<cell_size_x<<" "<<cell_size_y<<" "<<cell_size_z);
  sim_log("* Real cell sizes (in Debyes):    "<<hx/uth_e<<" "<<hy/uth_e<<" "<<hz/uth_e);
  sim_log("* Lx, Ly, Lz =                    "<<Lx<<" "<<Ly<<" "<<Lz);
  sim_log("* nx, ny, nz =                    "<<nx<<" "<<ny<<" "<<nz);
  sim_log("* Charge/macro electron =         "<<q_e);
  sim_log("* Average particles/processor:    "<<N_e/nproc());
  sim_log("* Average particles/cell:         "<<nppc);
  sim_log("* Omega_0, Omega_pe:              "<<omega<<" "<<1);
  sim_log("* Plasma density, ne/nc:          "<<n_e<<" "<<n_e_over_n_crit);
  sim_log("* Vac wavelength (nm):            "<<vacuum_wavelength*1e7);
  sim_log("* I_laser (W/cm^2):               "<<laser_intensity/1e7);
  sim_log("* T_e, T_i (eV)                   "<<t_e<<" "<<t_i);
  sim_log("* m_e, m_H, m_He                  "<<"1 "<<mime_H<<" "<<mime_He);
  sim_log("* Radiation damping:              "<<damp);
  sim_log("* Fraction of courant limit:      "<<cfl_req);
  sim_log("* vthe/c:                         "<<uth_e);
  sim_log("* vthi_H /c:                      "<<uthi_H);
  sim_log("* vthe_He/c:                      "<<uthi_He);
  sim_log("* emax at entrance:               "<<e0);
  sim_log("* emax at waist:                  "<<e0/(waist/width));
  sim_log("* Poynting interval:              "<<poynting_interval); 
  sim_log("* field interval:                 "<<field_interval); 
  sim_log("* restart interval:               "<<restart_interval); 
  sim_log("* num vacuum edge grids:          "<<iv_thick);
  sim_log("* width, waist, xfocus:           "<<width<<" "<<waist<<" "<<xfocus); 
  sim_log("* ycenter, zcenter, mask:         "<<ycenter<<" "<<zcenter<<" "<<mask); 
  sim_log("* field_interval                  "<<field_interval); 
  sim_log("* velocity_interval               "<<velocity_interval); 
  sim_log("* restart_interval:               "<<restart_interval); 
  sim_log("* quota_check_interval:           "<<quota_check_interval); 
  sim_log("* write_poynting_sum:             "<<(write_poynting_sum ? "Yes" : "No")); 
  sim_log("* write_poynting_faces:           "<<(write_poynting_faces? "Yes" : "No")); 
  sim_log("* write_eb_faces:                 "<<(write_eb_faces ? "Yes" : "No")); 
  sim_log("* write_backscatter_only:         "<<(write_backscatter_only ? "Yes" : "No")); 
  sim_log("*********************************");

  // Set up high level simulation parameters

  sim_log("Setting up high-level simulation parameters."); 
  num_step             = int(t_stop/(dt)); 

//??????????????????????????????????????????????????????????
  //status_interval      = 200; 
  status_interval      = 20; 
  sync_shared_interval = status_interval/1;
  clean_div_e_interval = status_interval/1;
  clean_div_b_interval = status_interval/10; 
  
  // Turn off some of the spam
  verbose = 1; 

  // For maxwellian reinjectoin, we need more than the default number of 
  // passes (3) through the boundary handler
  // Note:  We have to adjust sort intervals for maximum performance on Cell. 
  num_comm_round = 6;

  global->e0                   = e0; 
  global->omega                = omega; 
  global->field_interval       = field_interval; 
  global->particle_interval    = particle_interval;
  global->poynting_interval    = poynting_interval;
  global->restart_interval     = restart_interval;
  global->quota_check_interval = quota_check_interval;
  global->quota_sec            = quota_sec;
  global->rtoggle              = 0;
  global->load_particles       = load_particles;
  global->mobile_ions          = mobile_ions; 
  global->H_present            = H_present; 
  global->He_present           = He_present; 
  global->topology_x           = topology_x;  
  global->topology_y           = topology_y;  
  global->topology_z           = topology_z;  
  global->xfocus               = xfocus;  
  global->ycenter              = ycenter; 
  global->zcenter              = zcenter; 
  global->mask                 = mask; 
  global->waist                = waist; 
  global->width                = width; 
  global->lambda               = lambda; 

  global->write_poynting_data  = write_poynting_data;

  global->write_poynting_sum   = write_poynting_sum; 
  global->write_poynting_faces = write_poynting_faces; 
  global->write_eb_faces       = write_eb_faces;      
  global->write_backscatter_only = write_backscatter_only; 

  global->vthe                 = uth_e; 
  global->vthi_H               = uthi_H; 
  global->vthi_He              = uthi_He; 
  global->velocity_interval    = velocity_interval; 

  // Set up the species
  // Allow additional local particles in case of non-uniformity.

  sim_log("Setting up species."); 
  double max_local_np = 1.3*N_e/nproc(); 
//??????????????????????????????????????????????????????????????????????????????????????
  double max_local_nm = max_local_np/10.0;
  species_t * electron = define_species("electron", -1, max_local_np, max_local_nm, 15, 1); 

  // Start with two ion species.  We have option to go to Xe and Kr gas fills if 
  // we need a higher ion/electron macroparticle ratio.  

  species_t *ion_H, *ion_He; 
  if ( mobile_ions ) {
    if ( H_present  ) ion_H  = define_species("H",  Z_H /mime_H,  max_local_np, max_local_nm, 100, 1); 
    if ( He_present ) ion_He = define_species("He", Z_He/mime_He, max_local_np, max_local_nm, 100, 1); 
  }

  // Set up grid
  sim_log("Setting up computational grid."); 
  grid->dx       = hx; 
  grid->dy       = hy; 
  grid->dz       = hz; 
  grid->dt       = dt; 
  grid->cvac     = 1;
  grid->eps0     = 1; 
  grid->damp     = damp; 

  sim_log("Setting up absorbing mesh."); 
  define_absorbing_grid( 0,         -0.5*Ly,    -0.5*Lz,        // Low corner
                         Lx,         0.5*Ly,     0.5*Lz,        // High corner 
                         nx,         ny,         nz,            // Resolution
                         topology_x, topology_y, topology_z,    // Topology
                         reflect_particles );                   // Default particle boundary condition 

  // From grid/partition.c: used to determine which domains are on edge
# define RANK_TO_INDEX(rank,ix,iy,iz) BEGIN_PRIMITIVE {                   \
    int _ix, _iy, _iz;                                                    \
    _ix  = (rank);                        /* ix = ix+gpx*( iy+gpy*iz ) */ \
    _iy  = _ix/int(global->topology_x);   /* iy = iy+gpy*iz */            \
    _ix -= _iy*int(global->topology_x);   /* ix = ix */                   \
    _iz  = _iy/int(global->topology_y);   /* iz = iz */                   \
    _iy -= _iz*int(global->topology_y);   /* iy = iy */                   \
    (ix) = _ix;                                                           \
    (iy) = _iy;                                                           \
    (iz) = _iz;                                                           \
  } END_PRIMITIVE 

  sim_log("Overriding x boundaries to absorb fields."); 
  int ix, iy, iz;        // Domain location in mesh
  RANK_TO_INDEX( int(rank()), ix, iy, iz ); 

  // Set up Maxwellian reinjection B.C. 

  sim_log("Setting up Maxwellian reinjection boundary condition."); 

  maxwellian_reflux_t boundary_data[1];
  boundary_data->ut_perp[electron->id] = uth_e;     // vth_e in perp direction
  boundary_data->ut_para[electron->id] = uth_e;     // vth_e in para direction

  // Ions
  if ( mobile_ions ) { 
    if ( H_present ) { 
      boundary_data->ut_perp[ion_H->id] = uthi_H; 
      boundary_data->ut_para[ion_H->id] = uthi_H; 
    } 
    if ( He_present ) { 
      boundary_data->ut_perp[ion_He->id] = uthi_He; 
      boundary_data->ut_para[ion_He->id] = uthi_He; 
    } 
  }
  int maxwellian_reinjection = add_boundary( grid, maxwellian_reflux, boundary_data ); 
 
  if(maxwellian_reinjection == INVALID_BOUNDARY) {
  	ERROR(("Invalid Boundary Encountered"));
  } // if

  // Set up materials
  sim_log("Setting up materials."); 
  define_material( "vacuum", 1 );

  define_material( "impermeable_vacuum", 1 ); 

  finalize_field_advance( vacuum_field_advance); 
 
  // Paint the simulation volume with materials and boundary conditions
# define iv_region (   x<      hx*iv_thick || x>Lx  -hx*iv_thick  \
                    || y<-Ly/2+hy*iv_thick || y>Ly/2-hy*iv_thick  \
                    || z<-Lz/2+hz*iv_thick || z>Lz/2-hz*iv_thick ) /* all boundaries are i.v. */ 
  set_region_material( iv_region, "impermeable_vacuum",     "impermeable_vacuum" ); 
  set_region_bc( iv_region, maxwellian_reinjection, maxwellian_reinjection, maxwellian_reinjection ); 

  // Load particles 
  if ( load_particles ) {
    sim_log("Loading particles.");
    // Fast load of particles--don't bother fixing artificial domain correlations
    double xmin=grid->x0, xmax=grid->x1;
    double ymin=grid->y0, ymax=grid->y1;
    double zmin=grid->z0, zmax=grid->z1;
    repeat( N_e/(topology_x*topology_y*topology_z) ) {
      double x = uniform_rand( xmin, xmax );
      double y = uniform_rand( ymin, ymax );
      double z = uniform_rand( zmin, zmax );
      if ( iv_region ) continue;           // Particle fell in iv_region.  Don't load.
      inject_particle( electron, x, y, z,
                       maxwellian_rand(uth_e),
                       maxwellian_rand(uth_e),
                       maxwellian_rand(uth_e), q_e, 0, 0 );
      if ( mobile_ions ) {
        if ( H_present )  // Inject an H macroion on top of macroelectron
          inject_particle( ion_H, x, y, z, 
                           maxwellian_rand(uthi_H), 
                           maxwellian_rand(uthi_H), 
                           maxwellian_rand(uthi_H), qi_H, 0, 0 ); 
        if ( He_present ) // Inject an H macroion on top of macroelectron
          inject_particle( ion_He, x, y, z, 
                           maxwellian_rand(uthi_He), 
                           maxwellian_rand(uthi_He), 
                           maxwellian_rand(uthi_He), qi_He, 0, 0 ); 
      }
    }
  }


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

  global->hHdParams.format = band;
  sim_log ( "Hydrogen hydro output format : band" );

  global->hHedParams.format = band;
  sim_log ( "Helium hydro output format   : band" );

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
  // ????????????????????????????????????????????????????????????????????????
  int stride_x = 2, stride_y = 2, stride_z = 2; 
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
  // Hydrogen hydro

  // relative path to electron species data from global header
  sprintf(global->hHdParams.baseDir, "Hhydro");

  // base file name for fields output
  sprintf(global->hHdParams.baseFileName, "H_hydro");

  // set hydrogen hydro strides
  global->hHdParams.stride_x = stride_x;
  global->hHdParams.stride_y = stride_y;
  global->hHdParams.stride_z = stride_z;
  sim_log ( "Ion species x-stride " << global->hHdParams.stride_x );
  sim_log ( "Ion species y-stride " << global->hHdParams.stride_y );
  sim_log ( "Ion species z-stride " << global->hHdParams.stride_z );

  // add hydrogen hydro parameters to list
  global->outputParams.push_back(&global->hHdParams);

  //----------------------------------------------------------------------
  // Helium hydro

  // relative path to electron species data from global header
  sprintf(global->hHedParams.baseDir, "Hehydro");

  // base file name for fields output
  sprintf(global->hHedParams.baseFileName, "He_hydro");

  // set helium hydro strides
  global->hHedParams.stride_x = stride_x;
  global->hHedParams.stride_y = stride_y;
  global->hHedParams.stride_z = stride_z;
  sim_log ( "Ion species x-stride " << global->hHedParams.stride_x );
  sim_log ( "Ion species y-stride " << global->hHedParams.stride_y );
  sim_log ( "Ion species z-stride " << global->hHedParams.stride_z );

  // add helium hydro parameters to list
  global->outputParams.push_back(&global->hHedParams);

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
  global->hHdParams.output_variables(  current_density  | charge_density |
                                       momentum_density | ke_density |
                                       stress_tensor );
  global->hHedParams.output_variables( current_density  | charge_density |
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

  create_hydro_list(varlist, global->hHdParams);
  sim_log ( "Ion species variable list: " << varlist );

 /*------------------------------------------------------------------------*/

  sim_log("***Finished with user-specified initialization ***"); 

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
  if ( step%200==0 ) sim_log("Time step: "<<step); 

# define should_dump(x) \
  (global->x##_interval>0 && remainder(step,global->x##_interval)==0)

  // Do a mkdir at time t=0 to ensure we have all the directories we need
  if ( step==0 ) {
    dump_mkdir("rundata"); 
    dump_mkdir("field"); 
    dump_mkdir("ehydro"); 
    dump_mkdir("Hhydro"); 
    dump_mkdir("Hehydro"); 
    dump_mkdir("restart"); 
    dump_mkdir("particle"); 
    dump_mkdir("poynting"); 
    dump_mkdir("velocity"); 

    // Turn off rundata for now
    // dump_grid("rundata/grid");
    // dump_materials("rundata/materials");
    // dump_species("rundata/species");
    global_header("global", global->outputParams);

  }

  // Field and hydro data

  if ( should_dump(field) ) {
    field_dump( global->fdParams ); 

//?????????????????????????????????????????????????????????????????????????
#if 0  
    if ( global->load_particles ) {
      hydro_dump( "electron", global->hedParams );
      if ( global->mobile_ions ) {
        if ( global->H_present  ) hydro_dump( "H",  global->hHdParams );
        if ( global->He_present ) hydro_dump( "He", global->hHedParams );
      }
    }
#endif

  }

  // Particle dump data

  if ( should_dump(particle) && global->load_particles ) {
    dump_particles( "electron", "particle/eparticle" );
    if ( global->mobile_ions ) {
      if ( global->H_present  ) dump_particles( "H",  "particle/Hparticle" );
      if ( global->He_present ) dump_particles( "He", "particle/Heparticle" );
    }
  } 

//?????????????????????????????????????????????????????????????????????????
#if 1
  //------------------------------------------------------------------------------------
  // Old Poynting diagnostic (from GB run); only computes an integrated Poynting flux on
  // lower-x boundary. 
  //
  // Ponyting data - send to stdout for GB run
  // Write Poynting flux at left boundary
  // Poynting flux is defined positive if directed in the +x direction. 
  // TODO: It is assumed that Ponyting dumps are infrequent, so we can afford the 
  // mpi_allreduce() here.  This needs to be verified on RR hardware. 

  // If we use this stub, put the following in the globals block above: 
  // int write_poynting_data;     // Whether to write poynting data to file (or just stdout)
  
#define ALLOCATE(A,LEN,TYPE)                                             \
  if ( !((A)=(TYPE *)malloc((size_t)(LEN)*sizeof(TYPE))) ) ERROR(("Cannot allocate.")); 

  static float *pvec=NULL; 
  static double psum, gpsum; 
  FileIO       fileIO; 
  FileIOStatus status; 
  static char fname_poynting[]="poynting/poynting";
  static int stride, initted=0; 
  
  if ( !initted ) { 
    stride=(grid->ny-1)*(grid->nz-1); 
    ALLOCATE( pvec, stride, float ); 
    initted=1; 
  }

  if ( step>0 && should_dump(poynting) ) {
    int i, j, k, k1, k2, ix, iy, iz; 
    for ( i=0; i<stride; ++i ) pvec[i]=0;     // Initialize pvec to zero. 
    RANK_TO_INDEX( int(rank()), ix, iy, iz ); 
    if ( ix==0 ) {                            // Compute Poynting for domains on left of box
      for ( j=1; j<grid->ny; ++j ) {
        for ( k=1; k<grid->nz; ++k ) {
          k1 = INDEX_FORTRAN_3(1,j+1,k+1,0,grid->nx+1,0,grid->ny+1,0,grid->nz+1);
          k2 = INDEX_FORTRAN_3(2,j+1,k+1,0,grid->nx+1,0,grid->ny+1,0,grid->nz+1);
          pvec[(j-1)*(grid->nz-1)+k-1] = (  field[k2].ey*0.5*(field[k1].cbz+field[k2].cbz)
                                          - field[k2].ez*0.5*(field[k1].cby+field[k2].cby) )
                                          / (grid->cvac*grid->cvac*global->e0*global->e0);
        }
      }
    }                                         // Leave pvec = zero in mp_allsum_d for interior

    // Sum poynting flux on surface
    for ( i=0, psum=0; i<stride; ++i ) psum+=pvec[i]; 
    // Sum over all surfaces
    mp_allsum_d(&psum, &gpsum, 1, grid->mp);  
    // Divide by number of mesh points summed over
    gpsum /= stride*global->topology_y*global->topology_z; 
    
    if ( rank()==0 && global->write_poynting_data ) { 
      sim_log("Starting to write Poynting data."); 
      status=fileIO.open( fname_poynting, (step==global->poynting_interval ? io_write : io_read_write) );
      if ( status==fail ) ERROR(("Could not open file."));
      fileIO.seek( uint64_t((step/global->poynting_interval-1)*sizeof(double)), SEEK_SET );
      fileIO.write( &gpsum, 1 );
      fileIO.close(); 
      sim_log("Finished writing Poynting data."); 
    } 
    sim_log("** step = "<<step<<" Poynting = "<<gpsum);  // Dump data to stdout
  }

#endif 

#if 0
  //----------------------------------------------------------------------------
  // New Poynting diagnostic.  Lin needs the raw E and B fields at the boundary
  // in order to perform digital filtering to extract the SRS component from the 
  // SBS.  We use random access binary writes with a stride of length:
  // 
  // stride =   2* grid->ny*grid->nz * global->topology_y*global->topology_z  // lower, upper x faces
  //          + 2* grid->nz*grid->nx * global->topology_z*global->topology_x  // lower, upper y faces
  //          + 2* grid->nx*grid->ny * global->topology_x*global->topology_y; // lower, upper z faces
  // 
  // On the x faces, e1 = ey, e2 = ez, cb1 = cby, cb2 = cbz
  // On the y faces, e1 = ez, e2 = ex, cb1 = cbz, cb2 = cbx
  // On the z faces, e1 = ex, e2 = ey, cb1 = cbx, cb2 = cby
  //
  // We also write 6-element arrays of integrated poynting flux on each face:
  // 
  // vals = {lower x, upper x, lower y, upper y, lower z, upper z}
  // 
  // Note:  This diagnostic assumes uniform domains.
  // 
  // Also note:  Poynting flux in a given direction is defined as the projection
  // of E x B along the unit vector in that direction.  
  //---------------------------------------------------------------------------- 
  
#define ALLOCATE(A,LEN,TYPE)                                             \
  if ( !((A)=(TYPE *)malloc((size_t)(LEN)*sizeof(TYPE))) ) ERROR(("Cannot allocate.")); 

  BEGIN_PRIMITIVE { 
    static double *pvec =NULL, *e1vec =NULL, *e2vec =NULL, *cb1vec =NULL, *cb2vec =NULL; 
    static double *gpvec=NULL, *ge1vec=NULL, *ge2vec=NULL, *gcb1vec=NULL, *gcb2vec=NULL; 
    static double *psum, *gpsum, norm; 
    FileIO       fileIO; 
    FileIOStatus status; 
    static char fname_poynting_sum[]="poynting/poynting_sum",
                fname_poynting[]    ="poynting/poynting"    ,
                fname_e1[]          ="poynting/e1"          ,
                fname_e2[]          ="poynting/e2"          ,
                fname_cb1[]         ="poynting/cb1"         ,
                fname_cb2[]         ="poynting/cb2"         ;
    static uint64_t stride;  // Force tmp variable in seek() to be uint64_t and not int!
    static int sum_stride, initted=0; 
    int ncells_yz = int(grid->ny*grid->nz*global->topology_y*global->topology_z);
    int ncells_zx = int(grid->nz*grid->nx*global->topology_z*global->topology_x);
    int ncells_xy = int(grid->nx*grid->ny*global->topology_x*global->topology_y);

   
    if ( !initted ) { 
      if ( global->write_backscatter_only ) {
        stride     = uint64_t(ncells_yz); // x faces
        sum_stride = 1; 
      } else {
        stride     = uint64_t(  2*( ncells_yz + ncells_zx + ncells_xy ) ); 
        sum_stride = 6;
      }
      ALLOCATE( psum, sum_stride, double ); ALLOCATE( gpsum, sum_stride, double ); 
      ALLOCATE( pvec,   stride, double ); ALLOCATE( gpvec,   stride, double ); 
      ALLOCATE( e1vec,  stride, double ); ALLOCATE( ge1vec,  stride, double ); 
      ALLOCATE( e2vec,  stride, double ); ALLOCATE( ge2vec,  stride, double ); 
      ALLOCATE( cb1vec, stride, double ); ALLOCATE( gcb1vec, stride, double ); 
      ALLOCATE( cb2vec, stride, double ); ALLOCATE( gcb2vec, stride, double ); 
      norm = 1.0 / (grid->cvac*grid->cvac*global->e0*global->e0);
      initted=1; 
    }
 
    // FIXME:  On real 3D problem, measure overhead from MPI collectives
    // FIXME:  Rewrite using permutation-symmetric macros
    // FIXME:  Don't we have to do something special for mp on Roadrunner? 


    // Note:  Will dump core if we dump poynting by mistake on time t=0  

    if ( step>0 && should_dump(poynting) ) {
      uint64_t ii;  // To shut the compiler up.
      int i, j, k, k1, k2, ix, iy, iz, skip, index; 

      // Initialize arrays to zero
      for ( ii=0; ii<stride; ++ii ) { 
        pvec[ii]    = 0;  
        e1vec[ii]   = 0; 
        e2vec[ii]   = 0; 
        cb1vec[ii]  = 0;
        cb2vec[ii]  = 0;
        gpvec[ii]   = 0; 
        ge1vec[ii]  = 0;
        ge2vec[ii]  = 0;
        gcb1vec[ii] = 0;
        gcb2vec[ii] = 0;
      }
      RANK_TO_INDEX( int(rank()), ix, iy, iz );  // Get position of domain in global topology

      skip=0;
 
      // Lower x face
      if ( ix==0 ) {
        for ( j=1; j<=grid->ny; ++j ) {
          for ( k=1; k<=grid->nz; ++k ) {
            float e1, e2, cb1, cb2;
            // In output, the 2D surface arrays A[j,k] are FORTRAN indexed: 
            // The j quantity varyies fastest, k, slowest. 
            index = int(  ((iy*grid->ny) + j-1)  
                        + ((iz*grid->nz) + k-1) * (grid->ny*global->topology_y) 
                        + skip); 
            k1  = INDEX_FORTRAN_3(1,j+1,k+1,0,grid->nx+1,0,grid->ny+1,0,grid->nz+1);
            k2  = INDEX_FORTRAN_3(2,j+1,k+1,0,grid->nx+1,0,grid->ny+1,0,grid->nz+1);
            e1  = field[k2].ey;
            e2  = field[k2].ez;
            cb1 = 0.5*(field[k1].cby+field[k2].cby); 
            cb2 = 0.5*(field[k1].cbz+field[k2].cbz);
            pvec[index]   = ( e1*cb2-e2*cb1 )*norm; 
            e1vec[index]  = e1;         
            e2vec[index]  = e2;
            cb1vec[index] = cb1;
            cb2vec[index] = cb2;
          }
        }
      }

      if ( global->write_backscatter_only==0 ) {  

        skip+=ncells_yz;
  
        // Upper x face
        if ( ix==global->topology_x-1 ) {
          for ( j=1; j<=grid->ny; ++j ) {
            for ( k=1; k<=grid->nz; ++k ) {
              float e1, e2, cb1, cb2;
              index = int(  ((iy*grid->ny) + j-1)  
                          + ((iz*grid->nz) + k-1) * (grid->ny*global->topology_y) 
                          + skip); 
              k1  = INDEX_FORTRAN_3(grid->nx-1,j+1,k+1,0,grid->nx+1,0,grid->ny+1,0,grid->nz+1);
              k2  = INDEX_FORTRAN_3(grid->nx  ,j+1,k+1,0,grid->nx+1,0,grid->ny+1,0,grid->nz+1);
              e1  = field[k2].ey;
              e2  = field[k2].ez;
              cb1 = 0.5*(field[k1].cby+field[k2].cby); 
              cb2 = 0.5*(field[k1].cbz+field[k2].cbz);
              pvec[index]   = ( e1*cb2-e2*cb1 )*norm; 
              e1vec[index]  = e1;         
              e2vec[index]  = e2;
              cb1vec[index] = cb1;
              cb2vec[index] = cb2;
            }
          }
        }                       
        skip+=ncells_yz;

        // Lower y face
        if ( iy==0 ) {
          for ( j=1; j<=grid->nz; ++j ) {
            for ( k=1; k<=grid->nx; ++k ) {
              float e1, e2, cb1, cb2;
              index = int(  ((iz*grid->nz) + j-1)  
                          + ((ix*grid->nx) + k-1) * (grid->nz*global->topology_z) 
                          + skip); 
              k1  = INDEX_FORTRAN_3(k+1,1,j+1,0,grid->nx+1,0,grid->ny+1,0,grid->nz+1);
              k2  = INDEX_FORTRAN_3(k+1,2,j+1,0,grid->nx+1,0,grid->ny+1,0,grid->nz+1);
              e1  = field[k2].ez;
              e2  = field[k2].ex;
              cb1 = 0.5*(field[k1].cbz+field[k2].cbz); 
              cb2 = 0.5*(field[k1].cbx+field[k2].cbx);
              pvec[index]   = ( e1*cb2-e2*cb1 )*norm; 
              e1vec[index]  = e1;         
              e2vec[index]  = e2;
              cb1vec[index] = cb1;
              cb2vec[index] = cb2;
            }
          }
        }                      
        skip+=ncells_zx;
  
        // Upper y face
        if ( iy==global->topology_y-1 ) {  
          for ( j=1; j<=grid->nz; ++j ) {
            for ( k=1; k<=grid->nx; ++k ) {
              float e1, e2, cb1, cb2;
              index = int(  ((iz*grid->nz) + j-1)  
                          + ((ix*grid->nx) + k-1) * (grid->nz*global->topology_z) 
                          + skip); 
              k1  = INDEX_FORTRAN_3(k+1,grid->ny-1,j+1,0,grid->nx+1,0,grid->ny+1,0,grid->nz+1);
              k2  = INDEX_FORTRAN_3(k+1,grid->ny  ,j+1,0,grid->nx+1,0,grid->ny+1,0,grid->nz+1);
              e1  = field[k2].ez;
              e2  = field[k2].ex;
              cb1 = 0.5*(field[k1].cbz+field[k2].cbz); 
              cb2 = 0.5*(field[k1].cbx+field[k2].cbx);
              pvec[index]   = ( e1*cb2-e2*cb1 )*norm; 
              e1vec[index]  = e1;         
              e2vec[index]  = e2;
              cb1vec[index] = cb1;
              cb2vec[index] = cb2;
            }
          }
        }                         
        skip+=ncells_zx;
  
        // Lower z face
        if ( iz==0 ) {
          for ( j=1; j<=grid->nx; ++j ) {
            for ( k=1; k<=grid->ny; ++k ) {
              float e1, e2, cb1, cb2;
              index = int(  ((ix*grid->nx) + j-1)  
                          + ((iy*grid->ny) + k-1) * (grid->nx*global->topology_x) 
                          + skip); 
              k1  = INDEX_FORTRAN_3(j+1,k+1,1,0,grid->nx+1,0,grid->ny+1,0,grid->nz+1);
              k2  = INDEX_FORTRAN_3(j+1,k+1,2,0,grid->nx+1,0,grid->ny+1,0,grid->nz+1);
              e1  = field[k2].ex;
              e2  = field[k2].ey;
              cb1 = 0.5*(field[k1].cbx+field[k2].cbx); 
              cb2 = 0.5*(field[k1].cby+field[k2].cby);
              pvec[index]   = ( e1*cb2-e2*cb1 )*norm; 
              e1vec[index]  = e1;         
              e2vec[index]  = e2;
              cb1vec[index] = cb1;
              cb2vec[index] = cb2;
            }
          }
        }                      
        skip+=ncells_xy;
  
        // Upper z face
        if ( iz==global->topology_z-1 ) {  
          for ( j=1; j<=grid->nx; ++j ) {
            for ( k=1; k<=grid->ny; ++k ) {
              float e1, e2, cb1, cb2;
              index = int(  ((ix*grid->nx) + j-1)  
                          + ((iy*grid->ny) + k-1) * (grid->nx*global->topology_x) 
                          + skip); 
              k1  = INDEX_FORTRAN_3(j+1,k+1,grid->nz-1,0,grid->nx+1,0,grid->ny+1,0,grid->nz+1);
              k2  = INDEX_FORTRAN_3(j+1,k+1,grid->nz  ,0,grid->nx+1,0,grid->ny+1,0,grid->nz+1);
              e1  = field[k2].ex;
              e2  = field[k2].ey;
              cb1 = 0.5*(field[k1].cbx+field[k2].cbx); 
              cb2 = 0.5*(field[k1].cby+field[k2].cby);
              pvec[index]   = ( e1*cb2-e2*cb1 )*norm; 
              e1vec[index]  = e1;         
              e2vec[index]  = e2;
              cb1vec[index] = cb1;
              cb2vec[index] = cb2;
            }
          }
        }                         
      }
 
      if ( global->write_poynting_sum ) {  
        // Sum poynting flux on surface
        skip=0;
  
        // Lower x face
        for ( i=0, psum[0]=0; i<ncells_yz; ++i ) psum[0]+=pvec[i+skip]; 
        if ( global->write_backscatter_only==0 ) {  

          // Upper x face
          skip+=ncells_yz;
          for ( i=0, psum[1]=0; i<ncells_yz; ++i ) psum[1]+=pvec[i+skip]; 

          // Lower y face
          skip+=ncells_yz;
          for ( i=0, psum[2]=0; i<ncells_zx; ++i ) psum[2]+=pvec[i+skip]; 
   
          // Upper y face
          skip+=ncells_zx;
          for ( i=0, psum[3]=0; i<ncells_zx; ++i ) psum[3]+=pvec[i+skip]; 

          // Lower z face
          skip+=ncells_zx;
          for ( i=0, psum[4]=0; i<ncells_xy; ++i ) psum[4]+=pvec[i+skip]; 

          // Upper z face
          skip+=ncells_xy;
          for ( i=0, psum[5]=0; i<ncells_xy; ++i ) psum[5]+=pvec[i+skip]; 
        }

  
        // Sum over all surfaces
        mp_allsum_d(psum, gpsum, sum_stride, grid->mp);  
        // Divide by number of mesh points summed over
        gpsum[0] /= ncells_yz; 
        if ( global->write_backscatter_only==0 ) {  
          gpsum[1] /= ncells_yz; 
          gpsum[2] /= ncells_zx; 
          gpsum[3] /= ncells_zx; 
          gpsum[4] /= ncells_xy; 
          gpsum[5] /= ncells_xy; 
        }
  
        // Write summed Poynting data
        if ( rank()==0 && global->write_poynting_sum ) { 
          status=fileIO.open( fname_poynting_sum, 
                              (step==global->poynting_interval ? io_write : io_read_write) );
          if ( status==fail ) ERROR(("Could not open file."));
          fileIO.seek( uint64_t(sum_stride*(step/global->poynting_interval-1)*sizeof(double)), 
                       SEEK_SET );
          fileIO.write( gpsum, sum_stride );
          fileIO.close(); 
        } 
      }
      sim_log("** step = "<<step<<" Lower x Poynting = "<<gpsum[0]);  // Dump data to stdout
      if ( global->write_backscatter_only==0 ) {  
        sim_log("**        "<<step<<" Upper x Poynting = "<<gpsum[1]);  // Dump data to stdout
        sim_log("**        "<<step<<" Lower y Poynting = "<<gpsum[2]);  // Dump data to stdout
        sim_log("**        "<<step<<" Upper y Poynting = "<<gpsum[3]);  // Dump data to stdout
        sim_log("**        "<<step<<" Lower z Poynting = "<<gpsum[4]);  // Dump data to stdout
        sim_log("**        "<<step<<" Upper z Poynting = "<<gpsum[5]);  // Dump data to stdout
      }
 
      // FIXME:  Is the paranoia of explicit casts inside fileIO.seek() necessary?  

      if ( global->write_poynting_faces ) {  
        // Sum across all processors to get quantities on each surface, then write from proc 0 
        mp_allsum_d(pvec,   gpvec,   stride, grid->mp);  
        if ( rank()==0 ) { 
          status=fileIO.open( fname_poynting, 
                              (step==global->poynting_interval ? io_write : io_read_write) );
          if ( status==fail ) ERROR(("Could not open file."));
          fileIO.seek( stride * uint64_t(step/global->poynting_interval-1)
                              * uint64_t(sizeof(double)), 
                       SEEK_SET );
          fileIO.write( gpvec, stride );
          fileIO.close(); 
        }
      } 

      if ( global->write_eb_faces ) {  
        // Sum across all processors to get quantities on each surface, then write from proc 0 
        mp_allsum_d(e1vec,  ge1vec,  stride, grid->mp);  
        mp_allsum_d(e2vec,  ge2vec,  stride, grid->mp);  
        mp_allsum_d(cb1vec, gcb1vec, stride, grid->mp);  
        mp_allsum_d(cb2vec, gcb2vec, stride, grid->mp);  

        if ( rank()==0 ) { 
          // Write e1 face data
          status=fileIO.open( fname_e1, 
                              (step==global->poynting_interval ? io_write : io_read_write) );
          if ( status==fail ) ERROR(("Could not open file."));
          fileIO.seek( stride * uint64_t(step/global->poynting_interval-1)
                              * uint64_t(sizeof(double)), 
                       SEEK_SET );
          fileIO.write( ge1vec, stride );
          fileIO.close(); 
    
          // Write e2 data
          status=fileIO.open( fname_e2, 
                              (step==global->poynting_interval ? io_write : io_read_write) );
          if ( status==fail ) ERROR(("Could not open file."));
          fileIO.seek( stride * uint64_t(step/global->poynting_interval-1)
                              * uint64_t(sizeof(double)), 
                       SEEK_SET );
          fileIO.write( ge2vec, stride );
          fileIO.close(); 
    
          // Write cb1 data
          status=fileIO.open( fname_cb1, 
                              (step==global->poynting_interval ? io_write : io_read_write) );
          if ( status==fail ) ERROR(("Could not open file."));
          fileIO.seek( stride * uint64_t(step/global->poynting_interval-1)
                              * uint64_t(sizeof(double)), 
                       SEEK_SET );
          fileIO.write( gcb1vec, stride );
          fileIO.close(); 
    
          // Write cb2 Poynting data
          status=fileIO.open( fname_cb2, 
                              (step==global->poynting_interval ? io_write : io_read_write) );
          if ( status==fail ) ERROR(("Could not open file."));
          fileIO.seek( stride * uint64_t(step/global->poynting_interval-1)
                              * uint64_t(sizeof(double)), 
                       SEEK_SET );
          fileIO.write( gcb2vec, stride );
          fileIO.close(); 
        }
      }
    } // if
  } END_PRIMITIVE; 
#endif // switch for poynting diagnostic 


// How to implement a spatial mask on particles. 

# if 0
// Example of how to sum over only those particles in a cylinder with 
// radius PRMAX cells

/* Set PZMASK to an appropriate value */
// PRMAX is the number of cells each side from the center line z = 0
//#define PRMAX 100
//#define PMASK (pz*pz+py*py<PRMAX*PRMAX)

// Set mask to sum over all particles
#define PMASK 1

//FIXME: Use memset to speed up the zeroing of the distribution
#define PREPARE_VELOCITY_SPACE_DATA(SUFF, NAME)                           \
    {                                                                     \
      species_t *sp;                                                      \
      for (int i=0; i<NVX; ++i)                                           \
        for (int j=0; j<NVZ; ++j)                                         \
          f_##SUFF[i*NVZ+j]=0;                                            \
      sp = find_species_name(NAME, species_list);                         \
      for (int ip=0; ip<sp->np; ip++) {                                   \
        particle_t *p=&sp->p[ip];                                         \
        /* Lots of stuff commented because PMASK only has pz */           \
        int nxp2=grid->nx+2;                                              \
        int nyp2=grid->ny+2;                                              \
        /* int nzp2=grid->nz+2;    */                                     \
        /* Turn index i into separate ix, iy, iz indices */               \
        int iz = p->i/(nxp2*nyp2);                                        \
        int iy = (p->i-iz*nxp2*nyp2)/nxp2;                                \
        /* int ix = p->i-nxp2*(iy+nyp2*iz); */                            \
        /* Compute real particle position from relative coords and grid data */ \
        /* double px = grid->x0+((ix-1)+(p->dx+1)*0.5)*grid->dx; */       \
        double py = grid->y0+((iy-1)+(p->dy+1)*0.5)*grid->dy;             \
        double pz = grid->z0+((iz-1)+(p->dz+1)*0.5)*grid->dz;             \
        float invgamma=1/sqrt(1+p->ux*p->ux+p->uy*p->uy+p->uz*p->uz);     \
        float vx=p->ux*grid->cvac*invgamma;                               \
        float vz=p->uz*grid->cvac*invgamma;                               \
        long ivx=long(vx/dvx_##SUFF+(NVX/2));                             \
        long ivz=long(vz/dvz_##SUFF+(NVZ/2));                             \
        if ( abs(ivx)<NVX && abs(ivz)<NVZ && PMASK ) f_##SUFF[ivx*NVZ+ivz]+=p->q;  \
      }                                                                   \
    }
#endif 

#if 1
  // -------------------------------------------------------------------------
  // Diagnostic to write a 2d vx, vz binned velocity distribution 
  // Note that we have converted, using the relativistic gammas, from momentum
  // (which the code uses) to velocity prior to binning. 
  //
#define NVX 100
#define NVZ 100

//FIXME: Use memset to speed up the zeroing of the distribution
#define PREPARE_VELOCITY_SPACE_DATA(SUFF, NAME)                       \
    {                                                                 \
      species_t *sp;                                                  \
      for (int i=0; i<NVX; ++i)                                       \
        for (int j=0; j<NVZ; ++j)                                     \
          f_##SUFF[i*NVZ+j]=0;                                        \
      sp = find_species_name(NAME, species_list);                     \
      for (int ip=0; ip<sp->np; ip++) {                               \
        particle_t *p=&sp->p[ip];                                     \
        float invgamma=1/sqrt(1+p->ux*p->ux+p->uy*p->uy+p->uz*p->uz); \
        float vx=p->ux*grid->cvac*invgamma;                           \
        float vz=p->uz*grid->cvac*invgamma;                           \
        long ivx=long(vx/dvx_##SUFF+(NVX/2));                         \
        long ivz=long(vz/dvz_##SUFF+(NVZ/2));                         \
        if ( abs(ivx)<NVX && abs(ivz)<NVZ )                           \
          f_##SUFF[ivx*NVZ+ivz]+=p->q;                                \
      }                                                               \
    }

#define DUMP_VELOCITY_DATA(SUFF,LEN,HEADER_SIZE)                    \
    {                                                               \
      status = fileIO_##SUFF.open( fname_##SUFF,                    \
                              (step==global->velocity_interval ? io_write : io_read_write) ); \
      if ( status==fail ) ERROR(("Could not open file."));          \
      fileIO_##SUFF.seek( uint64_t( HEADER_SIZE +                   \
                                    (step/global->velocity_interval-1) \
                                    * LEN * sizeof(float)) ,        \
                          SEEK_SET );                               \
      fileIO_##SUFF.write( f_##SUFF, LEN );                         \
      fileIO_##SUFF.close();                                        \
    }

#define INCLUDE_VELOCITY_HEADER 0
#if INCLUDE_VELOCITY_HEADER 
#  define VELOCITY_HEADER_SIZE (2*sizeof(int)+2*sizeof(float))
#  define WRITE_VELOCITY_HEADER(SUFF)                               \
    {                                                               \
      int nvpts[2] = { NVX, NVZ };                                  \
      float dv[2];                                                  \
      dv[0] = dvx_##SUFF; dv[1] = dvz_##SUFF;                       \
      status = fileIO_##SUFF.open( fname_##SUFF, io_write );        \
      if ( status==fail ) ERROR(("Could not open file."));          \
      fileIO_##SUFF.write( &nvpts, 2 );                             \
      fileIO_##SUFF.write( &dv,    2 );                             \
      fileIO_##SUFF.close();                                        \
    }
#else
#  define VELOCITY_HEADER_SIZE 0
#  define WRITE_VELOCITY_HEADER(SUFF) 
#endif

  BEGIN_PRIMITIVE {
    float f_e[NVX*NVZ], f_He[NVX*NVZ], f_H[NVX*NVZ];
    float vmax_e  = 10*global->vthe,    dvx_e,  dvz_e;
    float vmax_H  = 10*global->vthi_H,  dvx_H,  dvz_H;
    float vmax_He = 10*global->vthi_He, dvx_He, dvz_He;
    //FILE *fp_e, *fp_H, *fp_He;
    FileIO       fileIO_e, fileIO_H, fileIO_He; 
    FileIOStatus status; 

    static char fname_e[256], fname_H[256], fname_He[256];
    dvx_e  = dvz_e  = 2*vmax_e /NVX;
    dvx_H  = dvz_H  = 2*vmax_H /NVX;
    dvx_He = dvz_He = 2*vmax_He/NVX;
    sprintf( fname_e,  "velocity/velocity_e.%i",  (int)rank() );
    sprintf( fname_H,  "velocity/velocity_H.%i",  (int)rank() );
    sprintf( fname_He, "velocity/velocity_He.%i", (int)rank() );

    // Prepare header (if necessary)
    if ( !step ) {
      WRITE_VELOCITY_HEADER(e);
      if ( global->mobile_ions ) {
        if ( global->H_present )  WRITE_VELOCITY_HEADER(H);
        if ( global->He_present ) WRITE_VELOCITY_HEADER(He);
      }
    }
    
    // Bin particle data and write to file. 
    if (step && !(step%global->velocity_interval) ) {
      PREPARE_VELOCITY_SPACE_DATA(e, "electron");
      if ( global->mobile_ions ) {
        if ( global->H_present )  PREPARE_VELOCITY_SPACE_DATA(H, "H");
        if ( global->He_present ) PREPARE_VELOCITY_SPACE_DATA(He, "He");
      }
      DUMP_VELOCITY_DATA(e, NVX*NVZ, VELOCITY_HEADER_SIZE);
      if ( global->mobile_ions ) {
        if ( global->H_present )  DUMP_VELOCITY_DATA(H,  NVX*NVZ, VELOCITY_HEADER_SIZE);
        if ( global->He_present ) DUMP_VELOCITY_DATA(He, NVX*NVZ, VELOCITY_HEADER_SIZE);
      }
    }
  } END_PRIMITIVE; 
#endif 


  // You SHALL use periodic checkpoint/restart on Roadrunner unless the file 
  // system proves to be so flaky that commenting it out as below is necessary to 
  // do any work!

# if 0
  // Restart dump 

  if ( step>0 && should_dump(restart) ) {
    static const char * restart_fbase[2] = { "restart/restart0", "restart/restart1" }; 

    double dumpstart = mp_elapsed(grid->mp);
    dump_restart( restart_fbase[global->rtoggle], 0 ); 
    double dumpelapsed = mp_elapsed(grid->mp) - dumpstart;

	sim_log("Restart duration "<<dumpelapsed);

    global->rtoggle^=1; 
  } 
# endif 

  if ( step>0 && should_dump(restart) ) {
    static const char * restart_fbase[2] = { "restart/restart0", "restart/restart1" };
    double dumpstart = mp_elapsed(grid->mp);

    // Employ turnstiles to partially serialize the writes
    // NUM_TURNSTILES is define above
    begin_turnstile(NUM_TURNSTILES);
    dump_restart( restart_fbase[global->rtoggle], 0 );
    end_turnstile;

    double dumpelapsed = mp_elapsed(grid->mp) - dumpstart;

    sim_log("Restart duration "<<dumpelapsed);

    global->rtoggle^=1;
  }


#if 0
  if ( step == 30470 ) {
    dump_restart( "restart/restart", 0 ); 
    sim_log( "Restart dump restart completed." ); 
  }
# endif 

#if 0
  if ( step>0 && global->quota_check_interval && (step%global->quota_check_interval)==0 ) { 
    if ( mp_elapsed(grid->mp) > global->quota_sec ) {
      dump_restart( "restart/restart", 0 ); 
      sim_log( "Restart dump restart completed." ); 
      sim_log( "Allowed runtime exceeded for this job.  Terminating." ); 
      mp_barrier( grid->mp ); // Just to be safe
      mp_finalize( grid->mp ); 
      exit(0); 
    }
  }
# endif 

  if ( step>0 && global->quota_check_interval && (step%global->quota_check_interval)==0 ) {
    if ( mp_elapsed(grid->mp) > global->quota_sec ) {

      // Employ turnstiles to partially serialize the writes
      // NUM_TURNSTILES is define above
      begin_turnstile(NUM_TURNSTILES);
      dump_restart( "restart/restart", 0 );
      end_turnstile;

      sim_log( "Restart dump restart completed." );
      sim_log( "Allowed runtime exceeded for this job.  Terminating." );
      mp_barrier( grid->mp ); // Just to be safe
      mp_finalize( grid->mp );
      exit(0);
    }
  }

} 
begin_field_injection { 
  // Inject a light wave from lhs boundary with E aligned along y
  // Use scalar diffraction theory for the Gaussian beam source.  (This is approximate). 

  // For quiet startup (i.e., so that we don't propagate a delta-function noise
  // pulse at time t=0) we multiply by a constant phase term exp(i phi) where: 
  //   phi = k*global->xfocus+atan(h)    (3d) 

  // Inject from the left a field of the form ey = e0 sin( omega t )

# define DY    ( grid->y0 + (iy-0.5)*grid->dy - global->ycenter )
# define DZ    ( grid->z0 + (iz-1  )*grid->dz - global->zcenter )
# define R2    ( DY*DY + DZ*DZ )                                   
# define PHASE ( global->omega*t + h*R2/(global->width*global->width) )
# define MASK  ( R2<=pow(global->mask*global->width,2) ? 1 : 0 )

  if ( grid->x0==0 ) {               // Node is on left boundary
    double alpha      = grid->cvac*grid->dt/grid->dx;
    double emax_coeff = (4/(1+alpha))*global->omega*grid->dt*global->e0;
    double prefactor  = emax_coeff*sqrt(2/M_PI); 
    double t          = grid->dt*step; 

    // Compute Rayleigh length in c/wpe
    double rl         = M_PI*global->waist*global->waist/global->lambda; 

    double pulse_shape_factor = 1;
    float pulse_length        = 70;  // units of 1/wpe
    float sin_t_tau           = sin(0.5*t*M_PI/pulse_length);
    pulse_shape_factor        = ( t<pulse_length ? sin_t_tau : 1 );
    double h                  = global->xfocus/rl;   // Distance / Rayleigh length

    // Loop over all Ey values on left edge of this node
    for ( int iz=1; iz<=grid->nz+1; ++iz ) 
      for ( int iy=1; iy<=grid->ny; ++iy )  
        field(1,iy,iz).ey += prefactor 
                             * cos(PHASE) 
                             * exp(-R2/(global->width*global->width)) 
                             * MASK * pulse_shape_factor; 
  }
}


begin_particle_injection {
  // No particle injection for this simulation
}


begin_current_injection {
  // No current injection for this simulation
}

begin_particle_collisions {
  // No particle collisions for this simulation
}

