/* Modified tracer header file. 
 * 
 * Changes to legacy version: 
 * - Now a proper header file with safeties for multiple inclusion
 * - Added user-configurable number of fields and a facility for 
 *   adding user-defined fields to tracer particle and trajectory
 *   dumps
 * - Added status checks on files to make more graceful error
 *   reporting when read/write files can't be opened 
 * - Removed first-pass buffered write macros
 * - Added loop to advance_tracers boundary_p calls to reflect 
 *   num_comm_round semantics
 * - Added hydro accumulation macros + switch to turn these on 
 *   as necessary 
 * - Added undefs of several of the helper macros local to this
 *   file. 
 * - Tracers changed to dump checkpoint/restart into 
 *   restart/restart[tag] where tag = 0, 1, 2
 * - We now use fseek (-like) file writing of trajectory files 
 *   to ensure we are checkpoint/restart safe without introducing
 *   glitches into the traj files
 *   
 * FIXME:
 * - Add typesafing of macro parameters for safety. (Obviously,
 *   this will be limited - as C/C++ macros are just a lexical 
 *   substitution, any exceptions would occur at runtime.)  
 *
 * Mods added by Brian Albright, XTD-PRI, 4 Aug. 2016 
 *--------------------------------------------------------------*/ 

// Old comments: 
/* This file contains all the necessary routines for tagging
 * and advancing tracer particles in VPIC. The methedology here
 * is to copy all tagged particles into a new tracer species.
 * This tracer species does not back react on the fields/plasma
 * so that the evolution of this tracer species is identical
 * to that of the actual simulation particles. The reason for
 * copying is that we can then use the q subfield in the
 * particle structure as a tag for the tracers, but by altering
 * q we cannot accuratly accumulate moments and thus there
 * can be no back-reaction and so we need to copy. The particle
 * push is unaffected since this depends on q/m which is
 * unmodified and specified external to the particle structure.
 *
 *
 *********************************************************************
 * Example of how these routines should be called form the input deck:
 *********************************************************************
 * begin_globals{
 *      ...
 *      species_t * tracers_list ;
 *      ...
 * }
 *
 * begin_initilization{
 *  ...
 *      define_species(,...);       // Because I'm lazy, tracer species
 *      ...                             // must be defined after all other
 *      ...                 // species !!
 *      tag_tracer(...);        // Call to tag/create a new tracer particle
 *      ...
 *      hijack_tracers(...);        // MUST BE CALLED AFTER DEFINING TRACERS
 *  ...
 * }
 *
 * begin_diagnostics{
 *      ...
 *      tag_tracer(...);                // Call to tag/create a new tracer particle
 *      ...
 *      dump_tracers(...);              // Call to dump tracers when required
 *      ...
 *      dump_tracer_restart();          // Call to dump tracers to a restart file...
 *      r
 *      q
 *      u
 *      ...
 * }
 *
 * begin_particle_injection{
 *      advance_tracers();              // MUST BE CALLED AT EVERY TIME STEP
 *      ...
 *  }
 *
 */

//--------------------------------------------------------------
#ifndef __TRACER_INTERP_HXX__
#define __TRACER_INTERP_HXX__

//--------------------------------------------------------------
// General purpose allocation macro
#ifndef ALLOCATE
#define ALLOCATE(A,LEN,TYPE)                                                  \
  if ( !((A)=(TYPE *)malloc((size_t)(LEN)*sizeof(TYPE))) )                    \
    ERROR(("Cannot allocate."));
#endif // ALLOCATE

//--------------------------------------------------------------
// Turn on some useful output for debugging
#ifndef TRACER_VERBOSE
#define TRACER_VERBOSE 0
#endif // TRACER_VERBOSE

//--------------------------------------------------------------
// Default number of fields written per particle

#ifndef TRACER_NUM_FIELDS_BASE
#define TRACER_NUM_FIELDS_BASE (18) 
#endif // TRACER_NUM_FIELDS_BASE

//--------------------------------------------------------------
// Allow user to override using the following syntax: 
//
// #undef  TRACER_NUM_ADDED_FIELDS 
// #define TRACER_NUM_ADDED_FIELDS (N) 
//
// where N >= 0 

#ifndef TRACER_NUM_ADDED_FIELDS
#define TRACER_NUM_ADDED_FIELDS (0) 
#endif // TRACER_NUM_ADDED_FIELDS

//--------------------------------------------------------------
// User can add additional fields to write per particle. This is 
// potentially useful if data such as, e.g., kinetic energy, Lorentz
// factor, or hydro dump data are desired. Syntax is as follows: 
//
// #undef TRACER_USER_DEFINED_DATA
// #define TRACER_USER_DEFINED_DATA \
//         pout[index + TRACER_NUM_FIELDS_BASE + 0]                           = ... ; \
//         pout[index + TRACER_NUM_FIELDS_BASE + 1]                           = ... ; \
//         ...
//         pout[index + TRACER_NUM_FIELDS_BASE + TRACER_NUM_ADDED_FIELDS - 1] = ... ; 
// 

#ifndef TRACER_USER_DEFINED_DATA
#define TRACER_USER_DEFINED_DATA
#endif // TRACER_USER_DEFINED_DATA

//--------------------------------------------------------------
// Redefine TRACER_DO_ACCUM_HYDRO to 1 in the input deck if we need to accumulate 
// hydro moments before writing trajectory output. Since this involves a pass 
// through all the particles in the system as well as synchronization (which 
// hits MPI), don't do this step unless we need these data! 

#ifndef TRACER_DO_ACCUM_HYDRO
#define TRACER_DO_ACCUM_HYDRO 0
#endif // #ifndef TRACER_ACCUM_HYDRO

//--------------------------------------------------------------
// Declare variables and allocate space for hydro arrays. TAG 
// identifies which hydro species type is used; ID is the 
// species name (type const char[]).  Note: because namespace needs i
// to be declared outside this block, there's no BEGIN_PRIMITIVE


#ifndef TRACER_HYDRO_SETUP
#define TRACER_HYDRO_SETUP(TAG,ID)                                            \
  static int TAG##_initted = 0;                                               \
  static species_t * TAG##_species;                                           \
  static hydro_array_t * ALIGNED(128) TAG##_hydro;                                  \
  if ( ! TAG##_initted ) {                                                    \
    TAG##_initted = 1;                                                        \
    TAG##_species = find_species(ID);                                         \
    TAG##_hydro   = new_hydro_array( grid );                                        \
  } // if initted                                  
#endif // TRACER_HYDRO_SETUP 

//--------------------------------------------------------------
// Macro to accumulate hydro data. TAG identifies which 
// hydro moments before writing trajectory output. Since this involves a pass 
// through all the particles in the system as well as synchronization (which 
// hits MPI), don't do this step unless we need these data! 

#ifndef TRACER_ACCUM_HYDRO       
#define TRACER_ACCUM_HYDRO(TAG)                                               \
  BEGIN_PRIMITIVE {                                                           \
    clear_hydro_array( TAG##_hydro );                                         \
    accumulate_hydro_p( TAG##_hydro,                                          \
                        TAG##_species,                                     \
                        interpolator_array);                                         \
    synchronize_hydro_array( TAG##_hydro);                                   \
  } END_PRIMITIVE                                     
#endif // TRACER_ACCUM_HYDRO 

//--------------------------------------------------------------
// Arguments: p is the particle to copy
//            tracer is the tracer species
//            addr is (float *) pointer (when cast to (int *)) to tag value
//
// Warning:  inject_particle_raw doesn't check for free storage 
// space, so we're silently assuming all nodes have enough room to
// hold all tracers
//
// Note: inject_particle_raw() does not test for out of bounds error - 
//       it's the caller's responsibility to do this. 

#ifndef tag_tracer
#define tag_tracer(p, tracer, q)                                              \
  BEGIN_PRIMITIVE {                                                           \
    inject_particle_raw( tracer, p->dx, p->dy, p->dz, p->i,                   \
                         p->ux, p->uy, p->uz, q );                            \
  } END_PRIMITIVE
#endif // tag_tracer

//--------------------------------------------------------------
// BJA - This makes use of VPIC species semantics. Provided tracers
// are defined last so they're at the head of the species_list linked
// list, we can snip them from the list and evolve them separately
// so they are not accumulated by the field solve. However, this
// puts the onus on the user to handle checkpoint/restart of 
// species in the input deck.
//
// FIXME: It would be cleaner and less error prone if this 
//        this functionality were enabled in the VPIC source. 
// 

#ifndef hijack_tracers
#define hijack_tracers( num_tracer_species )                                  \
  BEGIN_PRIMITIVE{                                                            \
    species_t *s         = species_list ;                                     \
    global->tracers_list = species_list ;                                     \
                                                                              \
    /* advance num_tracer_species-1 elements down species_list */             \
    for(int i = 1; i < num_tracer_species ; ++i) s = s->next;                 \
                                                                              \
    /* set species_list to point to the first physical */                     \
    /* particle species */                                                    \
    species_list = s->next;                                                   \
    s->next      = NULL ;                                                     \
  } END_PRIMITIVE
#endif // hijack_tracers

//--------------------------------------------------------------
// advance_p takes care of all local moves, but to resolve cross-domain
// moves and boundary interactions we need to call boundary_p as well.
//
// Note: We've configured this to pass through the boundary handler
//       num_comm_round times instead of assuming three passes maximum.
//       This is done in preparation for being able to process tracers
//       through randomizing boundary handlers such as maxwellian_reflux
//
// BJA - Note: randomizing boundary handlers (e.g., maxwellian_reflux) 
//       cannot be used with tracer particles presently since tracers' 
//       trajectories will be randomized upon interaction with the 
//       boundary in a non-deterministic fashion. FIXME: in the 
//       new head version, have each particle generate its own 
//       private RNG sequence
//
// BJA - Added dummy field because in boundary_p() if a tracer hits an 
//       absorbing boundary, it will accumulate rho_b in the field array 
//
//       FIXME: watch performance: if bad, consider aligned malloc instead
//       of ALLOCATE (which uses regular malloc)

// NOTE-RFB: We likely don't need this temp summation here, as surely s-> nm is 0?!
// New:
//int temp = s->nm;
//advance_p(s, dummy_a, interpolator_array);	\
//s->nm += temp;
// Old:
//s->nm += advance_p(s->p, s->np, s->q, s->pm, s->max_nm,		\
//             dummy_a, interpolator_array, grid) ;			\



#ifndef advance_tracers
#define advance_tracers() BEGIN_PRIMITIVE {                                   \
  static int a_initted = 0;                                                   \
  static accumulator_array_t *dummy_a;                                              \
  static field_array_t       *dummy_f;                                              \
  if (a_initted == 0) {                                                       \
    dummy_a = new_accumulator_array(grid);                                         \
    dummy_f = new_standard_field_array(grid,material_list,0);  \
  if(step()) read_tracer_restart(global->rtoggle);			\
    a_initted = 1;                                                            \
  } /* if */                                                                  \
                                                                              \
  species_t *s = global->tracers_list ;                                       \
  while( s ) {                             \
    int temp = s->nm; \
    advance_p(s, dummy_a, interpolator_array);	\
    s->nm += temp; \
    for ( int npass=0; npass<num_comm_round; ++npass ) {                      \
      boundary_p(particle_bc_list,s, dummy_f, dummy_a); \
 } /* for */                                                               \
    if (step()%s->sort_interval == 0) sort_p(s);			\
    s = s->next ;                                                             \
  } /* while */                                                               \
} END_PRIMITIVE
#endif // advance_tracers

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}
template <typename T> int is_negative(T val) {
    return (val < T(0));
}

//--------------------------------------------------------------
// Symbols used in macros below - Note: need to undef at the end 
// of include file so we don't pollute namespace 

#define nxg (grid->nx + 2)
#define nyg (grid->ny + 2)
#define nzg (grid->nz + 2)
#define i0 (ii%nxg)
#define j0 ((ii/nxg)%nyg)
#define k0 (ii/(nxg*nyg))
#define tracer_x ((i0 + (dx0-1)*0.5) * grid->dx + grid->x0)
#define tracer_y ((j0 + (dy0-1)*0.5) * grid->dy + grid->y0)
#define tracer_z ((k0 + (dz0-1)*0.5) * grid->dz + grid->z0)

//--------------------------------------------------------------
// This dump routine converts particle positions to global
// rather than local coordinates. User defined custom data
// fields are enabled as described above. 
//
// Note: we also need to store the local coords because conversion 
// from rel. to abs. position loses precision and this can lead to 
// irreproducible particle orbits in the second pass. 

#ifndef dump_tracers
#define dump_tracers(fbase)                                                   \
  BEGIN_PRIMITIVE {                                                           \
    char dname[256], fname[256] ;                                             \
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
    FileIO       fh;                                                          \
    FileIOStatus status;                                                      \
    sprintf(dname, "%s/T.%d", fbase, step );                                  \
    if ( rank()==0 ) dump_mkdir(dname);                                       \
    mp_barrier( grid->mp ); /* prevent race condition */                      \
                                                                              \
    nvar = TRACER_NUM_FIELDS_BASE + TRACER_NUM_ADDED_FIELDS;                  \
    *local_rank = rank();                                                     \
                                                                              \
    while( s ) {                                                              \
      if ( s->np > 0 ) {                                                      \
        pout = new float[s->np*nvar];                                         \
        for ( p=s->p, j=0; j<s->np; j++, p++ ) {                              \
          int index = j*nvar;                                                 \
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
          TRACER_USER_DEFINED_DATA  ;                                         \
        }                                                                     \
        sprintf(fname, "%s/%s.%d", dname , s->name, (int)rank());             \
        status = fh.open(fname, io_write);                                    \
        if ( status == fail ) ERROR(("Could not open file %s", fname));       \
        WRITE(int,    s->np ,         fh ) ;                                  \
        WRITE(int,    nvar,           fh ) ;                                  \
        WRITE(float,  s->q_m ,        fh ) ;                                  \
        WRITE(double, step*grid->dt , fh ) ;                                  \
        WRITE(float,  1.0 ,           fh ) ;                                  \
        WRITE(double, 1.0 ,           fh ) ;                                  \
        fh.write(pout, nvar*s->np);                                           \
        fh.close();                                                           \
        delete [] pout;                                                       \
      }                                                                       \
      s = s->next;                                                            \
    }                                                                         \
  } END_PRIMITIVE
#endif // dump_tracers


//--------------------------------------------------------------
// BJA - is this even used? 
// FIXME: remove if not 

#ifndef tracer_output
#define tracer_output(fbase)                                                  \
  BEGIN_PRIMITIVE {                                                           \
    FileIO       fh;                                                          \
    FileIOStatus status;                                                      \
    sprintf(dname, "%s/T.%d", fbase, step );                                  \
    if ( rank()==0 ) dump_mkdir(dname);                                       \
    mp_barrier( grid->mp ); /* prevent race condition */                      \
    isp = 0;                                                                  \
    while ( s ) {                                                             \
      sprintf(fname, "%s/%s.%d", dname , s->name, (int)rank());               \
      status = fh.open(fname, io_write);                                      \
      if ( status == fail ) ERROR(("Could not open file %s", fname));         \
      WRITE(float, s->q_m, fh);                                               \
      WRITE(int,     nvar, fh);                                               \
      WRITE(int,      tri, fh);                                               \
      fh.write(np + isp*tri, tri);                                            \
      fh.write(t, tri);                                                       \
      fh.write(pout + 2*isp*ntr*nvar*tri, ntraj_point[isp]*nvar);             \
      fh.close();                                                             \
      s = s->next;                                                            \
      isp++;                                                                  \
    }                                                                         \
  } END_PRIMITIVE
#endif // tracer_output

//--------------------------------------------------------------
// In input deck, need to execute this dump_tracer_restart macro 
// every time we do checkpoint/restart. 

#ifndef dump_tracer_restart
#define dump_tracer_restart(RESTART_TAG)                                      \
  BEGIN_PRIMITIVE {                                                           \
    species_t *s = global->tracers_list ;                                     \
    char fname[256] ;                                                         \
    int numspecies = 0;                                                       \
    FileIO       f;                                                           \
    FileIOStatus status;                                                      \
    sprintf(fname, "restart/restart%d/tracer%d.%d", RESTART_TAG,              \
            global->particle_tracing, (int)rank() ) ;                         \
    status = f.open(fname, io_write);                                         \
    if ( status == fail ) ERROR(("Could not open file %s", fname));           \
                                                                              \
    while( s ) {                                                              \
      ++numspecies;                                                           \
      s = s->next;                                                            \
    }                                                                         \
    WRITE(int, numspecies, f);                                                \
    s = global->tracers_list;                                                 \
    while( s ){                                                               \
      WRITE_STRING(s->name, f);                                               \
      WRITE(species_id, s->id, f);                                            \
      WRITE(float, s->q, f);                                                \
      WRITE(float, s->m, f);                                                \
      WRITE(int, s->max_np, f);                                               \
      WRITE(int, s->max_nm, f);                                               \
      WRITE(int, s->sort_interval, f);                                        \
      WRITE(int, s->sort_out_of_place, f);                                    \
      WRITE(int, s->np, f);                                                   \
      if (s->np > 0 ) f.write(s->p, s->np);                                   \
      s = s->next;                                                            \
    }                                                                         \
    f.close();                                                                \
  } END_PRIMITIVE                                                            
#endif // dump_tracer_restart
                                                                             
//--------------------------------------------------------------
// In input deck, need to execute this read_tracer_restart macro 
// before advancing tracers in begin_particle_injection block 
//
// FIXME: Need to test that tag = 2 (code runs to completion & 
//        writes restart dumps) restarts properly

#ifndef read_tracer_restart
#define read_tracer_restart(RESTART_TAG)                                      \
  BEGIN_PRIMITIVE {                                                           \
    char fname[256] ;                                                         \
    FileIO f;                                                                 \
    FileIOStatus status;                                                      \
    int namelen, maxnp, maxnm, sort ;                                         \
    int sorttype, np, numspecies, i ;                                         \
    float q,m,qm ;                                                                \
    species_id id;                                                            \
    char *name;                                                               \
    species_t * s;                                                            \
                                                                              \
    global->tracers_list = NULL ;                                             \
    sprintf(fname, "restart/restart%d/tracer%d.%d", (RESTART_TAG),            \
            global->particle_tracing, (int)rank() ) ;                         \
    status = f.open(fname, io_read);                                          \
    if ( status == fail ) ERROR(("Could not open file %s", fname));           \
                                                                              \
    READ(int, numspecies, f);                                                 \
    for ( i=0 ; i < numspecies ; ++i ) {                                      \
      READ(int, namelen, f) ;                                                 \
      name = (char *)malloc(namelen+1) ;                                      \
      f.read(name, namelen) ;                                                 \
      name[namelen] = '\0' ;                                                  \
      READ(species_id,    id,   f);                                           \
      READ(float,     q,   f);                                               \
      READ(float,     m,   f);                                               \
      READ(int,       maxnp,    f);                                           \
      READ(int,       maxnm,    f);                                           \
      READ(int,       sort,     f);                                           \
      READ(int,       sorttype, f);                                           \
      READ(int,       np,   f);                                               \
      s = append_species( species(name, q, m , maxnp, maxnm, sort, sorttype, grid), &(global->tracers_list)); \
      s->id = id ;                                                            \
      s->np = np ;                                                            \
      if ( np > 0 ) f.read(s->p, np);                                         \
    }                                                                         \
    f.close();                                                                \
  } END_PRIMITIVE
#endif // read_tracer_restart

//--------------------------------------------------------------
// dump tracer by particle trajectory - requires global->tracer2_interval
// to be defined
//
// Note: cleaned up directory-making syntax to avoid race condition and to 
//       prevent spamming mkdir from each proc each dump
//
// Note: Need to work out a buffered write strategy for this to improve 
//       performance. 


#ifndef dump_traj
#define dump_traj(fbase)                                                      \
  BEGIN_PRIMITIVE {                                                           \
    char dname[256], fname[256] ;                                             \
    species_t *s = global->tracers_list ;                                     \
    float ex, ey, ez, bx, by, bz;                                             \
    float dx0, dy0, dz0;                                                      \
    float ux, uy, uz, q;                                                      \
    uint32_t tag[1];                                                          \
    int ii, n, nvar;                                                          \
    int32_t local_rank[1];                                                    \
    const particle_t     * ALIGNED(32) p;                                     \
    const particle_t     * ALIGNED(32) p0;                                    \
    const interpolator_t * ALIGNED(16) f;                                     \
    const grid_t * g;                                                         \
    FileIO       fh;                                                          \
    FileIOStatus status;                                                      \
                                                                              \
    sprintf(dname, "%s", fbase );                                             \
                                                                              \
    if ( step==0 ) {                                                          \
      if ( rank()==0 ) {                                                      \
        dump_mkdir(dname);                                                    \
      }                                                                       \
      mp_barrier( grid->mp ); /* prevent race condition */                    \
    }                                                                         \
                                                                              \
    nvar = TRACER_NUM_FIELDS_BASE + TRACER_NUM_ADDED_FIELDS;                  \
    float pout[TRACER_NUM_FIELDS_BASE + TRACER_NUM_ADDED_FIELDS];             \
    g = grid;                                                                 \
    *local_rank = rank();                                                     \
                                                                              \
    while( s ) {                                                              \
      if ( TRACER_VERBOSE !=0 )                                               \
        sim_log( "Writing trajectories for: "<<s->name );                     \
      n = s->np;                                                              \
      if ( n > 0 ) {                                                          \
        p0 = s->p;                                                            \
        for ( p=p0; n; n--, p++ ) {                                           \
          dx0 = p->dx;                                                        \
          dy0 = p->dy;                                                        \
          dz0 = p->dz;                                                        \
          ii  = p->i;                                                         \
          ux  = p->ux;                                                        \
          uy  = p->uy;                                                        \
          uz  = p->uz;                                                        \
          q   = p->q;                                                         \
          uint32_t tag = *reinterpret_cast<uint32_t*>(&q);                    \
          f = &interpolator(ii);                                              \
          if (tag != 0) {                                                     \
            int index=0;                                                      \
            ex = f->ex + dy0*f->dexdy                                         \
                       + dz0*(f->dexdz+dy0*f->d2exdydz);                      \
            ey = f->ey + dz0*f->deydz                                         \
                       + dx0*(f->deydx+dz0*f->d2eydzdx);                      \
            ez = f->ez + dx0*f->dezdx                                         \
                       + dy0*(f->dezdy+dx0*f->d2ezdxdy);                      \
            bx = f->cbx + dx0*f->dcbxdx;                                      \
            by = f->cby + dy0*f->dcbydy;                                      \
            bz = f->cbz + dz0*f->dcbzdz;                                      \
            pout[0]  = step*grid->dt ;                                        \
            pout[1]  = (float) tracer_x ;                                     \
            pout[2]  = (float) tracer_y ;                                     \
            pout[3]  = (float) tracer_z ;                                     \
            pout[4]  = ux;                                                    \
            pout[5]  = uy;                                                    \
            pout[6]  = uz;                                                    \
            pout[7]  = ex;                                                    \
            pout[8]  = ey;                                                    \
            pout[9]  = ez;                                                    \
            pout[10] = bx;                                                    \
            pout[11] = by;                                                    \
            pout[12] = bz;                                                    \
            pout[13] = p->dx;                                                 \
            pout[14] = p->dy;                                                 \
            pout[15] = p->dz;                                                 \
            pout[16] = *reinterpret_cast<float*>( local_rank );               \
            pout[17] = *reinterpret_cast<float*>( &ii );                      \
            TRACER_USER_DEFINED_DATA;                                         \
                                                                              \
            sprintf(fname, "%s/%s.%i", dname , s->name, tag);                 \
            status = fh.open( fname,                                          \
                              ( step>0 ? io_read_write : io_write ) );        \
            if ( status == fail ) {                                           \
                ERROR(("Could not open file %s", fname));                     \
            }                                                                 \
            fh.seek( uint64_t( (step/global->tracer2_interval)                \
                               *nvar*sizeof(float) ), SEEK_SET );             \
            fh.write(pout,nvar);                                              \
            fh.close();                                                       \
          } else {                                                            \
            sim_log_local("Mangled tag: "<<tag);                              \
          } /* if */                                                          \
        } /* for */                                                           \
      } /* if */                                                              \
      s = s->next;                                                            \
    }                                                                         \
  } END_PRIMITIVE
#endif // dump_traj

//--------------------------------------------------------------
// FIXME: clean up the syntax to make input decks simpler and less
//        error prone

#ifndef read_tracer
#define read_tracer(arr,npart,fname)                                          \
  BEGIN_PRIMITIVE {                                                           \
    float        qm, const1;                                                  \
    double       timestep, const2;                                            \
    FileIO       fp;                                                          \
    FileIOStatus status;                                                      \
                                                                              \
    status = fp.open( fname, io_read );                                       \
    if ( status == fail )                                                     \
      ERROR(("Could not open file %s", fname));                               \
    fp.read( &npart,    1 );                                                  \
    fp.read( &qm,       1 );                                                  \
    fp.read( &timestep, 1 );                                                  \
    fp.read( &const1,   1 );                                                  \
    fp.read( &const2,   1 );                                                  \
    if ( TRACER_VERBOSE != 0 ) {                                              \
      sim_log_local("npart:                        "<<npart);                 \
      sim_log_local("qm:                           "<<qm);                    \
      sim_log_local("time:                         "<<timestep);              \
      sim_log_local("sanity check (should be 1 1): "<<const1<<" "<<const2);   \
    }                                                                         \
    if ( const1 != 1.0 || const2 != 1.0 ) {                                   \
      ERROR(("Mangled header: const1 = %f, const2 = %f\n", const1, const2));  \
      exit(1);                                                                \
    }                                                                         \
    ALLOCATE( arr, npart*9, float );                                          \
    if ( npart > 0 ) {                                                        \
      fp.read( arr, npart*9 );                                                \
      sim_log_local( npart << " loaded from " << fname );                     \
      if ( TRACER_VERBOSE != 0 ) {                                            \
        for (int tmpval=0; tmpval<npart*9; ++tmpval ) {                       \
          sim_log(tmpval<<" "<<arr[tmpval]);                                  \
        }                                                                     \
      }                                                                       \
    }                                                                         \
  } END_PRIMITIVE 
#endif // read_tracer



/****************************************************************
 ******************* Copied from dumpmacros.h *******************
 ***************************************************************/

#ifndef WRITE
#define WRITE(type,value,fileIO)                                              \
  BEGIN_PRIMITIVE {                                                           \
    type __WRITE_tmp = (type)(value);                                         \
    fileIO.write( &__WRITE_tmp, 1 );                                          \
  } END_PRIMITIVE
#endif // WRITE

// Note: strlen does not include the terminating NULL

#ifndef WRITE_STRING 
#define WRITE_STRING(string,fileIO)                                           \
  BEGIN_PRIMITIVE {                                                           \
    int __WRITE_STRING_len = 0;                                               \
    if( string!=NULL ) __WRITE_STRING_len = strlen(string);                   \
    fileIO.write( &__WRITE_STRING_len, 1 );                                   \
    if( __WRITE_STRING_len>0 ) {                                              \
        fileIO.write( string, __WRITE_STRING_len );                           \
    }                                                                         \
  } END_PRIMITIVE
#endif // WRITE_STRING 

#ifndef READ
#define READ(type,value,fileIO)                                               \
  BEGIN_PRIMITIVE {                                                           \
    type __READ_tmp;                                                          \
    fileIO.read(&__READ_tmp, 1 );                                             \
    (value) = __READ_tmp;                                                     \
  } END_PRIMITIVE
#endif // READ


#endif // __TRACER_HXX__

