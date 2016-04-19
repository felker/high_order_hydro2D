//!   \file main.cpp
//    \brief High Order 2D Hydro code
//======================================================================================

// C/C++ headers
#include <stdint.h>   // int64_t
#include <cstdio>     // sscanf()
#include <cstdlib>    // strtol
#include <ctime>      // clock(), CLOCKS_PER_SEC, clock_t
#include <exception>  // exception
#include <iomanip>    // setprecision()
#include <iostream>   // cout, endl
#include <new>        // bad_alloc
#include <string>     // string

// Athena++ headers
#include "athena.hpp"
#include "globals.hpp"
#include "mesh.hpp"
#include "parameter_input.hpp" 
#include "outputs/outputs.hpp"
#include "outputs/wrapper.hpp"
#include "utils/utils.hpp"

#ifdef OPENMP_PARALLEL
#include <omp.h>
#endif

//--------------------------------------------------------------------------------------
//!   \fn int main(int argc, char *argv[]) 
//    \brief Athena++ main program

int main(int argc, char *argv[])
{
  char *input_filename;
  char *prundir = NULL;
  int res_flag=0;   // set to 1 if -r        argument is on cmdline
  int narg_flag=0;  // set to 1 if -n        argument is on cmdline
  int iarg_flag=0;  // set to 1 if -i <file> argument is on cmdline
  int mesh_flag=0;  // set to <nproc> if -m <nproc> argument is on cmdline
  int walltime_flag=0; // wall time flag
  Real wtlim;
  int ncstart=0;

//--- Step 2. --------------------------------------------------------------------------
// Check for command line options and respond.

  for (int i=1; i<argc; i++) {

    // If argv[i] is a 2 character string of the form "-?" then:
    if(*argv[i] == '-'  && *(argv[i]+1) != '\0' && *(argv[i]+2) == '\0'){
      switch(*(argv[i]+1)) {
      case 'i':                      // -i <input_filename>
        ++i;
        if(res_flag==0) input_filename = argv[i];
        iarg_flag = 1;
      break;
      case 'r':                      // -r <restart_file>
        res_flag = 1;
        input_filename = argv[++i];
        break;
      case 'd':                      // -d <run_directory>
        prundir = argv[++i];
        break;
      case 'n':
        narg_flag = 1;
        break;
      case 'm':
        mesh_flag = std::strtol(argv[++i],NULL,10);
        break;
      case 't':
        int wth, wtm, wts;
        std::sscanf(argv[++i],"%d:%d:%d",&wth,&wtm,&wts);
        walltime_flag=1;
        wtlim=(Real)(wth*3600+wtm*60+wts);
        break;
      case 'c':
        if(Globals::my_rank==0) ShowConfig();
        return(0);
      break;
      case 'h':
      default:
        if(Globals::my_rank==0) {
          std::cout<<"Athena++ "<< athena_version <<std::endl;
          std::cout<<"Usage: "<<argv[0]<<" [options] [block/par=value ...]"<<std::endl;
          std::cout<<"Options:" << std::endl;
          std::cout<<"  -i <file>       specify input file [athinput]"<<std::endl;
          std::cout<<"  -r <file>       restart with this file"<<std::endl;
          std::cout<<"  -d <directory>  specify run dir [current dir]"<<std::endl;
          std::cout<<"  -n              parse input file and quit"<<std::endl;
          std::cout<<"  -c              show configuration and quit"<<std::endl;
          std::cout<<"  -m <nproc>      output mesh structure and quit"<<std::endl;
          std::cout<<"  -t hh:mm:ss     wall time limit for final output" << std::endl;
          std::cout<<"  -h              this help"<<std::endl;
          ShowConfig();
        }
        return(0);
        break;
      }
    } // else if argv[i] not of form "-?" ignore it here (tested in ModifyFromCmdline)
  }

// Note steps 3-6 are protected by a simple error handler
//--- Step 3. --------------------------------------------------------------------------
// Construct object to store input parameters, then parse input file and command line.
// With MPI, the input is read by every process in parallel using MPI-IO.

  ParameterInput *pinput;
  IOWrapper infile;
  try {
    pinput = new ParameterInput;
    infile.Open(input_filename,WRAPPER_READ_MODE);
    pinput->LoadFromFile(infile);
    pinput->ModifyFromCmdline(argc,argv);
  	// leave the input file open (for restarting)
  } 
  catch(std::bad_alloc& ba) {
    std::cout << "### FATAL ERROR in main" << std::endl
              << "memory allocation failed initializing class ParameterInput: " 
              << ba.what() << std::endl;
    infile.Close();
    return(0);
  }
  catch(std::exception const& ex) {
    std::cout << ex.what() << std::endl;  // prints diagnostic message  
    infile.Close();
    return(0);
  }

  // Dump input parameters and quit if code was run with -n option.
  if (narg_flag){
    if(Globals::my_rank==0) pinput->ParameterDump(std::cout);
    infile.Close();
    return(0);
  }

//--- Step 4. --------------------------------------------------------------------------
// Construct and initialize Mesh

  Mesh *pmesh;
  try {
    if(res_flag==0)
      pmesh = new Mesh(pinput, mesh_flag);
    else { 
      pmesh = new Mesh(pinput, infile, mesh_flag);
      ncstart=pmesh->ncycle;
    }
    infile.Close(); // close the input file here
  }
  catch(std::bad_alloc& ba) {
    std::cout << "### FATAL ERROR in main" << std::endl
              << "memory allocation failed initializing class Mesh: " 
              << ba.what() << std::endl;
    return(0);
  }
  catch(std::exception const& ex) {
    std::cout << ex.what() << std::endl;  // prints diagnostic message
    return(0);
  }

  // Quit if -m was on cmdline.  This option builds and outputs mesh structure.
  if (mesh_flag>0){
    return(0);
  }

//--- Step 5. --------------------------------------------------------------------------
// Set initial conditions by calling problem generator, or reading restart file

  try {
    pmesh->Initialize(res_flag, pinput);
  } 
  catch(std::bad_alloc& ba) {
    std::cout << "### FATAL ERROR in main" << std::endl << "memory allocation failed "
              << "in problem generator " << ba.what() << std::endl;
    return(0);
  }
  catch(std::exception const& ex) {
    std::cout << ex.what() << std::endl;  // prints diagnostic message
    return(0);
  }

//--- Step 6. --------------------------------------------------------------------------
// Change to run directory, initialize outputs object, and make output of ICs

  Outputs *pouts;
  try {
    ChangeRunDir(prundir);
    pouts = new Outputs(pmesh, pinput);
    if(res_flag==0) pouts->MakeOutputs(pmesh,pinput);
  } 
  catch(std::bad_alloc& ba) {
    std::cout << "### FATAL ERROR in main" << std::endl
              << "memory allocation failed setting initial conditions: " 
              << ba.what() << std::endl;
    return(0);
  }
  catch(std::exception const& ex) {
    std::cout << ex.what() << std::endl;  // prints diagnostic message
    return(0);
  }

//--- Step 9. === START OF MAIN INTEGRATION LOOP =======================================
// For performance, there is no error handler protecting this step (except outputs)


  if(walltime_flag==1)
    WallTimeLimit::InitWTLimit();

  if(Globals::my_rank==0) {
    std::cout<<std::endl<<"Setup complete, entering main loop..."<<std::endl<<std::endl;
  }

  clock_t tstart = clock();
#ifdef OPENMP_PARALLEL
  double omp_start_time = omp_get_wtime();
#endif

  while ((pmesh->time < pmesh->tlim) && 
         (pmesh->nlim < 0 || pmesh->ncycle < pmesh->nlim)){

    if(Globals::my_rank==0) {
      std::cout << "cycle=" << pmesh->ncycle << std::scientific <<std::setprecision(14)
                << " time=" << pmesh->time << " dt=" << pmesh->dt <<std::endl;
    }

    pmesh->UpdateOneStep();

    pmesh->ncycle++;
    pmesh->time += pmesh->dt;

//kyle: first ref to AMR. DELETED

    pmesh->NewTimeStep();

    try {
      pouts->MakeOutputs(pmesh,pinput);
    } 
    catch(std::bad_alloc& ba) {
      std::cout << "### FATAL ERROR in main" << std::endl
                << "memory allocation failed during output: " << ba.what() <<std::endl;
      return(0);
    }
    catch(std::exception const& ex) {
      std::cout << ex.what() << std::endl;  // prints diagnostic message
      return(0);
    }

    // check the wall time limit
    if(walltime_flag==1) {
      if(Globals::my_rank==0) {
        clock_t tnow = clock();
        Real wtnow = (Real)(tnow-tstart)/(Real)CLOCKS_PER_SEC;
        if(wtnow > wtlim && (pmesh->nlim-pmesh->ncycle>2 || pmesh->nlim<0)
                         && (pmesh->tlim-pmesh->time>3.0*pmesh->dt)) {
          walltime_flag=2;
          pmesh->nlim=pmesh->ncycle+2;
          WallTimeLimit::SendWTLimit(pmesh->nlim);
        }
      }
      else {
        if(WallTimeLimit::TestWTLimit(pmesh->nlim)==true)
          walltime_flag=2;
      }
    }

  } // END OF MAIN INTEGRATION LOOP ====================================================

  pmesh->UserWorkAfterLoop(pinput);

  if(walltime_flag==2) { // hit the wall time limit
    try {
      pouts->MakeOutputs(pmesh,pinput,true);
    } 
    catch(std::bad_alloc& ba) {
      std::cout << "### FATAL ERROR in main" << std::endl
                << "memory allocation failed during output: " << ba.what() <<std::endl;
      return(0);
    }
    catch(std::exception const& ex) {
      std::cout << ex.what() << std::endl;  // prints diagnostic message
      return(0);
    }
  }
  // print diagnostic messages
  if(Globals::my_rank==0) {
    std::cout << "cycle=" << pmesh->ncycle << " time=" << pmesh->time
              << " dt=" << pmesh->dt << std::endl;
    if (pmesh->ncycle == pmesh->nlim) {
      if(walltime_flag==2)
        std::cout << std::endl << "Terminating on wall-time limit" << std::endl;
      else 
        std::cout << std::endl << "Terminating on cycle limit" << std::endl;
    } else {
      std::cout << std::endl << "Terminating on time limit" << std::endl;
    }
    std::cout << "time=" << pmesh->time << " cycle=" << pmesh->ncycle << std::endl;
    std::cout << "tlim=" << pmesh->tlim << " nlim=" << pmesh->nlim << std::endl;

    // Calculate and print the zone-cycles/cpu-second and wall-second
#ifdef OPENMP_PARALLEL
    double omp_time = omp_get_wtime() - omp_start_time;;
#endif
    clock_t tstop = clock();
    float cpu_time = (tstop>tstart ? (float)(tstop-tstart) : 1.0)/(float)CLOCKS_PER_SEC;
    int64_t zones = pmesh->GetTotalCells();
    float zc_cpus = (float)(zones*(pmesh->ncycle-ncstart))/cpu_time;

    std::cout << std::endl << "cpu time used  = " << cpu_time << std::endl;
    std::cout << "zone-cycles/cpu_second = " << zc_cpus << std::endl;
#ifdef OPENMP_PARALLEL
    float zc_omps = (float)(zones*(pmesh->ncycle-ncstart))/omp_time;
    std::cout << std::endl << "omp wtime used = " << omp_time << std::endl;
    std::cout << "zone-cycles/omp_wsecond = " << zc_omps << std::endl;
#endif
  }

  delete pinput;
  delete pmesh;
  delete pouts;

  return(0); 
}