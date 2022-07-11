#include "header.h"
#include "createFolders.h"
#include "setup.h"
#include "definition.h"
#include "setupGrid.h"
#include "functions_PIC.h"
#include "functions_FDTD.h"
#include "functions_HLLC.h"
#include "output.h"
#include <iostream>
#include <fstream>

int main
(int argc, char *argv[])
{
  int my_id, p;
  my_id=0;p=0;
  int ierr=0;
  int istop_final; // error message from hypre solver
  int ilower[2],iupper[2]; // required for hypre domain

  int nt=0;
  double t;

    

  GridPoints ***gc = NULL;     // Eulerian
  FieldPoints ***field = NULL; // Eulerian
  FieldPoints ***field_copy = NULL; // Eulerian
  vector<ThreeVec> xi, pi;     // particles (Ions)
  vector<double>   qi;         // particles (Ions)
  vector<ThreeVec> xe, pe;     // particles (Ions)
  vector<double>   qe;         // particles (Ions)
  vector<double>   Vpar;         // particles (Ions)
  vector<double>   Vperp;         // particles (Ions)
  vector<ThreeVec> x_old;
  vector<ThreeVec> x_new;
  vector<ThreeVec> E_p;        // Efield on particles
  vector<ThreeVec> B_p;        // Bfield on particles
        
  int nparts_i_local;
  int nparts_i_total;
    
  int nparts_e_local;
  int nparts_e_total;
   // averaging global parameters (Boeuf)
   // ====================================
   // (xisize,xesize,anode_i,anode_e,cathode_i,cathode_e)
   int numlocal[6];   // local
   int numglobal[6];  // global
    
    //domain decomposition information
    bool notEast, notWest;
    int peast,pwest;
    bool notNorth, notSouth;
    int  psouth,pnorth;

  // ofstream outtime,outspace,outrho,outne,outni,outey;
  ofstream outt; // (for 1D azimuthal)
  ofstream outte, outfe, outfi;
  ofstream outtec; // for 2D defined in output.h (each output = nt)
  ofstream outtecfluid;
  ofstream outtecfield; // for field 
  ofstream outglobal; // for boeuf
 

  // start MPI
  // =========
  //ierr = MPI_Init(&argc,&argv);
  //ierr = MPI_Comm_size( MPI_COMM_WORLD, &p);
  //ierr = MPI_Comm_rank( MPI_COMM_WORLD, &my_id);
  //MPI_Barrier(MPI_COMM_WORLD);
  printf("%d my_id, %d procs, %d error message\n",my_id,p,ierr);
  // Root should have MPI_COMM_WORLD
  //MPI_Bcast(&p,1,MPI_INT,0,MPI_COMM_WORLD);
  //if(my_id == 0 ) printf("\n\n main # procs = %d \n\n",p);
  setProperties();
  // create folders
  // ==============
  createFolders();
  //MPI_Barrier(MPI_COMM_WORLD);
  printf("create folder complete \n");

  srand(time(NULL));

  // ============
  // Initialize
  // ============
  initialize();

  // check processor information
  // ===========================
  //checkProcs(my_id, p);
  //assignProcs(my_id, p);
  //distributeProcs(MPI_COMM_WORLD,my_id, p);
  getCellSize();
  getCoordinate();
    
  // Initialize (continue)
  // =====================
  initializeArray(gc,field);
  initializeArray(gc,field_copy);
    
  cout << " === Initialize Done (header.h, header_mpi.h) " <<endl;
    
  //FIELDS & PARTICLE SETUP
  // Setup arrays (restart included here)
  // ====================================
    setup(gc,field,field_copy,xi,pi,qi,xe,pe,qe,x_old,x_new,E_p,B_p,Vpar,Vperp);
    setup_fields(gc,field,field_copy);
    //backward
    //leapfrog(xi,pi,x_old,x_new,qi,0);
    //forward
    //leapfrog(xi,pi,x_old,x_new,qi,1);
    //moments_VDF(gc,xi,pi,qi,0,1,1);
    //current_density_zigzag(gc,x_old,x_new);
    
    cout << " === Setup Done (setup.h) " <<endl;
    cout << endl << " Before Iteration " <<endl;
    cout << " nmax: " << nmax << endl;
    cout << " nrec: " << nrec << endl << endl;

  // =====================
  // Global quantities
  // =====================

    outt.open("./global.dat"); //
 
    int nttmp0 = initC;

    int idt_final = 10;
    double _timer_start, _timer_check;
    double _timer[idt_final];
    for(int idt=0;idt<idt_final;++idt)
       _timer[idt] = 0.0;

   //_timer_start = MPI_Wtime();
   if(restart == 0) 
   {
	t = 0;
   }
   else
   { 
	t = initC*dt; 
   }

     cout << endl << "  Iteration START " <<endl << endl;
     cout << " initC is " << initC << " (if restart not zero); t is " << t << endl <<endl;

  // ==========
  // Iteration
  // ==========
  for(nt=initC;nt<=initC+nmax;++nt) // added nmax for output
  {
      //Start the loop checking the number of particles
      nparts_i_total = (int)xi.size();
      nparts_e_total = (int)xe.size();
      cout << "Time passed:" << dt*nt << " Timestep:"<< nt << " proton count: "<< nparts_i_total << " electron count: "<< nparts_e_total << endl ;
      //output ot global quantities
      Temporal_diagnostics(nt,dt*nt,outt,xi,pi,qi,xe,pe,qe,field);
      //output time
      if(nt%nout==0)
      {
         //initial conditions
         if(nt==0)
         {
             if(nt%nrec==0)
             {
                cout << " entering particle writing ... " <<endl;
                write_particles(nt,xi,pi,qi,E_p,B_p,Vpar,Vperp,0);
                write_particles(nt,xe,pe,qe,E_p,B_p,Vpar,Vperp,1);
                cout << " wrote particle  ... [Time: " << endl;
                output_node(nt,t,outtecfield,field);
                 
             }
            // velocity move (half timestep): v^{n}-> v^{n+1/2}
            boris(field,xi,pi,E_p,B_p,0,0,0.5); //protons
            boris(field,xi,pi,E_p,B_p,1,0,0.5); //electrons

            nodes_to_interc(field,field_copy, solver_Maxwell);
            EMfield_halftime(gc,field,field_copy,0,solver_Maxwell); //from n to n-1/2
            interc_to_nodes(field,field_copy,solver_Maxwell);//field_copy = fields on nodes
         }
         else
         {
             if(nt%nrec==0)
             {
                 //velocity move (half timestep): v^{n+1/2}-> v^{n}
                 boris(field_copy,xi,pi,E_p,B_p,0,0,-0.5);
                 boris(field_copy,xi,pi,E_p,B_p,1,0,-0.5);
                 
                 cout << " entering particle writing ... " <<endl;
                 write_particles(nt,xi,pi,qi,E_p,B_p,Vpar,Vperp,0);
                 write_particles(nt,xe,pe,qe,E_p,B_p,Vpar,Vperp,1);
                 cout << " wrote particle  ... "<< endl;
                 //velocity move (half timestep): v^{n}-> v^{n+1/2}
                 boris(field_copy,xi,pi,E_p,B_p,0,0,0.5);
                 boris(field_copy,xi,pi,E_p,B_p,1,0,0.5);
                 output_node(nt,t,outtecfield,field);
             }
         }
      }
      
        EMfield(field,solver_Maxwell);
        interc_to_nodes(field,field_copy,solver_Maxwell);
        // velocity update (one timestep): v^{n+1/2}->v^{n+3/2}
          // ===================================================
        current_reset(field);
        //proton evolution and current deposition
        boris(field_copy,xi,pi,E_p,B_p,0,0,1.0);
        // Advance position: r^{n} --> r^{n+1}
        // ================
        leapfrog(xi,pi,x_old,x_new,qi,1);
        current_zigzag(field,x_old,x_new,xi,0);
        //electrons evolution and current deposition
        boris(field_copy,xe,pe,E_p,B_p,1,0,1.0);
        // Advance position: r^{n} --> r^{n+1}
        // ================
        leapfrog(xe,pe,x_old,x_new,qe,1);
        current_zigzag(field,x_old,x_new,xe,1);

        //moments_VDF(gc,xi,pi,qi,0,1,1);
        //check conservation of number of particles
        if((int)xi.size() != (int)pi.size() || (int)xi.size() != (int)qi.size())
        cout << "ion size error" <<endl;
      if((int)xe.size() != (int)pe.size() || (int)xe.size() != (int)qe.size())
      cout << "ion size error" <<endl;
        t+=dt;

  }
  // ===============
  // After iteration
  // ===============
    outfi.close();
    cout << nmax << " " <<nt << endl;
   return 0;
}

