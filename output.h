//==========================================
// Write phase space 
//==========================================
void write_particles(
int nt,
 vector<ThreeVec> &xx, vector<ThreeVec> &px, vector<double> &qx, vector<ThreeVec> &E_p, vector<ThreeVec> &B_p,vector<double> &Vpar,vector<double> &Vperp,int iflag)
{
     int my_id=0;
      int Nptmp;
      ofstream outf;
      char filename[80];

      int nparticle_local;
      int nparticle_total;

      if(iflag==0)       
        sprintf(filename,"./output/ion/parts_%08d_%02d.dat",nt,my_id);
      else if(iflag==1) 
     sprintf(filename,"./output/electron/parts_%08d_%02d.dat",nt,my_id);
      else                cout << "Not defined iflag in write-particles" <<endl; 

        outf.open(filename);
        for(int np=0;np<(int)xx.size();np++)
         outf << xx[np].getX() << " \t"<< xx[np].getY() << " \t" <<  px[np].getX() << " \t" << px[np].getY()  << " \t" << px[np].getZ() << " \t"  << E_p[np].getX() << " \t" << E_p[np].getY() << " \t" << E_p[np].getZ() << " \t" <<  B_p[np].getX() << " \t" << B_p[np].getY()  << " \t" << B_p[np].getZ() << " \t" << Vpar[np] << " \t" << Vperp[np] << endl;
        outf.close(); 
}
//==========================================
// Trace particles 
//==========================================
void write_traceparts
(int my_id, int p, 
 double t, ofstream &out,
 vector<ThreeVec> &xx, vector<ThreeVec> &px, int iflag)
{
     int Np2;
     if(iflag==0)      Np2 = 2;
     else if(iflag==1) Np2 = 2;
      // cout << (int) xx.size() <<endl; 

     for(int ip=Np2;ip<Np2+1;ip++)
           out<< t<<" \t"<<xx[ip].getX() << " \t"<<xx[ip].getY() << " \t"<< px[ip].getX() <<" \t"<<px[ip].getY()<<" \t"<<endl;
}

//==========================================
// Output field
//==========================================
void output_cell
(int my_id, int p,
 int nt, double ttmp, ofstream &outtecfield,
 GridPoints ***xgc)
{
     // Output in root
     int my_id_tmp;
     if(p!=0) my_id_tmp = 0;
     else     my_id_tmp = 0;

     if(my_id==my_id_tmp)
     {
        char filename[80];
        sprintf(filename,"./output/ex_c_%08d.dat",nt);
        outtecfield.open(filename);
        for(int iy=0;iy<Ny;iy++)
        {
            for(int ix=0;ix<Nx;ix++) 
            {
               outtecfield << xgc[ix][iy][0].ex << " \t" << endl; 
            }
        }
        outtecfield.close();
     }

     if(my_id==my_id_tmp)
     {
        char filename[80];
        sprintf(filename,"./output/ey_c_%08d.dat",nt);
        outtecfield.open(filename);
        for(int iy=0;iy<Ny;iy++)
        {
            for(int ix=0;ix<Nx;ix++)
            {
               outtecfield << xgc[ix][iy][0].ey << " \t" << endl;
            }
        }
        outtecfield.close();
     }
     if(my_id==my_id_tmp)
     {
        char filename[80];
        sprintf(filename,"./output/ez_c_%08d.dat",nt);
        outtecfield.open(filename);
        for(int iy=0;iy<Ny;iy++)
        {
            for(int ix=0;ix<Nx;ix++)
            {
               outtecfield << xgc[ix][iy][0].ez << " \t" << endl;
            }
        }
        outtecfield.close();
     }

     if(my_id==my_id_tmp)
     {  
        char filename[80];
        sprintf(filename,"./output/bx_c_%08d.dat",nt);
        outtecfield.open(filename);
        for(int iy=0;iy<Ny;iy++)
        {   
            for(int ix=0;ix<Nx;ix++)
            {  
               outtecfield << xgc[ix][iy][0].bx << " \t" << endl;
            }
        }
        outtecfield.close();
     }
     
     if(my_id==my_id_tmp)
     {  
        char filename[80];
        sprintf(filename,"./output/by_c_%08d.dat",nt);
        outtecfield.open(filename);
        for(int iy=0;iy<Ny;iy++)
        {   
            for(int ix=0;ix<Nx;ix++)
            {  
               outtecfield << xgc[ix][iy][0].by << " \t" << endl;
            }
        }
        outtecfield.close();
     }
     if(my_id==my_id_tmp)
     {  
        char filename[80];
        sprintf(filename,"./output/bz_c_%08d.dat",nt);
        outtecfield.open(filename);
        for(int iy=0;iy<Ny;iy++)
        {   
            for(int ix=0;ix<Nx;ix++)
            {  
               outtecfield << xgc[ix][iy][0].bz << " \t" << endl;
            }
        }
        outtecfield.close();
     }
    if(my_id==my_id_tmp)
        {
           char filename[80];
           sprintf(filename,"./output/jx_c_%08d.dat",nt);
           outtecfield.open(filename);
           for(int iy=0;iy<Ny;iy++)
           {
               for(int ix=0;ix<Nx;ix++)
               {
                  outtecfield << xgc[ix][iy][0].jx << " \t" << endl;
               }
           }
           outtecfield.close();
        }
        
        if(my_id==my_id_tmp)
        {
           char filename[80];
           sprintf(filename,"./output/jy_c_%08d.dat",nt);
           outtecfield.open(filename);
           for(int iy=0;iy<Ny;iy++)
           {
               for(int ix=0;ix<Nx;ix++)
               {
                  outtecfield << xgc[ix][iy][0].jy << " \t" << endl;
               }
           }
           outtecfield.close();
        }
        if(my_id==my_id_tmp)
        {
           char filename[80];
           sprintf(filename,"./output/jz_c_%08d.dat",nt);
           outtecfield.open(filename);
           for(int iy=0;iy<Ny;iy++)
           {
               for(int ix=0;ix<Nx;ix++)
               {
                  outtecfield << xgc[ix][iy][0].jz << " \t" << endl;
               }
           }
           outtecfield.close();
        }
}

void output_node
(int nt, double ttmp, ofstream &outtecfield,
 FieldPoints ***field)
{
     // Output in root
     int my_id_tmp= 0;
     int my_id=0;

     if(my_id==my_id_tmp)
     {
        char filename[80];
        sprintf(filename,"./output/ex_n_%08d.dat",nt);
        outtecfield.open(filename);
        for(int iy=0;iy<Ny+1;iy++)
        {
            for(int ix=0;ix<Nx+1;ix++)
            {
               outtecfield << field[ix][iy][0].ex << " \t" << endl;
            }
        }
        outtecfield.close();
     }

     if(my_id==my_id_tmp)
     {
        char filename[80];
        sprintf(filename,"./output/ey_n_%08d.dat",nt);
        outtecfield.open(filename);
        for(int iy=0;iy<Ny+1;iy++)
        {
            for(int ix=0;ix<Nx+1;ix++)
            {
               outtecfield << field[ix][iy][0].ey << " \t" << endl;
            }
        }
        outtecfield.close();
     }
     if(my_id==my_id_tmp)
     {
        char filename[80];
        sprintf(filename,"./output/ez_n_%08d.dat",nt);
        outtecfield.open(filename);
        for(int iy=0;iy<Ny+1;iy++)
        {
            for(int ix=0;ix<Nx+1;ix++)
            {
               outtecfield << field[ix][iy][0].ez << " \t" << endl;
            }
        }
        outtecfield.close();
     }

     if(my_id==my_id_tmp)
     {
        char filename[80];
        sprintf(filename,"./output/bx_n_%08d.dat",nt);
        outtecfield.open(filename);
        for(int iy=0;iy<Ny+1;iy++)
        {
            for(int ix=0;ix<Nx+1;ix++)
            {
               outtecfield << field[ix][iy][0].bx << " \t" << endl;
            }
        }
        outtecfield.close();
     }
     
     if(my_id==my_id_tmp)
     {
        char filename[80];
        sprintf(filename,"./output/by_n_%08d.dat",nt);
        outtecfield.open(filename);
        for(int iy=0;iy<Ny+1;iy++)
        {
            for(int ix=0;ix<Nx+1;ix++)
            {
               outtecfield << field[ix][iy][0].by << " \t" << endl;
            }
        }
        outtecfield.close();
     }
     if(my_id==my_id_tmp)
     {
        char filename[80];
        sprintf(filename,"./output/bz_n_%08d.dat",nt);
        outtecfield.open(filename);
        for(int iy=0;iy<Ny+1;iy++)
        {
            for(int ix=0;ix<Nx+1;ix++)
            {
               outtecfield << field[ix][iy][0].bz << " \t" << endl;
            }
        }
        outtecfield.close();
     }

    if(my_id==my_id_tmp)
    {
       char filename[80];
       sprintf(filename,"./output/jx_n_%08d.dat",nt);
       outtecfield.open(filename);
       for(int iy=0;iy<Ny+1;iy++)
       {
           for(int ix=0;ix<Nx+1;ix++)
           {
              outtecfield << field[ix][iy][0].jx << " \t" << endl;
           }
       }
       outtecfield.close();
    }
    
    if(my_id==my_id_tmp)
    {
       char filename[80];
       sprintf(filename,"./output/jy_n_%08d.dat",nt);
       outtecfield.open(filename);
       for(int iy=0;iy<Ny+1;iy++)
       {
           for(int ix=0;ix<Nx+1;ix++)
           {
              outtecfield << field[ix][iy][0].jy << " \t" << endl;
           }
       }
       outtecfield.close();
    }
    if(my_id==my_id_tmp)
    {
       char filename[80];
       sprintf(filename,"./output/jz_n_%08d.dat",nt);
       outtecfield.open(filename);
       for(int iy=0;iy<Ny+1;iy++)
       {
           for(int ix=0;ix<Nx+1;ix++)
           {
              outtecfield << field[ix][iy][0].jz << " \t" << endl;
           }
       }
       outtecfield.close();
    }

}


//==========================================
// Output particle moments on the grid
//==========================================
void output_fluid
(int my_id, int p,
 int nt, double ttmp,ofstream &outtecfluid,
 GridPoints ***gc)
{
     // Need to Reduce to root
     int num_save = 12;
     double *dtemp_local;
     double *dtemp_total;
     dtemp_local = (double*)malloc(Nx*Ny*num_save*sizeof(double));
     dtemp_total = (double*)malloc(Nx*Ny*num_save*sizeof(double));

     int in =0;
    
     for(int iy=0;iy<Ny;iy++)
     {
         for(int ix=0;ix<Nx;ix++)
         {
             dtemp_local[in*Nx*Ny+iy*Nx+ix] = gc[ix][iy][0].den_i;
         }
     }
     ++in;
     for(int iy=0;iy<Ny;iy++)
     {
         for(int ix=0;ix<Nx;ix++)
         {
             dtemp_local[in*Nx*Ny+iy*Nx+ix] = gc[ix][iy][0].ux_i;
             
         }
         
     }
    ++in;
     for(int iy=0;iy<Ny;iy++)
     {
         for(int ix=0;ix<Nx;ix++)
         {
             dtemp_local[in*Nx*Ny+iy*Nx+ix] = gc[ix][iy][0].uy_i;
         }
     }
    ++in;
     for(int iy=0;iy<Ny;iy++)
     {
         for(int ix=0;ix<Nx;ix++)
         {
             dtemp_local[in*Nx*Ny+iy*Nx+ix] = gc[ix][iy][0].uz_i;
             
         }
     }
     ++in;
     for(int iy=0;iy<Ny;iy++)
     {
         for(int ix=0;ix<Nx;ix++)
         {
             dtemp_local[in*Nx*Ny+iy*Nx+ix] = gc[ix][iy][0].pxx_i;
         }
     }
     ++in;
     for(int iy=0;iy<Ny;iy++)
     {
         for(int ix=0;ix<Nx;ix++)
         {
             dtemp_local[in*Nx*Ny+iy*Nx+ix] = gc[ix][iy][0].pxy_i;
         }
     }
     ++in;
     for(int iy=0;iy<Ny;iy++)
     {
         for(int ix=0;ix<Nx;ix++)
         {
             dtemp_local[in*Nx*Ny+iy*Nx+ix] = gc[ix][iy][0].pyy_i;
         }
     }
     ++in;
    for(int iy=0;iy<Ny;iy++)
    {
        for(int ix=0;ix<Nx;ix++)
        {
            dtemp_local[in*Nx*Ny+iy*Nx+ix] = gc[ix][iy][0].pxz_i;
        }
    }
    ++in;
    for(int iy=0;iy<Ny;iy++)
    {
        for(int ix=0;ix<Nx;ix++)
        {
            dtemp_local[in*Nx*Ny+iy*Nx+ix] = gc[ix][iy][0].pyz_i;
        }
    }
     ++in;
    for(int iy=0;iy<Ny;iy++)
    {
        for(int ix=0;ix<Nx;ix++)
        {
            dtemp_local[in*Nx*Ny+iy*Nx+ix] = gc[ix][iy][0].pzz_i;
        }
    }
     ++in;
    for(int iy=0;iy<Ny;iy++)
    {
        for(int ix=0;ix<Nx;ix++)
        {
            dtemp_local[in*Nx*Ny+iy*Nx+ix] = gc[ix][iy][0].fxvx;
        }
    }
    ++in;
    for(int iy=0;iy<Ny;iy++)
    {
        for(int ix=0;ix<Nx;ix++)
        {
            dtemp_local[in*Nx*Ny+iy*Nx+ix] = gc[ix][iy][0].fvpape;
        }
    }
    ++in;

     // check consistency
     if(in!=num_save) cout << " wrong num_save " << num_save << " != " << in <<endl;

     // Reduce to root
     int ntotal = Nx*Ny*num_save;
     
    
    // Output in root
    int my_id_tmp;
    if(p!=0) my_id_tmp = 0;
    else     my_id_tmp = 0;
    
    in = 0;

    if(my_id==my_id_tmp)
    {
       char filename[80];
       sprintf(filename,"./output/ion/Density_%08d.dat",nt);
       outtecfluid.open(filename);
       for(int iy=0;iy<Ny;iy++)
       {
           for(int ix=0;ix<Nx;ix++)
           {
              outtecfluid <<  dtemp_total[in*Nx*Ny+iy*Nx+ix] << " \t" << endl;
           }
       }
       outtecfluid.close();
        ++in;
    }
    if(my_id==my_id_tmp)
    {
       char filename[80];
       sprintf(filename,"./output/ion/Ux_%08d.dat",nt);
       outtecfluid.open(filename);
       for(int iy=0;iy<Ny;iy++)
       {
           for(int ix=0;ix<Nx;ix++)
           {
              outtecfluid <<  dtemp_total[in*Nx*Ny+iy*Nx+ix] << " \t" << endl;
           }
       }
       outtecfluid.close();
        ++in;
    }
    if(my_id==my_id_tmp)
    {
       char filename[80];
       sprintf(filename,"./output/ion/Uy_%08d.dat",nt);
       outtecfluid.open(filename);
       for(int iy=0;iy<Ny;iy++)
       {
           for(int ix=0;ix<Nx;ix++)
           {
              outtecfluid << dtemp_total[in*Nx*Ny+iy*Nx+ix]  << " \t" << endl;
           }
       }
       outtecfluid.close();
      ++in;
    }
    if(my_id==my_id_tmp)
    {
       char filename[80];
       sprintf(filename,"./output/ion/Uz_%08d.dat",nt);
       outtecfluid.open(filename);
       for(int iy=0;iy<Ny;iy++)
       {
           for(int ix=0;ix<Nx;ix++)
           {
              outtecfluid << dtemp_total[in*Nx*Ny+iy*Nx+ix]  << " \t" << endl;
           }
       }
       outtecfluid.close();
       ++in;
    }
    if(my_id==my_id_tmp)
    {
       char filename[80];
       sprintf(filename,"./output/ion/Pxx_%08d.dat",nt);
       outtecfluid.open(filename);
       for(int iy=0;iy<Ny;iy++)
       {
           for(int ix=0;ix<Nx;ix++)
           {
              outtecfluid << dtemp_total[in*Nx*Ny+iy*Nx+ix] << " \t" << endl;
           }
       }
       outtecfluid.close();
       ++in;
    }
    if(my_id==my_id_tmp)
    {
       char filename[80];
       sprintf(filename,"./output/ion/Pxy_%08d.dat",nt);
       outtecfluid.open(filename);
       for(int iy=0;iy<Ny;iy++)
       {
           for(int ix=0;ix<Nx;ix++)
           {
              outtecfluid << dtemp_total[in*Nx*Ny+iy*Nx+ix] << " \t" << endl;
           }
       }
       outtecfluid.close();
       ++in;
    }
    if(my_id==my_id_tmp)
    {
       char filename[80];
       sprintf(filename,"./output/ion/Pyy_%08d.dat",nt);
       outtecfluid.open(filename);
       for(int iy=0;iy<Ny;iy++)
       {
           for(int ix=0;ix<Nx;ix++)
           {
              outtecfluid << dtemp_total[in*Nx*Ny+iy*Nx+ix] << " \t" << endl;
           }
       }
       outtecfluid.close();
       ++in;
    }
    if(my_id==my_id_tmp)
      {
         char filename[80];
         sprintf(filename,"./output/ion/Pxz_%08d.dat",nt);
         outtecfluid.open(filename);
         for(int iy=0;iy<Ny;iy++)
         {
             for(int ix=0;ix<Nx;ix++)
             {
                outtecfluid << dtemp_total[in*Nx*Ny+iy*Nx+ix] << " \t" << endl;
             }
         }
         outtecfluid.close();
         ++in;
      }
      if(my_id==my_id_tmp)
      {
         char filename[80];
         sprintf(filename,"./output/ion/Pyz_%08d.dat",nt);
         outtecfluid.open(filename);
         for(int iy=0;iy<Ny;iy++)
         {
             for(int ix=0;ix<Nx;ix++)
             {
                outtecfluid << dtemp_total[in*Nx*Ny+iy*Nx+ix] << " \t" << endl;
             }
         }
         outtecfluid.close();
         ++in;
      }
      if(my_id==my_id_tmp)
      {
         char filename[80];
         sprintf(filename,"./output/ion/Pzz_%08d.dat",nt);
         outtecfluid.open(filename);
         for(int iy=0;iy<Ny;iy++)
         {
             for(int ix=0;ix<Nx;ix++)
             {
                outtecfluid << dtemp_total[in*Nx*Ny+iy*Nx+ix] << " \t" << endl;
             }
         }
         outtecfluid.close();
         ++in;
      }
    if(my_id==my_id_tmp)
       {
          char filename[80];
          sprintf(filename,"./output/ion/Fxvx_%08d.dat",nt);
          outtecfluid.open(filename);
          for(int iy=0;iy<Ny;iy++)
          {
              for(int ix=0;ix<Nx;ix++)
              {
                 outtecfluid << dtemp_total[in*Nx*Ny+iy*Nx+ix] << " \t" << endl;
              }
          }
          outtecfluid.close();
         ++in;
       }
    if(my_id==my_id_tmp)
    {
       char filename[80];
       sprintf(filename,"./output/ion/Fvpape_%08d.dat",nt);
       outtecfluid.open(filename);
       for(int iy=0;iy<Ny;iy++)
       {
           for(int ix=0;ix<Nx;ix++)
           {
              outtecfluid << dtemp_total[in*Nx*Ny+iy*Nx+ix] << " \t" << endl;
           }
       }
       outtecfluid.close();
      ++in;
    }
    
     free(dtemp_local);
     free(dtemp_total);
}

//==========================================
// Particle -> Fluid 
//==========================================
void Temporal_diagnostics
(int nt, double ttmp,ofstream &outt,
 vector<ThreeVec> &xi, vector<ThreeVec> &pi, vector<double> &qi,vector<ThreeVec> &xe, vector<ThreeVec> &pe, vector<double> &qe,
FieldPoints ***field)
{
  double *ux_tot;
  double *dtemp_local;
  double *dtemp_total;
  double *dT_local;
  double *dT_total;

  // Gather space charge 
  // =====================
  int num_save = 2;
  dtemp_local = (double*)malloc(num_save*sizeof(double));
  dtemp_total = (double*)malloc(num_save*sizeof(double));
  dT_local = (double*)malloc(1*sizeof(double)); // for temperature
  dT_total = (double*)malloc(1*sizeof(double)); // for temperature

  for(int in=0; in<num_save; in++)
  {
     dtemp_local[in] = 0.0;
     dtemp_total[in] = 0.0;
  }
  
  // IONS
  for(int ip=0;ip<(int)xi.size();ip++)
  {
      dtemp_total[0]  += 0.5*mass_i/ech*(pow(pi[ip].getX(),2.0)+pow(pi[ip].getY(),2.0)+pow(pi[ip].getZ(),2.0))*qi[ip]/((double)Nx*Ny*Nz);
      dtemp_total[1]  += 0.5*mass_e/ech*(pow(pe[ip].getX(),2.0)+pow(pe[ip].getY(),2.0)+pow(pe[ip].getZ(),2.0))*qe[ip]/((double)Nx*Ny*Nz);
  }
   // Just get it to root


       double eenergy = 0.0;
       double benergy = 0.0;
        
       for(int ix=0;ix<Nx+1;ix++)
       {
           for(int iy=0;iy<Ny+1;iy++)
           {
               eenergy += 0.5*eps*(  pow(field[ix][iy][0].ex,2.0) + pow(field[ix][iy][0].ey,2.0) + pow(field[ix][iy][0].ez,2.0) )/((double)(Nx+1)*(Ny+1));
           }
        }
        for(int ix=0;ix<Nx+1;ix++)
        {
            for(int iy=0;iy<Ny+1;iy++)
            {
                benergy += 0.5*eps*(  pow(field[ix][iy][0].bx,2.0) + pow(field[ix][iy][0].by,2.0) + pow(field[ix][iy][0].bz,2.0) )/((double)(Nx+1)*(Ny+1));
            }
         }

       // Write!
        // Write!
       outt<< ttmp<< " \t";
       for(int in=0; in<num_save; in++){outt<< dtemp_total[in] << " \t";}
       outt<< eenergy<< " \t";
       outt<< benergy<<endl;

    free(dT_local);
    free(dT_total);
}

