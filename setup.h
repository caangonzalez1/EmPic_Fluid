#include <dirent.h>
#include <errno.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

//==========================================
//  Setup all initial conditions 
//==========================================
void setup_fields (GridPoints ***xgc,
FieldPoints ***field,
FieldPoints ***field_copy)
{


        for(int ix=0;ix<Nx;ix++)
            {
                for(int iy=0;iy<Ny;iy++)
                {
                    for(int iz=0;iz<Nz;iz++)
                    {
                        xgc[ix][iy][iz].x = xmin + dx*0.5 + double(ix)*dx;
                        xgc[ix][iy][iz].y = ymin + dy*0.5 + double(iy)*dy;
                        xgc[ix][iy][iz].z = zmin + dz*0.5 + double(iz)*dz;
                    }
                }
            }

        // Initialize cells
        for(int ix=0;ix<Nx+1;ix++)
        {
            for(int iy=0;iy<Ny+1;iy++)
            {
                for(int iz=0;iz<Nz+1;iz++)
                {
                        field[ix][iy][iz].x = xmin + double(ix)*dx;
                        field[ix][iy][iz].y = ymin + double(iy)*dy;
                        field[ix][iy][iz].z = zmin + double(iz)*dz;
                        field[ix][iy][iz].ex = 0.0;
                        field[ix][iy][iz].ey = 0.0;
                        field[ix][iy][iz].ez = 0.0;
                        field[ix][iy][iz].bx = 1.0;
                        field[ix][iy][iz].by = 0.0;
                        field[ix][iy][iz].bz = 0.0;
                        field[ix][iy][iz].jx = 0.0;
                        field[ix][iy][iz].jy = 0.0;
                        field[ix][iy][iz].jz = 0.0;
                        field_copy[ix][iy][iz].ex = field[ix][iy][iz].ex;
                        field_copy[ix][iy][iz].ey = field[ix][iy][iz].ey;
                        field_copy[ix][iy][iz].ez = field[ix][iy][iz].ez;
                        field_copy[ix][iy][iz].bx = field[ix][iy][iz].bx;
                        field_copy[ix][iy][iz].by = field[ix][iy][iz].by;
                        field_copy[ix][iy][iz].bz = field[ix][iy][iz].bz;
                        field_copy[ix][iy][iz].jx = field[ix][iy][iz].jx;
                        field_copy[ix][iy][iz].jy = field[ix][iy][iz].jy;
                        field_copy[ix][iy][iz].jz = field[ix][iy][iz].jz;
                }
            }
        }

    /*
    //CC initial condition
        for(int iy=0;iy<Ny;iy++)
        {
         for(int ix=0;ix<Nx;ix++)
          {
              xgc[ix][iy][0].ex = 0.0;
              xgc[ix][iy][0].ez = 0.0;
              xgc[ix][iy][0].bx = 1.0;
              xgc[ix][iy][0].bz = 0.0;
              
              if(xgc[ix][iy][0].x <= 0.5)
              {
                  xgc[ix][iy][0].ey = 1.0;
                  xgc[ix][iy][0].by = -0.75;
              }
              else
              {
                  xgc[ix][iy][0].ey = -1.0;
                  xgc[ix][iy][0].by = 0.75;
              }
            }
        }
     */
    /*
        double tmp_far;
        //nodes initial condition
        for(int iy=0;iy<Ny+1;iy++)
        {
         for(int ix=0;ix<Nx+1;ix++)
           {
               tmp_far = pow((field[ix][iy][0].x),2) + pow((field[ix][iy][0].y),2);
               if(tmp_far==0)
               {
                   field[ix][iy][0].bz=0.0;
                   field[ix][iy][0].ex=0.0;
                   field[ix][iy][0].ey=0.0;
               }
               else
               {
                   field[ix][iy][0].bz = sqrt(tmp_far);
                   field[ix][iy][0].ex = (0.01*(field[ix][iy][0].x))/pow(tmp_far,1.5);
                   field[ix][iy][0].ey = (0.01*(field[ix][iy][0].y))/pow(tmp_far,1.5);
               }
           }
         }
    */
    
}

void setup
(GridPoints ***xgc,
 FieldPoints ***field,
 FieldPoints ***field_copy,
 vector<ThreeVec> &xi, vector<ThreeVec> &pi,vector<double> &qi,vector<ThreeVec> &xe, vector<ThreeVec> &pe,vector<double> &qe,vector<ThreeVec> &x_old,vector<ThreeVec> &x_new,vector<ThreeVec> &E_p, vector<ThreeVec> &B_p,vector<double> &Vpar,vector<double> &Vperp)
{
    double dxp;
    double dyp;
    double del_x,kp,xa;
    double ranf;
    double xcosine;
    double tmp_far;
    ThreeVec tmp;
    double ranf1,ranf2;
    double vel_0, del_v_x,del_v_y,del_v_z; 


    // No-restart
    if (restart == 0)
    {
	initC = 0;
 
  	for(int ip=0;ip<Npi;ip++)
    	{
              qi.push_back(partweight);
              ranf = ((double) rand()) / ((double) RAND_MAX);
              tmp.setX(xmin + (xmax-xmin-1.0e-15)*ranf);
              ranf1 = ((double) rand()) / ((double) RAND_MAX);
              tmp.setY(ymin + (ymax-ymin-1.0e-15)*ranf1);
              ranf2 = ((double) rand()) / ((double) RAND_MAX);
              tmp.setZ(zmin+(zmax-zmin)*ranf2);

              xi.push_back(tmp);

              vel_0 = gaussian(vth_i);
              del_v_x = vel_0;
              vel_0 = gaussian_speed(vth_i);
              ranf = ((double) rand() / (RAND_MAX));
              del_v_y = vel_0 * sin(2.0*PI0*ranf);
              del_v_z = vel_0 * cos(2.0*PI0*ranf);
              tmp.setX(del_v_x);
              tmp.setY(del_v_y);
              tmp.setZ(del_v_z);
              pi.push_back(tmp);
              
//For output porpuses
              tmp.setX(0.0);
              tmp.setY(0.0);
              tmp.setZ(0.0);
              x_old.push_back(tmp);
              x_new.push_back(tmp);
              E_p.push_back(tmp);
              B_p.push_back(tmp);
              Vpar.push_back(0.0);
              Vperp.push_back(0.0);
        }
        
        for(int ip=0;ip<Npe;ip++)
            {
                  qe.push_back(partweight);
                  ranf = ((double) rand()) / ((double) RAND_MAX);
                  tmp.setX(xmin + (xmax-xmin-1.0e-15)*ranf);
                  ranf1 = ((double) rand()) / ((double) RAND_MAX);
                  tmp.setY(ymin + (ymax-ymin-1.0e-15)*ranf1);
                  ranf2 = ((double) rand()) / ((double) RAND_MAX);
                  tmp.setZ(zmin+(zmax-zmin)*ranf2);
                  xe.push_back(tmp);


                  vel_0 = gaussian(vth_e);
                  del_v_x = vel_0;
                  vel_0 = gaussian_speed(vth_e);
                  ranf = ((double) rand() / (RAND_MAX));
                  del_v_y = vel_0 * sin(2.0*PI0*ranf);
                  del_v_z = vel_0 * cos(2.0*PI0*ranf);
                  if(ip==int(Npe/2)){del_v_x=4.0*vth_e;}
                  else{del_v_x=-4.0*vth_e;}
                  tmp.setX(del_v_x);
                  tmp.setY(del_v_y);
                  tmp.setZ(del_v_z);
                  pe.push_back(tmp);
            }
        cout << "parts total (ppc) :" << (int)pi.size() << " " << (double)((int)pi.size())/((double)Nx*Ny*Nz) <<endl;
        cout << " particle weight " << partweight <<endl;

    }
    // Restart function
    else
    { 
	cout << " -- loading from files -- " << endl;
    }

}
void read_initial_fields
(int my_id ,int p,
 GridPoints ***xgc,int nt)
{
    double *array_tmp;
    int ntt;
    ntt = nt/nout;
    
    ifstream fin;
    char filename[80];
    array_tmp =  (double*)malloc(((Nx*Ny)+1)*sizeof(double));
    double ddd;
    int counter;
    
 //   if(my_id==0)
    {
        sprintf(filename,"./output/Ex_%d.out",ntt);
        fin.open(filename);
        counter = 0;
        while (counter < Nx*Ny)
        {
            fin >> ddd;
    
            array_tmp[counter] = ddd;
            counter++;
        }
        fin.close();
           //share the array_tmp to all the processors
        //MPI_Bcast(&array_tmp,(Nx*Ny)+1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    }
   //     MPI_Barrier(MPI_COMM_WORLD);
        counter=0;
        for(int ix=0;ix<Nx;ix++)
        {
            for(int iy=0;iy<Ny;iy++)
            {
                xgc[ix][iy][0].ex = array_tmp[counter];
                counter++;
            }
        }
    
//      if(my_id==0)
      {
          sprintf(filename,"./output/Ey_%d.out",ntt);
          fin.open(filename);
          counter = 0;
          while (counter < Nx*Ny)
          {
              fin >> ddd;
      
              array_tmp[counter] = ddd;
              counter++;
          }
          fin.close();
             //share the array_tmp to all the processors
          //MPI_Bcast(&array_tmp,(Nx*Ny)+1,MPI_DOUBLE,0,MPI_COMM_WORLD);
     }
     
     //share the array_tmp to all the processors
 //   MPI_Barrier(MPI_COMM_WORLD);
    counter=0;
    for(int ix=0;ix<Nx;ix++)
     {
         for(int iy=0;iy<Ny;iy++)
         {
             xgc[ix][iy][0].ey = array_tmp[counter];
             counter++;
         }
    }
    
 //     if(my_id==0)
       {
           sprintf(filename,"./output/Ez_%d.out",ntt);
           fin.open(filename);
           counter = 0;
           while (counter < Nx*Ny)
           {
               fin >> ddd;
       
               array_tmp[counter] = ddd;
               counter++;
           }
           fin.close();
              //share the array_tmp to all the processors
           //MPI_Bcast(&array_tmp,(Nx*Ny)+1,MPI_DOUBLE,0,MPI_COMM_WORLD);
       }
      
      //share the array_tmp to all the processors
 //   MPI_Barrier(MPI_COMM_WORLD);
    counter=0;
    for(int ix=0;ix<Nx;ix++)
     {
         for(int iy=0;iy<Ny;iy++)
         {
             xgc[ix][iy][0].ez = array_tmp[counter];
             counter++;
         }
    }
    
 //   if(my_id==0)
       {
           sprintf(filename,"./output/Bx_%d.out",ntt);
           fin.open(filename);
           counter = 0;
           while (counter < Nx*Ny)
           {
               fin >> ddd;
       
               array_tmp[counter] = ddd;
               counter++;
           }
           fin.close();
              //share the array_tmp to all the processors
           //MPI_Bcast(&array_tmp,(Nx*Ny)+1,MPI_DOUBLE,0,MPI_COMM_WORLD);
       }
      
      //share the array_tmp to all the processors
  //  MPI_Barrier(MPI_COMM_WORLD);
    counter=0;
    for(int ix=0;ix<Nx;ix++)
     {
         for(int iy=0;iy<Ny;iy++)
         {
             xgc[ix][iy][0].bx = array_tmp[counter];
             counter++;
         }
    }
    
//    if(my_id==0)
    {
        sprintf(filename,"./output/By_%d.out",ntt);
        fin.open(filename);
        counter = 0;
          while (counter < Nx*Ny)
          {
              fin >> ddd;
      
              array_tmp[counter] = ddd;
              counter++;
          }
          fin.close();
           //share the array_tmp to all the processors
        //MPI_Bcast(&array_tmp,(Nx*Ny)+1,MPI_DOUBLE,0,MPI_COMM_WORLD);
      }
     
     //share the array_tmp to all the processors
  //  MPI_Barrier(MPI_COMM_WORLD);
    counter=0;
    for(int ix=0;ix<Nx;ix++)
     {
         for(int iy=0;iy<Ny;iy++)
         {
             xgc[ix][iy][0].by = array_tmp[counter];
             counter++;
         }
    }

 //      if(my_id==0)
       {
           
           sprintf(filename,"./output/Bz_%d.out",ntt);
           fin.open(filename);
           counter = 0;
           while (counter < Nx*Ny)
           {
               fin >> ddd;
       
               array_tmp[counter] = ddd;
               counter++;
           }
           fin.close();
              //share the array_tmp to all the processors
           //MPI_Bcast(&array_tmp,(Nx*Ny)+1,MPI_DOUBLE,0,MPI_COMM_WORLD);
       }
      
      //share the array_tmp to all the processors
  // MPI_Barrier(MPI_COMM_WORLD);
    counter=0;
    for(int ix=0;ix<Nx;ix++)
     {
         for(int iy=0;iy<Ny;iy++)
         {
             //cout << array_tmp[counter];
             //cout << endl;
             xgc[ix][iy][0].bz = array_tmp[counter];
             counter++;
         }
    }
free(array_tmp);
}
