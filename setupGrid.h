void getCellSize()
{
    nlocal.x = Nx;
    nlocal.y = Ny;
    nlocal.vxe = Nz;
    nlocal_nodes.x = nlocal.x+1;
    nlocal_nodes.y = nlocal.y+1;
    nlocal.vxe = nlocal.vxe+1;

  
   // Start cell in global
   nstart.x = 0;
   nstart.y = 0;
   nstart.vxe = 0;


   // Total size
   ntot.x   = nlocal.x   +2*nbuffer;
   ntot.y =  nlocal.x +2*nbuffer;
   ntot.vxe = nlocal.x +2*nbuffer;
  
   ntot_nodes.x = nlocal_nodes.x +2*nbuffer;
   ntot_nodes.y = nlocal_nodes.y +2*nbuffer;

printf("---------------------------------------\n");
printf("nlocal.x=%d,nlocal.y=%d\n",nlocal.x,nlocal.y);
printf("nlocal_nodes.x=%d,nlocal_nodes.y=%d\n",nlocal_nodes.x,nlocal_nodes.y);
printf("ntot.x=%d, ntot.y=%d\n",ntot.x,ntot.y);
printf("ntot_nodes.x=%d,ntot_nodes.y=%d\n",ntot_nodes.x,ntot_nodes.y);
printf("---------------------------------------\n");

}

void getCoordinate()
{
   ind_x = (int*)malloc(ntot.x*sizeof(int));
   ind_y = (int*)malloc(ntot.y*sizeof(int));
   ind_vxe = (int*)malloc(ntot.vxe*sizeof(int));

   x = (double*)malloc(ntot.x*sizeof(double));
   y = (double*)malloc(ntot.y*sizeof(double));
   vxe = (double*)malloc(ntot.vxe*sizeof(double));

   // set index notation the same as vdf_i and vdf_e
   // takes the buffer cells into account 
 for(int i=0; i<ntot.x; i++)
 {
   ind_x[i] = nstart.x + (i-nbuffer);
   x[i] = xmin + dd.x*ind_x[i] + dd.x/2.0;
 }

 for(int i=0; i<ntot.y; i++)
 {
   ind_y[i] = nstart.y + (i-nbuffer);
   y[i] = ymin + dd.y*ind_y[i] + dd.y/2.0;
 }
   for(int i=0; i<ntot.vxe; i++)
   {
     ind_vxe[i] = nstart.vxe + (i-nbuffer);
     vxe[i] = zmin + dd.vxe*ind_vxe[i] + dd.vxe/2.0;
   }

      cout << "x,y,vxe: " << ntot.x << " "<< ntot.y << " "<< ntot.vxe <<endl;

   // global index
   x_global = (double*)malloc(nx_global*sizeof(double));
   for(int i=0; i<nx_global;i++)
   x_global[i] = xmin +dd.x*(double)(i);
}
