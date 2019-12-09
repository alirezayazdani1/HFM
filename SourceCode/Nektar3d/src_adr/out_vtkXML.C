/* created by Leopold Grinberg October 2006  */
/* each CPU writes data in "vtu" format      */
/* zero rank CPU creates a legend file - *.pvd file */

#ifdef NEK2VTK



#include <nektar.h>
#include <stdio.h>
#include <string.h>

#include <stdlib.h>
#include <unistd.h>
#include <sys/param.h>

/* include PDIO and Catamount header files */

#ifdef PDIO

#include <errno.h>
extern "C" {
#include "pdio.h"
}
#include "catamount/cnos_mpi_os.h"

#endif



// added by James 10/26/2006
#include <string>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <vector>
using std::vector;
using std::ostringstream;
using std::string;
using std::copy;
using std::fill_n;
// end added by James

#include <sys/stat.h>
#include <sys/types.h>
#include <fcntl.h>

#if !( defined (__blrts__) || defined (__bg__) )
int mkdir(const char *pathname, mode_t mode);
#endif

void write_vtk_file(Domain *omega);
void write_vtk_file_XML(Domain *omega);

/* pass PDIO data structure to the .pvd file writer */
#ifdef PDIO
void write_pvd_file(int Ndumps, Domain *omega, pdio_info_t info);
#else
void write_pvd_file(int Ndumps, Domain *omega);
#endif
void write_pvtu_file(int ndump);

static char  parent_name[256];
static char  fname_short[256];

void write_vtk_file(Domain *omega){

   FILE *vert_id_file;
   char vert_id_file_name[512];
   int iel,i,j,k;
   Coord X;
   Element *E;
   static int *vert_ID_array;

   static int max_solve_one  = 0, 
              min_solve_zero = 0, 
              max_solve_zero = 0,
              gap            = 0,
              Nvert_total    = 0,
              vert_ID_array_length = 0;

   static int FLAG_INIT = 0;

   sprintf(vert_id_file_name,"file_%d.vtk",mynode());

//   sprintf(vert_id_file_name,"%s_%d.vtk",omega->name,mynode());
   vert_id_file = fopen(vert_id_file_name,"w");

   if (FLAG_INIT == 0) {

     for(iel=0;iel<omega->U->nel;iel++) {
        E = omega->U->flist[iel];
        for (i = 0; i < E->Nverts; i++){
            if (E->vert[i].solve == 0){
              min_solve_zero = E->vert[i].gid;
              break;
            }
        }
        if (min_solve_zero != 0)
            break;
      }

      for(iel=0;iel<omega->U->nel;iel++) {
        E = omega->U->flist[iel];
        for (i = 0; i < E->Nverts; i++){
            if (E->vert[i].solve == 0){
                max_solve_zero = ( (E->vert[i].gid > max_solve_zero) ? E->vert[i].gid : max_solve_zero);
                min_solve_zero = ( (E->vert[i].gid < min_solve_zero) ? E->vert[i].gid : min_solve_zero);
            }
            else
                max_solve_one = ( (E->vert[i].gid > max_solve_one) ? E->vert[i].gid : max_solve_one);
        }
      }
      gap = min_solve_zero-max_solve_one-1;
      Nvert_total = max_solve_zero-gap+1;

      /* count length of vert_ID_array */
      vert_ID_array_length = 0;
      for(iel = 0; iel < omega->U->nel; iel++) {
        E = omega->U->flist[iel];
        vert_ID_array_length  += E->Nverts;
      }

      vert_ID_array = new int[vert_ID_array_length];
      
      k = 0;      
      for(iel = 0; iel < omega->U->nel; iel++) {
        E = omega->U->flist[iel];
        for (i = 0; i < E->Nverts; i++){
          if (E->vert[i].solve == 0)
             j = E->vert[i].gid - gap;
          else
             j = E->vert[i].gid;
          vert_ID_array[k] = j;
          k++;
        }
      } 
      FLAG_INIT = 1;
    } 

    /* allocate memory for coordinates of verices */

    double *XYZ,*x,*y,*z;
    XYZ = new double[Nvert_total*3];
    x = XYZ;
    y = &XYZ[Nvert_total];
    z = &XYZ[Nvert_total*2];

    double *UVWP, *Uvel, *Vvel, *Wvel, *Pres;
    UVWP = new double[Nvert_total*4];
    Uvel = UVWP;
    Vvel = &UVWP[Nvert_total];
    Wvel = &UVWP[Nvert_total*2];
    Pres = &UVWP[Nvert_total*3];


    /* collect coordinates of vertices                           */
    /* gid = "semi-global" ID of Vertex                          */
    /* "semi-global" - since it is global in this partition only */

    k = 0;
    for(iel = 0; iel < omega->U->nel; iel++) {
      E = omega->U->flist[iel];
      for (i = 0; i < E->Nverts; i++){
          j = vert_ID_array[k];
          k++;
          x[j] = E->vert[i].x;
          y[j] = E->vert[i].y;
          z[j] = E->vert[i].z;
      }
    }

   /*  write out geometry data */

   fprintf(vert_id_file,"# vtk DataFile Version 2.0 \n");
   fprintf(vert_id_file,"Unstructured Grid Example \n");
   fprintf(vert_id_file,"ASCII\n");

   fprintf(vert_id_file,"DATASET UNSTRUCTURED_GRID\n");
   fprintf(vert_id_file,"POINTS %d double \n",Nvert_total);

   for (j = 0; j < Nvert_total; j++)
     fprintf(vert_id_file," %.8f  %.8f  %.8f \n",x[j],y[j],z[j]);

    /*  map cells */

    fprintf(vert_id_file,"CELLS %d %d\n",omega->U->nel,omega->U->nel*5);

    k = 0;
    for(iel = 0; iel < omega->U->nel; iel++) {
      E = omega->U->flist[iel];
	  fprintf(vert_id_file,"4 ");
      for (i = 0; i < E->Nverts; i++){
        j = vert_ID_array[k];
        k++;
	fprintf(vert_id_file,"%d  ",j);
      }
      fprintf(vert_id_file," \n");
    }

    fprintf(vert_id_file,"CELL_TYPES %d \n",omega->U->nel);
    for (iel = 0; iel < omega->U->nel; ++iel)
      fprintf(vert_id_file,"%d \n",10);

    /* collect UVWP values */
    k = 0;
    for(iel = 0; iel < omega->U->nel; iel++) {
      E = omega->U->flist[iel];
      for (i = 0; i < E->Nverts; i++){
        j = vert_ID_array[k];
        k++;
        Uvel[j] = E->vert[i].hj[0];
        Vvel[j] = omega->V->flist[iel]->vert[i].hj[0];
        Wvel[j] = omega->W->flist[iel]->vert[i].hj[0];
        Pres[j] = omega->P->flist[iel]->vert[i].hj[0];
     }
   } 

  /* dump data scalar - pressure */

    fprintf(vert_id_file,"POINT_DATA %d \n",Nvert_total);
    fprintf(vert_id_file,"SCALARS pressure double \n");
    fprintf(vert_id_file,"LOOKUP_TABLE default \n");
    
    for (i = 0; i < Nvert_total; i++)
        fprintf(vert_id_file,"%.8f \n",Pres[i]);

    fprintf(vert_id_file,"VECTORS velocity double \n");
    /* dump data - vector Uvel,Vvel,Wvel */
    for (i = 0; i < Nvert_total; i++)
      fprintf(vert_id_file,"%.8f %.8f %.8f \n",Uvel[i],Vvel[i],Wvel[i]);

    fclose (vert_id_file);
 
    delete[] XYZ;
    delete[] UVWP;
}

/* write file in XML format. unstructured grid. Tetarahedral elements. binary output */
void write_vtk_file_XML(Domain *omega){

   FILE *vert_id_file;
   char vert_id_file_name[BUFSIZ], path_name[BUFSIZ];
   int iel,i,j,k;
   Coord X;
   Element *E;
   static int *vert_ID_array;

   static int max_solve_one  = 0, 
              min_solve_zero = 0, 
              max_solve_zero = 0,
              gap            = 0,
              Nvert_total    = 0,
              vert_ID_array_length = 0;

   static int FLAG_INIT = 0;
   static int FLAG_DUMP_INDEX = 0;
   static int SHIFT_DUMP_INDEX = iparam("NVTKSTART");
   static int rank;

   /* PDIO variables */

#ifdef PDIO
   size_t              bufsize = 1024*4000;
   size_t              offset = 0;
   ssize_t             rc;
   static pdio_info_t  info;
#endif
//   mode_t mask;
   umask(022);

   if (FLAG_INIT == 0) {

     /* PDIO configuration - pdio_init is called once */

#ifdef PDIO
     if (iparam("RMTHST")) {
       rank = cnos_get_rank();
       memset(&info, 0, sizeof(info));
       snprintf(info.host, MAXHOSTNAMELEN, omega->remote_hostname);

       info.port         = iparam("NPORT");   /* the default is 5001*/
       info.minTCPchunk  = 1024*1000;         /* so that it just sends *every* time */
       info.maxTCPchunk  = 2*bufsize;         /* throw in some extra to accomodate headers, etc. */
       info.TCPRingSize  = 0x100000*100;      /* 100 MB for the whole ring */
       info.writesize    = bufsize;           /* let the PDIO daemon know how much data is coming; MAXIMUM write size in bytes */

       if (0!=pdio_init(&info)){
         fprintf(stderr, "%d: failed to initialize PDIO\n", rank);
         return;
       }
     }
#else
     rank = mynode();
#endif

     for(iel=0;iel<omega->U->nel;iel++) {
        E = omega->U->flist[iel];
        for (i = 0; i < E->Nverts; i++){
            if (E->vert[i].solve == 0){
              min_solve_zero = E->vert[i].gid;
              break;
            }
        }
        if (min_solve_zero != 0)
            break;
      }

      for(iel=0;iel<omega->U->nel;iel++) {
        E = omega->U->flist[iel];
        for (i = 0; i < E->Nverts; i++){
            if (E->vert[i].solve == 0){
                max_solve_zero = ( (E->vert[i].gid > max_solve_zero) ? E->vert[i].gid : max_solve_zero );
                min_solve_zero = ( (E->vert[i].gid < min_solve_zero) ? E->vert[i].gid : min_solve_zero );
            }
            else
                max_solve_one = ( (E->vert[i].gid > max_solve_one) ? E->vert[i].gid : max_solve_one );
        }
      }
      gap = min_solve_zero-max_solve_one-1;
      Nvert_total = max_solve_zero-gap+1;

      /* count length of vert_ID_array */
      vert_ID_array_length = 0;
      for(iel = 0; iel < omega->U->nel; iel++) {
        E = omega->U->flist[iel];
        vert_ID_array_length  += E->Nverts;
      }

      vert_ID_array = new int[vert_ID_array_length];
      
      k = 0;      
      for(iel = 0; iel < omega->U->nel; iel++) {
        E = omega->U->flist[iel];
        for (i = 0; i < E->Nverts; i++){
          if (E->vert[i].solve == 0)
             j = E->vert[i].gid - gap;
          else
             j = E->vert[i].gid;
          vert_ID_array[k] = j;
          k++;
        }
      }

      char *pch;
      char temp_fname[BUFSIZ];
      sprintf(temp_fname,"%s",omega->name);
      pch = strtok(temp_fname,"/");
      while (pch != NULL){
        memset(fname_short,'\0',256*sizeof(char));
        sprintf(fname_short,"%s",pch);
        pch = strtok(NULL,"/");
      }

      /* create the parent directory */
			sprintf(parent_name,"%s-vtk",omega->name);
			ROOTONLY{
        i = mkdir(parent_name, S_IRWXU | S_IRWXG | S_IRWXO );
        if (i != 0 )
         fprintf(stderr,"Error in creating new directory \"%s\" \n ",parent_name);
			}

      /* create PVD file */
      ROOTONLY{
        j = iparam("NSTEPS")/iparam("VTKSTEPS")+1;
#ifdef PDIO
        write_pvd_file(j, omega, info);
#else
       //write_pvd_file(j, omega);
#endif
      }

      /* create directories for *vtu files */
      ROOTONLY{   
        sprintf(path_name,"%s/%s_time%d",parent_name,fname_short,FLAG_DUMP_INDEX+SHIFT_DUMP_INDEX);
        i = mkdir(path_name, S_IRWXU | S_IRWXG | S_IRWXO );
        if (i != 0 )
         fprintf(stderr,"Error in creating new directory \"%s\" \n ",path_name);
      }
      gsync ();

     FLAG_INIT = 1;
    } 

    /* create pvtu file*/
		ROOTONLY
	    write_pvtu_file(FLAG_DUMP_INDEX+SHIFT_DUMP_INDEX);	// Alireza: pvtu file is preferred

    /* create directory for the next dump , this is needed in order to avoid synchronization */
    ROOTONLY{
      sprintf(path_name,"%s/%s_time%d",parent_name,fname_short,FLAG_DUMP_INDEX+SHIFT_DUMP_INDEX+1);
      i = mkdir(path_name, S_IRWXU | S_IRWXG | S_IRWXO);
      if (i != 0 )
        fprintf(stderr,"Error in creating new directory \"%s\" \n ",path_name);
    }

    /* let last rank to tar previously created directory */    
    if ( rank == (numnodes()-1)  && FLAG_DUMP_INDEX > 0  ){
       char sys_call_string[BUFSIZ];
        /* note that index is not shifted by one */
       sprintf(sys_call_string,"tar -cf  %s/%s_time%d.tar %s/%s_time%d",parent_name,fname_short,
							FLAG_DUMP_INDEX+SHIFT_DUMP_INDEX-1,parent_name,fname_short,FLAG_DUMP_INDEX+SHIFT_DUMP_INDEX-1);
       system(sys_call_string);
    }

    /* set name for *vtu file */
    sprintf(vert_id_file_name,"%s/%s_time%d/%s_%d.vtu",parent_name,fname_short,FLAG_DUMP_INDEX+SHIFT_DUMP_INDEX,fname_short,mynode());
 
    /* allocate memory for coordinates of vertices */
//    float *XYZ;
//    XYZ = new float[Nvert_total*3];
    vector<float> XYZ;
    XYZ.resize(Nvert_total*3);
  
//    float *UVWP, *Pres;
//    UVWP = new float[Nvert_total*4];
//    Pres = &UVWP[Nvert_total*3];
		vector<float> UVW, Pres;
		UVW.resize(Nvert_total*3);
		Pres.resize(Nvert_total);

#ifdef ADR
		static int Nspec = iparam("NSPECS");
		float *Con;
    Con = new float[Nvert_total*Nspec];
#endif

    /* collect coordinates of vertices                           */
    /* gid = "semi-global" ID of Vertex                          */
    /* "semi-global" - since it is global in this partition only */

    k = 0;
    for(iel = 0; iel < omega->U->nel; iel++) {
      E = omega->U->flist[iel];
      for (i = 0; i < E->Nverts; i++){
          j = vert_ID_array[k];
          k++;
          XYZ[j*3] = (float) E->vert[i].x;
          XYZ[j*3+1] = (float) E->vert[i].y;
          XYZ[j*3+2] = (float) E->vert[i].z;
      }
    }

    /* collect UVWP values */
    k = 0;
    for(iel = 0; iel < omega->U->nel; iel++) {
      E = omega->U->flist[iel];
      for (i = 0; i < E->Nverts; i++){
        j = vert_ID_array[k];
        k++;
        UVW[j*3]   = (float) E->vert[i].hj[0];
        UVW[j*3+1] = (float) omega->V->flist[iel]->vert[i].hj[0];
        UVW[j*3+2] = (float) omega->W->flist[iel]->vert[i].hj[0];
        Pres[j] = (float) omega->P->flist[iel]->vert[i].hj[0];
#ifdef ADR
				for (int l = 0; l < Nspec; l++)
					Con[j+l] = (float) omega->T[l]->flist[iel]->vert[i].hj[0];
#endif
      }
    }

    ///////////////////////////////////////////////////////////////
    // handle the header part (not data) of the xml file
    // this part should eventually be written so that the string 
    // is created only once.
    ///////////////////////////////////////////////////////////////

  /*  write header */

  ostringstream xml;

  xml << string("<?xml version=\"1.0\"?>\n")
#if (defined (__bg__) || defined (__blrts__) ) 
      << string("<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"BigEndian\">\n")
#else
      << string("<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n")
#endif
      << string("<UnstructuredGrid>\n")
      << string("<Piece NumberOfPoints=\"")
      << Nvert_total << string("\" NumberOfCells=\"")
      << omega->U->nel << string("\">\n")

  /* start section data at points  */
      << string("<PointData Scalars=\"Pressure\" Vectors=\"Velocity\">\n ")

  /* write scalar data - Pressure */ 
      << string("<DataArray type=\"Float32\" Name=\"Pressure\" NumberOfComponents=\"1\" format=\"appended\" offset=\"0\"/>\n")

  /* write Vector data - Velocity */
      << string("<DataArray type=\"Float32\" Name=\"Velocity\" NumberOfComponents=\"3\" format=\"appended\" offset=\"")
      << (sizeof(float)*Nvert_total+sizeof(int)) << string("\"/>\n")

  /* end section data at points  */
      << string("</PointData>\n")

  /*  write geometry data - coordinates */
      << string("<Points>\n")
      << string("<DataArray type=\"Float32\" Name=\"xyz\" NumberOfComponents=\"3\" format=\"appended\" offset=\"")
      << (sizeof(float)*Nvert_total*4+sizeof(int)*2) << string("\"/>\n")
      << string("</Points>\n")

  /*  map cells */

      << string("<Cells>\n")
      << string("<DataArray type=\"Int32\" Name=\"connectivity\" NumberOfComponents=\"1\" format=\"appended\" offset=\"")
      << (sizeof(float)*Nvert_total*7+sizeof(int)*3) << string("\"/>\n")
  		
 
      << string("<DataArray type=\"Int32\" Name=\"offsets\" format=\"appended\" offset=\"")
      << (sizeof(float)*Nvert_total*7+sizeof(int)*vert_ID_array_length+sizeof(int)*4) << string("\"/>\n")
  

      << string("<DataArray type=\"UInt8\" Name=\"types\"  format=\"appended\" offset=\"")
      << (sizeof(float)*Nvert_total*7+sizeof(int)*(vert_ID_array_length+omega->U->nel)+sizeof(int)*5)
      << string("\"/>\n")
      << string("</Cells>\n")
			
      << string("</Piece>\n")
      << string("</UnstructuredGrid>\n")

  /* write data */
  /* format: 1. write length of DataArray in bytes (j) 2. DataArray   */
  
      << string("<AppendedData encoding=\"raw\">\n_");
			
   string xml_start(xml.str());
   int xml_start_length = xml_start.length();
			
	//////////////////////////////////////
	// handle the end of the xml file
	//////////////////////////////////////

   string xml_end("</AppendedData>\n</VTKFile>\n");
   int xml_end_length = xml_end.length();
	
  //fprintf(vert_id_file,"</AppendedData> \n");
  /* end of write data */

  //fprintf(vert_id_file,"</VTKFile> \n");
  
  //////////////////////////////////////////////////////////////////////
  // handle the "raw" binary section of the xml file
  // (note: techinically a valid xml document is NOT allowed to have 
  // unencoded ie raw binary data)
  //////////////////////////////////////////////////////////////////////

  // the total size of the binary section is the sum of the values of 
  // "size" in the comment blocks preceding the section of code
  int binary_total_size = 4 + sizeof(float)*Nvert_total +
			  4 + sizeof(float)*Nvert_total*3 +
			  4 + sizeof(float)*Nvert_total*3 +
			  4 + sizeof(int)*vert_ID_array_length +
		 	  4 + sizeof(int)*omega->U->nel +
			  4 + sizeof(unsigned char)*omega->U->nel;
													
  // the total length of the file in bytes
  int total_length = xml_start_length + binary_total_size + xml_end_length;

  // the xml file
  vector<char> buffer(total_length);
	
  // copy the strings first - the header and the tail as well
  vector<char>::iterator iter = copy( xml_start.begin(), xml_start.end(), buffer.begin() );
  copy( xml_end.begin(), xml_end.end(), buffer.begin() + xml_start_length + binary_total_size );
	
  // size = 4 + sizeof(float)*Nvert_total
  // j is an int and sizeof(int) == 4 even for x86-64 systems
  char* c;
  j = sizeof(float)*Nvert_total;
  c = reinterpret_cast<char*>(&j);
  iter = copy( c, c + 4, iter );

  // write Scalar (pressure) 
  c = reinterpret_cast<char*>(&Pres[0]);
  iter = copy( c, c + sizeof(float)*Nvert_total, iter );
  
  // size = 4 + sizeof(float)*Nvert_total*3
  j = sizeof(float)*Nvert_total*3;
  c = reinterpret_cast<char*>(&j);
  iter = copy( c, c + 4, iter );

  c = reinterpret_cast<char*>(&UVW[0]);
  iter = copy( c, c + sizeof(float)*Nvert_total*3, iter );
 
  // size = 4 + sizeof(float)*Nvert_total*3
  j = sizeof(float)*Nvert_total*3;
  c = reinterpret_cast<char*>(&j);
  iter = copy( c, c + 4, iter );
  //fwrite(&j,4,1,vert_id_file);
  
  c = reinterpret_cast<char*>(&XYZ[0]);
  iter = copy( c, c + sizeof(float)*Nvert_total*3, iter );
  
  // size = 4 + sizeof(int)*vert_ID_array_length
  j = sizeof(int)*vert_ID_array_length;
  c = reinterpret_cast<char*>(&j);
  iter = copy( c, c + 4, iter );
  //fwrite(&j,4,1,vert_id_file);
  
  k = 0;
  for(iel = 0; iel < omega->U->nel; iel++) {
    E = omega->U->flist[iel];
    for (i = 0; i < E->Nverts; i++){
      j = vert_ID_array[k];
      k++;
      c = reinterpret_cast<char*>(&j);
  		iter = copy( c, c + 4, iter );
      //fwrite(&j,sizeof(int),1,vert_id_file);
    }
  }
   
  // size = 4 + sizeof(int)*omega->U->nel
  j = sizeof(int)*omega->U->nel;
  c = reinterpret_cast<char*>(&j);
  iter = copy( c, c + 4, iter );
  //fwrite(&j,4,1,vert_id_file);

  for(iel = 0; iel < omega->U->nel; iel++){
    j = 4*(iel+1);
    c = reinterpret_cast<char*>(&j);
  	iter = copy( c, c + 4, iter );
    //fwrite(&j,sizeof(int),1,vert_id_file);
  }
  
  // size = 4 + sizeof(unsigned char)*omega->U->nel
  j = sizeof(unsigned char)*omega->U->nel;
  c = reinterpret_cast<char*>(&j);
  iter = copy( c, c + 4, iter );
 
  // define a type of elements - for Tets use type=10
  char cj = 10;
  fill_n( iter, omega->U->nel, cj );

  ////////////////////////
  // Write file :)
  ////////////////////////

#ifdef PDIO 
  /* PDIO write */
  if (iparam("RMTHST")) {
    snprintf(info.rfile_name, MAXPATHLEN, "%s", vert_id_file_name);
    rc = pdio_write(&buffer[0], buffer.size(), offset, &info);
    if (rc != buffer.size()) {
      fprintf(stdout, "%d: pdio_write failed (%ld) writting file locally  \n", rank, rc);
    }
  }
  else{
    fprintf(stdout,"else: RMTHST = %d \n",iparam("RMTHST"));
      vert_id_file = fopen(vert_id_file_name,"wb");
      fwrite( &buffer[0], 1, buffer.size(), vert_id_file );
      fclose(vert_id_file);
      i  = open(vert_id_file_name, O_RDWR);
      fchmod(i,  S_IRWXU | S_IRWXG | S_IROTH | S_IWOTH );
      close(i);
  }
  /* PDIO finish - this is optional */
  if (iparam("RMTHST")) {
    if ( ( FLAG_DUMP_INDEX+1 ) * iparam("VTKSTEPS") >= iparam("NSTEPS") ) {
      if (0!=pdio_fini()) {
	fprintf(stderr, "%d: failed to finalize PDIO\n", rank);
	return;
      }
    }
  }
#else
      vert_id_file = fopen(vert_id_file_name,"wb");
      fwrite( &buffer[0], 1, buffer.size(), vert_id_file );
      fclose(vert_id_file);
      i  = open(vert_id_file_name, O_RDWR);
      fchmod(i,  S_IRWXU | S_IRWXG | S_IROTH | S_IWOTH );
      close(i);
#endif

  FLAG_DUMP_INDEX++;

//  delete[] XYZ;
//  delete[] UVWP;
	XYZ.clear();
	UVW.clear();
	Pres.clear();
#ifdef ADR
	delete[] Con;
#endif

}

#ifdef PDIO
void write_pvd_file(int Ndumps, Domain *omega, pdio_info_t info) {
#else
void write_pvd_file(int Ndumps, Domain *omega){
#endif

  int i,timestep, Ncpu, VTKSTEPS = iparam("VTKSTEPS");
  int SHIFT_DUMP_INDEX = iparam("NVTKSTART");   
  double t  = dparam("t"),
         dt = dparam("DT"),
         timevalue = 0;
         
  FILE *pvd_file;
  char pvd_file_name[BUFSIZ];

  sprintf(pvd_file_name,"%s.pvd",omega->name);

  pvd_file = fopen(pvd_file_name,"w");

  ostringstream xml;
  int xml_start_length,xml_end_length,xml_core_length,xml_total_length;

  /*  create header */
  xml << string("<?xml version=\"1.0\"?>\n")
#if (defined (__bg__) || defined (__blrts__) )
      << string("<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"BigEndian\"> \n")
#else
      << string("<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\"> \n")
#endif

      << string("  <Collection> \n");

  string xml_start(xml.str());
  xml_start_length = xml_start.length();

  /* create tail */
  string xml_end("  </Collection>\n</VTKFile>\n");
  xml_end_length = xml_end.length();

  /* create core */
  string *xml_core;
  xml_core = new string[Ndumps*numnodes()];

  i=0;
  xml_core_length = 0;

  for (timestep = 0; timestep < Ndumps; timestep++){

    timevalue = t+(dt*VTKSTEPS)*timestep; 

    for (Ncpu = 0; Ncpu < numnodes(); Ncpu++){

      ostringstream xml_temp;
#if 1
      xml_temp    << string("    <DataSet timestep=\"")
                  << timestep+SHIFT_DUMP_INDEX  << string("\" timevalue=\"")
                  << timevalue   << string("\" part=\"")
                  << Ncpu        << string("\" file=\"")
                  << omega->name << string("_time")
                  << timestep+SHIFT_DUMP_INDEX  << string("/")
                  << fname_short << string("_")
                  << Ncpu        << string(".vtu\"/> \n");
#else
      xml_temp    << string("    <DataSet timestep=\"")
                  << timestep+SHIFT_DUMP_INDEX  << string("\" timevalue=\"")
                  << timevalue   << string("\" part=\"")
                  << Ncpu        << string("\" file=\"")
                  << omega->name << string("_time")
                  << timestep+SHIFT_DUMP_INDEX  << string("/file_")
                  << Ncpu        << string(".vtu\"/> \n");
#endif

      xml_core[i] = xml_temp.str();
      xml_core_length += xml_core[i].length();
     
      i++;
    }
  }

  /* combine header, core and tail  */

  xml_total_length = xml_start_length+xml_end_length+xml_core_length;

  vector<char> buffer(xml_total_length);
  copy( xml_start.begin(), xml_start.end(), buffer.begin());

  int offset = xml_start_length;
  for (i = 0; i < Ndumps*numnodes(); i++){
     copy( xml_core[i].begin(), xml_core[i].end(), buffer.begin()+offset);
     offset += xml_core[i].length();
  }
  copy( xml_end.begin(), xml_end.end(), buffer.begin()+offset);

  /* output */
  fwrite( &buffer[0], 1, buffer.size(), pvd_file );
  fclose(pvd_file);
  i  = open(pvd_file_name, O_RDWR);
  fchmod(i,  S_IRWXU | S_IRWXG | S_IROTH | S_IWOTH );


#ifdef PDIO
  if (iparam("RMTHST")) {
    int      rank = 1;
    ssize_t  rc;
    offset = 0; 
    snprintf(info.rfile_name, MAXPATHLEN, "%s", pvd_file_name);
    rc = pdio_write(&buffer[0], buffer.size(), offset, &info);
    if (rc != buffer.size()) {
      fprintf(stderr, "%d: pdio_write failed (%ld)\n", rank, rc);
      delete[] xml_core;
      return;
    }
  }
#endif

  delete[] xml_core;

}

void write_pvtu_file(int ndump){

  int i,timestep, Ncpu, VTKSTEPS = iparam("VTKSTEPS");
  int SHIFT_DUMP_INDEX = iparam("NVTKSTART");
  double t  = dparam("t"),
         dt = dparam("DT"),
         timevalue = 0;

  FILE *pvtu_file;
  char pvtu_file_name[BUFSIZ];

  sprintf(pvtu_file_name,"%s/%s_time%d.pvtu",parent_name,fname_short,ndump);

  pvtu_file = fopen(pvtu_file_name,"w");

  ostringstream xml;
  int xml_start_length,xml_end_length,xml_core_length,xml_total_length;

  xml << string("<?xml version=\"1.0\"?>\n")
#if (defined (__bg__) || defined (__blrts__) ) 
      << string("<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"BigEndian\">\n")
#else
      << string("<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n")
#endif
      << string("<PUnstructuredGrid GhostLevel=\"0\">\n")

  /*  write geometry data - coordinates */
      << string("<PPoints>\n")
      << string("  <PDataArray type=\"Float32\" Name=\"xyz\" NumberOfComponents=\"3\"/>\n")
      << string("</PPoints>\n")

  /* start section data at points  */
      << string("<PPointData Scalars=\"Pressure\" Vectors=\"Velocity\">\n")

  /* write scalar data - Pressure */
      << string("  <PDataArray type=\"Float32\" Name=\"Pressure\" NumberOfComponents=\"1\"/>\n")

  /* write Vector data - Velocity */
      << string("  <PDataArray type=\"Float32\" Name=\"Velocity\" NumberOfComponents=\"3\"/>\n")

  /* end section data at points  */
      << string("</PPointData>\n")

  /*  map cells */

      << string("<PCells>\n")
      << string("  <PDataArray type=\"Int32\" Name=\"connectivity\" NumberOfComponents=\"1\"/>\n")

      << string("  <PDataArray type=\"Int32\" Name=\"offsets\" NumberOfComponents=\"1\"/>\n")

      << string("  <PDataArray type=\"UInt8\" Name=\"types\" NumberOfComponents=\"1\"/>\n")
      << string("</PCells>\n");

  string xml_start(xml.str());
  xml_start_length = xml_start.length();

  /* create tail */
  string xml_end("</PUnstructuredGrid>\n</VTKFile>\n");
  xml_end_length = xml_end.length();

  /* create core */
	xml.str("");
  for (Ncpu = 0; Ncpu < numnodes(); Ncpu++){
    xml	  << string("<Piece Source=\"")
					<< parent_name << string("/")
					<< fname_short << string("_time")
          << ndump  << string("/")
          << fname_short << string("_")
          << Ncpu        << string(".vtu\"/>\n");
  }
  string xml_core(xml.str());
  xml_core_length = xml_core.length();

  /* combine header, core and tail  */

  xml_total_length = xml_start_length+xml_end_length+xml_core_length;

  vector<char> buffer(xml_total_length);
  copy( xml_start.begin(), xml_start.end(), buffer.begin());

  int offset = xml_start_length;
  copy( xml_core.begin(), xml_core.end(), buffer.begin()+offset);
  offset += xml_core_length;
  copy( xml_end.begin(), xml_end.end(), buffer.begin()+offset);

  /* output */
  fwrite( &buffer[0], 1, buffer.size(), pvtu_file );
  fclose(pvtu_file);
  i  = open(pvtu_file_name, O_RDWR);
  fchmod(i,  S_IRWXU | S_IRWXG | S_IROTH | S_IWOTH );

}


#endif //end of ifdef NEK2VTK  
