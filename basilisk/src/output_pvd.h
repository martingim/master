/*
This function writes one XML VTK file per PID process of type unstructured grid
(*.vtu) which can be read using Paraview. File stores scalar and vector fields
defined at the center points. Results are recorded on ASCII format. If one writes
one *.vts file per PID process this function may be combined with
output_pvtu_ascii() above to read in parallel. Tested in (quad- and oct-)trees
using MPI. Also works with solids (when not using MPI).
*/
//+ zb[-1,0]*zb[-1,0] + zb[0,-1]*zb[0,-1]+ zb[0,1]*zb[0,1]+ 
#define IN_DOMAIN zb[]*zb[] > 1e-6



@def foreach_xorder()
OMP_PARALLEL() {
    int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
    Point point = { 0 };
    point.level = depth(); point.n = 1 << point.level;
    int _k;
    OMP(omp for schedule(static))
        for (_k = GHOSTS; _k < point.n + GHOSTS; _k++) {
#if dimension == 1
            point.i = _k;
#endif
#if dimension > 1
            point.j = _k;
            for (point.i = GHOSTS; point.i < point.n + GHOSTS; point.i++)
                {
#endif
                    POINT_VARIABLES
                        @
                        @def end_foreach_xorder()
#if dimension > 1
                }
#endif
}
}
@


@def foreach_face_xorder()
OMP_PARALLEL() {
    int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
    Point point = { 0 };
    point.level = depth(); point.n = 1 << point.level;
    int _k;
    OMP(omp for schedule(static))
        for (_k = GHOSTS; _k <= point.n + GHOSTS; _k++) {
#if dimension == 1
            point.i = _k;
#endif

#if dimension > 1
            point.j = _k;
            for (point.i = GHOSTS; point.i <= point.n + GHOSTS; point.i++)
                {
#endif
                    POINT_VARIABLES
                        @
                        @def end_foreach_face_xorder()
#if dimension > 1
                }
#endif
        }
}
@

@def foreach_vertex_xorder()
foreach_face_xorder() {
    x -= Delta / 2.;
#if dimension > 1  
    y -= Delta / 2.;
#endif
    @
    @define end_foreach_vertex_xorder() } end_foreach_face_xorder()





@def foreach_xorder()
OMP_PARALLEL() {
    int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
    Point point = { 0 };
    point.level = depth(); point.n = 1 << point.level;
    int _k;
    OMP(omp for schedule(static))
        for (_k = GHOSTS; _k < point.n + GHOSTS; _k++) {
#if dimension == 1
            point.i = _k;
#endif
#if dimension > 1
            point.j = _k;
            for (point.i = GHOSTS; point.i < point.n + GHOSTS; point.i++)
                {
#endif
                    POINT_VARIABLES
                        @
                        @def end_foreach_xorder()
#if dimension > 1
                }
#endif
}
}
@


@def foreach_face_xorder()
OMP_PARALLEL() {
    int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
    Point point = { 0 };
    point.level = depth(); point.n = 1 << point.level;
    int _k;
    OMP(omp for schedule(static))
        for (_k = GHOSTS; _k <= point.n + GHOSTS; _k++) {
#if dimension == 1
            point.i = _k;
#endif

#if dimension > 1
            point.j = _k;
            for (point.i = GHOSTS; point.i <= point.n + GHOSTS; point.i++)
                {
#endif
                    POINT_VARIABLES
                        @
                        @def end_foreach_face_xorder()
#if dimension > 1
                }
#endif
        }
}
@

@def foreach_vertex_xorder()
foreach_face_xorder() {
    x -= Delta / 2.;
#if dimension > 1  
    y -= Delta / 2.;
#endif
    @
    @define end_foreach_vertex_xorder() } end_foreach_face_xorder()



@def foreach_xorder_boundary()
OMP_PARALLEL() {
int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
Point point = { 0 };
point.level = depth(); point.n = 1 << point.level;
int _k;
OMP(omp for schedule(static))
    for (_k = 0; _k < point.n + 2*GHOSTS; _k++) {
#if dimension == 1
                point.i = _k;
#endif
#if dimension > 1
                point.j = _k;
                for (point.i = 0; point.i < point.n + 2*GHOSTS; point.i++)
                {
#endif
                    POINT_VARIABLES
                        @
                        @def end_foreach_xorder_boundary()
#if dimension > 1
                }
#endif
            }
    }
    @


@def foreach_face_xorder_boundary()
OMP_PARALLEL() {
int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
Point point = { 0 };
point.level = depth(); point.n = 1 << point.level;
int _k;
OMP(omp for schedule(static))
    for (_k = 0; _k <= point.n + 2*GHOSTS; _k++) {
#if dimension == 1
                point.i = _k;
#endif

#if dimension > 1
                point.j = _k;
                for (point.i = 0; point.i <= point.n + 2*GHOSTS+1; point.i++)
                {
#endif
                    POINT_VARIABLES
                        @
                        @def end_foreach_face_xorder_boundary()
#if dimension > 1
                }
#endif
            }
    }
    @

        @def foreach_vertex_xorder_boundary()
        foreach_face_xorder_boundary() {
        x -= Delta / 2.;
#if dimension > 1  
        y -= Delta / 2.;
#endif
        @
        @define end_foreach_vertex_xorder_boundary() } end_foreach_face_xorder_boundary()

    
double amax(double a, double b){
  if (a < b){
    return b;
  } else {
    return a;
  }
}
double dt_start = 0.;
// Define alternative traversal order to match with the x ordering of vts (basilisk uses z ordering by default)

/*
Writes a structured vtk file (.vts) of a single layer.
*/
void output_vts_ascii_single_layer(FILE* fp, scalar* list, int layerid, bool displace_z, bool outputbreaking )
{

    
    int nthreads_ = omp_get_num_threads();
    omp_set_num_threads(1); 
    
        // MULTIGRID
    
 fputs("<?xml version=\"1.0\"?>\n<VTKFile type=\"StructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n", fp);

#if dimension == 1
    fprintf(fp, "\t <StructuredGrid WholeExtent=\"%d %d %d %d 0 0\">\n", 0, N, 0, 0);
    fprintf(fp, "\t\t <Piece Extent=\"%d %d %d %d  0 0\">\n", 0, N, 0, nl);
#endif

#if dimension == 2
    fprintf(fp, "\t <StructuredGrid WholeExtent=\"%d %d %d %d %d %d\">\n", 0, N, 0, N, 0, 0);
    fprintf(fp, "\t\t <Piece Extent=\"%d %d %d %d %d %d\">\n", 0, N, 0, N, 0, nl);
#endif

    // Loop over velocity data and store kinematics in cell vector stucture
    int layertemp = _layer;
    _layer = layerid;
    fputs("\t\t\t <CellData Scalars=\"scalars\">\n", fp);
    fprintf(fp, "\t\t\t\t <DataArray type=\"Float64\" NumberOfComponents=\"3\" Name=\"Velocity\" format=\"ascii\">\n");
    foreach() {
#if dimension == 1
        fprintf(fp, "%g %g 0.\n", u.x[], w[]);
#endif
#if dimension == 2
        fprintf(fp, "%g %g %g\n", u.x[], u.y[], w[]);
#endif
    }
    fputs("\t\t\t\t </DataArray>\n", fp);


    // loop over all scalars in scalarlist and store values as cell data
    for (scalar s in list) {
        if (strcmp(s.name, "eta") == 0) {
            fprintf(fp, "\t\t\t\t <DataArray type=\"Float64\" Name=\"%s\" format=\"ascii\">\n", s.name);
            foreach() {
                if (h[] > dry) {
                    fprintf(fp, "%g\n", val(s));
                }
                else {
                    fprintf(fp, "nan\n");
                }
            }
            fputs("\t\t\t\t </DataArray>\n", fp);
        }
        else {
            fprintf(fp, "\t\t\t\t <DataArray type=\"Float64\" Name=\"%s\" format=\"ascii\">\n", s.name);
            foreach() {
                fprintf(fp, "%g\n", val(s));
            }
            fputs("\t\t\t\t </DataArray>\n", fp);
        }
        
    }

   

    fputs("\t\t\t </CellData>\n", fp);


    // Coordinates 
    fputs("\t\t\t <Points>\n", fp);
#if dimension == 1
    fputs("\t\t\t\t <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n", fp);

    if (displace_z) {
        scalar z_trans = list[0];
        // set z coordinate to the value of the first scalar field in the provided list.
        double hh;
        foreach_vertex(serial) {
            hh = (z_trans[] + z_trans[-1]) / 2.;
            fprintf(fp, "%12.4f %12.4f 0.\n", x, hh);

        }
    }
    else {
        // set z coordinate to 0.
        foreach_vertex() {
            fprintf(fp, "%12.4f 0. 0.\n", x);
        }
    }
#endif
#if dimension == 2
    fputs("\t\t\t\t <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n", fp);
    
    if (displace_z) {
        scalar z_trans = list[0];
        // first do the bottom coordinates stored in zb
        double hh;
        foreach_vertex(serial) {
            hh = (z_trans[] + z_trans[-1] + z_trans[0, -1] + z_trans[-1, -1]) / 4.;
            fprintf(fp, "%12.4f %12.4f %12.4f\n", x, y, hh);
        }
    }
    else {
        // first do the bottom coordinates stored in zb
        foreach_vertex() {
            fprintf(fp, "%12.4f %12.4f 0.0\n", x, y);
        }
    }
#endif

    fputs("\t\t\t\t </DataArray>\n", fp);
    fputs("\t\t\t </Points>\n", fp);
    fputs("\t\t </Piece>\n", fp);

    // write time value
    fprintf(fp, "\t\t <FieldData> \n");
    fprintf(fp, "\t\t\t <DataArray type = \"Float64\" Name = \"%s\" NumberOfTuples = \"1\" format = \"ascii\" RangeMin = \"%.3f\" RangeMax = \"%.3f\"> \n", "TimeValue", t + dt_start, t + dt_start);
    fprintf(fp, "\t\t\t %.3f \n", t + dt_start);
    fprintf(fp, "\t\t\t </DataArray > \n");
    fprintf(fp, "\t\t </FieldData> \n");
    fputs("\t </StructuredGrid>\n", fp);
    fputs("</VTKFile>\n", fp);
    fflush(fp);

    _layer = layertemp;
    omp_set_num_threads(nthreads_);

}


void output_vts_ascii_all_layers(FILE* fp, scalar* list, int Nx, int Ny=0)
{
  int nthreads_ = omp_get_max_threads();
  omp_set_num_threads(1);
    
  // MULTIGRID
    
  fputs("<?xml version=\"1.0\"?>\n<VTKFile type=\"StructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n", fp);

  #if dimension == 1
    fprintf(fp, "\t <StructuredGrid WholeExtent=\"%d %d %d %d 0 0\">\n", 0, Nx, 0, nl);
    fprintf(fp, "\t\t <Piece Extent=\"%d %d %d %d  0 0\">\n", 0, Nx, 0, nl);
  #endif

  #if dimension == 2
    fprintf(fp, "\t <StructuredGrid WholeExtent=\"%d %d %d %d %d %d\">\n", 0, Ny, 0, Nx, 0, nl);
    fprintf(fp, "\t\t <Piece Extent=\"%d %d %d %d %d %d\">\n", 0, Ny, 0, Nx, 0, nl);
  #endif

  // Loop over velocity data and store kinematics in cell vector stucture
  fputs("\t\t\t <CellData Scalars=\"scalars\">\n", fp);
  fprintf(fp, "\t\t\t\t <DataArray type=\"Float64\" NumberOfComponents=\"3\" Name=\"Velocity\" format=\"ascii\">\n");
  for(int i = nl-1; i >= 0; i--){
    foreach() {
      #if dimension == 1
        fprintf(fp, "%g %g 0.\n", u.x[0,0,i], w[0,0,i]);
      #endif
      #if dimension == 2
        if (h[] > dry){
          fprintf(fp, "%g %g %g\n", u.x[0,0,i], u.y[0,0,i], w[0,0,i]);
        } else {
          fprintf(fp, "nan nan nan\n");
        }
      #endif
    }
  }
  // Dummy text to get correct amount of layers
  foreach() {
    if (h[] > dry){
      fprintf(fp, "0 0 0.\n");
    } else {
      fprintf(fp, "nan nan nan\n");
    }
  }
  fputs("\t\t\t\t </DataArray>\n", fp);


  // loop over all scalars in scalarlist and store values as cell data
  for (scalar s in list) {
    if (strcmp(s.name, "eta") == 0) {
      fprintf(fp, "\t\t\t\t <DataArray type=\"Float64\" Name=\"%s\" format=\"ascii\">\n", s.name);
        for(int i = nl-1; i >= 0; i--){
          foreach() {
            if (h[] > dry) {
              fprintf(fp, "%g\n", s[0,0,i]);
            } else {
              fprintf(fp, "nan\n");
            }
          }
        }
      fputs("\t\t\t\t </DataArray>\n", fp);
    } else {
      fprintf(fp, "\t\t\t\t <DataArray type=\"Float64\" Name=\"%s\" format=\"ascii\">\n", s.name);
      for(int i = nl-1; i >= 0; i--){
        foreach() {
          if (h[] > dry){
            fprintf(fp, "%g\n", s[0,0,i]);
          } else {
            fprintf(fp, "nan\n");
          }
        }
      }
      foreach(){
        if (h[] > dry){
          fprintf(fp, "0\n");
        } else {
          fprintf(fp, "nan\n");
        }
      }
      fputs("\t\t\t\t </DataArray>\n", fp);
    }   
  }   

    fputs("\t\t\t </CellData>\n", fp);
    // Coordinates 
    fputs("\t\t\t <Points>\n", fp);
#if dimension == 1
  fputs("\t\t\t\t <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n", fp);
  double* zcorr = (double *)malloc((Nx+1)*sizeof(double));
  // for (int j = 0; j <= N; j++){
  //   zcorr[j] = eta[j,0,0];
  // }
  int mi = 0;
  foreach_vertex (serial){
    zcorr[mi] = eta[0,0,0];
    //fprintf(stdout, "%f \n", eta[0,0,0]); 
    mi++;
  }
  int k;
  for(int i = nl-1; i >= 0; i--){
    k = 0;
    foreach_vertex(serial) {
        fprintf(fp, "%12.4f %12.4f 0.\n", x,zcorr[k]);
        if (h[0,0,k] < 1e-3){
          zcorr[k] = zcorr[k] - h[-1,-1,i];
        } else{ 
          
          // if(x>0.37){
            
          // }
          zcorr[k] = zcorr[k] - h[0,0,i];
        }
        k++;
    }
  }
  k = 0;
  foreach_vertex(serial){
    fprintf(fp, "%12.4f %12.4f 0.\n", x, zcorr[k]);
    k++;
  }
  free(zcorr);
#endif
#if dimension == 2
  fputs("\t\t\t\t <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n", fp);
  double** zcorr = (double**)malloc(nl*sizeof(double*));
  int j = 0;
  zcorr[0] = (double *)malloc((Ny+1)*(Nx+1)*sizeof(double));
  foreach_vertex(serial){
      if (IN_DOMAIN){
        zcorr[0][j] = zb[] + h[0,0,0];
      } else {
        zcorr[0][j] = zb[-1,0] + h[-1,0,0];
      }
      j++;
    }
    j = 0;
  for (int i = 1; i < nl; i++) {
    zcorr[i] = (double *)malloc((Ny+1)*(Nx+1)*sizeof(double));
    foreach_vertex(serial){
      if (IN_DOMAIN){
        zcorr[i][j] = zcorr[i-1][j] + h[0,0,i];
      } else {
        zcorr[i][j] = zb[-1,0] + h[-1,0,i];
      }
      j++;
    }
    j = 0;
  }
  for(int i = nl-1; i >= 0; i--){
    j = 0;
    foreach_vertex(serial) {
      if (IN_DOMAIN){
        fprintf(fp, "%12.4f %12.4f %12.4f\n", x,y, zcorr[i][j]);
        /*
        if (h[0,0,i] < 1e-3){
          double corr = max(h[-1,0,i], h[1,0,i]);
          corr = max(corr, h[0,1,i]);
          corr = max(corr, h[0,-1,i]);
          zcorr[i][j] = zcorr[j] -corr;
        } else {
          zcorr[j] = zcorr[j] - h[0,0,i];
        }
        */
      } else {
        fprintf(fp, "nan nan nan\n");
      }
      j++;
    }
  }
  j = 0;

  foreach_vertex(serial){
    if (IN_DOMAIN){
      fprintf(fp, "%12.4f %12.4f %12.4f\n", x,y, zb[]);
    } else {
      fprintf(fp, "nan nan nan\n");
    }
    j++;
  }
  for (int i = 0; i < nl; i++){
    free(zcorr[i]);
  }
  free(zcorr);
#endif

    fputs("\t\t\t\t </DataArray>\n", fp);
    fputs("\t\t\t </Points>\n", fp);
    fputs("\t\t </Piece>\n", fp);

    // write time value
    fprintf(fp, "\t\t <FieldData> \n");
    fprintf(fp, "\t\t\t <DataArray type = \"Float64\" Name = \"%s\" NumberOfTuples = \"1\" format = \"ascii\" RangeMin = \"%.3f\" RangeMax = \"%.3f\"> \n", "TimeValue", t + dt_start, t + dt_start);
    fprintf(fp, "\t\t\t %.3f \n", t + dt_start);
    fprintf(fp, "\t\t\t </DataArray > \n");
    fprintf(fp, "\t\t </FieldData> \n");
    fputs("\t </StructuredGrid>\n", fp);
    fputs("</VTKFile>\n", fp);
    fflush(fp);
    omp_set_num_threads(nthreads_);

}