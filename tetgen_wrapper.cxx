#include "tetgen.h"

extern "C" int main_tetgen_wrapper(char *command, int npts, double *x, int nquad
			, int ntri, int *icontag, int nhole, double *xh

				   , int *nptsf, double *xf, int *ntet
				   , int *tetcon, int *neigh, int *nbntri, int *bntri) 
{

  //local vars
  //
  // declare and construct (init) io
  tetgenio io, out;
  // Define a pointer of facet.
  tetgenio::facet *f = NULL;
  // Define a pointer of polygon.
  tetgenio::polygon *p = NULL; 
  int indx;

  // All indices start from 1.
  io.firstnumber = 1;
  io.mesh_dim = 3;   // must be 3.

  // fill points
  io.numberofpoints = npts;
  io.pointlist = new REAL[ io.mesh_dim * io.numberofpoints];
  for(int i = 0; i < io.numberofpoints; i++)
    for(int j = 0; j < 3; j++)
      io.pointlist[3*i+j] = x[3*i+j];

  io.numberofpointattributes = 0;
  io.pointattributelist = NULL;
  // io.numberofaddpoints = 0;
  // io.addpointlist = NULL;
  io.pointmarkerlist = NULL;


  // allocating the faces
  io.numberoffacets = nquad+ntri;
  io.facetlist = new tetgenio::facet[io.numberoffacets];
  io.facetmarkerlist = new int[io.numberoffacets];
  indx = 0;

  // add quads if any
  for(int i = 0; i < nquad; i++)
    {

      f = &io.facetlist[i];
      f->numberofpolygons = 1;
      f->polygonlist = new tetgenio::polygon[1];
      f->numberofholes = 0;
      f->holelist = NULL;
      p = &f->polygonlist[0];
      p->numberofvertices = 4;
      // Allocate memory for vertices.
      p->vertexlist = new int[4];
      // fill vertices and tags
      for(int j =0; j < 4; j++)
	{
	  p->vertexlist[j] = icontag[indx];
	  indx++;
	}

      // grab the tag and store
      io.facetmarkerlist[i] = icontag[indx];
      indx++;
    }

  // add triangles
  for(int i = nquad; i < io.numberoffacets; i++)
    {

      f = &io.facetlist[i];
      f->numberofpolygons = 1;
      f->polygonlist = new tetgenio::polygon[1];
      f->numberofholes = 0;
      f->holelist = NULL;
      p = &f->polygonlist[0];
      p->numberofvertices = 3;
      // Allocate memory for vertices.
      p->vertexlist = new int[3];
      // fill vertices and tags
      for(int j =0; j < 3; j++)
	{
	  p->vertexlist[j] = icontag[indx];
	  indx++;
	}

      // grab the tag and store
      io.facetmarkerlist[i] = icontag[indx];
      indx++;
    }


  // adding holes
  io.numberofholes = nhole;
  io.holelist = NULL;
  if ( io.numberofholes > 0)
    {

      io.holelist = new REAL[io.mesh_dim * io.numberofholes];

      for(int i = 0; i < io.numberofholes; i++)
	for(int j = 0; j < 3; j++)
	  io.holelist[3*i+j] = xh[3*i+j];
    }

  io.numberofregions = 0;
  io.regionlist = NULL;

  io.numberoffacetconstraints = 0;
  io.facetconstraintlist = NULL;

  io.numberofsegmentconstraints = 0;
  io.segmentconstraintlist = NULL;

  io.numberoftrifaces = 0;
  io.trifacelist = NULL;
  io.trifacemarkerlist = NULL;

  io.numberofedges = 0;
  io.edgelist = NULL;
  io.edgemarkerlist = NULL;


  // io.save_nodes("alltri_in");
  // io.save_poly("alltri_in");

  // perform tetmeshing ...
  tetrahedralize(command, &io, &out);

  // out.save_nodes("dumped");
  // out.save_elements("dumped");
  // out.save_faces("dumped");

  // preparing the outputs
  //
  //

  *nptsf = out.numberofpoints;

  // bullet proofing ...
  if ( (*nptsf) <= 0 )
    {
      printf("\n no output nodes generated in tetmeshing! stop \n");
      exit(0);
    }
  // copy output nodes (initial + Steiner nodes)
  for (int i = 0; i < (*nptsf); i++)
    for(int j = 0; j < 3; j++)
      xf[3*i+j] = out.pointlist[3*i+j];

  // collect generated tets
  *ntet = out.numberoftetrahedra;
  // bullet proofing ...
  if ( (*ntet) <= 0 )
    {
      printf("\n no output tetrahedrals generated in tetmeshing! stop \n");
      exit(0);
    }

  // gather tet vertices number and neighbors
  for (int i = 0; i < (*ntet); i++)
    for (int j =0; j < 4; j++)
      {
	tetcon[4*i+j] = out.tetrahedronlist[4*i+j];
	neigh[4*i+j] = out.neighborlist[4*i+j]; 
      }

  // gather boundary triangle elements info
  *nbntri = out.numberoftrifaces;
  // bullet proofing ...
  if ( (*nbntri) <= 0 )
    {
      printf("\n no boundary triangles generated in tetmeshing! stop \n");
      exit(0);
    }

  // pack boundary triangles info
  // bntri = [ vertex1, vertex2, vertex3, face marker, adjacent tet1, adjacent tet2]
  for (int i = 0; i < (*nbntri); i++)
    {
      // vertices
      for (int j =0; j < 3; j++)
	bntri[6*i+j] = out.trifacelist[3*i+j];
      // // face marker
      // for (int j =3; j < 4; j++)
      // 	bntri[6*i+j] = out.trifacemarkerlist[i];
      // // adjacent tets 1, 2
      // for (int j =4; j < 6; j++)
      // 	bntri[6*i+j] = out.adjtetlist[2*i+j-4];
    }


  //done here
  return 0;
}
