//
// An API for OpenCASCADE library
// for performing CAD functionality
// such as talking to IGES and STEP models
// and surface parametrization and interpolation
// used in variety of Finite Element software.
//
//
// to compile and run successfully, first put opencascade shared 
// object library files *.so into LD_LIBRARY_PATH so OS can 
// locate them. For example the following will do it:
//
// $>> export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:<<path to Ocascade \lib>>
//
// then compile and link using:
// 
// $>> g++ -Wall -DTEST ocas_hooks.cxx -I Ocascade_inc/ -L Ocascade_lib/lib*
//
// to get the tester program executable. <ake sure that always 
// have the opencascade header file .hxx in the path. 
// I used -I to include the path for the linker.
//
// written by : A. Ghasemi
// version 1.01
// --------------------------
// Sample Compilation String:
// --------------------------
// g++ -Wall -DTEST ocas_hooks.cxx 
// -I ~/Desktop/fortran_opencascade_wrapper/opencascade-6.9.0/install/inc/ 
// -L ~/Desktop/fortran_opencascade_wrapper/opencascade-6.9.0/install
//                                    /lin64/gcc/lib/lib*
//
#include "IGESControl_Reader.hxx" 
#include "TColStd_HSequenceOfTransient.hxx" 
#include "TopoDS_Shape.hxx" 
#include "Interface_Static.hxx"
#include "TopExp_Explorer.hxx"
#include "TopoDS.hxx"
#include "Geom_Surface.hxx"
#include "BRep_Tool.hxx"
#include "ShapeAnalysis_Surface.hxx"
#include "STEPControl_Reader.hxx"
#include <BRepMesh_IncrementalMesh.hxx>
#include "Poly_Triangulation.hxx"
//#include <omp.h>
#include "GeomAPI_ProjectPointOnSurf.hxx"

//prototypes
int init_bn_boxs(void);
int find_cad_faces_bounding_boxes(void);
int pt_in_box(int ii, const gp_Pnt& tpt);

// macros
#define MY_MIN(X, Y) (((X) < (Y)) ? (X) : (Y))
#define MY_MAX(X, Y) (((X) > (Y)) ? (X) : (Y))
#define TRI_QUERY_FAST
#define USE_SURF_ANALY
#define ECHO_SURF_TRIANG_TO_FILE //also includes bn_boxes

//
// Statically accessible vars and storage
//
TopoDS_Shape sh_static;
TopExp_Explorer anExp_static;
int CURRENT_CAD_FACE_TAG = -1;
gp_Pnt2d UVprev_stored;
int CURRENT_TOPODS_FACES = -1;

// The following is used to store bounding
// box for each TopoDS face of CAD model.
#define MAX_TOPODS_FACES 10000
struct topods_face_bn_box 
{
Standard_Real xmin;
Standard_Real xmax;
Standard_Real ymin;
Standard_Real ymax;
Standard_Real zmin;
Standard_Real zmax;
};
struct topods_face_bn_box bn_boxs[MAX_TOPODS_FACES];
#define BOX_MARGIN_FACTOR 0.1

// reads STEP database and finally initializes
// the static topological classes Shape and Explorer.
// This is only done once to increase performance
// of repeated surface parametrization queries.
//
extern "C" int init_STEP(const char* fname)
{
  // Read the CAD database in STEP format
  STEPControl_Reader mySTEPReader; 
  Standard_Integer nSTEPFaces, nTransFaces; 
  mySTEPReader.ReadFile (fname);

  // info
  //Show the general Information
  mySTEPReader.PrintCheckLoad(Standard_False, IFSelect_GeneralInfo);
  mySTEPReader.PrintCheckLoad(Standard_False, IFSelect_ItemsByEntity); 
  std::cout << "The tolerance defined in STEP file is : " 
       << Interface_Static::RVal("read.precision.val") << endl;

  //selects all STEP faces in the file and puts them into a list,  
  Handle(TColStd_HSequenceOfTransient) myList =  
    mySTEPReader.GiveList("step-faces");
 
  //Transfer myList, 
  nSTEPFaces = myList->Length();  
  nTransFaces = mySTEPReader.TransferList(myList); 
  std::cout << "STEP Faces: " << nSTEPFaces << ", Transferred: " 
	    << nTransFaces << endl; 
  //and obtains the results in an OCCT shape. 
  sh_static = mySTEPReader.OneShape(); 

  anExp_static.Init(sh_static, TopAbs_FACE);

  find_cad_faces_bounding_boxes();

  // done here!
  return 0;
}

// reads IGES database and finally initializes
// the static topological classes Shape and Explorer.
// This is only done once to increase performance
// of repeated surface parametrization queries.
//
extern "C" int init_IGES(const char* fname)
{
  // Read the CAD database in IGES format
  IGESControl_Reader myIgesReader; 
  Standard_Integer nIgesFaces, nTransFaces; 
  myIgesReader.ReadFile (fname);

  // info
  //Show the general Information
  myIgesReader.PrintCheckLoad(Standard_False, IFSelect_GeneralInfo);
  myIgesReader.PrintCheckLoad(Standard_False, IFSelect_ItemsByEntity); 
  std::cout << "The tolerance defined in IGES file is : " 
       << Interface_Static::RVal("read.precision.val") << endl;

  //selects all IGES faces in the file and puts them into a list,  
  Handle(TColStd_HSequenceOfTransient) myList =  
    myIgesReader.GiveList("iges-faces");
 
  //Transfer myList, 
  nIgesFaces = myList->Length();  
  nTransFaces = myIgesReader.TransferList(myList); 
  std::cout << "IGES Faces: " << nIgesFaces << ", Transferred: " 
	    << nTransFaces << endl; 
  //and obtains the results in an OCCT shape. 
  sh_static = myIgesReader.OneShape(); 

  anExp_static.Init(sh_static, TopAbs_FACE);

  find_cad_faces_bounding_boxes();

  // done here!
  return 0;
}

// Cleans the static memory vars 
//
extern "C" int clean_statics(void)
{
  // clean shape objects
  sh_static.Nullify();

  // Cleans up the explorer object
  anExp_static.Clear();

  // done here!
  return 0;
}

#ifndef TRI_QUERY_FAST

// finds the local parametric coords "uv[]" of the given physical
// point "XYZ" on a CAD database in IGES format.
//
// npts = number of points in the query
// *pts = [x1, y1, z1, x2, y2, z2, ....] coords of the points
// *found = [CAD Face Number = Found,
//          -1              = Not Found] 
//          for each point (a vector of npts length)
// *uv = [u1, v1, u2, v2, u3, v3, ....] surface parameters of found points
// tol = the given tolerance to match the CAD face
//
extern "C" int find_pts_on_database(int npts, double *pts
			 , int *found, double *uv, double tol)
{

  // Now, Start querying ...
  gp_Pnt pt_samp, pt_samp2;
  gp_Pnt2d pt_uv;
  double dist, min_dist, min_u, min_v;
  int ii, min_ii;

  for (int indx = 0; indx < npts; indx++)
    {
      printf("indx = %d out of %d, percent completed %3.2f \n"
	     , indx, npts, ((float)indx/(float)npts)*100.);
      // HARD Reset!
      anExp_static.ReInit();

      // take one point from input array
      pt_samp.SetX(pts[indx*3]);
      pt_samp.SetY(pts[indx*3+1]);
      pt_samp.SetZ(pts[indx*3+2]);

      ii = 0; //Now, explore what face in database this point lies on
      min_dist = 1.0e14;

      for(;(anExp_static.More() && (1));anExp_static.Next()){
	ii = ii + 1;

	if ( !pt_in_box(ii, pt_samp) ) continue;

	const TopoDS_Face& anFace = TopoDS::Face(anExp_static.Current());
	// get face as surface
	const Handle(Geom_Surface) &surface = BRep_Tool::Surface(anFace);
	ShapeAnalysis_Surface sas(surface);

	// find parameters uv of that point on the
	// given surface with the given tolerance
	pt_uv = sas.ValueOfUV(pt_samp, tol);
  
	//reevaluate the 3D coordinate of the point using surface parametrization 
	pt_samp2 = sas.Value(pt_uv.X(),pt_uv.Y());
	dist = sqrt( (pt_samp2.X() - pt_samp.X()) * (pt_samp2.X() - pt_samp.X()) 
		   + (pt_samp2.Y() - pt_samp.Y()) * (pt_samp2.Y() - pt_samp.Y()) 
		   + (pt_samp2.Z() - pt_samp.Z()) * (pt_samp2.Z() - pt_samp.Z()));

	// update closest point and its face tag
	if (dist < min_dist)
	  {
	    min_dist = dist;
	    min_ii = ii;
	    min_u = pt_uv.X();
	    min_v = pt_uv.Y();
	  }

	// std::cout << "dist = " << dist << endl;

      } //next face in database

      if ( min_dist <= tol ) // found! 
	{
	  found[indx]  = min_ii;
	  uv[indx*2]   = min_u;
	  uv[indx*2+1] = min_v;
	  // OK, since we found it then break the loop
	  // break;
	  printf("\n foundddd! \n");
	}
      else //NOT found!
	{

	  found[indx] = -1;

	}


    } //next point in the query

  // done here!
  return 0;
}

#else //use the FAST version

// finds the local parametric coords "uv[]" of the given physical
// point "XYZ" on a CAD database in IGES or STEP format.
//
// npts = number of points in the query
// *pts = [x1, y1, z1, x2, y2, z2, ....] coords of the points
// *found = [CAD Face Number = Found,
//          -1              = Not Found] 
//          for each point (a vector of npts length)
// *uv = [u1, v1, u2, v2, u3, v3, ....] surface parameters of found points
// tol = the given tolerance to match the CAD face
//
extern "C" int find_pts_on_database(int npts, double *pts
				    , int *found, double *uv, double tol)
{

  // local vars
  gp_Pnt pt_samp, pt_samp2;
  gp_Pnt2d pt_uv;
  double dist;
  int ii, indx;
  double *min_dist = NULL;

  // init.
  min_dist = (double *)malloc( npts * sizeof(double) );
  for (indx = 0; indx < npts; indx++)
    {
      found[indx] = -1;
      min_dist[indx] = 10. * tol; //1.e14;
    }

  // HARD Reset!
  anExp_static.ReInit();

  ii = 0; //Now, explore and loop over faces
  for(;(anExp_static.More() && (1));anExp_static.Next()){
    ii = ii + 1;

    const TopoDS_Face& anFace = TopoDS::Face(anExp_static.Current());
    // get face as surface
    const Handle(Geom_Surface) &surface = BRep_Tool::Surface(anFace);

#ifdef USE_SURF_ANALY
    ShapeAnalysis_Surface sas(surface);
#else

    GeomAPI_ProjectPointOnSurf projpnt;
#endif

// #pragma omp parallel for shared(found, uv, tol, pts, npts, anExp_static, ii, min_dist) private(pt_samp, pt_samp2, pt_uv, dist) firstprivate(sas) schedule(auto)
    for (indx = 0; indx < npts; indx++)
      {

    // take one point from input array
    pt_samp.SetX(pts[indx*3]);
    pt_samp.SetY(pts[indx*3+1]);
    pt_samp.SetZ(pts[indx*3+2]);

    // skip if out-of-bound for "ii"th CAD face
    if ( !pt_in_box(ii, pt_samp) ) continue;

#ifdef USE_SURF_ANALY
    // find parameters uv of that point on the
    // given surface with the given tolerance
    pt_uv = sas.ValueOfUV(pt_samp, tol);
  
    //reevaluate the 3D coordinate of the point using surface parametrization 
    pt_samp2 = sas.Value(pt_uv.X(),pt_uv.Y());

#else

    // init again the orthogonal projection object
    projpnt.Init(pt_samp, surface);
    if ( projpnt.NbPoints() ) 
      { 
	pt_samp2 = projpnt.NearestPoint();
	double au, av; //the u- and v-coordinates of the projected point
	projpnt.LowerDistanceParameters(au, av); //get the nearest projection
	pt_uv = gp_Pnt2d(au, av); //equivalent 2d description of pnt on surf
      }
    else
      {
	printf("Warning : in GeomAPI_ProjectPointOnSurf class, Init() yields zero NbPoints()!");  
      }

#endif

    dist = sqrt( (pt_samp2.X() - pt_samp.X()) * (pt_samp2.X() - pt_samp.X()) 
      + (pt_samp2.Y() - pt_samp.Y()) * (pt_samp2.Y() - pt_samp.Y()) 
      + (pt_samp2.Z() - pt_samp.Z()) * (pt_samp2.Z() - pt_samp.Z()));

    if ( dist <= min_dist[indx] ) // found! 
      {
	min_dist[indx] = dist;
	found[indx]  = ii;
	uv[indx*2]   = pt_uv.X();
	uv[indx*2+1] = pt_uv.Y();
      }
    else //NOT found!
      {

      }

      } //next node in the list


  } //next CAD face

  //cleanups;
  free(min_dist);

  // done here!
  return 0;
}

#endif

//
// assuming point "A" with parametric coords "uv"
// is somewhere on face#CAD_face, then the following
// will return the *xyz=[x,y,z] of that point
// in the global space.
//
extern "C" int uv2xyz(int CAD_face, double *uv, double *xyz)
{

  gp_Pnt pt_samp;
  anExp_static.ReInit();

  // go to that face sequntially 
  for(int ii = 1; ii < CAD_face; ii++)
    {
      anExp_static.Next();
    }

  const TopoDS_Face& anFace = TopoDS::Face(anExp_static.Current());
  // get face as surface
  const Handle(Geom_Surface) &surface = BRep_Tool::Surface(anFace);
  ShapeAnalysis_Surface sas(surface);

  //evaluate the 3D coordinate of the point using surface parametrization 
  pt_samp = sas.Value(uv[0],uv[1]);

  xyz[0] = pt_samp.X();
  xyz[1] = pt_samp.Y();
  xyz[2] = pt_samp.Z();

  // done here!
  return 0;

}

int export_bn_boxes_tecplot(const char *fname)
{

  int i, j;
  ofstream myFile;

  printf(" writing bounding boxes to %s ...\n", fname); 
  myFile.open (fname);

  //bullet proofing
  if (CURRENT_TOPODS_FACES == -1)
    {
      printf("\n\n\n bounding boxes first need to be computed! stop \n\n\n");
      exit(0);
    }

  //write the header
  myFile << "title = \"bounding boxes here\" " << endl;
  myFile << "variables = \"x\", \"y\", \"z\" " << endl;
  //continue writing header of tecplot file ...
  myFile << "zone n = " << (8*CURRENT_TOPODS_FACES) << " , e = " 
	 << CURRENT_TOPODS_FACES << " , f = fepoint, et = brick " << endl;

  // first write the coords
  for( i = 1; i <= CURRENT_TOPODS_FACES; i++)
    {
      myFile << bn_boxs[i].xmin << " " << bn_boxs[i].ymin << " " << bn_boxs[i].zmin << endl;
      myFile << bn_boxs[i].xmax << " " << bn_boxs[i].ymin << " " << bn_boxs[i].zmin << endl;
      myFile << bn_boxs[i].xmax << " " << bn_boxs[i].ymax << " " << bn_boxs[i].zmin << endl;
      myFile << bn_boxs[i].xmin << " " << bn_boxs[i].ymax << " " << bn_boxs[i].zmin << endl;

      myFile << bn_boxs[i].xmin << " " << bn_boxs[i].ymin << " " << bn_boxs[i].zmax << endl;
      myFile << bn_boxs[i].xmax << " " << bn_boxs[i].ymin << " " << bn_boxs[i].zmax << endl;
      myFile << bn_boxs[i].xmax << " " << bn_boxs[i].ymax << " " << bn_boxs[i].zmax << endl;
      myFile << bn_boxs[i].xmin << " " << bn_boxs[i].ymax << " " << bn_boxs[i].zmax << endl;

    }

  // then, write the connectivity
  for( i = 1; i <= CURRENT_TOPODS_FACES; i++)
    {
      for ( j = 1; j <= 8; j++)
	myFile << (8*(i-1) + j) << " "; 

      myFile << endl;
    }

  //close the outputfile
  myFile.close();
  printf(" done writing to %s! \n", fname); 

  // done here!
  return 0;
}

// finds the bounding box of each TopoDS_Face
// and stores it in static array bn_boxs[]
//  
int find_cad_faces_bounding_boxes(void)
{

  // local vars
  int ii, jj;
  int nCellFace, nNodeFace;
  gp_Pnt vert[3];
  Standard_Integer pts[3];

#ifdef ECHO_SURF_TRIANG_TO_FILE
  // this file output is only for debug
  // and/or visualization
  ofstream myFile;
  // printf(" writing openCASCADE incremental face mesh to opencascade_faces.m \n"); 
  // myFile.open ("opencascade_faces.m");
  // myFile << "tris = [" << endl;
  printf(" writing openCASCADE incremental face mesh to opencascade_faces.tec \n"); 
  myFile.open ("opencascade_faces.tec");
#endif

  // HARD Reset!
  anExp_static.ReInit();

  // init bounding boxes
  init_bn_boxs();

  // mesh the shape to get the bounds from
  // triangulation <NOT VERY ACCURATE>
  BRepMesh_IncrementalMesh(sh_static,100.0);

  ii = 0; //found this!!!

  // loop over faces
  for(;(anExp_static.More() && (1));anExp_static.Next()){
    ii = ii + 1;

    const TopoDS_Face& anFace = TopoDS::Face(anExp_static.Current());

    // Get triangulation
    TopLoc_Location L;
    Handle (Poly_Triangulation) facing = BRep_Tool::Triangulation(anFace,L);
    const Poly_Array1OfTriangle & triangles = facing->Triangles();
    const TColgp_Array1OfPnt & nodes = facing->Nodes();

#ifdef ECHO_SURF_TRIANG_TO_FILE
    nCellFace = triangles.Upper() - triangles.Lower() + 1;
    nNodeFace = nodes.Upper() - nodes.Lower() + 1;

    myFile << "title = \"tris\" " << endl;
    myFile << "variables = \"x\", \"y\", \"z\" " << endl;
    //continue writing header of tecplot file ...
    myFile << "zone n = " << nNodeFace << " , e = " 
           << nCellFace << " , f = fepoint, et = triangle " << endl;

    if (!facing.IsNull()) //if that face is triangulated
      {
	// write face triangulation nodes to file
	for ( int i=1; i <= nNodeFace; i++ )
	  {
	      myFile << nodes(i).X() << " " << nodes(i).Y() 
		     << " " << nodes(i).Z() << endl;
          }

	for ( int i=facing->NbTriangles(); i >= 1; --i )
	  {
	    //get one of the triangles on that face
	    Poly_Triangle triangle = triangles(i);
	    // get the node number of that triangle
	    triangle.Get(pts[0], pts[1], pts[2]);
	    myFile << pts[0] << " " << pts[1] << " " << pts[2] << endl;
          }
       }
    else
       {
         myFile << " No triangulation for this CAD face! " << endl;
       }
#endif

    if (!facing.IsNull()) //if that face is triangulated
      {
	for ( int i=facing->NbTriangles(); i >= 1; --i )
	  {
	    //get one of the triangles on that face
	    Poly_Triangle triangle = triangles(i);
	    // get the node number of that triangle
	    triangle.Get(pts[0], pts[1], pts[2]);

	    for (jj = 0; jj < 3; jj++)
	      {
		vert[jj] = nodes(pts[jj]); //store vertices for visualization

		//
		// Now, update the bounding box of each face
		//
		// xmin
		bn_boxs[ii].xmin = MY_MIN(bn_boxs[ii].xmin, vert[jj].X());
		// xmax
		bn_boxs[ii].xmax = MY_MAX(bn_boxs[ii].xmax, vert[jj].X());
		// ymin
		bn_boxs[ii].ymin = MY_MIN(bn_boxs[ii].ymin, vert[jj].Y());
		// ymax
		bn_boxs[ii].ymax = MY_MAX(bn_boxs[ii].ymax, vert[jj].Y());
		// zmin
		bn_boxs[ii].zmin = MY_MIN(bn_boxs[ii].zmin, vert[jj].Z());
		// zmax
		bn_boxs[ii].zmax = MY_MAX(bn_boxs[ii].zmax, vert[jj].Z());
	      }

// #ifdef ECHO_SURF_TRIANG_TO_FILE
// 	    // write face triangulation to file
// 	    // three vertices
// 	    for (jj = 0; jj < 3; jj++)
// 	      myFile << vert[jj].X() << " " << vert[jj].Y() 
// 		     << " " << vert[jj].Z() << ";" << endl;
// 	    // and the last one repeated
// 	    jj = 0;
// 	    myFile << vert[jj].X() << " " << vert[jj].Y() 
// 		   << " " << vert[jj].Z() << ";" << endl;
// #endif

	  }
      }

  }

  // save the number of CAD faces
  CURRENT_TOPODS_FACES = ii;

#ifdef ECHO_SURF_TRIANG_TO_FILE
  // finalize the debug file
  // myFile << "];" << endl;
  myFile.close();
  // printf(" done writing opencascade_faces.m! \n"); 
  printf(" done writing opencascade_faces.tec! \n"); 

  //export bounding boxes
  export_bn_boxes_tecplot("bn_boxes_tec");

#endif

  // done here!
  return 0;
}

int init_bn_boxs(void)
{

  //local vars
  int i;

  for( i = 0; i < MAX_TOPODS_FACES; i++)
    {
      bn_boxs[i].xmin = 1.0e15;
      bn_boxs[i].xmax = -1.0e15;
      bn_boxs[i].ymin = 1.0e15;
      bn_boxs[i].ymax = -1.0e15;
      bn_boxs[i].zmin = 1.0e15;
      bn_boxs[i].zmax = -1.0e15;
    }

  // done here!
  return 0;

}

// return 1 if this point "tpt" is inside
// or on the box number ii (one-based)
//
int pt_in_box(int ii, const gp_Pnt& tpt)
{
  double dx = BOX_MARGIN_FACTOR * (bn_boxs[ii].xmax - bn_boxs[ii].xmin);
  double dy = BOX_MARGIN_FACTOR * (bn_boxs[ii].ymax - bn_boxs[ii].ymin);
  double dz = BOX_MARGIN_FACTOR * (bn_boxs[ii].zmax - bn_boxs[ii].zmin);

  if ( (tpt.X() >= (bn_boxs[ii].xmin - dx)) && 
       (tpt.X() <= (bn_boxs[ii].xmax + dx)) &&
       (tpt.Y() >= (bn_boxs[ii].ymin - dy)) &&
       (tpt.Y() <= (bn_boxs[ii].ymax + dy)) &&
       (tpt.Z() >= (bn_boxs[ii].zmin - dz)) &&
       (tpt.Z() <= (bn_boxs[ii].zmax + dz)) ) 
    return 1;
  else
    return 0;

  // done here!

}

extern "C" int xyz2uv(int CAD_face, double *xyz, double *uv, double tol)
{

  gp_Pnt pt_samp;
  gp_Pnt2d pt_uv;
  anExp_static.ReInit();

  // go to that face sequntially 
  for(int ii = 1; ii < CAD_face; ii++)
    {
      anExp_static.Next();
    }

  const TopoDS_Face& anFace = TopoDS::Face(anExp_static.Current());
  // get face as surface
  const Handle(Geom_Surface) &surface = BRep_Tool::Surface(anFace);
  ShapeAnalysis_Surface sas(surface);

  // fill the point object
  pt_samp.SetX(xyz[0]);
  pt_samp.SetY(xyz[1]);
  pt_samp.SetZ(xyz[2]);

  // find parameters uv of that point on the
  // given surface with the given tolerance
  pt_uv = sas.ValueOfUV(pt_samp, tol);
  uv[0] = pt_uv.X();
  uv[1] = pt_uv.Y();
  
  // done here!
  return 0;
}

extern "C" int xyz2close_xyz(int CAD_face, double *xyz, double *close_xyz, double tol)
{

  gp_Pnt pt_samp, pt_samp2;
  gp_Pnt2d pt_uv;
  anExp_static.ReInit();

  // go to that face sequntially 
  for(int ii = 1; ii < CAD_face; ii++)
    {
      anExp_static.Next();
    }

  const TopoDS_Face& anFace = TopoDS::Face(anExp_static.Current());
  // get face as surface
  const Handle(Geom_Surface) &surface = BRep_Tool::Surface(anFace);
  ShapeAnalysis_Surface sas(surface);

  // fill the point object
  pt_samp.SetX(xyz[0]);
  pt_samp.SetY(xyz[1]);
  pt_samp.SetZ(xyz[2]);

  if ( CAD_face != CURRENT_CAD_FACE_TAG) //New UV required
    {
      // find parameters uv of that point on the
      // given surface with the given tolerance
      pt_uv = sas.ValueOfUV(pt_samp, tol);

      // store this CAD face for next query to be fast
      CURRENT_CAD_FACE_TAG = CAD_face;
    }
  else
    {
      // pt_uv = sas.NextValueOfUV (UVprev_stored, pt_samp, tol, 1.e-10);
      // if not working, uncomment the following and comment the above
      pt_uv = sas.ValueOfUV(pt_samp, tol);
    }

  // store current UV
  UVprev_stored = pt_uv;

  //reevaluate the 3D coordinate of the point using surface parametrization 
  pt_samp2 = sas.Value(pt_uv.X(),pt_uv.Y());

  // return the closest point (orthogonal projection)
  close_xyz[0] = pt_samp2.X();
  close_xyz[1] = pt_samp2.Y();
  close_xyz[2] = pt_samp2.Z();

  // done here!
  return 0;
}

#ifdef TEST  
int main(int argc, char *argv[])
{

  // local vars
  const char *fname_step = "./crude_samples/store.step";
  const char *fname_iges = "./crude_samples/store.iges";
  // chose a point in 3D space for query
  int npts = 1;
  double *pts = new double[3];
  pts[0]= 0.373333;   pts[1]= 0.433333;   pts[2]= 0.118;
  int *found = new int[1];
  double *uv = new double[2];
  double tol = 1.e-1;

  // Perform a query ...
  // Using IGES CAD model
  init_IGES(fname_iges);
  find_pts_on_database(npts, pts, found, uv, tol);
  find_pts_on_database(npts, pts, found, uv, tol);
  find_pts_on_database(npts, pts, found, uv, tol);
  clean_statics();

  //show output
  std::cout << " ===========  in IGES file: =========== " << endl;
  std::cout << "found[0] = " << found[0] << endl;
  std::cout << "(u, v) = (" << uv[0] << ", " << uv[1] << ") " << endl;

  // Using STEP CAD model
  init_STEP(fname_step);
  find_pts_on_database(npts, pts, found, uv, tol);
  find_pts_on_database(npts, pts, found, uv, tol);
  find_pts_on_database(npts, pts, found, uv, tol);

  //show output
  std::cout << " ===========  in STEP file: =========== " << endl;
  std::cout << "found[0] = " << found[0] << endl;
  std::cout << "(u, v) = (" << uv[0] << ", " << uv[1] << ") " << endl;

  // Now, test snapping function uv2xyz
  pts[0] = 0.; pts[1] = 0.; pts[2] = 0.;
  uv2xyz(found[0], uv, pts);
  std::cout << "snapped point = (" << pts[0] << ", " 
	    << pts[1] << ", " << pts[2] << ")" << endl;  

  // final cleanups
  clean_statics();

  //done here! 
  return 0;
} 

#endif
