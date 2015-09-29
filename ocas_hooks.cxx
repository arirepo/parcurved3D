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

//prototypes
int init_bn_boxs(void);
int find_cad_faces_bounding_boxes(void);

// macros
#define MY_MIN(X, Y) (((X) < (Y)) ? (X) : (Y))
#define MY_MAX(X, Y) (((X) > (Y)) ? (X) : (Y))

//
// Statically accessible vars and storage
//
TopoDS_Shape sh_static;
TopExp_Explorer anExp_static;

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

// finds the bounding box of each TopoDS_Face
// and stores it in static array bn_boxs[]
//  
int find_cad_faces_bounding_boxes(void)
{

  // local vars
  int ii, jj;
  gp_Pnt vert[3];
  Standard_Integer pts[3];

  // this file output is only for debug
  // and/or visualization
  ofstream myFile;
  printf(" writing openCASCADE incremental face mesh to opencascade_faces.m \n"); 
  myFile.open ("opencascade_faces.m");
  myFile << "tris = [" << endl;

  // HARD Reset!
  anExp_static.ReInit();

  // init bounding boxes
  init_bn_boxs();

  // mesh the shape to get the bounds from
  // triangulation <NOT VERY ACCURATE>
  BRepMesh_IncrementalMesh(sh_static,100.0);

  // loop over faces
  for(;(anExp_static.More() && (1));anExp_static.Next()){
    ii = ii + 1;
    const TopoDS_Face& anFace = TopoDS::Face(anExp_static.Current());

    // Get triangulation
    TopLoc_Location L;
    Handle (Poly_Triangulation) facing = BRep_Tool::Triangulation(anFace,L);
    const Poly_Array1OfTriangle & triangles = facing->Triangles();
    const TColgp_Array1OfPnt & nodes = facing->Nodes();
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

	    // write face triangulation to file
	    // three vertices
	    for (jj = 0; jj < 3; jj++)
	      myFile << vert[jj].X() << " " << vert[jj].Y() 
		     << " " << vert[jj].Z() << ";" << endl;
	    // and the last one repeated
	    jj = 0;
	    myFile << vert[jj].X() << " " << vert[jj].Y() 
		   << " " << vert[jj].Z() << ";" << endl;

	  }
      }


  }

  // finalize the debug file
  myFile << "];" << endl;
  myFile.close();

  printf(" done writing opencascade_faces.m! \n"); 

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
