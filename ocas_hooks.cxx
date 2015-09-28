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

TopoDS_Shape sh_static;
TopExp_Explorer anExp_static;

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
