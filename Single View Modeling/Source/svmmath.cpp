/***************************************************************
 * CS4670/5670, Fall 2012 Project 4
 * File to be modified #1:
 * svmmath.cpp
 *		a routine for intersecting >2 lines (for vanishing point
 *		computation);
 *		routines for computing the homography for the reference
 *		plane and arbitrary polygons
 **************************************************************/

#pragma warning(disable : 4996)

#include "svmmath.h"
#include "jacob.h"
#include "vec.h"
#include <cstring>
#include <cstdio>
#include <assert.h>
#include <iostream>

#include "Eigen/Core"
#include "MinEig.h"

using namespace Eigen;
using namespace std;

//
// TODO 1: BestFitIntersect()
//		Given lines, the list of 3 or more lines to be intersected,
//		find the best fit intersection point.
//		See http://www-2.cs.cmu.edu/~ph/869/www/notes/vanishing.txt.
//
SVMPoint BestFitIntersect(const std::list<SVMLine> &lines, int imgWidth, int imgHeight)
{
  // check
    if (lines.size() < 2)
    //here should check if the num is below 3
 // if (lines.size() < 3)
	{
      fprintf(stderr, "Not enough lines to compute the best fit.");
      abort();
    }

  SVMPoint bestfit;
  list<SVMLine>::const_iterator iter;

  // To accumulate stuff
  typedef Matrix<double, Dynamic, 3, RowMajor> Matrix3;

  int numLines = (int) lines.size();
  // M to store sum of matrix
  Matrix3 M = Matrix3::Zero(3,3);

  // Deal with numerical conditioning
  int offsW = imgWidth/2;
  int offsH = imgHeight/2;
  int w = (offsH + offsW)/2;

  // temp parameters
  Vec3d temp = Vec3d(0,0,0);
  SVMPoint *p1,*p2;

  int count=0;
  iter = lines.begin();
  while(iter != lines.end()){
	  p1 = iter->pnt1;
	  p2 = iter->pnt2;

	  // numerical conditioning
	  Vec3d nc1 = Vec3d(p1->u-offsW,p1->v-offsH,w);
	  Vec3d nc2 = Vec3d(p2->u-offsW,p2->v-offsH,w);
//	  Vec3d nc1 = Vec3d(p1->u,p1->v,p1->w);
//	  Vec3d nc2 = Vec3d(p2->u,p2->v,p2->w);

	  // cross product
	  Vec3d cp = cross(nc1,nc2);

	  M(0,0) += cp[0] * cp[0];  
	  M(0,1) += cp[0] * cp[1];
	  M(0,2) += cp[0] * cp[2];
	  M(1,1) += cp[1] * cp[1];
	  M(1,2) += cp[1] * cp[2];
	  M(2,2) += cp[2] * cp[2];
	  M(1,0) += cp[1] * cp[0];
	  M(2,0) += cp[2] * cp[0];
	  M(2,1) += cp[2] * cp[1];

	  iter++;
	  count++;
  }

  double eigenVector[3];
  double eigenValue;
  MinEig(M,eigenValue,eigenVector);

  // return the coordinate to 3rd = 1
  bestfit.u = (eigenVector[0]/eigenVector[2])*w + offsW;
  bestfit.v = (eigenVector[1]/eigenVector[2])*w + offsH;
//  bestfit.u = (eigenVector[0]/eigenVector[2]);
//  bestfit.v = (eigenVector[1]/eigenVector[2]);
  bestfit.w = 1;



  return bestfit;
}


//
// TODO 2: ConvertToPlaneCoordinate()
//		Given a plane defined by points, converts their coordinates into
//		a plane coordinate of your choise.
//              See the pdf titled "Homography from Polygon in R^3 to Image Plane",
//              whose link can be found from the project page.
//
//      The final divisors you apply to the u and v coordinates should be saved uScale and vScale
//
void ConvertToPlaneCoordinate(const vector<SVMPoint>& points, vector<Vec3d>& basisPts, double &uScale, double &vScale)
{
  int numPoints = points.size();
  
  // firstly, choose the first and second 
  Vec4d r = Vec4d(points[0].X,points[0].Y,points[0].Z,points[0].W);
  Vec4d p = Vec4d(points[1].X,points[1].Y,points[1].Z,points[1].W);

  // and get the e_x
  Vec4d ex = p-r;
  ex.normalize();

  // then we need to find a point with angle between ex and itself is 90, whose cos is 0
  // and for the worst situation, the angle is 0, they are colinear, cos equals 1
  double bestAngle = 1.0;
  double angle = 1.0;

  // store q
  Vec4d rq;

  for(int count=2;count<numPoints;count++){
	  Vec4d temp = Vec4d(points[count].X,points[count].Y,points[count].Z,points[count].W);
	  Vec4d vTemp = temp - r;

	  // normalize too
	  vTemp.normalize();

	  // the cos is their dot product
	  angle = vTemp * ex;
	  angle = abs(angle);

	  if(angle < bestAngle){
		  bestAngle = angle;
		  rq = vTemp;
	  }
  }

  // obtain e_y
  Vec4d ey = rq - (rq * ex) * ex;
  ey.normalize();

  // pushback and obtain max and min
  double uMax = 0 , uMin = FLT_MAX , vMax = 0 , vMin = FLT_MAX;

  for(int i=0;i<numPoints;i++){
	  Vec4d vec = Vec4d(points[i].X,points[i].Y,points[i].Z,points[i].W);
	  Vec4d a = vec - r;

	  Vec3d basis = Vec3d(a * ex, a * ey , 1);
	  basisPts.push_back(basis);

	  if(basis[0] > uMax) { uMax = basis[0]; }
	  else if (basis[0] < uMin) { uMin = basis[0]; }

	  if(basis[1] > vMax) { vMax = basis[1]; }
	  else if(basis[1] < vMin) { vMin = basis[1]; }
  }

  uScale = uMax - uMin;
  vScale = vMax - vMin;

}



//
// TODO 3: ComputeHomography()
//		Computes the homography H from the plane specified by "points" to the image plane,
//		and its inverse Hinv.
//		If the plane is the reference plane (isRefPlane == true), don't convert the
//		coordinate system to the plane. Only do this for polygon patches where
//		texture mapping is necessary.
//		Coordinate system conversion is to be implemented in a separate routine
//		ConvertToPlaneCoordinate.
//		For more detailed explaination, see the pdf titled
//              "Homography from Polygon in R^3 to Image Plane", whose link can be found from
//              the project page.
//
void ComputeHomography(CTransform3x3 &H, CTransform3x3 &Hinv, const vector<SVMPoint> &points, vector<Vec3d> &basisPts, bool isRefPlane)
{
  int i;
  int numPoints = (int) points.size();
  assert( numPoints >= 4 );

  basisPts.clear();
  if (isRefPlane) // reference plane
    {
      for (i=0; i < numPoints; i++)
        {
          Vec3d tmp = Vec3d(points[i].X, points[i].Y, points[i].W); // was Z, not W
          basisPts.push_back(tmp);
        }
    }
  else // arbitrary polygon
    {
      double uScale, vScale; // unused in this function
      ConvertToPlaneCoordinate(points, basisPts, uScale, vScale);
    }

  // A: 2n x 9 matrix where n is the number of points on the plane
  //    as discussed in lecture
  int numRows = 2 * numPoints;
  const int numCols = 9;

  typedef Matrix<double, Dynamic, 9, RowMajor> MatrixType;
  MatrixType A = MatrixType::Zero(numRows, numCols);

  /******** BEGIN TODO ********/
  for (i = 0; i < numPoints; i++)
	{
		A(2*i,0) = basisPts[i][0];
		A(2*i,1) = basisPts[i][1];
		A(2*i,2) = 1;
		A(2*i,3) = 0;
		A(2*i,4) = 0;
		A(2*i,5) = 0;
		A(2*i,6) = - points[i].u * basisPts[i][0];
		A(2*i,7) = - points[i].u * basisPts[i][1];
		A(2*i,8) = - points[i].u;

		A(2*i+1,0) = 0;
		A(2*i+1,1) = 0;
		A(2*i+1,2) = 0;
		A(2*i+1,3) = basisPts[i][0];
		A(2*i+1,4) = basisPts[i][1];
		A(2*i+1,5) = 1;
		A(2*i+1,6) = - points[i].v * basisPts[i][0];
		A(2*i+1,7) = - points[i].v * basisPts[i][1];
		A(2*i+1,8) = - points[i].v;
	}


  double eval, h[9];
  MinEig(A, eval, h);

  H[0][0] = h[0]; 
  H[0][1] = h[1]; 
  H[0][2] = h[2]; 

  H[1][0] = h[3];
  H[1][1] = h[4]; 
  H[1][2] = h[5];
  
  H[2][0] = h[6];
  H[2][1] = h[7];
  H[2][2] = h[8];


  /******** END TODO ********/

  // compute inverse of H
  if (H.Determinant() == 0)
    fl_alert("Computed homography matrix is uninvertible \n");
  else
    Hinv = H.Inverse();

  int ii;
  printf("\nH=[\n");
  for (ii=0; ii<3; ii++)
    printf("%e\t%e\t%e;\n", H[ii][0]/H[2][2], H[ii][1]/H[2][2], H[ii][2]/H[2][2]);
  printf("]\nHinv=[\n");

  for (ii=0; ii<3; ii++)
    printf("%e\t%e\t%e;\n", Hinv[ii][0]/Hinv[2][2], Hinv[ii][1]/Hinv[2][2], Hinv[ii][2]/Hinv[2][2]);
  printf("]\n\n");
}
