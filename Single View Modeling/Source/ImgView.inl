/***************************************************************
 * CS4670/5670, Fall 2012 Project 4
 * File to be modified #2:
 * ImgView.inl (included from ImgView.cpp)
 *		contains routines for computing the 3D position of points
 ***************************************************************/

//
// TODO 4: sameXY()
//		Computes the 3D position of newPoint using knownPoint
//		that has the same X and Y coordinate, i.e. is directly
//		below or above newPoint.
//		See lecture slide on measuring heights.
//
// HINT1: make sure to dehomogenize points when necessary
// HINT2: there is a degeneracy that you should look out for involving points already in line with the reference
// HINT3: make sure to get the sign of the result right, i.e. whether it is above or below ground
void ImgView::sameXY()
{
  if (pntSelStack.size() < 2)
    {
      fl_alert("Not enough points on the stack.");
      return;
    }
  
  SVMPoint &newPoint = *pntSelStack[pntSelStack.size() - 1];
  SVMPoint &knownPoint = *pntSelStack[pntSelStack.size() - 2];
  
  if( !knownPoint.known() )
    {
      fl_alert("Can't compute relative values for unknown point.");
      return;
    }
  
  if( refPointOffPlane == NULL )
    {
      fl_alert("Need to specify the reference height first.");
      return;
    }
  

  /******** BEGIN TODO ********/
 
  // for now, three points known, ref,new,basis
  // want to calculate the height, need to find intersection with horizon
  Vec3d vx = Vec3d(xVanish.u,xVanish.v,xVanish.w);
  Vec3d vy = Vec3d(yVanish.u,yVanish.v,yVanish.w);
  Vec3d vz = Vec3d(zVanish.u,zVanish.v,zVanish.w);

  // in order to get the horizon, we cross product the two vanishing points
  // the cross product means the normal and also the line formed by two points
  Vec3d horizon = cross(vx,vy);

  // now we need to find t, according to the known point and homogeneous T for ref
  Mat3d matH = Mat3d(H[0][0],H[0][1],H[0][2],H[1][0],H[1][1],H[1][2],H[2][0],H[2][1],H[2][2]);
  // the transformation has two reasons, first,those points are unknown, second, known point may not on the ground
  // just for H T, so simply set to this, neither a 2d pos nor 3d
  Vec3d pos3dOfBottom0 = Vec3d(knownPoint.X,knownPoint.Y,1);
  Vec3d Pos3dOfBottom = Vec3d(refPointOffPlane->X,refPointOffPlane->Y,1);
  
  // reminder,the left product of mat has already been modified in vec.h, if want to use another, please go and do
  Vec3d b0 = matH*pos3dOfBottom0;
  Vec3d b = matH*Pos3dOfBottom;

  // find the vanishing point in horizon
  Vec3d bb0 = cross(b,b0);
  Vec3d v = cross(bb0,horizon);

  // find t
  Vec3d nul = Vec3d(0,0,0);
  Vec3d t;
  Vec3d t0 = Vec3d(newPoint.u/newPoint.w,newPoint.v/newPoint.w,1);
  Vec3d r = Vec3d(refPointOffPlane->u/refPointOffPlane->w,refPointOffPlane->v/refPointOffPlane->w,1);
  if (v == nul){
	  // if the ref in the same xy of bottom
	  // because we finally use the vertical line through ref to calculate H and R, need to put in there
	  t = t0 - b0 + b;
  }
  else {
	  // if not same
	  Vec3d vt0 = cross(v,t0);
	  Vec3d rb = cross(r,b);
	  t = cross(vt0,rb);
  }

  // here we come to the cross ratio
  // we need to use t,b,vz,r    before using them, normalize first, till their w to 1
  t /= t[2];
  b /= b[2];
  r /= r[2];
  vz /= vz[2];

  double crossRatio = ((t-b).length() * (vz-r).length()) / ((r-b).length() * (vz-t).length());
  double height = crossRatio * referenceHeight;

  // now we get height, but I don't know it should be negative or positive
  // I use the cos angle

  // this operation is because the result from homogeneous T has very little value, which refers to the low value of w
  b0 /= b0[2];
  b /= b[2];

  // here you can also choose r to replace vz, and make sure r is positive in z
  // and notice that here use the angle ,cos. > or < 0 etc
  if( ((t0-b0) * (r-b)) < 0 ) {
	  height *= -1;
  }

  
  
  newPoint.X = knownPoint.X;
  newPoint.Y = knownPoint.Y;
  newPoint.Z = height;
  newPoint.W = 1;
  /******** END TODO ********/
 
  newPoint.known(true);
  
  printf( "Calculated new coordinates for point: (%e, %e, %e)\n", newPoint.X, newPoint.Y, newPoint.Z );
  
  redraw();
}



//
// TODO 5: sameZPlane()
//		Compute the 3D position of newPoint using knownPoint
//		that lies on the same plane and whose 3D position is known.
//		See the man on the box lecture slide.
//		If newPoint is on the reference plane (Z==0), use homography (this->H, or simply H) directly.
//
// HINT: For this function, you will only need to use the three vanishing points and the reference homography 
//       (in addition to the known 3D location of knownPoint, and the 2D location of newPoint)
void ImgView::sameZPlane()
{
  if (pntSelStack.size() < 2)
    {
      fl_alert("Not enough points on the stack.");
      return;
    }
  
  SVMPoint &newPoint = *pntSelStack[pntSelStack.size() - 1];
  SVMPoint &knownPoint = *pntSelStack[pntSelStack.size() - 2];
  
  if( !knownPoint.known() )
    {
      fl_alert("Can't compute relative values for unknown point.");
      return;
    }
  
  /******** BEGIN TODO ********/
  
  // store the inv and normal H T in mat
  Vec3d b0 = Vec3d(newPoint.u,newPoint.v,newPoint.w);
  Mat3d matH = Mat3d(H[0][0],H[0][1],H[0][2],H[1][0],H[1][1],H[1][2],H[2][0],H[2][1],H[2][2]);
  Mat3d invMatH=Mat3d(Hinv[0][0],Hinv[0][1],Hinv[0][2],Hinv[1][0],Hinv[1][1],Hinv[1][2],Hinv[2][0],Hinv[2][1],Hinv[2][2]);


  if (knownPoint.Z != 0){
	  Vec3d t1 = Vec3d(knownPoint.u,knownPoint.v,knownPoint.w);
	  Vec3d m0 = Vec3d(newPoint.u,newPoint.v,newPoint.w);

	  Vec3d vx = Vec3d(xVanish.u,xVanish.v,xVanish.w);
	  Vec3d vy = Vec3d(yVanish.u,yVanish.v,yVanish.w);
	  Vec3d vz = Vec3d(zVanish.u,zVanish.v,zVanish.w);

	  // according to the hint, using vz instead of t0 to find b0
	  Vec3d b0t0 = cross (m0,vz);
	  Vec3d horizon = cross(vx,vy);
	  Vec3d t1m0 = cross(t1,m0);
	  Vec3d v = cross(t1m0,horizon);

	  Vec3d pos3dOfb1 = Vec3d(knownPoint.X,knownPoint.Y,1);
	  Vec3d b1 = matH*pos3dOfb1; 

	  Vec3d nul = Vec3d(0,0,0);
	  if(v == nul){
		  b0 = cross(b1,b0t0);
	  }
	  else {
		  Vec3d b1v = cross(b1,v);
		  b0 = cross(b1v,b0t0);
	  }
  }

//  b0 /= b0[2];
  b0 = invMatH*b0;
  b0 /= b0[2];

  newPoint.X = b0[0];
  newPoint.Y = b0[1];
  newPoint.Z = knownPoint.Z;
  newPoint.W = 1;


  /******** END TODO ********/

  newPoint.known(true);
 
  printf( "Calculated new coordinates for point: (%e, %e, %e)\n", newPoint.X, newPoint.Y, newPoint.Z );
 
  redraw();
}

