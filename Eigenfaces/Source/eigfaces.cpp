
/////////////////////////////////////////////////////////////////////////////////////////////////
//	Project 4: Eigenfaces                                                                      //
//  CSE 455 Winter 2003                                                                        //
//	Copyright (c) 2003 University of Washington Department of Computer Science and Engineering //
//                                                                                             //
//  File: eigfaces.cpp                                                                         //
//	Author: David Laurence Dewey                                                               //
//	Contact: ddewey@cs.washington.edu                                                          //
//           http://www.cs.washington.edu/homes/ddewey/                                        //
//                                                                                             //
/////////////////////////////////////////////////////////////////////////////////////////////////



#include "stdafx.h"


EigFaces::EigFaces()
:
Faces()
{
	//empty
}

EigFaces::EigFaces(int count, int width, int height)
:
Faces(count, width, height)
{
	//empty
}

void EigFaces::projectFace(const Face& face, Vector& coefficients) const
{
	if (face.getWidth()!=width || face.getHeight()!=height) {
		throw Error("Project: Face to project has different dimensions");
	}

	coefficients.resize(getSize());
	// ----------- TODO #2: compute the coefficients for the face and store in coefficients.

	int eigenCount = getSize();
	double coe = 0.0;
	for (int count=0;count<eigenCount;count++){
		coe = 0;
		for(int y=0;y<height;y++){
			for(int x=0;x<width;x++){
				coe += (*this)[count].pixel(x,y,0) * (face.pixel(x,y,0)-average_face.pixel(x,y,0));
			}
		}
		coefficients[count] = coe;
	}
}

void EigFaces::constructFace(const Vector& coefficients, Face& result) const
{	
	// ----------- TODO #3: construct a face given the coefficients
	int eigenCount = getSize();

   	for(int y=0;y<height;y++){
	   	for(int x=0;x<width;x++){
			result.pixel(x,y,0) = average_face.pixel(x,y,0);
			for(int count=0;count<eigenCount;count++){
				result.pixel(x,y,0) += (*this)[count].pixel(x,y,0) * coefficients[count];
			}
		}
	}
}

bool EigFaces::isFace(const Face& face, double max_reconstructed_mse, double& mse) const
{
	// ----------- TODO #4: Determine if an image is a face and return true if it is. Return the actual
	// MSE you calculated for the determination in mse
	// Be sure to test this method out with some face images and some non face images
	// to verify it is working correctly.
	Face result(width,height);
	Vector coefficients;

	double accumulatedError = 0;
	mse = 0;
	int n = width*height;

	projectFace(face,coefficients);
	constructFace(coefficients,result);

	for(int y=0;y<height;y++){
		for(int x=0;x<width;x++){
			accumulatedError = face.pixel(x,y,0) - result.pixel(x,y,0);
			mse += accumulatedError*accumulatedError;
		}
	}
	mse /= n;

	if(mse<max_reconstructed_mse){
		return true;
	}
	else{
		return false;
	}
}

bool EigFaces::verifyFace(const Face& face, const Vector& user_coefficients, double max_coefficients_mse, double& mse) const
{
	// ----------- TODO #5 : Determine if face is the same user give the user's coefficients.
	// return the MSE you calculated for the determination in mse.
	Vector coefficients;
	mse = 0;
	int eigenSize = getSize();

	projectFace(face,coefficients);

	for(int i=0;i<eigenSize;i++){
		mse += (user_coefficients[i]-coefficients[i])*(user_coefficients[i]-coefficients[i]);
	}

	mse /= eigenSize;
	if(mse<=max_coefficients_mse){
		return true;
	}
	else{
		return false;
	}
}

void EigFaces::recognizeFace(const Face& face, Users& users) const
{
	// ----------- TODO #6: Sort the users by closeness of match to the face
	Vector coefficients;
	double mse = 0;
	int eigenSize = getSize();

	projectFace(face,coefficients);

	for(int count=0;count<users.getSize();count++){
		mse = 0;
		for(int i=0;i<eigenSize;i++){
			mse += (users[count][i]-coefficients[i])*(users[count][i]-coefficients[i]);
    	}
		mse /= eigenSize;
		users[count].setMse(mse);
	}

	users.sort();

}

void EigFaces::findFace(const Image& img, double min_scale, double max_scale, double step, int n, bool crop, Image& result) const
{
	// ----------- TODO #7: Find the faces in Image. Search image scales from min_scale to max_scale inclusive,
	// stepping by step in between. Find the best n faces that do not overlap each other. If crop is true,
	// n is one and you should return the cropped original img in result. The result must be identical
	// to the original besides being cropped. It cannot be scaled and it must be full color. If crop is
	// false, draw green boxes (use r=100, g=255, b=100) around the n faces found. The result must be
	// identical to the original image except for the addition of the boxes.
	
	Image scaledImage;
	Face subOriginalFace;
	Face subFace;

	int imgWidth = img.getWidth();
	int imgHeight = img.getHeight();

	int cropWidth = width;
	int cropheight = 0;

	double localMse = 0;
	double xdistance;
	double ydistance;
	double distance;
	double faceToface = sqrt(width*width+height*height);

	std::list<FacePosition> facepos;
	std::list<FacePosition>::iterator it;

	// the outmost loop to scale the image and check the face
	for(double scale = min_scale;scale <= max_scale;scale += step){
		std::cout<< "Analysis about [" << scale << "] of image engaged." << std::endl;

		scaledImage.resize((int)(imgWidth*scale),(int)(imgHeight*scale));

		img.resample(scaledImage);

		// ATTENTION:: maybe just loop till x=getWidth-width, y= getHeight-height??
		for(int y=0;y<(scaledImage.getHeight()-cropWidth);y++){
			for(int x=0;x<(scaledImage.getWidth()-cropheight);x++){
		
				cropheight = (int) cropWidth*10/7;

				subOriginalFace.resize(cropWidth,cropheight);
				subOriginalFace.subimage(x,x+cropWidth-1,y,y+cropheight-1,scaledImage,false);
				subFace.resize(width,height);
				subOriginalFace.resample(subFace);

				// the algorithm is always fooled by low variance area, so enhanced them
				bool isface = isFace(subFace,1000.0,localMse);
				if(isface){
					
					FacePosition localPos = FacePosition();
					localPos.x = x;
					localPos.y = y;
					localPos.scale = scale;
					localPos.error = localMse;

					it = facepos.begin();
					if(it == facepos.end()){
						//empty? push back!
						facepos.push_back(localPos);
					}
					else{
						//not empty? suchs we need to go through the list
						//if there is a point around this one and similar,decide to use which
						bool overlap = false;
						for(it=facepos.begin();it!=facepos.end();it++){
							xdistance = ((*it).x/(*it).scale) - (localPos.x/localPos.scale);
							ydistance = ((*it).y/(*it).scale) - (localPos.y/localPos.scale);
							xdistance *= xdistance;
							ydistance *= ydistance;
							//suppose that the face is circular, so use SSD to determine whether overlap
							distance = sqrt(xdistance+ydistance);
							if(distance <= faceToface){
								overlap = true;
								if((*it).error>localPos.error){
									(*it).x = localPos.x;
									(*it).y = localPos.y;
									(*it).scale = localPos.scale;
									(*it).error = localPos.error;

									facepos.sort();
								}
								break;
							}
						}
						//not anyone overlay?? push back!
						if(!overlap){
							facepos.push_back(localPos);
							facepos.sort();
						}
					}
				}//end of the dealing with subface like a face
				
			}
		}//end of loop for scaled image 

		
	}//end of a loop for specific scale

	//now we only need n face candidate
	facepos.sort();
	while(facepos.size()>n){
		facepos.pop_back();
	}
	// now we check what's in the list
	for(it=facepos.begin();it!=facepos.end();it++){
		std::cout << "subFace from["<<(int)((*it).x/(*it).scale)<<","<<(int)((*it).y/(*it).scale)<<"] has mse: " << (*it).error << std::endl;
	}

	if(!crop){
		result.resize(img.getWidth(),img.getHeight(),img.getColors());
		img.resample(result);

		it = facepos.begin();

		int lux,luy,rlx,rly;

		for(it=facepos.begin();it!=facepos.end();it++){
			lux = (int)((*it).x/(*it).scale);
			luy = (int)((*it).y/(*it).scale);
			rlx = (int)(((*it).x+cropWidth)/((*it).scale));
			rly = (int)(((*it).y+cropheight)/((*it).scale));

			result.line(lux,luy,rlx,luy,100,256,100);
			result.line(rlx,luy,rlx,rly,100,256,100);
			result.line(lux,luy,lux,rly,100,256,100);
			result.line(lux,rly,rlx,rly,100,256,100);
		}
	}
	else{
		FacePosition win = facepos.front();
		img.crop(win.x/win.scale,
			win.y/win.scale,
			(win.x+cropWidth)/win.scale,
			(win.y+cropheight)/win.scale,result);
	}
}

void EigFaces::morphFaces(const Face& face1, const Face& face2, double distance, Face& result) const
{
	// TODO (extra credit): MORPH along *distance* fraction of the vector from face1 to face2 by
	// interpolating between the coefficients for the two faces and reconstructing the result.
	// For example, distance 0.0 will approximate the first, while distance 1.0 will approximate the second.
	// Negative distances are ok two.

	Vector coe1;
	Vector coe2;
	Vector coe;

	int eigenCount = getSize();
	coe.resize(eigenCount);
	coe1.resize(eigenCount);
	coe2.resize(eigenCount);

	projectFace(face1,coe1);
	projectFace(face2,coe2);

	// I don't understand the negative,so I assume the input distance is between 0~1

	for(int count=0;count<eigenCount;count++){
		coe[count] = coe1[count]*distance + coe2[count]*(1-distance);
	}

	result.resize(width,height);
	constructFace(coe,result);

}

const Face& EigFaces::getAverage() const
{
	return average_face;
}

void EigFaces::setAverage(const Face& average)
{
	average_face=average;
}



