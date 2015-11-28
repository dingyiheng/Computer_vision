#include <assert.h>
#include <math.h>
#include <cmath>
#include <FL/Fl.H>
#include <FL/Fl_Image.H>
#include "features.h"
#include "ImageLib/FileIO.h"

#define PI 3.14159265358979323846

// Compute features of an image.
bool computeFeatures(CFloatImage &image, FeatureSet &features, int featureType) {
	// TODO: Instead of calling dummyComputeFeatures, write your own
	// feature computation routines and call them here.
	switch (featureType) {
	case 1:
		dummyComputeFeatures(image, features);
		break;
	case 2:
		ComputeHarrisFeatures(image, features);
		break;
	default:
		return false;
	}

	// This is just to make sure the IDs are assigned in order, because
	// the ID gets used to index into the feature array.
	for (unsigned int i=0; i<features.size(); i++) {
		features[i].id = i+1;
	}

	return true;
}

// Perform a query on the database.  This simply runs matchFeatures on
// each image in the database, and returns the feature set of the best
// matching image.
bool performQuery(const FeatureSet &f, const ImageDatabase &db, int &bestIndex, vector<FeatureMatch> &bestMatches, double &bestScore, int matchType) {
	// Here's a nice low number.
	bestScore = -1e100;

	vector<FeatureMatch> tempMatches;
	double tempScore;

	for (unsigned int i=0; i<db.size(); i++) {
		if (!matchFeatures(f, db[i].features, tempMatches, tempScore, matchType)) {
			return false;
		}

		if (tempScore > bestScore) {
			bestIndex = i;
			bestScore = tempScore;
			bestMatches = tempMatches;
		}
	}

	return true;
}

// Match one feature set with another.
bool matchFeatures(const FeatureSet &f1, const FeatureSet &f2, vector<FeatureMatch> &matches, double &totalScore, int matchType) {
	// TODO: We have given you the ssd matching function, you must write your own
	// feature matching function for the ratio test.
	
	printf("\nMatching features.......\n");

	switch (matchType) {
	case 1:
		//case 1,2 for match 5*5 window , normal and simplest descriptor
		ssdMatchFeatures(f1, f2, matches, totalScore);
		return true;
	case 2:
		ratioMatchFeatures(f1, f2, matches, totalScore);
		return true;
	default:
		return false;
	}
}

// Evaluate a match using a ground truth homography.  This computes the
// average SSD distance between the matched feature points and
// the actual transformed positions.
double evaluateMatch(const FeatureSet &f1, const FeatureSet &f2, const vector<FeatureMatch> &matches, double h[9]) {
	double d = 0;
	int n = 0;

	double xNew;
	double yNew;

    unsigned int num_matches = matches.size();
	for (unsigned int i=0; i<num_matches; i++) {
		int id1 = matches[i].id1;
        int id2 = matches[i].id2;
        applyHomography(f1[id1-1].x, f1[id1-1].y, xNew, yNew, h);
		d += sqrt(pow(xNew-f2[id2-1].x,2)+pow(yNew-f2[id2-1].y,2));
		n++;
	}	

	return d / n;
}

void addRocData(const FeatureSet &f1, const FeatureSet &f2, const vector<FeatureMatch> &matches, double h[9],vector<bool> &isMatch,double threshold,double &maxD) {
	double d = 0;

	double xNew;
	double yNew;

    unsigned int num_matches = matches.size();
	for (unsigned int i=0; i<num_matches; i++) {
		int id1 = matches[i].id1;
        int id2 = matches[i].id2;
		applyHomography(f1[id1-1].x, f1[id1-1].y, xNew, yNew, h);

		// Ignore unmatched points.  There might be a better way to
		// handle this.
		d = sqrt(pow(xNew-f2[id2-1].x,2)+pow(yNew-f2[id2-1].y,2));
		if (d<=threshold)
		{
			isMatch.push_back(1);
		}
		else
		{
			isMatch.push_back(0);
		}

		if (matches[i].score>maxD)
			maxD=matches[i].score;
	}	
}

vector<ROCPoint> computeRocCurve(vector<FeatureMatch> &matches,vector<bool> &isMatch,vector<double> &thresholds)
{
	vector<ROCPoint> dataPoints;

	for (int i=0; i < (int)thresholds.size();i++)
	{
		//printf("Checking threshold: %lf.\r\n",thresholds[i]);
		int tp=0;
		int actualCorrect=0;
		int fp=0;
		int actualError=0;
		int total=0;

        int num_matches = (int) matches.size();
		for (int j=0;j < num_matches;j++)
		{
			if (isMatch[j])
			{
				actualCorrect++;
				if (matches[j].score<thresholds[i])
				{
					tp++;
				}
			}
			else
			{
				actualError++;
				if (matches[j].score<thresholds[i])
				{
					fp++;
				}
            }
			
			total++;
		}

		ROCPoint newPoint;
		//printf("newPoints: %lf,%lf",newPoint.trueRate,newPoint.falseRate);
		newPoint.trueRate=(double(tp)/actualCorrect);
		newPoint.falseRate=(double(fp)/actualError);
		//printf("newPoints: %lf,%lf",newPoint.trueRate,newPoint.falseRate);

		dataPoints.push_back(newPoint);
	}

	return dataPoints;
}






//TO DO---------------------------------------------------------------------
//Loop through the image to compute the harris corner values as described in class
// srcImage:  grayscale of original image
// harrisImage:  populate the harris values per pixel in this image
void computeHarrisValues(CFloatImage &srcImage, CFloatImage &harrisImage, CFloatImage &sobelImage, CFloatImage &harrisMatrixSumImage)
{
    
	int w = srcImage.Shape().width;
    int h = srcImage.Shape().height;
	/*
	ofstream image("image.txt");
	ofstream sobel("sobel.txt");
	ofstream harris("harris.txt");
	ofstream harris2("harris2.txt");
	int count =0;
	*/
	
	CFloatImage harrisMatrixImage(srcImage.Shape().width,srcImage.Shape().height,4);
	

	//gaussian mask with sigma 1.25
	double gaussianMask[5][5]={{0.0085,0.0223,0.0307,0.0223,0.0085},
	                           {0.0223,0.0583,0.0802,0.0583,0.0223},
	                           {0.0307,0.0802,0.1105,0.0802,0.0307},
	                           {0.0223,0.0583,0.0802,0.0583,0.0223},
	                           {0.0085,0.0223,0.0307,0.0223,0.0085}};



	//first , to calculate the gradient of each pixel
	//because the sobel operator needs 3*3 matrix, so the margin will be 0 instead
	//also the harris matrix image
	//maybe this part should be deleted

	harrisMatrixImage.ClearPixels();
	harrisMatrixSumImage.ClearPixels();



	//the rest of the srcimage should apply sobel operator to generate the sobelmatrix 
	//and then ,generate the single harris matrix , not the sum!!
    for (int y = 1; y < h-1; y++) {
        for (int x = 1; x < w-1; x++) {
           
			// TODO:  Compute the harris score for 'srcImage' at this pixel and store in 'harrisImage'.  See the project
            //   page for pointers on how to do this
			//0 channel for gradient in x direction ,and 1 channel for y
			sobelImage.Pixel(x,y,0) = (srcImage.Pixel(x+1,y-1,0)+2*srcImage.Pixel(x+1,y,0)+srcImage.Pixel(x+1,y+1,0))-(srcImage.Pixel(x-1,y-1,0)+2*srcImage.Pixel(x-1,y,0)+srcImage.Pixel(x-1,y+1,0));
			sobelImage.Pixel(x,y,1) = (srcImage.Pixel(x+1,y+1,0)+2*srcImage.Pixel(x,y+1,0)+srcImage.Pixel(x-1,y+1,0))-(srcImage.Pixel(x+1,y-1,0)+2*srcImage.Pixel(x,y-1,0)+srcImage.Pixel(x-1,y-1,0));
    
			harrisMatrixImage.Pixel(x,y,0) =  sobelImage.Pixel(x,y,0) * sobelImage.Pixel(x,y,0);
			harrisMatrixImage.Pixel(x,y,1) =  sobelImage.Pixel(x,y,0) * sobelImage.Pixel(x,y,1);
			harrisMatrixImage.Pixel(x,y,2) =  sobelImage.Pixel(x,y,1) * sobelImage.Pixel(x,y,0);
			harrisMatrixImage.Pixel(x,y,3) =  sobelImage.Pixel(x,y,1) * sobelImage.Pixel(x,y,1);
			/*
			image <<  "["<< x<< "," <<y<<"]" << srcImage.Pixel(x,y,0) << "   ";
			sobel << "["<< x<< "," <<y<<"]" << sobelImage.Pixel(x,y,0) << "," << sobelImage.Pixel(x,y,1) << "   ";
			harris << "["<< x<< "," <<y<<"]" << harrisMatrixImage.Pixel(x,y,0) << ","  <<harrisMatrixImage.Pixel(x,y,1) << "," << harrisMatrixImage.Pixel(x,y,2) << "," <<harrisMatrixImage.Pixel(x,y,3) << "   ";
			count ++;
			if (count >= w-2 ){
				count =0;
				image << endl << endl << endl;
				sobel << endl << endl << endl;
				harris << endl << endl << endl;
			}
			*/

		}
    }

	
//	count =0;

	//after the gradient,we should calculate the c for each pixel, within this process, we apply a 5*5 gaussian matrix to weight
	//and three pixel width in total of the margin will be ignored
	for (int y = 3; y < h-3; y++) {
        for (int x = 3; x < w-3; x++) {
			harrisMatrixSumImage.Pixel(x,y,0) = gaussianMask[0][0]*harrisMatrixImage.Pixel(x-2,y-2,0) + gaussianMask[0][1]*harrisMatrixImage.Pixel(x-1,y-2,0) + gaussianMask[0][2]*harrisMatrixImage.Pixel(x,y-2,0) + gaussianMask[0][3]*harrisMatrixImage.Pixel(x+1,y-2,0) + gaussianMask[0][4]*harrisMatrixImage.Pixel(x+2,y-2,0) +
				                                gaussianMask[1][0]*harrisMatrixImage.Pixel(x-2,y-1,0) + gaussianMask[1][1]*harrisMatrixImage.Pixel(x-1,y-1,0) + gaussianMask[1][2]*harrisMatrixImage.Pixel(x,y-1,0) + gaussianMask[1][3]*harrisMatrixImage.Pixel(x+1,y-1,0) + gaussianMask[1][4]*harrisMatrixImage.Pixel(x+2,y-1,0) +
												gaussianMask[2][0]*harrisMatrixImage.Pixel(x-2,y,0)   + gaussianMask[2][1]*harrisMatrixImage.Pixel(x-1,y,0)   + gaussianMask[2][2]*harrisMatrixImage.Pixel(x,y,0)   + gaussianMask[2][3]*harrisMatrixImage.Pixel(x+1,y,0)   + gaussianMask[2][4]*harrisMatrixImage.Pixel(x+2,y,0)   +
												gaussianMask[3][0]*harrisMatrixImage.Pixel(x-2,y+1,0) + gaussianMask[3][1]*harrisMatrixImage.Pixel(x-1,y+1,0) + gaussianMask[3][2]*harrisMatrixImage.Pixel(x,y+1,0) + gaussianMask[3][3]*harrisMatrixImage.Pixel(x+1,y+1,0) + gaussianMask[3][4]*harrisMatrixImage.Pixel(x+2,y+1,0) +
												gaussianMask[4][0]*harrisMatrixImage.Pixel(x-2,y+2,0) + gaussianMask[4][1]*harrisMatrixImage.Pixel(x-1,y+2,0) + gaussianMask[4][2]*harrisMatrixImage.Pixel(x,y+2,0) + gaussianMask[4][3]*harrisMatrixImage.Pixel(x+1,y+2,0) + gaussianMask[4][4]*harrisMatrixImage.Pixel(x+2,y+2,0);
			
			harrisMatrixSumImage.Pixel(x,y,1) = gaussianMask[0][0]*harrisMatrixImage.Pixel(x-2,y-2,1) + gaussianMask[0][1]*harrisMatrixImage.Pixel(x-1,y-2,1) + gaussianMask[0][2]*harrisMatrixImage.Pixel(x,y-2,1) + gaussianMask[0][3]*harrisMatrixImage.Pixel(x+1,y-2,1) + gaussianMask[0][4]*harrisMatrixImage.Pixel(x+2,y-2,1) +
				                                gaussianMask[1][0]*harrisMatrixImage.Pixel(x-2,y-1,1) + gaussianMask[1][1]*harrisMatrixImage.Pixel(x-1,y-1,1) + gaussianMask[1][2]*harrisMatrixImage.Pixel(x,y-1,1) + gaussianMask[1][3]*harrisMatrixImage.Pixel(x+1,y-1,1) + gaussianMask[1][4]*harrisMatrixImage.Pixel(x+2,y-1,1) +
												gaussianMask[2][0]*harrisMatrixImage.Pixel(x-2,y,1)   + gaussianMask[2][1]*harrisMatrixImage.Pixel(x-1,y,1)   + gaussianMask[2][2]*harrisMatrixImage.Pixel(x,y,1)   + gaussianMask[2][3]*harrisMatrixImage.Pixel(x+1,y,1)   + gaussianMask[2][4]*harrisMatrixImage.Pixel(x+2,y,1)   +
												gaussianMask[3][0]*harrisMatrixImage.Pixel(x-2,y+1,1) + gaussianMask[3][1]*harrisMatrixImage.Pixel(x-1,y+1,1) + gaussianMask[3][2]*harrisMatrixImage.Pixel(x,y+1,1) + gaussianMask[3][3]*harrisMatrixImage.Pixel(x+1,y+1,1) + gaussianMask[3][4]*harrisMatrixImage.Pixel(x+2,y+1,1) +
												gaussianMask[4][0]*harrisMatrixImage.Pixel(x-2,y+2,1) + gaussianMask[4][1]*harrisMatrixImage.Pixel(x-1,y+2,1) + gaussianMask[4][2]*harrisMatrixImage.Pixel(x,y+2,1) + gaussianMask[4][3]*harrisMatrixImage.Pixel(x+1,y+2,1) + gaussianMask[4][4]*harrisMatrixImage.Pixel(x+2,y+2,1);

			harrisMatrixSumImage.Pixel(x,y,2) = gaussianMask[0][0]*harrisMatrixImage.Pixel(x-2,y-2,2) + gaussianMask[0][1]*harrisMatrixImage.Pixel(x-1,y-2,2) + gaussianMask[0][2]*harrisMatrixImage.Pixel(x,y-2,2) + gaussianMask[0][3]*harrisMatrixImage.Pixel(x+1,y-2,2) + gaussianMask[0][4]*harrisMatrixImage.Pixel(x+2,y-2,2) +
				                                gaussianMask[1][0]*harrisMatrixImage.Pixel(x-2,y-1,2) + gaussianMask[1][1]*harrisMatrixImage.Pixel(x-1,y-1,2) + gaussianMask[1][2]*harrisMatrixImage.Pixel(x,y-1,2) + gaussianMask[1][3]*harrisMatrixImage.Pixel(x+1,y-1,2) + gaussianMask[1][4]*harrisMatrixImage.Pixel(x+2,y-1,2) +
												gaussianMask[2][0]*harrisMatrixImage.Pixel(x-2,y,2)   + gaussianMask[2][1]*harrisMatrixImage.Pixel(x-1,y,2)   + gaussianMask[2][2]*harrisMatrixImage.Pixel(x,y,2)   + gaussianMask[2][3]*harrisMatrixImage.Pixel(x+1,y,2)   + gaussianMask[2][4]*harrisMatrixImage.Pixel(x+2,y,2)   +
												gaussianMask[3][0]*harrisMatrixImage.Pixel(x-2,y+1,2) + gaussianMask[3][1]*harrisMatrixImage.Pixel(x-1,y+1,2) + gaussianMask[3][2]*harrisMatrixImage.Pixel(x,y+1,2) + gaussianMask[3][3]*harrisMatrixImage.Pixel(x+1,y+1,2) + gaussianMask[3][4]*harrisMatrixImage.Pixel(x+2,y+1,2) +
												gaussianMask[4][0]*harrisMatrixImage.Pixel(x-2,y+2,2) + gaussianMask[4][1]*harrisMatrixImage.Pixel(x-1,y+2,2) + gaussianMask[4][2]*harrisMatrixImage.Pixel(x,y+2,2) + gaussianMask[4][3]*harrisMatrixImage.Pixel(x+1,y+2,2) + gaussianMask[4][4]*harrisMatrixImage.Pixel(x+2,y+2,2);

			harrisMatrixSumImage.Pixel(x,y,3) = gaussianMask[0][0]*harrisMatrixImage.Pixel(x-2,y-2,3) + gaussianMask[0][1]*harrisMatrixImage.Pixel(x-1,y-2,3) + gaussianMask[0][2]*harrisMatrixImage.Pixel(x,y-2,3) + gaussianMask[0][3]*harrisMatrixImage.Pixel(x+1,y-2,3) + gaussianMask[0][4]*harrisMatrixImage.Pixel(x+2,y-2,3) +
				                                gaussianMask[1][0]*harrisMatrixImage.Pixel(x-2,y-1,3) + gaussianMask[1][1]*harrisMatrixImage.Pixel(x-1,y-1,3) + gaussianMask[1][2]*harrisMatrixImage.Pixel(x,y-1,3) + gaussianMask[1][3]*harrisMatrixImage.Pixel(x+1,y-1,3) + gaussianMask[1][4]*harrisMatrixImage.Pixel(x+2,y-1,3) +
												gaussianMask[2][0]*harrisMatrixImage.Pixel(x-2,y,3)   + gaussianMask[2][1]*harrisMatrixImage.Pixel(x-1,y,3)   + gaussianMask[2][2]*harrisMatrixImage.Pixel(x,y,3)   + gaussianMask[2][3]*harrisMatrixImage.Pixel(x+1,y,3)   + gaussianMask[2][4]*harrisMatrixImage.Pixel(x+2,y,3)   +
												gaussianMask[3][0]*harrisMatrixImage.Pixel(x-2,y+1,3) + gaussianMask[3][1]*harrisMatrixImage.Pixel(x-1,y+1,3) + gaussianMask[3][2]*harrisMatrixImage.Pixel(x,y+1,3) + gaussianMask[3][3]*harrisMatrixImage.Pixel(x+1,y+1,3) + gaussianMask[3][4]*harrisMatrixImage.Pixel(x+2,y+1,3) +
												gaussianMask[4][0]*harrisMatrixImage.Pixel(x-2,y+2,3) + gaussianMask[4][1]*harrisMatrixImage.Pixel(x-1,y+2,3) + gaussianMask[4][2]*harrisMatrixImage.Pixel(x,y+2,3) + gaussianMask[4][3]*harrisMatrixImage.Pixel(x+1,y+2,3) + gaussianMask[4][4]*harrisMatrixImage.Pixel(x+2,y+2,3);
	        	
			//it's possible that the denominator is zero. so the value will become 1.#IND000000
			if ( (harrisMatrixSumImage.Pixel(x,y,0) + harrisMatrixSumImage.Pixel(x,y,3))==0 ){
				harrisImage.Pixel(x,y,0) = 0;
			}
			else {
				harrisImage.Pixel(x,y,0) = (harrisMatrixSumImage.Pixel(x,y,0) * harrisMatrixSumImage.Pixel(x,y,3) - harrisMatrixSumImage.Pixel(x,y,1) * harrisMatrixSumImage.Pixel(x,y,2)) / (harrisMatrixSumImage.Pixel(x,y,0) + harrisMatrixSumImage.Pixel(x,y,3));
//	            harrisImage.Pixel(x,y,0) +=harrisImage.Pixel(x,y,0);
			}
			/*
			harris2 << "["<< x<< "," <<y<<"]" << harrisMatrixSumImage.Pixel(x,y,0) <<","<<harrisMatrixSumImage.Pixel(x,y,1)<<","<<harrisMatrixSumImage.Pixel(x,y,2)<<","<<harrisMatrixSumImage.Pixel(x,y,3)<<"," << harrisImage.Pixel(x,y,0) << "   ";
			count ++;
			if(count >= w-6){
				count =0 ;
				harris2 << endl << endl << endl;
			}
			*/

		}
	}

	/*
	image.close();
	sobel.close();
	harris.close();
	harris2.close();
	*/


}



// TO DO---------------------------------------------------------------------
// Loop through the harrisImage to threshold and compute the local maxima in a neighborhood
// srcImage:  image with Harris values
// destImage: Assign 1 to a pixel if it is above a threshold and is the local maximum in 3x3 window, 0 otherwise.
//    You'll need to find a good threshold to use.
void computeLocalMaxima(CFloatImage &srcImage,CByteImage &destImage)
{
	//the threshold to filter out the unsatisfied feature point
	//the point with the bigger harris operator will be better
	//so higher the threshold will lower down the number of features,vice versa
	double threshold = 0.2;

	int w = srcImage.Shape().width;
	int h = srcImage.Shape().height;
	/*
	double sum = 0.0;
	for (int y = 0; y < h; y++) {
        for (int x = 0; x < w; x++) {
			sum += srcImage.Pixel(x,y,0);
		}
	}
	double mean = sum /(w*h);
	threshold = mean*1.618;
	*/

	int count = 0;
	while (count<200 || count >1000){
		if(count <200){
			threshold = threshold /2 ;
		}
		else if (count >1000){
			threshold = threshold *1.618;
		}
		count =0;
		for (int y = 0; y < h; y++) {
			for (int x = 0; x < w; x++) {
	     		if (srcImage.Pixel(x,y,0) <= threshold){
		     		destImage.Pixel(x,y,0) = 0;
	     		}
     			else {
		    		double highest = srcImage.Pixel(x,y,0);

		    		destImage.Pixel(x,y,0) = 1;
		    		for (int i1=y-1; i1<=y+1 ;i1++){
		     			for(int i2=x-1; i2<=x+1 ;i2++){
	     					if ( srcImage.Pixel(i2,i1,0) > highest ){
					    		destImage.Pixel(x,y,0) = 0;
				    		}
			    		}
			    	}
					if(destImage.Pixel(x,y,0) == 1){
						count++;
					}
	  	    	}
	    	}
    	}
	}
	


}

double computePN(CFloatImage &srcImage,double xP , double yP, double xN , double yN , double centerX , double centerY)
{
	double xPXFloor = floor(xP);
	double xPXCeil = xPXFloor + 1;
	double xPYFloor = floor(yP);
	double xPYCeil = xPYFloor + 1;

	double xNXFloor = floor(xN);
	double xNXCeil = xNXFloor + 1;
	double xNYFloor = floor(yN);
	double xNYCeil = xNYFloor + 1;

	double grayXP = srcImage.Pixel(xPXFloor,xPYFloor,0) * ( xPXCeil - centerX ) * (xPYCeil - centerY)
		+ srcImage.Pixel(xPXCeil,xPYFloor,0) * (centerX - xPXFloor) * (xPYCeil - centerY)
		+ srcImage.Pixel(xPXFloor,xPYCeil,0) * (xPXCeil - centerX) * (centerY - xPYFloor)
		+ srcImage.Pixel(xPXCeil,xPYCeil,0) * (centerX - xPXFloor) * (centerY - xPYFloor);

	double grayXN = srcImage.Pixel(xNXFloor,xNYFloor,0) * ( xNXCeil - centerX ) * (xNYCeil - centerY)
		+ srcImage.Pixel(xNXCeil,xNYFloor,0) * (centerX - xNXFloor) * (xNYCeil - centerY)
		+ srcImage.Pixel(xNXFloor,xNYCeil,0) * (xNXCeil - centerX) * (centerY - xNYFloor)
		+ srcImage.Pixel(xNXCeil,xNYCeil,0) * (centerX - xNXFloor) * (centerY - xNYFloor);

	return (grayXP - grayXN);
}

void ComputeHarrisFeatures(CFloatImage &image, FeatureSet &features)
{
	//Create grayscale image used for Harris detection
	CFloatImage grayImage=ConvertToGray(image);

	//Create image to store Harris values

	CFloatImage sobel(image.Shape().width,image.Shape().height,2);
	CFloatImage harrisMatrixSumImage(image.Shape().width,image.Shape().height,4);
	CFloatImage harrisImage(image.Shape().width,image.Shape().height,1);

	//Create image to store local maximum harris values as 1, other pixels 0
	CByteImage harrisMaxImage(image.Shape().width,image.Shape().height,1);

	sobel.ClearPixels();
	harrisImage.ClearPixels();
	harrisMaxImage.ClearPixels();

	
	//compute Harris values puts harris values at each pixel position in harrisImage. 
	//You'll need to implement this function.
    computeHarrisValues(grayImage, harrisImage, sobel,harrisMatrixSumImage);
	
	// Threshold the harris image and compute local maxima.  You'll need to implement this function.
	computeLocalMaxima(harrisImage,harrisMaxImage);
    
    // Prints out the harris image for debugging purposes
	CByteImage tmp(harrisImage.Shape());
	convertToByteImage(harrisImage, tmp);
    WriteFile(tmp, "harris.tga");
    

	// TO DO--------------------------------------------------------------------
	//Loop through feature points in harrisMaxImage and create feature descriptor 
	//for each point above a threshold

//	ofstream descriptor("histogram.txt");

	int id = 0;

    for (int y=0;y<harrisMaxImage.Shape().height;y++) {
		for (int x=0;x<harrisMaxImage.Shape().width;x++) {
		
			// Skip over non-maxima
            if (harrisMaxImage.Pixel(x, y, 0) == 0)
                continue;

            //TO DO---------------------------------------------------------------------
		    // Fill in feature with descriptor data here. 

			
           
#if 0
			//First ,self-implemented simple descrptor
			//for SSD ,we need a 5*5 matrix and modify the data in FeatureSet.h
			//don't worry about the v will be in the margin and program will crash
			//actually, this part has been done in computeharrisFeature
			Feature f;

			f.type = 2;
			id++;
		    f.id =id;
			f.x = x;
		    f.y = y;

			for (int v=y-2; v<=y+2 ;v++){
				for(int u=x-2; u<=x+2 ;u++){
					//the original image is rgb
					double r = image.Pixel(u,v,0);
			        double g = image.Pixel(u,v,1);
			        double b = image.Pixel(u,v,2);

					f.data.push_back(r);
					f.data.push_back(g);
					f.data.push_back(b);
				}
			}
#endif


			//Second, using the sift descriptor
			//circumstance of calculating the sift but with out 16*16 enough pixel exists
	/*		if ( x <13 || x >harrisMaxImage.Shape().width-14 || y <13 || y>harrisMaxImage.Shape().height-14){
				continue;
			}
			*/
			if ( x <46 || x >harrisMaxImage.Shape().width-47 || y <46|| y>harrisMaxImage.Shape().height-47){
				continue;
			}

			 Feature f;

			f.type = 2;
			id++;
		    f.id =id;
			f.x = x;
		    f.y = y;

			//Descriptor with rotation invariance
			//Assuming that the original angle is upperward ,which is +90 degree to x+
			//and eigenvector has a angleRadians, substract the +90 will generate the angle to perform rotation matrix
			//+thelta for counter_clockwise, - for clockwise



			double a = harrisMatrixSumImage.Pixel(x,y,0);
			double b = harrisMatrixSumImage.Pixel(x,y,1);
			double c = harrisMatrixSumImage.Pixel(x,y,2);
			double d = harrisMatrixSumImage.Pixel(x,y,3);

			double T = a + d;
			double D = a*d - b*c;
            
			double L1 = T/2 + sqrt((T*T)/4-D);
			double L2 = T/2 - sqrt((T*T)/4-D);
			
			double vectorX = 0.0;
			double vectorY = 0.0;

			if (c != 0){
				vectorX = L1 -d;
				vectorY = c;
			}
			else if ( b != 0){
				vectorX = b;
				vectorY = L1 - a;
			}
			else if (b==0 && c==0){
				vectorX = 1;
				vectorY = 0;
			}

			f.angleRadians = atan2 (vectorY , vectorX);

//			descriptor << "["<< x<< "," <<y<<"]"  << " id: " << f.id  <<endl;
//			descriptor << "[Eigenvector: x,y]" << vectorX << ","<< vectorY << "  [angleOfEigenVector:]" << f.angleRadians <<endl;

			for(int iy1=0; iy1<4 ; iy1++){
				for(int ix1=0; ix1<4 ; ix1++){
					double histogram[8] = {0,0,0,0,0,0,0,0};
					for(int iy2=0; iy2<16 ; iy2++){
						for(int ix2=0; ix2<16 ; ix2++){
							double angleToRotate = f.angleRadians - (PI/2);

							//coordiantion transformation,twice
							int xLocal = x-32+(16*ix1)+ix2;
						    int yLocal = y-32+(16*iy1)+iy2;
							double xMiddle = cos(angleToRotate)*(xLocal-x) - sin(angleToRotate)*(y-yLocal);
							double yMiddle = sin(angleToRotate)*(xLocal-x) + cos(angleToRotate)*(y-yLocal);
							double xNew = xMiddle + x;
							double yNew = y - yMiddle;

							double xNewPositiveX = cos(angleToRotate)*(xLocal+1-x) - sin(angleToRotate)*(y-yLocal) + x;
							double xNewPositiveY = y - (sin(angleToRotate)*(xLocal+1-x) + cos(angleToRotate)*(y-yLocal));
							double xNewNegativeX = cos(angleToRotate)*(xLocal-1-x) - sin(angleToRotate)*(y-yLocal) + x;
							double xNewNegativeY = y - (sin(angleToRotate)*(xLocal-1-x) + cos(angleToRotate)*(y-yLocal));
							double yNewPositiveX = cos(angleToRotate)*(xLocal-x) - sin(angleToRotate)*(y-1-yLocal) + x;
							double yNewPositiveY = y-(sin(angleToRotate)*(xLocal-x) + cos(angleToRotate)*(y-1-yLocal));
							double yNewNegativeX = cos(angleToRotate)*(xLocal-x) - sin(angleToRotate)*(y+1-yLocal) + x;
							double yNewNegativeY = y-(sin(angleToRotate)*(xLocal-x) + cos(angleToRotate)*(y+1-yLocal));
#if 0
							//force into integer
							double xFloor = floor(xNew) ;
							double xCeil  = ceil(xNew); 
							double yFloor = floor(yNew) ;
							double yCeil  = ceil(yNew) ;

							//intepolation
							//There might be some problem , the matrix rotated maybe comtains the index overflow the original matrix
							//option: using sobel to interpolate, or using gray image to interpolate

							double Ix = sobel.Pixel(xFloor,yFloor,0) * (xCeil - xNew) * (yCeil - yNew) 
								+sobel.Pixel(xCeil,yFloor,0) * (xNew - xFloor) * (yCeil - yNew)
								+sobel.Pixel(xFloor,yCeil,0) * (xCeil - xNew) * (yNew - yFloor)
								+sobel.Pixel(xCeil,yCeil,0) * (xNew - xFloor) * (yNew - yFloor);

							double Iy = sobel.Pixel(xFloor,yFloor,1) * (xCeil - xNew) * (yCeil - yNew) 
								+sobel.Pixel(xCeil,yFloor,1) * (xNew - xFloor) * (yCeil - yNew)
								+sobel.Pixel(xFloor,yCeil,1) * (xCeil - xNew) * (yNew - yFloor)
								+sobel.Pixel(xCeil,yCeil,1) * (xNew - xFloor) * (yNew - yFloor);
#endif

							
							double Ix = computePN( grayImage , xNewPositiveX , xNewPositiveY , xNewNegativeX , xNewNegativeY , xNew ,yNew );
							double Iy = computePN( grayImage , yNewPositiveX , yNewPositiveY , yNewNegativeX , yNewNegativeY , xNew ,yNew );
/*
							if ( Ix == 0 && Iy == 0){
								continue;
							}
*/
							double module = sqrt(Ix*Ix + Iy*Iy);
							double angleToX = atan2 (Iy , Ix);
							//we need to record the corresponding angle between the IXIY vector rotated to eigenvector
							double angle = f.angleRadians - angleToX;

							if(angle <= (-1)*PI){
								angle = (-1)*angle -PI;
							}
							else if (angle > PI){
								angle = PI - angle;
							}

//							descriptor << "[xNew,yNew]" << xNew << "," << yNew << "    " << " [Ix,Iy] " << Ix <<"," << Iy 
//								<< " [angle] " << angle  <<" [module] "<<module<<endl; 
								
							if     (angle < ((-3)*PI)/4)  {histogram[0] += module;}
							else if(angle < ((-1)*PI)/2)  {histogram[1] += module;}
							else if(angle < ((-1)*PI)/4)  {histogram[2] += module;}
							else if(angle < 0)            {histogram[3] += module;}
							else if(angle < PI/4)         {histogram[4] += module;}
							else if(angle < PI/2)         {histogram[5] += module;}
							else if(angle < (3*PI)/4)     {histogram[6] += module;}
							else                          {histogram[7] += module;}
							  
						}
					}
					
//					descriptor << "--------------------------------------------------" << endl;
					
					for (int num=0; num<8 ;num++){
				     	f.data.push_back(histogram[num]);
//				        descriptor << histogram[num] << "    ";
				    }
//			     	descriptor << endl;

				}
			}

//			descriptor << "####################################################" << endl <<endl <<endl;

#if 0
			for(int iy1=0; iy1<4 ; iy1++){
				for(int ix1=0; ix1<4 ; ix1++){
					double histogram[8] = {0,0,0,0,0,0,0,0};
					for(int iy2=0; iy2<4 ; iy2++){
						for(int ix2=0; ix2<4 ; ix2++){
							double angle = atan2 ( sobel.Pixel(x-7+(4*ix1)+ix2,y-7+(4*iy1)+iy2,1) , sobel.Pixel(x-7+(4*ix1)+ix2,y-7+(4*iy1)+iy2,0) );
							if     (angle < ((-3)*PI)/4)  {histogram[0]++;}
							else if(angle < ((-1)*PI)/2)  {histogram[1]++;}
							else if(angle < ((-1)*PI)/4)  {histogram[2]++;}
							else if(angle < 0)            {histogram[3]++;}
							else if(angle < PI/4)         {histogram[4]++;}
							else if(angle < PI/2)         {histogram[5]++;}
							else if(angle < (3*PI)/4)     {histogram[6]++;}
							else                          {histogram[7]++;}
						}
					}

					for (int num=0; num<8 ;num++){
				     	f.data.push_back(histogram[num]);
				    }


				}
			}

#endif
         
            features.push_back(f);

        }
	}
//	descriptor.close();


}


// Compute silly example features.  This doesn't do anything
// meaningful.
void dummyComputeFeatures(CFloatImage &image, FeatureSet &features) {
	//Create grayscale image used for Harris detection
	CFloatImage grayImage=ConvertToGray(image);

	//Create image to store Harris values

	CFloatImage sobel(image.Shape().width,image.Shape().height,2);
	CFloatImage harrisMatrixSumImage(image.Shape().width,image.Shape().height,4);
	CFloatImage harrisImage(image.Shape().width,image.Shape().height,1);

	//Create image to store local maximum harris values as 1, other pixels 0
	CByteImage harrisMaxImage(image.Shape().width,image.Shape().height,1);

	sobel.ClearPixels();
	harrisImage.ClearPixels();
	harrisMaxImage.ClearPixels();

	
	//compute Harris values puts harris values at each pixel position in harrisImage. 
	//You'll need to implement this function.
    computeHarrisValues(grayImage, harrisImage, sobel,harrisMatrixSumImage);
	
	// Threshold the harris image and compute local maxima.  You'll need to implement this function.
	computeLocalMaxima(harrisImage,harrisMaxImage);
    
    // Prints out the harris image for debugging purposes
	CByteImage tmp(harrisImage.Shape());
	convertToByteImage(harrisImage, tmp);
    WriteFile(tmp, "harris.tga");
    

	// TO DO--------------------------------------------------------------------
	//Loop through feature points in harrisMaxImage and create feature descriptor 
	//for each point above a threshold

//	ofstream descriptor("histogram.txt");

	int id = 0;

    for (int y=0;y<harrisMaxImage.Shape().height;y++) {
		for (int x=0;x<harrisMaxImage.Shape().width;x++) {
		
			// Skip over non-maxima
            if (harrisMaxImage.Pixel(x, y, 0) == 0)
                continue;

			Feature f;

			f.type = 2;
			id++;
		    f.id =id;
			f.x = x;
		    f.y = y;

			for (int v=y-2; v<=y+2 ;v++){
				for(int u=x-2; u<=x+2 ;u++){
					//the original image is rgb
					double r = image.Pixel(u,v,0);
			        double g = image.Pixel(u,v,1);
			        double b = image.Pixel(u,v,2);

					f.data.push_back(r);
					f.data.push_back(g);
					f.data.push_back(b);
				}
			}
		 features.push_back(f);
		}
	}
}



// Perform simple feature matching.  This just uses the SSD
// distance between two feature vectors, and matches a feature in the
// first image with the closest feature in the second image.  It can
// match multiple features in the first image to the same feature in
// the second image.
void ssdMatchFeatures(const FeatureSet &f1, const FeatureSet &f2, vector<FeatureMatch> &matches, double &totalScore) {
	int m = f1.size();
	int n = f2.size();

	while (!matches.empty()){
		matches.pop_back();
	}
	totalScore = 0;

	double d;
	double dBest;
	int idBest;

	for (int i=0; i<m; i++) {
		dBest = 1e100;
		idBest = 0;

		for (int j=0; j<n; j++) {
			d = distanceSSD(f1[i].data, f2[j].data);

			if (d < dBest) {
				dBest = d;
				idBest = f2[j].id;
			}
		}

		FeatureMatch matched;

		matched.id1 = f1[i].id;
		matched.id2 = idBest;
		matched.score = dBest;

		matches.push_back(matched);

		totalScore += matched.score;
	}
}

// TODO: Write this function to perform ratio feature matching.  
// This just uses the ratio of the SSD distance of the two best matches as the score
// and matches a feature in the first image with the closest feature in the second image.
// It can match multiple features in the first image to the same feature in
// the second image.  (See class notes for more information, and the sshMatchFeatures function above as a reference)
void ratioMatchFeatures(const FeatureSet &f1, const FeatureSet &f2, vector<FeatureMatch> &matches, double &totalScore) 
{
	double threshold = 0.9;

	int m = f1.size();
	int n = f2.size();

	while (!matches.empty()){
		matches.pop_back();
	}
	totalScore = 0;

	double d;
	double dBest,dSecond;
	int idBest,idSecond;

	for (int i=0; i<m; i++) {
		dBest = 1e100;
		idBest = 0;
		dSecond = 1e100;
		dSecond = 0;

		for (int j=0; j<n; j++) {
			d = distanceSSD(f1[i].data, f2[j].data);

			if (d < dBest) {
				dSecond = dBest;
				idSecond = idBest;
				dBest = d;
				idBest = f2[j].id;
			}
			else if(d < dSecond){
				dSecond = d;
				idSecond = f2[j].id;
			}
		}
		//after look trough all the features , go to the ratio
		double ratio = dBest / dSecond;
		if(ratio < threshold){
			 FeatureMatch matched;

		     matched.id1 = f1[i].id;
		     matched.id2 = idBest;
		     matched.score = dBest;

		     matches.push_back(matched);

		     totalScore += matched.score;
		}
	}
    
}





// Convert Fl_Image to CFloatImage.
bool convertImage(const Fl_Image *image, CFloatImage &convertedImage) {
	if (image == NULL) {
		return false;
	}

	// Let's not handle indexed color images.
	if (image->count() != 1) {
		return false;
	}

	int w = image->w();
	int h = image->h();
	int d = image->d();

	// Get the image data.
	const char *const *data = image->data();

	int index = 0;

	for (int y=0; y<h; y++) {
		for (int x=0; x<w; x++) {
			if (d < 3) {
				// If there are fewer than 3 channels, just use the
				// first one for all colors.
				convertedImage.Pixel(x,y,0) = ((uchar) data[0][index]) / 255.0f;
				convertedImage.Pixel(x,y,1) = ((uchar) data[0][index]) / 255.0f;
				convertedImage.Pixel(x,y,2) = ((uchar) data[0][index]) / 255.0f;
			}
			else {
				// Otherwise, use the first 3.
				convertedImage.Pixel(x,y,0) = ((uchar) data[0][index]) / 255.0f;
				convertedImage.Pixel(x,y,1) = ((uchar) data[0][index+1]) / 255.0f;
				convertedImage.Pixel(x,y,2) = ((uchar) data[0][index+2]) / 255.0f;
			}

			index += d;
		}
	}
	
	return true;
}

// Convert CFloatImage to CByteImage.
void convertToByteImage(CFloatImage &floatImage, CByteImage &byteImage) {
	CShape sh = floatImage.Shape();

    assert(floatImage.Shape().nBands == byteImage.Shape().nBands);
	for (int y=0; y<sh.height; y++) {
		for (int x=0; x<sh.width; x++) {
			for (int c=0; c<sh.nBands; c++) {
				float value = floor(255*floatImage.Pixel(x,y,c) + 0.5f);

				if (value < byteImage.MinVal()) {
					value = byteImage.MinVal();
				}
				else if (value > byteImage.MaxVal()) {
					value = byteImage.MaxVal();
				}

				// We have to flip the image and reverse the color
				// channels to get it to come out right.  How silly!
				byteImage.Pixel(x,sh.height-y-1,sh.nBands-c-1) = (uchar) value;
			}
		}
	}
}

// Compute SSD distance between two vectors.
double distanceSSD(const vector<double> &v1, const vector<double> &v2) {
	int m = v1.size();
	int n = v2.size();

	if (m != n) {
		// Here's a big number.
		return 1e100;
	}

	double dist = 0;

	for (int i=0; i<m; i++) {
		dist += pow(v1[i]-v2[i], 2);
	}

	
	return sqrt(dist);
}



// Transform point by homography.
void applyHomography(double x, double y, double &xNew, double &yNew, double h[9]) {
	double d = h[6]*x + h[7]*y + h[8];

	xNew = (h[0]*x + h[1]*y + h[2]) / d;
	yNew = (h[3]*x + h[4]*y + h[5]) / d;
}

// Compute AUC given a ROC curve
double computeAUC(vector<ROCPoint> &results)
{
	double auc=0;
	double xdiff,ydiff;
	for (int i = 1; i < (int) results.size(); i++)
    {
        //fprintf(stream,"%lf\t%lf\t%lf\n",thresholdList[i],results[i].falseRate,results[i].trueRate);
		xdiff=(results[i].falseRate-results[i-1].falseRate);
		ydiff=(results[i].trueRate-results[i-1].trueRate);
		auc=auc+xdiff*results[i-1].trueRate+xdiff*ydiff/2;
    	    
    }
	return auc;
}

