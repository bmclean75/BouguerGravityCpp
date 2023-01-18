/***************************************************************************
 *   Copyright (C) 2021 Ben F McLean                                       *
 *   drbenmclean@gmail.com                                                 *
 *                                                                         *
 ***************************************************************************/
#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <time.h>
#include <sstream>
#include <fstream>
#include <cstdio> //for stdin stdout
#include <string>
#include <cmath> //for sin()

//FGL_Numerics and the following files are unavailable due to licensing.
#include "FGL_Numerics.h" 
#include "FGL_Commandline.h"

#include "FGL_ValuePoint.h"
#include "FGL_ValuePointSet.h"
#include "FGL_Kriging2.h"
#include "FGL_Point.h"
#include "FGL_Vector.h"
#include "FGL_Plane.h"


//g++ -O3 -Wall -fpermissive -I /home/ben/code/archive/FGL_Numerics -I /home/ben/code/archive/FGL_Programming -I /home/ben/code/archive/FGL_Geometrics BouguerGrav.cpp -o BouguerGrav

//BouguerGrav LLEOset=rawgrav.txt headerlines=4 gravfactor=0.00001

FGL_ValuePointSet<double> lattitudeCorrection(const FGL_ValuePointSet<double> &rawGrav);
FGL_ValuePointSet<double> freeAirCorrection(const FGL_ValuePointSet<double> &rawGrav);
FGL_ValuePointSet<double> simpleBouguerCorrection(const FGL_ValuePointSet<double> &rawGrav);
	

using namespace std;

int main(int argc, char* argv[]){
	
	cout << endl;
	cout << "**********************************************************" << endl;
	cout << "       extract Bouguer gravity from raw gravity data"       << endl;
	cout << "**********************************************************" << endl << endl;
	
	cout << "       use gravity in mGal"       << endl;
	cout << "       specify --Isogal65 or --isogal84"       << endl;
	
	//raw grav data in Potsdam datum (Isogal 65), IGSN71 datum (Isogal 84), or AAGD07
	//https://www.publish.csiro.au/eg/pdf/EG17094
	
	//read commandline params
	FGL_Commandline commandline(argc,argv);
	
	//read in raw gravity file
	//************************
	
	//create a Lattitude Longitude Elevation Observation file for raw gravity data
	FGL_ValuePointSet<double> rawGravity;

	int headerlines = 0;
	if(commandline.check_option("headerlines")){
		headerlines = stoi(commandline.get_option("headerlines"));
	}
	
	//cout << "number of header lines = " << headerlines << endl;
	rawGravity.readValuePointSet(commandline.get_option("LLEOset"),headerlines);
	cout << "read in " << rawGravity.size() << " raw gravity datapoints..." << endl;
	
	//assume raw grav in mGal, convert to m/s^2
	for(int i=0;i<rawGravity.size();i++){
		rawGravity(i).setValue(rawGravity(i).value() * 0.00001);
	}
	
	//output station map... scale z to km for visual balance with lat/long
	rawGravity.writeValuePointSet_vtkpoints("mydata_initialgrav",1.0,1.0,-0.001,1.0);
	
	//check if data uses Isogal65 or Isogal85
	if(commandline.check_flag("isogal65")){
		for(int i=0;i<rawGravity.size();i++){
			//convert to isogal84, https://www.publish.csiro.au/eg/pdf/EG17094
			//cout << "convert from " << rawGravity(i).value() << " to ";
			rawGravity(i).setValue(9.7967188+1.00053*(rawGravity(i).value()-9.7968574));
			//cout << rawGravity(i).value() << endl;
		}
	}
	
	//create correction for lattitude effects
	FGL_ValuePointSet<double> lattCorr = lattitudeCorrection(rawGravity);
	FGL_ValuePointSet<double> freeAirCorr = freeAirCorrection(rawGravity);
	FGL_ValuePointSet<double> simpleBougCorr = simpleBouguerCorrection(rawGravity);
	
	//file out corrections
	lattCorr.writeValuePointSet_vtkpoints("mydata_lattCorr",1.0,1.0,-0.001,1.0);
	freeAirCorr.writeValuePointSet_vtkpoints("mydata_freeAirCorr",1.0,1.0,-0.001,1.0);
	simpleBougCorr.writeValuePointSet_vtkpoints("mydata_simpleBougCorr",1.0,1.0,-0.001,1.0);
	
	//apply corrections
	FGL_ValuePointSet<double> BouguerGrav = rawGravity; //get coordinates into the BouguerGrav ValuePointSet
	for(int i=0;i<rawGravity.size();i++){
		BouguerGrav(i).setValue(rawGravity[i].value() + lattCorr[i].value() + freeAirCorr[i].value() + simpleBougCorr[i].value());
	}
	
	BouguerGrav.writeValuePointSet_vtkpoints("BouguerGrav",1.0,1.0,-0.001,1.0);
	
	//kriging data
	//************
	//kriging currently uses vector<FGL_ValuePoint<double> >, eventually write for vps...
	
	//check no double-up points in (x,y) that will cause A to be singular!
	//(consider deleting near points as well, within a given radius, so matrix not close to singular = more inaccurate)
	cout << "Kriging dataset: " << endl;
	cout << "	erasing any doubled (x,y) locations to ensure matrix invertability..." << endl;
	int nErased = 0;
	for(int i=0;i<BouguerGrav.size();i++){
		for(int j=i+1;j<BouguerGrav.size();j++){ //start at point 2
			if(BouguerGrav[i].cx()==BouguerGrav[j].cx()
				&&BouguerGrav[i].cy()==BouguerGrav[j].cy()){
				nErased++;
				//cout << "erase j = " << j << " as same loc as i = " << i << endl;
				BouguerGrav.erase(j);
				j--;
			}
		}
	}
	cout << "	erased " << nErased << " datapoints, " << BouguerGrav.size() << " datapoints remaining..." << endl;
	
	vector <FGL_ValuePoint<double> > datapoints;
	for(int i=0;i<BouguerGrav.size();i++){
		datapoints.push_back(BouguerGrav(i));
	}
	
	//find min and max x/y extent of points
	cout << "	determine the range of the survey..." << endl;
	double  minCx, minCy, maxCx, maxCy;
	BouguerGrav.valuePointSet_extent(minCx,minCy,maxCx,maxCy);
	cout << "	survey ranges from (" << minCx << "," << minCy << ") to (" << maxCx << "," << maxCy << ")" << endl;
	double rangex = maxCx - minCx;
	double rangey = maxCy - minCy;
	//extend border by 50% each direction
	double ox = minCx - rangex/2.0;
	double oy = minCy - rangey/2.0;
	double tx = maxCx + rangex/2.0;
	double ty = maxCy + rangey/2.0;
	cout << "	survey will be gridded from (" << ox << "," << oy << ") to (" << tx << "," << ty << ")" << endl;
	//choose nx ny
	int nx = 100;
	int ny = 100;
	//calc dx dy
	double dx = (tx-ox)/nx;
	double dy = (ty-oy)/ny;
	//make a target array for kriged data
	FGL_Array3D<double> BouguerKriged(nx,ny,1,dx,dy,1.0,ox,oy,0.0);
	//krige datapoints
	cout << "	perform universal Kriging in the x/y plane..." << endl;
	universalkrige_allpoints_xy(datapoints,BouguerKriged); //can also include maxpairs in operands
	vtkout_STRUCTUREDPOINTS(BouguerKriged,"bouger_output");
	
	cout << "processing complete!" << endl;
  	return EXIT_SUCCESS;
	
}


FGL_ValuePointSet<double> lattitudeCorrection(const FGL_ValuePointSet<double> &rawGrav){
	
	FGL_ValuePointSet<double> corr; //set of corrections to be applied //could make and use copy constructor
	corr = rawGrav;//get coordinates into corr, requires copy constructor or operator
	
	//int nPoints = rawGrav.size();
	for(int i=0;i<rawGrav.size();i++){
		
		//https://en.wikipedia.org/wiki/Theoretical_gravity
		//International gravity formula 1980 series expansion
		//note lattitude is rawGrav.cx()
		double c1 = 5.2790414e-3;
		double c2 = 2.32718e-5;
		double c3 = 1.262e-7;
		double c4 = 7.0e-10;
		double latitudedeg = rawGrav[i].cx()*M_PI/180.0;
		double sin2theta = pow(sin(latitudedeg),2.0);
		double sin4theta = pow(sin(latitudedeg),4.0);
		double sin6theta = pow(sin(latitudedeg),6.0);
		double sin8theta = pow(sin(latitudedeg),8.0);
		double correction = -9.780327 * (1.0 + c1*sin2theta + c2*sin4theta  + c3*sin6theta  + c4*sin8theta);
		corr(i).setValue(correction);
		
		//cout << "correction applied = " << (1.0 + c1*sin2theta + c2*sin4theta  + c3*sin6theta  + c4*sin8theta) << endl;
		//cout << "rawGrav[i].cx() = " << rawGrav[i].cx() << endl;
	}
	
	return corr;
	
};

FGL_ValuePointSet<double> freeAirCorrection(const FGL_ValuePointSet<double> &rawGrav){
	
	FGL_ValuePointSet<double> corr; //set of corrections to be applied //could make and use copy constructor
	corr = rawGrav;//get coordinates into corr, requires copy constructor or operator
	
	for(int i=0;i<rawGrav.size();i++){
		
		//see https://en.wikipedia.org/wiki/Free-air_gravity_anomaly
		//note altitude is rawGrav.cz()
		//for a simple correction of elevation from sea level to altitude cz
		corr(i).setValue(0.000003086*rawGrav[i].cz());
		//see also https://en.wikipedia.org/wiki/Theoretical_gravity
		
		
	}
	
	return corr;
	
};


FGL_ValuePointSet<double> simpleBouguerCorrection(const FGL_ValuePointSet<double> &rawGrav){
	
	FGL_ValuePointSet<double> corr; //set of corrections to be applied //could make and use copy constructor
	corr = rawGrav;//get coordinates into corr, requires copy constructor or operator
	
	for(int i=0;i<rawGrav.size();i++){
		
		//see https://en.wikipedia.org/wiki/Bouguer_anomaly
		//note altitude is rawGrav.cz()
		//for a simple correction of a Bouguer plate of thickness from sea level to altitude cz
		double G = 6.67430e-11; //https://en.wikipedia.org/wiki/Gravitational_constant
		double rho = 2670.0; //kg/m^3
		corr(i).setValue(-2.0*M_PI*G*rho*rawGrav[i].cz());
	}
	
	return corr; //needs copy constructor?
	
};
