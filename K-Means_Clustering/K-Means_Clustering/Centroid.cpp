//
//  Centroid.cpp
//  K-Means_Clustering
//
//  Created by Martin Gregory Sendrowicz on 5/16/23.
//

#include "Centroid.h"

//....................................................................................................
Centroid::Centroid( double x , double y , int lbl ): xcoord{x}, ycoord(y), label(lbl){
    //std::cout << "Centroid default constructor called" << std::endl ;
}
//....................................................................................................
void Centroid::showCentroid(){
    std::cout <<"c"<< label<<" is at x="<< xcoord <<" y="<< ycoord << std::endl ;
}//...................................................................................................
