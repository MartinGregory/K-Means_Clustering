//
//  Point.cpp
//  K-Means_Clustering
//
//  Created by Martin Gregory Sendrowicz on 5/15/23.
//

#include "Point.h"

//................................................................................................
// -1 is just a diagnostic value (a debugging aid)
Point::Point( double x , double y , int lbl , double dist ): xcoord{x}, ycoord{y}, label{lbl}, distance{dist}{
    //std::cout << "Point default constructor called" << std::endl ;
}
//................................................................................................
void Point::showPoint(){
    std::cout << xcoord <<" : "<< ycoord <<" : "<< label <<" : "<< distance << std::endl ;
}//...............................................................................................

