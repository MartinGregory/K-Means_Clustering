//
//  Centroid.hpp
//  K-Means_Clustering
//
//  Created by Martin Gregory Sendrowicz on 5/16/23.
//

#ifndef _Centroid_H_
#define _Centroid_H_

#include <istream>
#include <iostream>

//||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
class Centroid{
    
    friend class KMeans ;
      
public : // data members
      double xcoord ;       // Cartesian x-coordinate
      double ycoord ;       // Cartesian y-coordinate
      int label ;           // label denotes the cluster
      
public :
    /* Default Constructor : -1 is just a diagnistic value (debugging aid) */
    Centroid( double x=-1.0 , double y=-1.0 , int lbl=-1 ) ;
      
    // getters (inlined)
    int getX() { return xcoord ; }
    int getY() { return ycoord ; }
    int getLabel() { return label ; }
      
    // setters (inlined)
    void setX( int x ) { this->xcoord = x ; }
    void setY( int y ) { this->ycoord = y ; }
    void setLbl( int lbl ) { this->label = lbl ; }
    
    void showCentroid() ;
};//||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||


#endif /* _Centroid_H_ */
