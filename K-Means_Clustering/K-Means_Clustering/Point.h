//
//  Point.h
//  K-Means_Clustering
//
//  Created by Martin Gregory Sendrowicz on 5/15/23.
//

#ifndef Point_h
#define Point_h

#include <limits>       // __DBL_MAX__ = std::numeric_limits<double>::max() ;
#include <iostream>

//||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
class Point {
    
    friend class KMeans ;
    
    private : // data members
    double xcoord ;     // point's Cartesian x-coordinate
    double ycoord ;     // point's Cartesian y-coordinate
    int label ;         // label is assigned according to the cluster to which the point belongs (varies dynamically)
    
    double distance ;   /* distance to its own cluster centroid. Initially set ALL points to largest
                        possible distance: this->distance = std::numeric_limits<double>::max() ;  */
    
    public :
    //................................................................................................
    // -1 is just a diagnostic value (a debugging aid)
    Point( double x=-1 , double y=-1 , int lbl=1 , double dist=__DBL_MAX__ ) ;   // default constructor
    //................................................................................................
    
    // getters (inlined)
    int getX() { return xcoord ; }
    int getY() { return ycoord ; }
    int getLabel() { return label ; }
    int getDist() { return distance ; }

    // setters (inlined)
    void setX( int x ) { this->xcoord = x ; }
    void setY( int y ) { this->ycoord = y ; }
    void setLbl( int lbl ) { this->label = lbl ; }
    void setDist( int dist ) { this->distance = dist ; }
    
    void showPoint();
};//||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

#endif /* Point_h */
