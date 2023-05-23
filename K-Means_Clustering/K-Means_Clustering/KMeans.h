//
//  KMeans.h
//  K-Means_Clustering
//
//  Created by Martin Gregory Sendrowicz on 5/15/23.
//

#ifndef _KMeans_H_
#define _KMeans_H_


#include "Point.h"
#include "Centroid.h"

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>      // std::stringstream iss( line ) ;
#include <memory>       // std::unique_ptr<T>
#include <cmath>        // sqrt()
#include <math.h>       // M_PI
#define _USE_MATH_DEFINES

#include <cstdlib>      // rand()
#include <time.h>       // time()


//||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
class KMeans {
    
    int _K ;        // number of clusters <- ask the user for the number of clusters
    int num_pts ;   // total number of ALL data points
    
    std::unique_ptr< Point[] > point_set ;
    /*Point* point_set ;   pointer to the 1-D array of Point Class to be dynamically allocated
                           at run-time  based on the number of given data points; initially set all
                           Point's 'distance' data member to DBL_MAX */
    
    // to be initialized from image header
    int num_rows ;       // number of rows
    int num_cols ;       // number of columns
    int min_val ;        // minimum data point value i.e. binary 0
    int max_val ;        // maximum data point value i.e. binary 1
    
    std::unique_ptr< std::unique_ptr< int[] >[] > data ;
    /* int** data ;          a 2D array, size of num_rows by num_cols i.e. this is the given
                             to be clustered data */

    
    std::unique_ptr< Centroid[] > _Kcentroids ;
    /*Centroid* _Kcentroids ;  to be dynamically allocated in class constructor, size = K */
    
    int change ;     // for tracking the label changes, initialize to 0
    
public :
    KMeans( int K=1 ) ;     // default constructor
    
    void extractPts( std::ifstream & in_file , std::ofstream & out_file ) ;
    void show_data() ;
    void loadPointSet( std::ifstream & in_ps ) ;
    void showPointSet() ;
    
    void randomSeedCentroids() ;
    void intelliSeed( double delta ) ;
    void kMeansClustering( std::ofstream & out_iter ) ;
    
    void computeCentroids() ;
    void distanceMinLabel( Point & point ) ;

    double computeDist( Point p1, Centroid p2) ;
    void writePtSet( std::ofstream & output ) ;
    void point2Image () ;
    void printImage( std::ofstream & out_file, int iteration) ;
private :
    bool compare_dbl( double a , double b ) ;
    // helper methods for my "Intelligent Seed" algorithm
    int find_pz_quadrant( Point & pz , Centroid & c0 ) ;
    void adjust_c0pz_to_zero( Point & pz , Centroid & c0 , int pz_quad ) ;
    double find_theta( Point & pz , int pz_quad , double dist_z ) ;
    void place_c1( double r , int pz_quad , Centroid & c1 , double theta  ) ;
    double find_c1_angle( int pz_quad , double theta  ) ;
    double place_ci_centroid( Centroid & c2 , int pz_quad , double c1_angle , double beta , double r ) ;
    void adjust_coords( Centroid & c0 , Centroid & ci , int ci_angle ) ;
    
};//||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

#endif /* KMeans_h */




