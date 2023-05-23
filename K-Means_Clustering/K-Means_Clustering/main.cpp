//
//  main.cpp
//  K-Means_Clustering
//
//  Created by Martin Gregory Sendrowicz on 5/15/23.
//

#include <iostream>
#include <fstream>
//#include <limits>


#include "KMeans.h"

//....................................................................................................
template < typename T >
void is_open( T & stream , std::string file_name ){
    
    if( !stream ){
        std::cerr << "Error openeing the file:"<< file_name << std::endl ;
        exit( EXIT_FAILURE ) ;
    }
}//...................................................................................................

//||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
int main( int argc , const char * argv[] ) {
    
    // rand() is called in KMeans::seedCentroids()
    std::srand( static_cast<unsigned int>( std::time(0) )) ;

    if( argc != 3 ){
        std::cerr << "This program requires exactly 2 command line arguments:\n" ;
        std::cerr << "(1) input file name: e.g. input.txt\n" ;
        std::cerr << "(2) output file name: e.g. output.txt\n" ;
        return 1 ;
    }
    
    std::ifstream input  { argv[ 1 ] } ; is_open( input  , argv[ 1 ] ) ;
    std::ofstream output { argv[ 2 ] } ; is_open( output , argv[ 2 ] ) ;
    
    //std::cout << argv[ 1 ] << " : " << argv[ 2 ] << std::endl ;
    
    std::cin.clear() ;
    //std::cin.ignore( std::numeric_limits<std::streamsize>::max(),'\n' ) ;
    
    int k = 3 ;
    std::cout << "Enter K:" ;
    std::cin >> k ;
    
    std::unique_ptr<KMeans> kmpt { new KMeans( k ) } ;
    
    // Step-1: Extract data points into pointSet.txt
    std::ofstream out_ps { "pointSet.txt" } ; is_open( out_ps , "pointSet.txt" ) ;
    (*kmpt).extractPts( input , out_ps ) ;
    input.close() ;
    out_ps.close() ; //(*kmpt).show_data() ;
    
    std::ifstream in_ps  { "pointSet.txt" } ; is_open( in_ps  , "pointSet.txt" ) ;
    
    // Step-2: Load pointSet.txt into 'point_set' data member array
    (*kmpt).loadPointSet( in_ps ) ;
    in_ps.close() ; //(*kmpt).showPointSet() ;
    
    
    std::ofstream out_iter { "iterations.txt" } ; is_open( out_iter , "iterations.txt" ) ;
    (*kmpt).kMeansClustering( out_iter ) ;
    out_iter.close() ;
    
    (*kmpt).writePtSet( output ) ;
    output.close() ;
    
    return 0;
}//|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
