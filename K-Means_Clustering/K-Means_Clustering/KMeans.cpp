//
//  KMeans.cpp
//  K-Means_Clustering
//
//  Created by Martin Gregory Sendrowicz on 5/15/23.
//

#include "KMeans.h"

//....................................................................................................
/* Default Constructor : -1 is just a diagnistic value (debugging aid) */
KMeans::KMeans( int K ) : _K{K}, num_pts{0}, num_rows{-1}, num_cols{-1}, min_val{-1}, max_val{-1}, change{-1} {
    (*this)._Kcentroids = std::unique_ptr< Centroid[] >( new Centroid[ K ] ) ;
    for( int i=1 ; i<=K ; i++ )
        (*this)._Kcentroids[ i-1 ] = Centroid( i,i,i ) ;
    //std::cout << "KMeans default constructor called" << std::endl ;
}
//....................................................................................................
/* Algorithm - steps for extracting points from a binary data image:
 step-0:    infile ← open input.txt image file (DONE in main())
            rows,cols,min,max ← read image_header i.e. data parameters
            data ← dynamically allocate a 2-D array of size num_rows × num_cols.
 
 step-1:    read the image L→R & T→B (i.e. Left-to-Right , Top-to-Bottom)
 step-2:    p(x,y) ← read point from image     i.e. data point p with coordinates x,y

            if( p(x,y).value >= 1 )
                output  x,y         i.e. output x,y coordinates to PointSet out-file
                num_pts++           count how many points there are
 
 step-3:    repeat steps 1 to 2 until eof (‘end of file’)
 step-4:    point_set ← dynamically allocate a 1-D array of Points of size num_pts. */
void KMeans::extractPts( std::ifstream & input , std::ofstream & output ) {
  
    // step-0:
    input >> (*this).num_rows >> (*this).num_cols >> (*this).min_val >> (*this).max_val ;

    (*this).data = std::unique_ptr< std::unique_ptr< int[] >[] >( new std::unique_ptr<int[]>[ (*this).num_rows ] ) ;
    
    for( int row=0 ; row<(*this).num_rows ; row++ ){
        (*this).data[row] = std::unique_ptr< int[] >( new int[ (*this).num_cols ] ) ;
        
        for( int col=0 ; col<(*this).num_cols ; col++ )
            (*this).data[row][col] = 0 ;        //at first let's initialize 'data' with 0s
    }
    
    // step-1:
    int value=0 ; int row=0 ; int col=0 ;
    while( !input.eof() ){
        
        input >> value ;                        // stream extraction operator >> stops on white spaces
        //std::cout << value << " " ;           // ergo we are extracting one data point at a time
        
    // step-2:
        if( value > 0 ){                        // only collect the non-zero data points
            output << row <<" "<< col << std::endl ;
            (*this).num_pts++ ;
        }
        (*this).data[row][col] = value ;        // mirror the actual data into 'data'
        
        ++col ;
        if( col == (*this).num_cols ){
            ++row ; col=0 ;
            if( row==(*this).num_rows ) break ;
        }
    }//step-3: //std::cout << row <<" : "<< col << std::endl ;
    
    // step-4:
    (*this).point_set = std::unique_ptr< Point[] >( new Point[ (*this).num_pts ] ) ;
    for(int p=0 ; p<(*this).num_pts ; ++p)
        this->point_set[ p ] = Point() ;
}
//....................................................................................................
void KMeans::show_data(){
    
    for( int row=0 ; row<(*this).num_rows ; row++ ){
        for( int col=0 ; col<(*this).num_cols ; col++ )
            std::cout << (*this).data[row][col] << " " ;
        std::cout << std::endl ;
    }
}
//....................................................................................................
/* Read each Point p from pointSet.txt file into 'point_set' data member array. Initially set each
   Point's 'distance' data member to DBL_MAX. */
void KMeans::loadPointSet( std::ifstream & in_ps ) {

    std::cout <<"Number of Points:"<< (*this).num_pts << "\n" ;
    
    int row=0 ;
    while( row <(*this).num_pts ) {                  //OR while( !in_ps.eof() ) OR while( in_ps ) {
          
        in_ps >> (*this).point_set[ row ].xcoord ;
        in_ps >> (*this).point_set[ row ].ycoord ;
        point_set[ row ].distance = __DBL_MAX__ ;
        
        //std::cout <<"row:"<< row << " = " ; (*this).point_set[ row ].showPoint() ;
        row++ ;
    }
}
//....................................................................................................
void KMeans::showPointSet(){
    for( int i=0 ; i<(*this).num_pts ; ++i )
        (*this).point_set[ i ].showPoint() ;
}
//....................................................................................................
void KMeans::randomSeedCentroids(){

    /* Genereate a sequence of random numbers within ranges 0-to-'num_rows' and 0-to-'num_cols'*/
    int x=0,y=0 ;
    for( int i=0 ; i<(*this)._K ; ++i ){
        
again:
        x = rand() % (*this).num_rows ;
        y = rand() % (*this).num_cols ;
        
        if( i>0 ){
            for( int j=0 ; j<i ; ++j ){
                if( (x == (*this)._Kcentroids[ j ].xcoord) && (y == (*this)._Kcentroids[ j ].ycoord) )
                    goto again ;
            }
        }
        (*this)._Kcentroids[ i ].xcoord = double(x) ;
        (*this)._Kcentroids[ i ].ycoord = double(y) ;
        (*this)._Kcentroids[ i ].showCentroid() ;
    }
}
//....................................................................................................
bool KMeans::compare_dbl( double a, double b ){
    double epsilon = 0.0000001f ;     // for float use 0.01f
    //double epsilon = std::numeric_limits<double>::epsilon() ;
    //std::cout << "epsilon: "<< epsilon << std::endl ;   // on my system epsilon = 2.22045e-16
    return( std::fabs( a-b ) > epsilon )? true : false ;
}
//....................................................................................................
int KMeans::find_pz_quadrant( Point & pz, Centroid & c0 ){
    
    // FIX this for real value EPSILON &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    if( (pz.xcoord < c0.xcoord) && (pz.ycoord > c0.ycoord) )
        return 1 ;
    if( (pz.xcoord < c0.xcoord) && (pz.ycoord < c0.ycoord) )
        return 2 ;
    if( (pz.xcoord > c0.xcoord) && (pz.ycoord < c0.ycoord) )
        return 3 ;
    if( (pz.xcoord > c0.xcoord) && (pz.ycoord > c0.ycoord) )
        return 4 ;
    
    std::cerr<<"Function KMeans::find_pz_quadrant failed to find pz quadrant\n" ;
    exit( EXIT_FAILURE ) ;
}
//....................................................................................................
void KMeans::adjust_c0pz_to_zero( Point & pz, Centroid & c0, int pz_quad ){
    
    switch (pz_quad){
        case 1:
            pz.xcoord = std::abs( c0.xcoord - pz.xcoord ) ;
            pz.ycoord = std::abs( pz.ycoord - c0.ycoord ) ;
            break ;
        case 2:
            pz.xcoord = (-1)*( c0.xcoord - pz.xcoord ) ;
            pz.ycoord = std::abs( c0.ycoord - pz.ycoord ) ;
            break ;
        case 3:
            pz.xcoord = (-1)*( pz.xcoord - c0.xcoord ) ;
            pz.ycoord = (-1)*( c0.ycoord - pz.ycoord ) ;
            break ;
        case 4:
            pz.xcoord = std::abs( pz.xcoord - c0.xcoord ) ;
            pz.ycoord = (-1)*( pz.ycoord - c0.ycoord ) ;
            break ;
        default:
            std::cerr<<"Function KMeans::adjust_c0pz_to_zero failed!\n" ;
            exit( EXIT_FAILURE ) ;
    }
}
//....................................................................................................
double KMeans::find_theta( Point & pz, int pz_quad , double dist_z ){
    
    switch (pz_quad){
        case 1:
            return std::asinf( pz.ycoord / dist_z ) ;
        case 2:
            return std::asinf( pz.ycoord / dist_z ) ;
        case 3:
            return std::asinf( ((-1)*pz.ycoord) / dist_z ) ;
        case 4:
            return std::asinf( ((-1)*pz.ycoord) / dist_z ) ;
        default:
            std::cerr<<"Function KMeans::find_theta failed to find theta!\n" ;
            exit( EXIT_FAILURE ) ;
    }
}
//....................................................................................................
void KMeans::place_c1( double r , int pz_quad , Centroid & c1 , double theta  ){
    
    switch (pz_quad){
        case 1:
            c1.xcoord = r * std::cosf( theta ) ;
            c1.ycoord = r * std::sinf( theta ) ;
            break ;
        case 2:
            c1.xcoord = r * ((-1) * std::cosf( theta )) ;
            c1.ycoord = r * std::sinf( theta ) ;
            break ;
        case 3:
            c1.xcoord = r * ((-1) * std::cosf( theta )) ;
            c1.ycoord = r * ((-1) * std::sinf( theta )) ;
            break ;
        case 4:
            c1.xcoord = r * std::cosf( theta ) ;
            c1.ycoord = r * ((-1) * std::sinf( theta )) ;
            break ;
        default:
            std::cerr<<"Function KMeans::place_c1 failed to place c1!\n" ;
            exit( EXIT_FAILURE ) ;
    }
}
//....................................................................................................
double KMeans::find_c1_angle( int pz_quad , double theta  ){
    
    switch (pz_quad){
        case 1:
            return theta ;
        case 2:
            return M_PI - theta ;
        case 3:
            return M_PI + theta ;
        case 4:
            return (2*M_PI) - theta ;
        default:
            std::cerr<<"Function KMeans::find_c1_angle failed!\n" ;
            exit( EXIT_FAILURE ) ;
    }
}
//....................................................................................................
/* Naming notation:
        c1 = ventroid c_i
        c2 = centroid c_i+1
        c1_angle = the angle of centroid c_i
        c2_angle = the angle of centroid c_i+1
        c2_quad = the quadrant in which centroid c_i+1 is located
        c2_theta = is the theta angle of centroid c_i+1 (i.e. the angle required to calculate the
                   c_i+1's x,y-coordinates)      */
double KMeans::place_ci_centroid( Centroid & c2 , int pz_quad , double c1_angle , double beta , double r ){
    
    // find the c2's angle i.e. 'c2_angle' (in radians):
    double c2_angle = std::fmod( (c1_angle + beta) , (2*M_PI) ) ;

    /* determine the Quadrant for c2 using 'c2_angle' */
    /* calculate c2's theta-θ angle (i.e. the angle required to calculate the c2's x,y-coordinates)*/
    /* finally calculate the c2's x,y-coordinates: */
    
    /* WARNING!  switch case with number ranges (i.e. k ... m) is NOT part of the standard C or C++.
    It is provided as as extension of the GNU C compiler -- so compile using GNU */
    int c2_angle_int_deg = static_cast<int>( c2_angle * (180.0/M_PI) ) ; // in degrees
    double c2_theta = 0.0 ;
    
    switch ( c2_angle_int_deg ) {
        case 0 ... 89:              // Quadrant 1
            c2_theta = c2_angle ;
            c2.xcoord = r * std::cosf( c2_theta ) ;
            c2.ycoord = r * std::sinf( c2_theta ) ;
            return c2_angle ;
            
        case 90 :
            c2.xcoord = 0.0 ; c2.ycoord = r ; return c2_angle ;
            
        case 91 ... 179:            // Quadrant 2
            c2_theta = M_PI - c2_angle ;
            c2.xcoord = r * ((-1) * std::cosf( c2_theta )) ;
            c2.ycoord = r * std::sinf( c2_theta ) ;
            return c2_angle ;
            
        case 180:
            c2.xcoord = (-1)*r ; c2.ycoord = 0.0 ; return c2_angle ;
            
        case 181 ... 269:           // Quadrant 3
            c2_theta = c2_angle - M_PI ;
            c2.xcoord = r * ((-1) * std::cosf( c2_theta )) ;
            c2.ycoord = r * ((-1) * std::sinf( c2_theta )) ;
            return c2_angle ;
            
        case 270:
            c2.xcoord = 0.0 ; c2.ycoord = (-1)*r ; return c2_angle ;
            
        case 271 ... 359:           // Quadrant 4
            c2_theta = (2*M_PI) - c2_angle ;
            c2.xcoord = r * std::cosf( c2_theta ) ;
            c2.ycoord = r * ((-1) * std::sinf( c2_theta )) ;
            return c2_angle ;
            
        case 360:
            c2.xcoord = r ; c2.ycoord = 0.0 ; return c2_angle ;
            
        default:
            std::cerr<<"Function KMeans::place_ci_centroid failed to determine the Quadrant!\n" ;
            exit( EXIT_FAILURE ) ;
    }
}
//....................................................................................................
void swap( double & x , double & y  ){
    double temp = x ;
    x = y ;
    y = temp ;
}
//....................................................................................................
void KMeans::adjust_coords( Centroid & c0 , Centroid & ci , int ci_angle ){
    
    double temp = 0.0 ;
    switch ( ci_angle ) {
        case 0 ... 89:              // Quadrant 1
            temp = ci.xcoord ;
            ci.xcoord = c0.xcoord - ci.ycoord ;
            ci.ycoord = c0.ycoord + temp ;
            return ;
            
        case 90 :
            swap( ci.xcoord , ci.ycoord  ) ; return ;
            
        case 91 ... 179 :            // Quadrant 2
            temp = ci.xcoord ;
            ci.xcoord = c0.xcoord - ci.ycoord ;
            ci.ycoord = c0.ycoord - (-1*temp) ;
            return ;
            
        case 180:
            swap( ci.xcoord , ci.ycoord  ) ; return ;
            
        case 181 ... 269 :           // Quadrant 3
            temp = ci.xcoord ;
            ci.xcoord = c0.xcoord + (-1*ci.ycoord) ;
            ci.ycoord = c0.ycoord + temp ;
            return ;
            
        case 270 :
            swap( ci.xcoord , ci.ycoord  ) ; return ;
            
        case 271 ... 359 :           // Quadrant 4
            temp = ci.xcoord ;
            ci.xcoord = c0.xcoord + (-1*ci.ycoord)  ;
            ci.ycoord = c0.ycoord + temp ;
            return ;
            
        case 360 :
            swap( ci.xcoord , ci.ycoord  ) ; return ;
            
        default:
            std::cerr << "Function KMeans::adjust_coords failed to adjust the coordinates!\n" ;
            exit( EXIT_FAILURE ) ;
    }
}
//....................................................................................................
/*                          "IntelliSeed" alg by Martin Gregory Sendrowicz
                            ©Copyright by Martin Gregory Sendrowicz 2023
                                        All Rights Reserved                  */
void KMeans::intelliSeed( double delta ){
    
    /* STEP-1:    consider ALL the data-points and find their collective centroid c0 */
    int sum_x=0, sum_y=0 ;
    for( int p=0 ; p < (*this).num_pts ; ++p ){
        
        sum_x += (*this).point_set[ p ].xcoord ;
        sum_y += (*this).point_set[ p ].ycoord ;
    }
    
    double coord_x = static_cast<double>(sum_x) / static_cast<double>((*this).num_pts) ;
    double coord_y = static_cast<double>(sum_y) / static_cast<double>((*this).num_pts) ;
    
    Centroid centroid( coord_x , coord_y , (*this)._K+1 ) ;
    
    std::cout << "\nThe Collective Centroid is: " ; centroid.showCentroid() ;
    
    /* STEP-2:      find point pz such that pz is the farthest point away from 'centroid c0'. Let’s
    call this distance 'dist_z' and let 'p_idx' be the index of the farthest Point in 'point_set' */
    double dist_z=0.0 ; int p_idx=0 ;
    
    double temp=0.0 ;
    for( int p=0 ; p < (*this).num_pts ; ++p ){
        
        temp = computeDist( (*this).point_set[ p ] , centroid ) ;
        
        if( compare_dbl( dist_z,temp) )
            if( dist_z < temp ){
                dist_z = temp ;
                p_idx = p ;
            }
    }
    std::cout << "The farthest Point is: " ; (*this).point_set[ p_idx ].showPoint() ;
    std::cout << "The dist between farthest Point and Centroid is: "<< dist_z <<std::endl ;
    
    /* STEP-3"  use c0 as the center of the circle with radius r = δz  where δ = e.g. 0.5 */
    double r = delta * dist_z ;
    std::cout << "radius r is: " << r << std::endl ;
    
    /* STEP-4   place 1st centroid c1 at δ z distance away c0 at the intersection of the
    circumference of the circle and the segment joining c0 and pz*/
    
    /* Step-4.1
    Perform the  Coordinate System Transform  from ‘row,col-coordinate’ Image Space into ‘x,y-coordinate’
    Cartesian Space. To adjust our image’s ‘row,col-coordinate’ system to be congruent with the
    ‘x-y-coordinate’ system we need to flip the coordinates; i.e.  row → y  and  col → x */
    Centroid c0 ;
    c0.xcoord = centroid.ycoord ;
    c0.ycoord = centroid.xcoord ;
    c0.label = 0 ;
    
    Point pz ;
    pz.xcoord = (*this).point_set[ p_idx ].ycoord ;
    pz.ycoord = (*this).point_set[ p_idx].xcoord ;
    
    /* Step-4.2 Find pz’s Quadrant (there are 4 cases) */
    int pz_quad = find_pz_quadrant( pz,c0 ) ;
    std::cout << "pz is in quadrant: " << pz_quad << std::endl ;
    
    /* Step-4.2
    Adjust the coordinates of c0 so that c0 from position (h,k) translates to position: (c0h=0, c0k=0).
    Then, based on this translation, find the corresponding position for pz. */
    adjust_c0pz_to_zero( pz, c0, pz_quad ) ;
    std::cout << "pz is @ location: " << pz.xcoord<<" : "<<  pz.ycoord << std::endl ;
    
    /* Step-4.3
    Find the angle θ. Using Unit Circle Symmetry, consider point pz’ located in Quadrant I. The
    coordinates of such pz’ must be mirrored across the x-axis, if pz is in Quadrant II. OR if pz is
    in Quadrant IV, then the mirroring would be across y-axis; and if pz is in Quadrant III, the
    mirroring would be across the origin.
    Also, note that, to find θ, we only need the y-coordinate.
    *Exception: pz may happen to be in Quadrant I—in that case, just calculate angle θ . */
    double theta = find_theta( pz , pz_quad , dist_z ) ;
    std::cout << "Theta is: " << theta * (180.0/M_PI) << std::endl ;     // show in degrees°
    
    /* Step-4.4
    Find placement for c1 centroid . To place c1, we must scale the dist_z using parameter δ (delta).
    Note that c1 is always in the same Quadrant as pz (since c1 is based on pz). */
    Centroid c1 ;
    c1.label = 1 ;
    place_c1( r , pz_quad , c1 , theta  ) ;
    std::cout << "Centroid c1 is at: " << c1.xcoord<<" : "<< c1.ycoord << std::endl ;
    
    /* Step-4.5
    Find the central angle β  by dividing the circumference into k equal parts (2π by k) or (360° by k)
    — based on which we’ll place the other k–1 centroids. */
    double beta = static_cast<double>(2*M_PI) / static_cast<double>((*this)._K) ;
    std::cout << "Cental Angle Beta is: " << beta * (180.0/M_PI) << std::endl ;
    
    /* Step-4.6
    Find the c1’s angle c1° */
    double c1_angle = find_c1_angle( pz_quad , theta  ) ;
    std::cout << "Centroid c1 angle is: " << c1_angle * (180.0/M_PI) << std::endl ;
    
    /* Step-5.1
    Find the remaining k–1 centroids. Starting from centroid c1’s quadrant, determine the quadrant
    in which to place centroids ci, where 2 ≤ i ≤ k) (i.e. ALL the remaining centroids) */
    std::vector<Centroid> centroids ; centroids.emplace_back( c1 ) ;
    std::vector<double> angles ; angles.emplace_back( c1_angle ) ;
    std::vector<int> angles_deg ; angles_deg.emplace_back( static_cast<int>( c1_angle * (180.0/M_PI)) ) ;
    
    double c2_angle = c1_angle ;
    int i = 2 ;
    while( i <= (*this)._K ) {
        Centroid c2 ;
        c2.label = i ;
        
        c2_angle = place_ci_centroid( c2 , pz_quad , c2_angle , beta , r ) ;
        
        angles.emplace_back( c2_angle ) ;
        angles_deg.emplace_back( static_cast<int>( c2_angle * (180.0/M_PI)) ) ;
        
        std::cout << "Centroid c"<<i<<" is at: " << c2.xcoord <<" : "<<c2.ycoord << std::endl ;
        
        centroids.emplace_back( c2 ) ;
        i++ ;
    }
    /* Step-5.2
    Adjust the coordinates  of c1, c2, c3,…,ck so that the origin is back in the Top-Left corner
    (i.e. perform transform to the original image’s coordinate system).
    Note that:      y-coordinates correspond to rows
                    x-coordinates correspond to cols   */
    i = 0 ;
    while( i < (*this)._K ) {
        adjust_coords( centroid , centroids.at(i) , angles_deg.at(i) ) ;
        
        (*this)._Kcentroids[ i ] = centroids.at(i) ;
        
        std::cout << "Centroid c"<<i+1<<" is at: " << centroids.at(i).xcoord <<" : "<<
                  centroids.at(i).ycoord << std::endl ;
        i++ ;
    }
    std::cout << std::endl ;
}
//....................................................................................................
void KMeans::kMeansClustering( std::ofstream & out_iter ) {

    int iteration = 0 ;
    while ( true ) {
        
        /* Iterate over the given 'point-set' and relabel each 'data' point. The 'point-set' contains
        the updated (in the previous iterations) labels */
        point2Image() ;

        /* Output the result of each iteration to 'iterations.txt' file. This is great for visual
        confirmaction */
        printImage ( out_iter , iteration ) ;
        
        /* reset the 'change's of the previous iteration to see if this one will also produce any
        changes -- otherwise you are done. */
        this->change = 0 ;
        
        /* Seed the Centroids -- random seeding leaves the correct clustering to chance and at times
        may fail -- e.g. having the same cluster labeled heterogeneously between two centroids. */
        if( iteration == 0 ){
            //randomSeedCentroids() ;
            intelliSeed( 0.5 ) ;
        } else
            computeCentroids() ;
        
        /* Relabel the data points based on their closest centroid. */
        int index = 0 ;
        while( index < this->num_pts ){

            distanceMinLabel( point_set[ index ] ) ;
            index++ ;
        }
        iteration++ ;
        /* repeat until NO more changes occur */
        if( this->change == 0 ) break ;
    }
}
//....................................................................................................
/* Iterates via the entire 'point_set' array and computes each of the K centroids. The computed
   centroids are stored in '_Kcentroids' array'; where each consecutive centroid's label increments
   from 1 to K.*/
void KMeans::computeCentroids() {
    
    /*STEP 0:   Dynamically allocate and initialize ALL the necessary data structures that are required
                to compute centroids */
    std::unique_ptr< double[] > sumX { new double[ this->_K ] } ;
    std::unique_ptr< double[] > sumY { new double[ this->_K ] } ;
    std::unique_ptr< int[] > totalPt { new int[ this->_K ] } ;
    
    for( int i=0 ; i<(*this)._K ; ++i ){
        sumX[ i ] = 0.0 ;
        sumY[ i ] = 0.0 ;
        totalPt[ i ] = 0 ;
    }
    //STEP 1:   iterate via ALL the Points in 'point_set' and extract the necessary information
    int index=0, label=0 ;
    while( index < (*this).num_pts ){

        label = (*this).point_set[ index ].label ;      // get the given Point’s cluster label
        sumX[ label-1 ] +=  (*this).point_set[ index ].xcoord ;
        sumY[ label-1 ] +=  (*this).point_set[ index ].ycoord ;
        totalPt[ label-1 ]++ ;      /* for each of the k labels, capture how many points have
                                       the given label */
        index++ ;
    }
    //STEP 2:   compute each centroid
    label = 0 ; double x_before = 0.0 ; double y_before = 0.0 ;
    while ( label < this->_K ) {
        
        x_before = (*this)._Kcentroids[ label ].xcoord ;
        y_before = (*this)._Kcentroids[ label ].ycoord ;

        (*this)._Kcentroids[ label ].xcoord = sumX[ label ] / static_cast<double>(totalPt[ label ])  ;
        (*this)._Kcentroids[ label ].ycoord = sumY[ label ] / static_cast<double>(totalPt[ label ])  ;
        
        if( compare_dbl( x_before , (*this)._Kcentroids[ label ].xcoord ) )
            this->change++ ;
        if( compare_dbl( y_before , (*this)._Kcentroids[ label ].ycoord ) )
            this->change++ ;
        
        (*this)._Kcentroids[ label ].showCentroid() ;

        label++ ;
    }
}
//....................................................................................................
/* Computes the distance from the given Point p to each of the K centroids. Then checks if point's
   label needs to be adjusted (or not) */
void KMeans::distanceMinLabel ( Point & point ) {
    
    /*STEP 0:   'min_dist' i.e. distance to the nearest centroid
                'min_lbl' is the label that matches the label of the nearest centroid */
    double min_dist = __DBL_MAX__ ;
    int min_lbl = 0 ;
    
    //STEP 1:   compute k distances to each of the k centroids
    int label=1 ; double dist=0.0 ;
    while ( label <= (*this)._K ) {
        
        dist = computeDist( point , _Kcentroids[ label-1 ] ) ;
        
    //STEP 2:   find the min distance and its corresponding label. I.e. if(mid_dist > dist)...
        if( compare_dbl( min_dist, dist) )
            if( min_dist > dist) {
                min_dist = dist ;
                min_lbl = label ;
            }
        label++ ;
    }
    /* Compare the 'min_dist' against the current point's distance and see if it needs to be adjusted.
    If the adjustement occurs, increment the 'change' data member so that another iteration will be
    in order. I.e. if(point.distance > min_dist)... */
    /* Important Update: sometimes the distance to a different centroid (which has already moved away) could
    have been smaller, now the distance to that very centroid has increased and a new centroid has
    moved in. Unfortunately the distance to that new centroid is larger than the previously set min.
    In that scenario the labels will NOT update accordingly. Solution: check if the previously set
    label is different that the new min_lbl. If it is, update to the new larger min.*/
//    if( compare_dbl( point.distance, min_dist) )
//        if( point.distance > min_dist || point.label != min_lbl ){
//            point.distance = min_dist ;
//            point.label = min_lbl ;
//            this->change++ ;
//        }
    // OR even simpler ...
     if( point.label != min_lbl ){
             point.distance = min_dist ;
             point.label = min_lbl ;
             this->change++ ;
     }
}
//...................................................................................↑↑↑ helper ↑↑↑...
/* Computes the Euclidean distance from point p1 to the given centroid p2. */
double KMeans::computeDist( Point p1, Centroid p2) {
    
    double x = std::powf( p1.xcoord - p2.xcoord,2 ) ;
    double y = std::powf( p1.ycoord - p2.ycoord,2 ) ;
    return std::sqrtf( x + y ) ;
}
//....................................................................................................
/* Reads each Point in 'point_set' and writes the corresponding point's label into 'data' 2D array. */
void KMeans::point2Image () {
      
    for(int i=0 ; i<(*this).num_pts ; ++i)
        this->data  [ static_cast<int>((*this).point_set[i].xcoord) ]
                    [ static_cast<int>((*this).point_set[i].ycoord) ] = (*this).point_set[ i ].label ;
}
//....................................................................................................
/* Outputs each iteration to show how the data point's labels adjust according to the traveling
   centroids. The centroids are demarcated as 'X'. */
void KMeans::printImage ( std::ofstream & out_iter, int iteration){
    
    out_iter << "\n*** Result of iteration "<< iteration <<" ***" ;
    if( iteration==1 )
        out_iter << " Centroids placed by IntelliSeed Algorithm ***" ;
    
    std::unique_ptr< std::unique_ptr< int[] >[] > temp ;
    temp = std::unique_ptr< std::unique_ptr< int[] >[] >( new std::unique_ptr<int[]>[ (*this).num_rows ] ) ;
    
    for( int row=0 ; row<(*this).num_rows ; row++ ){
        temp[row] = std::unique_ptr< int[] >( new int[ (*this).num_cols ] ) ;
        for( int col=0 ; col<(*this).num_cols ; col++ )
            temp[row][col] = data[row][col] ;
    }
    int x=0 , y=0 , l=(*this)._K+1 ;
    for( int k=0 ; k<(*this)._K ; ++k ){
        x = (*this)._Kcentroids[ k ].xcoord ;
        y = (*this)._Kcentroids[ k ].ycoord ;
        temp[x][y] = l++  ;                 // distinctly demarcate centroid locations as K+1 label
    }
    for( int r=0 ; r<this->num_rows ; ++r ) {
        out_iter << std::endl ;
        
        for( int c=0 ; c<this->num_cols ; ++c ) {
            if( temp[r][c]==0 )
                out_iter<< "." << " " ;
            else {
                if( temp[r][c] < (*this)._K+1 )
                    out_iter << temp[r][c] << " " ;
                else
                    out_iter << "X" << " " ;    // show where the centroids are
            }
        }
    }
}
//....................................................................................................
/* Writes each Point in 'point_set' with its label and xy location into the final output.txt file.
   This step is for future processing -- i.e. the data is now clustered and can be used for other
   projects .*/
void KMeans::writePtSet( std::ofstream & output ) {

    output << "Number of Data Points: "<< num_pts << std::endl ;
    //output << "Rows:" << num_rows << " Cols:" << num_cols << std::endl ;

    for( int i=0 ; i<(*this).num_pts ; ++i ){
        output << point_set[i].xcoord <<" "<< point_set[i].ycoord <<" "<< point_set[i].label << std::endl ;
    }
}
//....................................................................................................

