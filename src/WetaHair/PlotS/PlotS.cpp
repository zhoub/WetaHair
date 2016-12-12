/*
 * The MIT License (MIT)
 * Copyright (c) 2011-2016 Bo Zhou<bo.schwarzstein@gmail.com>
 * Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
 * The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#include <omp.h>

#include <cmath>
#include <iostream>
#include <vector>

#include <boost/timer/timer.hpp>

#include <WetaHair/Implementation/Origin.hpp>

#include "rgbe.h"

int main(int iArgc, char * aArgv[])
{
    -- iArgc , ++ aArgv ;
    if ( iArgc != 5 )
    {
        return EXIT_FAILURE ;
    }

    int iP = atoi( aArgv[0] ) ;
    double lfBeta = atof( aArgv[1] ) * boost::math::constants::degree<double>() ;
    double lfEta = atof( aArgv[2] ) ;
    double aMuA[3] = { 0.42f , 0.70f , 1.37f } ;

    int iSize = atoi( aArgv[3] ) ;
    FILE * pOutput = fopen( aArgv[4] , "wb" ) ;
    if ( ! pOutput )
    {
        return EXIT_FAILURE ;
    }
    RGBE_WriteHeader( pOutput , iSize , iSize , NULL ) ;

    //
    std::vector<float> vPixels;
    vPixels.resize( iSize * iSize * 3 ) ;

    double lfDelta = boost::math::constants::two_pi<double>() / iSize ;

    boost::timer::auto_cpu_timer cTimer ;
    #pragma omp parallel for
    for ( int y = 0 ; y < iSize ; ++ y )
    {
        double lfThetaD = - boost::math::constants::pi<double>() + lfDelta * ( y + 0.5 ) ;
        double lfEtaPrime = EtaPrime( lfEta , lfThetaD ) ;
        for ( int x = 0 ; x < iSize ; ++ x )
        {
            double aResult[3] = { 0.0 , 0.0 , 0.0 } ;
            if ( iP == 0 )
            {
                // $ \acos \in [ 0 , \pi ] $
                double lfPhi = - boost::math::constants::pi<double>() + lfDelta * ( x + 0.5 ) ;
                double lfAngle = cos( lfDelta * 0.5 * ( x + 0.5 ) ) ;
                N( 3 , aResult, aMuA , iP , lfBeta , lfPhi , lfEta , lfEtaPrime , lfAngle ) ;
            }
            else
            {
                // $\phi \in [ - \pi , pi ]$
                double lfPhi = - boost::math::constants::pi<double>() + lfDelta * ( x + 0.5 ) ;
                N( 3 , aResult, aMuA , iP , lfBeta , lfPhi , lfEta , lfEtaPrime , lfThetaD ) ;
            }
            int z = ( x + y * iSize ) * 3 ;
            vPixels[z + 0] = aResult[0] ;
            vPixels[z + 1] = aResult[1] ;
            vPixels[z + 2] = aResult[2] ;
        }
    }

    RGBE_WritePixels( pOutput , reinterpret_cast<float *>( & vPixels[0] ) , vPixels.size() );
    fclose( pOutput ) ;

    return EXIT_SUCCESS ;
}
