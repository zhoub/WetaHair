/*
 * The MIT License (MIT)
 * Copyright (c) 2011-2016 Bo Zhou<bo.schwarzstein@gmail.com>
 * Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
 * The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#include <iostream>
#include <cmath>

#include <WetaHair/Implementation/Origin.hpp>

int main(int iArgc, char * aArgv[])
{
    double lfBeta = 8.1 * boost::math::constants::degree<double>() ;

    double aThetaI[5] = { - 1.0 , - 1.1  , - 1.2 , - 1.3 , - 1.4 } ;

    static const double lfDelta = boost::math::constants::half_pi<double>() / 1024.0 ;
    for ( int i = 0 ; i < 5 ; ++ i )
    {
        double lfThetaI = aThetaI[i] ;
        for ( double lfThetaR = 0.0 ; lfThetaR <= boost::math::constants::half_pi<double>() ; lfThetaR += lfDelta )
        {
            std::cout << std::fixed << M( lfBeta , lfThetaI ,  lfThetaR ) << ',' ;
        }
        std::cout << std::endl ;
    }

    return EXIT_SUCCESS ;
}
