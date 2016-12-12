/*
 * The MIT License (MIT)
 * Copyright (c) 2011-2016 Bo Zhou<bo.schwarzstein@gmail.com>
 * Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
 * The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#pragma once

#include <cmath>

#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions.hpp>

static inline double G( double lfBeta , double lfX)
{
    double lfResult = exp( - lfX * lfX / ( 2.0 * lfBeta * lfBeta ) )
                    / ( sqrt( boost::math::constants::two_pi<double>() ) * lfBeta ) ;
    return lfResult ;
}

///////////////////////////////////////////////////////////////////////////////////////////////////

static inline double M( double lfBeta , double lfThetaI , double lfThetaR )
{
    double lfV = lfBeta * lfBeta ;
    double lfResult = 1 / sinh( 1.0 / lfV ) / ( 2.0 * lfV )
                      * exp( sin(- lfThetaI ) * sin( lfThetaR ) / lfV )
                      * boost::math::cyl_bessel_i( 0 , cos( - lfThetaI ) * cos( lfThetaR ) / lfV ) ;
    return lfResult ;
}

///////////////////////////////////////////////////////////////////////////////////////////////////

static inline double F( double lfEta , double lfX )
{
    double lfR = pow( ( 1.0 - lfEta ) / ( 1.0 + lfEta) , 2.0 ) ;
    double lfResult = lfR + ( 1.0 - lfR ) * pow( 1.0 - cos(lfX) , 5.0 ) ;
    return lfResult ;
}

static inline double T( double lfMuA , double lfGammaT )
{
    double lfCosGammaT = cos( lfGammaT ) ;
    double lfMuAPrime = lfMuA / lfCosGammaT ;
    double lfResult = exp( - 2.0 * lfMuAPrime * ( 1.0 + cos( 2.0 * lfGammaT ) ) ) ;
    return lfResult ;
}

static inline double GammaI( double lfH )
{
    double lfResult = asin( lfH ) ;
    return lfResult ;
}

static inline double GammaT( double lfH , double lfEtaPrime )
{
    double lfResult = asin( lfH / lfEtaPrime ) ;
    return lfResult ;
}

static inline double EtaPrime( double lfEta , double lfThetaD )
{
    double lfSinThetaD = sin( lfThetaD ) ;
    double lfCosThetaD = cos( lfThetaD ) ;
    double lfResult = sqrt( lfEta * lfEta - lfSinThetaD * lfSinThetaD )
                      / lfCosThetaD ;
    return lfResult ;
}

static inline double Phi( int iP , double lfH , double lfEtaPrime )
{
    double lfResult = 2.0 * iP * GammaT( lfH , lfEtaPrime )
                    - 2.0 * GammaI( lfH )
                    + iP * boost::math::constants::pi<double>() ;
    return lfResult ;
}

static inline double D( double lfBeta , double lfPhi , const int iK = 32)
{
    double lfResult = 0.0 ;
    for ( int k = - iK ; k < iK ; ++ k )
    {
        lfResult += G( lfBeta , lfPhi - boost::math::constants::two_pi<double>() * k ) ;
    }
    return lfResult;
}

static inline double A( int iP , double lfH , double lfEta , double lfEtaPrime , double lfAngle , double lfMuA )
{
    double lfResult = 0.0 ;
    if ( iP == 0 )
    {
        lfResult = F( lfEta , 0.5 * acos( lfAngle ) ) ;
    }
    else
    {
        double lfF = F( lfEta , acos( cos( lfAngle ) * cos( asin( lfH ) ) ) ) ;
        double lfGammaT = GammaT( lfH , lfEtaPrime ) ;
        lfResult = pow( 1.0 - lfF , 2.0 ) * pow( lfF , iP - 1.0 ) * T( lfMuA , lfGammaT ) ;
    }
    return lfResult ;
}

static void N( int iLength , double * aResult , const double * aMuA , int iP , double lfBeta , double lfPhi , double lfEta , double lfEtaPrime , double lfAngle )
{
    static const double aAbscissas[32] = { 0.0243502926634244325089558 , 0.0729931217877990394495429 , 0.1214628192961205544703765 , 0.1696444204239928180373136 , 0.2174236437400070841496487 , 0.2646871622087674163739642 , 0.3113228719902109561575127 , 0.3572201583376681159504426 , 0.4022701579639916036957668 , 0.4463660172534640879849477 , 0.4894031457070529574785263 , 0.5312794640198945456580139 , 0.5718956462026340342838781 , 0.6111553551723932502488530 , 0.6489654712546573398577612 , 0.6852363130542332425635584 , 0.7198818501716108268489402 , 0.7528199072605318966118638 , 0.7839723589433414076102205 , 0.8132653151227975597419233 , 0.8406292962525803627516915 , 0.8659993981540928197607834 , 0.8893154459951141058534040 , 0.9105221370785028057563807 , 0.9295691721319395758214902 , 0.9464113748584028160624815 , 0.9610087996520537189186141 , 0.9733268277899109637418535 , 0.9833362538846259569312993 , 0.9910133714767443207393824 , 0.9963401167719552793469245 , 0.9993050417357721394569056 } ;
    static const double aWeights[32]   = { 0.0486909570091397203833654 , 0.0485754674415034269347991 , 0.0483447622348029571697695 , 0.0479993885964583077281262 , 0.0475401657148303086622822 , 0.0469681828162100173253263 , 0.0462847965813144172959532 , 0.0454916279274181444797710 , 0.0445905581637565630601347 , 0.0435837245293234533768279 , 0.0424735151236535890073398 , 0.0412625632426235286101563 , 0.0399537411327203413866569 , 0.0385501531786156291289625 , 0.0370551285402400460404151 , 0.0354722132568823838106931 , 0.0338051618371416093915655 , 0.0320579283548515535854675 , 0.0302346570724024788679741 , 0.0283396726142594832275113 , 0.0263774697150546586716918 , 0.0243527025687108733381776 , 0.0222701738083832541592983 , 0.0201348231535302093723403 , 0.0179517157756973430850453 , 0.0157260304760247193219660 , 0.0134630478967186425980608 , 0.0111681394601311288185905 , 0.0088467598263639477230309 , 0.0065044579689783628561174 , 0.0041470332605624676352875 , 0.0017832807216964329472961 } ;

    for ( int i = 0 ; i < 32 ; ++ i )
    {
        double lfH = aAbscissas[i] ;
        double lfWeight = aWeights[i] ;
        for ( int j = 0 ; j < iLength ; ++ j )
        {
            aResult[j] += lfWeight * ( A( iP ,   lfH , lfEta , lfEtaPrime , lfAngle , aMuA[j] ) * D( lfBeta , lfPhi - Phi( iP ,   lfH , lfEtaPrime ) )
                                     + A( iP , - lfH , lfEta , lfEtaPrime , lfAngle , aMuA[j] ) * D( lfBeta , lfPhi - Phi( iP , - lfH , lfEtaPrime ) ) ) ;
        }

    }

    for ( int i = 0 ; i < iLength ; ++ i )
    {
        aResult[i] *= 0.5 ;
    }
}
