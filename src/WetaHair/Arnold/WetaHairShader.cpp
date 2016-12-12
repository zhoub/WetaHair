#include "WetaHair.hpp"

AI_SHADER_NODE_EXPORT_METHODS(WetaHairMethods);

enum WETA_HAIR_PARAMETERS
{
    MU_A ,
    ETA ,

    WEIGHT_R ,
    BETA_R ,

    WEIGHT_TT ,
    BETA_TT ,

    WEIGHT_TRT ,
    BETA_TRT ,

    WEIGHT_TRRT ,
    BETA_TRRT
};

node_parameters
{
    AiParameterRGB( "MuA" ,  0.1f ,  0.2f , 0.8f ) ;

    AiParameterFlt( "Eta" ,  1.55f );

    //
    AiParameterFlt( "WeightR" , 0.0f ) ;

    AiParameterFlt( "BetaR" , 8.0f );

    //
    AiParameterFlt( "WeightTT" , 0.0f ) ;

    AiParameterFlt( "BetaTT" , 8.0f );

    //
    AiParameterFlt( "WeightTRT" , 0.0f ) ;

    AiParameterFlt( "BetaTRT", 8.0f );

    //
    AiParameterFlt( "WeightTRRT" , 0.0f ) ;

    AiParameterFlt( "BetaTRRT" , 8.0f );

    //
    AiMetaDataSetInt( mds , NULL , "maya.id" , 0x20131122 );
	AiMetaDataSetStr( mds , NULL , "maya.classification" , "shader/surface" );
	AiMetaDataSetStr( mds , NULL , "maya.output_name" , "outColor" );
	AiMetaDataSetStr( mds , NULL , "maya.output_shortname" , "out" );
}

node_initialize
{
}

node_update
{
}

node_finish
{
}

node_loader
{
    switch ( i )
    {
        case 0 :
        {
            node->methods = WetaHairMethods ;
            node->output_type = AI_TYPE_RGB ;
            node->name = "WetaHair" ;
            node->node_type = AI_NODE_SHADER ;
            strcpy( node->version , AI_VERSION ) ;
            
            return TRUE;
        }
    }

    return FALSE;
}

shader_evaluate
{
    if( sg->Rt & AI_RAY_SHADOW )
	{
		return;
	}

    //
    AtColor cMuA = AiShaderEvalParamRGB( MU_A ) ;

    AtFloat fEta = AiShaderEvalParamFlt( ETA ) ;

    AtFloat fWeightR = AiShaderEvalParamFlt( WEIGHT_R ) ;
    AtFloat fBetaR = AI_DTOR * AiShaderEvalParamFlt( BETA_R ) ;

    AtFloat fWeightTT = AiShaderEvalParamFlt( WEIGHT_TT ) ;
    AtFloat fBetaTT = AI_DTOR * AiShaderEvalParamFlt( BETA_TT ) ;

    AtFloat fWeightTRT = AiShaderEvalParamFlt( WEIGHT_TRT ) ;
    AtFloat fBetaTRT  = AI_DTOR * AiShaderEvalParamFlt( BETA_TRT ) ;

    AtFloat fWeightTRRT = AiShaderEvalParamFlt( WEIGHT_TRRT ) ;
    AtFloat fBetaTRRT = AI_DTOR * AiShaderEvalParamFlt( BETA_TRRT ) ;

    //
    AtVector vV = AiV3Normalize( - sg->Rd ) ;

    AtVector vT = AiV3Normalize( sg->dPdv ) ;
    AtVector vX = AiV3Normalize( AiV3Cross( vT , vV ) ) ;
    AtVector vY = AiV3Normalize( AiV3Cross( vX , vT ) ) ;

    AtFloat fThetaR = AI_PIOVER2 - acosf( AiV3Dot( vV , vT ) ) ;
    AtFloat fPhiR   = atan2f( AiV3Dot( vV , vX ) , AiV3Dot( vV , vY ) ) ;

    //
    AtColor cDirectR = AI_RGB_BLACK , cDirectTT = AI_RGB_BLACK , cDirectTRT = AI_RGB_BLACK , cDirectTRRT = AI_RGB_BLACK ;

    sg->fhemi = FALSE;
    AiLightsPrepare( sg ) ;
    while ( AiLightsGetSample( sg ) )
    {
        if ( AiLightGetAffectSpecular( sg->Lp ) )
        {
            AtVector vL = AiV3Normalize( - sg->Ld ) ;
            AtFloat fThetaI = AI_PIOVER2 - acosf( AiV3Dot( vL , vT ) ) ;
            AtFloat fPhiI = atan2f( AiV3Dot( vL , vX ) , AiV3Dot( vL , vY ) ) ;

            AtFloat fThetaD = ( fThetaR - fThetaI ) * 0.5f ;

            AtFloat fPhi = fPhiR - fPhiI ;
            if ( fPhi > AI_PI )
            {
                fPhi -= AI_PITIMES2 ;
            }
            if ( fPhi < - AI_PI ) 
            {
                fPhi += AI_PITIMES2 ;
            }

            if ( fWeightR > 0.0f )
            {
                cDirectR += sg->Li
                          * M( fBetaR , fThetaI , fThetaR )
                          * N( 0 , fBetaR , fPhi , fEta , EtaPrime( fEta , fThetaD ) , AiV3Dot( vV , vL ) , cMuA ) ;
            }

            if ( fWeightTT > 0.0f )
            {
                cDirectTT += sg->Li
                           * M( fBetaTT , fThetaI , fThetaR )
                           * N( 1 , fBetaTT , fPhi , fEta , EtaPrime( fEta , fThetaD ) , fThetaD , cMuA ) ;
            }

            if ( fWeightTRT > 0.0f )
            {
                cDirectTRT += sg->Li
                            * M( fBetaTRT , fThetaI , fThetaR )
                            * N( 2 , fBetaTRT , fPhi , fEta , EtaPrime( fEta , fThetaD ) , fThetaD , cMuA ) ;
            }
 
            if ( fWeightTRRT > 0.0f )
            {
                cDirectTRRT += sg->Li
                             * M( fBetaTRT , fThetaI , fThetaR )
                             * N( 3 , fBetaTRRT , fPhi , fEta , EtaPrime( fEta , fThetaD ) , fThetaD , cMuA ) ;
            }
        }
    }

    sg->out.RGB = fWeightR * cDirectR
                + fWeightTT * cDirectTT
                + fWeightTRT * cDirectTRT
                + fWeightTRRT * cDirectTRRT ;
}
