function hdilim = hdi(samplevec,credmass)
% hdilim = hdi(samplevec,credmass)
% Computes highest density interval from a sample of representative values,
% estimated as shortest credible interval.
%
% INPUTS:
%     samplevec = vector of representative values from a probability distribution.
%     credmass  = scalar between 0 and 1, indicating the mass within the credible
%                 interval that is to be estimated, e.g. 0.90.
% OUTPUT:
%     hdilim    = vector containing the limits of the hdi
%
% HDI implementation based on original R code HDIofMCMC from John K. Kruschke:
% https://github.com/boboppie/kruschke-doing_bayesian_data_analysis/blob/master/1e/HDIofMCMC.R
% ------------------------------------------
% Copyright (C) Guillaume Rousselet 2015

% GAR - University of Glasgow - 16 Dec 2015

sortedPts = sort( samplevec );
ciIdxInc = floor( credmass * length( sortedPts ) );
nCIs = length( sortedPts ) - ciIdxInc;
ciWidth = zeros(nCIs,1);
for ci = 1:nCIs
    ciWidth(ci) = sortedPts(ci + ciIdxInc) - sortedPts(ci);
end

HDImin = sortedPts( find(ciWidth == min(ciWidth),1) );
HDImax = sortedPts( find(ciWidth == min(ciWidth),1)  + ciIdxInc);
hdilim = [HDImin HDImax];

% original R code from
% https://github.com/boboppie/kruschke-doing_bayesian_data_analysis/tree/master/1e
% HDIofMCMC = function( sampleVec , credMass=0.95 ) {
%     # Computes highest density interval from a sample of representative values,
%     #   estimated as shortest credible interval.
%     # Arguments:
%     #   sampleVec
%     #     is a vector of representative values from a probability distribution.
%     #   credMass
%     #     is a scalar between 0 and 1, indicating the mass within the credible
%     #     interval that is to be estimated.
%     # Value:
%     #   HDIlim is a vector containing the limits of the HDI
%     sortedPts = sort( sampleVec )
%     ciIdxInc = floor( credMass * length( sortedPts ) )
%     nCIs = length( sortedPts ) - ciIdxInc
%     ciWidth = rep( 0 , nCIs )
%     for ( i in 1:nCIs ) {
%         ciWidth[ i ] = sortedPts[ i + ciIdxInc ] - sortedPts[ i ]
%     }
%     HDImin = sortedPts[ which.min( ciWidth ) ]
%     HDImax = sortedPts[ which.min( ciWidth ) + ciIdxInc ]
%     HDIlim = c( HDImin , HDImax )
%     return( HDIlim )
% }