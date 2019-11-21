//Inmersed Boundary Wettability (Based on Hourge, P. 2004)

Nalpha = fvc::grad(alpha1)/mag(fvc::grad(alpha1));

Nwall = fvc::grad(B)/mag(fvc::grad(B));

thetaI = acos( (Nalpha & Nwall)/(mag(Nalpha)*mag(Nwall)) );

NCorr = Nwall*(cos(thetaEq)-cos(thetaI)*cos(thetaI-thetaEq))/(1-cos(thetaI)^2)
      + Nalpha*(cos(thetaI - thetaEq)-cos(thetaI)*cos(thetaEq))/(1-cos(thetaI)^2);
