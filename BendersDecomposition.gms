*The input file is included first.
$Include C:\Users\Emily\Dropbox\master\gitfiles\master\eurdata.gms

scalar epsi convergence tolerance /1e-3/ ;
scalar epsi2 alternative convergence tolerance /1e-4/ ;

*Scalars For Timing
scalar
starttime
elapsed        Total Wall time
SPstart
SPel           SP wall time
MPstart
MPel           MP wall time
rm             CPU time MP
rs             CPU time SP
;

option profile =3;
option profiletol =0.00001;

set iteration maximum number of iterations /1*1000 / ;
set theIter(iteration) dynamic set of conducted iterations ;
alias (theIter, theIterI);

parameters
xeItPar(m,n,e)                   parameter of transformation expansion values at present iteration
feItPar(m,n,n,e)                 parameter of flow expansion values at present iteration
xePar(m,n,e,iteration)           parameter of transformation expansion values at each iteration
fePar(m,n,n,e,iteration)         parameter of flow expansion values at each iteration
lambdaF(m,n,n,e,iteration)       dual to constraint fixing flow expansions
lambdaX(m,n,e,iteration)         dual to constraint fixing transformation expansions
qsPar(m,n,p,d,e,iteration)       parameter of qs in MP (vi - based)
qpPar(m,n,p,e,iteration)         parameter of qp in MP (vi - based)
fPar(m,n,n,e,iteration)          parameter of f in MP (vi - based)
xPar(m,n,e,e,iteration)          parameter of x in MP (vi - based)
alphaPar(m,n,p,e,iteration)      parameter of alpha in MP (VI - based)
CG(iteration)                    convergence gap
CG2                              parameter keeps track of previous CG during iterations
upsilonPar(m,n,n,e,iteration)    market clareing dual for f
phiPar(m,n,e,e,iteration)        market clareing dual for x
iotaPar(iteration)               parameter for iota to be
pbp                              parameter forn scenario specific probability in SPs
dcp                              parameter forn scenario specific disc. rate in SPs
intp(n,d,e)                      parameter forn scenario specific int in SPs
slpp(n,d,e)                      parameter forn scenario specific slp in SPs
lp(n,e,e)                        parameter forn scenario specific efficiency rate in SPs
maxProdVp(n,p,e)                 parameter forn scenario specific maximum priduction capacity in SPs
fmip(n,n,e)                      parameter forn scenario specific sum of earlier flow expansions in SPs
xmip(n,e)                        parameter forn scenario specific sum of earlier transformation expansions in SPs

;

positive variables

qp2(n,p,e)              quantity produced
qs2(n,p,d,e)            quantity supplied
qt2(n,n,p,e)            quantity transported
qc2(n,p,e,e)            quantity transformed
f2(n,n,e)               flow (sum of transported quantities)
x2(n,e,e)               quantity transformed

fe(m,n,n,e)             investments in flow capacity
xe(m,n,e)               investments in transformation output capacity

delta(m,n,n,e)          dual to maxExF
zeta(m,n,e)             dual to maxExX
eta                     dual to thetaCon
iota(iteration)         dual to newcut

alpha2(n,p,e)           dual to maxProd
gamma2(n,n,e)           dual to maxFlow
epsilon2(n,e)           dual to maxTrans
;

variables
theta                   approximation of SP solution in MP objective
zMP                     optimal value of MP objective
zMP0                    optimal value of initial MP objective
zSP                     optimal value of SP objective
z                       optimal value of total
beta(m,n,p,e)           dual to massBalance
beta2(n,p,e)            dual to massBalance
upsilon(m,n,n,e)        dual to eq marketClearingF (flow) and unit price of transportation in SP
phi(m,n,e,e)            dual to eq marketClearingX
upsilon2(n,n,e)         dual to eq marketClearingF (flow) and unit price of transportation in SP
phi2(n,e,e)             dual to eq marketClearingX
;

Equations

convexityCon            convexity constraint in MP

newcut3                 cuts added in MP VI-based

maxProd2(n,p,e)         production constraint for producer
massBalance2(n,p,e)     mass balance constraint for producer

maxFlow2(n,n,e)         flow constraint for transporter in SP
maxExF(m,n,n,e)         expansion constraint (upper limit) for transporter in MP
maxTrans2(n,ei)         transformation constraint for transformer in SP
maxExX(m,n,e)           expansion constraint for transformer in MP

c12(n,p,d,e)            comp to qs2
c22(n,p,e)              comp to qp2
c32(n,ni,p,e)           comp to qt2
c42(n,p,ei,e)           comp to qc2
c52(n,ni,e)             comp to f2
c6m(m,n,ni,e)           comp to fe
c72(n,ei,e)             comp to x2
c8m(m,n,e)              comp to xe

marketClearingF2(n,n,e) links the producer and transporter problem by assuring concistiency among transported quantities in SP
marketClearingX2(n,e,e) links the producer and transformer problem by assuring concistiency among transformed quantities in SP

;

c12(n,p,d,e)$ spno(n,p)..       - pbp*dcp*(intp(n,d,e) - slpp(n,d,e)*qs2(n,p,d,e) - slpp(n,d,e)*sum(pi,qs2(n,pi,d,e))) - beta2(n,p,e) =G= 0;
*the following option for c12 can be used for monopolies
*c12(m,n,p,d,e)..                  - pb(m)*dc(m)*(int(m,n,d,e) - slp(m,n,d,e)*qs(m,n,p,d,e)) - beta(m,n,p,e) =G= 0;

c22(n,p,e)..                     pbp*dcp*(2*k1*qp2(n,p,e) + k2) + alpha2(n,p,e) + beta2(n,p,e) =G= 0;
c32(n,ni,p,e)..                  pbp*dcp*(upsilon2(n,ni,e)) - beta2(n,p,e) + beta2(ni,p,e) =G= 0;
c42(n,p,ei,e)..                  pbp*dcp*(phi2(n,ei,e)) + lp(n,ei,e)*beta2(n,p,e) - beta2(n,p,ei) =G= 0;
c52(n,ni,e)..                    pbp*dcp*(2*k3*f2(n,ni,e) - upsilon2(n,ni,e)) + gamma2(n,ni,e) =G= 0;
c6m(m,n,ni,e)..                  pb(m)*dc(m)*k4 + delta(m,n,ni,e) - sum(theIter, iota(theIter)*sum(mi $ ancestor(mi,m), pb(m)*dc(m)*lambdaF(mi,n,ni,e,theIter) ) ) =G= 0;
c72(n,ei,e)..                    pbp*dcp*(2*k5*x2(n,ei,e) - phi2(n,ei,e)) + lp(n,ei,e)*epsilon2(n,e) =G= 0;
c8m(m,n,e)..                     pb(m)*dc(m)*k6 + zeta(m,n,e) - sum(theIter, iota(theIter)*sum(mi $ ancestor (mi,m), pb(m)*dc(m)*lambdaX(mi,n,e,theIter) ) ) =G= 0;

maxProd2(n,p,e)..                maxProdVp(n,p,e)   - qp2(n,p,e) =G= 0;
massBalance2(n,p,e)..            - sum(d,qs2(n,p,d,e)) - sum( ni $ (not (ord(ni) = ord(n))), qt2(n,ni,p,e)) - sum( ei $ (not (ord(ei) = ord(e))), qc2(n,p,e,ei)) + qp2(n,p,e) + sum( ni $ (not (ord(ni) = ord(n))) , qt2(ni,n,p,e)) + sum( ei $ (not (ord(ei) = ord(e))), lp(n,ei,e)*qc2(n,p,ei,e)) =E=  0 ;

maxFlow2(n,ni,e)..               maxFlowV(n,ni,e) + fmip(n,ni,e) - f2(n,ni,e) =G=  0 ;
maxExF(m,n,ni,e)..               maxExFV(n,ni,e) - fe(m,n,ni,e) =G= 0 ;

maxTrans2(n,e)..                 maxTransV(n,e) + xmip(n,e) - sum(ei  $ (not (ord(e) = ord(ei))) , lp(n,ei,e)*x2(n,ei,e)) =G=0 ;
maxExX(m,n,e)..                  maxExXV(n,e) - xe(m,n,e) =G= 0 ;

marketClearingF2(n,ni,e)..       - f2(n,ni,e) + sum( p, qt2(n,ni,p,e)) =E= 0 ;
marketClearingX2(n,e,ei)..       - x2(n,e,ei) + sum( p, qc2(n,p,e,ei)) =E= 0 ;

convexityCon.. 1 - sum(theIter, iota(theIter)) =E= 0 ;

newcut3(theIter)..
theta
+ sum((n,p,e,m),  alphaPar(m,n,p,e,theIter)* maxProdV(m,n,p,e) )
+ sum((n,ni,e,m), lambdaF(m,n,ni,e,theIter)*(maxFlowV(n,ni,e) + sum(mi $ ancestor(m,mi), fe(mi,n,ni,e))))
+ sum((n,e,m),    lambdaX(m,n,e,theIter)*(maxTransV(n,e)      + sum(mi $ ancestor(m,mi), xe(mi,n,e))))
- sum((n,e,p,d,theIterI,m), pb(m)*dc(m)*(int(m,n,d,e) - slp(m,n,d,e)*( qsPar(m,n,p,d,e,theIter)+ sum(pi,qsPar(m,n,pi,d,e,theIter)) ))*iota(theIterI)*qsPar(m,n,p,d,e,theIterI)  )
+ sum((n,p,e,theIterI,m),   pb(m)*dc(m)*(2*k1*qpPar(m,n,p,e,theIter) + k2)*qpPar(m,n,p,e,theIterI)*iota(theIterI) )
+ sum((n,ni,e,theIterI,m),  pb(m)*dc(m)*(k3*2*fPar(m,n,ni,e,theIter) - upsilonPar(m,n,n,e,theIter))*fPar(m,n,ni,e,theIterI)*iota(theIterI)  )
+ sum((n,e,ei,theIterI,m),  pb(m)*dc(m)*(k5*2*xPar(m,n,e,ei,theIter) - phiPar(m,n,e,e,theIterI)   )*xPar(m,n,e,ei,theIterI)*iota(theIterI)  )
 =G= 0 ;

*Subproblem
model SPC2 /c12.qs2, c22.qp2, c32.qt2, c42.qc2, c52.f2, c72.x2, maxProd2.alpha2, massBalance2.beta2, maxFlow2.gamma2, maxTrans2.epsilon2, marketClearingF2.upsilon2, marketClearingX2.phi2 /;

*Master Problem
model MPVI /maxExF.delta, maxExX.zeta, c6m.fe, c8m.xe, newcut3.iota, convexityCon.theta  / ;

scalar converged number will be set equal to 1 if convergence criterion is satisfied /0/;

*Time measurements
starttime = jnow*24*3600;
MPel = 0.0000000;
SPel = 0.0000000;

*INITIALIZATION OF VARIABLES
theIter(iteration) = no ;
iotaPar(iteration) =0;
iota.l(iteration) = 0;
xeItPar(m,n,e) = 0 ;
feItPar(m,n,ni,e) = 0 ;
CG2 = 1000000;
rm = 0;
rs = 0;

*------------START ITERATING--------------
loop(iteration$(not converged) ,
*------------SOVE SPs--------------
         loop(m,
                 pbp = pb(m);
                 dcp = dc(m);
                 intp(n,d,e) = int(m,n,d,e);
                 slpp(n,d,e) = slp(m,n,d,e);
                 lp(n,ei,e) = l(m,n,ei,e);
                 fmip(n,ni,e) = sum( mi $ ancestor(m,mi), feItPar(mi,n,ni,e));
                 xmip(n,e) =  sum(mi $ ancestor(m,mi), xeItPar(mi,n,e));
                 maxProdVp(n,p,e) = maxProdV(m,n,p,e);

                 SPstart = jnow*24*3600;

                 solve SPC2 using MCP ;

                 rs = rs + SPC2.resUsd  ;
                 SPel=SPel + (jnow*24*3600-SPstart);

                 lambdaF(m,n,ni,e,iteration) = maxFlow2.m(n,ni,e);
                 lambdaX(m,n,e,iteration) = maxTrans2.m(n,e);
                 qsPar(m,n,p,d,e,iteration) = qs2.l(n,p,d,e);
                 qpPar(m,n,p,e,iteration) = qp2.l(n,p,e);
                 fPar(m,n,ni,e,iteration) = f2.l(n,ni,e);
                 xPar(m,n,e,ei,iteration) = x2.l(n,e,ei);
                 alphaPar(m,n,p,e,iteration) = alpha2.l(n,p,e);
                 upsilonPar(m,n,n,e,iteration) = upsilon2.l(n,n,e);
                 phiPar(m,n,e,e,iteration) = phi2.l(n,e,e);
         );
*------------COMPUTE CG--------------
         if(ord(iteration)>1,
*                  hs= jnow*24*3600;
                  CG(iteration) = sum(m, pb(m)*dc(m)*(
                  -sum((n,e,p,d), ((int(m,n,d,e) - slp(m,n,d,e)*( qsPar(m,n,p,d,e,iteration))- slp(m,n,d,e)*sum(pi,qsPar(m,n,pi,d,e,iteration)) ) - (int(m,n,d,e) - slp(m,n,d,e)*( sum( theIterI, qsPar(m,n,p,d,e,theIterI)*iotaPar(theIterI))) - sum(( theIterI,pi),qsPar(m,n,pi,d,e,theIterI)*iotaPar(theIterI))))*(sum(theIterI,qsPar(m,n,p,d,e,theIterI)*iotaPar(theIterI))))
                  +sum((n,e,p)  , ((2*k1*qpPar(m,n,p,e,iteration) + k2) - (2*k1*(sum(theIterI,(qpPar(m,n,p,e,theIterI)*iotaPar(theIterI)))) +k2))*(sum(theIterI,(qpPar(m,n,p,e,theIterI)*iotaPar(theIterI)))))
                  +sum((n,ni,e) , ((k3*2*fPar(m,n,ni,e,iteration) - upsilonPar(m,n,ni,e,iteration)) - ((sum(theIterI,(k3*2*fPar(m,n,ni,e,theIterI) - upsilonPar(m,n,ni,e,theIterI))*iotaPar(theIterI)))))*(sum(theIterI,fPar(m,n,ni,e,theIterI)*iotaPar(theIterI))))
                  +sum((n,e,ei) , ((k5*2*xPar(m,n,e,ei,iteration) - phiPar(m,n,e,ei,iteration)) -     ((sum(theIterI,(k5*2*xPar(m,n,e,ei,theIterI) - phiPar(m,n,e,ei,theIterI))*iotaPar(theIterI)))))*(sum(theIterI,xPar(m,n,e,ei,theIterI)*iotaPar(theIterI))))
                  +sum((n,ni,e), (lambdaF(m,n,ni,e,iteration) - sum(theIterI,(lambdaF(m,n,ni,e,theIterI)*iota.l(theIterI))))*(maxFlowV(n,ni,e) + sum(mi $ ancestor(m,mi),feItPar(mi,n,ni,e)))   )
                  +sum((n,e,ei), (lambdaX(m,n,e,iteration) - sum(theIterI,(lambdaX(m,n,e,theIterI)*iotaPar(theIterI))))*(maxTransV(n,e) + sum(mi $ ancestor(m,mi),xeItPar(mi,n,e)))  )
                  +sum((n,p,e),  (alphaPar(m,n,p,e,iteration) - sum(theIterI,(alphaPar(m,n,p,e,theIterI)*iotaPar(theIterI))))*maxProdV(m,n,p,e) )
                  )) ;
*------------CONVERGENCE CHECK--------------
                  if (((CG(iteration))<epsi and (CG(iteration))> -1* epsi),
                          converged = 1;
                          elapsed = (jnow*24*3600 - starttime);
                  );

                  CG2 = CG(iteration);
                 display CG;
         );

         theIter(iteration) = yes  ;

*------------SOLVE MP--------------
       MPstart=jnow*24*3600;
       solve MPVI using mcp;
       if (converged < 1,
*        running time computations. the last MP is not included in this sum, as convergence is already detected
         MPel=MPel +(jnow*24*3600-MPstart);
       );
       rm = rm + MPVI.resUsd;

*------------PREPARE INFO TO BE PASSED TO SP + ALTERNETIVE CONVERGENCE CHECK--------------
       loop(theIter,
           iotaPar(theIter)=iota.l(theIter);
       );

       if( converged <> 1 and ord(iteration)>4 and smax((m,n,e), (xe.l(m,n,e) - xeItPar(m,n,e))) < epsi2 and smin((m,n,e), (xe.l(m,n,e) - xeItPar(m,n,e)) > -1*epsi2) and (smax((m,n,ni,e), (fe.l(m,n,ni,e) - feItPar(m,n,ni,e))) < epsi2 and smin((m,n,ni,e), (fe.l(m,n,ni,e) - feItPar(m,n,ni,e))) > -1*epsi2) ,
           converged = 1;
           elapsed = (jnow*24*3600 - starttime);
       );
       xeItPar(m,n,e) = xe.l(m,n,e)  ;
       feItPar(m,n,ni,e) = fe.l(m,n,ni,e)  ;
);
*------------DISPLAY CONVERGENCE RESULTS--------------
display CG, elapsed, theIter, MPel, SPel, rm, rs;


