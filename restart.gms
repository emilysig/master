
*$Include C:\Users\Kevin\Desktop\Dropbox\master\gams\dataSmallDetMulti3d.gms
*$Include C:\Users\Kevin\Desktop\Dropbox\master\gams\m1.gms
*$Include C:\Users\Emily\Dropbox\master\gams\dataSmallDetMulti3d.gms
$Include C:\Users\Emily\Dropbox\master\gams\m1.gms
*$Include C:\Users\Emily\Dropbox\master\gams\prosjekt.gms

*Full NLP
*SP MCP
*MP LP    -> MCP
*SP not profit maximization






*TO DO  Make an MPC0 with complete independence from iteration/theIter!!!







scalar UB upper bound /INF/ ;
scalar LB lower bound /-INF/ ;
scalar epsi convergence tolerance /1e-4/ ;
scalar thetadown initial value or guess for alpha /-4000/ ;


set iteration maximum number of iterations /1*20/ ;
set theIter(iteration) dynamic set of conducted iterations ;
set theIterA(iteration) dynamic set of conducted iterations -1 ;

alias (theIter, theIterI);

parameters
xeItPar(m,n,e)                   parameter of transformation expansion values at present iteration
feItPar(m,n,n,e)                 parameter of flow expansion values at present iteration
xePar(m,n,e,iteration)           parameter of transformation expansion values at each iteration
fePar(m,n,n,e,iteration)         parameter of flow expansion values at each iteration
lambdaF(m,n,n,e,iteration)       dual to constraint fixing flow expansions
lambdaX(m,n,e,iteration)         dual to constraint fixing transformation expansions
lambdaParF(m,n,n,e,iteration)    dual to constraint fixing flow expansions
lambdaParX(m,n,e,iteration)      dual to constraint fixing transformation expansions
sj(iteration)                    solution of SP at each iteration
sjPar
qsPar(m,n,p,d,e,iteration)       parameter of qs in MP (vi - based)
qpPar(m,n,p,e,iteration)         parameter of qp in MP (vi - based)
fPar(m,n,n,e,iteration)          parameter of f in MP (vi - based)
xPar(m,n,e,e,iteration)          parameter of x in MP (vi - based)
alphaPar(m,n,p,e,iteration)      parameter of alpha in MP (VI - based)
CG(iteration)                    convergence gap
;


parameter oI(iteration) iteration count test;

parameter rep_arci(*,*,*,*,*),rep_trfi(*,*,*,*),rep_iter(*,*)
         rep_arcd(*,*,*,*,*),rep_trfd(*,*,*,*);

positive variables
qp(m,n,p,e)              quantity produced
qs(m,n,p,d,e)            quantity supplied
qt(m,n,n,p,e)            quantity transported
qc(m,n,p,e,e)            quantity transformed
f(m,n,n,e)               flow (sum of transported quantities)
fe(m,n,n,e)              flow capacity expansions
x(m,n,e,e)               transformation
xe(m,n,e)                transformation expansions

alpha(m,n,p,e)           dual to maxProd
gamma(m,n,n,e)           dual to maxFlow
delta(m,n,n,e)           dual to maxExF
epsilon(m,n,e)           dual to maxTrans
zeta(m,n,e)              dual to maxExX
eta                      dual to thetaCon
iota(iteration)          dual to newcut
;

variables
theta                    approximation of SP solution in MP objective
zMP                      optimal value of MP objective
zMP0                     optimal value of initial MP objective
zSP                      optimal value of SP objective
z                        optimal value of total
beta(m,n,p,e)            dual to massBalance
*zeta(m,n,p,e)            dual to eq massBalance in SP
upsilon(m,n,n,e)         dual to eq marketClearingF (flow) and unit price of transportation in SP
phi(m,n,e,e)             dual to eq marketClearingX
;


Equations

MPobj                    MP objective function
SPobj                    SP objective function
obj                      objective for total

thetacon                 lower bound on alpha in MP
convexityCon             convexity constraint in MP

newcut                   cuts added in MP
newcut2                  cuts added in MP
newcut3                  cuts added in MP VI-based

maxProd(m,n,p,e)         production constraint for producer
massBalance(m,n,p,e)     mass balance constraint for producer

maxFlowOrg(m,n,n,e)      flow constraint for transporter in Original
maxFlow(m,n,n,e)         flow constraint for transporter in SP
maxFlowS(m,n,n,e)        constraint in SP binding flow expansions to value found in MP
maxExF(m,n,n,e)          expansion constraint (upper limit) for transporter in MP

maxTransOrg(m,n,ei)      transformation constraint for transformer in Original
maxTrans(m,n,ei)         transformation constraint for transformer in SP
maxTransS(m,n,e)         constraint in SP binding transformation expansions to value found in MP
maxExX(m,n,e)            expansion constraint for transformer in MP

c1(m,n,p,d,e)
c2(m,n,p,e)
c3(m,n,ni,p,e)
c4(m,n,p,ei,e)
c5(m,n,ni,e)
c6(m,n,ni,e)
c6m(m,n,ni,e)
c7(m,n,ei,e)
c8(m,n,e)
c8m(m,n,e)

marketClearingF(m,n,n,e) links the producer and transporter problem by assuring concistiency among transported quantities in SP
marketClearingX(m,n,e,e) links the producer and transformer problem by assuring concistiency among transformed quantities in SP

;

*c1(m,n,p,d,e)$ spno(n,p)..      - pb(m)*dc(m)*(int(m,n,d,e) - slp(m,n,d,e)*qs(m,n,p,d,e) - slp(m,n,d,e)*sum(pi,qs(m,n,pi,d,e))) - beta(m,n,p,e) =G= 0;
c1(m,n,p,d,e)..                  - pb(m)*dc(m)*(int(m,n,d,e) - slp(m,n,d,e)*qs(m,n,p,d,e)) - beta(m,n,p,e) =G= 0;

c2(m,n,p,e)..                    pb(m)*dc(m)*(2*k1*qp(m,n,p,e) + k2) + alpha(m,n,p,e) + beta(m,n,p,e) =G= 0;
c3(m,n,ni,p,e)..                 pb(m)*dc(m)*(upsilon(m,n,ni,e)) - beta(m,n,p,e) + beta(m,ni,p,e) =G= 0;
c4(m,n,p,ei,e)..                 pb(m)*dc(m)*(phi(m,n,ei,e)) + l(m,n,ei,e)*beta(m,n,p,e) - beta(m,n,p,ei) =G= 0;
c5(m,n,ni,e)..                   pb(m)*dc(m)*(2*k3*f(m,n,ni,e) - upsilon(m,n,ni,e)) + gamma(m,n,ni,e) =G= 0;
c6m(m,n,ni,e)..                  pb(m)*dc(m)*k4 + delta(m,n,ni,e) - sum(theIter, iota(theIter)*sum(mi $ ancestor (mi,m), pb(m)*dc(m)*lambdaF(mi,n,ni,e,theIter) ) ) =G= 0;
c6(m,n,ni,e)..                   pb(m)*dc(m)*k4 + delta(m,n,ni,e) - sum(mi $ ancestor (mi,m), gamma(mi,n,ni,e) ) =G= 0;
c7(m,n,ei,e)..                   pb(m)*dc(m)*(2*k5*x(m,n,ei,e) - phi(m,n,ei,e)) + l(m,n,ei,e)*epsilon(m,n,e) =G= 0;
c8m(m,n,e)..                     pb(m)*dc(m)*k6 + zeta(m,n,e) - sum(theIter, iota(theIter)*sum(mi $ ancestor (mi,m), pb(m)*dc(m)*lambdaX(mi,n,e,theIter) ) ) =G= 0;
*c8(m,n,e)..                     pb(m)*dc(m)*k6 + zeta(m,n,e) - sum(mi $ ancestor (m,mi), epsilon(mi,n,e))  =G= 0;
c8(m,n,e)..                      pb(m)*dc(m)*k6 + zeta(m,n,e) - sum(mi $ ancestor (mi,m), epsilon(mi,n,e))  =G= 0;

maxProd(m,n,p,e)..               maxProdV(m,n,p,e) - qp(m,n,p,e) =G= 0;
massBalance(m,n,p,e)..           - sum(d,qs(m,n,p,d,e)) - sum( ni $ (not (ord(ni) = ord(n))), qt(m,n,ni,p,e)) - sum( ei $ (not (ord(ei) = ord(e))), qc(m,n,p,e,ei)) + qp(m,n,p,e) + sum( ni $ (not (ord(ni) = ord(n))) , qt(m,ni,n,p,e)) + sum( ei $ (not (ord(ei) = ord(e))), l(m,n,ei,e)*qc(m,n,p,ei,e)) =E=  0 ;

maxFlowOrg(m,n,ni,e)..           maxFlowV(n,ni,e) + sum( mi $ ancestor(m,mi), fe(mi,n,ni,e)) - f(m,n,ni,e) =G=  0 ;
maxFlow(m,n,ni,e)..              maxFlowV(n,ni,e) + sum( mi $ ancestor(m,mi), feItPar(mi,n,ni,e)) - f(m,n,ni,e) =G=  0 ;
maxExF(m,n,ni,e)..               maxExFV(n,ni,e) - fe(m,n,ni,e) =G= 0 ;

maxTransOrg(m,n,e)..             maxTransV(n,e) + sum(mi $ ancestor(m,mi), xe(mi,n,e)) - sum(ei  $ (not (ord(e) = ord(ei))) , l(m,n,ei,e)*x(m,n,ei,e)) =G=0 ;
maxTrans(m,n,e)..                maxTransV(n,e) + sum(mi $ ancestor(m,mi), xeItPar(mi,n,e)) - sum(ei  $ (not (ord(e) = ord(ei))) , l(m,n,ei,e)*x(m,n,ei,e)) =G=0 ;
maxExX(m,n,e)..                  maxExXV(n,e) - xe(m,n,e) =G= 0 ;

marketClearingF(m,n,ni,e)..      - f(m,n,ni,e) + sum( p, qt(m,n,ni,p,e)) =E= 0 ;
marketClearingX(m,n,e,ei)..      - x(m,n,e,ei) + sum( p, qc(m,n,p,e,ei)) =E= 0 ;


*SPobj.. zSP =E= - (sum((m,e,n), pb(m)*dc(m)*(
*sum(p $ spno(n,p), sum(d, (int(m,n,d,e) - slp(m,n,d,e)*(sum( pi, qs(m,n,pi,d,e))))*qs(m,n,p,d,e))
* - k1*sqr(qp(m,n,p,e)) - k2*qp(m,n,p,e))
*+ sum(d,0.5*slp(m,n,d,e)*sqr(sum(p, qs(m,n,p,d,e))))
*- sum(ni,k5*sqr(sum(p,qt(m,n,ni,p,e))))
*- sum(ei,k6*sqr(sum(p,qc(m,n,p,ei,e))))
*)) );
SPobj.. zSP =E= - (sum((m,e,n), pb(m)*dc(m)*(
sum(p $ spno(n,p), sum(d, (int(m,n,d,e) - slp(m,n,d,e)*(sum( pi, qs(m,n,pi,d,e))))*qs(m,n,p,d,e))
 - sqr(qp(m,n,p,e)) - qp(m,n,p,e))
+ sum(d,0.5*slp(m,n,d,e)*sqr(sum(p, qs(m,n,p,d,e))))
- sum(ni,k3*sqr(f(m,n,ni,e)))
- sum(ei,k5*sqr(x(m,n,ei,e)))
)) ) ;

MPobj.. zMP =E= sum((m,n), pb(m)*dc(m)* (
  sum((e,ni), k4*fe(m,n,ni,e))
+ sum(ei, k6*xe(m,n,ei))
)) + theta ;

thetacon.. - thetadown + theta =G= 0 ;
convexityCon.. 1 - sum(theIter, iota(theIter)) - eta =E= 0 ;

newcut(theIter).. theta - sj(theIter) +
sum( (m,n,ni,ei), lambdaParF(m,n,ni,ei,theIter)*sum(mi $ ancestor(m,mi ), (fe(mi,n,ni,ei) - fePar(mi,n,ni,ei,theIter)))) +
sum( (m,n,ei), lambdaParX(m,n,ei,theIter)*(sum(mi $ ancestor(m,mi),(xe(mi,n,ei) -  xePar(mi,n,ei,theIter))))) =G= 0 ;

newcut2(theIter).. theta - sj(theIter) -
sum( (m,n,ni,e), pb(m)*dc(m)*lambdaParF(m,n,ni,e,theIter)*sum(mi $ ancestor(m,mi ), (fe(mi,n,ni,e) - fePar(mi,n,ni,e,theIter)))) -
sum( (m,n,ei), pb(m)*dc(m)*lambdaParX(m,n,ei,theIter)*(sum(mi $ ancestor(m,mi),(xe(mi,n,ei) -  xePar(mi,n,ei,theIter))))) =G= 0 ;

newcut3(theIter)..
theta
+ sum((n,p,e,m),  alphaPar(m,n,p,e,theIter)* maxProdV(m,n,p,e) )
+ sum((n,ni,e,m), lambdaParF(m,n,ni,e,theIter)*(maxFlowV(n,ni,e) + sum(mi $ ancestor(m,mi), fe(mi,n,ni,e))))
+ sum((n,e,m),    lambdaParX(m,n,e,theIter)*(maxTransV(n,e)      + sum(mi $ ancestor(m,mi), xe(mi,n,e))))
- sum((n,e,p,d,theIterI,m), pb(m)*dc(m)*(int(m,n,d,e) - slp(m,n,d,e)*( qsPar(m,n,p,d,e,theIter)))*iota(theIterI)*qsPar(m,n,p,d,e,theIterI)  )
+ sum((n,p,e,theIterI,m),   pb(m)*dc(m)*(2*k1*qpPar(m,n,p,e,theIter) + k2)*qpPar(m,n,p,e,theIterI)*iota(theIterI) )
+ sum((n,ni,e,theIterI,m),  pb(m)*dc(m)*k3*2*fPar(m,n,ni,e,theIter)      *fPar(m,n,ni,e,theIterI)*iota(theIterI)  )
+ sum((n,e,ei,theIterI,m),  pb(m)*dc(m)*k5*2*xPar(m,n,e,ei,theIter)      *xPar(m,n,e,ei,theIterI)*iota(theIterI)  )
 =G= 0 ;

*newcut3(theIter)..
*  theta
*+ sum((n,p,e,m), alphaPar(m,n,p,e,theIter) * maxProdV(m,n,p,e))
*+ sum((n,ni,e,m), lambdaF(m,n,ni,e,theIter)*(maxFlowV(n,ni,e) + sum(mi $ ancestor(m,mi), fe(mi,n,ni,e))))
*+ sum((n,e,m),    lambdaX(m,n,e,theIter)   *(maxTransV(n,e)   + sum(mi $ ancestor(m,mi), xe(mi,n,e))))
*- sum((n,e,p,d,theIterI,m) $ (oI(theIterI) <= oI(theIter)),  pb(m)*dc(m)*(int(m,n,d,e) - slp(m,n,d,e)*(qsPar(m,n,p,d,e,theIter) + sum(pi,qsPar(m,n,pi,d,e,theIter)) ))*iota(theIterI)*qsPar(m,n,p,d,e,theIterI)   )
*+ sum((n,e,p,theIterI,m)   $ (oI(theIterI) <= oI(theIter)),  pb(m)*dc(m)*(k1*2*qpPar(m,n,p,e,theIter) + k2)*qpPar(m,n,p,e,theIterI)*iota(theIterI) )
*+ sum((n,ni,e,theIterI,m)  $ (oI(theIterI) <= oI(theIter)),  pb(m)*dc(m)*(k3*2*fPar(m,n,ni,e,theIter))*fPar(m,n,ni,e,theIterI)*iota(theIterI)  )
*+ sum((n,e,ei,theIterI,m)  $ (oI(theIterI) <= oI(theIter)),  pb(m)*dc(m)*(k5*2*xPar(m,n,ei,e,theIter))*xPar(m,n,ei,e,theIterI)*iota(theIterI)  )
*=G= 0 ;

obj.. z =E=  - (sum((m,e,n), pb(m)*dc(m)*(
sum(p $ spno(n,p), sum(d,(int(m,n,d,e) - slp(m,n,d,e)*(sum( pi, qs(m,n,pi,d,e))))*qs(m,n,p,d,e))
 - k1*sqr(qp(m,n,p,e)) - k2*qp(m,n,p,e))
+ sum(d,0.5*slp(m,n,d,e)*sqr(sum(p, qs(m,n,p,d,e))))
- sum(ni,k3*sqr(f(m,n,ni,e)) + k4*fe(m,n,ni,e))
- sum(ei,k5*sqr(x(m,n,ei,e)))
- k6*xe(m,n,e)
)) )
;


model org /obj, maxProd, massBalance, maxFlowOrg, maxTransOrg, maxExX, maxExF, marketClearingF, marketClearingX /;
model SP definition of SP /SPobj, maxProd, massBalance, maxFlow, maxTrans, marketClearingF, marketClearingX /;
model MP definition of MP /MPobj, maxExX, maxExF, thetacon, newcut/;
model MP0 definition of initial MP /MPobj, maxExX, maxExF, thetacon/;

*orgC is wrong....
model orgC /c1.qs, c2.qp, c3.qt, c4.qc, c5.f, c6.fe, c7.x, c8.xe, maxProd.alpha , massBalance.beta , maxFlowOrg.gamma , maxTransOrg.epsilon, maxExF.delta, maxExX.zeta, marketClearingF.upsilon, marketClearingX.phi  /;
model SPC /c1.qs, c2.qp, c3.qt, c4.qc, c5.f, c7.x, maxProd.alpha, massBalance.beta, maxFlow.gamma, maxTrans.epsilon, marketClearingF.upsilon, marketClearingX.phi /;
*model MPC /c6.fe, c8.xe,  maxExF.delta, maxExX.zeta /  ;
model MPC /c6m.fe, c8m.xe, newcut.iota, maxExF.delta, maxExX.zeta, thetacon.eta, convexityCon.theta / ;
*model MPC0 /c6.fe, c8.xe, , maxExF.zeta, maxExX.delta, thetacon.eta, convexityCon.theta /

model MPVI /maxExF.delta, maxExX.zeta, c6m.fe, c8m.xe, newcut3.iota, thetacon.eta, convexityCon.theta  / ;

option NLP = minos;
*solve org using nlp minimizing z;
solve orgC using MCP;

rep_arci(m,n,ni,e,'full_mcp') =fe.l(m,n,ni,e);
rep_arcd(m,n,ni,e,'full_mcp') =gamma.l(m,n,ni,e);
rep_trfi(m,n,e,'full_mcp')    = xe.l(m,n,e);
rep_trfd(m,n,e,'full_mcp')=     epsilon.l(m,n,e);
rep_iter(iteration,'OPT')    = 7777;

scalar converged number will be set equal to 1 if convergence criterion is satisfied /0/;
scalar iterationCount to keep track of iterations /0/;

*------------SOLVE INTIAL MP----------
Solve MP0 using nlp minimizing zMP ;

LB = zMP.l ;
theIter(iteration) = no ;
theIterA(iteration) = no ;

*------------START ITERATING--------------
loop(iteration$(not converged) ,

         iterationCount = iterationCount +  1;

         xePar(m,n,e,iteration) = xe.l(m,n,e)  ;
         fePar(m,n,ni,e,iteration) = fe.l(m,n,ni,e)  ;

         xeItPar(m,n,e) = xe.l(m,n,e)  ;
         feItPar(m,n,ni,e) = fe.l(m,n,ni,e)  ;

         solve SPC using MCP ;
sj(iteration) =  - (sum((m,e,n), pb(m)*dc(m)*(
sum(p $ spno(n,p), sum(d, (int(m,n,d,e) - slp(m,n,d,e)*(sum( pi, qs.l(m,n,pi,d,e))))*qs.l(m,n,p,d,e))
 - sqr(qp.l(m,n,p,e)) - qp.l(m,n,p,e))
+ sum(d,0.5*slp(m,n,d,e)*sqr(sum(p, qs.l(m,n,p,d,e))))
- sum(ni,k3*sqr(f.l(m,n,ni,e)))
- sum(ei,k5*sqr(x.l(m,n,ei,e)))
)) );

sjPar =  - (sum((m,e,n), pb(m)*dc(m)*(
sum(p $ spno(n,p), sum(d, (int(m,n,d,e) - slp(m,n,d,e)*(sum( pi, qs.l(m,n,pi,d,e))))*qs.l(m,n,p,d,e))
 - sqr(qp.l(m,n,p,e)) - qp.l(m,n,p,e))
+ sum(d,0.5*slp(m,n,d,e)*sqr(sum(p, qs.l(m,n,p,d,e))))
- sum(ni,k3*sqr(f.l(m,n,ni,e)))
- sum(ei,k5*sqr(x.l(m,n,ei,e)))
)) );
*         solve SP minimizing zSP using NLP;

         lambdaParF(m,n,ni,e,iteration) = maxFlow.m(m,n,ni,e);
         lambdaParX(m,n,e,iteration) = maxTrans.m(m,n,e);
         lambdaF(m,n,ni,e,iteration) = maxFlow.m(m,n,ni,e);
         lambdaX(m,n,e,iteration) = maxTrans.m(m,n,e);
         qsPar(m,n,p,d,e,iteration) = qs.l(m,n,p,d,e);
         qpPar(m,n,p,e,iteration) = qp.l(m,n,p,e);
         fPar(m,n,ni,e,iteration) = f.l(m,n,ni,e);
         xPar(m,n,e,ei,iteration) = x.l(m,n,e,ei);
         alphaPar(m,n,p,e,iteration) = alpha.l(m,n,p,e);

*CG(theIter) = sum(m, pb(m)*dc(m)*(
*-sum((n,e,p,d), ((int(m,n,d,e) - slp(m,n,d,e)*( qsPar(m,n,p,d,e,theIter))) - (int(m,n,d,e) - slp(m,n,d,e)*( sum( theIterI, qsPar(m,n,p,d,e,theIterI)*iota(theIterI)))))*(sum(theIterI,qsPar(m,n,p,d,e,theIterI)*iota(theIterI))))
*-sum((n,e,p)  , ((2*k1*qpPar(m,n,p,e,theIter) + k2) - (2*k1*(sum(theIterI(qpPar(m,n,p,e,theIterI)*iota(theIterI)))) +k2))*(sum(theIterI(qpPar(m,n,p,e,theIterI)*iota(theIterI)))))
*-sum((n,ni,e) , ((k3*2*fPar(m,n,ni,e,theIter)) - (k3*2*(sum(theIterI,fPar(m,n,ni,e,theIterI)*iota(theIterI)))))*(sum(theIterI,fPar(m,n,ni,e,theIterI)*iota(theIterI))))
*-sum((n,e,ei) , ((k5*2*xPar(m,n,e,ei,theIter)) - (k5*2*(sum(theIterI,xPar(m,n,e,ei,theIterI)*iota(theIterI)))))*(sum(theIterI,xPar(m,n,e,ei,theIterI)*iota(theIterI))))
*+sum((n,ni,e), (lambdaF(m,n,ni,e,theIter) - sum(theIterI,(lambdaF(m,n,ni,e,theIterI)*iota(theIterI))))*(maxFlowV(n,ni,e) + sum(mi $ ancestor(m,mi),fe(mi,n,ni,e)))   )
*+sum((n,e,ei), (lambdaX(m,n,e,ei,theIter) - sum(theIterI,(lambdaX(m,n,e,ei,theIterI)*iota(theIterI))))*(maxTransV(n,e,ei) + sum(mi $ ancestor(m,mi),xe(mi,n,e,ei)))  )
*+sum((n,p,e),  (alphaPar(m,n,p,e,theIter) - sum(theIterI,(alphaPar(m,n,p,e,theIterI)*iota(theIterI))))*maxProdV(m,n,p,e) )
*))
*;

          theIter(iteration) = yes ;
*          display CG(iteration);
*         display lambdaParX, lambdaParF, sjPar;

*         UB = min (sj(iteration) + zMP.l - theta.l ,UB);
*         display UB, LB ;
*         if ( abs(LB - UB) < epsi ,
*                 converged = 1 ;
*         );
*        solve MP minimizing zMP using nlp ;
*         solve MPC using mcp ;
        solve MPVI using mcp;

if(iterationCount>2,
CG(iteration) = sum(m, pb(m)*dc(m)*(
-sum((n,e,p,d), ((int(m,n,d,e) - slp(m,n,d,e)*( qsPar(m,n,p,d,e,iteration))) - (int(m,n,d,e) - slp(m,n,d,e)*( sum( theIterI, qsPar(m,n,p,d,e,theIterI)*iota.l(theIterI)))))*(sum(theIterI,qsPar(m,n,p,d,e,theIterI)*iota.l(theIterI))))
-sum((n,e,p)  , ((2*k1*qpPar(m,n,p,e,iteration) + k2) - (2*k1*(sum(theIterI,(qpPar(m,n,p,e,theIterI)*iota.l(theIterI)))) +k2))*(sum(theIterI,(qpPar(m,n,p,e,theIterI)*iota.l(theIterI)))))
-sum((n,ni,e) , ((k3*2*fPar(m,n,ni,e,iteration)) - (k3*2*(sum(theIterI,fPar(m,n,ni,e,theIterI)*iota.l(theIterI)))))*(sum(theIterI,fPar(m,n,ni,e,theIterI)*iota.l(theIterI))))
-sum((n,e,ei) , ((k5*2*xPar(m,n,e,ei,iteration)) - (k5*2*(sum(theIterI,xPar(m,n,e,ei,theIterI)*iota.l(theIterI)))))*(sum(theIterI,xPar(m,n,e,ei,theIterI)*iota.l(theIterI))))
+sum((n,ni,e), (lambdaF(m,n,ni,e,iteration) - sum(theIterI,(lambdaF(m,n,ni,e,theIterI)*iota.l(theIterI))))*(maxFlowV(n,ni,e) + sum(mi $ ancestor(m,mi),fe.l(mi,n,ni,e)))   )
+sum((n,e,ei), (lambdaX(m,n,e,iteration) - sum(theIterI,(lambdaX(m,n,e,theIterI)*iota.l(theIterI))))*(maxTransV(n,e) + sum(mi $ ancestor(m,mi),xe.l(mi,n,e)))  )
+sum((n,p,e),  (alphaPar(m,n,p,e,iteration) - sum(theIterI,(alphaPar(m,n,p,e,theIterI)*iota.l(theIterI))))*maxProdV(m,n,p,e) )
))
;

display CG;
);


*         zMP.l = sum((m,n), pb(m)*dc(m)* (sum((e,ni), k3*fe.l(m,n,ni,e))
*+ sum(ei, k4*xe.l(m,n,ei)))  ) + theta.l ;
*         solve MP minimizing zMP using nlp ;

         theIterA(iteration) = yes ;


rep_arci(m,n,ni,e,iteration) = fe.l(m,n,ni,e);
rep_arcd(m,n,ni,e,iteration) = lambdaF(m,n,ni,e,iteration);
rep_trfi(m,n,e,iteration)    = xe.l(m,n,e);
rep_trfd(m,n,e,iteration)    = lambdaX(m,n,e,iteration);
rep_iter(iteration,'LB')     = LB;
rep_iter(iteration,'UB')     = UB;

);

option rep_arci:4:4:1, rep_trfi:4:3:1, rep_iter:4:1:1
       rep_arcd:4:4:1, rep_trfd:4:3:1;

display rep_arci, rep_trfi, rep_arcd, rep_trfd, rep_iter;

