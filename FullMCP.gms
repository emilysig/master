*The input file is included first.
$Include C:\Users\Emily\Dropbox\master\gitfiles\master\eurdata.gms

*TIME scalars
scalar
starttime
elapsed
ORGstart
ORGel
r
;
starttime = jnow*24*3600;


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
;

variables
beta(m,n,p,e)            dual to massBalance
upsilon(m,n,n,e)         dual to eq marketClearingF (flow) and unit price of transportation in SP
phi(m,n,e,e)             dual to eq marketClearingX
;


Equations

maxProd(m,n,p,e)         production constraint for producer
massBalance(m,n,p,e)     mass balance constraint for producer

maxFlowOrg(m,n,n,e)      flow constraint for transporter in Original
maxExF(m,n,n,e)          expansion constraint (upper limit) for transporter in MP

maxTransOrg(m,n,ei)      transformation constraint for transformer in Original
maxExX(m,n,e)            expansion constraint for transformer in MP

c1(m,n,p,d,e)            comp to qs
c2(m,n,p,e)              comp to qp
c3(m,n,ni,p,e)           comp to qt
c4(m,n,p,ei,e)           comp to qc
c5(m,n,ni,e)             comp to f
c6(m,n,ni,e)             comp to fe
c7(m,n,ei,e)             comp to x
c8(m,n,e)                comp to xe

marketClearingF(m,n,n,e) links the producer and transporter problem by assuring concistiency among transported quantities in SP
marketClearingX(m,n,e,e) links the producer and transformer problem by assuring concistiency among transformed quantities in SP

;

c1(m,n,p,d,e)$ spno(n,p)..       - pb(m)*dc(m)*(int(m,n,d,e) - slp(m,n,d,e)*qs(m,n,p,d,e) - slp(m,n,d,e)*sum(pi,qs(m,n,pi,d,e))) - beta(m,n,p,e) =G= 0;
*The second c1 can be used for monopoly
*c1(m,n,p,d,e)..                  - pb(m)*dc(m)*(int(m,n,d,e) - slp(m,n,d,e)*qs(m,n,p,d,e)) - beta(m,n,p,e) =G= 0;

c2(m,n,p,e)..                    pb(m)*dc(m)*(2*k1*qp(m,n,p,e) + k2) + alpha(m,n,p,e) + beta(m,n,p,e) =G= 0;
c3(m,n,ni,p,e)..                 pb(m)*dc(m)*(upsilon(m,n,ni,e)) - beta(m,n,p,e) + beta(m,ni,p,e) =G= 0;
c4(m,n,p,ei,e)..                 pb(m)*dc(m)*(phi(m,n,ei,e)) + l(m,n,ei,e)*beta(m,n,p,e) - beta(m,n,p,ei) =G= 0;
c5(m,n,ni,e)..                   pb(m)*dc(m)*(2*k3*f(m,n,ni,e) - upsilon(m,n,ni,e)) + gamma(m,n,ni,e) =G= 0;
c6(m,n,ni,e)..                   pb(m)*dc(m)*k4 + delta(m,n,ni,e) - sum(mi $ ancestor (mi,m), gamma(mi,n,ni,e) ) =G= 0;
c7(m,n,ei,e)..                   pb(m)*dc(m)*(2*k5*x(m,n,ei,e) - phi(m,n,ei,e)) + l(m,n,ei,e)*epsilon(m,n,e) =G= 0;
c8(m,n,e)..                      pb(m)*dc(m)*k6 + zeta(m,n,e) - sum(mi $ ancestor (mi,m), epsilon(mi,n,e))  =G= 0;

maxProd(m,n,p,e)..               maxProdV(m,n,p,e) - qp(m,n,p,e) =G= 0;
massBalance(m,n,p,e)..           - sum(d,qs(m,n,p,d,e)) - sum( ni $ (not (ord(ni) = ord(n))), qt(m,n,ni,p,e)) - sum( ei $ (not (ord(ei) = ord(e))), qc(m,n,p,e,ei)) + qp(m,n,p,e) + sum( ni $ (not (ord(ni) = ord(n))) , qt(m,ni,n,p,e)) + sum( ei $ (not (ord(ei) = ord(e))), l(m,n,ei,e)*qc(m,n,p,ei,e)) =E=  0 ;

maxFlowOrg(m,n,ni,e)..           maxFlowV(n,ni,e) + sum( mi $ ancestor(m,mi), fe(mi,n,ni,e)) - f(m,n,ni,e) =G=  0 ;
maxExF(m,n,ni,e)..               maxExFV(n,ni,e) - fe(m,n,ni,e) =G= 0 ;

maxTransOrg(m,n,e)..             maxTransV(n,e) + sum(mi $ ancestor(m,mi), xe(mi,n,e)) - sum(ei  $ (not (ord(e) = ord(ei))) , l(m,n,ei,e)*x(m,n,ei,e)) =G=0 ;
maxExX(m,n,e)..                  maxExXV(n,e) - xe(m,n,e) =G= 0 ;

marketClearingF(m,n,ni,e)..      - f(m,n,ni,e) + sum( p, qt(m,n,ni,p,e)) =E= 0 ;
marketClearingX(m,n,e,ei)..      - x(m,n,e,ei) + sum( p, qc(m,n,p,e,ei)) =E= 0 ;

model orgC /c1.qs, c2.qp, c3.qt, c4.qc, c5.f, c6.fe, c7.x, c8.xe, maxProd.alpha , massBalance.beta , maxFlowOrg.gamma , maxTransOrg.epsilon, maxExF.delta, maxExX.zeta, marketClearingF.upsilon, marketClearingX.phi  /;

ORGstart = jnow;
solve orgC using MCP;
r = orgC.resUsd;
ORGel = jnow*24*3600 - ORGstart;
elapsed = jnow*24*3600 - starttime;

display ORGel, elapsed, r;
