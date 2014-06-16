*This dataset creates an example with 16 scenarios and 2 timestages
*The parameter int varies from a 50% reduction to a 5 times increase in stage 2 relative to stage 1.
*The probabilities of each scenario are normally distributed with expected value 2.75.



sets
t time stages /1*2/
m scenario tree nodes /1*17/
nm number of scenarios /1*16/
p producers /1*2/
n nodes /1*2/
*a arcs /0/

*m1(m) /1/
*stage(t) /1/
e types of energy natural gas as 1 and electricity as 2/1*2/
et(e) /2/
ef(e) /1/
d demand sectors /1*2/
;

alias
(p,pi)
(n,ni)
(t,ti)
(m,mi)
(e,ei)
;

scalars
k1  quadratic production cost coeff  /1 /
k2  linear production cost coeff     /1 /
k3  quad expansion cost coeff F     /1 /
k4  quad expansion cost coeff X     /1 /
k5  linear expansion cost coeff F     /1 /
k6  linear expansion cost coeff X     /1 /;

parameter
slp(m,n,d,e) slope of demand function;
slp(m,n,d,e) = 1;

table maxFlowV(n,n,e)        maximum flow on each arc
      1     2
1.2   15   15
2.1   0     0
;


table maxTransV(n,e)        maximum transformation output at each node
    1       2
1   0       0.1
2   0       0.1
;

table maxExXV(n,e)        maximum transformation output capacity expansions at each node at a time
   1     2
1  0     10
2  0     10
;

table maxExFV(n,n,e)        maximum transformation output capacity expansions at each node at a time
     1     2
1.1  0     0
1.2  30    30
2.1  0     0
2.2  0     0
;


set table arcs(n,n)      map of arcs connecting nodes
    1    2
1   0    1
2   0    0
;

set table spno(n,p)      tabel of nodes where a supplier can supply
    1    2
1   1    1
2   1    1
;

table maxProdV(m,n,p,e)
                 1       2
1.1.1            30      0
1.1.2            30      0
1.2.1            30      30
1.2.2            30      30
2.1.1            30      0
2.1.2            30      0
2.2.1            30      30
2.2.2            30      30
3.1.1            30      0
3.1.2            30      0
3.2.1            30      30
3.2.2            30      30
4.1.1            30      0
4.1.2            30      0
4.2.1            30      30
4.2.2            30      30
5.1.1            30      0
5.1.2            30      0
5.2.1            30      30
5.2.2            30      30
6.1.1            30      0
6.1.2            30      0
6.2.1            30      30
6.2.2            30      30
7.1.1            30      0
7.1.2            30      0
7.2.1            30      30
7.2.2            30      30
8.1.1            30      0
8.1.2            30      0
8.2.1            30      30
8.2.2            30      30
9.1.1            30      0
9.1.2            30      0
9.2.1            30      30
9.2.2            30      30
10.1.1           30      0
10.1.2           30      0
10.2.1           30      30
10.2.2           30      30
11.1.1           30      0
11.1.2           30      0
11.2.1           30      30
11.2.2           30      30
12.1.1           30      0
12.1.2           30      0
12.2.1           30      30
12.2.2           30      30
13.1.1           30      0
13.1.2           30      0
13.2.1           30      30
13.2.2           30      30
14.1.1           30      0
14.1.2           30      0
14.2.1           30      30
14.2.2           30      30
15.1.1           30      0
15.1.2           30      0
15.2.1           30      30
15.2.2           30      30
16.1.1           30      0
16.1.2           30      0
16.2.1           30      30
16.2.2           30      30
17.1.1           30      0
17.1.2           30      0
17.2.1           30      30
17.2.2           30      30
;

table l(m,n,e,e)
         1       2
1.1.1    1       0.4
1.1.2    1       1
1.2.1    1       0.4
1.2.2    1       1
2.1.1    1       0.4
2.1.2    1       1
2.2.1    1       0.4
2.2.2    1       1
3.1.1    1       0.4
3.1.2    1       1
3.2.1    1       0.4
3.2.2    1       1
4.1.1    1       0.4
4.1.2    1       1
4.2.1    1       0.4
4.2.2    1       1
5.1.1    1       0.4
5.1.2    1       1
5.2.1    1       0.4
5.2.2    1       1
6.1.1    1       0.4
6.1.2    1       1
6.2.1    1       0.4
6.2.2    1       1
7.1.1    1       0.4
7.1.2    1       1
7.2.1    1       0.4
7.2.2    1       1
8.1.1    1       0.4
8.1.2    1       1
8.2.1    1       0.4
8.2.2    1       1
9.1.1    1       0.4
9.1.2    1       1
9.2.1    1       0.4
9.2.2    1       1
10.1.1   1       0.4
10.1.2   1       1
10.2.1   1       0.4
10.2.2   1       1
11.1.1   1       0.4
11.1.2   1       1
11.2.1   1       0.4
11.2.2   1       1
12.1.1   1       0.4
12.1.2   1       1
12.2.1   1       0.4
12.2.2   1       1
13.1.1   1       0.4
13.1.2   1       1
13.2.1   1       0.4
13.2.2   1       1
14.1.1   1       0.4
14.1.2   1       1
14.2.1   1       0.4
14.2.2   1       1
15.1.1   1       0.4
15.1.2   1       1
15.2.1   1       0.4
15.2.2   1       1
16.1.1   1       0.4
16.1.2   1       1
16.2.1   1       0.4
16.2.2   1       1
17.1.1   1       0.4
17.1.2   1       1
17.2.1   1       0.4
17.2.2   1       1
;

table int(m,n,d,e)
         1               2
1.1.1    10.000000       10.000000
1.1.2    25.000000       30.000000
1.2.1    20.000000       10.000000
1.2.2    65.000000       80.000000
2.1.1    5.000000        5.000000
2.1.2    12.500000       15.000000
2.2.1    10.000000       5.000000
2.2.2    32.500000       40.000000
3.1.1    8.000000        8.000000
3.1.2    20.000000       24.000000
3.2.1    16.000000       8.000000
3.2.2    52.000000       64.000000
4.1.1    11.000000       11.000000
4.1.2    27.500000       33.000000
4.2.1    22.000000       11.000000
4.2.2    71.500000       88.000000
5.1.1    14.000000       14.000000
5.1.2    35.000000       42.000000
5.2.1    28.000000       14.000000
5.2.2    91.000000       112.000000
6.1.1    17.000000       17.000000
6.1.2    42.500000       51.000000
6.2.1    34.000000       17.000000
6.2.2    110.500000      136.000000
7.1.1    20.000000       20.000000
7.1.2    50.000000       60.000000
7.2.1    40.000000       20.000000
7.2.2    130.000000      160.000000
8.1.1    23.000000       23.000000
8.1.2    57.500000       69.000000
8.2.1    46.000000       23.000000
8.2.2    149.500000      184.000000
9.1.1    26.000000       26.000000
9.1.2    65.000000       78.000000
9.2.1    52.000000       26.000000
9.2.2    169.000000      208.000000
10.1.1   29.000000       29.000000
10.1.2   72.500000       87.000000
10.2.1   58.000000       29.000000
10.2.2   188.500000      232.000000
11.1.1   32.000000       32.000000
11.1.2   80.000000       96.000000
11.2.1   64.000000       32.000000
11.2.2   208.000000      256.000000
12.1.1   35.000000       35.000000
12.1.2   87.500000       105.000000
12.2.1   70.000000       35.000000
12.2.2   227.500000      280.000000
13.1.1   38.000000       38.000000
13.1.2   95.000000       114.000000
13.2.1   76.000000       38.000000
13.2.2   247.000000      304.000000
14.1.1   41.000000       41.000000
14.1.2   102.500000      123.000000
14.2.1   82.000000       41.000000
14.2.2   266.500000      328.000000
15.1.1   44.000000       44.000000
15.1.2   110.000000      132.000000
15.2.1   88.000000       44.000000
15.2.2   286.000000      352.000000
16.1.1   47.000000       47.000000
16.1.2   117.500000      141.000000
16.2.1   94.000000       47.000000
16.2.2   305.500000      376.000000
17.1.1   50.000000       50.000000
17.1.2   125.000000      150.000000
17.2.1   100.000000      50.000000
17.2.2   325.000000      400.000000
;

table scenarios(nm,t)
         1       2
1        1       2
2        1       3
3        1       4
4        1       5
5        1       6
6        1       7
7        1       8
8        1       9
9        1       10
10       1       11
11       1       12
12       1       13
13       1       14
14       1       15
15       1       16
16       1       17
;

table ancestor(m,m)
         1       2
1        0       0
2        1       0
3        1       0
4        1       0
5        1       0
6        1       0
7        1       0
8        1       0
9        1       0
10       1       0
11       1       0
12       1       0
13       1       0
14       1       0
15       1       0
16       1       0
17       1       0
;

parameters
dc(m) discount rate      /1  1.000000,2  0.980000,3  0.980000,4  0.980000,5  0.980000,6  0.980000,7  0.980000,8  0.980000,9  0.980000,10  0.980000,11  0.980000,12  0.980000,13  0.980000,14  0.980000,15  0.980000,16  0.980000,17  0.980000/
pb(m) probability of each scenario node /1  1.000000,2  0.026693,3  0.036351,4  0.047367,5  0.059057,6  0.070455,7  0.080426,8  0.087845,9  0.091807,10  0.091807,11  0.087845,12  0.080426,13  0.070455,14  0.059057,15  0.047367,16  0.036351,17  0.026693/
 ;


