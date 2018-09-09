 /*-------------------------------------------------------------------*
  * MC1 SAS    Monte Carlo study for ISI conference                   *
  *-------------------------------------------------------------------*
  *  Author:  Stef van Buuren, TNO Leiden <BUUREN@NIPG.TNO.NL>        *
  * Created:  Aug 18, 1993                                            *
  * Revised:  Aug 20, 1993                                            *
  * Version:  0.2                                                     *
  *    Code:  SAS/IML 6.07                                            *
  *     See:  1 van Buuren, S., van Rijckevorsel, J.L.A. &            *
  *             Rubin, D.B. (1993). "Multiple imputation by splines"  *
  *             Submitted to ISI 1993 Meeting in Florence             *
  *-------------------------------------------------------------------*/

proc iml;
    reset log noprint fw=5 linesize=85 nocenter;

/* POLY
 * Polynomial regression
 * IN: y, x1 (n x 1), p (>0, number of polynomial terms)
 * OUT: b (1 x p) sd (1 x p)
 */
start poly(b, covb, y, x1, p);
    x = repeat(x1, 1, p);
    pow = repeat((1:p), nrow(x), 1);
    x = x ## pow;
    xpx = t(x) * x;
    xpy = t(x) * y;
    xpxi = inv(xpx);
    b = xpxi * xpy;
    yhat = x * b;
    resid = y - yhat;
    sse = t(resid) * resid;
    dfe = nrow(y) - p;
    mse = sse / dfe;
    covb = xpxi # mse;
    b = t(b);
finish poly;


    /* Simulation parameters */
    s = 1000;
    n = 300;
    q0 = {1 2};
    beta0 = {1 3};
    sigma0 = {3}; *sigma0 = {0.01 0.5 1 3};
    p0 = {0.3}; *p0 = {0.1 0.2 0.3};
    m = 5;
    t = 25;
    k = 3;
    alpha = 0.95;

    do bi = 1 to ncol(beta0);
      beta = beta0[bi];
      do si = 1 to ncol(sigma0);
        sigma = sigma0[si];
        do pi = 1 to ncol(p0);
          p = p0[pi];
          do qi = 1 to ncol(q0);
            q = q0[qi];
            c_hits = j(1, ncol(q0), 0);
            c_beta0 = c_hits;
            c_beta0[q] = beta;
            do r = 1 to s;

              x = normal(j(n, 1, 0));
              e = normal(j(n, 1, 0)) # sigma;
              y = beta # x ## q + e;

              call poly(b, covb, y, x, ncol(q0));
              sd = t(sqrt(vecdiag(covb)));
              low = b - 1.96 # sd;
              up  = b + 1.96 # sd;
              hit = (low <= c_beta0 & c_beta0 <= up);
              c_hits = c_hits + hit;
              *print c_beta0 b sd hit;

            end;
            print beta sigma p q s c_hits;

          end;
        end;
      end;
    end;
quit;
run;




/*   NSS5.SAS
 *   Fitting logistic regression models to the moscow data.
 *   Date: 6/08/93 SvB
 */

proc format;
    value efff  0 = 'observed' 1 = 'none';
    value ishf  0 = 'success'  1 = 'fail';
    value dibsf 0 = 'none';
    value np    0 = 'none' 1 = 'present';
    value imp  -1 = 'Original' 0 = '5 x Imp';
*proc printto; run;
options linesize=80 pagesize=66 pageno=1;
proc printto file='sas.out' new; run;
libname nsslib '/users/buuren/mos/nss/lib';

title 'Model fitting - LISTWISE DELETION';
data temp;
    set nsslib.data;
    drg1 = .; drg2 = .;
    if drg = 1 then do; drg1 = 0; drg2 = 0; end;
    if drg = 2 then do; drg1 = 1; drg2 = 0; end;
    if drg = 3 then do; drg1 = 0; drg2 = 1; end;

proc logistic;
    model ish = drg1 drg2 dibs dns pok st
    /covb;
run;

title 'Model selection - COMPLETED DATA';
data temp;
    set nsslib.data2;
    if 1 <= imp <= 5;
    drg1 = 0; drg2 = 0;
    if drg = 2 then drg1 = 1;
    if drg = 3 then drg2 = 1;

proc logistic covout outest=parms;
    by imp;
    model ish = drg1 drg2 dibs dns pok st
    /covb;
run;

proc printto; run;


/* Pool results */

*data means;
*    set parms;
*    if _TYPE_ = 'PARMS';
*data vars;
*    set parms;
*    if _TYPE = 'COV';

proc iml;
    reset log nocenter;
    m = 5;
    k = 7;
    use work.parms;
    read all var _num_ into qdak where(_type_='PARMS');
    read all var _num_ into u where(_type_='COV');
    qdak = qdak[,2:k+1];
    u = u[,2:k+1];
    print 'Parameter estimates', qdak;
    print 'Covariance matrix of estimates', u;


/* pool */
    qdak = shape(qdak, m, k);
    qbar = qdak[:,];
    ubar = (u[1:7,] + u[8:14,] + u[15:21,] + u[22:28,] + u[29:35,]) / 5;
    print 'Average covariance matrix Ubar',  ubar;

    d = qdak - repeat(qbar, m, 1);
    b = j(k, k, 0);
    do l = 1 to m;
        b = b + t(d[l,]) * d[l,];
    end;
    b = b / (m-1);
    print 'Between imputation covariance matrix', b;

    t = ubar + (1 + 1/m) # b;
    print 'Total covariance matrix', t;


/* 95% CI */
    b = t(vecdiag(b));
    ubar = t(vecdiag(ubar));
    t = t(vecdiag(t));
    r = (1 + 1/m) # b / ubar;
    v = (m - 1) # ((1 + 1/r) ## 2);
    se = sqrt(t);

    low = qbar - 1.96 # sqrt(t);
    up = qbar + 1.96 # sqrt(t);

    f  = (qbar ## 2) / t;
    p = 1 - probf(f, 1, v);

    gamma = (r + (2/(v+3))) / (r + 1);

    rows = {est, se, sig, low, up, gamma};
    mat =  qbar // se // p // low // up // gamma;
    print 'POOLING RESULTS', mat[rowname=rows];

quit;
run;
