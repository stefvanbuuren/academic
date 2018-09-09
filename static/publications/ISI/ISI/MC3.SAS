 /*-------------------------------------------------------------------*
  * MC3 SAS    Monte Carlo study for TRIM program                     *
  *-------------------------------------------------------------------*
  *  Author:  Stef van Buuren, TNO Leiden <BUUREN@NIPG.TNO.NL>        *
  * Created:  Aug 18, 1993                                            *
  * Revised:  Sept 20, 1993                                           *
  * Version:  0.3                                                     *
  *    Code:  SAS/IML 6.07                                            *
  *     See:  1 van Buuren, S., van Rijckevorsel, J.L.A. &            *
  *             Rubin, D.B. (1993). "Multiple imputation by splines"  *
  *             Submitted to ISI 1993 Meeting in Florence             *
  *-------------------------------------------------------------------*/

proc iml;
    reset log noprint fw=5 linesize=85 nocenter;
    reset storage = trim;
    load module=_all_;


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

start main;
    errcode = {" if errors >= 0 then do;",
               "    errors = errors + 1;",
               "    if errors > 100 then abort;",
               "    else do; call push(errcode); resume; end;",
               " end;" };
    call push(errcode);
    errors = 0;

    /* Simulation parameters */
    s = 100;
    n = 300;
    q0 = {1 2};
    beta0 = {1 3};
    sigma0 = {3 1 0.5 0.01};
    p0 = {0.1 0.2 0.3};
    m = 5;
    maxit = 20;
    k = 2; /* n of parameters */

    template = j(1, 18, 0);
    create results from template;
    do bi = 1 to ncol(beta0);
      beta = beta0[bi];
      do si = 1 to ncol(sigma0);
        sigma = sigma0[si];
        do pi = 1 to ncol(p0);
          p = p0[pi];
          do qi = 1 to ncol(q0);
            q = q0[qi];
            a_hits = j(1, ncol(q0), 0);
            b_hits = j(1, ncol(q0), 0);
            c_hits = j(1, ncol(q0), 0);
            hits = j(1, ncol(q0), 0);
            c_beta0 = c_hits;
            c_beta0[q] = beta;
            do r = 1 to s;
              free result;
              result = result || bi || si || pi || qi || r;

              /* Generate data */
              x = normal(j(n, 1, 0));
              e = normal(j(n, 1, 0)) # sigma;
              y = beta # x ## q + e;

              /* Complete data analysis for control */
              call poly(b, covb, y, x, ncol(q0));
              sd = t(sqrt(vecdiag(covb)));
              low = b - 1.96 # sd;
              up  = b + 1.96 # sd;
              hit = (low <= c_beta0 & c_beta0 <= up);
              c_hits = c_hits + hit;
              *print b;
              *print covb;
              *  print low up;
              result = result || b || low || up;

              /* Make p percent missing data (MCAR) */
              id1 = uniform(j(n, 1, 0)) < p;
              id2 = floor(1 + 2 # uniform(j(n, 1, 0)));
              data = y || x;
              do ii = 1 to n;
                if id1[ii,] = 1 then data[ii, id2[ii,]] = .;
              end;
              r_in = choose(data = ., 0, 1);
              nval = sum(r_in);
              nmis = nrow(data) * ncol(data) - nval;
              *print nmis;
              *datar = data || r_in || id1 || id2;
              *print datar;

              /* Multiply impute the data */
              do anal = 1 to 2;
                if anal = 1 then knots = {1 1}; /* Rubin-Schafer */
                if anal = 2 then knots = {1 3}; /* spline        */
                type  = {1 1};
                impx = j(nmis, m, 0);
                model = 1;
                meth = 1;
                mns = 0;
                call trim(impx, mns, data, r_in, knots, m, maxit,
                          type, meth, model);

                /* Do m complete-data analysis */
                qdak = j(m, 2, 0);
                free u;
                do i = 1 to m;
                    call subimp(data, r_in, impx, i);
                    call poly(b, covb, data[,1], data[,2], ncol(q0));
                    qdak[i,] = b;
                    u = u // covb;
                end;
               * print anal;
               * print qdak;
               * print u;

                /* Pool parameters */
                qbar = qdak[:,];
                ubar = (u[1:2,]+u[3:4,]+u[5:6,]+u[7:8,]+u[9:10,]) / m;
                *print 'Average covariance matrix Ubar',  ubar;
                d = qdak - repeat(qbar, m, 1);
                b = j(k, k, 0);
                do l = 1 to m;
                  b = b + t(d[l,]) * d[l,];
                end;
                b = b / (m-1);
                *print 'Between imputation covariance matrix', b;
                t = ubar + (1 + 1/m) # b;
                *print 'Total covariance matrix', t;


                /* Compute MI 95% CI */
                b = t(vecdiag(b));
                ubar = t(vecdiag(ubar));
                se = t(sqrt(vecdiag(t)));
                low = qbar - 1.96 # se;
                up  = qbar + 1.96 # se;
                *print low up;

                /* Compute coverages */
                hit = (low <= c_beta0 & c_beta0 <= up);
                if anal = 1 then a_hits = a_hits + hit;
                if anal = 2 then b_hits = b_hits + hit;
                resx = result || anal || b || low || up;
                *print resx;
                append from resx;
              end;

            end;
            print p q sigma c_beta0 a_hits b_hits c_hits;
            save;
          end;
        end;
      end;
    end;
finish main;

call main;
errors = -1;

quit;
