proc sort data=mc.data3;
    by power sigma0 prmis method;
data test;
   set mc.data3;
   if _n_ < 5000;
 run;
 /*-------------------------------------------------------------------*
  * MC& SAS    Monte Carlo study for TRIM program                     *
  *-------------------------------------------------------------------*
  *  Author:  Stef van Buuren, TNO Leiden <BUUREN@NIPG.TNO.NL>        *
  * Created:  Aug 18, 1993                                            *
  * Revised:  Oct, 5, 1993                                            *
  * Version:  0.1                                                     *
  *    Code:  SAS/IML 6.07                                            *
  *     Analyze the result                       *
  *-------------------------------------------------------------------*/
libname mc '/users/buuren/trim';

data

proc means data=mc.data3 n mean maxdec=2 out;
    class power sigma0 methods;
    var hita hitb;
    output out=mc.stats n
               mean(a0)=mu_a mean(b0)=mu_b mean(a)=ave_a mean(b)=ave_b

run;
