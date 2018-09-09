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

*proc sort data=mc.data3;
*    by method power prmis sigma0;

proc means data=mc.data3 n mean maxdec=2;
    class power sigma0 methods;
    var hita;
run;
proc means data=mc.data3 n mean maxdec=2;
    class power sigma0 methods;
    var hitb;
run;
