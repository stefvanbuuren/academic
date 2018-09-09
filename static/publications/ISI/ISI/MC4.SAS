 /*-------------------------------------------------------------------*
  * MC4 SAS    Monte Carlo study for TRIM program                     *
  *-------------------------------------------------------------------*
  *  Author:  Stef van Buuren, TNO Leiden <BUUREN@NIPG.TNO.NL>        *
  * Created:  Aug 18, 1993                                            *
  * Revised:  Sept 20, 1993                                           *
  * Version:  0.3                                                     *
  *    Code:  SAS/IML 6.07                                            *
  *     Analyze the results dataset                                   *
  *-------------------------------------------------------------------*/
libname mc '/users/buuren/trim';

/* first read the data */
data mc.data;
    keep par a0 b0 sigma0 prmis power repli method a alow aup hita b blow bup hitb;

    set mc.results;
    par = .;
    if col1 = 1 then par = 1;
    if col1 = 2 then par = 3;
    sigma0 = .;
    if col2 = 1 then sigma0 = 3;
    if col2 = 2 then sigma0 = 1;
    if col2 = 3 then sigma0 = 0.5;
    if col2 = 4 then sigma0 = 0.01;
    prmis = col3 / 10;
    power = col4;
    if power = 1 then do;
       a0 = par;
       b0 = 0;
    end;
    if power = 2 then do;
       a0 = 0;
       b0 = par;
    end;
    repli = col5;
    type = col12;

    if type = 1 then do;
        method = 1;
        alow = col8;
        a    = col6;
        aup  = col10;
        hita = (alow <= a0 & a0 <= aup);
        blow = col9;
        b    = col7;
        bup  = col11;
        hitb = (blow <= b0 & b0 <= bup);
        output;
    end;

    method = type + 1;
    alow = col15;
    a    = (col15 + col17) / 2;
    aup  = col17;
    hita = (alow <= a0 & a0 <= aup);
    blow = col16;
    b    = (col16 + col18) / 2;
    bup  = col18;
    hitb = (blow <= b0 & b0 <= bup);
    output;
run;

*proc sort data=mc.data;
*    by method power prmis sigma0;

proc means data=mc.data n mean maxdec=2;
    class power sigma0 method;
    var hita;
run;
proc means data=mc.data n mean maxdec=2;
    class power sigma0 method;
    var hitb;
run;
