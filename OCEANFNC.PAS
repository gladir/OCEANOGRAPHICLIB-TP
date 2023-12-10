{ @author: Sylvain Maltais (support@gladir.com)
  @created: 2023
  @website(https://www.gladir.com/CODER/OCEANOGRAPHICLIB-TP)
  @abstract(Target: Turbo Pascal 7, Free Pascal 3.2)
}

Unit OCEANFNC;

INTERFACE

Function Alpha_Wrt_T_Exact(Sa,T,P:Real):Real;
Function Beta_Const_T_Exact(Sa,T,P:Real):Real;
Function Cp_T_Exact(Sa,T,P:Real):Real;
Function Dilution_Coefficient_T_Exact(Sa,T,P:Real):Real;
Function Enthalpy_T_Exact(Sa,T,P:Real):Real;
Function Entropy_From_Pt(Sa,Pt:Real):Real;
Function Entropy_From_T(Sa,T,P:Real):Real;
Function Entropy_Part(Sa,T,P:Real):Real;
Function Kappa_T_Exact(Sa,T,P:Real):Real;
Function Pt_From_T(Sa,T,P,P_Ref:Real):Real;

IMPLEMENTATION

Uses GIBBSFNC;

Const
 gsw_cp0:Real=3991.86795711963;
 gsw_t0:Real=273.15;
 gsw_sso:Real=35.16504;

Var
 gsw_ups:Real;

{@description: Cette fonction permet de calculer le coefficient de
               dilatation thermique de l'eau de mer par rapport … la
               temp‚rature in situ.
 @param: Sa: Salinit‚ Absolue (g/kg)
         T: Temp‚rature in situ (degr‚ Celsius)
         p: Pression de la mer (dbar)
}

Function Alpha_Wrt_T_Exact(Sa,T,P:Real):Real;
Var
 n0,n1:Integer;
Begin
 n0:=0;
 n1:=1;
 Alpha_Wrt_T_Exact:=(Gibbs(n0,n1,n1,Sa,T,P)/Gibbs(n0,n0,n1,Sa,T,P));
End;

{@description: Cette fonction permet de calculer le coefficient de
               contraction saline (haline) de l'eau de mer … temp‚rature
               constante in situ.
 @param: Sa: Salinit‚ Absolue (g/kg)
         T: Temp‚rature in situ (degr‚ Celsius)
         p: Pression de la mer (dbar)
}

Function Beta_Const_T_Exact(Sa,T,P:Real):Real;
Var
 n0,n1:Integer;
Begin
 n0:=0;
 n1:=1;
 Beta_Const_T_Exact:=-Gibbs(n1,n0,n1,Sa,T,P)/Gibbs(n0,n0,n1,Sa,T,P);
End;

{@description: Cette fonction permet de calculer la capacit‚ thermique
               isobare de l'eau de mer.
}

Function Cp_T_Exact(Sa,T,P:Real):Real;
Var
 n0,n2:Integer;
Begin
 n0:=0;
 n2:=2;
 Cp_T_Exact:=-(T+273.15e0)*Gibbs(n0,n2,n0,Sa,T,P);
End;

{@description: Cette fonction permet de calculer le coefficient de
               dilution de l'eau de mer. Le coefficient de dilution de
               l'eau de mer est d‚fini comme la salinit‚ absolue
               multipli‚e par la d‚riv‚e seconde de la fonction de Gibbs
               par rapport … la salinit‚ absolue.
}

Function Dilution_Coefficient_T_Exact(Sa,T,P:Real):Real;
Var
 g08,x,x2,y,z:Real;
Begin
 x2:=SFac*Sa;
 x:=Sqrt(x2);
 y:=T*0.025;
 z:=P*1E-4;
 g08:=2.0*(8103.20462414788+
        y*(2175.341332000392+
        y*(-274.2290036817964+
        y*(197.4670779425016+
        y*(-68.5590309679152+9.98788038278032*y)))-90.6734234051316*z)+
        1.5*x*(-5458.34205214835-980.14153344888*y+
        (4.0/3.0)*x*(2247.60742726704-340.1237483177863*1.25*x+220.542973797483*y)+
        180.142097805543*z)+z*(-219.1676534131548+
        (-16.32775915649044 - 120.7020447884644*z)*z));
 g08:=x2*g08+x*(-7296.43987145382+z*(598.378809221703+
      z*(-156.8822727844005+(204.1334828179377-10.23755797323846*z)*z))+
      y*(-1480.222530425046+
      z*(-525.876123559641+(249.57717834054571-88.449193048287*z)*z)+
      y*(-129.1994027934126+z*(1149.174198007428+
      z*(-162.5751787551336+76.9195462169742*z))+
      y*(-30.0682112585625-1380.9597954037708*z+
      y*(2.626801985426835+703.695562834065*z)))))+
      11625.62913253464+1702.453469893412*y;
 Dilution_Coefficient_T_Exact:=0.25*Sfac*g08;
End;

{@description: Cette fonction permet de calculer l'enthalpie sp‚cifique
               de l'eau de mer.
 @param: Sa: Salinit‚ Absolue (g/kg)
         T: Temp‚rature in situ (degr‚ Celsius)
         p: Pression de la mer (dbar)
}

Function Enthalpy_T_Exact(Sa,T,P:Real):Real;
Var
 n0,n1:Integer;
Begin
 n0:=0;
 n1:=1;
 Enthalpy_T_Exact:=(Gibbs(n0,n0,n0,sa,t,p)-(t+Gsw_T0)*Gibbs(n0,n1,n0,sa,t,p));
End;

{@description: Cette fonction permet de calculer l'entropie sp‚cifique de
               l'eau de mer par temp‚rature potentielle.
}

Function Entropy_From_Pt(Sa,Pt:Real):Real;
Var
 n0,n1:Integer;
 pr0:Real;
Begin
 n0:=0;
 n1:=1;
 pr0:=0.0;
 Entropy_From_Pt:=-Gibbs(n0,n1,n0,Sa,Pt,pr0);
End;

{@description: Cette fonction permet de calculer l'entropie sp‚cifique de
               l'eau de mer.
}

Function Entropy_From_T(Sa,T,P:Real):Real;
Var
 n0,n1:Integer;
Begin
 n0:=0;
 n1:=1;
 Entropy_From_T:=-Gibbs(n0,n1,n0,Sa,T,P);
End;

{@description: Cette fonction permet de calculer l'entropie moins les
               termes d‚pendeat uniquement de SA.
}

Function Entropy_Part(Sa,T,P:Real):Real;
Var
 x2,x,y,z,g03,g08:Real;
Begin
 x2:=Sfac*sa;
 x:=Sqrt(x2);
 y:=t*0.025;
 z:=p*1e-4;
 g03:=z*(-270.983805184062+z*(776.153611613101 + z*(-196.51255088122 +
     (28.9796526294175-2.13290083518327*z)*z)))+
     y*(-24715.571866078+z*(2910.0729080936+
     z*(-1513.116771538718+z*(546.959324647056+
     z*(-111.1208127634436+8.68841343834394*z))))+
     y*(2210.2236124548363+z*(-2017.52334943521+
     z*(1498.081172457456+z*(-718.6359919632359+
       (146.4037555781616-4.9892131862671505*z)*z)))+
     y*(-592.743745734632+z*(1591.873781627888+
     z*(-1207.261522487504+(608.785486935364-105.4993508931208*z)*z))+
     y*(290.12956292128547+z*(-973.091553087975+
     z*(602.603274510125+z*(-276.361526170076+32.40953340386105*z)))+
     y*(-113.90630790850321+y*(21.35571525415769-67.41756835751434*z)+
     z*(381.06836198507096+z*(-133.7383902842754+49.023632509086724*z)))))));
 g08:=x2*(z*(729.116529735046+z*(-343.956902961561+z*(124.687671116248+
      z*(-31.656964386073+7.04658803315449*z))))+
      x*(x*(y*(-137.1145018408982+y*(148.10030845687618+
      y*(-68.5590309679152+12.4848504784754*y)))-
         22.6683558512829*z)+z*(-175.292041186547+
        (83.1923927801819-29.483064349429*z)*z)+
      y*(-86.1329351956084+z*(766.116132004952+
      z*(-108.3834525034224+51.2796974779828*z))+
      y*(-30.0682112585625-1380.9597954037708*z+
      y*(3.50240264723578+938.26075044542*z))))+
      y*(1760.062705994408+y*(-675.802947790203+
      y*(365.7041791005036+y*(-108.30162043765552+12.78101825083098*y)+
      z*(-1190.914967948748+(298.904564555024-145.9491676006352*z)*z))+
      z*(2082.7344423998043+z*(-614.668925894709+
        (340.685093521782-33.3848202979239*z)*z)))+
      z*(-1721.528607567954+z*(674.819060538734+
      z*(-356.629112415276+(88.4080716616-15.84003094423364*z)*z)))));
 Entropy_Part:=-(g03+g08)*0.025;
End;

{@description: Cette fonction permet d'indiquer la compressibilit‚
               isentropique de l'eau de mer.
}

Function Kappa_T_Exact(Sa,T,P:Real):Real;
Var
 n0,n1,n2:Integer;
 g_tt,g_tp:Real;
Begin
 n0:=0;
 n1:=1;
 n2:=2;
 g_tt:=Gibbs(n0,n2,n0,sa,t,p);
 g_tp:=Gibbs(n0,n1,n1,sa,t,p);
 Kappa_T_Exact:=((g_tp*g_tp-g_tt*Gibbs(n0,n0,n2,sa,t,p))/(Gibbs(n0,n0,n1,sa,t,p)*g_tt));
End;

{@description: Cette fonction permet de calculer la temp‚rature potentielle
               de l'eau de mer … partir de la temp‚rature in situ.
 @param: Sa:    Salinit‚ Absolue (g/kg)
         T:     Temp‚rature in situ (degr‚ Celsius)
         p:     Pression de la mer (dbar)
         p_ref: Pression de mer de r‚f‚rence (dbar)
}

Function Pt_From_T(Sa,T,P,P_Ref:Real):Real;
Var
 n0,n2,no_iter:Integer;
 s1,pt,ptm,pt_old,dentropy,dentropy_dt,true_entropy_part:Real;
Begin
 n0:=0;
 n2:=2;
 s1:=sa/gsw_ups;
 pt:=t+(p-p_ref)*(8.65483913395442e-6-s1*1.41636299744881e-6-
     (p+p_ref)*7.38286467135737e-9+t*(-8.38241357039698e-6+
     s1*2.83933368585534e-8+t*1.77803965218656e-8+(p+p_ref)*1.71155619208233e-10));
 dentropy_dt:=gsw_cp0/((gsw_t0+pt)*(1.0-0.05*(1.0-sa/gsw_sso)));
 true_entropy_part:=Entropy_part(sa,t,p);
 For no_iter:=1 to 2 do Begin
  pt_old:=pt;
  dentropy:=Entropy_part(sa,pt_old,p_ref) - true_entropy_part;
  pt:=pt_old-dentropy/dentropy_dt;
  ptm:=0.5*(pt+pt_old);
  dentropy_dt:=-Gibbs(n0,n2,n0,sa,ptm,p_ref);
  pt:=pt_old - dentropy/dentropy_dt;
 End;
 Pt_From_T:=Pt;
End;

BEGIN
 gsw_ups:=gsw_sso/35.0;
END.