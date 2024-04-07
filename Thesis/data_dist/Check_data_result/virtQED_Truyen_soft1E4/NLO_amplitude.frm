*#-
* choose: 1 (counter terms), 2(self), 3(vert el), 4(vert mu), 5(box s), 6(box u)
#define termNLO "4"
*
*************************                 L.D Truyen 12/7/2020               *************************************************
*
*******************************************************************************************************************************
*************************************************   DECLARE   ***********************************************************************
*******************************************************************************************************************************
*  Parameters and Mandelstam variables
Symbols a,t,s,u,e,mme,mmu,pi;
Symbols EL,ME,MM,ME2,MM2;
Symbols dZe1,dZwf1,dZwf2;
*  Particles momentum
vector k,k1,p,p1,q;
vector k1m,pm,qm;
*  Dirac spinor
function U,V;
Cfunction Ubar,Vbar,
*  Particles mass
m;
*********************
*  N-point function *
*********************
*  4-point function
Cfunction Dget, Dval, D0, D1, D2, D3 , D00, D11, D22, D33, D12, D13, D23;
symbol dd0,dd1,dd2,dd3,dd00,dd11,dd22,dd33,dd12,dd13,dd23,id1,id2;
* box s:
set Q1: qm,k1m,p1,qm,qm,k1m,k1m,p1,p1;
set P1: qm,k1m,p1,k1m,p1,p1,qm,qm,k1m;
* box u:
set Q2: qm,k1m,pm,qm,qm,k1m,k1m,pm,pm;
set P2: qm,k1m,pm,k1m,pm,pm,qm,qm,k1m;

set Di: D1,D2,D3;
set Dii: D11,D22,D33,D12,D13,D23,D12,D13,D23;
set D: D0,D1,D2,D3,D00,D11,D22,D33,D12,D13,D23;
set dd: dd0,dd1,dd2,dd3,dd00,dd11,dd22,dd33,dd12,dd13,dd23;
*  3-point function
Cfunction Cget,Cval,C0,C1,C2,C00,C11,C22,C12;
symbol cc0,cc1,cc2,cc00,cc11,cc22,cc12,ic1,ic2;
set Q3: k,k1,k,k1;
set P3: k,k1,k1,k;
set Q4: p,p1,p,p1;
set P4: p,p1,p1,p;
set Ci: C1,C2;
set Cii: C11,C22,C12,C12;
set C: C0,C1,C2,C00,C11,C22,C12;
set cc: cc0,cc1,cc2,cc00,cc11,cc22,cc12;
*  2-point function
Cfunction Bget, Bval, B0, B1, B00, B11, dB0, dB1;
symbol bb0,bb1,bb00,bb11,dbb0,dbb1,ib1,ib2,ib3,ib4,ib5,ib6;
set B: B0,B1,B00,B11,dB0,dB1;
set bb: bb0,bb1,bb00,bb11,dbb0,dbb1;
*  1-point function
Cfunction Aget, Aval, A0;
symbol aa0,ia1,ia2;
*********************
*  Tensor index
indices i,i1,i2,j,j1,j2,mu,nu,mu1,mu2,nu1,nu2,rho,sigma;
************************************
******   Counter parameter   *******
local [d_e]=e^2/(8*pi^2)*(-2/9+2/3*mme^2*dB0(0,mme^2,mme^2)+B1(0,mme^2,mme^2)/3+B0(0,mme^2,mme^2)/2+2/3*mmu^2*dB0(0,mmu^2,mmu^2)+B1(0,mmu^2,mmu^2)/3+B0(0,mmu^2,mmu^2)/2);
local [d_pe]=4*e^2*mme/(16*pi^2)*(+1/(4*mme)-B1(mme^2,0,mme^2)/(2*mme)-B0(mme^2,0,mme^2)/(2*mme)-mme*dB1(mme^2,0,mme^2)+mme*dB0(mme^2,0,mme^2));
local [d_pu]=4*e^2*mmu/(16*pi^2)*(+1/(4*mmu)-B1(mmu^2,0,mmu^2)/(2*mmu)-B0(mmu^2,0,mmu^2)/(2*mmu)-mmu*dB1(mmu^2,0,mmu^2)+mmu*dB0(mmu^2,0,mmu^2));
*
*******************************************************************************************************************************
**************************************  Feynman Amplitude   ************************************************************************
*******************************************************************************************************************************
*   NLO_Amplitude
********************
local [M_NLO]=
********************
***[box
#if ( `termNLO' == 5 )
*  Box1 Amplitude
+i_*e^4/(16*pi^2)*(-Ubar(1,k1)*g_(1,mu,i,nu)*U(1,k)*Ubar(2,p1)*g_(2,mu,j,nu)*U(2,p)*(sum_(a,1,9,Q1[a](i)*P1[a](j)*Dii[a](t,k.k,s,p1.p1,k1.k1,p.p,0,0,mme^2,mmu^2))+d_(i,j)*D00(t,k.k,s,p1.p1,k1.k1,p.p,0,0,mme^2,mmu^2))+
4*Ubar(1,k1)*g_(1,mu)*U(1,k)*Ubar(2,p1)*g_(2,mu)*U(2,p)*k1.p1*D0(t,k.k,s,p1.p1,k1.k1,p.p,0,0,mme^2,mmu^2)+
Ubar(1,k1)*g_(1,mu)*U(1,k)*Ubar(2,p1)*g_(2,i,j,mu)*U(2,p)*(2*k1(i)*sum_(a,1,3,Q1[a](j)*Di[a](t,k.k,s,p1.p1,k1.k1,p.p,0,0,mme^2,mmu^2)))-
Ubar(1,k1)*g_(1,i,j,mu)*U(1,k)*Ubar(2,p1)*g_(2,mu)*U(2,p)*(2*p1(i)*sum_(a,1,3,Q1[a](j)*Di[a](t,k.k,s,p1.p1,k1.k1,p.p,0,0,mme^2,mmu^2))))
#endif
#if ( `termNLO' == 6 )
*********************
*   Box2 Amplitude
+i_*e^4/(16*pi^2)*(+Ubar(1,k1)*g_(1,mu,i,nu)*U(1,k)*Ubar(2,p1)*g_(2,nu,j,mu)*U(2,p)*(sum_(a,1,9,Q2[a](i)*P2[a](j)*Dii[a](t,k.k,u,p.p,k1.k1,p1.p1,0,0,mme^2,mmu^2))+d_(i,j)*D00(t,k.k,u,p.p,k1.k1,p1.p1,0,0,mme^2,mmu^2))+
4*Ubar(1,k1)*g_(1,mu)*U(1,k)*Ubar(2,p1)*g_(2,mu)*U(2,p)*k1.p*D0(t,k.k,u,p.p,k1.k1,p1.p1,0,0,mme^2,mmu^2)-
Ubar(1,k1)*g_(1,mu)*U(1,k)*Ubar(2,p1)*g_(2,mu,j,i)*U(2,p)*(2*k1(i)*sum_(a,1,3,Q2[a](j)*Di[a](t,k.k,u,p.p,k1.k1,p1.p1,0,0,mme^2,mmu^2)))-
Ubar(1,k1)*g_(1,i,j,mu)*U(1,k)*Ubar(2,p1)*g_(2,mu)*U(2,p)*(2*p(i)*sum_(a,1,3,Q2[a](j)*Di[a](t,k.k,u,p.p,k1.k1,p1.p1,0,0,mme^2,mmu^2))))
#endif
***]box
***[self
#if ( `termNLO' == 2 )
*************************
*   Vaccum1 Amplitude
*************************
-i_*e^4/(4*pi^2*t^2)*Ubar(1,k1)*g_(1,mu)*U(1,k)*Ubar(2,p1)*g_(2,nu)*U(2,p)*
(2*d_(mu,nu)*B00(q.q,mme^2,mme^2)+2*q(mu)*q(nu)*B11(q.q,mme^2,mme^2)+2*q(mu)*q(nu)*B1(q.q,mme^2,mme^2)+d_(mu,nu)*(1/2*(q.q)*B0(q.q,mme^2,mme^2)-A0(mme^2))
*************************
*   Vacuum2 Amplitude 
+2*d_(mu,nu)*B00(q.q,mmu^2,mmu^2)+2*q(mu)*q(nu)*B11(q.q,mmu^2,mmu^2)+2*q(mu)*q(nu)*B1(q.q,mmu^2,mmu^2)+d_(mu,nu)*(1/2*(q.q)*B0(q.q,mmu^2,mmu^2)-A0(mmu^2)))
#endif
***]self
***[vert
#if ( `termNLO' == 3 )
*************************
*   Vertex1 Amplitude
*************************
+i_*e^4/(t*16*pi^2)*Ubar(1,k1)*(-2*g_(1,mu)-2*g_(1,i,mu,j)*(sum_(a,1,4,Q3[a](i)*P3[a](j)*Cii[a](k.k,t,k1.k1,0,mme^2,mme^2))+d_(i,j)*C00(k.k,t,k1.k1,0,mme^2,mme^2)+
k(i)*sum_(a,1,2,Q3[a](j)*Ci[a](k.k,t,k1.k1,0,mme^2,mme^2))+k(j)*sum_(a,1,2,Q3[a](i)*Ci[a](k.k,t,k1.k1,0,mme^2,mme^2))+k(i)*k(j)*C0(k.k,t,k1.k1,0,mme^2,mme^2)+
q(j)*(sum_(a,1,2,Q3[a](i)*Ci[a](k.k,t,k1.k1,0,mme^2,mme^2))+k(i)*C0(k.k,t,k1.k1,0,mme^2,mme^2)))+
4*mme*q(mu)*C0(k.k,t,k1.k1,0,mme^2,mme^2)-2*mme^2*g_(1,mu)*C0(k.k,t,k1.k1,0,mme^2,mme^2)+8*mme*(sum_(a,1,2,Q3[a](mu)*Ci[a](k.k,t,k1.k1,0,mme^2,mme^2))+k(mu)*C0(k.k,t,k1.k1,0,mme^2,mme^2)))*U(1,k)*Ubar(2,p1)*g_(2,mu)*U(2,p)
#endif
#if ( `termNLO' == 4 )
*************************
*   Vertex2 Amplitude
*************************
+i_*e^4/(t*16*pi^2)*Ubar(2,p1)*(-2*g_(2,mu)-2*g_(2,i,mu,j)*(sum_(a,1,4,Q4[a](i)*P4[a](j)*Cii[a](p.p,t,p1.p1,0,mmu^2,mmu^2))+d_(i,j)*C00(p.p,t,p1.p1,0,mmu^2,mmu^2)+
p(i)*sum_(a,1,2,Q4[a](j)*Ci[a](p.p,t,p1.p1,0,mmu^2,mmu^2))+p(j)*sum_(a,1,2,Q4[a](i)*Ci[a](p.p,t,p1.p1,0,mmu^2,mmu^2))+p(i)*p(j)*C0(p.p,t,p1.p1,0,mmu^2,mmu^2)-
q(j)*(sum_(a,1,2,Q4[a](i)*Ci[a](p.p,t,p1.p1,0,mmu^2,mmu^2))+p(i)*C0(p.p,t,p1.p1,0,mmu^2,mmu^2)))-
4*mmu*q(mu)*C0(p.p,t,p1.p1,0,mmu^2,mmu^2)-2*mmu^2*g_(2,mu)*C0(p.p,t,p1.p1,0,mmu^2,mmu^2)+8*mmu*(sum_(a,1,2,Q4[a](mu)*Ci[a](p.p,t,p1.p1,0,mmu^2,mmu^2))+p(mu)*C0(p.p,t,p1.p1,0,mmu^2,mmu^2)))*U(2,p)*Ubar(1,k1)*g_(1,mu)*U(1,k)
#endif
***]vert
***[counter-terms
#if ( `termNLO' == 1 )
*  Counter-term Amplitude
*************************
*+i_*e^2/t*Ubar(1,k1)*g_(1,mu)*U(1,k)*Ubar(2,p1)*g_(2,mu)*U(2,p)*(2*dZe1+dZwf1+dZwf2)
+i_*e^2/t*Ubar(1,k1)*g_(1,mu)*U(1,k)*Ubar(2,p1)*g_(2,mu)*U(2,p)*(2*[d_e]+[d_pe]+[d_pu])
#endif
***]counter-terms
;
*************************
*  LO Amplitude conjugate
*************************
local [M0*] = -i_*e^2*
*     electron line
Ubar(1,k)*g_(1,rho)*U(1,k1)*
*     photon propagator
d_(rho,sigma)/t*
*     muon line
Ubar(2,p)*g_(2,sigma)*U(2,p1);
*************************
*  Square Amplitude
*************************
local [M_NLO^2]=2*[M0*]*[M_NLO];
id k1m = -k1;
id pm = -p;
id qm = -q;
bracket e,t;
print;
.sort 
*************************
*     identities        *
*************************
*  Dirac equation  
id U(i?,k?)*Ubar(i?,k?)=(g_(i,k)+m(k));
id V(i?,k?)*Vbar(i?,k?)=(g_(i,k)-m(k));
***************************************
******   Trace calculation   **********
Trace4,1;
Trace4,2;
*
********************************************
*      Mandelstam variable transform      **
id k?{k1,p}.q = t/2;
id k?{k,p1}.q = -t/2;
id q.q = t;
id k1.k = mme^2-t/2;
id p1.p = mmu^2-t/2;
id k.p = (s-mmu^2-mme^2)/2;
id k1.p1 = (s-mmu^2-mme^2)/2;
id k1.p = (mmu^2+mme^2-u)/2;
id k.p1 = (mmu^2+mme^2-u)/2;
id only s=2*mmu^2+2*mme^2-u-t;
*  mass identification
id k?.k?=m(k)*m(k);
id m(k?{k,k1})=mme;
id m(k?{p,p1})=mmu;
*********************************************
****  Numerical calculation conventions  ****
id D0?D[t](?a) = Dval(dd[t],Dget(?a));
id D0?C[t](?a) = Cval(cc[t],Cget(?a));
id B0?B[t](?a) = Bval(bb[t],Bget(?a));
*id A0(?a) = Aval(aa0,Aget(?a));
argument;
id Dget(t,k.k,s,p1.p1,k1.k1,p.p,0,0,mme^2,mmu^2) = id1;
id Dget(t,k.k,u,p.p,k1.k1,p1.p1,0,0,mme^2,mmu^2) = id2;
id Cget(k.k,t,k1.k1,0,mme^2,mme^2) = ic1;
id Cget(p.p,t,p1.p1,0,mmu^2,mmu^2) = ic2;
id Bget(q.q,mme^2,mme^2) = ib1;
id Bget(q.q,mmu^2,mmu^2) = ib2;
id Bget(mme^2,0,mme^2) = ib3;
id Bget(mmu^2,0,mmu^2) = ib4;
id Bget(0,mme^2,mme^2) = ib5;
id Bget(0,mmu^2,mmu^2) = ib6;
id Aget(mme^2) = ia1;
id Aget(mmu^2) = ia2;
endargument;
.sort
* to FormCalc convention:
id e=EL;
id mme=ME;
id mmu=MM;
argument;
id mme^2=ME2;
id mmu^2=MM2;
endargument;
*
*****   export result  *******
*
bracket e,t,u,pi,dZe1,dZwf1,dZwf2;
*format mathematica;
format doublefortran;
*print +s [M_NLO^2];
print [M_NLO^2];
print [d_pe];
print [d_pu];
print [d_e];
.end  

* convention:
*(q.q,mme^2,mme^2) > (q,mme,mme)
*(k.k,t,k1.k1,0,mme^2,mme^2) > (k,k1,0,mme,mme)
*(p.p,t,p1.p1,0,mmu^2,mmu^2) > (p,p1,0,mmu,mmu)
*(t,k.k,s,p1.p1,k1.k1,p.p,0,0,mme^2,mmu^2) > (-q,-k1,p1,0,0,mme,mmu)
*(t,k.k,u,p.p,k1.k1,p1.p1,0,0,mme^2,mmu^2) > (-q,-k1,-p,0,0,mme,mmu)

id Dget(t,k.k,s,p1.p1,k1.k1,p.p,0,0,mme^2,mmu^2) = id1;
id Dget(t,k.k,u,p.p,k1.k1,p1.p1,0,0,mme^2,mmu^2) = id2;
id Cget(k.k,t,k1.k1,0,mme^2,mme^2) = ic1;
id Cget(p.p,t,p1.p1,0,mmu^2,mmu^2) = ic2;
id Bget(q.q,mme^2,mme^2) = ib1;
id Bget(q.q,mmu^2,mmu^2) = ib2;
id Bget(mme^2,0,mme^2) = ib3;
id Bget(mmu^2,0,mmu^2) = ib4;
id Bget(0,mme^2,mme^2) = ib5;
id Bget(0,mmu^2,mmu^2) = ib6;
id Aget(mme^2) = ia1;
id Aget(mmu^2) = ia2;
