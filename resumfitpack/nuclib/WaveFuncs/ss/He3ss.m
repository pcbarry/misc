(* :Title: He3ss

   :Author:	Sergey Kulagin

   :Package Version: May 29, 2003
   :Updated:         June 27, 2007

   :Summary:

Spectral function of He-3.
Based on Sch\"ultze-Sauer Fortran code ss_sumreg.f.
This package reads the data file of Sch\"ultze-Sauer spectral function
[R. Sch\"ultze and P. Sauer, Phys. Rev. C48 (1993) 38]
and returns the interpolation functions for different channels.
C               The spectral function is based on the three-body    *
C               wave function calculated on a finite grid with      *
C               40*30 mesh points for the (p,q=p_N)-momenta. It is  *
C               again calculated on a finite grid with (30*30) mesh *
C               points for the momentum/energy mesh (p_N,E).        * 

*)

BeginPackage["ss`"]

ssprotected = Unprotect[setHe3ss,dss,css]

setHe3ss::usage="Executing setHe3ss[locpath] sets up SS spectral function of He-3.
The argument locpath points to the location of data file He3ss.dat of values of the spectral function
and the package He3ss.m " 

dss::usage="dss[n,p] the deuteron pole contribution to the function f_n (n=0,1,2)
normalized per 1 proton."

css::usage="css[n,nuc,e,p] the contribution from the two-body continuum
to the function f_n (n=0,1,2) for proton (nuc=1) and neutron (nuc=2) normalized
per 1 particle."

pmin::usage= "Minimal momentum in the momentum mesh."
pmax::usage= "Maximum momentum in the momentum mesh."
emin0::usage= "Minimum energy of isospin 0 state in the energy mesh (deuteron binding energy)."
emin1::usage= "Minimum energy of isospin 1 state in the energy mesh."
emax::usage= "Maximum energy in the energy mesh."

(* Messages: *)

dss::arg1 = "Illegal spectral function type `1` "
css::arg1 = "Illegal spectral function type `1` "
css::arg2 = "Illegal nucleon state `1` "


Begin["`Private`"]

(* Interpolation order *)
intord = 3;

setHe3ss[He3path_String]:=(
Print["Setting up the SS spectral function of He-3..."];

(*
Read spectral function from data file. Note the following rules:

  The spectral function is given for the pair-isospin T=0 and T=1.
  The two-body breakup is given in the first columne, i.e. R0T0(1,J) etc.
  R0T0(2..NERG,J) contains the three-body breakup.
  The normalization is chosen such that the proton and  neutron
  contribution is given by:

      R0P = (F0T0 + F0T1)*4.*PI/3.
      R0N = 2.*F0T1*4.*PI/3.

  with R0P = 2/3
  and  R0N = 1/3

  Since the spectral function is given on a finite mesh the sumrule is
  not exactly fullfilled !
     The mesh was chosen for a particular kinematic situation.
     The momentum mesh is given for Gaussian points.
*)
(* ENERGIES ARE IN MEV AND MOMENTA IN FM^-1 *)

nmesh = 30;
r0t0 = Table[0, {i, 1, nmesh}, {j, 1, nmesh}];
r1t0 = r0t0; r2t0 = r0t0;
r0t1 = r0t0; r1t1 = r0t0; r2t1 = r0t0;

(* ener   is the energy mesh points in MeV.
   nerg   is the number of energy mesh points.
   qmesh  are the momentum mesh points in Fm^-1.
   nq     is the number of momentum mesh points.
   qweigh are the Gassian weights (not used in this code).
   r0t0, r1t0, r2t0 are the points of the corresponding spectral functions
                    with isospin 0.
   r0t1, r1t1, r2t1 are the points of the corresponding spectral functions
                    with isospin 1.
*)

(* BEGIN reading the spectral function from data file *)

rawdat = OpenRead[He3path<>"He3ss.dat"];
nerg = Read[rawdat, Number];
nq = Read[rawdat, Number];
ener = ReadList[rawdat, Number, nerg];
qmesh = ReadList[rawdat, Number, nq];
qweigh = ReadList[rawdat, Number, nq];

Do[(
      Do[ (
          r0t0[[i, j]] = Read[rawdat, Number];
          r1t0[[i, j]] = Read[rawdat, Number];
          r2t0[[i, j]] = Read[rawdat, Number]
          )
        , {j, 1, nq}])
    , {i, 1, nerg}];

Do[(
      Do[ (
          r0t1[[i, j]] = Read[rawdat, Number];
          r1t1[[i, j]] = Read[rawdat, Number];
          r2t1[[i, j]] = Read[rawdat, Number]
          )
        , {j, 1, nq}])
    , {i, 1, nerg}];

Close[rawdat];

(* END reading the spectral function from data file *)
(*
        Interpolation of spectral functions
Since the spectral functions change for a few orders of magnitude,
technically it is more convenient to interpolate the Log of the spectral
functions.
*)

(* Isospin 0 *)
(*
The deuteron in the final state (2-body breakup) corresponds to the 1st
row in r0t0, r1t0, and r2t0.
We treat this state separately.
*)

deut[0] = Table[{qmesh[[j]], r0t0[[1, j]] //Log}, {j, 1, nq}];
deut[1] = Table[{qmesh[[j]], r1t0[[1, j]] //Log}, {j, 1, nq}];
deut[2] = Table[{qmesh[[j]], r2t0[[1, j]] //Log}, {j, 1, nq}];

Do[ dint[in,q_] = Interpolation[deut[in], InterpolationOrder -> intord][q], {in,0,2}];

Clear[in,q];
(* Deuteron pole contribution in He-3 normalized per 1 proton *)

dss[in_, q_]:= If[ MemberQ[{0,1,2},in], ((dint[in,q] //Exp)// Re) /2, 
                                        (Message[dss::arg1, in]; Abort[]) ];

(* Two-nucleon continuum in the final state (3-body breakup): *)

(*
Do[ con[in,0] = Flatten[
Table[{ener[[i]], qmesh[[j]], ssval[in,0][[i, j]]}, {i, 2, nerg}, {j,1, nq}],
1], {in,0,2}];
*)

con[0,0] = Table[{ ener[[i]], qmesh[[j]], r0t0[[i, j]] //Log}, {i, 2, nerg}, {j,1, nq}];
con[1,0] = Table[{ ener[[i]], qmesh[[j]], r1t0[[i, j]] //Log}, {i, 2, nerg}, {j,1, nq}];
con[2,0] = Table[{ ener[[i]], qmesh[[j]], r2t0[[i, j]] //Log}, {i, 2, nerg}, {j,1, nq}];

con[0,0] = Flatten[con[0,0], 1];
con[1,0] = Flatten[con[1,0], 1];
con[2,0] = Flatten[con[2,0], 1];

Clear[e,q];
Do[ cint[in,0, e_,q_] = Interpolation[con[in,0], InterpolationOrder -> intord][e,q], {in,0,2}];

Clear[in];
css0[in_,e_,q_] := Which[
     MemberQ[{0,1,2},in], (cint[in,0,e,q] //Exp)// Re,
     True,                (Message[css::arg1, in]; Abort[])];


(* Isospin 1 *)
(*
The deuteron (two-body breakup) does not contribute in this channel.
Two-nucleon continuum in the final state (3-body breakup):
*)

(*
Do[ con[in,1] = Flatten[
           Table[{ener[[i]], qmesh[[j]], ssval[in,1][[i, j]] }, 
       {i, 2, nerg}, {j,1, nq}],
                      1], 
{in,0,2}];
*)

con[0,1] = Table[{ ener[[i]], qmesh[[j]], r0t1[[i, j]] //Log}, {i, 2, nerg}, {j,1, nq}];
con[1,1] = Table[{ ener[[i]], qmesh[[j]], r1t1[[i, j]] //Log}, {i, 2, nerg}, {j,1, nq}];
con[2,1] = Table[{ ener[[i]], qmesh[[j]], r2t1[[i, j]] //Log}, {i, 2, nerg}, {j,1, nq}];

con[0,1] = Flatten[con[0,1], 1];
con[1,1] = Flatten[con[1,1], 1];
con[2,1] = Flatten[con[2,1], 1];

Clear[e,q];
Do[ cint[in,1, e_,q_] = Interpolation[con[in,1], InterpolationOrder -> intord][e,q], {in,0,2}];

Clear[in,inuc];
css1[in_,e_,q_] := Which[
     MemberQ[{0,1,2},in], (cint[in,1,e,q] //Exp)// Re,
     True,                (Message[css::arg1, in]; Abort[])];

(* Two-nucleon continuum contribution to He-3 spectral function
   normalized per 1 particle 
*)

css[in_,inuc_, e_,q_] := Which[ inuc==1, (css0[in,e,q] + css1[in,e,q]) /2, (* per 1 proton *)
                                inuc==2, 2*css1[in,e,q],                   (* neutron *)
                                True, (Message[css::arg2, inuc]; Abort[])];

(*
The proton and neutron momentum distributions mdp[p] and mdn[p].
These distributions are given by the energy integrals of the spectral
functions with definite isospin:

mdp[p] = 4*Pi*(md0[p]+md1[p])
mdn[p] = 4*Pi*(2*md1[p])

The distributions are normalized as

\int_pmin^pmax dp p^2 mdp(p) = 2
                      mdn(p) = 1
*)

(*   The edge points of the momentum and energy mesh   *)

pmin = First[qmesh];
pmax = Last[qmesh]; (* pmax=6.5; *)
emin0 = ener[[1]];
emin1 = ener[[2]];
emax = Last[ener];

Clear[r0t0,r1t0,r2t0];
Clear[r0t1,r1t1,r2t1];
Clear[con,deut];

Print["*** Done"]
)

End[ ]
Protect[Evaluate[ssprotected]]
EndPackage[ ]

Print["If no error messages, the SS spectral function was successfully loaded.
For more information type ? setHe3ss, dss, css, pmin, pmax, emin0, emin1, emax"];

(* End of He3ss *)
