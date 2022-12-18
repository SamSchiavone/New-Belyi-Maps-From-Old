ChangeDirectory("~/github/BelyiDB");
load "config.m";
AttachSpec("~/github/Belyi/Code/spec");
AttachSpec("~/github/Gm-Reduce/spec");

function FiberProduct(phi1, phi2)
  // Form fiber product of phi1 and phi2
  f1 := PlaneModel(phi1);
  f2 := PlaneModel(phi2);
  K1 := BaseRing(Parent(f1));
  R<t,u,v> := PolynomialRing(K1,3);
  f1 := Evaluate(f1, [t,u]);
  f2 := Evaluate(f2, [t,v]);
  return Curve(Spec(R), [f1,f2]);
end function;

function FiberProductPlane(phi1, phi2 : ind := 1)
  C := FiberProduct(phi1, phi2);
  fs := DefiningEquations(C);
  assert #fs eq 2;
  R<t,u,v> := Parent(fs[1]);
  K1 := BaseRing(R);
  f3 := Resultant(fs[1], fs[2], R.ind);
  printf "resultant = %o\n", f3;
  S<X,Y> := PolynomialRing(K1,2);
  ev := [S.i : i in [1..ind-1]] cat [S!0] cat [S.i : i in [ind..Rank(S)]];
  f3 := Evaluate(f3, ev);
  return Curve(Spec(S), f3);
end function;

function FindCommonUniverse(evs)
  RRoo := ExtendedReals();
  QQ := Rationals();
  evs2 := [];
    for el in evs do     
    if el cmpeq Infinity() then
      Append(~evs2, RRoo!(el));
    else       
      Append(~evs2, RRoo!(QQ!el));
    end if;
  end for;
  return evs2;
end function;

function VerifyBelyiRamification(phi)
  print "Computing ramification points";
  dpts, dmults := Support(Divisor(Differential(phi)));
  print "Computing ramification values";
  evs := [* Evaluate(phi, RepresentativePoint(el)) : el in dpts *];
  evs := Seqset(FindCommonUniverse(evs));
  printf "Ramified above %o\n", evs;
  return #evs eq 3;
end function;

lab1 := LMFDBLabelToBelyiDBLabel("5T5-5_2.1.1.1_3.2-a" : dot_m := true);
s1 := BelyiDBRead(lab1);
phi1 := s1`BelyiDBBelyiMaps[1];
X1 := s1`BelyiDBBelyiCurves[1];
K1 := BaseField(X1);
//X1 := Curve(ProjectiveSpace(PolynomialRing(K1, 2)));
//KX1<x> := FunctionField(X1);
//phi1 := KX1!(3/2*x^5/(x^5 - 5/2*x^4 - 5/4*x^3 + 45/8*x^2 - 27/8));
lab2 := LMFDBLabelToBelyiDBLabel("3T2-3_2.1_2.1-a" : dot_m := true);
s2 := BelyiDBRead(lab2);
phi2 := s2`BelyiDBBelyiMaps[1];
X2 := s2`BelyiDBBelyiCurves[1];

C3 := FiberProduct(phi1,phi2);
C3;
D3 := FiberProductPlane(phi1, phi2 : ind := 1);
//D3;
//KC3 := FunctionField(C3);
//phi := KC3.1;
//KD3 := FunctionField(D3);
//phi := KD3.1;
VerifyBelyiRamification(phi);



names := BelyiDBFilenames(5);
for name in names do
  s := BelyiDBRead(name);
  phi := s`BelyiDBBelyiMaps[1];
  C := FiberProduct(phi, phi2);
  print C;
  print Genus(C);
end for;

/*
  PlaneModel(phi1);
  f1 := $1;
  PlaneModel(phi2);
  f2 := $1;
  R<t,u,v> := PolynomialRing(K1,3);
  Evaluate(f1, [t,u]);
  f1 := $1;
  Evaluate(f2, [t,v]);
  f2 := $1;
  Resultant(f1, f2, t);
  f3 := $1;
  Parent(f3);
  S<X,Y> := PolynomialRing(K1,2);
  f3 := Evaluate(f3, [0, X, Y]);
  C3 := Curve(Spec(S), f3);
  C3;
  Genus(C3);
  IsIrreducible(f3);
  D3 := Curve(Spec(R), [f1,f2]);
  D3;
  Genus(D3);
*/
