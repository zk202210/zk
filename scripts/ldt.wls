#!/usr/bin/env -S wolframscript

argv = Rest @ $ScriptCommandLine;

(* arguments: *)
logInput := ToExpression[argv[[1]]];
logRho := ToExpression[argv[[2]]];
logF := ToExpression[argv[[3]]];
s := ToExpression[argv[[4]]];

logL = logInput + logRho;
L = 2^logL;
rho = 2^(- logRho);
F = 2^logF;
r = logInput;

(* distance ref: 2018/828 Aurora p.38 *)
d = Min[{1 - rho, (1 - rho) / 2,  (1 - rho/2)/3}];


(* error ref: TR17-134 *)
error0 := 3 * L / F + (1 - Min[{d, (1-3*rho-2/Sqrt[L])/4}])^s;

(* optimized error*)
error0 := 2 * L / F + (r/2-1) * Sqrt[L]/F + (1 - Min[{d, (1-3*rho-2/Sqrt[L])/4}])^s;


J[e_,x_] := 1 - Sqrt[1 - x(1-e)];
(* error ref: BKS18 https://www.math.toronto.edu/swastik/fri.pdf *)
error1[e_] := 2*logL/e^3/F + (1 + e * logL - Min[{d, J[e, J[e, 1 - rho]]}])^s;


(* error ref: 2020/654 *)
error2[m_] := (m + 0.5)^7 * L^2 / 2 / rho^(3/2) / F + (2*m + 1) * (L + 1) * r / Sqrt[rho] / F + (Sqrt[rho] * (1 + 1 / (2 * m)))^s;


(* conj error ref: 2018/828 Aurora p.57 *)
error3[e_] := 2 * logL^2 / e / F + (1 - d * (1 - e))^s;


Print["bits of security"];

a := FindArgMin[{error3[e], 0 < e < 1}, e][[1]];
Print["Conjecture:\t", a,"\t", N[-Log[2, error3[a]]]];

Print["TR17-134:\t\t\t\t", N[-Log[2, error0]]];

a := FindArgMin[{error1[e], 0 < e < 1}, e][[1]];
Print["BKS18:\t\t", a, "\t", N[-Log[2, error1[a]]]];

a := FindArgMin[{error2[m], m > 3}, m][[1]];
Print["2020/654:\t", a,"\t", N[-Log[2, error2[a]]]];
