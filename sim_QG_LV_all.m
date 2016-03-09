(* This Mathematica 10.1 script integrates the quantitative genetic Lotka-Volterra equations and saves the results in a text file. The file name to which the data is saved is stored in the variable "outfile"; be sure to change this as needed. *)

ClearAll["Global`*"]; (* Clear all variables and functions *)

outfile = "output_QG_LV_simulation.txt"; (* Path/name of file where data should be stored *)

extinctionthreshold = 0.001; (* Extinction threshold: species below this density will be considered extinct *)
alleethreshold = 10.0^(-20); (* Threshold density where weak Allee effect kicks in *)
alpha[sigmai_, sigmaj_, omega_][mui_, muj_] := omega/(E^((mui - muj)^2 / (2*(sigmai^2 + sigmaj^2) + omega^2)) * Sqrt[2*sigmai^2 + 2*sigmaj^2 + omega^2]); (* Competition coefficient between species i (trait mean: mui; trait standard deviation: sigmai) and species j (trait mean: muj; trait standard deviation: sigmaj), with competition width omega *)
beta[sigmai_, sigmaj_, omega_][mui_, muj_] := (2*(muj - mui) * sigmai^2 * omega) / (E^((mui - muj)^2 / (2*(sigmai^2 + sigmaj^2) + omega^2)) * (2*sigmai^2 + 2*sigmaj^2 + omega^2)^(3/2)); (* Selective pressure on species i from species j *)
vars[t_] := Join[Table[n[i][t], {i, 1, S}], Table[mu[i][t], {i, 1, S}]]; (* The dynamical variables of the model, concatenated in a single vector vars[t]: n[i][t] is the density, mu[i][t] is the mean trait of species i at time t *)
cutoff[n_] := Piecewise[{{3*n^2/alleethreshold^2 - 2*n^3/alleethreshold^3, 0 <= n <= alleethreshold}, {1, n > alleethreshold}}]; (* This function is equal to 1 for n greater than the Allee threshold, to 0 for n negative, and varies smoothly as a cubic function in between. One can think of this as a smoothed step function; its role is to prevent species abundances that fall too low to reach negative values due to numerical error *)
dndt[i_, t_] := n[i][t]*(b[sigma[[i]], theta][mu[i][t]] - Sum[alpha[sigma[[i]], sigma[[j]], omega][mu[i][t], mu[j][t]] * n[j][t], {j, 1, S}]) * cutoff[n[i][t]]; (* Right hand side of equation governing the densities. Notice that it is multiplied by cutoff[n[i][t]] to prevent the abundances from getting negative due to numerical error. This may also be thought of as an effective weak Allee effect reducing the growth rate of very small populations *)
dmudt[i_, t_, hsquare_] := hsquare*(bbar[sigma[[i]], theta][mu[i][t]] - Sum[beta[sigma[[i]], sigma[[j]], omega][mu[i][t], mu[j][t]] * n[j][t], {j, 1, S}]); (* Right hand side of equation governing the trait means *)
rhs[t_, hsquare_] := Join[Table[dndt[i, t], {i, 1, S}], Table[dmudt[i, t, hsquare], {i, 1, S}]]; (* Join the right hand sides of the density- and trait mean equations into a single vector rhs[t,hsquare], where hsquare is the heritability *)
eqs := Thread[D[vars[t], t] == rhs[t, hsquare]]; (* The list of dynamical equations *)
initialcond := Join[Table[n[i][0] == ninit[[i]], {i, 1, S}], Table[mu[i][0] == muinit[[i]], {i, 1, S}]]; (* Initial conditions *)

S = 51; (* Number of species *)
theta = 1/2; (* Width of intrinsic growth function b[sigma,theta][mu] *)
tmax = 10^(10); (* Amount of time to simulate equations for *)
simsets = 100; (* Number of simulation sets to run *)
lowvar = {0.01, 0.05}; (* Lower and upper limits of the (uniform) distribution of intraspecific standard deviations, for the case of low levels of variation *)
mixedvar = {0.01, 0.3}; (* Same for the case of mixed levels of intraspecific variation *)
highvar = {0.1, 0.05}; (* Same for the case of high levels of intraspecific variation *)

SeedRandom[12345]; (* Seed the random number generator, so results are exactly replicable *)
ninit = Table[1, {i, 1, S}]; (* Initial densities: all of them are set to 1 *)
muinittab = RandomReal[{-theta, theta}, {simsets, S}]; (* Pre-generate a table of [simsets] sets of initial conditions for the trait means *)
sigmatab = Table[{Table[0, {S}], Sort[RandomReal[lowvar, S]], Sort[RandomReal[mixedvar, S]],Sort[RandomReal[highvar, S]]}, {simsets}]; (* Pre-generate a table of [simsets] sets of intraspecific standard deviations. Each of the [simsets] sets has four columns: one for zero, one for low, one for mixed, and one for high intraspecific variance. Each of these columns then has as many entries as the number of species S *)

str = OpenWrite[outfile]; (* Create output file *)
WriteString[str, "svar hsquare omega bshape replicate"]; (* Create header line, with the names of the table's columns *)
Do[WriteString[str, " sigma"<>TextString[i]], {i, 1, S}];
Do[WriteString[str, " n"<>TextString[i]], {i, 1, S}];
Do[WriteString[str, " mu"<>TextString[i]], {i, 1, S}];
WriteString[str, "\n"];

Do[
  muinit = muinittab[[rep]]; (* Set initial trait means *)
  Do[
    Do[
      sigma = sigmatab[[rep, svar]]; (* Set intraspecific standard deviations *)
      Clear[b, bbar]; (* Clear intrinsic growth and trait-averaged intrinsic growth functions *)
      If[bshape == 1, (* For bshape==1 the intrinsic growth function is rectangular *)
        If[svar > 1, (* Use these definitions whenever species have intraspecific variation (svar > 1) *)
          b[sigma_, theta_][mu_] := (Erf[(theta - mu) / (Sqrt[2] * sigma)] + Erf[(theta + mu) / (Sqrt[2] * sigma)]) / 2;
          bbar[sigma_, theta_][mu_] := (((1 - E^((2*theta*mu) / sigma^2)) * sigma) / (E^((theta + mu)^2 / (2*sigma^2)) * Sqrt[2*Pi]));
        ];
        If[svar == 1, (* Use these definitions when species have no intraspecific variation (svar = 1) *)
          b[sigma_, theta_][mu_] := Piecewise[{{1, -theta <= mu <= theta}}];
          bbar[sigma_, theta_][mu_] := 0;
        ];
      ];
      If[bshape == 2, (* For bshape==2 the intrinsic growth function is quadratic *)
        If[svar > 1, (* Use these definitions whenever species have intraspecific variation (svar > 1) *)
          b[sigma_, theta_][mu_] := (E^((theta - mu)^2 / (2*sigma^2)) * Sqrt[2/Pi] * (theta - mu + E^((2*theta*mu) / sigma^2) * (theta + mu)) * sigma + E^((theta^2 + mu^2) / sigma^2) * (theta^2 - mu^2 - sigma^2) * (Erf[(theta - mu) / (Sqrt[2]*sigma)] + Erf[(theta + mu) / (Sqrt[2]*sigma)])) / (2*E^((theta^2 + mu^2) / sigma^2) * theta^2);
          bbar[sigma_, theta_][mu_] := -((sigma^2 * ((E^((theta - mu)^2 / (2*sigma^2)) - E^((theta + mu)^2 / (2*sigma^2))) * Sqrt[2/Pi] * sigma + E^((theta^2 + mu^2) / sigma^2) * mu * (Erf[(theta - mu) / (Sqrt[2]*sigma)] + Erf[(theta + mu) / (Sqrt[2]*sigma)]))) / (E^((theta^2 + mu^2) / sigma^2) * theta^2));
        ];
        If[svar == 1, (* Use these definitions when species have no intraspecific variation (svar = 1) *)
          b[sigma_, theta_][mu_] := Piecewise[{{1 - mu^2 / theta^2, -theta <= mu <= theta}}];
          bbar[sigma_, theta_][mu_] := 0;
        ];
      ];
      If[bshape == 3, (* For bshape==3 the intrinsic growth function is triangular *)
        If[svar > 1, (* Use these definitions whenever species have intraspecific variation (svar > 1) *)
          b[sigma_, theta_][mu_] := (-(((-1 + E^((2 * theta * mu) / sigma^2)) * Sqrt[2/Pi] * sigma) / E^((theta + mu)^2 / (2*sigma^2))) + (theta + mu) * (Erf[(theta - mu) / (Sqrt[2]*sigma)] + Erf[(theta + mu) / (Sqrt[2] * sigma)])) / (4 * theta);
          bbar[sigma_, theta_][mu_] := -(sigma / (E^((theta - mu)^2 / (2*sigma^2)) * Sqrt[2*Pi])) + (sigma^2 * (Erf[(theta - mu) / (Sqrt[2]*sigma)] + Erf[(theta + mu) / (Sqrt[2]*sigma)])) / (4 * theta);
        ];
        If[svar == 1, (* Use these definitions when species have no intraspecific variation (svar = 1) *)
          b[sigma_, theta_][mu_] := Piecewise[{{(mu + theta) / (2 * theta), -theta <= mu <= theta}}];
          bbar[sigma_, theta_][mu_] := 0;
        ];
      ];
      Do[
        If[Not[(hsquare > 0) && (svar == 1)], (* Do not solve the biologically impossible case where species have positive heritability but no intraspecific variability (svar = 1) *)
          sol = NDSolve[Join[eqs, initialcond], vars[t], {t, 0, tmax}, MaxSteps -> Infinity]; (* Integrate equations *)
          nf = Table[n[i][t] /. Flatten[sol] /. t -> tmax, {i, 1, S}]; (* Final densities *)
          muf = Table[mu[i][t] /. Flatten[sol] /. t -> tmax, {i, 1, S}]; (* Final trait means *)
          Do[If[nf[[i]] < extinctionthreshold, nf[[i]] = 0], {i, 1, S}]; (* Consider species with density less than extinctionthreshold extinct *)
          (* Write results to file *)
          WriteString[str, TextString[svar - 1] <> " "];
          WriteString[str, TextString[hsquare] <> " "];
          WriteString[str, TextString[omega] <> " "];
          WriteString[str, TextString[bshape] <> " "];
          WriteString[str, TextString[rep - 1]];
          Do[WriteString[str, " " <> TextString[sigma[[i]]]], {i, 1, S}];
          Do[WriteString[str, " " <> TextString[nf[[i]]]], {i, 1, S}];
          Do[WriteString[str, " " <> TextString[muf[[i]]]], {i, 1, S}];
          WriteString[str, "\n"];
        ];
        ,{hsquare, {0, 0.1, 0.5}}, {omega, {0.1, 0.15, 0.2}}
      ];
      ,{svar, {1, 2, 3, 4}}
    ];
    ,{bshape, {1, 2, 3}}
  ]; (* Iterate through all parameter combinations within a simulation set *)
  ,{rep, Range[1, simsets]} (* Repeat for all simulation sets *)
];
Close[str]; (* Close connection to output file *)
