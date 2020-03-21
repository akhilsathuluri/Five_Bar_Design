(* ::Package:: *)

BeginPackage["mathutils`"]
(***********************************************************)
(* Defining the zero-tolerance and imaginary tolerance. *)
(***********************************************************)
zeroTol = 10^-15;
imgTol = 10^-15;
(* Protect the symbol zeroTol, imgTol *)

Protect[zeroTol]; 
Protect[imgTol];
Off[Set::wrsym]  (* Just turn off the "Protect" messages *)!
(***********************************************************)

mathutils::usage ="'mathutils' package provides some utilities which come handy in 
manipulating algebraic or trignometric equations. 

The symbols 'zeroTol' and 'imgTol' are protected.

It contains the following functions :

1) realQ[a, \[Epsilon]]
			Returns 'True', if 'a' is real valued within the imaginary-tolerance '\[Epsilon]' and 'False' otherwise.

2) getCoefficient[exp, \[Theta]]  
			Extracts the coefficients from a linear expression 'exp' in the sine and cosine of one variable '\[Theta]'.

3) negativeQ[x]  
			Returns 'True', if 'x' is a negative real number and 'False' otherwise.

4) zeroQ[x, \[Epsilon]] 
			Returns 'True', if 'x' is less than zero-tolerance '\[Epsilon]' and 'False' otherwise.  
            
5) solveLinTrig1[exp, \[Theta]] 
			Solves the trignometric equation 'exp' == 0 linear in sine and cosine of unknown variable '\[Theta]'.
           
6) solveLinTrig2[{exp1, exp2}, \[Theta]] 
			Solves two simultaneous trignometric equations 'exp1'==0 & 'exp2'==0 linear in sines and cosines of unknown variable '\[Theta]' and returns the eliminant. 

7) solveQuad[exp, x] 
			Solves a quadratic equation 'exp' == 0 in the variable 'x'.

8) isCircle[f, {x,y}]
			Checks if a polynomial equation 'f' == 0 represents a circle in the Cartesian coordinates 'x' and 'y'.

9) isLine[f, {x,y}]
			Checks if a polynomial equation 'f' == 0 in the Cartesian coordinates 'x' and 'y' would represent a line. 

10) intersectCircleLine[fc, fl, {x,y}]
			Finds the points of intersection between a circle 'fc' and a line 'fl' defined in the Cartesian coordinates 'x' and 'y'.

11) intersectCircles[fc1, fc2, {x,y}]
			Finds the points of intersection between two circles 'fc1' and 'fc2' in the Cartesian coordinates 'x' and 'y'.

12) rot2D[\[Theta]]
			Returns the 2-D rotation matrix corresponding to the angle '\[Theta]' CCW.

13) rrDyad[{x1,y1},{x2,y2},{l1,l2},{\[Theta]1,\[Theta]2}]
			Solves the RR-dyad with known points {x1,y1},{x2,y2} and link lengths l1,l2 respectively and returns the solutions for \[Theta]1,\[Theta]2.

14) rpDyad[{{x1,y1},l},{{x2,y2},\[Alpha]},{d,\[Theta]}]
			Solves the RP-dyad with prismatic joint along {x2,y2} with slope angle \[Alpha] and link length l with known point {x1,y1}.

P.S : Type '?<function_name>' for more details."


realQ::usage                     = " realQ[a, \[Epsilon]:\!\(\*SuperscriptBox[\(10\), \(-15\)]\)] : Returns 'True', if 'a' is real valued within an imaginary-tolerace of '\[Epsilon]' and 'False' otherwise. Default tolerance is assigned at \!\(\*SuperscriptBox[\(10\), \(-15\)]\)."

getCoefficient::usage            = " getCoefficient[exp, \[Theta]] : Returns the list containing coefficients of Cos[\[Theta]], Sin[\[Theta]] & constant in the trignometric expression of form: a Cos[\[Theta]]+b Sin[\[Theta]]+c, where a, b, c are real numbers. "

negativeQ::usage                 = " negativeQ[x] : Returns 'True', if 'x' is a negative real number and 'False' otherwise."

zeroQ::usage                     = " zeroQ[x, \[Epsilon]:\!\(\*SuperscriptBox[\(10\), \(-15\)]\)] : Returns 'True', if 'x' is less than specified zero-tolerance '\[Epsilon]' and 'False' otherwise. Default tolerance is assigned at \!\(\*SuperscriptBox[\(10\), \(-15\)]\)."

solveLinTrig1::usage             = " solveLinTrig1[exp, \[Theta], singularityCond:0] : Returns a list of rule(s) containing the solutions for '\[Theta]' that satisfy the trignometric equation of the form: a Cos[\[Theta]]+b Sin[\[Theta]]+c == 0 , where a, b, c are real numbers, with an optional argument that takes 0 or 1, to return the singularity condition as well."
                
solveLinTrig2::usage             = " solveLinTrig2[{exp1,exp2}, \[Theta], returnCS:0] : Returns the solution for '\[Theta]' and the eliminant that satisfies the trignometric equations of the form: a1*Cos[\[Theta]]+b1*Sin[\[Theta]]+c1 == 0 & a2*Cos[\[Theta]]+b2*Sin[\[Theta]]+c2 == 0, where a1, a2, b1, b2, c1, c2 are real numbers, with an optional argument that takes 0 or 1 to return sine and cosine functions of '\[Theta]'."

solveQuad::usage                 = " solveQuad[f, x] : Solves a quadratic equation 'f' == 0 in variable x."

isCircle::usage                  = " isCircle[f, {x,y}] : Checks if the polynomial equation 'f' == 0 represents a circle in the Cartesian coordinates x and y. "

isLine::usage                    = " isLine[f, {x,y}] : Checks if the polynomial 'f' == 0 in the Cartesian coordinates 'x' and 'y' represents a line."

intersectCircleLine::usage       = " intersectCircleLine[fc,fl,{x,y}] : Finds the points of intersection between circle 'fc' and line 'fl' expressed in the Cartesian coordinates x and y. The circle can be represented as an equation 'fc' == 0 or in centre-radius form as {{h,k},r}. Similarly, the line can be represented as an equation 'fl' == 0 or in two-point form {{x1,y1},{x2,y2}}."

intersectCircles::usage          = " intersectCircles[fc1,fc2,{x,y}] : Finds the points of intersection between two circles 'fc1' and 'fc2' expressed in the Cartesian coordinates x and y. The circles can be represented as an equation 'fc' == 0 or in centre-radius form as {{h,k},r}."

rot2D::usage                     = " rot2d[\[Theta]] : Returns the 2-D rotation matrix corresponding the angle \[Theta] CCW."

rrDyad::usage                     = "rrDyad[{x1,y1},{x2,y2},{l1,l2},{\[Theta]1,\[Theta]2}] : Solves the RR-dyad with known points {x1,y1},{x2,y2} and link lengths l1,l2 respectively and returns the solutions for \[Theta]1,\[Theta]2 measured CCW from +ve x-axis at known points."

rpDyad::usage                     = "rpDyad[{{x1,y1},l},{{x2,y2},\[Alpha]},{d,\[Theta]}] : Solves the RP-dyad with prismatic joint along {x2,y2} with slope angle \[Alpha] and link length l with known point {x1,y1}, and returns the solutions for d,\[Theta]."



Begin["Private`"]

(*****************************************************)
(* Returns 'True' if the value of 'a' is Real,   *)
(* with the imaginary-tolerance \[Epsilon].               *)
(*****************************************************)
realQ[a_, \[Epsilon]_:imgTol]:= Module[{},
If[NumericQ[a] && Element[Chop[N[a],\[Epsilon]],Reals],
	Return[True],
	Return[False]
];
];

(********************************************************)
(* Returns the list containing coefficients of Cos[\[Theta]], *)
(* Sin[\[Theta]] & constant in the trignometric expression *)
(* of form: a*Cos[\[Theta]]+b*Sin[\[Theta]]+c *)
(********************************************************)
getCoefficient[\[Eta]_,\[Theta]_]:=
Module[{a,b,c},
	{a,b}=Coefficient[\[Eta],{Cos[\[Theta]],Sin[\[Theta]]}];
	c=TrigExpand[Simplify[\[Eta]-{a,b}.{Cos[\[Theta]],Sin[\[Theta]]}]];

(*Print warning if the coefficients are complex*)
	If[((!(realQ[a]) && NumericQ[a]) || (!(realQ[b]) && NumericQ[b]) || (!(realQ[c]) && NumericQ[c])),
		Print[StringForm["getCoefficient::warning: The coefficients of sine and cosine functions in the equation `` are not real valued.",\[Eta]]];
      ];
	Return[{a,b,c}];
];

(********************************************************)
(* Returns 'True', if 'x' is a negative real number  *)
(* else returns False. *)
(********************************************************)
negativeQ[a_] := Module[{},
     If[realQ[a], 
          If[Chop[N[a],imgTol] < 0, 
             Return[True];
            ];
      ];
     Return[False];
];
(************************************************************)
(* Returns 'True', if absolute value of 'x' is less  *)
(* than zero-tolerance \[Epsilon].Default value of \[Epsilon] is 10^(-15) *)
(************************************************************)
zeroQ[x_, \[Epsilon]_: zeroTol] := Module[{},
If[NumericQ[x] && Abs[x] < \[Epsilon],
	Return[True],
	Return[False]
]
];




(************************************************************)
(* Solves one equation '\[Eta]' linear in the trignometric    *)
(* functions Cos[\[Theta]] and Sin[\[Theta]] for '\[Theta]', with an optional*)
(* argument to return the singularity condition as well. *)
(************************************************************)
solveLinTrig1[\[Eta]_, \[Theta]_, singularityCond_:0] :=
  Module[{a, b, c, \[Psi], \[Phi]},
    
	(* Check for the optional argument to be either 1 or 0. *)
	If[!IntersectingQ[{singularityCond},{0,1}],
		Print[StringForm["solveLinTrig1::error: The optional argument 'SingularityCond' should be either 1 or 0."]];
		Return[{{}}];
	  ];	

   (* Get the coefficients of cosine and sine of the angle '\[Theta]' *)
   {a,b,c}=getCoefficient[\[Eta],\[Theta]];
   (* Check that the coefficients of Cos[\[Theta]] and Sin[\[Theta]] are both not zero *)
    If[ zeroQ[(a^2 + b^2)],
		If[ zeroQ[c]
			,(* a\[Equal]0, b\[Equal]0, c\[Equal]0 *)
			Print[StringForm["solveLinTrig1::warning: Trivial equation. The variable `` can take any value in ``.",\[Theta], \
			Interval[{0,2\[Pi]}]]];
			Return[{{\[Theta]-> Interval[{0,2\[Pi]}]}}];,
			(* Else :: a\[Equal]0, b\[Equal]0, c\[NotEqual]0 *)
			Print[StringForm["solveLinTrig1::error: The input equation is inconsistent. No solution exists."]];
			Return[{{}}];
		];
    ];
   
   (* Solve for the \[Psi]=ArcTan[a,b] *)
   \[Psi] = ArcTan[a, b];
   
   (* Solve for \[Phi]=ArcCos[-C/Sqrt[a^2+b^2]] *)
   \[Phi] = ArcCos[-c/Sqrt[a^2 + b^2]];

   (*Root characterization*)
	If[zeroQ[a^2+b^2-c^2],
		Print[StringForm["solveLinTrig1::warning: Singularity occurs. The roots of `` are real and repeated.",\[Theta]]];
  	];

	If[negativeQ[a^2+b^2-c^2],
		Print[StringForm["solveLinTrig1::warning: The roots of `` are complex.",\[Theta]]];
  	];

	(*Check for the optional argument to return the singularity condition as well.*)
    If[singularityCond==0,
		(* Return the solution \[Psi]+\[Phi], \[Psi]-\[Phi] *)
		Return[{{\[Theta] -> \[Psi] - \[Phi]}, {\[Theta] -> \[Psi] + \[Phi]}}];,
		Return[{{{\[Theta] -> \[Psi] - \[Phi]}, {\[Theta] -> \[Psi] + \[Phi]}},a^2+b^2-c^2}]
	]
	
];

(*****************************************************)
(* Solves two equations '\[Eta]1' and '\[Eta]2 linear in the         *)
(* trignometric functions Cos[\[Theta]1], Sin[\[Theta]1], Cos[\[Theta]2], and  *) 
(* Sin[\[Theta]2] for \[Theta]1 and \[Theta]2 .                                *)
(*****************************************************)
solveLinTrig2[{\[Eta]1_, \[Eta]2_}, \[Theta]_, returnCS_: 0] := 
  Module[{a1, b1, c1, a2, b2, c2, detab, detbc, detac, elim, 
    sol\[Theta]},
    
	(* Check for the optional argument to be either 1 or 0. *)
	If[!IntersectingQ[{returnCS},{0,1}],
		Print[StringForm["solveLinTrig2::error: The optional argument 'returnCS' should be either 1 or 0."]];
		Return[{{}}];
	  ];
	
   {a1, b1, c1} = getCoefficient[\[Eta]1, \[Theta]];
   {a2, b2, c2} = getCoefficient[\[Eta]2, \[Theta]];
   
   detab = (a1*b2 - b1*a2);
   detbc = (b1*c2 - b2*c1);
   detac = (a2*c1 - a1*c2);
   
   If[zeroQ[detab],
	 If[zeroQ[detbc] && zeroQ[detac],
      Print[StringForm["solveLinTrig2::error: The input equations, cannot be solved uniquely for ``, as they are linearly dependent in the sine and cosine functions of ``. Try using solveLinTrig1!",\[Theta],\[Theta]]],
	  Print[StringForm["solveLinTrig2::error: The input equations, cannot be solved for ``, as they are inconsistent in the sine and cosine functions of ``.",\[Theta],\[Theta]]]
		];
	  Return[{{}}];
     ];
   
   elim = TrigExpand[Simplify[detac^2 + detbc^2 - detab^2]];
   
   sol\[Theta] = {{\[Theta] -> ArcTan[detbc/detab, detac/detab]},elim};
   

   (* Check for the optional argument to return the substitutions for cosine and sine of \[Theta]1. *)
	If[returnCS!=0,
		sol\[Theta] = Append[sol\[Theta], Solve[{\[Eta]1, \[Eta]2}=={0,0},{Cos[\[Theta]],Sin[\[Theta]]}]];
	  ];

   Return[sol\[Theta]];
];


(********************************************************)
(* Solves a quadratic equation 'f'\[Equal]0 in variable 'x' *)
(********************************************************)
solveQuad[f_, x_] := 
Module[{n, a, b, c, \[CapitalDelta]},
   
	n = Exponent[f, x];
   
   (* If n \[GreaterEqual]  3, then print an error, and return {} *)
     If[n >= 3,
    Print[StringForm["solveQuad::error: Equation `` == 0 is not a quadratic equation in ``. ", f, x]];
    Return[{{}}];
    ];
   
   (* If n \[Equal] 1, then solve the linear equation, 
   print a warning, and return the solution *)
   If[n == 1,
    (* a*x +b = 0 *)
    {a, b} = Reverse[CoefficientList[f, x]];
    Print[
     StringForm[
      "solveQuad::warning: Equation `` == 0 is a linear equation in ``.", f, x]];
    Return[{{x -> -b/a}}];
    ];
   
   (* If n \[Equal] 0, then f is independent of x. Print error, return {}.  *)
   If[n == 0,
		If[NumericQ[f],
		   Print[StringForm["solveQuad::error: Equation `` == 0 is inconsistent.", f]],
           Print[StringForm["solveQuad::error: Equation `` == 0 is independent of ``.", f, x]]
		  ];
    Return[{{}}];
    ];
   
   (* If n \[Equal] -\[Infinity], then f \[Equal] 0. Print error, return NULL. *)
   If[n == -\[Infinity],
    Print[StringForm["solveQuad::error: Equation `` \[Equal] 0 is trivial with respect to  ``.", f, x]];
    Return[{{x-> Interval[{-\[Infinity],\[Infinity]}]}}];
    ];
   
   (* n \[Equal] 2, then proceed normally *)
   If[n == 2,
    {a, b, c} = Reverse[CoefficientList[f, x]];
    \[CapitalDelta] = b^2 - 4 a*c;
    
    (*Print warning if the coefficients are complex*)
	If[((!(realQ[a]) && NumericQ[a]) || (!(realQ[b]) && NumericQ[b]) || (!(realQ[c]) && NumericQ[c])),
		Print[StringForm["solveQuad::warning: The coefficients of the quadratic equation `` are not real valued.",f]];
      ];

    (* If \[CapitalDelta] = 0 or < 0, then simply print warnings, 
    but return the solutions anyway -- complex roots are allowed, 
    and may be processed appropriately later. *)    
    If[(negativeQ[\[CapitalDelta]]), 
     Print[StringForm["solveQuad::warning: The equation `` == 0 has complex conjugate roots for ``.", f, x]];];
    
    If[(realQ[\[CapitalDelta]] && zeroQ[\[CapitalDelta]]), 
     Print[StringForm["solveQuad::warning: The equation `` == 0 has repeated roots for ``.", f, x]];];
    
		Return[{{x -> (-b + Sqrt[\[CapitalDelta]])/(2 a)},{ x -> (-b - Sqrt[\[CapitalDelta]])/(2 a)}}];
    ];

];


(******************************************************)
(* Checks if the polynomial equation 'f'\[Equal]0 represents     *)
(* a circle in the variables 'x' and 'y'                   *)
(*****************************************************)
isCircle[f_,{x_,y_}]:=
Module[{fc,cfList,cfListNl,ax,ay,rsq},

(* Chopping the coefficients with default zero-tolerance of 10^-15 *)
fc = Chop[N[f],zeroTol];

(*Check for the input 'fc' to be a polynomial in 'x' and 'y'. *)
If[!(PolynomialQ[fc,{x,y}]),
	Print[StringForm["isCircle::error: The equation `` == 0 is not a polynomial in `` and ``.",fc,x,y]];
	Return[False];
  ];

cfList = CoefficientList[fc, {x,y}];

(*Check if the coefficients are real*)
If[ MatrixQ[cfList,NumericQ],
    If[!MatrixQ[cfList,realQ],
       Print[StringForm["isCircle::error: The equation `` == 0 does not represent a circle in `` and ``. The coefficients are not real valued. ",fc,x,y]];
       Return[False];
      ];
  ];

(* Check for the polynomial 'fc' to have a total degree of 2. *)
If[Dimensions[cfList]!={3,3},
	Print[StringForm["isCircle::error: The equation `` == 0 does not represent a circle in `` and ``. The degree of the equation is not 2 in `` or ``. ",fc,x,y,x,y]];
	Return[False];
  ];

cfListNl = cfList[[2;;3]][[All,2;;3]];
(*Check for the polynomial 'fc', not to include the terms 'xy', 'x^2y', and 'xy^2'*)
If[!MatrixQ[cfListNl,zeroQ],
	Print[StringForm["isCircle::error: The equation `` == 0 does not represent a circle in `` and ``. The equation contains the terms of the form \!\(\*SuperscriptBox[\(``\), \(\(m\)\(\\\ \)\)]\)\!\(\*SuperscriptBox[\(``\), \(n\)]\) .",fc,x,y,x,y]];
	Return[False];
  ];

{ax, ay}={cfList[[3,1]],cfList[[1,3]]};
(* Check that the coefficients of the x^2 and y^2 terms are the same. *)
If[ !zeroQ[ax - ay], 
	Print[StringForm["isCircle::error: The equation `` == 0 does not represent a circle in `` and ``. The coefficients of ``^2 and ``^2 are not the same. ",fc,x,y,x,y]];
	Return[False];
  ];

rsq = (cfList[[1,2]]/(2 ax))^2 + (cfList[[2,1]]/(2 ax))^2 - cfList[[1,1]]/(ax);
(* Check for the positiveness of the square of the radius of the circle. *)
If[negativeQ[rsq], 
	Print[StringForm["isCircle::error: The radius of the circle `` == 0 is negative. ",fc]];
    Return[False];
  ];

If[zeroQ[rsq], 
	Print[StringForm["isCircle::warning: The radius of the circle `` == 0 is zero. It is a point circle.",fc]];
  ];

Return[True];
];


(*****************************************************)
(* Check for a polynomial 'f' in 'x', and 'y'   *)
(* to represent a line.                            *)
(****************************************************)
isLine[f_,{x_,y_}]:=
Module[{fl,cfList,exp},

(*Chopping the coefficients with default zero-tolerance of 10^-15*)
fl = Chop[N[f],zeroTol];

(* Check for the input 'fl' to be a polynomial in 'x' and 'y' *)
If[!PolynomialQ[fl,{x,y}],
	Print[StringForm["isLine::error: The equation `` == 0 is not a polynomial in Cartesian coordinates `` and ``.",fl,x,y]];
	Return[False];
];

cfList=CoefficientList[f, {x,y}];

(*Check if the coefficients are real*)
If[ MatrixQ[cfList,NumericQ],
    If[!MatrixQ[cfList,realQ],
      Print[StringForm["isLine::error: The equation `` == 0 does not represent a line in Cartesian coordinates `` and ``. The coefficients are not real valued.",fl,x,y]];
       Return[False];
      ];
  ];

exp = Exponent[fl,{x,y}];
(* Check for the linearity of the polynomial 'fl' in 'x' and 'y'. *)
If[!IntersectingQ[{exp.{1,1}},{2,1}],
	Print[StringForm["isLine::error: The equation `` == 0 does not represent a line in Cartesian coordinates `` and ``. The degree of the equation is not 1 in `` or ``. ",fl,x,y,x,y]];
	Return[False];
];

(* Check for non-linear term 'xy' in the polynomial 'fl'.*)
If[!zeroQ[Coefficient[fl,x y]], 
	Print[StringForm["isLine::error: The equation `` == 0 does not represent a line in Cartesian coordinates `` and ``. The equation has the non-linear term ```` .",fl,x,y,x,y]];
	Return[False];
];

Return[True];
];


(**********************************************)
(* Finds the points of intersection between a line *)
(*'fl'\[Equal]0 and a circle 'fc'\[Equal]0 in the Cartesian coordinates 'x' *)
(* and 'y'                                         *)
(**********************************************)
intersectCircleLine[fc_,fl_,{x_,y_}]:=
Module[{al,bl,cl,var,solvar,circquad,circsol},

(*Check for the input equation 'fc' to be a circle.*)
If[!isCircle[fc,{x,y}],
	Print[StringForm["intersectionCircleLine::error: 'isCircle' returned 'False' for the circle check."]];
	Return[{}];
  ];

(*Check for the input equation 'fl' to be a line.*)
If[!isLine[fl,{x,y}],
	Print[StringForm["intersectionCircleLine::error: 'isLine' returned 'False' for the line check."]];
	Return[{}];
  ];

(* Get the coefficients in the line equation*)
{al,bl}=Coefficient[fl,{x,y}];
 cl=Simplify[fl-{al,bl}.{x,y}];

(* Check if any of the coefficients of x or y in the line equation are zero.*)
(* Solve for one of the variable in the line equation*)
If[zeroQ[al],
	var=x;
	solvar={y->((-al/bl)*x-(cl/bl))};,
	(* else solve for x *)
	var=y;
	solvar={x->((-bl/al)*y-(cl/al))};
  ];

(* Substitute the variable in the circle equation *)
circquad=fc/.solvar;

circsol=solveQuad[circquad,var];

(* Return the results *)
Return[Map[Union[solvar/.#,#]&,circsol]];

];


(**********************************************)
(* Finds the points of intersection between *)
(* a line in two-point form and a circle, represented *)
(* in the centre-radius format, in the Cartesian *) 
(* coordinates 'x' and 'y'                   *)
(**********************************************)
intersectCircleLine[{{xc1_,yc1_},r1_},{{x1_,y1_},{x2_,y2_}},{x_,y_}]:=
Module[{circ1,line1},
circ1=(x-xc1)^2+(y-yc1)^2-r1^2;
line1=(x2-x1)(y-y1)-(y2-y1)(x-x1);

Return[intersectCircleLine[circ1,line1,{x,y}]];
];


(**********************************************)
(* Finds the points of intersection between *)
(* a line "fl==0" and a circle, represented *)
(* in the centre-radius format, in the Cartesian *) 
(* coordinates 'x' and 'y'                   *)
(**********************************************)
intersectCircleLine[{{xc1_,yc1_},r1_},fl_,{x_,y_}]:=
Module[{circ1},
circ1=(x-xc1)^2+(y-yc1)^2-r1^2;

Return[intersectCircleLine[circ1,fl,{x,y}]];
];


(**********************************************)
(* Finds the points of intersection between *)
(* a line in two-point form and a circle "fc==0" in the Cartesian *) 
(* coordinates 'x' and 'y'                   *)
(**********************************************)
intersectCircleLine[fc_,{{x1_,y1_},{x2_,y2_}},{x_,y_}]:=
Module[{line1},
line1=(x2-x1)(y-y1)-(y2-y1)(x-x1);

Return[intersectCircleLine[fc,line1,{x,y}]];
];


(********************************************************)
(* Find the points of intersection of two circles    *)
(* 'fc1'\[Equal]0 and 'fc2'\[Equal]0 defined in the variables 'x' *)
(* and 'y' *)
(********************************************************)
intersectCircles[fc1_,fc2_,{x_,y_}]:=
Module[{a1,a2,circ1,circ2,line,al,bl,cl, sol},

(* Verify that these are indeed circles else return NULL *)
If[!isCircle[fc1,{x,y}],
  Print[StringForm["intersectionCircleLine::error: 'isCircle' returned 'False' for the equation `` == 0 circle check.",fc1]];
  Return[{{}}];
  ];
If[!isCircle[fc2,{x,y}],
  Print[StringForm["intersectionCircleLine::error: 'isCircle' returned 'False' for the equation `` == 0 circle check.",fc2]];
  Return[{{}}];
  ];

a1=Coefficient[fc1,x^2];
a2=Coefficient[fc2,x^2];

(* If the circles are scaled, then convert them to the generic form of the equation. *)
circ1=fc1/a1;
circ2=fc2/a2;


(* Subtract the two circle equations to get a line equation. *)
line=circ1-circ2;

(* Get the coefficients in the line equation*)
{al,bl}=Coefficient[line,{x,y}];
cl=Simplify[line-({al,bl}.{x,y})];

(* Check for conditions for colinearity/concurrence of the circles. *)
If[zeroQ[al] && zeroQ[bl],
	If[!zeroQ[cl],
		Print["intersectCircle::error: The circles are concentric."];
		Return[{{}}];,
		(* else *)
		Print["intersectCircle::error: The circles are colinear."];
		Return[{{}}];
	  ];
  ];

(* Solve the circle-line intersection problem using the module above *)
sol=intersectCircleLine[fc1,line,{x,y}];

(* Return the results *)
Return[sol];
];



(********************************************************)
(* Find the points of intersection of two circles    *)
(* represented using the centre and radius defined   *)
(* in the variables 'x' and 'y'                       *)
(********************************************************)
intersectCircles[{{xc1_,yc1_},r1_},fc2_,{x_,y_}]:=
Module[{circ1},
circ1=(x-xc1)^2+(y-yc1)^2-r1^2;

Return[intersectCircles[circ1,fc2,{x,y}]];


];


(********************************************************)
(* Find the points of intersection of two circles    *)
(* represented using the centre and radius defined   *)
(* in the variables 'x' and 'y'                       *)
(********************************************************)
intersectCircles[fc1_,{{xc2_,yc2_},r2_},{x_,y_}]:=
Module[{circ2},
circ2=(x-xc2)^2+(y-yc2)^2-r2^2;

Return[intersectCircles[fc1,circ2,{x,y}]];


];


(********************************************************)
(* Find the points of intersection of two circles    *)
(* represented using the centre and radius defined   *)
(* in the variables 'x' and 'y'                       *)
(********************************************************)
intersectCircles[{{xc1_,yc1_},r1_},{{xc2_,yc2_},r2_},{x_,y_}]:=
Module[{circ1,circ2},
circ1=(x-xc1)^2+(y-yc1)^2-r1^2;
circ2=(x-xc2)^2+(y-yc2)^2-r2^2;

Return[intersectCircles[circ1,circ2,{x,y}]];


];


(********************************************************)
(* Returns the 2D rotation matrix corresponding the *)
(* angle \[Theta]                                            *)
(* CCW rotation is considered positive.              *)
(********************************************************)
rot2D[\[Theta]_]:= {{Cos[\[Theta]],-Sin[\[Theta]]},{Sin[\[Theta]],Cos[\[Theta]]}};


(*********************************************************)
(* Solves the RRdyad module and returns the solutions *)
(* for the passive variables. All angles are measured *)
(* CCW to +ve x-axis from the known points            *)
(*********************************************************)
rrDyad[{x1_,y1_},{x2_,y2_},{l1_,l2_},{\[Theta]1_,\[Theta]2_}]:=
Module[{\[Eta],sol\[Theta]2,elim\[Theta]2,sol\[Theta]1,sol},

\[Eta]={x1,y1}+l1*{Cos[\[Theta]1],Sin[\[Theta]1]}-l2*{Cos[\[Theta]2],Sin[\[Theta]2]}-{x2,y2};
{sol\[Theta]2,elim\[Theta]2}=solveLinTrig2[\[Eta],\[Theta]2];
sol\[Theta]1=solveLinTrig1[TrigExpand[elim\[Theta]2],\[Theta]1];
sol=((#)\[Union](sol\[Theta]2/.(#)))&/@sol\[Theta]1;

Return[sol];
];


(***********************************************************)
(* Solves the RPdyad module and returns the solutions   *)
(* for the passive variables d,\[Theta]. \[Theta] is measured CCW to *)
(* +ve x-axis from the known point {x1,y1}, while d    *)
(* is measured along the line from known point {x2,y2} *)
(***********************************************************)
rpDyad[{{x1_,y1_},l_},{{x2_,y2_},\[Alpha]_},{d_,\[Theta]_}]:=
Module[{\[Eta],sol\[Theta],elim\[Theta],sold,sol},

\[Eta]={x1,y1}+l*{Cos[\[Theta]],Sin[\[Theta]]}-d*{Cos[\[Alpha]],Sin[\[Alpha]]}-{x2,y2};
{sol\[Theta],elim\[Theta]}=solveLinTrig2[\[Eta],\[Theta]];
sold=solveQuad[elim\[Theta],d];
sol=((#)\[Union](sol\[Theta]/.(#)))&/@sold;

Return[sol];
];


End[]
EndPackage[]
