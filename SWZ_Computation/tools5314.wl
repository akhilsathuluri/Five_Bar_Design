(* ::Package:: *)

BeginPackage["tools5314`",{"tools`"}]

tools5314::usage = "'tools5314' package is a bundle of all the functions developed as a part of ED5314 course"


CCInt::usage = "p is the distance between the centers of the circles in canonical form and r1, r2 are their respective radii."
FindDeg::usage = "To find the degree of a given polynomial"
Circularity::usage = "To find the circularity of a given polynomial"
fCShalftan::usage = "To find the solutions of a polynomial of the form f(Cos[\[Theta]],Sin[\[Theta]]) using half tangent identity"
fCSeuler::usage = "To find the solutions of a polynomial of the form f(Cos[\[Theta]],Sin[\[Theta]]) using euler's identity"
fCScossin::usage=="To find the solutions of a polynomial of the form f(Cos[\[Theta]],Sin[\[Theta]]) using Cos-Sin identity"
PlanarIntDiv::usage = "Function to solve intersection problem of any two planar curves but enter poly with higher order first"
LinkIntersect::usage = "Given link parameters for a given pose, the module returns if the links intersect or not"
CCIntA::usage = "Gives the common area between two intersecting circles. Returns 0 if not."
WSS1::usage = "Given the parameters l,r,b of a symmetric architecture 5 bar mechanism, it returns the workspace area"
ResidualSylv::usage = "Given two planar curves and the module gives the resultant using Sylvesters method"


Begin["Private`"]


(*Function to find ouut the degree of the polynomial*)
(*Make sure that only one variabe is recognised in the poly given as input and also that no other constant is given recognized as a variable*)
FindDeg[poly_]:=
Module[{\[FormalX]},
(*Convert each term to \[FormalX] times the term and then find the highest power of \[FormalX] to find the degree of the polynomial*)
Return[Max@Exponent[MonomialList[poly]/.Thread[Variables[poly]->\[FormalX]],\[FormalX]]]
];


(*Function handling all the seven cases of two circle intersection problem*)
CCInt[p_, r1_, r2_]:=
Module[{x,y,C1,C2,L12,solx,temp,soly,solution},
	If[p==0&&r1==r2,
		Print[StringForm["CCInt::warning: Both the circles overlap and have infinite solutions"]];
		Return[{{}}];
      ];
	If[p==0&&r1!= r2,
		Print[StringForm["CCInt::warning: Both the circles are concentric and may not have a real solution"]];
		Return[{{}}];
	  ];
	If[p!=0,
		C1 = x^2+y^2-r1^2;
		C2 = (x-p)^2+y^2-r2^2;
		L12 = Simplify[C1-C2];
		solx = Solve[L12==0,x];
		temp = C1/.solx;
		soly = Solve[temp==0,y];
		solution = {{solx[[1]][[1]],soly[[1]][[1]]},{solx[[1]][[1]],soly[[2]][[1]]}};
		Return[{x,y}/.solution];
	  ];
];


(*Function to find the circularity of a given polynomial*)
Circularity[poly_]:=
Module[{X,Y,w,deg,temp,homogenpoly,highestdegterms,pos,circularity},
(*First verify if its a planar curve*)
		If[Dimensions[Variables[poly]][[1]]!=2,Print[StringForm["Circularity::warning: The number of variables is greater than 2. Please input a planar curve"]];
			Return[{{}}];
           ];
		If[Dimensions[Variables[poly]][[1]]== 2,
(*Using FinDeg function to find the degree of the polynomial*)
			deg  = FindDeg[poly];	
(*Homogenizing the polynomial, taking w's common and setting w\[Rule]0 to get the homogenous equation at w=0*)
			temp = Simplify[(w^deg)*(poly/.{Variables[poly][[1]]->X/w,Variables[poly][[2]]->Y/w})];
			homogenpoly = temp/.w->0;
(*Extracting the terms with the highest degree by introducong the \[FormalX] terms along with each variable and then finding the coefficient of \[FormalX]^deg*)
			highestdegterms = Coefficient[homogenpoly/.Thread[Variables[homogenpoly]->\[FormalX]*Variables[homogenpoly]],\[FormalX]^deg];
(*Using Position and FactorList functions to find the position of the terms X^2+Y^2 in the factors list of the given polynomial and hence finding the power of the x^2+y^2 terms*)
(*Return the power = circularity else if there is no X^2+Y^2 factor then return 0*)
			pos = Position[FactorList[highestdegterms],X^2+Y^2];
				If[pos!= {},
					circularity= FactorList[highestdegterms][[pos[[1]][[1]]]][[pos[[1]][[2]]+1]];
					Return[circularity];
                  ];
				If[pos=={},
					circularity=0;
					Return[circularity];
                  ];
         ];
];


fCShalftan[{poly_,\[Theta]_}]:=
Module[{t,temp,tempt,sol, \[Theta]sol},
temp = Together[poly/.{Cos[\[Theta]]->((1-Tan[\[Theta]/2]^2)/(1+Tan[\[Theta]/2]^2)),Sin[\[Theta]]->2*Tan[\[Theta]/2]/(1+Tan[\[Theta]/2]^2)}];
tempt = temp/.Tan[\[Theta]/2]->t;
(*Insert code for verifying the highest degree coefficient\[NotEqual]0. Where m,n are powers of Cos and Sin in f*)
(*If[Coefficient[Numerator[temp],Tan[\[Theta]/2]^(2*m+n)]\[Equal]0, Print[StringForm["fSChalftan::warning: The coefficient of the highest degree is zero";]*)
(*Numerically solving the polynomial to obtain roots*)
sol = NSolve[Numerator[tempt]==0,t];
\[Theta]sol = 2*ArcTan[t]/.sol//N;
Return[\[Theta]sol];
];


fCScossin[{poly_,\[Theta]_}]:=
Module[{remaind,cosval,temp,temp2.sinval,\[Theta]sol,t},
remaind = PolynomialRemainder[poly,(Cos[\[Theta]])^2+(Sin[\[Theta]])^2-1,Cos[\[Theta]]];
cosval = Solve[remaind==0,Cos[\[Theta]]];
temp = poly/.cosval;
temp2 = Together[poly/.Sin[\[Theta]]->t];
sinval = NSolve[Numerator[temp2]==0,t]
\[Theta]sol = ArcTan[(Sin[\[Theta]]/.sinval)/(Cos[\[Theta]]/.cosval/.sinval)];
Return[\[Theta]sol];
];



fCSeuler[{poly_,\[Theta]_}]:=
Module[{t, temp, tempt, sol, \[Theta]sol},
temp = Together[poly/.{Cos[\[Theta]]->(t+1/t)/2,Sin[\[Theta]]->(t-1/t)/(2*I)}];
(*Now we have a polynomial in e^(\[ImaginaryI]\[Theta]) and insert code to solve them assuming *)
(*Insert flags and warnings to handle singularities*)
sol = NSolve[Numerator[temp]==0,t];
\[Theta]sol = ArcTan[(t-1/t)/(I*(t+1/t))]/.sol//N;
Return[\[Theta]sol];
];


(*Function to solve intersection problem of any two planar curves but enter poly with higher order first*)
PlanarIntDiv[poly1_,poly2_,var1_,var2_]:=
Module[{temp,tempvar1sol,poly1var2,solvar2val,solvar1val,n,m,solvec,v1,v2,xvalid,yvalid,ysol,xsol},
(*Finding the reminder of the polynomial division between both the polynomials*)
(*Polynomial division is preferred in this case as we are dividing with a circle we get a linar reminder*)
temp = PolynomialRemainder[poly1,poly2,var1];
xvalid = {};yvalid = {};
If[temp==0,
Print[StringForm["PlanarDivInt::warning::The polynomials have a common factor i.e. a subset of the curves might be overlapping"]];
Return[{}];
];
tempvar1sol = Solve[temp==0,var1];
poly1var2 = poly1/.tempvar1sol;
solvar2val = NSolve[poly1var2==0,var2];
solvar1val = tempvar1sol/.solvar2val;
xsol = var1/.solvar1val;
xsol = Flatten[xsol];
ysol = var2/.solvar2val;
ysol = Flatten[ysol];
For[n=1,n<=Dimensions[solvar2val][[1]],n=n+1,
For[m=1,m<=Dimensions[solvar1val][[1]],m=m+1,
v1 = Chop[poly1/.{var2->ysol[[n]],var1->xsol[[m]]}];
v2 = Chop[poly2/.{var2->ysol[[n]],var1->xsol[[m]]}];
If[v1==0&&v2== 0 ,
AppendTo[yvalid,ysol[[n]]];
AppendTo[xvalid,xsol[[m]]];
(*solvec = Join[solvec,{solvar1val[[m]],solvar2val[[n]]}];*);
];
];
];
Return[{xvalid,yvalid}];
];


LinkIntersect[{x1_,y1_},{x2_,y2_},{x3_,y3_},{x4_,y4_},w1_,w2_,fos_]:=
Module[{a1,b1,a2,b2,pc,\[Alpha],\[Beta],R1,z,z1,R2,z2,e1,e2,param,x,y,\[Theta],p11,p12,p21,p22,flag},
(*Now a1,b1,a2,b2 are the properties of the bounding boxes. Assuming the length is along x axis*)
a1 = Sqrt[(x1-x2)^2+(y1-y2)^2]/2;
b1 = fos*w1/2;
a2 = Sqrt[(x3-x4)^2+(y3-y4)^2]/2;
b2 = fos*w2/2;
(*The global coordinates are assumed to be along one of the links with origin at the mid-point of the link. Now we'll need to transform the other ellipse into this coordinates. \[Alpha] is the angle of inclination of the first ellipse*)
p11 = {x1,y1};
p12 = {x2,y2};
p21 = {x3,y3};
p22 = {x4,y4};
pc = 0.5*(p11+p12);
\[Alpha] = ArcTan[(p12[[2]]-p11[[2]])/(p12[[1]]-p11[[1]])];
\[Beta] = ArcTan[(p22[[2]]-p21[[2]])/(p22[[1]]-p21[[1]])];

R1 = {{Cos[\[Beta]],-Sin[\[Beta]]},{Sin[\[Beta]],Cos[\[Beta]]}};
z = {x-0.5*(x3+x4),y-0.5*(y3+y4)}+pc;
z1 = Inverse[R1].z;
R2 = {{Cos[\[Alpha]],-Sin[\[Alpha]]},{Sin[\[Alpha]],Cos[\[Alpha]]}};
z2 = Inverse[R2].z1;
(*Defining the first ellipse in its local coordinates*)
e1 = ((x/a1)^4+(y/b1)^4)-1;
(*Therefore the second ellipse is of the form*)
e2 = ((z2[[1]]/a2)^4+(z2[[2]]/b2)^4)-1;
(*Now we have defined the bounding boxes for the two links all that is left is to check if they are intersecting. The method followed is to define the points on one of the ellipses parametrically with respect to some loop variable \[Theta] and then descretize \[Theta] from 0 to 2\[Pi] and check if any point is inside the second ellipse. It is obvious that the result would depend on the fineness of descritization so instead of making the discretization smaller we increase the bounding box widths with some factor of safety associated with it.*)
param = {a1*(Cos[\[Theta]]^2)^(1/4),b1*(Sin[\[Theta]]^2)^(1/4)};
(*Loop running *)
For[\[Theta] =0, \[Theta]<=2\[Pi],\[Theta] = \[Theta]+0.01,
If[(e2/.{x->param[[1]],y->param[[2]]})<0\[Or](e2/.{x->-param[[1]],y->param[[2]]})<0\[Or](e2/.{x->param[[1]],y->-param[[2]]})<0\[Or](e2/.{x->-param[[1]],y->-param[[2]]})<0,
flag=1;
Break[];
];
flag=0;
];
If[flag==1,Return[Print[StringForm["You just saved a world from robot link collison by finding one!!!"]]]
];
If[flag==0,Return[Print[StringForm["Awesome!!! The given pose has no interferences."]]]
]
];


CCIntA[r1_,r2_,b_]:=
Module[{\[Alpha],\[Beta],eq1,eq2,a1,a2,\[Alpha]angle,\[Beta]angle},
	If[r1+r2>=b\[And]Abs[r1-r2]<b,
		eq1= r1*Cos[\[Alpha]]+r2*Cos[\[Beta]]-b;
		eq2=  r1*Sin[\[Alpha]]-r2*Sin[\[Beta]];
		\[Beta]angle = eangle[{eq1,eq2},\[Alpha],\[Beta]];
		\[Alpha]angle = eangle[{eq1,eq2},\[Beta],\[Alpha]];
	a1 = r1^2*(\[Alpha]angle[[2]][[1]])-r1*Cos[(\[Alpha]angle[[2]][[1]])]*r1*Sin[(\[Alpha]angle[[2]][[1]])];
	a2 = r2^2*(\[Beta]angle[[2]][[1]])-r2*Cos[(\[Beta]angle[[2]][[1]])]*r2*Sin[(\[Beta]angle[[2]][[1]])];
		Return[(a1+a2)];
		];
If[Abs[r1-r2]>=b,
Print[StringForm["CCint::warning::One of the circle is completely inside the other"]];
Return[\[Pi]*(Min[r1,r2])^2];
];
If[r1+r2<b,Print[StringForm["CCint::warning::The circles dont intersect"]];
		Return[0];
	];
];


WSS1[l_,r_,b_] :=
Module[{a1,a2,a3},
If[2*(l+r)==b,Print[StringForm["WSS1::warning: The links form a straight chain. Not really a manipulator. The area is 0"]];
Return[0];
];
If[b==0,Print[StringForm["WSS1::warning:: The loop closes where it began!"]];
Return[N[(\[Pi]*((l+r)^2-l^2))]];
];
If[b-2*l<= 0,
a1 = CCIntA[l+r,l+r,b];
a2 = 2*CCIntA[l,l+r,b];
a3 = CCIntA[l,l,b];
Print[StringForm["WSS1::info:: The full workspace is a union of two non intersecting areas. The total are is given."]];
Return[N[(a1-a2+a3)]];
];
If[b-2*(r+l)<0\[And]b-(l+r)>l,
Print[StringForm["WSS1::info:: The full workspace is just the intersection of two circles of radii l+r"]];
	Return[N[CCIntA[l+r,l+r,b]]];
];
If[b-2*(l+r)<= 0\[And]b-(l+r)<l\[And]b-2*l>=0,
a1 = CCIntA[l+r,l+r,b];
a2 = 2*CCIntA[l+r,l,b];
Print[StringForm["WSS1::info:: The full workspace is just the intersection of two circles of radii l+r minus the twice the circle circle intersection of radii l and l+r"]];
Return[N[(a1-a2)]];
];
];


ResidualSylv[f1_,g1_,x_]:=
Module[{n1,m1,n, m,v,tempmat1,i,j,Mat1,f,g},
n1 = Exponent[f1,x];
m1 = Exponent[g1,x];
If[n1>= m1,f = f1;g = g1;];

If[m1>n1,f = g1;g = f1;];

n = Exponent[f,x];
m = Exponent[g,x];

tempmat1 = {f,g};
v=1;
While[v<= n-1,
 AppendTo[tempmat1,(x^v)*g];
v++;
];
v=1;
While[v<= m-1,
 AppendTo[tempmat1,(x^v)*f];
v++;
];
Mat1 = ConstantArray[0,{m+n,m+n}];
For[i=1,i<=m+n,i++,
For[j=1,j<m+n,j++,
Mat1[[i,j]] = Coefficient[tempmat1[[i]],x^(m+n-j)];
];
Mat1[[i,m+n]] = (tempmat1[[i]])/.x->0;
];
(*Return[{Det[Mat1]}];*)
Return[{Det[Mat1],Mat1}];
];


End[]
EndPackage[]
