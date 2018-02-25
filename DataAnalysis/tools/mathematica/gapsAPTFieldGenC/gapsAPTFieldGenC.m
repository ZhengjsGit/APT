(* ::Package:: *)

BeginPackage["gapsAPTFieldGenC`"]


(* ::Section:: *)
(*Introduction*)


(* ::Text:: *)
(*The "gapsAPTFieldGenC"  package is a sub-project of the APT code. It helps users to generate electromagnetic field c-source functions. Users can write down the analyitical expressions of field configurations via this package, which will automatically generate the standard c-source for "GAPT_APT_Field" prototype. Through the  inputed potential fields, this package also generates all the functions for different order tensors which is defined in APT.*)


(* ::Section:: *)
(*Function massages*)


(* ::Text:: *)
(*Functions for users:*)


gapsAPTDefST::usage="gapsAPTDefST[corlist,cortype]
	Define the symbols of spacetime coordinates used in the expressions of field functions.
	
	corlist is a 4D list containing symbols for 4D spacetime coordinate; 
	cortype is a string telling the type of spatial coordinate. 
	cortype can be set as \"Cartesian\", or \"Cylindrical\".
	
	Example: gapsAPTDefST[{t,R,Theta,Z},\"Cylindrical\"];"
gapsAPTDefName::usage="gapsAPTDefName[name]
	Define the name of the new electromagnetic field function.
	
	name is a string of the name of field function.

	Example: gapsAPTDefName[\"Tokamak\"];"
gapsAPTDefParameters::usage="gapsAPTDefParameters[paralist]
	Define the extra parameters needed by the electromagnetic function.
	Notice that the scale paramters for electric and magnetic field strength must be 
	\"B0\" and \"E0\". Any other parameter except these two should be defined here.
	The listed parameters should be single float number.

	paralist is a list of strings containing all the extra parameters use in new function.
	
	Example: gapsAPTDefParameters[{\"q\",\"a\",\"R0\"}];"
gapsAPTGenAllFiles::usage="gapsAPTGenAllFiles[A4,maxorder]
	Generate all the tensor (with maximum order of maxorder) by use of A4, the 4D 
	representation of potential field.
	This function outputs two files: 1, name.c, the c source file the field function
	based on the standard of APT project; 2, ADD_Field, the command file used in the 
	ADD script of APT.

	A4 is a 4D list of 4-potential field which is composed by the parameters, the 
	coordinate symbols define by users, and the deflaut parameters B0 and E0.
	maxorder is an integer of the maximum order tensors supported by the new function.

	Example: 
		Atokamak={0,\!\(\*FractionBox[\(B0\\\ R0\\\ Z\), \(2\\\ R\)]\),\!\(\*FractionBox[\(\(-B0\)\\\ \((\*SuperscriptBox[\(Z\), \(2\)] + \*SuperscriptBox[\((R - R0)\), \(2\)])\)\), \(2\\\ q\\\ R\)]\) -\!\(\*FractionBox[\(E0\\\ R0\), \(R\)]\) t, \!\(\*FractionBox[\(\(-B0\)\\\ R0\), \(2\)]\) Log[\!\(\*FractionBox[\(R\), \(R0\)]\)]};
		gapsAPTGenAllFiles[Atokamak,3];
"


(*Hiden functions used for debuging*)
(*
gapsAPTReplaceST::usage="1"
gapsAPTTransCorCYLD::usage=="1"
gapsAPTTransCor::usage=="1"
gapsAPTA2EB::usage=="1"
gapsAPTTensorA::usage=="1"
gapsAPTGradA::usage=="1"
gapsAPTGenCFile::usage="1"
gapsAPTIdx1D2ND::usage="1"
ToMatlab::usage = "ToMatlab[expr], transform mathematica code into matlab form"
*)


(* ::Section:: *)
(*Main package*)


Begin["`Private`"]


(* ::Subsection::Closed:: *)
(*ToMatlab package*)


(*http://library.wolfram.com/infocenter/MathSource/577/*)


ToMatlab[e_] := foldlines[ToMatlabaux[e] <> ";\n"]

ToMatlab[e_, name_] :=
    ToMatlabaux[name] <> "=" <> foldlines[ToMatlabaux[e] <> ";\n"]   /; 
		(!NumberQ[name])

ToMatlab[e_, margin_Integer] :=
    Block[{s},
	SetMargin[margin];
	s = foldlines[ToMatlabaux[e] <> ";\n"];
	RestoreMargin[];
	s]

ToMatlab[e_, name_, margin_Integer] :=
    Block[{s},
	SetMargin[margin];
	s = ToMatlabaux[name] <> "=" <> foldlines[ToMatlabaux[e] <> ";\n"];
	RestoreMargin[];
	s]


(*** Numbers and strings *****************************************************)

ToMatlabaux[s_String] := s

ToMatlabaux[n_Integer] :=
    If[n >= 0, ToString[n], "(" <> ToString[n] <> ")"]

(*ToMatlabaux[r_Rational] := 
    "(" <> ToMatlabaux[Numerator[r]] <> "/" <>
           ToMatlabaux[Denominator[r]] <> ")"*)

ToMatlabaux[r_Rational] := 
    "(" <> ToString[Numerator[r]] <> "/" <>
           ToString[Denominator[r]] <> ")"


ToMatlabaux[r_Real] := 
    Block[{a = MantissaExponent[r]},
        If[r >= 0,
            ToString[N[a[[1]],18]] <> "E" <> ToString[a[[2]]],
            "(" <> ToString[N[a[[1]],18]] <> "E" <> ToString[a[[2]]] <> ")"]]

ToMatlabaux[I] := "sqrt(-1)";

ToMatlabaux[c_Complex] :=
    "(" <>
    If[Re[c] === 0,
        "",
        ToMatlabaux[Re[c]] <> "+"] <>
    If[Im[c] === 1,
        "sqrt(-1)",
        "sqrt(-1)*" <> ToMatlabaux[Im[c]] ] <> ")"



(*** Lists, vectors and matrices *********************************************)


numberMatrixQ[m_] := MatrixQ[m] && (And @@ Map[numberListQ,m])

numberListQ[l_] := ListQ[l] && (And @@ Map[NumberQ,l])


numbermatrixToMatlab[m_] :=
    Block[{i, s=""}, 
	For[i=1, i<=Length[m], i++,
	    s = s <> numbermatrixrow[m[[i]]];    
	    If[i < Length[m], s = s <> ";"]];
	s]

numbermatrixrow[l_] :=
    Block[{i, s=""},
	For[i=1, i<=Length[l], i++, 
	    s = s <> ToMatlabaux[l[[i]]];
	    If[i < Length[l], s = s <> ","]];
	s]


ToMatlabaux[l_List /; MatrixQ[l]] :=
    If[numberMatrixQ[l],
	"[" <> numbermatrixToMatlab[l] <> "]",
	"[" <> matrixToMatlab[l] <> "]"]

matrixToMatlab[m_] :=
    If[Length[m] === 1, 
        ToMatlabargs[m[[1]]],
        ToMatlabargs[m[[1]]] <> ";" <>
            matrixToMatlab[ argslistdrop[m] ] ]

ToMatlabaux[l_List] := "[" <> ToMatlabargs[l] <> "]"



(*** Symbols *****************************************************************)

ToMatlabaux[Colon] = ":"
ToMatlabaux[Abs] = "abs"
ToMatlabaux[Min] = "min"
ToMatlabaux[Max] = "max"
ToMatlabaux[Sin] = "sin"
ToMatlabaux[Cos] = "cos"
ToMatlabaux[Tan] = "tan"
ToMatlabaux[Cot] = "cot"
ToMatlabaux[Csc] = "csc"
ToMatlabaux[Sec] = "sec"
ToMatlabaux[ArcSin] = "asin"
ToMatlabaux[ArcCos] = "acos"
ToMatlabaux[ArcTan] = "atan"
ToMatlabaux[ArcCot] = "acot"
ToMatlabaux[ArcCsc] = "acsc"
ToMatlabaux[ArcSec] = "asec"
ToMatlabaux[Sinh] := "sinh"
ToMatlabaux[Cosh] := "cosh"
ToMatlabaux[Tanh] := "tanh"
ToMatlabaux[Coth] := "coth"
ToMatlabaux[Csch] := "csch"
ToMatlabaux[Sech] := "sech"
ToMatlabaux[ArcSinh] := "asinh"
ToMatlabaux[ArcCosh] := "acosh"
ToMatlabaux[ArcTanh] := "atanh"
ToMatlabaux[ArcCoth] := "acoth"
ToMatlabaux[ArcCsch] := "acsch"
ToMatlabaux[ArcSech] := "asech"
ToMatlabaux[Log] := "log"
ToMatlabaux[Exp] := "exp"
ToMatlabaux[MatrixExp] := "expm"
ToMatlabaux[Pi] := "pi"
ToMatlabaux[E] := "exp(1)"
ToMatlabaux[True] := "1"
ToMatlabaux[False] := "0"

ToMatlabaux[e_Symbol] := ToString[e]



(*** Relational operators ****************************************************)

ToMatlabaux[e_ /; Head[e] === Equal] :=
    ToMatlabrelop[ argslist[e], "=="]
ToMatlabaux[e_ /; Head[e] === Unequal] :=
    ToMatlabrelop[ argslist[e], "~="]
ToMatlabaux[e_ /; Head[e] === Less] :=
    ToMatlabrelop[ argslist[e], "<"]
ToMatlabaux[e_ /; Head[e] === Greater] :=
    ToMatlabrelop[ argslist[e], ">"]
ToMatlabaux[e_ /; Head[e] === LessEqual] :=
    ToMatlabrelop[ argslist[e], "<="]
ToMatlabaux[e_ /; Head[e] === GreaterEqual] :=
    ToMatlabrelop[ argslist[e], ">="]
ToMatlabaux[e_ /; Head[e] === And] :=
    ToMatlabrelop[ argslist[e], "&"]
ToMatlabaux[e_ /; Head[e] === Or] :=
    ToMatlabrelop[ argslist[e], "|"]
ToMatlabaux[e_ /; Head[e] === Not] :=
    "~(" <> ToMatlabaux[e[[1]]] <> ")"

ToMatlabrelop[e_, o_] :=
    If[Length[e] === 1, 
        "(" <> ToMatlabaux[e[[1]]] <> ")",
        "(" <> ToMatlabaux[e[[1]]] <> ")" <> o <>
         ToMatlabrelop[ argslistdrop[e], o] ]

relopQ[e_] := MemberQ[{Equal, Unequal, Less, Greater, LessEqual,
    GreaterEqual, And, Or, Not}, Head[e]]



(*** Addition, multiplication and powers *************************************)

ToMatlabaux[e_ /; Head[e] === Plus] :=
    If[relopQ[e[[1]]],
        "(" <> ToMatlabaux[e[[1]]] <> ")",
        ToMatlabaux[e[[1]]] ] <>
    "+" <>
        If[Length[e] === 2,
            If[relopQ[e[[2]]],
                "(" <> ToMatlabaux[e[[2]]] <> ")",
                ToMatlabaux[e[[2]]] ],
            ToMatlabaux[ dropfirst[e] ]]

ToMatlabaux[e_ /; Head[e] === Times] :=
    If[Head[e[[1]]] === Plus,
        "(" <> ToMatlabaux[e[[1]]] <> ")",
        ToMatlabaux[e[[1]]] ] <>
    ".*" <>
        If[Length[e] === 2,
            If[Head[e[[2]]] === Plus,
                "(" <> ToMatlabaux[e[[2]]] <> ")",
                ToMatlabaux[e[[2]]] ],
            ToMatlabaux[ dropfirst[e] ]]

ToMatlabaux[e_ /; Head[e] === Power] :=
    If[Head[e[[1]]] === Plus || Head[e[[1]]] === Times || Head[e[[1]]] === Power,
        "(" <> ToMatlabaux[e[[1]]] <> ")",
        ToMatlabaux[e[[1]]] ] <>
    ".^" <>
        If[Length[e] === 2,
            If[Head[e[[2]]] === Plus || Head[e[[2]]] === Times || Head[e[[2]]] === Power,
                "(" <> ToMatlabaux[e[[2]]] <> ")",
                ToMatlabaux[e[[2]]] ],
            ToMatlabaux[ dropfirst[e] ]]



(*** Special cases of functions **********************************************)

ToMatlabaux[Rule[_,r_]] := ToMatlabaux[r]

ToMatlabaux[Log[10, z_]] := "log10(" <> ToMatlabaux[z] <> ")"
ToMatlabaux[Log[b_, z_]] :=
    "log(" <> ToMatlabaux[z] <> ")./log(" <> ToMatlabaux[b] <> ")"

ToMatlabaux[Power[e_, 1/2]] := "sqrt(" <> ToMatlabaux[e] <> ")"
ToMatlabaux[Power[E, z_]] := "exp(" <> ToMatlabaux[z] <> ")"

ToMatlabaux[If[test_, t_, f_]] :=
    Block[{teststr = ToMatlabaux[test]},
        "((" <> teststr <> ").*(" <> ToMatlabaux[t] <> ")+(~("
             <> teststr <> ")).*(" <> ToMatlabaux[f] <> "))"]

ToMatlabaux[e__ /; (Head[e] === Max || Head[e] == Min)] :=
    ToMatlabaux[Head[e]] <> "(" <>
        If[ Length[e] === 2,
            ToMatlabargs[e] <> ")",
            ToMatlabaux[e[[1]]] <> "," <> ToMatlabaux[dropfirst[e]] <> ")"]

ToMatlabaux[Colon[a_,b_]] :=
    "((" <> ToMatlabaux[a] <> "):(" <> ToMatlabaux[b] <> "))"
ToMatlabaux[Colon[a_,b_,c_]] :=
    "((" <> ToMatlabaux[a] <> "):(" <> ToMatlabaux[b] <> 
        "):(" <> ToMatlabaux[c] <> "))"



(*** General functions *******************************************************)

ToMatlabaux[e_] :=
    ToMatlabaux[Head[e]] <> "(" <>
        ToMatlabargs[ argslist[e] ] <> ")"

ToMatlabargs[e_] :=
    If[Length[e] === 1, 
        ToMatlabaux[e[[1]]],
        ToMatlabaux[e[[1]]] <> "," <>
            ToMatlabargs[ argslistdrop[e] ] ]



(*** Argument lists **********************************************************)

(*** argslist returns a List of the arguments ***)
argslist[e_] :=
    Block[{ARGSLISTINDEX}, Table[ e[[ARGSLISTINDEX]],
        {ARGSLISTINDEX, 1, Length[e]}]]

(*** argslistdrop returns a List of all arguments except the first one ***)
argslistdrop[e_] :=
    Block[{ARGSLISTINDEX}, Table[ e[[ARGSLISTINDEX]], 
        {ARGSLISTINDEX, 2, Length[e]}]]

(*** dropfirst is like argslistdrop but retains the original Head ***)
dropfirst[e_] :=
    e[[ Block[{i}, Table[i, {i,2,Length[e]}]] ]]



(*** Folding long lines ******************************************************)

MARGIN = 66
MARGINS = {}

SetMargin[m_] := (MARGINS = Prepend[MARGINS, MARGIN]; MARGIN = m; MARGINS)

RestoreMargin[] := 
    If[Length[MARGINS] > 0,
	MARGIN = MARGINS[[1]];
	MARGINS = Drop[MARGINS, 1]]		

foldlines[s_String] :=
    Block[{cut, sin=s, sout=""},
	While[StringLength[sin] >= MARGIN, 
	    cut = findcut[sin];
	    If[cut > 0,		
		sout = sout <> StringTake[sin,cut] <> " ...\n  ";
		sin = StringDrop[sin,cut],
		(* else *)
		sout = sout <> StringTake[sin,MARGIN];
		sin = StringDrop[sin,MARGIN]]];
	sout <> sin]

findcut[s_String] :=
    Block[{i=MARGIN}, 
        While[i > 0 &&
              !MemberQ[{";", ",", "(", ")", "+", "*", " "}, StringTake[s,{i}]],
            i--];
        i]


(* ::Subsection:: *)
(*Global variables*)


ClearAll[gapsST];ClearAll[gapsName];ClearAll[gapsParameters];ClearAll[gapsCorType];
gapsST;gapsName;gapsParameters;gapsCorType;


(* ::Subsection:: *)
(*Definition functions*)


gapsAPTDefST[ST_,type_]:=Module[{},gapsST=ST;gapsCorType=type];
gapsAPTDefName[x_]:=Module[{},gapsName=x];
gapsAPTDefParameters[x_]:=Module[{},gapsParameters=x];


(* ::Subsection:: *)
(*Derivation functions*)


gapsAPTReplaceST[A_]:=A/.{gapsST[[1]]->t,gapsST[[2]]->x1,gapsST[[3]]->x2,gapsST[[4]]->x3};


gapsAPTIdx1D2ND[idx1D_,dim_,order_]:=Module[
{i,idxND,J1=idx1D,J2=0},
	idxND={};
	For[i=1,i<order,i++,
		J2=Mod[J1,dim];
		idxND=Append[idxND,J2];
		J1=Floor[J1/dim];
	];
	Flatten[Append[idxND,J1]]
]


gapsAPTTransCorCYLD[A_]:=Module[
{M={{xx/Sqrt[xx^2+yy^2],-yy/Sqrt[xx^2+yy^2],0},{yy/Sqrt[xx^2+yy^2],xx/Sqrt[xx^2+yy^2],0},{0,0,1}},AA,Amid},
	Amid=A/.{t->tt,x1->Sqrt[xx^2+yy^2],x2->ArcTan[xx,yy],x3->zz};
	AA=Amid[[2;;4]];
	AA=M.AA;
	Simplify[{Amid[[1]],AA[[1]],AA[[2]],AA[[3]]}]
];


gapsAPTTransCor[A_,type_]:=Switch[
type,
	"Cartesian",
		Simplify[A/.{t->tt,x1->xx,x2->yy,x3->zz}],
	"Cylindrical",
		gapsAPTTransCorCYLD[A],
	_,
		Print["Wrong Type of inputted coordinates!"];
];


gapsAPTA2EB[A4_]:=Module[
{E,B,A=A4[[2;;4]],phi=A4[[1]]},
	B=Simplify[Curl[A,{xx,yy,zz}]];
	E=Simplify[-Grad[phi,{xx,yy,zz}]-D[A,{tt}]];
	Join[E,B]
];


gapsAPTGradA[A4_,order_]:=Module[
{i,Dim=4,out=A4,idx,out1D={}},
	For[i=1,i<order,i++,
		out=Grad[out,{tt,xx,yy,zz}]
	];
	out=Simplify[out];
	For[i=0,i<4^order,i++,
		idx=gapsAPTIdx1D2ND[i,4,order];
		out1D=Append[out1D,Extract[out,idx+1]];
	];
	Flatten[out1D]
];


gapsAPTTensorA[A4_,order_]:=Which[
		order==-1,
			gapsAPTA2EB[A4],
		order>=1,	
			gapsAPTGradA[A4,order],
		order==0 || order<-1,
			Print["Wrong order!"];
	];


(* ::Subsection:: *)
(*Generation functions*)


gapsAPTCForm2C[exp_]:=Module[
{cform},
	cform=ToString[CForm[Simplify[exp]]];
	cform=StringReplace[cform,"gapsAPTFieldGenC_Private_tt"->"tt"];
	cform=StringReplace[cform,"gapsAPTFieldGenC_Private_xx"->"xx"];
	cform=StringReplace[cform,"gapsAPTFieldGenC_Private_yy"->"yy"];
	cform=StringReplace[cform,"gapsAPTFieldGenC_Private_zz"->"zz"];
	cform=StringReplace[cform,"Power"->"pow"];
	cform=StringReplace[cform,"Sin"->"sin"];
	cform=StringReplace[cform,"Cos"->"cos"];
	cform=StringReplace[cform,"Exp"->"exp"];
	cform=StringReplace[cform,"Log"->"log"];
	cform=StringReplace[cform,"Sqrt"->"sqrt"];
	cform
]


gapsAPTMForm2M[exp_]:=Module[
{mform},
	mform=ToMatlab[Simplify[exp],1000];
	mform=StringReplace[mform,"gapsAPTFieldGenC`Private`tt"->"tt"];
	mform=StringReplace[mform,"gapsAPTFieldGenC`Private`xx"->"xx"];
	mform=StringReplace[mform,"gapsAPTFieldGenC`Private`yy"->"yy"];
	mform=StringReplace[mform,"gapsAPTFieldGenC`Private`zz"->"zz"];
	mform
]


gapsAPTGenCFileEB[file_,A_]:=Module[
{EB,i,EBstr},
	EB=gapsAPTTensorA[A,-1];
	(*Write part of E*)
	WriteString[file,"\tif(-1 == Order)\n\t{\n\t\tif(pInputs->EMField_Cal_E)\n\t\t{\n"];
	Do[
		EBstr=gapsAPTCForm2C[EB[[i]]];
		WriteString[file,"\t\t\tpTensor["<>ToString[i-1]<>"] = "<>EBstr<>";\n"],
		{i,1,3}
	];
	WriteString[file,"\t\t}\n\n"];
	(*Write part of B*)
	WriteString[file,"\t\tif(pInputs->EMField_Cal_B)\n\t\t{\n"];
	Do[
		EBstr=gapsAPTCForm2C[EB[[i]]];
		WriteString[file,"\t\t\tpTensor["<>ToString[i-1]<>"] = "<>EBstr<>";\n"],
		{i,4,6}
	];
	WriteString[file,"\t\t}\n\t}\n"];
]


gapsAPTGenMFileEB[file_,A_]:=Module[
{EB,i,EBstr},
	EB=gapsAPTTensorA[A,-1];
	(*Write part of E*)
	WriteString[file,"\tif -1 == Order\n\t\tTensor=zeros(6,1);\n\t\tif Inputs.EMField_Cal_E\n"];
	Do[
		EBstr=gapsAPTMForm2M[EB[[i]]];
		WriteString[file,"\t\t\tTensor("<>ToString[i]<>") = "<>EBstr],
		{i,1,3}
	];
	WriteString[file,"\t\tend\n"];
	(*Write part of B*)
	WriteString[file,"\t\tif Inputs.EMField_Cal_B\n"];
	Do[
		EBstr=gapsAPTMForm2M[EB[[i]]];
		WriteString[file,"\t\t\tTensor("<>ToString[i]<>") = "<>EBstr],
		{i,4,6}
	];
	WriteString[file,"\t\tend\n\tend\n"];
]


gapsAPTGenCFileAllATensor[file_,A_,maxorder_]:=Module[
{tensor,iorder,num,i,tensorStr},
	Do[
		tensor=gapsAPTTensorA[A,iorder];
		num=4^iorder;
		WriteString[file,"\tif("<>ToString[iorder]<>" == Order)\n\t{\n"];
		Do[
			tensorStr=gapsAPTCForm2C[tensor[[i]]];
			WriteString[file,"\t\tpTensor["<>ToString[i-1]<>"] = "<>tensorStr<>";\n"];
		,{i,1,num}];
		WriteString[file,"\t}\n"];
	,{iorder,1,maxorder}];
]


gapsAPTGenMFileAllATensor[file_,A_,maxorder_]:=Module[
{tensor,iorder,num,i,tensorStr},
	Do[
		tensor=gapsAPTTensorA[A,iorder];
		num=4^iorder;
		WriteString[file,"\tif "<>ToString[iorder]<>" == Order\n"];
		WriteString[file,"\t\tTensor=zeros("<>ToString[num]<>",1);\n"];
		Do[
			tensorStr=gapsAPTMForm2M[tensor[[i]]];
			WriteString[file,"\t\tTensor("<>ToString[i]<>") = "<>tensorStr];
		,{i,1,num}];
		WriteString[file,"\tend\n"];
	,{iorder,1,maxorder}];
]


gapsAPTGenCFile[A4_,maxorder_]:=Module[
{A,filename,funcname,file},
	A=gapsAPTTransCor[gapsAPTReplaceST[A4],gapsCorType];
	filename="./"<>ToString[gapsName]<>".c";
	funcname="GAPS_APT_Field_"<>ToString[gapsName];
	file=OpenWrite[filename];
	WriteString[file,"#include \"APT_AllHeaders.h\"\n\n"];
	WriteString[file,"int "<>funcname<>"(double *pTensor,double *pSpaceTime4,int Order,Gaps_IO_InputsContainer *pInputs)\n"];
	WriteString[file,"{\n"];
	WriteString[file,"\tint MaxOrder = "<>ToString[maxorder]<>";\n"];
	WriteString[file,"\tdouble tt=pSpaceTime4[0],xx=pSpaceTime4[1],yy=pSpaceTime4[2],zz=pSpaceTime4[3];\n"];
	WriteString[file,#]&/@Map["\tdouble "<>ToString[#]<>" = pInputs->EMField_"<>ToString[gapsName]<>"_"<>ToString[#]<>";\n"&,gapsParameters];
	WriteString[file,"\tdouble B0=pInputs->EMField_B0;\n"];
	WriteString[file,"\tdouble E0=pInputs->EMField_E0;\n"];
	Which[
		maxorder==0||maxorder<-1,
			Print["Wrong order!"],
		maxorder==-1,
			gapsAPTGenCFileEB[file,A],
		maxorder>=1,
			Module[{},
				gapsAPTGenCFileEB[file,A];
				gapsAPTGenCFileAllATensor[file,A,maxorder];
			]
	];
	WriteString[file,"\tif(MaxOrder<Order)\n\t{\n"];
	WriteString[file,"\t\tfprintf(stderr,\"ERROR: In function "<>funcname<>". This field function does NOT support tensors order larger than %d.\\n\",MaxOrder);\n"];
	WriteString[file,"\t}\n"];
	WriteString[file,"\treturn 0;\n"];
	WriteString[file,"}"];
	Close[file];
];


gapsAPTGenMFile[A4_,maxorder_]:=Module[
{A,filename,funcname,file},
	A=gapsAPTTransCor[gapsAPTReplaceST[A4],gapsCorType];
	funcname="GAPS_APT_Field_"<>ToString[gapsName];
	filename="./"<>funcname<>".m";
	file=OpenWrite[filename];
	WriteString[file,"function Tensor="<>funcname<>"(SpaceTime4,Order,Inputs)\n"];
	WriteString[file,"\tMaxOrder = "<>ToString[maxorder]<>";\n"];
	WriteString[file,"\ttt=SpaceTime4(1,:);xx=SpaceTime4(2,:);yy=SpaceTime4(3,:);zz=SpaceTime4(4,:);\n"];
	WriteString[file,#]&/@Map["\t"<>ToString[#]<>" = Inputs.EMField_"<>ToString[gapsName]<>"_"<>ToString[#]<>";\n"&,gapsParameters];
	WriteString[file,"\tB0=Inputs.EMField_B0;\n"];
	WriteString[file,"\tE0=Inputs.EMField_E0;\n"];
	Which[
		maxorder==0||maxorder<-1,
			Print["Wrong order!"],
		maxorder==-1,
			gapsAPTGenMFileEB[file,A],
		maxorder>=1,
			Module[{},
				gapsAPTGenMFileEB[file,A];
				gapsAPTGenMFileAllATensor[file,A,maxorder];
			]
	];
	WriteString[file,"\tif MaxOrder<Order\n"];
	WriteString[file,"\t\terror('ERROR: In function "<>funcname<>". This field function does NOT support tensors order larger than %d.\\n',MaxOrder);\n"];
	WriteString[file,"\tend\n"];
	WriteString[file,"end"];
	Close[file];
];


gapsAPTGenADDFile[maxorder_]:=Module[
{filename,file,paralist,parameters},
	filename="./ADD_Field";
	file=OpenWrite[filename];
	WriteString[file,#]&/@Map["Add_Inputs\t"<>"double\t"<>"1\tEMField_"<>ToString[gapsName]<>"_"<>ToString[#]<>"\t\"Add description:\"\n"&,gapsParameters];
	parameters="";
	paralist=Map["EMField_"<>ToString[gapsName]<>"_"<>ToString[#]&,gapsParameters];
	For[i=0,i<Length[paralist],i++,parameters=parameters<>paralist[[i+1]]<>";"];
	parameters="\""<>parameters<>"\"";
	WriteString[file,"Add_EMField    "<>ToString[gapsName]<>"    "<>parameters<>"    \"MaxOrder:"<>ToString[maxorder]<>"\"    \"Add description: \"\n"];
	Close[file];
];


gapsAPTGenAllFiles[A4_,maxorder_]:=Module[
{},
	gapsAPTGenCFile[A4,maxorder];
	gapsAPTGenMFile[A4,maxorder];
	gapsAPTGenADDFile[maxorder];
]


(* ::Subsection:: *)
(*End package*)


End[]
EndPackage[]
