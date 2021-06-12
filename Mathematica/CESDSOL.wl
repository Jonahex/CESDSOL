(* ::Package:: *)

BeginPackage["CESDSOL`"]; 


Begin["Private`"]; 


ToRegionIndex[Left[coord_], coords_] :=
    2 (Position[coords, coord][[1, 1]] - 1) + 1; 

ToRegionIndex[Right[coord_], coords_] :=
    2 (Position[coords, coord][[1, 1]] - 1) + 2; 


TimeDerivativeField[field_, time_, order_]:=Symbol[ToString[field] <> "$" <> ToString[time] <> ToString[order]];


ToCpp[expr_] :=
    StringReplace[ToString[CForm[expr]], {WordBoundary ~~ "Power" ~~ 
        WordBoundary -> "pow", WordBoundary ~~ "Sin" ~~ WordBoundary -> "sin",
         WordBoundary ~~ "Cos" ~~ WordBoundary -> "cos", WordBoundary ~~ "Tan"
         ~~ WordBoundary -> "tan", WordBoundary ~~ "Sec" ~~ WordBoundary -> "sec",
         WordBoundary ~~ "Csc" ~~ WordBoundary -> "csc", WordBoundary ~~ "Cot"
         ~~ WordBoundary -> "cot", WordBoundary ~~ "\[Pi]" ~~ WordBoundary -> "Pi", WordBoundary ~~ "Sqrt" ~~ WordBoundary -> "sqrt"
        }]; 


WordPattern[word_] :=
    WordBoundary ~~ ToString[word] ~~ WordBoundary; 


Quoted[s_] :=
    "\"" <> s <> "\""; 


PrintGridArray[{points__}] :=
    ToString[N[{points}]]; 

PrintGridArray[CESDSOL`UniformRange[min_, max_, n_, opts:OptionsPattern[
    ]]] :=
    "{MakeUniformRange<double>(" <> ToString[min] <> ", " <> ToString[
        max] <> ", " <> ToString[n] <> "), " <> ToString[OptionValue[UniformRange,
         opts, CESDSOL`Period]] <> "}"; 


Options[CESDSOL`UniformRange] = {CESDSOL`Period -> "std::nullopt"}; 


PrintGrid[CESDSOL`DirectProductGrid[grids__]] :=
    Module[{},
        "auto grid = std::make_shared<DirectProductGrid<" <> ToString[
            Length[{grids}]] <> ", double>>(" <> StringRiffle["SingleDimensionalGrid<double>("
             <> PrintGridArray[#] <> ")"& /@ {grids}, ", "] <> ");\n"
    ]; 

PrintGrid[grid_String] :=
    grid; 


PrintDiscretization[dim_, CESDSOL`StructuredFiniteDifferenceDiscretization[
    stencilSize_]] :=
    "std::unique_ptr<Discretization<"<>ToString[dim]<>">> discretization = std::make_unique<StructuredFiniteDifferenceDiscretization<" <>
         ToString[dim] <> ">>(" <> ToString[stencilSize] <> ")"; 


End[]; 


StationaryProblem[ceqsRaw_, fields_, coords_, grid_, opts:OptionsPattern[
    ]] :=
    Module[{dim, ceqsCount, fieldsCount, params, paramsCount, deqsRaw,
         deqsCount, vars, varsCount, regionCount, needsJacobian, fieldsTransformed,
         derivativeOperators, dOpsCount, rules, result, fieldDerivatives, problemName,
         localPIEs, localPIEExprs, localPIECount, globalPIEs, globalPIEExprs,
         globalPIECount, localVIEs, localVIEExprs, localVIECount, globalVIEs,
         globalVIEExprs, globalVIECount, localVDEs, localVDEExprs, localVDECount,
         globalVDEs, globalVDEExprs, globalVDECount, jacobianDerivators, lvdeJacobian,
         lvdeRules, ceqs, deqs, gvdeJacobian, gvdeRules, fieldDerivativesUnflatten,
         integrals, integralExprs, integralCount, integralJacobian, integralRules,
         gridName},
        dim = Length[coords];
        ceqsCount = Length[ceqsRaw];
        fieldsCount = Length[fields];
        params = OptionValue[Parameters];
        paramsCount = Length[params];
        deqsRaw = OptionValue[DiscreteEquations];
        deqsCount = Length[deqs];
        vars = OptionValue[Variables];
        varsCount = Length[vars];
        regionCount = 2 dim + 1;
        needsJacobian = OptionValue[CalculateJacobian];
        problemName = OptionValue[Name];
        localPIEs = #[[1]]& /@ OptionValue[LocalPIEs];
        localPIEExprs = #[[2]]& /@ OptionValue[LocalPIEs];
        localPIECount = Length[localPIEs];
        globalPIEs = #[[1]]& /@ OptionValue[GlobalPIEs];
        globalPIEExprs = #[[2]]& /@ OptionValue[GlobalPIEs];
        globalPIECount = Length[globalPIEs];
        localVIEs = #[[1]]& /@ OptionValue[LocalVIEs];
        localVIEExprs = #[[2]]& /@ OptionValue[LocalVIEs];
        localVIECount = Length[localVIEs];
        globalVIEs = #[[1]]& /@ OptionValue[GlobalVIEs];
        globalVIEExprs = #[[2]]& /@ OptionValue[GlobalVIEs];
        globalVIECount = Length[globalVIEs];
        localVDEs = #[[1]]& /@ OptionValue[LocalVDEs];
        localVDEExprs = #[[2]]& /@ OptionValue[LocalVDEs];
        localVDECount = Length[localVDEs];
        globalVDEs = #[[1]]& /@ OptionValue[GlobalVDEs];
        globalVDEExprs = #[[2]]& /@ OptionValue[GlobalVDEs];
        globalVDECount = Length[globalVDEs];
        integrals = #[[1]]& /@ OptionValue[Integrals];
        integralExprs = #[[2]]& /@ OptionValue[Integrals];
        integralCount = Length[integrals];
        fieldsTransformed = (# @@ coords&) /@ ceqs;
        derivativeOperators = DeleteDuplicates[Cases[{ceqsRaw}, (Derivative[
            ords__][#] @@ coords -> {ords}), Infinity]]& /@ fields;
        dOpsCount = Length[derivativeOperators];
        fieldDerivativesUnflatten = Thread /@ MapThread[List, {fields,
             derivativeOperators}];
        fieldDerivatives = Flatten[fieldDerivativesUnflatten, 1];
        jacobianDerivators = Join[MapIndexed[(#1 @@ coords -> {#2[[1]
            ] - 1, 0})&, fields], MapIndexed[(#1 -> {#2[[1]] + fieldsCount - 1, 0
            })&, vars], Flatten[MapIndexed[(((Derivative @@ #1[[2]])[#1[[1]]]) @@
             coords -> {#2[[1]] - 1, #2[[2]]})&, fieldDerivativesUnflatten, {2}],
             1]];
        lvdeJacobian = DeleteCases[Flatten[Outer[({#1[[1]], #2[[1]]} 
            -> D[#1[[2]], #2[[1]]])&, OptionValue[LocalVDEs], jacobianDerivators,
             1], 1], _ -> 0];
        integralJacobian = DeleteCases[Flatten[Outer[({#1[[1]], #2[[1
            ]]} -> D[#1[[2]], #2[[1]]])&, OptionValue[Integrals], jacobianDerivators,
             1], 1], _ -> 0];
        gvdeJacobian = DeleteCases[Flatten[Outer[({#1[[1]], #2} -> D[
            #1[[2]], #2])&, OptionValue[GlobalVDEs], vars], 1], _ -> 0];
        lvdeRules = (# @@ coords -> # @@ Cases[lvdeJacobian, Rule[{#,
             der_}, _] -> der])& /@ localVDEs;
        gvdeRules = (# -> # @@ Cases[gvdeJacobian, Rule[{#, der_}, _]
             -> der])& /@ globalVDEs;
        integralRules = (# -> # @@ Cases[integralJacobian, Rule[{#, der_
            }, _] -> der])& /@ integrals;
        {ceqs, deqs} = {ceqsRaw, deqsRaw} /. Join[lvdeRules, gvdeRules,
             integralRules];
        rules =
            Join[
                (Private`ToCpp[#[[2]]] -> "l.VDEValues[" <> ToString[
                    Position[localVDEs, Head[#[[1]]]][[1, 1]] - 1] <> "]")& /@ lvdeRules
                ,
                (Private`ToCpp[#[[2]]] -> "g.VDEValues[" <> ToString[
                    Position[globalVDEs, #[[1]]][[1, 1]] - 1] <> "]")& /@ gvdeRules
                ,
                (Private`ToCpp[#[[2]]] -> "g.ReductionValues[" <> ToString[
                    Position[integrals, #[[1]]][[1, 1]] - 1] <> "]")& /@ integralRules
                ,
                (Private`ToCpp[D[(#[[1, 1]]) @@ coords /. lvdeRules, 
                    #[[1, 2]]]] -> "l.LVDEJacobianComponentValues[" <> ToString[Position[
                    localVDEs, #[[1, 1]]][[1, 1]] - 1] <> ", " <> ToString[(#[[1, 2]] /. 
                    jacobianDerivators)[[1]]] <> ", " <> ToString[(#[[1, 2]] /. jacobianDerivators
                    )[[2]]] <> "]")& /@ lvdeJacobian
                ,
                (Private`ToCpp[D[(#[[1, 1]]) /. integralRules, #[[1, 
                    2]]]] -> "l.ReductionJacobianComponentValues[" <> ToString[Position[integrals,
                     #[[1, 1]]][[1, 1]] - 1] <> ", " <> ToString[(#[[1, 2]] /. jacobianDerivators
                    )[[1]]] <> ", " <> ToString[(#[[1, 2]] /. jacobianDerivators)[[2]]] <>
                     "]")& /@ integralJacobian
                ,
                (Private`ToCpp[D[#[[1, 1]] /. gvdeRules, #[[1, 2]]]] 
                    -> "g.GVDEJacobianComponentValues[" <> ToString[Position[globalVDEs, 
                    #[[1, 1]]][[1, 1]] - 1] <> ", " <> ToString[ceqsCount + Position[vars,
                     #[[1, 2]]][[1, 1]] - 1] <> "]")& /@ gvdeJacobian
                ,
                MapIndexed[(Private`WordPattern[#1] -> "l.Point[" <> 
                    ToString[#2[[1]] - 1] <> "]")&, coords]
                ,
                MapIndexed[(Private`ToCpp[#1 @@ coords] -> "l.FieldValues["
                     <> ToString[#2[[1]] - 1] <> "]")&, fields]
                ,
                MapIndexed[(Private`WordPattern[#1] -> "g.DiscreteVariables["
                     <> ToString[#2[[1]] - 1] <> "]")&, vars]
                ,
                MapIndexed[(Private`WordPattern[#1] -> "g.Parameters["
                     <> ToString[#2[[1]] - 1] <> "]")&, params]
                ,
                MapIndexed[(Private`ToCpp[#1 @@ coords] -> "l.PIEValues["
                     <> ToString[#2[[1]] - 1] <> "]")&, localPIEs]
                ,
                MapIndexed[(Private`WordPattern[#1] -> "g.PIEValues["
                     <> ToString[#2[[1]] - 1] <> "]")&, globalPIEs]
                ,
                MapIndexed[(Private`ToCpp[#1 @@ coords] -> "l.VIEValues["
                     <> ToString[#2[[1]] - 1] <> "]")&, localVIEs]
                ,
                MapIndexed[(Private`WordPattern[#1] -> "g.VIEValues["
                     <> ToString[#2[[1]] - 1] <> "]")&, globalVIEs]
                ,
                MapIndexed[(Private`ToCpp[#1 @@ coords] -> "l.VDEValues["
                     <> ToString[#2[[1]] - 1] <> "]")&, localVDEs]
                ,
                MapIndexed[(Private`WordPattern[#1] -> "g.VDEValues["
                     <> ToString[#2[[1]] - 1] <> "]")&, globalVDEs]
                ,
                MapIndexed[(Private`WordPattern[#1] -> "g.ReductionValues["
                     <> ToString[#2[[1]] - 1] <> "]")&, integrals]
                ,
                Module[{fieldIndex, derivativeIndex},
                    fieldIndex = Position[fields, #[[1]]][[1, 1]];
                    derivativeIndex = Position[derivativeOperators[[fieldIndex
                        ]], #[[2]]][[1, 1]];
                    ToString[CForm[((Derivative @@ #[[2]]) @@ {#[[1]]
                        }) @@ coords]] -> "l.DerivativeValues[" <> ToString[fieldIndex - 1] <>
                         "][" <> ToString[derivativeIndex - 1] <> "]"
                ]& /@ fieldDerivatives
            ];
        result = "auto descriptor = StationaryProblemDescriptor<" <> 
            ToString[dim] <> ", double, double>(GridDescriptor<" <> ToString[dim]
             <> ", double>(" <> ToString[regionCount] <> "), Array<Array<std::array<size_t, "
             <> ToString[dim] <> ">>>" <> ToString[derivativeOperators] <> ", " <>
             ToString[ceqsCount] <> ", " <> ToString[paramsCount] <> ", " <> ToString[
            deqsCount] <> ", " <> ToString[localPIECount] <> ", " <> ToString[globalPIECount
            ] <> ", " <> ToString[localVIECount] <> ", " <> ToString[globalVIECount
            ] <> ", " <> ToString[localVDECount] <> ", " <> ToString[globalVDECount
            ] <> ", " <> ToString[integralCount] <> ");\n";
        If[problemName != "",
            result = StringJoin[result, "descriptor.SetProblemName(" 
                <> Private`Quoted[problemName] <> ");\n"]
        ];
        result = StringJoin[result, MapIndexed["descriptor.SetParameterName("
             <> ToString[#2[[1]] - 1] <> ", " <> Private`Quoted[ToString[#1]] <> 
            ");\n"&, params], MapIndexed["descriptor.SetVariableName(" <> ToString[
            #2[[1]] - 1] <> ", " <> Private`Quoted[ToString[#1]] <> ");\n"&, Join[
            fields, vars]]];
        result = StringJoin[result, MapIndexed["descriptor.SetLocalParameterIndependentExpression("
             <> ToString[#2[[1]] - 1] <> ", [](const auto& l, const auto& g) {return "
             <> StringReplace[Private`ToCpp[#1], rules] <> ";});\n"&, localPIEExprs
            ], MapIndexed["descriptor.SetGlobalParameterIndependentExpression(" <>
             ToString[#2[[1]] - 1] <> ", [](const auto& g) {return " <> StringReplace[
            Private`ToCpp[#1], rules] <> ";});\n"&, globalPIEExprs], MapIndexed["descriptor.SetLocalVariableIndependentExpression("
             <> ToString[#2[[1]] - 1] <> ", [](const auto& l, const auto& g) {return "
             <> StringReplace[Private`ToCpp[#1], rules] <> ";});\n"&, localVIEExprs
            ], MapIndexed["descriptor.SetGlobalVariableIndependentExpression(" <>
             ToString[#2[[1]] - 1] <> ", [](const auto& g) {return " <> StringReplace[
            Private`ToCpp[#1], rules] <> ";});\n"&, globalVIEExprs], MapIndexed["descriptor.SetLocalVariableDependentExpression("
             <> ToString[#2[[1]] - 1] <> ", [](const auto& l, const auto& g) {return "
             <> StringReplace[Private`ToCpp[#1], rules] <> ";});\n"&, localVDEExprs
            ], MapIndexed["descriptor.SetGlobalVariableDependentExpression(" <> ToString[
            #2[[1]] - 1] <> ", [](const auto& g) {return " <> StringReplace[Private`ToCpp[
            #1], rules] <> ";});\n"&, globalVDEExprs], MapIndexed["descriptor.SetIntegrand("
             <> ToString[#2[[1]] - 1] <> ", [](const auto& g) {return " <> StringReplace[
            Private`ToCpp[#1], rules] <> ";});\n"&, integralExprs]];
        result =
            StringJoin[
                result
                ,
                MapIndexed[
                    Module[{eq, region},
                        If[Head[#1] === Rule,
                            eq = #1[[2]];
                            region = Private`ToRegionIndex[#1[[1]], coords
                                ]
                            ,
                            eq = #1;
                            region = 0
                        ];
                        "descriptor.SetContinuousEquation(" <> ToString[
                            #2[[1]] - 1] <> ", " <> ToString[region] <> ", [](const auto& l, const auto& g) {return "
                             <> StringReplace[Private`ToCpp[eq], rules] <> ";});\n"
                    ]&
                    ,
                    ceqs
                    ,
                    {2}
                ]
            ];
        result = StringJoin[result, MapIndexed["descriptor.SetDiscreteEquation("
             <> ToString[#2[[1]] - 1] <> ", [](const auto& globals) {return " <> 
            StringReplace[Private`ToCpp[#1], rules] <> ";});\n"&, deqs]];
        If[needsJacobian,
            result = StringJoin[result, "descriptor.SetLocalVariableDependentExpressionJacobianComponent("
                 <> ToString[Position[localVDEs, #[[1, 1]]][[1, 1]] - 1] <> ", " <> ToString[
                (#[[1, 2]] /. jacobianDerivators)[[1]]] <> ", " <> ToString[(#[[1, 2]
                ] /. jacobianDerivators)[[2]]] <> ", [](const auto& l, const auto& g) {return "
                 <> StringReplace[Private`ToCpp[#[[2]]], rules] <> ";});\n"& /@ lvdeJacobian
                ];
            result = StringJoin[result, "descriptor.SetIntegrandJacobianComponent("
                 <> ToString[Position[integrals, #[[1, 1]]][[1, 1]] - 1] <> ", " <> ToString[
                (#[[1, 2]] /. jacobianDerivators)[[1]]] <> ", " <> ToString[(#[[1, 2]
                ] /. jacobianDerivators)[[2]]] <> ", [](const auto& l, const auto& g) {return "
                 <> StringReplace[Private`ToCpp[#[[2]]], rules] <> ";});\n"& /@ integralJacobian
                ];
            result = StringJoin[result, "descriptor.SetGlobalVariableDependentExpressionJacobianComponent("
                 <> ToString[Position[globalVDEs, #[[1, 1]]][[1, 1]] - 1] <> ", " <> 
                ToString[ceqsCount + Position[vars, #[[1, 2]]][[1, 1]] - 1] <> ", [](const auto& l, const auto& g) {return "
                 <> StringReplace[Private`ToCpp[#[[2]]], rules] <> ";});\n"& /@ gvdeJacobian
                ];
            result =
                StringJoin[
                    result
                    ,
                    MapIndexed[
                        Module[{eq, region, fieldDers, derDers, varDers,
                             ceqIndex},
                            If[Head[#1] === Rule,
                                eq = #1[[2]];
                                region = Private`ToRegionIndex[#1[[1]
                                    ], coords]
                                ,
                                eq = #1;
                                region = 0
                            ];
                            ceqIndex = #2[[1]] - 1;
                            fieldDers = Simplify[D[eq, # @@ coords], 
                                Trig -> False]& /@ fields;
                            varDers = Simplify[D[eq, #], Trig -> False
                                ]& /@ vars;
                            derDers = Simplify[D[eq, ((Derivative @@ 
                                #[[2]]) @@ {#[[1]]}) @@ coords], Trig -> False]& /@ fieldDerivatives;
                                
                            StringJoin[
                                MapIndexed[
                                    If[!SameQ[#1, 0],
                                        "descriptor.SetJacobianComponent("
                                             <> ToString[ceqIndex] <> ", " <> ToString[#2[[1]] - 1] <> ", 0, " <>
                                             ToString[region] <> ", [](const auto& l, const auto& g) {return " <>
                                             StringReplace[Private`ToCpp[#1], rules] <> ";});\n"
                                        ,
                                        ""
                                    ]&
                                    ,
                                    fieldDers
                                ]
                                ,
                                MapIndexed[
                                    If[!SameQ[#1, 0],
                                        "descriptor.SetJacobianComponent("
                                             <> ToString[ceqIndex] <> ", " <> ToString[ceqsCount + #2[[1]] - 1] <>
                                             ", 0, 0, [](const auto& l, const auto& g) {return " <> StringReplace[
                                            Private`ToCpp[#1], rules] <> ";});\n"
                                        ,
                                        ""
                                    ]&
                                    ,
                                    varDers
                                ]
                                ,
                                MapIndexed[
                                    Module[{fieldIndex, derIndex},
                                        If[!SameQ[#1, 0],
                                            fieldIndex = Position[fields,
                                                 fieldDerivatives[[#2[[1]], 1]]][[1, 1]];
                                            derIndex = Position[derivativeOperators[[
                                                fieldIndex]], fieldDerivatives[[#2[[1]], 2]]][[1, 1]];
                                            "descriptor.SetJacobianComponent("
                                                 <> ToString[ceqIndex] <> ", " <> ToString[fieldIndex - 1] <> ", " <>
                                                 ToString[derIndex] <> ", " <> ToString[region] <> ", [](const auto& l, const auto& g) {return "
                                                 <> StringReplace[Private`ToCpp[#1], rules] <> ";});\n"
                                            ,
                                            ""
                                        ]
                                    ]&
                                    ,
                                    derDers
                                ]
                            ]
                        ]&
                        ,
                        ceqs
                        ,
                        {2}
                    ]
                ];
            result =
                StringJoin[
                    result
                    ,
                    MapIndexed[
                        Module[{eq, region, fieldDers, derDers, varDers,
                             deqIndex},
                            If[Head[#1] === Rule,
                                eq = #1[[2]];
                                region = Private`ToRegionIndex[#1[[1]
                                    ], coords]
                                ,
                                eq = #1;
                                region = 0
                            ];
                            deqIndex = ceqsCount + #2[[1]] - 1;
                            varDers = Simplify[D[eq, #]]& /@ vars;
                            StringJoin[
                                MapIndexed[
                                    If[!SameQ[#1, 0],
                                        "descriptor.SetJacobianComponent("
                                             <> ToString[deqIndex] <> ", " <> ToString[ceqsCount + #2[[1]] - 1] <>
                                             ", 0, 0, [](const auto& l, const auto& g) {return " <> StringReplace[
                                            Private`ToCpp[#1], rules] <> ";});\n"
                                        ,
                                        ""
                                    ]&
                                    ,
                                    varDers
                                ]
                            ]
                        ]&
                        ,
                        deqs
                    ]
                ];
            
        ];
        If[StringQ[grid],
            gridName = grid
            ,
            gridName = "grid";
            result = StringJoin[result, Private`PrintGrid[grid]]
        ];
        result = StringJoin[result, Private`PrintDiscretization[dim, 
            OptionValue[Discretization]], ";\n"];
        result = StringJoin[result, "auto problem = descriptor.MakeProblem("
             <> gridName <> ", std::move(discretization));\n"];
        result = StringJoin[result, "problem->AddLocalOutputExpression([](const auto& l, const auto& g){ return "
             <> StringReplace[Private`ToCpp[#[[2]]], rules] <> ";}, " <> Private`Quoted[
            #[[1]]] <> ");\n"& /@ OptionValue[LocalOutput]];
        result = StringJoin[result, "problem->AddGlobalOutputExpression([](const auto& g){ return "
             <> StringReplace[Private`ToCpp[#[[2]]], rules] <> ";}, " <> Private`Quoted[
            #[[1]]] <> ");\n"& /@ OptionValue[GlobalOutput]];
        result = StringJoin[result, "problem->AddIntegralOutputExpression([](const auto& l, const auto& g){ return "
             <> StringReplace[Private`ToCpp[#[[2]]], rules] <> ";}, " <> Private`Quoted[#[[
            1]]] <> ");\n"& /@ OptionValue[IntegralOutput]];
        result = StringJoin[result, "problem->AddPointOutputExpression("
             <> ToString[N[#[[2, 1]]]] <> ", [](const auto& l, const auto& g){ return "
             <> StringReplace[Private`ToCpp[#[[2, 2]]], rules] <> ";}, " <> Private`Quoted[
            #[[1]]] <> ");\n"& /@ OptionValue[PointOutput]];
        result
    ]; 

Options[StationaryProblem] = {CalculateJacobian -> True, Parameters ->
     {}, DiscreteEquations -> {}, Variables -> {}, Name -> "", LocalPIEs 
    -> {}, GlobalPIEs -> {}, LocalVIEs -> {}, GlobalVIEs -> {}, LocalVDEs
     -> {}, GlobalVDEs -> {}, Integrals -> {}, Discretization -> StructuredFiniteDifferenceDiscretization[
    5], LocalOutput -> {}, GlobalOutput -> {}, IntegralOutput -> {}, PointOutput
     -> {}}; 


ExplicitTransientProblem[ceqsRaw_, time_, spatialCoords_, grid_, opts
    :OptionsPattern[]] :=
    Module[{coords, dim, ceqsCount, fieldsCount, params, paramsCount,
         deqsRaw, deqsCount, vars, varsCount, regionCount, derivativeOperators,
         rules, result, fieldDerivatives, problemName, localPIEs, localPIEExprs,
         localPIECount, globalPIEs, globalPIEExprs, globalPIECount, localVIEs,
         localVIEExprs, localVIECount, globalVIEs, globalVIEExprs, globalVIECount,
         localVDEs, localVDEExprs, localVDECount, globalVDEs, globalVDEExprs,
         globalVDECount, ceqs, deqs, fieldDerivativesUnflatten, integrals, integralExprs,
         integralCount, gridName, fieldTopDerivativeOrders, fields, rawVars, 
        varTopDerivativeOrders, rawFields, rawRules},
        coords = Join[spatialCoords, {time}];
        dim = Length[spatialCoords];
        fieldTopDerivativeOrders = #[[1, 0, 0, dim + 1]]& /@ ceqsRaw;
        rawFields = #[[1, 0, 1]]& /@ ceqsRaw;
        fields = Flatten[MapIndexed[Join[{rawFields[[#2]]}, Table[Private`TimeDerivativeField[rawFields[[First[#2]]], time, i],{i, 1, #1-1}]]&,fieldTopDerivativeOrders]];
        fieldsCount = Length[fields];
        deqsRaw = OptionValue[DiscreteEquations];
        rawVars = #[[1, 0, 1]]& /@ deqsRaw;
        varTopDerivativeOrders = #[[1, 0, 0, 1]]& /@ deqsRaw;
        vars = Flatten[MapIndexed[Join[{rawVars[[#2]]}, Table[Private`TimeDerivativeField[rawVars[[First[#2]]], time, i],{i, 2, #1}]]&,varTopDerivativeOrders]];
        varsCount = Length[vars];
        rawRules = {Derivative[Repeated[0, {dim}], n_][f_] :> Private`TimeDerivativeField[f, time, n], Derivative[
                ord:Repeated[x_ /; x > 0, {dim}]
                ,
                n_;
                n > 0
            ][f_] :> Derivative[ord, 0][Private`TimeDerivativeField[f, time, n]], Derivative[$n_][$f_] :> Private`TimeDerivativeField[f, time, n]};
        ceqs = Flatten[MapIndexed[Join[Table[{Private`TimeDerivativeField[rawFields[[First[#2]]], time, i]@@coords},{i, fieldTopDerivativeOrders[[First[#2]]]-1}],{#1[[2]]}]&,ceqsRaw],1]/. rawRules;
        deqs = Flatten[MapIndexed[Join[Table[{Private`TimeDerivativeField[rawVars[[First[#2]]], time, i]@@coords},{i, varTopDerivativeOrders[[First[#2]]]-1}],{#1[[2]]}]&,deqsRaw],1]/. rawRules;
        ceqsCount = Length[ceqs];
        deqsCount = Length[deqs];
        params = OptionValue[Parameters];
        paramsCount = Length[params];
        regionCount = 2 dim + 1;
        problemName = OptionValue[Name];
        localPIEs = #[[1]]& /@ OptionValue[LocalPIEs];
        localPIEExprs = #[[2]]& /@ OptionValue[LocalPIEs];
        localPIECount = Length[localPIEs];
        globalPIEs = #[[1]]& /@ OptionValue[GlobalPIEs];
        globalPIEExprs = #[[2]]& /@ OptionValue[GlobalPIEs];
        globalPIECount = Length[globalPIEs];
        localVIEs = #[[1]]& /@ OptionValue[LocalVIEs];
        localVIEExprs = #[[2]]& /@ OptionValue[LocalVIEs];
        localVIECount = Length[localVIEs];
        globalVIEs = #[[1]]& /@ OptionValue[GlobalVIEs];
        globalVIEExprs = #[[2]]& /@ OptionValue[GlobalVIEs];
        globalVIECount = Length[globalVIEs];
        localVDEs = #[[1]]& /@ OptionValue[LocalVDEs];
        localVDEExprs = #[[2]]& /@ OptionValue[LocalVDEs];
        localVDECount = Length[localVDEs];
        globalVDEs = #[[1]]& /@ OptionValue[GlobalVDEs];
        globalVDEExprs = #[[2]]& /@ OptionValue[GlobalVDEs];
        globalVDECount = Length[globalVDEs];
        integrals = #[[1]]& /@ OptionValue[Integrals];
        integralExprs = #[[2]]& /@ OptionValue[Integrals];
        integralCount = Length[integrals];
        derivativeOperators = DeleteDuplicates[Cases[ceqs, (Derivative[
            ords__][#] @@ coords -> {ords}), Infinity]]& /@ fields;
        fieldDerivativesUnflatten = Thread /@ MapThread[List, {fields,
             derivativeOperators}];
        fieldDerivatives = Flatten[fieldDerivativesUnflatten, 1];
        rules =
            Join[
                MapIndexed[(Private`WordPattern[#1] -> "l.Point[" <> 
                    ToString[#2[[1]] - 1] <> "]")&, coords]
                ,
                MapIndexed[(Private`ToCpp[#1 @@ coords] -> "l.FieldValues["
                     <> ToString[#2[[1]] - 1] <> "]")&, fields]
                ,
                MapIndexed[(Private`ToCpp[#1[time]] -> "g.DiscreteVariables["
                     <> ToString[#2[[1]] - 1] <> "]")&, vars]
                ,
                MapIndexed[(Private`WordPattern[#1] -> "g.Parameters["
                     <> ToString[#2[[1]] - 1] <> "]")&, params]
                ,
                MapIndexed[(Private`ToCpp[#1 @@ coords] -> "l.PIEValues["
                     <> ToString[#2[[1]] - 1] <> "]")&, localPIEs]
                ,
                MapIndexed[(Private`WordPattern[#1] -> "g.PIEValues["
                     <> ToString[#2[[1]] - 1] <> "]")&, globalPIEs]
                ,
                MapIndexed[(Private`ToCpp[#1 @@ coords] -> "l.VIEValues["
                     <> ToString[#2[[1]] - 1] <> "]")&, localVIEs]
                ,
                MapIndexed[(Private`WordPattern[#1] -> "g.VIEValues["
                     <> ToString[#2[[1]] - 1] <> "]")&, globalVIEs]
                ,
                MapIndexed[(Private`ToCpp[#1 @@ coords] -> "l.VDEValues["
                     <> ToString[#2[[1]] - 1] <> "]")&, localVDEs]
                ,
                MapIndexed[(Private`WordPattern[#1] -> "g.VDEValues["
                     <> ToString[#2[[1]] - 1] <> "]")&, globalVDEs]
                ,
                MapIndexed[(Private`WordPattern[#1] -> "g.ReductionValues["
                     <> ToString[#2[[1]] - 1] <> "]")&, integrals]
                ,
                Module[{fieldIndex, derivativeIndex},
                    fieldIndex = Position[fields, #[[1]]][[1, 1]];
                    derivativeIndex = Position[derivativeOperators[[fieldIndex
                        ]], #[[2]]][[1, 1]];
                    ToString[CForm[((Derivative @@ #[[2]]) @@ {#[[1]]
                        }) @@ coords]] -> "l.DerivativeValues[" <> ToString[fieldIndex - 1] <>
                         "][" <> ToString[derivativeIndex - 1] <> "]"
                ]& /@ fieldDerivatives
            ];
        result = "auto descriptor = ExplicitTransientProblemDescriptor<"
             <> ToString[dim] <> ", double, double>(GridDescriptor<" <> ToString[
            dim] <> ", double>(" <> ToString[regionCount] <> "), Array<Array<std::array<size_t, "
             <> ToString[dim] <> ">>>" <> ToString[Map[Drop[#,-1]&,derivativeOperators,{2}]] <> ", " <>
             ToString[ceqsCount] <> ", " <> ToString[paramsCount] <> ", " <> ToString[
            deqsCount] <> ", " <> ToString[localPIECount] <> ", " <> ToString[globalPIECount
            ] <> ", " <> ToString[localVIECount] <> ", " <> ToString[globalVIECount
            ] <> ", " <> ToString[localVDECount] <> ", " <> ToString[globalVDECount
            ] <> ", " <> ToString[integralCount] <> ");\n";
        If[problemName != "",
            result = StringJoin[result, "descriptor.SetProblemName(" 
                <> Private`Quoted[problemName] <> ");\n"]
        ];
        result = StringJoin[result, MapIndexed["descriptor.SetParameterName("
             <> ToString[#2[[1]] - 1] <> ", " <> Private`Quoted[ToString[#1]] <> 
            ");\n"&, params], MapIndexed["descriptor.SetVariableName(" <> ToString[
            #2[[1]] - 1] <> ", " <> Private`Quoted[ToString[#1]] <> ");\n"&, Join[
            fields, vars]]];
        result = StringJoin[result, MapIndexed["descriptor.SetLocalParameterIndependentExpression("
             <> ToString[#2[[1]] - 1] <> ", [](const auto& l, const auto& g) {return "
             <> StringReplace[Private`ToCpp[#1], rules] <> ";});\n"&, localPIEExprs
            ], MapIndexed["descriptor.SetGlobalParameterIndependentExpression(" <>
             ToString[#2[[1]] - 1] <> ", [](const auto& g) {return " <> StringReplace[
            Private`ToCpp[#1], rules] <> ";});\n"&, globalPIEExprs], MapIndexed["descriptor.SetLocalVariableIndependentExpression("
             <> ToString[#2[[1]] - 1] <> ", [](const auto& l, const auto& g) {return "
             <> StringReplace[Private`ToCpp[#1], rules] <> ";});\n"&, localVIEExprs
            ], MapIndexed["descriptor.SetGlobalVariableIndependentExpression(" <>
             ToString[#2[[1]] - 1] <> ", [](const auto& g) {return " <> StringReplace[
            Private`ToCpp[#1], rules] <> ";});\n"&, globalVIEExprs], MapIndexed["descriptor.SetLocalVariableDependentExpression("
             <> ToString[#2[[1]] - 1] <> ", [](const auto& l, const auto& g) {return "
             <> StringReplace[Private`ToCpp[#1], rules] <> ";});\n"&, localVDEExprs
            ], MapIndexed["descriptor.SetGlobalVariableDependentExpression(" <> ToString[
            #2[[1]] - 1] <> ", [](const auto& g) {return " <> StringReplace[Private`ToCpp[
            #1], rules] <> ";});\n"&, globalVDEExprs], MapIndexed["descriptor.SetIntegrand("
             <> ToString[#2[[1]] - 1] <> ", [](const auto& g) {return " <> StringReplace[
            Private`ToCpp[#1], rules] <> ";});\n"&, integralExprs]];
        result =
            StringJoin[
                result
                ,
                MapIndexed[
                    Module[{eq, region},
                        If[Head[#1] === Rule,
                            eq = #1[[2]];
                            region = Private`ToRegionIndex[#1[[1]], coords
                                ]
                            ,
                            eq = #1;
                            region = 0
                        ];
                        "descriptor.SetContinuousEquation(" <> ToString[
                            #2[[1]] - 1] <> ", " <> ToString[region] <> ", [](const auto& l, const auto& g) {return "
                             <> StringReplace[Private`ToCpp[eq], rules] <> ";});\n"
                    ]&
                    ,
                    ceqs
                    ,
                    {2}
                ]
            ];
        result = StringJoin[result, MapIndexed["descriptor.SetDiscreteEquation("
             <> ToString[#2[[1]] - 1] <> ", [](const auto& globals) {return " <> 
            StringReplace[Private`ToCpp[#1], rules] <> ";});\n"&, deqs]];
        If[StringQ[grid],
            gridName = grid
            ,
            gridName = "grid";
            result = StringJoin[result, Private`PrintGrid[grid]]
        ];
        result = StringJoin[result, Private`PrintDiscretization[dim, 
            OptionValue[Discretization]], ";\n"];
        result = StringJoin[result, "auto problem = descriptor.MakeProblem("
             <> gridName <> ", std::move(discretization));\n"];
        result = StringJoin[result, "problem->AddLocalOutputExpression([](const auto& l, const auto& g){ return "
             <> StringReplace[Private`ToCpp[#[[2]]], rules] <> ";}, " <> Private`Quoted[
            #[[1]]] <> ");\n"& /@ OptionValue[LocalOutput]];
        result = StringJoin[result, "problem->AddGlobalOutputExpression([](const auto& g){ return "
             <> StringReplace[Private`ToCpp[#[[2]]], rules] <> ";}, " <> Private`Quoted[
            #[[1]]] <> ");\n"& /@ OptionValue[GlobalOutput]];
        result = StringJoin[result, "problem->AddIntegralOutputExpression([](const auto& l, const auto& g){ return "
             <> StringReplace[ToCpp[#[[2]]], rules] <> ";}, " <> Private`Quoted[#[[
            1]]] <> ");\n"& /@ OptionValue[IntegralOutput]];
        result = StringJoin[result, "problem->AddPointOutputExpression("
             <> ToString[N[#[[2, 1]]]] <> ", [](const auto& l, const auto& g){ return "
             <> StringReplace[Private`ToCpp[#[[2, 2]]], rules] <> ";}, " <> Private`Quoted[
            #[[1]]] <> ");\n"& /@ OptionValue[PointOutput]];
        result
    ]; 

Options[ExplicitTransientProblem] = {Parameters -> {}, DiscreteEquations
     -> {}, Name -> "", LocalPIEs -> {}, GlobalPIEs -> {}, LocalVIEs -> {
    }, GlobalVIEs -> {}, LocalVDEs -> {}, GlobalVDEs -> {}, Integrals -> 
    {}, Discretization -> StructuredFiniteDifferenceDiscretization[5], LocalOutput
     -> {}, GlobalOutput -> {}, IntegralOutput -> {}, PointOutput -> {}}; 


ReadSolution::wrongfile = "Can't read solution data from file `1`"; 

ReadSolution::wrongversion = "Unsupported solution file version `1`"; 

ReadSolution::wronggridtype = "Unsupported grid type `1`"; 

ReadSolution::wrongdimension = "Impossible dimension count `1`"; 

ReadSolution::wronggriddatatype = "Unsupported grid data type `1`"; 

ReadSolution::wrongproblemdatatype = "Unsupported problem data type `1`"; 

ReadSolution::wrongproblemtype = "Unsupported problem type `1`"; 

Begin["Private`"]; 

DataTypes = {"Integer8", "Integer16", "Integer32", "Integer64", "UnsignedInteger8",
     "UnsignedInteger16", "UnsignedInteger32", "UnsignedInteger64", "Real16",
     "Real32", "Real64"}; 

ReadDirectProductGrid[stream_, dimension_, dataType_] :=
    Module[{dimensionSizes, subgrids},
        dimensionSizes = BinaryReadList[stream, "UnsignedInteger64", 
            dimension];
        BinaryReadList[stream, DataTypes[[dataType + 1]], dimension];
            
        subgrids = BinaryReadList[stream, DataTypes[[dataType + 1]], 
            #]& /@ dimensionSizes;
        Flatten[Outer[List, Sequence @@ subgrids], dimension - 1]
    ]; 

ReadStationaryProblemData[stream_, grid_] :=
    Module[{dataType, fieldCount, varCount, paramCount, localOutputCount,
         globalOutputCount, reductionCount, pointOutputCount, parameters, dataTypeString,
         fields, vars, localOutput, globalOutput, reductions, pointOutput, gridSize
        },
        dataType = BinaryRead[stream, "UnsignedInteger64"];
        dataTypeString = DataTypes[[dataType + 1]];
        If[dataType >= Length[DataTypes],
            Message[ReadSolution::wrongproblemdatatype, dataType];
            Return[$Failed]
        ];
        gridSize = Length[grid];
        {fieldCount, varCount, paramCount, localOutputCount, globalOutputCount,
             reductionCount, pointOutputCount} = BinaryReadList[stream, "UnsignedInteger64",
             7];
        parameters = BinaryReadList[stream, dataTypeString, paramCount
            ];
        fields = Partition[BinaryReadList[stream, dataTypeString, fieldCount
             gridSize], gridSize];
        vars = BinaryReadList[stream, dataTypeString, varCount];
        localOutput = Partition[BinaryReadList[stream, dataTypeString,
             localOutputCount gridSize], gridSize];
        globalOutput = BinaryReadList[stream, dataTypeString, globalOutputCount
            ];
        reductions = BinaryReadList[stream, dataTypeString, reductionCount
            ];
        pointOutput = BinaryReadList[stream, dataTypeString, pointOutputCount
            ];
        fields = (ReadString[stream, FromCharacterCode[0]] -> Interpolation[
            Join[grid, List @ #\[Transpose], 2]])& /@ fields;
        vars = (ReadString[stream, FromCharacterCode[0]] -> #)& /@ vars
            ;
        parameters = (ReadString[stream, FromCharacterCode[0]] -> #)&
             /@ parameters;
        localOutput = (ReadString[stream, FromCharacterCode[0]] -> Interpolation[
            Join[grid, List @ #\[Transpose], 2]])& /@ localOutput;
        globalOutput = (ReadString[stream, FromCharacterCode[0]] -> #
            )& /@ globalOutput;
        reductions = (ReadString[stream, FromCharacterCode[0]] -> #)&
             /@ reductions;
        pointOutput = (ReadString[stream, FromCharacterCode[0]] -> #)
            & /@ pointOutput;
        <|"Parameters" -> parameters, "Fields" -> fields, "Variables"
             -> vars, "LocalOutput" -> localOutput, "GlobalOutput" -> globalOutput,
             "Reductions" -> reductions, "PointOutput" -> pointOutput|>
    ]; 

ReadTransientProblemData[stream_, grid_] :=
    Module[{dataType, fieldCount, varCount, paramCount, localOutputCount,
         globalOutputCount, reductionCount, pointOutputCount, parameters, dataTypeString,
         fields, vars, localOutput, globalOutput, reductions, pointOutput, timeCount,
         timeData, gridSize, time, currentGrid},
        dataType = BinaryRead[stream, "UnsignedInteger64"];
        dataTypeString = DataTypes[[dataType + 1]];
        If[dataType >= Length[DataTypes],
            Message[ReadSolution::wrongproblemdatatype, dataType];
            Return[$Failed]
        ];
        gridSize = Length[grid];
        {fieldCount, varCount, paramCount, localOutputCount, globalOutputCount,
             reductionCount, pointOutputCount, timeCount} = BinaryReadList[stream,
             "UnsignedInteger64", 8];
        parameters = BinaryReadList[stream, dataTypeString, paramCount
            ];
        fields = ConstantArray[{}, fieldCount];
        vars = ConstantArray[{}, varCount];
        localOutput = ConstantArray[{}, localOutputCount];
        globalOutput = ConstantArray[{}, globalOutputCount];
        reductions = ConstantArray[{}, reductionCount];
        pointOutput = ConstantArray[{}, pointOutputCount];
        Do[
            time = BinaryRead[stream, dataTypeString];
            currentGrid = Join[{time}, #]& /@ grid;
            fields = Join[fields, Join[currentGrid, List @ #\[Transpose], 2]& /@
                 Partition[BinaryReadList[stream, dataTypeString, fieldCount gridSize],gridSize], 2];
            If[varCount > 0, vars = Join[vars, {time, #}& /@ BinaryReadList[stream,
                 dataTypeString, varCount]]];
            If[localOutputCount > 0, localOutput = Join[localOutput, Join[currentGrid, List @ 
                #\[Transpose], 2]& /@ Partition[BinaryReadList[stream, dataTypeString, localOutputCount gridSize],gridSize
                ], 2]];
            If[globalOutputCount > 0, globalOutput = Join[globalOutput, {time, #}& /@ BinaryReadList[
                stream, dataTypeString, globalOutputCount]]];
            If[reductionCount > 0, reductions = Join[reductions, {time, #}& /@ BinaryReadList[
                stream, dataTypeString, reductionCount]]];
            If[pointOutputCount > 0, pointOutput = Join[pointOutput, {time, #}& /@ BinaryReadList[
                stream, dataTypeString, pointOutputCount]]];
            
            ,
            timeCount
        ];
        fields = (ReadString[stream, FromCharacterCode[0]] -> Interpolation[
            #])& /@ fields;
        vars = (ReadString[stream, FromCharacterCode[0]] -> Interpolation[
            #])& /@ vars;
        parameters = (ReadString[stream, FromCharacterCode[0]] -> #)& /@ parameters;
        localOutput = (ReadString[stream, FromCharacterCode[0]] -> Interpolation[
            #])& /@ localOutput;
        globalOutput = (ReadString[stream, FromCharacterCode[0]] -> Interpolation[
            #])& /@ globalOutput;
        reductions = (ReadString[stream, FromCharacterCode[0]] -> Interpolation[
            #])& /@ reductions;
        pointOutput = (ReadString[stream, FromCharacterCode[0]] -> Interpolation[
            #])& /@ pointOutput;
        <|"Parameters" -> parameters, "Fields" -> fields, "Variables"
             -> vars, "LocalOutput" -> localOutput, "GlobalOutput" -> globalOutput,
             "Reductions" -> reductions, "PointOutput" -> pointOutput|>
    ]; 

ReadGrid[stream_] :=
    Module[{gridType, dimension, dataType},
        {gridType, dimension, dataType} = BinaryReadList[stream, "UnsignedInteger64",
             3];
        If[dimension == 0,
            Message[ReadSolution::wrongdimension, dimension];
            Return[$Failed]
        ];
        If[dataType >= Length[DataTypes],
            Message[ReadSolution::wronggriddatatype, dataType];
            Return[$Failed]
        ];
        Skip[stream, Byte, 16];
        If[gridType == 1,
            ReadDirectProductGrid[stream, dimension, dataType]
            ,
            Message[ReadSolution::wronggridtype, gridType];
            $Failed
        ]
    ]; 

ReadProblemData[stream_, gridSize_] :=
    Module[{problemType},
        problemType = BinaryRead[stream, "UnsignedInteger32"];
        Skip[stream, Byte, 16];
        Switch[problemType,
        0,
            ReadStationaryProblemData[stream, gridSize],
                1
            ,
            ReadTransientProblemData[stream, gridSize],
                _
            ,
            Message[ReadSolution::wrongproblemtype, problemType];
            $Failed
        ]
    ]; 

End[]; 

ReadSolution[path_] :=
    Module[{stream, magic, version, grid, data},
        stream = OpenRead[path, BinaryFormat -> True];
        If[stream == $Failed,
            Return[$Failed]
        ];
        {magic, version} = BinaryReadList[stream, "UnsignedInteger32",
             2];
        If[magic != 1146307907,
            Message[ReadSolution::wrongfile, path];
            Return[$Failed]
        ];
        If[version != 1,
            Message[ReadSolution::wrongversion, version];
            Return[$Failed]
        ];
        grid = Private`ReadGrid[stream];
        If[grid == $Failed,
            Return[$Failed]
        ];
        data = Private`ReadProblemData[stream, grid];
        Close[stream];
        Append[data, Join[data["Fields"], data["Variables"], data["Parameters"
            ], data["LocalOutput"], data["GlobalOutput"], data["Reductions"], data[
            "PointOutput"]]]
    ]; 


EndPackage[]; 
