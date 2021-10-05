// Parametric CSC semi
// Erin Bachynski, October 2016 
GenieRules.Units.setInputUnit(Angle, "deg");

GenieRules.Tolerances.useTolerantModelling = true; 

//
//
Draft = 19.150000m; // Floater draft 
Diameter1 = 10m; // side column diameter  
Diameter2 = 12.800000m; // side column diameter
Diameter3 = 20m; // heave plate diameter
FB1 = 11m; // center column freeboard 
FB2 = 11m; // side column freeboard 
//
x1 = 0m; 	// center-center x distance for side column
y1 = 0m;	// center-center y distance for side column


x2 = 58.500000m; 	// center-center x distance for side column 
y2 = 67.550000m; 	// center-center x distance for side column 

hhp = 2m;

//
Mesh_Length_P = 0.8m; 
NumElWL = NumberOfElements(48); 
NumElWLC = NumberOfElements(12); 
NumEltop= NumberOfElements(12); 
dz = 2m; 
dx = 2m;
// 
elenWL = MeshDensity(Mesh_Length_P); 
Md_def = MeshDensity(Mesh_Length_P);

//Wet Surfaces
outer_wet = WetSurface();
// 
outer_shell = Set();


//Meshing rules
GenieRules.Meshing.elementType = mp1stOrder;
GenieRules.Meshing.superElementType = 1;
GenieRules.Meshing.autoSimplifyTopology = false;
GenieRules.Meshing.autoSplitPeriodicGeometry = false;
GenieRules.Meshing.preference(mpPreferRectangularMesh, false);
GenieRules.Meshing.preference(mpAllowTriangularElements, false);
GenieRules.Meshing.preference(mpPreferPointMassAsNodeMass, true);
GenieRules.Meshing.preference(mpUseDrillingElements, false);
GenieRules.Meshing.preference(mpUseEccentricHinges, true);
GenieRules.Meshing.eliminateInternalEdges = false;
GenieRules.Meshing.eliminateInternalVertices = true;
GenieRules.Meshing.preference(mpIncludeUnusedProperties, false);
GenieRules.Meshing.preference(mpUseLongLoadcaseNames, false);
GenieRules.Meshing.preference(mpUseLongSetNames, false);
GenieRules.Meshing.preference(mpUseLongPropertyNames, false);
GenieRules.Meshing.preference(mpMeshDensityRounded, true);
GenieRules.Meshing.scantlings = msGross;
GenieRules.Meshing.ignoreEccentricities = false;
GenieRules.Meshing.useCocentricBeams = false;
//GenieRules.Meshing.faceMeshStrategy = SesamQuadMesher;
GenieRules.Meshing.edgeMeshStrategy = LinearDistributionEdge;
GenieRules.Meshing.activate(mpMaxAngle, mpFail, true);
GenieRules.Meshing.setLimit(mpMaxAngle, mpFail, 179 deg);
GenieRules.Meshing.activate(mpMaxAngle, mpSplit, false);
GenieRules.Meshing.setLimit(mpMaxAngle, mpSplit, 165 deg);
GenieRules.Meshing.activate(mpMinAngle, mpFail, false);
GenieRules.Meshing.setLimit(mpMinAngle, mpFail, 1 deg);
GenieRules.Meshing.activate(mpMinAngle, mpSplit, false);
GenieRules.Meshing.setLimit(mpMinAngle, mpSplit, 15 deg);
GenieRules.Meshing.activate(mpMaxRelativeJacobi, mpFail, false);
GenieRules.Meshing.setLimit(mpMaxRelativeJacobi, mpFail, 10);
GenieRules.Meshing.activate(mpMaxRelativeJacobi, mpSplit, false);
GenieRules.Meshing.setLimit(mpMaxRelativeJacobi, mpSplit, 5);
GenieRules.Meshing.activate(mpMinNormalizedJacobi, mpFail, false);
GenieRules.Meshing.setLimit(mpMinNormalizedJacobi, mpFail, 0);
GenieRules.Meshing.activate(mpMinNormalizedJacobi, mpSplit, false);
GenieRules.Meshing.setLimit(mpMinNormalizedJacobi, mpSplit, 0.2);
GenieRules.Meshing.activate(mpMinEdge, false);
GenieRules.Meshing.setLimit(mpMinEdge, 0.1);
GenieRules.Meshing.activate(mpMaxChord, false);
GenieRules.Meshing.setLimit(mpMaxChord, 0.2);
GenieRules.Meshing.activate(mpMaxTwistAngle, mpFail, false);
GenieRules.Meshing.setLimit(mpMaxTwistAngle, mpFail, 30 deg);
GenieRules.Meshing.activate(mpMaxTwistAngle, mpSplit, false);
GenieRules.Meshing.setLimit(mpMaxTwistAngle, mpSplit, 10 deg);

//Tolerances Rules
GenieRules.Tolerances.angleTolerance = 2 deg;
GenieRules.Tolerances.pointTolerance = 0.01 m;
GenieRules.Tolerances.useTolerantModelling = true;

//Beam Creation Rules

//Beam Creation Rules
//***** LOAD MODELLING AND ANALYSIS *****//
LC1 = DummyHydroLoadCase();
LC1.setFemLoadcase(1);
LC1.designCondition(lcOperating);
LC1.wetSurface = outer_wet;
//Analyses
Meshing_panels = Analysis(true);
Meshing_panels.add(MeshActivity());
//******************************************************************* 
// Panels 
//******************************************************************* 
// 
// Define start/end points for guide curves 
// Column 1 
// 
Ps11 = Point(-x1,y1,FB2);
Ps12 = Point(-x1,y1,-Draft+hhp);
Ps13 = Point(-x1,y1,-Draft); 
// 
Ps21 = Point(x2,y2/2,FB2); 
Ps22 = Point(x2,y2/2,-Draft+hhp);  
Ps23 = Point(x2,y2/2,-Draft);
//
// 
// 
Pl2 = CreateShellCircularConeCylinder(Ps11,Diameter2/2, Ps12, Diameter2/2, -180,0);
Pl3 = CreateShellCircularConeCylinder(Ps12,Diameter3/2, Ps13, Diameter3/2, -180,0);
Pl4 = CreateShellCircularConeCylinder(Ps21,Diameter2/2, Ps22, Diameter2/2, 0, 360);
Pl5 = CreateShellCircularConeCylinder(Ps22,Diameter3/2, Ps23, Diameter3/2, 0, 360);
//
//
//
// full column
Curve_topS2 = CreateCircleFromPlaneAndRadius(Ps21,Vector3d(0,0,1),Diameter2/2); 
Curve_WLS2 = Curve_topS2.CopyTranslate(Vector3d(0,0,-FB2));
Curve_lowsS2 = Curve_topS2.CopyTranslate(Vector3d(0,0,-FB2-Draft+hhp));
Curve_lowbS2 = CreateCircleFromPlaneAndRadius(Ps22,Vector3d(0,0,1),Diameter3/2); 
Curve_botS2 = Curve_lowbS2.CopyTranslate(Vector3d(0,0,-hhp)); 
PlS2mid = CoverCurves(Curve_lowsS2,Curve_lowbS2);
PlS2bot = CoverCurves(Curve_botS2);
PlS2top = CoverCurves(Curve_topS2);
PlS2top.flipNormal(); 
//
// half column 
Ps113 = Point(-x1+Diameter2/2,0.0,FB2);
Ps114 = Point(-x1,Diameter2/2,FB2);
Ps115 = Point(-x1+Diameter2/2,0.0,FB2);
Curve_topS1 = GuideArcElliptic(Ps11, Ps113, Ps114, true); 
Curve_topS2 = Curve_topS1.copyRotate(Ps11,Vector3d(0,0,1),90*1); 
Curve_WLS11 = Curve_topS1.copyTranslate(Vector3d(0,0,-FB2));
Curve_WLS12 = Curve_topS2.copyTranslate(Vector3d(0,0,-FB2));
Curve_midS11 = Curve_topS1.copyTranslate(Vector3d(0,0,-FB2-Draft+hhp));
Curve_midS12 = Curve_topS2.copyTranslate(Vector3d(0,0,-FB2-Draft+hhp));
Ps116 = Point(-x1+Diameter3/2,0.0,-Draft+hhp);
Ps117 = Point(-x1,Diameter3/2,-Draft+hhp);
Ps118 = Point(-x1+Diameter3/2,0.0,-Draft+hhp);
Ps119 = Point(-x1,Diameter3/2,-Draft); 
Curve_midbS11 = GuideArcElliptic(Ps12, Ps116, Ps117, true); 
Curve_midbS12 = Curve_midbS11.copyRotate(Ps12,Vector3d(0,0,1),90*1); 
Curve_botS11 = Curve_midbS11.copyTranslate(Vector3d(0,0,-hhp)); 
Curve_botS12 = Curve_midbS12.copyTranslate(Vector3d(0,0,-hhp)); 
Line3 = GuideLine(Ps13,Ps119,1); 
PlS11 = CoverCurves(Curve_midS11,Curve_midbS11);
PlS12 = CoverCurves(Curve_midS12,Curve_midbS12);
PlS13 = CoverCurves(Curve_botS11,Curve_botS12);
PlS14 = CoverCurves(Curve_topS1,Curve_topS2);
PlS14.flipNormal();
FEdge1 = FeatureEdge(Line3);


Pl2.front.wetSurface = outer_wet;
Pl3.front.wetSurface = outer_wet;
Pl4.front.wetSurface = outer_wet;
Pl5.front.wetSurface = outer_wet;
PlS2mid.front.wetSurface = outer_wet;
PlS2bot.front.wetSurface = outer_wet;
PlS2top.front.wetSurface = outer_wet;
PlS11.front.wetSurface = outer_wet;
PlS12.front.wetSurface = outer_wet;
PlS13.front.wetSurface = outer_wet;
PlS14.front.wetSurface = outer_wet;

outer_shell.add(Pl2); 
outer_shell.add(Pl3); 
outer_shell.add(Pl4); 
outer_shell.add(Pl5); 
outer_shell.add(PlS2mid); 
outer_shell.add(PlS2bot); 
outer_shell.add(PlS2top); 
outer_shell.add(PlS11); 
outer_shell.add(PlS12); 
outer_shell.add(PlS13); 
outer_shell.add(PlS14); 

FEdge4 = FeatureEdge(Curve_WLS2);
FEdge5 = FeatureEdge(Curve_WLS11);
FEdge6 = FeatureEdge(Curve_WLS12);

FEdge4.numberOfElements = NumElWL;
FEdge5.numberOfElements = NumElWLC;
FEdge6.numberOfElements = NumElWLC;

LC1.generateAppliedLoads();


// Set mesh density
//Md_def = MeshDensity(Mesh_Length_P);
//Md_def.setDefault();

//  
Md_def.setDefault();


LC1.generateAppliedLoads();  
LC1.setFemLoadcase(1);  

Meshing_panels.execute();
ExportMeshFem().DoExport("T1.FEM");



