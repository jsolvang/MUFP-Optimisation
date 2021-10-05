//Exported using: HydroD V4.10-01 started 23-Feb-2021 14:02:41
// ***********************************************************
DirectionSet1 = DirectionSet(Array(0 deg,90 deg));
// ***********************************************************
FrequencySet1 = FrequencySet(FrequencyTypeFrequency, Array(0.01 rad/s,0.1 rad/s,0.2 rad/s, 0.25 rad/s, 0.275 rad/s, 0.3 rad/s, 0.325 rad/s, 0.35 rad/s, 0.4 rad/s, 0.45 rad/s, 0.475 rad/s, 0.5 rad/s, 0.525 rad/s, 0.55 rad/s, 0.575 rad/s, 0.6 rad/s, 0.65 rad/s, 0.7 rad/s,0.8 rad/s,0.9 rad/s,1 rad/s, 1 rad/s,1.025 rad/s,1.05 rad/s,1.075 rad/s,1.1 rad/s,1.125 rad/s,1.15 rad/s,1.175 rad/s,1.2 rad/s,1.225 rad/s,1.25 rad/s,1.275 rad/s,1.3 rad/s,1.325 rad/s,1.35 rad/s,1.375 rad/s,1.4 rad/s,1.425 rad/s,1.45 rad/s,1.475 rad/s,1.5 rad/s,1.525 rad/s,1.55 rad/s,1.575 rad/s,1.6 rad/s,1.625 rad/s,1.65 rad/s,1.675 rad/s,1.7 rad/s,1.725 rad/s,1.75 rad/s,1.775 rad/s,1.8 rad/s,1.825 rad/s,1.85 rad/s,1.875 rad/s,1.9 rad/s,1.925 rad/s,1.95 rad/s,1.975 rad/s,2 rad/s, 2.1 rad/s,2.2 rad/s,2.25 rad/s,2.3 rad/s, 2.35 rad/s, 2.4 rad/s, 2.45 rad/s, 2.5 rad/s, 2.55 rad/s, 2.6 rad/s, 2.65 rad/s, 2.7 rad/s,2.8 rad/s,2.9 rad/s,3 rad/s,3.05 rad/s,3.2 rad/s,3.4 rad/s,3.6 rad/s,3.8 rad/s,4 rad/s,4.5 rad/s,5 rad/s));
// ***********************************************************
Location1 = Location();
Location1.setDepth(300 m);
Location1.gravity = 9.80665 m/s^2;
Location1.air().density = 1.226 Kg/m^3;
Location1.air().kinematicViscosity = 1.462e-005 m^2/s;
Location1.water().density = 1025 Kg/m^3;
Location1.water().kinematicViscosity = 1.19e-006 m^2/s;
Location1.seabed().normaldirection = Vector3d(0 m,0 m,1 m);
// ***********************************************************
Condition1 = FrequencyDomain(Location1);
Condition1.waterSurface().directionSet = DirectionSet1;
Condition1.waterSurface().frequencySet = FrequencySet1;
Condition1.water().setNoCurrent();
// ***********************************************************
HydroModel1 = HydroModel(HydroModelFloating);
HydroModel1.setColumnStabilized(false);
HydroModel1.setBaselineZPos(0 m);
HydroModel1.setAPXPos(0 m);
HydroModel1.setFPXPos(100 m);
HydroModel1.clearReportMetaCenterRotationAxisAzim();
HydroModel1.addReportMetaCenterRotationAxisAzim(0 deg);
HydroModel1.addReportMetaCenterRotationAxisAzim(90 deg);
HydroModel1.clearReportHeelTrimCombinations();
HydroModel1.addReportHeelTrimCombination(0 deg, 0 deg);
HydroModel1.clearReportZWaterlines();
HydroModel1.addReportZWaterline(-15.5 m);
HydroModel1.addReportZWaterline(-15 m);
HydroModel1.addReportZWaterline(-14.5 m);
HydroModel1.addReportZWaterline(-14 m);
HydroModel1.addReportZWaterline(-13.5 m);
HydroModel1.addReportZWaterline(-13 m);
HydroModel1.addReportZWaterline(-12.5 m);
HydroModel1.addReportZWaterline(-12 m);
HydroModel1.addReportZWaterline(-11.5 m);
HydroModel1.addReportZWaterline(-11 m);
HydroModel1.addReportZWaterline(-10.5 m);
HydroModel1.addReportZWaterline(-10 m);
HydroModel1.addReportZWaterline(-9.5 m);
HydroModel1.addReportZWaterline(-9 m);
HydroModel1.addReportZWaterline(-8.5 m);
HydroModel1.addReportZWaterline(-8 m);
HydroModel1.addReportZWaterline(-7.5 m);
HydroModel1.addReportZWaterline(-7 m);
HydroModel1.addReportZWaterline(-6.5 m);
HydroModel1.addReportZWaterline(-6 m);
HydroModel1.addReportZWaterline(-5.5 m);
HydroModel1.addReportZWaterline(-5 m);
HydroModel1.addReportZWaterline(-4.5 m);
HydroModel1.addReportZWaterline(-4 m);
HydroModel1.addReportZWaterline(-3.5 m);
HydroModel1.addReportZWaterline(-3 m);
HydroModel1.addReportZWaterline(-2.5 m);
HydroModel1.addReportZWaterline(-2 m);
HydroModel1.addReportZWaterline(-1.5 m);
HydroModel1.addReportZWaterline(-1 m);
HydroModel1.addReportZWaterline(-0.5 m);
HydroModel1.addReportZWaterline(0 m);
HydroModel1.addReportZWaterline(0.5 m);
HydroModel1.addReportZWaterline(1 m);
HydroModel1.setKidTag("6bb094c6-003c-4928-98c7-3f4a53e6809f");
// ***********************************************************
PanelModel1 = PanelModel(HydroModel1, ElementEnumsFEMFile, "T1.FEM", true, false);
PanelModel1.setModelTranslation(Vector3d(-39.000000 m,0 m,0 m)); 
PanelModel1.regenerateGeometry();
PanelModel1.setKidTag("6319dba9-fb77-4504-901e-8d89d1352ad4");
// ***********************************************************
LoadingCondition1 = LoadingCondition(HydroModel1, 0 m, 0 deg, 0 deg);
LoadingCondition1.setByDraft(false);
LoadingCondition1.interpolateDampingMatrices(false);
LoadingCondition1.addDampingMatricesToWadam(true);
LoadingCondition1.setKidTag("f08fe34e-3ecb-424f-a62a-7eba96ed1091");
// ***********************************************************
MassModel1 = MassModel(LoadingCondition1, MassModelSpecified);
MassModel1.setUserMassCoordinateSystem(GlobalCoordinateSystem);
MassModel1.setTotalMass(8718160.888027 Kg); 
MassModel1.setCOG(Point(0 m,0 m,0.449800 m)); 
MassModel1.setRadiusGyration(Vector3d(41.117803 m,46.495319 m,43.906068 m)); 
MassModel1.setSpecificProductInertia(0.000000 m,-276.443168 m,0.000000 m); 
MassModel1.internalDynamics(false);
MassModel1.updateStiffnessWithFreeSurfaceEffect(true);
MassModel1.setKidTag("3a731dc1-3bee-4723-b31d-4ae6ba9a8cb0");
// ***********************************************************
WadamRun1 = WadamRun();
WadamRun1.useMultiBody(false);
WadamRun1.setHydroModel(HydroModel1);
WadamRun1.setLoadingCondition(LoadingCondition1);
WadamRun1.setEnvironmentData(Condition1);
WadamRun1.useSeaState(false);
WadamRun1.useFreeSurfaceIntegral(true);
WadamRun1.useThreeLayerFreeSurfaceModel(false);
WadamRun1.setFreeSurfaceAnnulusSize(500 m);
WadamRun1.dataCheck(false);
WadamRun1.wadamLite(false);
WadamRun1.setAnalysisType(WadamRunGlobalResponse);
WadamRun1.calculateDrift(false);
WadamRun1.setDriftForceType(WadamRunUnidirectionalWaves);
WadamRun1.driftByFarFieldIntegration(false);
WadamRun1.waveDriftDamping(false);
WadamRun1.setSolverType(WadamRunDirectSolver);
WadamRun1.setMaxMatrixDimension(50000);
WadamRun1.setSingularityType(WadamRunAnalyticalSingularity);
WadamRun1.setIntegrationType(WadamRunOneNodeGauss);
WadamRun1.setPanelDimensionType(WadamRunMaximumDiagonalPanelDimension);
WadamRun1.removeIrrFrequencies(false);
WadamRun1.setIncludeForwardSpeedEffect(false);
WadamRun1.setForwardSpeedX(0 m/s);
WadamRun1.setForwardSpeedY(0 m/s);
WadamRun1.saveTempWamitFiles(false);
WadamRun1.stopBeforePotenExecution(false);
WadamRun1.bypassPotenExecution(false);
WadamRun1.stopBeforeFirstForceExecution(false);
WadamRun1.bypassFirstForceExecution(false);
WadamRun1.useWadamMassCalculation(false);
WadamRun1.stopBeforeSecondForceExecution(false);
WadamRun1.bypassSecondForceExecution(false);
WadamRun1.useSaveRestart(false);
WadamRun1.setPrintType(WadamRunModelDataPrint);
WadamRun1.setResponseFileType(WadamRunSIFFormatted);
WadamRun1.calculateEigenvalues(true);
WadamRun1.defineGlobalResponseReferencePoint(false);
WadamRun1.setGlobalResponseReferencePoint(Point(0 m,0 m,0 m));
WadamRun1.setWaterlinePanelMethod(WadamRunNoPanelPressureAdjustment);
WadamRun1.setLoadTransferLinearInterpolationXmin(0);
WadamRun1.setLoadTransferLinearInterpolationXmax(0);
WadamRun1.sumFrequencyResults(false);
WadamRun1.differenceFrequencyResults(false);
WadamRun1.setToleranceWaterLine(5);
WadamRun1.setToleranceCOG(5);
WadamRun1.setCharacteristicLength(100 m);
WadamRun1.calculateRollDamping(false);
WadamRun1.specifyOutputDirectory(false, false);
WadamRun1.autoOwerwriteExistingResultFiles(false);
WadamRun1.moveResultFiles(false);
WadamRun1.useRunNameAsPrefix(false);
WadamRun1.setKidTag("15091209-2341-43a3-88cb-4ba1b421f615");
// ***********************************************************
WadamWizard1 = WadamWizard();
WadamWizard1.setModelConfiguration(WadamWorkflowPanelModel);
WadamWizard1.setPanelModelConfiguration(WadamWorkflowElementPanelModel);
WadamWizard1.timeDomain(false);
WadamWizard1.rollDamping = true;
WadamWizard1.stochasticRollDamping = false;
WadamWizard1.loadCrossection = true;
WadamWizard1.loadTransfer = true;
WadamWizard1.tankPressure = false;
WadamWizard1.pressurePanels = true;
WadamWizard1.offbodyPoints = true;
WadamWizard1.secondOrder = false;
WadamWizard1.dampingMatrix = false;
WadamWizard1.criticalDamping = false;
WadamWizard1.restoringMatrix = false;
WadamWizard1.addEditObject(DirectionSet1);
WadamWizard1.addEditObject(FrequencySet1);
WadamWizard1.addEditObject(Location1);
WadamWizard1.addEditObject(Condition1);
WadamWizard1.addEditObject(HydroModel1);
WadamWizard1.addEditObject(PanelModel1);
WadamWizard1.addEditObject(LoadingCondition1);
WadamWizard1.addEditObject(MassModel1);
WadamWizard1.addEditObject(WadamRun1);
WadamWizard1.setKidTag("a34fa70f-7639-40f7-b4a3-dc7cce1d4e21");
WadamRun1.execute();