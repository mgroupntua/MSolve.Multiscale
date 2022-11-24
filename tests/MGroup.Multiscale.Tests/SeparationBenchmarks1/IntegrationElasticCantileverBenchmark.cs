namespace MGroup.Multiscale.Tests.SeparationBenchmarks1
{
	using System.Collections.Generic;

	using MGroup.Constitutive.Structural;
	using MGroup.Constitutive.Structural.BoundaryConditions;
	using MGroup.Constitutive.Structural.Continuum;
	using MGroup.FEM.Structural.Continuum;
	using MGroup.LinearAlgebra.Matrices;

	using MGroup.MSolve.Discretization;
	using MGroup.MSolve.Discretization.Dofs;
	using MGroup.MSolve.Discretization.Entities;
	using MGroup.MSolve.MultiscaleAnalysis;
	using MGroup.MSolve.MultiscaleAnalysis.Interfaces;

	using MGroup.MSolve.Numerics.Integration.Quadratures;
	using MGroup.MSolve.Numerics.Interpolation;
	using MGroup.MSolve.Solution;
	//using MGroup.Multiscale.Interfaces;
	using MGroup.NumericalAnalyzers;
	using MGroup.NumericalAnalyzers.Discretization.NonLinear;
	using MGroup.NumericalAnalyzers.Logging;
	//using MGroup.Problems;
	using MGroup.Solvers.Direct;

	using MiMsolve.SolutionStrategies;

	public static class IntegrationElasticCantileverBenchmark 
	{
		
		public static TotalDisplacementsPerIterationLog RunExample()
		{
			
			Model model = new Model();
			int subdomainID = 1; model.SubdomainsDictionary.Add(subdomainID, new Subdomain(subdomainID));
			HexaCantileverBuilder_copyMS_222(model, 0.00219881744271988174427);

			

			// Solver
			var solverBuilder = new SkylineSolver.Factory();
			var algebraicModel = solverBuilder.BuildAlgebraicModel(model);
			ISolver solver = solverBuilder.BuildSolver(algebraicModel);

			// Problem type
			var provider = new ProblemStructural(model, algebraicModel, solver);

			

			var increments = 2;
			var childAnalyzerBuilder = new LoadControlAnalyzer.Builder(algebraicModel, solver, provider, increments);
			childAnalyzerBuilder.MaxIterationsPerIncrement = 100;
			childAnalyzerBuilder.NumIterationsForMatrixRebuild = 1;
			LoadControlAnalyzer childAnalyzer = childAnalyzerBuilder.Build();
			var parentAnalyzer = new StaticAnalyzer(algebraicModel,  provider, childAnalyzer);
			var watchDofs = new Dictionary<int, int[]>();
			watchDofs.Add(subdomainID, new int[5] { 0, 11, 23, 35, 47 });
			var log1 = new TotalDisplacementsPerIterationLog(
				new List<(INode node, IDofType dof)>()
				{
					(model.NodesDictionary[5], StructuralDof.TranslationX),
					(model.NodesDictionary[8], StructuralDof.TranslationZ),
					(model.NodesDictionary[12], StructuralDof.TranslationZ),
					(model.NodesDictionary[16], StructuralDof.TranslationZ),
					(model.NodesDictionary[20], StructuralDof.TranslationZ)
				}, algebraicModel
			);
			childAnalyzer.TotalDisplacementsPerIterationLog = log1;

			parentAnalyzer.Initialize();
			parentAnalyzer.Solve();


			return log1;
		}

		public static void HexaCantileverBuilder_copyMS_222(Model model, double load_value)
		{
			//Origin: ParadeigmataElegxwnBuilder.HexaCantileverBuilder_copyMS_222(Model model, double load_value)

			IRVEbuilder homogeneousRveBuilder1 = new HomogeneousRVEBuilderNonLinear();

			IContinuumMaterial3DDefGrad material1 = new MicrostructureDefGrad3D<SkylineMatrix>(homogeneousRveBuilder1, false, 1, new SkylineSolverPrefernce());

			double[,] nodeData = new double[,] { {-0.250000,-0.250000,-1.000000},
			{0.250000,-0.250000,-1.000000},
			{-0.250000,0.250000,-1.000000},
			{0.250000,0.250000,-1.000000},
			{-0.250000,-0.250000,-0.500000},
			{0.250000,-0.250000,-0.500000},
			{-0.250000,0.250000,-0.500000},
			{0.250000,0.250000,-0.500000},
			{-0.250000,-0.250000,0.000000},
			{0.250000,-0.250000,0.000000},
			{-0.250000,0.250000,0.000000},
			{0.250000,0.250000,0.000000},
			{-0.250000,-0.250000,0.500000},
			{0.250000,-0.250000,0.500000},
			{-0.250000,0.250000,0.500000},
			{0.250000,0.250000,0.500000},
			{-0.250000,-0.250000,1.000000},
			{0.250000,-0.250000,1.000000},
			{-0.250000,0.250000,1.000000},
			{0.250000,0.250000,1.000000}};

			int[,] elementData = new int[,] {{1,8,7,5,6,4,3,1,2},
			{2,12,11,9,10,8,7,5,6},
			{3,16,15,13,14,12,11,9,10},
			{4,20,19,17,18,16,15,13,14}, };

			// orismos shmeiwn
			for (int nNode = 0; nNode < nodeData.GetLength(0); nNode++)
			{
				model.NodesDictionary.Add(nNode + 1, new Node(id: nNode + 1, x: nodeData[nNode, 0], y: nodeData[nNode, 1], z: nodeData[nNode, 2]));

			}

			// orismos elements 
			IElementType e1;
			int subdomainID = 1;
			for (int nElement = 0; nElement < elementData.GetLength(0); nElement++)
			{
				List<Node> elementNodes = new List<Node>();
				for (int j = 0; j < 8; j++)
				{
					elementNodes.Add((Node)model.NodesDictionary[elementData[nElement, j + 1]]);
				}
				e1 = new ContinuumElement3DNonLinearDefGrad(elementNodes, material1, GaussLegendre3D.GetQuadratureWithOrder(2, 2, 2), InterpolationHexa8.UniqueInstance) // dixws to e. exoume sfalma enw sto beambuilding oxi//edw kaleitai me ena orisma to Hexa8
				{
					ID = nElement + 1,
					//ElementType = new ContinuumElement3DNonLinearDefGrad(material1, GaussLegendre3D.GetQuadratureWithOrder(2, 2, 2)) // dixws to e. exoume sfalma enw sto beambuilding oxi//edw kaleitai me ena orisma to Hexa8
				};
				
				model.ElementsDictionary.Add(e1.ID, e1);
				model.SubdomainsDictionary[subdomainID].Elements.Add(e1);
			}

			// constraint vashh opou z=-1
			var boundaryConditions = new List<NodalDisplacement>();
			for (int k = 1; k < 5; k++)
			{
				boundaryConditions.Add(new NodalDisplacement(model.NodesDictionary[k], StructuralDof.TranslationX, 0));
				boundaryConditions.Add(new NodalDisplacement(model.NodesDictionary[k], StructuralDof.TranslationY, 0));
				boundaryConditions.Add(new NodalDisplacement(model.NodesDictionary[k], StructuralDof.TranslationZ, 0));
				
			}

			// fortish korufhs
			//Load load1;
			var nodalLoads = new List<NodalLoad>();
			for (int k = 17; k < 21; k++)
			{
				nodalLoads.Add(new NodalLoad(model.NodesDictionary[k], StructuralDof.TranslationX, 1 * load_value));
				
			}
			var boundaryConditionsAndNodalLoads = new StructuralBoundaryConditionSet(boundaryConditions, nodalLoads);
			model.BoundaryConditions.Add(boundaryConditionsAndNodalLoads);
		}
	}
}
