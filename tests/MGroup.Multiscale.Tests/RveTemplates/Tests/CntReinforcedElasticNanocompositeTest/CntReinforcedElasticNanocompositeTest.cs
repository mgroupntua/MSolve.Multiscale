//namespace MGroup.Multiscale.Tests.RveTemplates.Tests.CntReinforcedElasticNanocompositeTest
//{
//	using System;
//	using System.Collections.Generic;
//	using System.Text;
//	using MGroup.Constitutive.Structural.BoundaryConditions;
//	using MGroup.Constitutive.Structural;
//	using MGroup.NumericalAnalyzers;
//	using MGroup.Solvers.Direct;
//	using MGroup.Stochastic;

//	using Xunit;
//	using System.Linq;
//	using MGroup.Constitutive.Structural.Continuum;
//	using MGroup.LinearAlgebra.Matrices;
//	using MGroup.MSolve.MultiscaleAnalysis;
//	using MiMsolve.SolutionStrategies;

//	public class CntReinforcedElasticNanocompositeTest
//	{
//		[Fact]
//		public static void CheckCntReinforcedElasticNanocompositeExample()
//		{
//			var rve = new CntReinforcedElasticNanocomposite(10);

//			IContinuumMaterial3D microstructure = new Microstructure3D<SkylineMatrix>(rve, false, 1, new SkylineSolverPrefernce());

//			double[,] consCheck1 = new double[6, 6];
//			for (int i1 = 0; i1 < 6; i1++) { for (int i2 = 0; i2 < 6; i2++) { consCheck1[i1, i2] = microstructure.ConstitutiveMatrix[i1, i2]; } }

//			double[] stresses = microstructure.UpdateConstitutiveMatrixAndEvaluateResponse(new double[9] { 1.05, 0, 0, 0, 0, 0, 0, 0, 0 });
//			double[] stressesCheck1 = stresses;
//			microstructure.CreateState();

//			stresses = microstructure.UpdateConstitutiveMatrixAndEvaluateResponse(new double[9] { 1.10, 0, 0, 0, 0, 0, 0, 0, 0 });
//			double[] stressesCheck2 = stresses;
//			var matrix1 = Matrix.CreateZero(3, 3); for (int i1 = 0; i1 < 3; i1++) { for (int i2 = 0; i2 < 3; i2++) { matrix1[i1, i2] = microstructure.ConstitutiveMatrix[i1, i2]; } }

//		}
//	}
//}
