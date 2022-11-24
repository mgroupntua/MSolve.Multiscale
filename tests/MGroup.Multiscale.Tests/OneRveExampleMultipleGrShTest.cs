using MGroup.LinearAlgebra.Commons;
using MGroup.LinearAlgebra.Matrices;
using MGroup.LinearAlgebra.Vectors;
//using MGroup.MSolve.Discretization.Commons;
using MGroup.MSolve.MultiscaleAnalysis.SupportiveClasses;
using MGroup.MSolve.Solution.LinearSystem;
using MGroup.Multiscale.SupportiveClasses;
using MGroup.Multiscale.Tests.ExampleModels;
using MGroup.Multiscale.Tests.SeparationBenchmarks2;
using MGroup.Solvers.AlgebraicModel;

using Xunit;

namespace MGroup.Multiscale.Tests
{
	public static class OneRveExampleMultipleGrShTest
	{
		[Fact]
		public static void CheckOneRveSerial()
		{
			(double[] stressesCheck3, double[] stressesCheck4, double[,] consCheck1, double[] uInitialFreeDOFs_state1, double[] uInitialFreeDOFs_state2) = OneRveExample.Check_Graphene_rve_serial();

			string results_file1 = "..\\..\\..\\MGroup.Multiscale.Tests\\InputFiles\\MicroStructureAndIntegrationOneRveMultipleGrSh\\uInitialFreeDOFs_state1.txt";
			string results_file2 = "..\\..\\..\\MGroup.Multiscale.Tests\\InputFiles\\MicroStructureAndIntegrationOneRveMultipleGrSh\\uInitialFreeDOFs_state2.txt";
			//string results_file3 = "..\\..\\..\\InputFiles\\MicroStructureAndIntegrationTestsB\\consCheck1.txt";
			string results_file4 = "..\\..\\..\\MGroup.Multiscale.Tests\\InputFiles\\MicroStructureAndIntegrationOneRveMultipleGrSh\\stressesCheck3.txt";
			string results_file5 = "..\\..\\..\\MGroup.Multiscale.Tests\\InputFiles\\MicroStructureAndIntegrationOneRveMultipleGrSh\\stressesCheck4.txt";
			double[] displacements1sIncrement = PrintUtilities.ReadVector(results_file1);
			double[] displacements2ndIncrement = PrintUtilities.ReadVector(results_file2);

			double[,] consCheck1Expected = new double[6, 6]
			{{8.247794441602281,5.040483382644040,5.045179342838760,-0.034573545680066,0.012873618640199,0.067413461733790},
			{5.040483382644040,7.758675250745090,5.083447516662590,-0.017660393516958,0.086264761000810,-0.001886483315119},
			{5.045179342838760,5.083447516662600,7.889514025249530,0.014993568822868,0.174547712576532,0.013639601528685},
			{-0.034573545680067,-0.017660393516956,0.014993568822868,1.404689076704550,0.023343385610862,0.099337624448147},
			{0.012873618640199,0.086264761000810,0.174547712576533,0.023343385610861,1.347276707954930,-0.002212957880199},
			{0.067413461733791,-0.001886483315119,0.013639601528686,0.099337624448147,-0.002212957880199,1.454060010268960} };

			double[] stressesCheck3Expected = PrintUtilities.ReadVector(results_file4);
			double[] stressesCheck4Expected = PrintUtilities.ReadVector(results_file5);

			Assert.True(AreDisplacementsSame(displacements1sIncrement, uInitialFreeDOFs_state1));
			Assert.True(AreDisplacementsSame(displacements2ndIncrement, uInitialFreeDOFs_state2));
			Assert.True(AreDisplacementsSame(consCheck1, consCheck1Expected));
			Assert.True(AreDisplacementsSame(stressesCheck3, stressesCheck3Expected));
			Assert.True(AreDisplacementsSame(stressesCheck4, stressesCheck4Expected));
		}

		
		public static void CheckOneRveParallel()
		{
			(int[] hexaPrint, int[] cohePrint, int[] shellPrint) = OneRveExample.Check_Graphene_rve_parallel();

			string results_file1 = "..\\..\\..\\MGroup.Multiscale.Tests\\RveTemplates\\Input\\RveGrShMultiple\\rve_no_1\\subdomainHexas.txt";
			string results_file2 = "..\\..\\..\\MGroup.Multiscale.Tests\\RveTemplates\\Input\\RveGrShMultiple\\rve_no_1\\subdomainCohesiveElements.txt";
			string results_file3 = "..\\..\\..\\MGroup.Multiscale.Tests\\RveTemplates\\Input\\RveGrShMultiple\\rve_no_1\\subdomainShellElements.txt";
			//string results_file4 = "..\\..\\..\\InputFiles\\MicroStructureAndIntegrationOneRveMultipleGrSh\\stressesCheck3.txt";
			//string results_file5 = "..\\..\\..\\InputFiles\\MicroStructureAndIntegrationOneRveMultipleGrSh\\stressesCheck4.txt";
			int[] hexaPrintExpected = PrintUtilities.ReadIntVector(results_file1);
			int[] cohePrintExpected = PrintUtilities.ReadIntVector(results_file2);
			int[] shellPrintExpected = PrintUtilities.ReadIntVector(results_file3);


			Assert.True(AreDisplacementsSame(hexaPrint, hexaPrintExpected));
			Assert.True(AreDisplacementsSame(cohePrint, cohePrintExpected));
			Assert.True(AreDisplacementsSame(shellPrint, shellPrintExpected));
		}

		public static bool AreDisplacementsSame(int[] expectedValues,
			int[] computedValues)
		{
			var comparer = new ValueComparer(1E-14);
			for (int i1 = 0; i1 < expectedValues.GetLength(0); i1++)
			{
				if (!comparer.AreEqual(expectedValues[i1], computedValues[i1]))
				{
					return false;
				}
			}
			return true;
		}

		public static bool AreDisplacementsSame(double[] expectedValues,
			double[] computedValues)
		{
			var comparer = new ValueComparer(1E-9);
			for (int i1 = 0; i1 < expectedValues.GetLength(0); i1++)
			{
				if (!comparer.AreEqual(expectedValues[i1], computedValues[i1]))
				{
					return false;
				}
			}
			return true;
		}

		public static bool AreDisplacementsSame(double[,] expectedValues,
			double[,] computedValues)
		{
			var comparer = new ValueComparer(1E-8);
			for (int i1 = 0; i1 < expectedValues.GetLength(0); i1++)
			{
				for (int i2 = 0; i2 < expectedValues.GetLength(1); i2++)
				{
					if (!comparer.AreEqual(expectedValues[i1, i2], computedValues[i1, i2]))
					{
						return false;
					}
				}
			}
			return true;
		}
	}
}
