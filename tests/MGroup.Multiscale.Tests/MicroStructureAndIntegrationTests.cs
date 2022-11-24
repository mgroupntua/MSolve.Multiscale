using MGroup.LinearAlgebra.Commons;
using MGroup.LinearAlgebra.Vectors;
using MGroup.MSolve.MultiscaleAnalysis.SupportiveClasses;
using MGroup.Multiscale.SupportiveClasses;
using MGroup.Multiscale.Tests.FEMpartB;

using Xunit;

namespace MGroup.Multiscale.Tests
{
	public static class MicroStructureAndIntegrationTests
	{
		[Fact]
		public static void CheckElasticMicrostructure()
		{
			(double[] stressesCheck1, double[] stressesCheck2, double[] stressesCheck3, double[] stressesCheck4, double[] uInitialFreeDOFs_state1, double[] uInitialFreeDOFs_state2) = SeparateCodeCheckingClass4.Check05bStressIntegrationObje_Integration();

			string results_file1 = "..\\..\\..\\MGroup.Multiscale.Tests\\InputFiles\\MicroStructureAndIntegrationTestsA\\U_sunol_micro_1.txt";
			string results_file2 = "..\\..\\..\\MGroup.Multiscale.Tests\\InputFiles\\MicroStructureAndIntegrationTestsA\\U_sunol_micro_2.txt";
			double[] displacements1sIncrement = PrintUtilities.ReadVector(results_file1);
			double[] displacements2ndIncrement = PrintUtilities.ReadVector(results_file2);

			Assert.True(ComparisonMetods.AreDisplacementsSame(displacements1sIncrement, uInitialFreeDOFs_state1));
			Assert.True(ComparisonMetods.AreDisplacementsSame(displacements2ndIncrement, uInitialFreeDOFs_state2));

		}

		[Fact]
		public static void CheckGrapheneMicrostructure()
		{
			(double[] stressesCheck3, double[] stressesCheck4, double[,] consCheck1, double[] uInitialFreeDOFs_state1, double[] uInitialFreeDOFs_state2) = SeparateCodeCheckingClass4.Check_Graphene_rve_Obje_Integration();

			string results_file1 = "..\\..\\..\\MGroup.Multiscale.Tests\\InputFiles\\MicroStructureAndIntegrationTestsB\\uInitialFreeDOFs_state1.txt";
			string results_file2 = "..\\..\\..\\MGroup.Multiscale.Tests\\InputFiles\\MicroStructureAndIntegrationTestsB\\uInitialFreeDOFs_state2.txt";
			//string results_file3 = "..\\..\\..\\InputFiles\\MicroStructureAndIntegrationTestsB\\consCheck1.txt";
			string results_file4 = "..\\..\\..\\MGroup.Multiscale.Tests\\InputFiles\\MicroStructureAndIntegrationTestsB\\stressesCheck3.txt";
			string results_file5 = "..\\..\\..\\MGroup.Multiscale.Tests\\InputFiles\\MicroStructureAndIntegrationTestsB\\stressesCheck4.txt";
			double[] displacements1sIncrement = PrintUtilities.ReadVector(results_file1);
			double[] displacements2ndIncrement = PrintUtilities.ReadVector(results_file2);

			double[,] consCheck1Expected = new double[6, 6] {{7.586076977452720,5.002676833389410,5.002022820741120,-0.013645880596011,-0.000827683355376,0.009516079959912},
															 {5.002676833389410,7.512154809921630,5.019436489046710,-0.004684943090660,0.015322738052621,-0.005624423827287},
															 {5.002022820741120,5.019436489046710,7.541834561367190,0.004377190727383,0.028459910927874,0.006755413641553},
															 {-0.013645880596012,-0.004684943090660,0.004377190727383,1.266215688675580,-0.001138059373760,0.018013447177758},
															 {-0.000827683355376,0.015322738052620,0.028459910927874,-0.001138059373760,1.270911173682930,-0.001860391453226},
															 {0.009516079959912,-0.005624423827287,0.006755413641553,0.018013447177758,-0.001860391453226,1.278204450851190},};

			double[] stressesCheck3Expected = PrintUtilities.ReadVector(results_file4);
			double[] stressesCheck4Expected = PrintUtilities.ReadVector(results_file5);

			Assert.True(ComparisonMetods.AreDisplacementsSame(displacements1sIncrement, uInitialFreeDOFs_state1, 1E-9));
			Assert.True(ComparisonMetods.AreDisplacementsSame(displacements2ndIncrement, uInitialFreeDOFs_state2, 1E-9));
			Assert.True(ComparisonMetods.AreDisplacementsSame(consCheck1, consCheck1Expected, 1E-8));
			Assert.True(ComparisonMetods.AreDisplacementsSame(stressesCheck3, stressesCheck3Expected, 1E-12));
			Assert.True(ComparisonMetods.AreDisplacementsSame(stressesCheck4, stressesCheck4Expected, 1E-11));//11

		}
	}
}
