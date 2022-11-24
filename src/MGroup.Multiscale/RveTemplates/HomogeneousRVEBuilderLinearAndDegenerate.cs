using System;
using System.Collections.Generic;

using MGroup.Constitutive.Structural;
//using MGroup.FEM.Entities;
using MGroup.MSolve.Discretization.Dofs;
using MGroup.MSolve.Discretization.Entities;
using MGroup.MSolve.MultiscaleAnalysis.Interfaces;
using MGroup.MSolve.MultiscaleAnalysis.SupportiveClasses;
//using MGroup.MSolve.Discretization.FreedomDegrees;
//using MGroup.MSolve.Discretization.Interfaces;
using MGroup.Multiscale.SupportiveClasses;

namespace MGroup.Multiscale
{
	/// <summary>
	/// Creates an elastic matrix rve for (3d generate to) 2D problems
	/// Authors Gerasimos Sotiropoulos
	/// </summary>
	public class HomogeneousRVEBuilderLinearAndDegenerate : IdegenerateRVEbuilder
	{
		private Tuple<rveMatrixParameters, grapheneSheetParameters> mpgp;
		private rveMatrixParameters mp;
		private grapheneSheetParameters gp;
		private string renumbering_vector_path;

		public HomogeneousRVEBuilderLinearAndDegenerate()
		{
			//TODOGerasimos
			// this.renumbering_vector_path=renumbering_vector_path,
			// this.subdiscr1=subdiscr1
		}
		public IRVEbuilder Clone(int a) => new HomogeneousRVEBuilderLinearAndDegenerate();

		public Tuple<Model, Dictionary<int, INode>, double> GetModelAndBoundaryNodes()
		{
			return Reference2RVEExample10_000withRenumbering_mono_hexa();
		}



		public Tuple<Model, Dictionary<int, INode>, double> Reference2RVEExample10_000withRenumbering_mono_hexa()
		{
			Model model = new Model();
			model.SubdomainsDictionary.Add(1, new Subdomain(1));

			Dictionary<int, INode> boundaryNodes = new Dictionary<int, INode>();
			// COPY APO: Reference2RVEExample100_000withRenumbering_mono_hexa
			double[,] Dq = new double[1, 1];
			//Tuple<rveMatrixParameters, grapheneSheetParameters> mpgp;
			//rveMatrixParameters mp;
			//grapheneSheetParameters gp;
			renumbering_vector_path = "..\\..\\..\\MGroup.Multiscale.Tests\\RveTemplates\\Input\\RveHomogeneous\\REF_new_total_numbering27.txt";
			string Fxk_p_komvoi_rve_path = @"C:\Users\turbo-x\Desktop\notes_elegxoi\REFERENCE_fe2_diafora_check\fe2_tax_me1_arxiko_chol_dixws_me1_OriginalRVEExampleChol_me_a1_REF2_10_000_renu_new_multiple_algorithms_check_stress_27hexa\Fxk_p_komvoi_rve.txt";
			int subdiscr1 = 1;
			int discr1 = 3;
			// int discr2 dn xrhsimopoieitai
			int discr3 = 3;
			int subdiscr1_shell = 7;
			int discr1_shell = 1;

			mpgp = FEMMeshBuilder.GetReferenceRveExampleParameters(subdiscr1, discr1, discr3, subdiscr1_shell, discr1_shell);
			mp = mpgp.Item1;
			gp = mpgp.Item2;
			double[][] ekk_xyz = new double[2][] { new double[] { 0, 0, 0 }, new double[] { 0.25 * 105, 0, 0.25 * 40 } };

			int graphene_sheets_number = 0; // 0 gra sheets afou exoume mono hexa
			o_x_parameters[] model_o_x_parameteroi = new o_x_parameters[graphene_sheets_number];


			FEMMeshBuilder.LinearHexaElementsOnlyRVEwithRenumbering_forMS_PeripheralNodes(model, mp, Dq, renumbering_vector_path, boundaryNodes);
			double volume = mp.L01 * mp.L02 * mp.L03;



			return new Tuple<Model, Dictionary<int, INode>, double>(model, boundaryNodes, volume);
		}

		public Dictionary<Node, IList<IStructuralDofType>> GetModelRigidBodyNodeConstraints(Model model)
		{
			return FEMMeshBuilder.GetConstraintsOfDegenerateRVEForNonSingularStiffnessMatrix_withRenumbering(model, mp.hexa1, mp.hexa2, mp.hexa3, renumbering_vector_path);
			//TODO:  Pithanws na epistrefetai apo GetModelAndBoundaryNodes ... AndConstraints.
		}


	}
	
}
