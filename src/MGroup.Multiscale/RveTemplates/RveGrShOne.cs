using System;
using System.Collections.Generic;
using System.Linq;
using MGroup.FEM.Embedding;
using MGroup.FEM.Entities;
using MGroup.FEM.Structural.Embedding;
using MGroup.MSolve.Discretization;
using MGroup.Multiscale.Interfaces;
using MGroup.Multiscale.SupportiveClasses;

namespace MGroup.Multiscale.RveTemplates
{
	/// <summary>
	/// Model builder that can be used to create graphene reinforced rves 
	/// Authors: Gerasimos Sotiropoulos
	/// </summary>
	public class RveGrShOne : IRVEbuilder //IdegenerateRVEbuilder
	{
		//GrapheneReinforcedRVEBuilderExample35fe2boundstiffHostTestPostData_forCheck
		//Origin  einai o updated GrapheneReinforcedRVEBuilderExample35fe2boundstiffHostTestPostData(  to idio einai) apo thn ekdosh ms_development_nl_elements_merge
		// modifications: discretization data kai paths wste
		//anti gia REFERENCE_kanonikh_gewmetria_fe2_post_dg\REF2_10__000_renu_new_multiple_algorithms_check_develop_gia_fe2_3grsh_4182dofs_multiple2
		//na trexoume REFERENCE_kanonikh_gewmetria_fe2_post_dg\REF2_10__000_renu_new_multiple_algorithms_check_develop_gia_fe2_3grsh_4182dofs_multiple2_forCHECK
		// pou exei th gewmetria tou REFERENCE_kanonikh_gewmetria_2\REF2_10__000_renu_new_multiple_algorithms_check_develop_1GrSh_correct_coh_CHECK_integration
		//dld tou "GrapheneReinforcedRVEBuilderCHECK" gia na elegxoume tis nees domes

		Tuple<rveMatrixParameters, grapheneSheetParameters> mpgp;
		rveMatrixParameters mp;
		grapheneSheetParameters gp;
		string renumbering_vector_path;
		int RVE_id;

		public RveGrShOne(int RVE_id)
		{
			this.RVE_id = RVE_id;
		}

		public IRVEbuilder Clone(int a) => new RveGrShOne(a);

		public Tuple<Model, Dictionary<int, INode>, double> GetModelAndBoundaryNodes()
		{
			return Reference2RVEExample10000withRenumberingwithInput_forMS();
		}

		private Tuple<Model, Dictionary<int, INode>, double> Reference2RVEExample10000withRenumberingwithInput_forMS()
		{
			Model model = new Model();
			model.SubdomainsDictionary.Add(1, new Subdomain(1));

			Dictionary<int, INode> boundaryNodes = new Dictionary<int, INode>();

			//Origin  public static void Reference2RVEExample10000withRenumberingwithInput(Model model)
			double[,] Dq;
			//Tuple<rveMatrixParameters, grapheneSheetParameters> mpgp;
			//rveMatrixParameters mp;
			//grapheneSheetParameters gp;
			var rve_id_data = RVE_id.ToString();

			renumbering_vector_path = "..\\..\\..\\MGroup.Multiscale.Tests\\RveTemplates\\Input\\RveGrShOne\\rve_no_{0}\\REF_new_total_numbering.txt";
			renumbering_vector_path = string.Format(renumbering_vector_path, rve_id_data);

			string Fxk_p_komvoi_rve_path = "..\\..\\..\\MGroup.Multiscale.Tests\\RveTemplates\\Input\\RveGrShOne\\rve_no_{0}\\Fxk_p_komvoi_rve.txt";
			Fxk_p_komvoi_rve_path = string.Format(Fxk_p_komvoi_rve_path, rve_id_data);

			string o_xsunol_input_path_gen = "..\\..\\..\\MGroup.Multiscale.Tests\\RveTemplates\\Input\\RveGrShOne\\rve_no_{0}\\o_xsunol_gs_";
			o_xsunol_input_path_gen = string.Format(o_xsunol_input_path_gen, rve_id_data);
			o_xsunol_input_path_gen = o_xsunol_input_path_gen + "{0}.txt";
			int subdiscr1 = 6;
			int discr1 = 1;
			// int discr2 dn xrhsimopoieitai
			int discr3 = 6;
			int subdiscr1_shell = 3;
			int discr1_shell = 1;
			//mpgp = FEMMeshBuilder.GetReferenceKanonikhGewmetriaRveExampleParametersStiffCase(subdiscr1, discr1, discr3, subdiscr1_shell, discr1_shell);
			mpgp = FEMMeshBuilder.GetReferenceKanonikhGewmetriaRveExampleParameters(subdiscr1, discr1, discr3, subdiscr1_shell, discr1_shell);
			mp = mpgp.Item1; //mp.hexa1 = 9; mp.hexa2 = 9; mp.hexa3 = 9;
			gp = mpgp.Item2;


			int graphene_sheets_number = 1;
			o_x_parameters[] model_o_x_parameteroi = new o_x_parameters[graphene_sheets_number];
			double[][] ekk_xyz = new double[graphene_sheets_number][];



			Dq = new double[9, 3 * (((mp.hexa1 + 1) * (mp.hexa2 + 1) * (mp.hexa3 + 1)) - ((mp.hexa1 - 1) * (mp.hexa2 - 1) * (mp.hexa3 - 1)))];
			FEMMeshBuilder.HexaElementsOnlyRVEwithRenumbering_forMS(model, mp, Dq, renumbering_vector_path, boundaryNodes);
			double volume = mp.L01 * mp.L02 * mp.L03;

			int hexaElementsNumber = model.ElementsDictionary.Count();

			//IEnumerable<Element> hostGroup = model.ElementsDictionary.Where(x => (x.Key < hexaElementsNumber + 1)).Select(kv => kv.Value);
			List<int> EmbeddedElementsIDs = new List<int>();
			int element_counter_after_Adding_sheet;
			element_counter_after_Adding_sheet = hexaElementsNumber; // initial value before adding first graphene sheet
			int shellElementsNumber;

			for (int j = 0; j < graphene_sheets_number; j++)
			{
				string file_no = (j + 1).ToString();
				string ox_sunol_input_path = string.Format(o_xsunol_input_path_gen, file_no);
				FEMMeshBuilder.AddGrapheneSheet_with_o_x_Input_withRenumbering(model, gp, ekk_xyz[j], model_o_x_parameteroi[j], renumbering_vector_path, ox_sunol_input_path);
				shellElementsNumber = (model.ElementsDictionary.Count() - element_counter_after_Adding_sheet) / 3; //tha xrhsimefsei
				for (int k = shellElementsNumber + element_counter_after_Adding_sheet + 1; k < model.ElementsDictionary.Count() + 1; k++)
				{
					EmbeddedElementsIDs.Add(model.ElementsDictionary[k].ID);
				}
				element_counter_after_Adding_sheet = model.ElementsDictionary.Count();

			}



			int[] EmbElementsIds = EmbeddedElementsIDs.ToArray();
			IEnumerable<Element> embdeddedGroup = model.ElementsDictionary.Where(x => (Array.IndexOf(EmbElementsIds, x.Key) > -1)).Select(kv => kv.Value); // dld einai null afth th stigmh
																																						   //var embeddedGrouping = new EmbeddedCohesiveGrouping(model, hostGroup, embdeddedGroup);

			//var CohesiveGroupings = new EmbeddedCohesiveGrouping[EmbElementsIds.GetLength(0)];

			var hostSubGroups = new Dictionary<int, IEnumerable<Element>>();
			for (int i1 = 0; i1 < EmbElementsIds.GetLength(0); i1++)
			{
				hostSubGroups.Add(EmbElementsIds[i1], FEMMeshBuilder.GetHostGroupForCohesiveElement(model.ElementsDictionary[EmbElementsIds[i1]], mp, model, renumbering_vector_path));
				//var embeddedGroup_i1 = new List<Element>(1) { model.ElementsDictionary[EmbElementsIds[i1]] };
				//CohesiveGroupings[i1] = new EmbeddedCohesiveGrouping(model, hostGroup_i1, embeddedGroup_i1);
			}

			var CohesiveGroupping = new EmbeddedCohesiveSubGrouping(model, hostSubGroups, embdeddedGroup);

			return new Tuple<Model, Dictionary<int, INode>, double>(model, boundaryNodes, volume);

		}

		// PROSOXH DEN ARKEI MONO TO PARAKATW NA GINEI UNCOMMENT WSTE NA GINEI IMPLEMENT TO IDegenerateRVEBuilder 
		//xreiazetai kai na xrhsimopoithei h katallhlh methodos tou femmeshbuilder gia to model and boundary nodes na dinei mono ta peripheral
		//public Dictionary<Node, IList<DOFType>> GetModelRigidBodyNodeConstraints(Model model)
		//{
		//    return FEMMeshBuilder.GetConstraintsOfDegenerateRVEForNonSingularStiffnessMatrix_withRenumbering(model, mp.hexa1, mp.hexa2, mp.hexa3, renumbering_vector_path);
		//    //TODO:  Pithanws na epistrefetai apo GetModelAndBoundaryNodes ... AndConstraints.
		//}

	}
}
