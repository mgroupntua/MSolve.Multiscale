



// TODO UNCOMMENT THIS TEST WHEN ELEMENT PULL REQUEST IS ACCEPTED


//using System;
//using System.Collections.Generic;
//using System.Diagnostics;
//using System.Globalization;
//using System.IO;
//using System.Linq;
//using System.Threading;

//using MGroup.Constitutive.Structural;
//using MGroup.Constitutive.Structural.BoundaryConditions;
//using MGroup.Constitutive.Structural.Cohesive;
//using MGroup.Constitutive.Structural.Continuum;
//using MGroup.Constitutive.Structural.Line;
//using MGroup.FEM.Structural.Continuum;
//using MGroup.FEM.Structural.Embedding;
//using MGroup.FEM.Structural.Line;
//using MGroup.MSolve.Discretization;
//using MGroup.MSolve.Discretization.Dofs;
//using MGroup.MSolve.Discretization.Entities;
//using MGroup.MSolve.MultiscaleAnalysis.Interfaces;
//using MGroup.MSolve.Numerics.Integration.Quadratures;
//using MGroup.MSolve.Numerics.Interpolation;
//using MGroup.Multiscale.SupportiveClasses;

//namespace MGroup.Stochastic
//{
//	public class CntReinforcedElasticNanocomposite : IRVEbuilder
//	{
//		int hexa1 = 3;
//		int hexa2 = 3;
//		int hexa3 = 3;

//		double L01 = 100;
//		double L02 = 100;
//		double L03 = 100;

//		IIsotropicContinuumMaterial3D matrixMaterial;
//		IIsotropicContinuumMaterial3D inclusionMaterial;
//		int hostElements { get; set; }
//		int embeddedElements { get; set; }
//		int hostNodes { get; set; }
//		int embeddedNodes { get; set; }

//		// cnt paramaters
//		IIsotropicContinuumMaterial3D CntMaterial;
//		int numberOfCnts;
//		int cntLength = 50;
//		// define mechanical properties
//		private double youngModulus = 1.0;//1.051e12; // 5490; // 
//		private double shearModulus = 1.0;//0.45e12; // 871; // 
//		double poissonRatio;  //2.15; // 0.034;
//		double area = 694.77; // 1776.65;  // CNT(20,20)-LinearEBE-TBT-L = 10nm
//		double inertiaY = 100.18; //1058.55;
//		double inertiaZ = 100.18; //1058.55;1058.55;
//		double torsionalInertia = 68.77; //496.38;
//		double effectiveAreaY;
//		double effectiveAreaZ;
//		int subdomainID = 0;

		
//		public double K_el { get; set; }
//		public double K_pl { get; set; }
//		public double T_max { get; set; }

//		private double[] constParameters;

//		public Model model { get; private set; }
//		List<INode> elementNodesClone;
//		List<INode> elementNodesBeam;

//		public bool readFromText { get; set; }

//		public CntReinforcedElasticNanocomposite(int numberOfCnts)
//		{
//			K_el = 10; K_pl = 1.0; T_max = 0.10;
//			this.matrixMaterial = new ElasticMaterial3D(youngModulus: 3.5, poissonRatio: 0.4);
//			constParameters = new double[3] { K_el, K_pl, T_max };
			

//			this.CntMaterial = new ElasticMaterial3D(youngModulus, (double)((youngModulus / (2 * shearModulus)) - 1));

//			double effectiveAreaY = area;
//			double effectiveAreaZ = area;

//			this.numberOfCnts = numberOfCnts;
//		}

//		public Tuple<Model, Dictionary<int, INode>, double> GetModelAndBoundaryNodes()
//		{
//			(int[] NodeIds, double[,] node_coords) = GetHexaRveNodesData();
//			int[,] elementConnectivity = GetHexaRveConnectivity();

//			model = new Model();
//			model.SubdomainsDictionary[0] = new Subdomain(0);
//			AddHexaElements(model, NodeIds, node_coords, elementConnectivity);

		
//			hostNodes = model.NodesDictionary.Count;
//			embeddedNodes = 2 * numberOfCnts;
//			hostElements = model.ElementsDictionary.Count;
//			embeddedElements = numberOfCnts;

//			(int[] cntNodeIds, double[,] cntNodeCoords, int[,] cntElementConnectivity) = GetCntBeamsNodesData(hostNodes, hostElements, readFromText);

//			AddCntBeamElements(model, cntNodeIds, cntNodeCoords, cntElementConnectivity);
//			AddCohesiveBeamElements(model, cntNodeIds, cntNodeCoords, cntElementConnectivity);
//			var embeddedGrouping = EmbeddedBeam3DGrouping.CreateCohesive(model, model.ElementsDictionary
//						.Where(x => x.Key < hostElements).Select(kv => kv.Value).ToArray(), model.ElementsDictionary.Where(x => x.Key >= hostElements + embeddedElements)
//						.Select(kv => kv.Value).ToArray(), true);

			

//			var boundaryNodesIds = GetBoundaryNodeIds();
//			boundaryNodesIds.Sort();
//			Dictionary<int, INode> boundaryNodes = boundaryNodesIds.ToDictionary(t => t, t => model.NodesDictionary[t]);
//			return new Tuple<Model, Dictionary<int, INode>, double>(model, boundaryNodes, L01 * L02 * L03);
//		}

		

//		private List<int> GetBoundaryNodeIds()
//		{
//			var boundaryNodes = new List<int>();
//			for (int k = 0; k < hexa3 + 1; k++)
//			{
//				var indexZeta = (hexa1 + 1) * (hexa2 + 1) * k;
//				if (k == 0 || k == hexa3)
//				{
//					for (int j = 0; j < hexa2 + 1; j++)
//					{
//						for (int i = 0; i < hexa1 + 1; i++)
//						{
//							boundaryNodes.Add(i + (hexa1 + 1) * j + indexZeta);
//						}
//					}
//				}
//				else
//				{
//					for (int i = 0; i < hexa1 + 1; i++)
//					{
//						boundaryNodes.Add(i + indexZeta);
//						boundaryNodes.Add(i + (hexa1 + 1) * hexa2 + indexZeta);
//					}

//					for (int j = 1; j < hexa2; j++)
//					{
//						boundaryNodes.Add(indexZeta + (hexa1 + 1) * j);
//						boundaryNodes.Add(indexZeta + (hexa1 + 1) * j + hexa1);
//					}
//				}
//			}
			
//			return boundaryNodes;
//		}

//		private void AddCntBeamElements(Model model, int[] cntNodeIds, double[,] cntNodeCoordinates, int[,] cntElementConnectivity)
//		{
//			var hostElements = model.ElementsDictionary.Count;
//			for (int i = 0; i < cntNodeIds.Length; i++)
//			{
//				var nodeID = cntNodeIds[i];
//				var coordX = cntNodeCoordinates[i, 0];
//				var coordY = cntNodeCoordinates[i, 1];
//				var coordZ = cntNodeCoordinates[i, 2];

//				model.NodesDictionary.Add(nodeID, new Node(id: nodeID, x: coordX, y: coordY, z: coordZ));
//			}

//			var beamSection =
//				new BeamSection3D(area, inertiaY, inertiaZ, torsionalInertia, effectiveAreaY, effectiveAreaZ);

//			for (int i = 0; i < cntElementConnectivity.GetLength(0); i++)
//			{
//				var elementNodes = new List<INode>
//				{
//					model.NodesDictionary[cntElementConnectivity[i, 0]],
//					model.NodesDictionary[cntElementConnectivity[i, 1]]
//				};
//				var cntMaterialClone = (ElasticMaterial3D)this.CntMaterial;
//				var CntMaterial = new ElasticMaterial3D(cntMaterialClone.YoungModulus, cntMaterialClone.PoissonRatio);
//				IElementType beam_1 = new Beam3DCorotationalQuaternion(elementNodes, CntMaterial.YoungModulus, CntMaterial.PoissonRatio, 7.85, beamSection);
//				beam_1.ID = i + hostElements;
				

//				model.ElementsDictionary.Add(beam_1.ID, beam_1);
//				model.SubdomainsDictionary[0].Elements.Add(beam_1);
//			}

//		}

//		private void AddCohesiveBeamElements(Model model, int[] cntNodeIds, double[,] cntNodeCoordinates, int[,] cntElementConnectivity)
//		{
//			// define mechanical properties

//			double mi = 8.0;
//			double ni = 8.0;
//			double thickness_CNT = 0.34;
//			double a = 0.241;
//			double diameter_CNT = (a / Math.PI) * Math.Sqrt(Math.Pow(ni, 2) + ni * mi + Math.Pow(mi, 2));
//			double radius_CNT = diameter_CNT / 2.0;
//			double radius_CNT_outer = radius_CNT + (thickness_CNT / 2);
//			double CntPerimeter = 2.0 * Math.PI * radius_CNT_outer;

//			for (int i = 0; i < cntNodeIds.Length; i++)
//			{
//				var nodeID = i + hostNodes + embeddedNodes;
//				var coordX = cntNodeCoordinates[i, 0];
//				var coordY = cntNodeCoordinates[i, 1];
//				var coordZ = cntNodeCoordinates[i, 2];
//				model.NodesDictionary.Add(nodeID, new Node(id: nodeID, x: coordX, y: coordY, z: coordZ));
//			}

//			// Create Cohesive Material
//			var cohesiveMaterial = new BondSlipMaterial(K_el, K_pl, 100.0, T_max, new double[2], new double[2], 1e-3);
//			/////////////var cohesiveMaterial = new BondSlipCohMatUniaxial(K_el, K_pl, 100.0, T_max, new double[2], new double[2], 1e-3);

//			// Create Beam3D Section
//			var beamSection = new BeamSection3D(area, inertiaY, inertiaZ, torsionalInertia, effectiveAreaY, effectiveAreaZ);

//			// element nodes
//			for (int i = 0; i < cntElementConnectivity.GetLength(0); i++)
//			{

//				int elementID = i + (hostElements + embeddedElements); // matrixElements + CNTelements
//				var node1 = model.NodesDictionary[cntElementConnectivity[i, 0] + embeddedNodes].ID;
//				var node2 = model.NodesDictionary[cntElementConnectivity[i, 1] + embeddedNodes].ID;  // matrixNodes + CNTnodes

//				elementNodesClone = new List<INode>();
//				elementNodesClone.Add(model.NodesDictionary[cntElementConnectivity[i, 0]]);
//				elementNodesClone.Add(model.NodesDictionary[cntElementConnectivity[i, 1]]);
//				// element nodes beam
//				elementNodesBeam = new List<INode>();
//				elementNodesBeam.Add(model.NodesDictionary[cntElementConnectivity[i, 0]]);
//				elementNodesBeam.Add(model.NodesDictionary[cntElementConnectivity[i, 1]]);

//				List<INode> nodeSet = new List<INode>()
//				{
//					model.NodesDictionary[node1 - embeddedNodes],
//					model.NodesDictionary[node2 - embeddedNodes],
//					model.NodesDictionary[node1],
//					model.NodesDictionary[node2]
//				};
//				IElementType cohesiveElement = new CohesiveBeam3DToBeam3D(nodeSet, cohesiveMaterial, GaussLegendre1D.GetQuadratureWithOrder(2), elementNodesBeam, elementNodesClone, matrixMaterial, 1, beamSection, CntPerimeter);
//				cohesiveElement.ID = elementID;
				

//				model.ElementsDictionary.Add(cohesiveElement.ID, cohesiveElement);
//				model.SubdomainsDictionary[0].Elements.Add(cohesiveElement);
//			}
//		}

//		private (int[] cntNodeIds, double[,] cntNodeCoordinates, int[,] cntElementConnectivity) GetCntBeamsNodesData(int hostNodes, int hostElements, bool readFromText)
//		{
//			if (readFromText == false)
//			{
				
//				return (new int[5], new double[5, 5], new int[5, 5]);
//			}
//			else
//			{
//				var cntNodeIds = new int[numberOfCnts * 2];
//				var cntNodeCoordinates = new double[numberOfCnts * 2, 3];
//				var cntElementConnectivity = new int[numberOfCnts, 2];

//				var SpecPath = "..\\..\\RveTemplates\\Input\\Output_Input Files";
				

//				string CNTgeometryFileName = "nodes.txt";
//				string CNTconnectivityFileName = "connectivity.txt";

//				string fileNameOnlyCNTgeometryFileName = Path.Combine(SpecPath, Path.GetFileNameWithoutExtension(CNTgeometryFileName));
//				string fileNameOnlyCNTconnectivityFileName = Path.Combine(SpecPath, Path.GetFileNameWithoutExtension(CNTconnectivityFileName));
//				string extension = Path.GetExtension(CNTgeometryFileName);
//				string extension_2 = Path.GetExtension(CNTconnectivityFileName);

//				string currentCNTconnectivityFileName = string.Format("{0}{1}", fileNameOnlyCNTconnectivityFileName, extension_2);

//				string currentCNTgeometryFileName = string.Format("{0}{1}", fileNameOnlyCNTgeometryFileName, extension);

				

//				int CNTNodes = File.ReadLines(currentCNTgeometryFileName).Count();
//				int CNTElems = File.ReadLines(currentCNTconnectivityFileName).Count();

//				// Geometry
//				using (TextReader reader = File.OpenText(currentCNTgeometryFileName))
//				{
//					string CultureName = Thread.CurrentThread.CurrentCulture.Name;
//					CultureInfo ci = new CultureInfo(CultureName);
//					ci.NumberFormat.NumberDecimalSeparator = ".";
//					Thread.CurrentThread.CurrentCulture = ci;
//					for (int i = 0; i < CNTNodes; i++)
//					{
//						string text = reader.ReadLine();
//						string[] bits = text.Split(' ');
//						int nodeID = int.Parse(bits[0]); // matrixNodes
//						double nodeX = double.Parse(bits[1]);
//						double nodeY = double.Parse(bits[2]);
//						double nodeZ = double.Parse(bits[3]);
//						cntNodeIds[i] = nodeID;
//						cntNodeCoordinates[i, 0] = nodeX; cntNodeCoordinates[i, 1] = nodeY; cntNodeCoordinates[i, 2] = nodeZ;
//					}
//				}
//				using (TextReader reader = File.OpenText(currentCNTconnectivityFileName))
//				{
//					for (int i = 0; i < CNTElems; i++)
//					{
//						string text = reader.ReadLine();
//						string[] bits = text.Split(' ');
//						int node1 = int.Parse(bits[0]); // matrixNodes + CNTnodes
//						int node2 = int.Parse(bits[1]); // matrixNodes + CNTnodes
//						cntElementConnectivity[i, 0] = node1; cntElementConnectivity[i, 1] = node2;
//					}
//				}
//				return (cntNodeIds, cntNodeCoordinates, cntElementConnectivity);
//			}
//		}

//		private void AddHexaElements(Model model, int[] nodeIds, double[,] node_coords, int[,] elementConnectivity)
//		{
//			for (int i1 = 0; i1 < nodeIds.Length; i1++)
//			{
//				int nodeID = nodeIds[i1];
//				double nodeCoordX = node_coords[i1, 0];
//				double nodeCoordY = node_coords[i1, 1];
//				double nodeCoordZ = node_coords[i1, 2];

//				model.NodesDictionary.Add(nodeID, new Node(id: nodeID, x: nodeCoordX, y: nodeCoordY, z: nodeCoordZ));
//			}

//			var renumbering = new int[] { 6, 7, 4, 5, 2, 3, 0, 1 };
//			for (int i1 = 0; i1 < elementConnectivity.GetLength(0); i1++)
//			{
//				List<INode> nodeSet = new List<INode>();
//				for (int j = 0; j < 8; j++)
//				{
//					int nodeID = elementConnectivity[i1, renumbering[j]];
//					nodeSet.Add(model.NodesDictionary[nodeID]);
//				}
//				var matrixMaterialClone = (ElasticMaterial3D)this.matrixMaterial;
				
//				var elementFactory = new Multiscale.SupportiveClasses.ContinuumElement3DFactory(new ElasticMaterial3D(matrixMaterialClone.YoungModulus, matrixMaterialClone.PoissonRatio), commonDynamicProperties: null);
//				IElementType e1 = elementFactory.CreateElement(CellType.Hexa8, nodeSet);
//				e1.ID = i1;

//				model.ElementsDictionary.Add(e1.ID, e1);
//				model.SubdomainsDictionary[subdomainID].Elements.Add(e1);
//			}
//		}



//		private int[,] GetHexaRveConnectivity()
//		{
//			int numberOfElements = hexa1 * hexa2 * hexa3;
//			int[,] elementConnectivity = new int[numberOfElements, 8];

//			int counterElement = 0;
//			for (int i = 0; i < hexa3; i++)
//			{
//				for (int j = 0; j < hexa2; j++)
//				{
//					for (int k = 0; k < hexa1; k++)
//					{
//						elementConnectivity[counterElement, 0] = k + (hexa1 + 1) * j + (hexa1 + 1) * (hexa2 + 1) * i;
//						elementConnectivity[counterElement, 1] = elementConnectivity[counterElement, 0] + 1;
//						elementConnectivity[counterElement, 2] = elementConnectivity[counterElement, 1] + (hexa1 + 1);
//						elementConnectivity[counterElement, 3] = elementConnectivity[counterElement, 2] - 1;

//						elementConnectivity[counterElement, 4] = k + (hexa1 + 1) * j + (hexa1 + 1) * (hexa2 + 1) * (i + 1);
//						elementConnectivity[counterElement, 5] = elementConnectivity[counterElement, 4] + 1;
//						elementConnectivity[counterElement, 6] = elementConnectivity[counterElement, 5] + (hexa1 + 1);
//						elementConnectivity[counterElement, 7] = elementConnectivity[counterElement, 6] - 1;

//						counterElement++;
//					}
//				}
//			}

//			return elementConnectivity;
//		}

//		private (int[] NodeIds, double[,] node_coords) GetHexaRveNodesData()
//		{
//			int numberOfNodes = (hexa1 + 1) * (hexa2 + 1) * (hexa3 + 1);
//			int[] NodeIds = new int[numberOfNodes];
//			double[,] node_coords = new double[numberOfNodes, 3];

//			int counterNode = 0;

//			for (int i = 0; i <= hexa3; i++)
//			{
//				for (int j = 0; j <= hexa2; j++)
//				{
//					for (int k = 0; k <= hexa1; k++)
//					{
//						NodeIds[counterNode] = counterNode;
//						node_coords[counterNode, 0] = k * (L01 / hexa1) - L01 * 0.5;
//						node_coords[counterNode, 1] = j * (L02 / hexa2) - L02 * 0.5;
//						node_coords[counterNode, 2] = i * (L03 / hexa3) - L03 * 0.5;
//						counterNode++;
//					}
//				}
//			}

//			return (NodeIds, node_coords);
//		}

//		public NodalDisplacement[] GetModelRigidBodyNodeConstraints(Model model)
//		{
//			NodalDisplacement[] rigidBodyNodeConstraints = new NodalDisplacement[9]
//			{
//				new NodalDisplacement(model.NodesDictionary[0], StructuralDof.TranslationX, amount: 0),
//				new NodalDisplacement(model.NodesDictionary[0], StructuralDof.TranslationY, amount: 0),
//				new NodalDisplacement(model.NodesDictionary[0], StructuralDof.TranslationZ, amount: 0),
//				new NodalDisplacement(model.NodesDictionary[1], StructuralDof.TranslationX, amount: 0),
//				new NodalDisplacement(model.NodesDictionary[1], StructuralDof.TranslationY, amount: 0),
//				new NodalDisplacement(model.NodesDictionary[1], StructuralDof.TranslationZ, amount: 0),
//				new NodalDisplacement(model.NodesDictionary[2], StructuralDof.TranslationX, amount: 0),
//				new NodalDisplacement(model.NodesDictionary[2], StructuralDof.TranslationY, amount: 0),
//				new NodalDisplacement(model.NodesDictionary[2], StructuralDof.TranslationZ, amount: 0),
//			};
			
//			return rigidBodyNodeConstraints;
//		}

//		private List<int> RigidNodes => new List<int>
//		{
//			{ 0},
//			{hexa1 },
//			{(hexa1+1)*hexa2 }
//		};

//		public IRVEbuilder Clone(int a)
//		{
			
//			var cnts = numberOfCnts;
//			return new CntReinforcedElasticNanocomposite(cnts);
//		}

//		public void UpdateCohesiveMaterial()
//		{
//			// define mechanical properties

//			double mi = 8.0;
//			double ni = 8.0;
//			double thickness_CNT = 0.34;
//			double a = 0.241;
//			double diameter_CNT = (a / Math.PI) * Math.Sqrt(Math.Pow(ni, 2) + ni * mi + Math.Pow(mi, 2));
//			double radius_CNT = diameter_CNT / 2.0;
//			double radius_CNT_outer = radius_CNT + (thickness_CNT / 2);
//			double CntPerimeter = 2.0 * Math.PI * radius_CNT_outer;

//			// Create Cohesive Material
//			//var cohesiveMaterial = new BondSlipCohMatUniaxial(K_el, K_pl, 100.0, T_max, new double[2], new double[2], 1e-3);
//			var cohesiveMaterial = new BondSlipMaterial(K_el, K_pl, 100.0, T_max, new double[2], new double[2], 1e-3);


//			// Create Beam3D Section
//			var beamSection = new BeamSection3D(area, inertiaY, inertiaZ, torsionalInertia, effectiveAreaY, effectiveAreaZ);

//			for (int i = 0; i < numberOfCnts; i++)
//			{
//				int elementID = i + (hostElements + embeddedElements);
//				model.ElementsDictionary[elementID] = new CohesiveBeam3DToBeam3D((List<INode>)model.ElementsDictionary[elementID].Nodes, (ICohesiveZoneMaterial)cohesiveMaterial, GaussLegendre1D.GetQuadratureWithOrder(2), elementNodesBeam,
//				elementNodesClone, matrixMaterial, 1, beamSection, CntPerimeter);
//			}
//		}
//	}
//}
