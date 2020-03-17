using System.Collections.Generic;
//using MGroup.Analyzers.Interfaces;
using MGroup.LinearAlgebra.Matrices;
using MGroup.MSolve.AnalysisWorkflow.Providers;
using MGroup.MSolve.Discretization;
using MGroup.MSolve.Solution.LinearSystems;
//using MGroup.MSolve.Discretization.Interfaces;
using MGroup.Solvers;
using MGroup.Solvers.LinearSystems;
using MGroup.NumericalAnalyzers;
using MGroup.MSolve.Solution;
using MGroup.Multiscale.InterfacesThermal;

namespace MGroup.Multiscale.Analyzers
{
	public class HomogenizationAnalyzer
	{
		private readonly IReadOnlyDictionary<int, ILinearSystem> linearSystems;
		private readonly IModel model;
		private readonly IStaticProvider provider;
		private readonly IReferenceVolumeElement rve;
		private readonly ISolver solver;

		public HomogenizationAnalyzer(IModel model, ISolver solver, IStaticProvider provider, 
			IReferenceVolumeElement rve)
		{
			this.model = model;
			this.linearSystems = solver.LinearSystems;
			this.solver = solver;
			this.provider = provider;
			this.rve = rve;
		}

		public Dictionary<int, IMatrix> EffectiveConstitutiveTensors { get; private set; }
		
		public void Initialize()
		{
			// The order in which the next initializations happen is very important.
			rve.ApplyBoundaryConditions();
			model.ConnectDataStructures();
			solver.OrderDofs(true);
			foreach (ILinearSystem linearSystem in linearSystems.Values)
			{
				linearSystem.Reset(); // Necessary to define the linear system's size 
			}
		}

		public void Solve()
		{
			var matricesFreeConstr = new Dictionary<int, IMatrixView>();
			var matricesConstrFree = new Dictionary<int, IMatrixView>();
			var matricesConstrConstr = new Dictionary<int, IMatrixView>();

			// Build all matrices for free and constrained dofs.
			foreach (ILinearSystem linearSystem in linearSystems.Values)
			{
				int id = linearSystem.Subdomain.ID;
				(IMatrixView matrixFreeFree, IMatrixView matrixFreeConstr, IMatrixView matrixConstrFree,
					IMatrixView matrixConstrConstr) = provider.CalculateSubMatrices(linearSystem.Subdomain);
				linearSystem.Matrix = matrixFreeFree;

				matricesFreeConstr[id] = matrixFreeConstr;
				matricesConstrFree[id] = matrixConstrFree;
				matricesConstrConstr[id] = matrixConstrConstr;
			}

			// Static condensation: Acond = Acc - Acf * inv(Aff) * Afc
			Dictionary<int, Matrix> invAffTimesAfc = solver.InverseSystemMatrixTimesOtherMatrix(matricesFreeConstr);
			var condensedMatrices = new Dictionary<int, IMatrix>();
			foreach (ILinearSystem linearSystem in linearSystems.Values)
			{
				int id = linearSystem.Subdomain.ID;
				IMatrix condensedMatrix = matricesConstrConstr[id].Subtract(
					matricesConstrFree[id].MultiplyRight(invAffTimesAfc[id]));
				condensedMatrices[id] = condensedMatrix;
			}

			// Calculate effective elasticity/conductivity/whatever tensor: C = 1 / V * (D * Acond * D^T)
			EffectiveConstitutiveTensors = new Dictionary<int, IMatrix>();
			foreach (ILinearSystem linearSystem in linearSystems.Values)
			{
				int id = linearSystem.Subdomain.ID;
				IMatrixView kinematicRelationsMatrix = rve.CalculateKinematicRelationsMatrix(linearSystem.Subdomain);
				Matrix effectiveTensor = kinematicRelationsMatrix.ThisTimesOtherTimesThisTranspose(condensedMatrices[id]);
				double rveVolume = rve.CalculateRveVolume();
				effectiveTensor.ScaleIntoThis(1.0 / rveVolume);
				EffectiveConstitutiveTensors[id] = effectiveTensor;
			}
		}
	}
}
