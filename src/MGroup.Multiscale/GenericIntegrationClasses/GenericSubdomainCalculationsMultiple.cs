using MGroup.LinearAlgebra.Matrices;
using MGroup.LinearAlgebra.Vectors;
using MGroup.MSolve.Discretization.Dofs;
using MGroup.MSolve.Discretization.Entities;
using MGroup.MSolve.Discretization.Providers;
using MGroup.MSolve.MultiscaleAnalysis.Interfaces;
using MGroup.MSolve.Solution;
using MGroup.MSolve.Solution.LinearSystem;
using MGroup.Solvers.AlgebraicModel;
using MGroup.Solvers.LinearSystem;
using System.Collections.Generic;
using System.Linq;

namespace MGroup.MSolve.MultiscaleAnalysis
{
    /// <summary>
    /// Supportive class  that implements nesessary integration methods associated with FE2 multiscale analysis 
    /// Authors: Gerasimos Sotiropoulos
    /// </summary>
    public static class GenericSubdomainCalculationsMultiple<TMatrix> where TMatrix : class, IMatrix
    {
        
        #region v2 methods
        public static double[][] CombineMultipleSubdomainsIntegrationVectorsIntoTotal(Dictionary<int, double[][]> VectorsSubdomains, IScaleTransitions scaleTransitions)
        {


            double[][] totalVectors = new double[scaleTransitions.MacroscaleVariableDimension()][]; //or VectorsSubdomains.getLength(0);
            int oneSubdomainID = VectorsSubdomains.Keys.ElementAt(0);
            for (int i1 = 0; i1 < scaleTransitions.MacroscaleVariableDimension(); i1++)
            {
                totalVectors[i1] = new double[VectorsSubdomains[oneSubdomainID][i1].GetLength(0)];
            }


            foreach (int subdomainID in VectorsSubdomains.Keys)
            {
                for (int i1 = 0; i1 < scaleTransitions.MacroscaleVariableDimension(); i1++)
                {
                    for (int i2 = 0; i2 < VectorsSubdomains[subdomainID][i1].GetLength(0); i2++)
                    {
                        totalVectors[i1][i2] += VectorsSubdomains[subdomainID][i1][i2];
                    }
                }
            }

            return totalVectors;

        }

        public static Dictionary<int, double[]> CalculateFppReactionsVectorSubdomains(Model model, IElementMatrixProvider elementProvider,
            IScaleTransitions scaleTransitions, Dictionary<int, INode> boundaryNodes, List<int> boundaryNodesIds, IGlobalVector solution,/* Dictionary<int, IVector> dSolution,*/
            Dictionary<int, Dictionary<IDofType, double>> initialConvergedBoundaryDisplacements, Dictionary<int, Dictionary<IDofType, double>> totalBoundaryDisplacements,
            int nIncrement, int totalIncrements, GlobalAlgebraicModel<TMatrix> globalAlgebraicModel)
        {
            Dictionary<int, double[]> FppReactionVectorSubdomains = new Dictionary<int, double[]>();

            
            FppReactionVectorSubdomains.Add(model.EnumerateSubdomains().First().ID, GenericSubdomainCalculations<TMatrix>.CalculateFppReactionsVector((Subdomain)model.EnumerateSubdomains().First(), elementProvider, scaleTransitions, boundaryNodes, 
                    boundaryNodesIds, solution, initialConvergedBoundaryDisplacements, totalBoundaryDisplacements, nIncrement, totalIncrements, globalAlgebraicModel));
            

            return FppReactionVectorSubdomains;

        }

        internal static double[] CombineMultipleSubdomainsStressesIntegrationVectorsIntoTotal(Dictionary<int, double[]> fppReactionVectorSubdomains)
        {

            double[] totalVector = new double[fppReactionVectorSubdomains.ElementAt(0).Value.GetLength(0)];

            foreach (int subdomainID in fppReactionVectorSubdomains.Keys)
            {
                for (int i1 = 0; i1 < totalVector.GetLength(0); i1++)
                {
                    totalVector[i1] += fppReactionVectorSubdomains[subdomainID][i1];
                }
            }

            return totalVector;

        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="KfpDqSubdomains"></param>
        /// <param name="model"></param>
        /// <param name="elementProvider"></param>
        /// <param name="scaleTransitions"></param>
        /// <param name="boundaryNodes"></param>
        /// <param name="solver">
        /// <paramref name="solver"/>.<see cref="ISolver.Initialize"/> must already have been called. Also the linear system matrices must already have been set.
        /// </param>
        public static IGlobalVector[] CalculateKffinverseKfpDqSubdomains(IGlobalVector[] KfpDqVectors, Model model, IElementMatrixProvider elementProvider,
            IScaleTransitions scaleTransitions, Dictionary<int, INode> boundaryNodes, ISolver solver)
        {
            // IReadOnlyDictionary<int, ILinearSystem> linearSystems = solver.LinearSystems;

            #region Creation of solution vectors structure
            IGlobalVector[] f2_vectorsSubdomains = new IGlobalVector[KfpDqVectors.Length];
            //foreach (int subdomainID in KfpDqSubdomains.Keys)
            //{
            //    f2_vectorsSubdomains.Add(subdomainID, new double[KfpDqSubdomains[subdomainID].GetLength(0)][]);
            //}
            #endregion

            //#region Creation of linear systems with no RHS (first RHS value can be assigned too )
            //ILinearSystem[] seclinearSystems = new ILinearSystem[linearSystems.Count];
            //int counter = 0;
            //foreach (ILinearSystem subdomain in linearSystems.Values)
            //{
            //    seclinearSystems[counter] = new SkylineLinearSystem(subdomain.ID, new double[KfpDqSubdomains[subdomain.ID][0].GetLength(0)]);
            //    seclinearSystems[counter].Matrix = subdomain.Matrix;
            //}
            //#endregion

            //#region creation of solver
            //var secSolver = new SolverSkyline(seclinearSystems[0]);
            //secSolver.Initialize();

            //#endregion

            #region Consecutively(for macroscaleVariableDimension times) Set proper right hand side. Solve. Copy solution in output vector 
            //int oneSubomainID = linearSystems.First().Value.Subdomain.ID;           //seclinearSystems[0].ID;
            for (int k = 0; k < scaleTransitions.MacroscaleVariableDimension(); k++) //KfpDqSubdomains[linearSystems[0].ID].GetLength(0)=Mac
            {
                #region Set proper RHS 
                //var globalRHS = new Vector(model.TotalDOFs); //TODO: uncoomment if globalRHS is needed for solver
                //foreach (ILinearSystem secSubdomain in linearSystems.Values)
                //{

                //    secSubdomain.RhsVector = Vector.CreateFromArray(KfpDqSubdomains[secSubdomain.Subdomain.ID][k], false);
                //    //secSubdomain.RhsVector = Vector.CreateFromArray(KfpDqSubdomains[secSubdomain.Subdomain.ID][k], true); Wste sigoura na mhn peiraxthei to double[]


                //    //mappings[seclinearSystems.Select((v, i) => new { System = v, Index = i }).First(x => x.System.ID == secSubdomain.ID).Index].SubdomainToGlobalVector(subdomainRHS.Data, globalRHS.Data);
                //    //TODO: uncoomment if globalRHS is needed for solver
                //}
                solver.LinearSystem.RhsVector.Clear();
                solver.LinearSystem.RhsVector.CopyFrom(KfpDqVectors[k]);
                #endregion

                #region Solve
                solver.Solve();
                #endregion

                #region Copy solution in output vector
                //foreach (ILinearSystem secSubdomain in linearSystems.Values)
                //{
                //    f2_vectorsSubdomains[secSubdomain.Subdomain.ID][k] = secSubdomain.Solution.CopyToArray();
                //}
                f2_vectorsSubdomains[k] = solver.LinearSystem.Solution.Copy();
                #endregion
            }
            #endregion

            return f2_vectorsSubdomains;
        }

        // TODO GEr maintain te following lines for future multi subdomain implemntation
         //public static Dictionary<int, double[][]> CalculateKpfKffinverseKfpDqSubdomains(Dictionary<int, double[][]> f2_vectorsSubdomains, Model model, IElementMatrixProvider elementProvider, IScaleTransitions scaleTransitions, Dictionary<int, Node> boundaryNodes)
        //{
        //    Dictionary<int, double[][]> f3_vectorsSubdomains = new Dictionary<int, double[][]>();

        //    foreach (Subdomain subdomain in model.SubdomainsDictionary.Values)
        //    {
        //        f3_vectorsSubdomains.Add(subdomain.ID, SubdomainCalculations.CalculateKpfKffinverseKfpDq(f2_vectorsSubdomains[subdomain.ID], subdomain, elementProvider, scaleTransitions, boundaryNodes));
        //    }
        //    return f3_vectorsSubdomains;
        //}
        #endregion
    }
}
