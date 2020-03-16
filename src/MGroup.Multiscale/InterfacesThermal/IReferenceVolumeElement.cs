using System;
using System.Collections.Generic;
using System.Text;

using MGroup.LinearAlgebra.Matrices;
using MGroup.MSolve.Discretization;
//using MGroup.MSolve.Discretization.Interfaces;

namespace MGroup.Multiscale.InterfacesThermal
{
	public interface IReferenceVolumeElement
	{
		void ApplyBoundaryConditions();
		IMatrixView CalculateKinematicRelationsMatrix(ISubdomain subdomain);
		double CalculateRveVolume();
	}
}



