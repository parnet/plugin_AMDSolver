/*
 * Copyright (c) 2011-2016:  G-CSC, Goethe University Frankfurt
 * Author:
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

#include "bridge/util.h"

// if possible, replace this with util_domain_dependent.h or
// util_algebra_dependent.h to speed up compilation time
#include "bridge/util_domain_algebra_dependent.h"

#include "iterative/iterative_solver.hpp"
#include "direct_solver/inversion.hpp"
#include "direct_solver/lu.hpp"
#include "direct_solver/qr.hpp"
#include "iterative/bicgstab.hpp"
#include "preconditioner/preconditioner.hpp"
#include "multigrid/ruge_stueben_amg.hpp"
#include "driver.hpp"
#include "iteration_control.hpp"

using namespace std;
using namespace ug::bridge;

namespace ug{
namespace AMDSolver{




/** 
 *  \defgroup plugin_template Plugin Template
 *  \ingroup plugins_experimental
 *  This is a template for new plugins.
 *  \{
 */

/**
 * Class exporting the functionality. All functionality that is to
 * be used in scripts or visualization must be registered here.
 */
struct Functionality
{

/**
 * Function called for the registration of Domain and Algebra dependent parts.
 * All Functions and Classes depending on both Domain and Algebra
 * are to be placed here when registering. The method is called for all
 * available Domain and Algebra types, based on the current build options.
 *
 * @param reg				registry
 * @param grp       		group for sorting of functionality
 */
template <typename TDomain, typename TAlgebra>
static void DomainAlgebra(Registry& reg, string grp)
{
	string suffix = GetDomainAlgebraSuffix<TDomain,TAlgebra>();
	string tag = GetDomainAlgebraTag<TDomain,TAlgebra>();

}

/**
 * Function called for the registration of Domain dependent parts.
 * All Functions and Classes depending on the Domain
 * are to be placed here when registering. The method is called for all
 * available Domain types, based on the current build options.
 *
 * @param reg				registry
 * @param parentGroup		group for sorting of functionality
 */
template <typename TDomain>
static void Domain(Registry& reg, string grp)
{
	string suffix = GetDomainSuffix<TDomain>();
	string tag = GetDomainTag<TDomain>();

}

/**
 * Function called for the registration of Dimension dependent parts.
 * All Functions and Classes depending on the Dimension
 * are to be placed here when registering. The method is called for all
 * available Dimension types, based on the current build options.
 *
 * @param reg				registry
 * @param parentGroup		group for sorting of functionality
 */
template <int dim>
static void Dimension(Registry& reg, string grp)
{
	string suffix = GetDimensionSuffix<dim>();
	string tag = GetDimensionTag<dim>();

}

/**
 * Function called for the registration of Algebra dependent parts.
 * All Functions and Classes depending on Algebra
 * are to be placed here when registering. The method is called for all
 * available Algebra types, based on the current build options.
 *
 * @param reg				registry
 * @param grp				group for sorting of functionality
 */
template <typename TAlgebra>
static void Algebra(Registry& reg, string grp)
{
	string suffix = GetAlgebraSuffix<TAlgebra>();
	string tag = GetAlgebraTag<TAlgebra>();

	// Direct solver methods
	{
		using T = AMD_LU<TAlgebra>;
		using TBase1 = IExternalSolver<TAlgebra>;
		using TBase2 = ILinearOperatorInverse<typename TAlgebra::vector_type> ;
		string name = string("AMD_LU").append(suffix);
		reg.add_class_<T,TBase1, TBase2>(name, grp)
			.add_constructor()

			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "AMD_LU", tag);
	}	{
		using T = AMD_Inversion<TAlgebra>;
		using TBase1 = IExternalSolver<TAlgebra>;
		using TBase2 = ILinearOperatorInverse<typename TAlgebra::vector_type>;
		string name = string("AMD_Inversion").append(suffix);
		reg.add_class_<T,TBase1, TBase2>(name, grp)
			.add_constructor()

			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "AMD_Inversion", tag);
	}	{
		using T = AMD_QR<TAlgebra>;
		using TBase1 = IExternalSolver<TAlgebra>;
		using TBase2 = ILinearOperatorInverse<typename TAlgebra::vector_type>;
		string name = string("AMD_QR").append(suffix);
		reg.add_class_<T,TBase1, TBase2>(name, grp)
			.add_constructor()

			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "AMD_QR", tag);
	}

	{
		using T = AMD_BiCGStab<TAlgebra>;
		using TBase1 = IExternalSolver<TAlgebra>;
		using AT = AMD_IterativeLinearSolver;
		string name = string("AMD_BiCGStab").append(suffix);
		reg.add_class_<T, AT, TBase1>(name, grp)
			.add_constructor()
			.add_method("set_preconditioner", &T::set_preconditioner, "", "", "")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "AMD_BiCGStab", tag);
	}

	{
		using T = AMD_RugeStuebenAMG<TAlgebra>;
		using TBase1 = IExternalSolver<TAlgebra>;
		using AT = AMD_IterativeLinearSolver;

		string name = string("AMD_RugeStuebenAMG").append(suffix);
		reg.add_class_<T,TBase1, AT>(name, grp)
			.add_constructor()
			.add_method("set_preconditioner", &T::set_preconditioner, "", "", "")
			.add_method("set_smoother", &T::set_smoother, "", "", "")
			.add_method("set_coarsest_level", &T::set_coarsest_level, "", "", "")
			.add_method("set_strength_threashold", &T::set_strength_threashold, "", "", "")
			.add_method("set_coarsening_strategy", &T::set_coarsening_strategy, "", "", "")
			.add_method("set_interpolation_type", &T::set_interpolation_type, "", "", "")
			.add_method("set_interpolation_ff1_limit", &T::set_interpolation_ff1_limit, "", "", "")
			.add_method("re_build_numeric", &T::re_build_numeric, "", "", "")
			.add_method("set_cycle", &T::set_cycle, "", "", "")
			.add_method("set_smooth_num", &T::set_cycle, "", "", "")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "AMD_RugeStuebenAMG", tag);
	}
}

/**
 * Function called for the registration of Domain and Algebra independent parts.
 * All Functions and Classes not depending on Domain and Algebra
 * are to be placed here when registering.
 *
 * @param reg				registry
 * @param grp		  		group for sorting of functionality
 */
static void Common(Registry& reg, string grp)
{

	{
		using T = ug::AMD_IterationControl;
		string name = string("AMD_IterationControl");
		reg.add_class_<T>(name, grp)
			.add_constructor()
			.add_method("set_verbose",&T::set_verbose)
			.add_method("init",&T::init)
			.add_method("init_tol",&T::init_tol)
			.add_method("init_min_iter",&T::init_min_iter)
			.add_method("init_max_iter",&T::init_max_iter)
			.add_method("set_residual_norm",&T::set_residual_norm)
			.set_construct_as_smart_pointer(true);
		// reg.add_class_to_group(name, "AMD_Solver", tag);
	}
	{
		using T = AMD_Solver;
		string name = string("AMD_Solver");
		reg.add_class_<T>(name, grp)
			.set_construct_as_smart_pointer(true);
		// reg.add_class_to_group(name, "AMD_Solver", tag);
	}
	{
		using T = AMD_Preconditioner;
		string name = string("AMD_Preconditioner");
		reg.add_class_<T>(name, grp)
			.set_construct_as_smart_pointer(true);
		// reg.add_class_to_group(name, "AMD_Preconditioner", tag);
	}	{
		using T = ug::AMD_IterativeLinearSolver;
		using TBase = ug::AMD_Solver;
		string name = string("AMD_IterativeLinearSolver");
		reg.add_class_<T,TBase>(name, grp)
			.add_method("set_conv_check",&T::set_conv_check)
			.set_construct_as_smart_pointer(true);
	}
	{
		using T = AMD_Jacobi;
		using TB = AMD_Preconditioner;
		string name = string("AMD_Jacobi");
		reg.add_class_<T, TB>(name, grp)
			.add_constructor()
			.set_construct_as_smart_pointer(true);
		//reg.add_class_to_group(name, "AMD_Jacobi");
	}
//	The code below shows how a simple function can be registered
	// reg.add_function("TemplateSampleFunction", &TemplateSampleFunction, grp,
	// 				 "", "", "Prints a short message");
}

}; // end Functionality

// end group plugin_template
/// \}

} // end namespace TemplatePlugin


/**
 * This function is called when the plugin is loaded.
 */
extern "C" void
InitUGPlugin_AMDSolver(Registry* reg, string grp)
{
	grp.append("AMDSolver");
	using Functionality = AMDSolver::Functionality;
	//if (!rocalution::is_initialized_rocalution()) {
		ug::rocalution_init();
	//}

	try{
		RegisterCommon<Functionality>(*reg,grp);
		RegisterDimensionDependent<Functionality>(*reg,grp);
		RegisterDomainDependent<Functionality>(*reg,grp);
		RegisterAlgebraDependent<Functionality>(*reg,grp);
		RegisterDomainAlgebraDependent<Functionality>(*reg,grp);
	}
	UG_REGISTRY_CATCH_THROW(grp);
}

}// namespace ug
