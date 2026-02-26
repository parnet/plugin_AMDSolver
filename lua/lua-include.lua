
-- Extend utilities
util = util or {}
util.amdsolver = util.amdsolver or {}

function util.amdsolver.create_conv_check(desc)
    desc.convCheck.absolute = desc.convCheck.absolute or 1e-14
    desc.convCheck.reduction = desc.convCheck.reduction or 1e-10
    desc.convCheck.divergence = desc.convCheck.divergence or 1e20

    desc.convCheck.iterations = desc.convCheck.iterations or 100
    desc.convCheck.iterations_min = desc.convCheck.iterations_min or 0
    desc.convCheck.verbose = desc.convCheck.verbose or False

    -- print(desc.convCheck)
    conv_check = AMD_IterationControl()
    conv_check:init(desc.convCheck.absolute,desc.convCheck.reduction, desc.convCheck.divergence,desc.convCheck.iterations_min,desc.convCheck.iterations)
    conv_check:set_verbose(desc.convCheck.verbose)
    return conv_check
end

function util.amdsolver.create_preconditioner(desc)
    precond_type = desc.precond.type

    if precond_type == "amg" then

        preconditioner = AMD_RugeStuebenAMG()
        preconditioner:set_coarsening_strategy("pmis")
        preconditioner:set_interpolation_type("extpi")
        preconditioner:set_cycle("V")

        return preconditioner

    elseif precond_type == "jac" then
        preconditioner = AMD_Jacobi();
        return preconditioner
    end
end


function util.amdsolver.create_solver(name, desc)
    preconditioner = nil

    if desc.convCheck ~= nil then
        conv_check = util.amdsolver.create_conv_check(desc)
    end

    if desc.precond ~= nil then
        preconditioner = util.amdsolver.create_preconditioner(desc)
    end

    if name == "amd_bicgstab" then
        solver = AMD_BiCGStab();
        solver:set_conv_check(conv_check)
        solver:set_preconditioner(preconditioner)
        linSolver = AgglomeratingSolver(solver);
        return linSolver
    elseif name == "amd_amg" then
        solver = AMD_RugeStuebenAMG()
        solver:set_conv_check(conv_check)
        solver:set_coarsening_strategy("pmis")
        solver:set_interpolation_type("extpi")
        solver:set_cycle("V")

        linSolver = AgglomeratingSolver(solver);

    elseif name == "amd_lu" then
        linSolver = AgglomeratingSolver(AMD_LU());
        return linSolver

    elseif name == "amd_qr" then
        linSolver = AgglomeratingSolver(AMD_QR());
        return linSolver

    elseif name == "amd_qr" then
        linSolver = AgglomeratingSolver(AMD_Inversion());
        return linSolver
    else
        print("unknown solver name")
    end
end