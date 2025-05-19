module ModuleLEPEE

using Distributions
using SpecialFunctions
using Interpolations
using Optim

include("ModuleLEPEE_aux.jl")
include("ModuleLEPEE_Rain.jl")
include("ModuleLEPEE_Wind.jl")

export MAT_m, fit_curve, coeffs, LTPrediction,Cr_dsd

mat_increment = 10^-6
mats = LinRange(1e-14,1e-3,100)

tyear=(365*24*60*60)


"""
    fit_curve(LT, R)

This function fits a power law curve to the LT data using non-linear least squares optimization.
It minimizes the squared differences between the observed LT values and the model defined as:
    
    LT ≈ x[1] * (R .^ x[2])
    
# Arguments
- LT: Array of observed dependent values.
- R: Array of independent variable values.

# Returns
An optimization result object containing the optimal parameters.
"""
function fit_curve(LT::Reals, R::Reals)
    objective = x -> sum((LT .- x[1] * (R .^ x[2])).^2)
    gradient!(g, x) = begin
        diff1 = -2 .* (LT .- x[1] * (R .^ x[2])).*(R .^ x[2])
        diff2 = -2 .* (LT .- x[1] * (R .^ x[2])).*(x[1] * (R .^ x[2]) .* log.(R))
        g[1] = sum(diff1)
        g[2] = sum(diff2)
    end
    result = optimize(objective, gradient!, [4e8, -5.7], BFGS(), Optim.Options(g_tol=1e-15, f_tol=1e-15, x_tol=1e-15))
    return result
end


"""
    coeffs(res, dsd::Function, wind::Function)

This function computes the coefficients necessary for LT prediction based on the optimization results and the provided droplet size distribution and wind functions. It calculates the exponent m, the wind coefficient (Cw), the rain coefficient (Cr), and the material coefficient M.

# Arguments
- res: The optimization result object from fit_curve.
- dsd: A function providing droplet size distribution corrections.
- wind: A function providing wind corrections.

# Returns
A named tuple with:
- m: The computed exponent.
- Cw: The wind correction factor.
- Cr: The rain correction factor.
- M: The final computed coefficient.
"""
function coeffs(res, dsd::Function, wind::Function)

    m = -res.minimizer[2] - 1

    Cw = wind(m);

    Cr = Cr_dsd(dsd)

    M = res.minimizer[1] * tyear * Cw * Cr / 8.9

    M = M ^ (1/m)

    return (m=m,Cw=Cw,Cr=Cr,M=M)

end

function LTPrediction(coeff, dsd::Function, wind::Function)
    A = 8.9 * coeff.M ^coeff.m /(tyear * Cr_dsd(dsd) * wind(coeff.m) );
    B = -coeff.m -1.0;
    res=(A=A,B=B);
    return(res);
end

end # module