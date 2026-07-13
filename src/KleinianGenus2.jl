module KleinianGenus2

import PolynomialRoots
import Random
using Polynomials
using LinearAlgebra
include("Miscellanea.jl")
include("CurveConstruction.jl")
include("KleinianFunctions.jl")
include("PeriodsCalculation.jl")
include("TripleOfDiscs.jl")
include("AbelMap.jl")

export Genus2WCurve, Genus2WCurveDivisor

end # module KleinianGenus2
