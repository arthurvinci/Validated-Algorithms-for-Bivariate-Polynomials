module ValidatedHensel

using AbstractAlgebra
using Nemo

################################################################################
#
#   Types
#
################################################################################

CFloat = RDF
CFloatElem = typeof(CFloat(0))
CBall = ArbField(64)
CBallElem = ArbFieldElem

export CFloat, CFloatElem
export CBall, CBallElem

################################################################################
#
#   Methods
#
################################################################################

include("math_operations.jl")
include("fast_div_rem.jl")
include("validated_fast_div_rem.jl")
include("hensel_lifting.jl")
include("validated_hensel_lifting.jl")


end # module ValidatedHensel
