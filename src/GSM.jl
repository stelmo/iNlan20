module GSM # Genome scale model code

# import packages
using BioSequences
using JSON
using ExcelReaders
using DataValues
using Statistics
using FASTX

include("parsetools.jl")
include("model.jl")
include("curaterxnsmets.jl")

end # module
