module HilbertSpace

using Combinatorics
using MLStyle



include("./quantum_statistics.jl")
export Quantum_Statistics, Fermionic, Bosonic

include("./quantum_state.jl")
export Abstract_State, Single_Particle_State
export sort_inplace_with_sign_flipped_fast!, sort_inplace_fast!


include("./hilbert_space.jl")
export Abstract_Hilbert_Space, Single_Particle_Hilbert_Space, Many_Particle_Hilbert_Space


end # module HilbertSpace
