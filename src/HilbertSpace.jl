module HilbertSpace

using Combinatorics
using MLStyle
using StaticArrays


abstract type Abstract_Hilbert_Space end
abstract type Abstract_State end


include("./quantum_statistics.jl") # where struct `Fermionic <: Quantum_Statistics` and `Bosonic <: Quantum_Statistics` are defined
include("./quantum_state.jl") # where struct `Single_Particle_State <: Abstract_State` and `Many_Particle_State <: Abstract_State` are defined, and methods like `sort()`, `sort_sgn()` are defined (comparison traits are also defined)
include("./hilbert_space.jl") # where struct `Finite_Dimensional_Single_Particle_Hilbert_Space <: Finite_Hilbert_Space` and `Finite_Dimensional_Many_Particle_Hilbert_Space <: Finite_Hilbert_Space` are defined


# export structs
export Quantum_Statistics, Fermionic, Bosonic
export Abstract_State, Single_Particle_State, Many_Particle_State

export Abstract_Hilbert_Space, Finite_Dimensional_Hilbert_Space, Finite_Dimensional_Single_Particle_Hilbert_Space, Finite_Dimensional_Many_Particle_Hilbert_Space


# export methods
export sort_inplace_with_sign_flipped_fast!



end # module HilbertSpace
