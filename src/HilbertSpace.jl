module HilbertSpace

using MLStyle

include("./quantum_statistics.jl") # where struct `Fermionic <: Quantum_Statistics` and `Bosonic <: Quantum_Statistics` are defined
include("./quantum_state.jl") # where struct `Single_Particle_State <: Abstract_State` and `Many_Particle_State <: Abstract_State` are defined
include("./hilbert_space.jl") # where struct `Finite_Single_Particle_Hilbert_Space <: Finite_Hilbert_Space` and `Finite_Many_Particle_Hilbert_Space <: Finite_Hilbert_Space` are defined

# export structs
export Quantum_Statistics, Fermionic, Bosonic
export Abstract_State, Single_Particle_State, Many_Particle_State
export Abstract_Hilbert_Space, Finite_Hilbert_Space, Finite_Single_Particle_Hilbert_Space, Finite_Many_Particle_Hilbert_Space

# export methods
export initialize_finite_single_particle_Hilbert_space, initialize_finite_fermionic_many_particle_Hilbert_space



end # module HilbertSpace
