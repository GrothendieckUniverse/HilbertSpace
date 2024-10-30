abstract type Abstract_State end

#=====================================================================================================================================#
#======================================================== Single-Particle State ======================================================#
#=====================================================================================================================================#
"""
General Single-Particle State `Single_Particle_State{S<:Quantum_Statistics}` 
---
with `S` labelling the quantum statistics of the single-particle state.
- Fields:
    - `statistics::S`: quantum statistics of the single-particle state
    - `dof_indices::Tuple`: tuple of indices to each physical degree of freedom. It can be either numerical (finite) or symbolic
"""
struct Single_Particle_State{S} <: Abstract_State where {S<:Quantum_Statistics}
    statistics::S
    dof_indices::Tuple # tuple of indices to each physical degree of freedom. It can be either numerical (finite) or symbolic
end
"add a show method for `Single_Particle_State`"
Base.show(io::IO, state::Single_Particle_State) = @match state.statistics begin
    Fermionic() => print(io, "|f_$(state.dof_indices)⟩")
    Bosonic() => print(io, "|b_$(state.dof_indices)⟩")
end


# add comparison between `Single_Particle_State`
Base.:(==)(state1::Single_Particle_State, state2::Single_Particle_State) = state1.statistics == state2.statistics && state1.dof_indices == state2.dof_indices
Base.:>(state1::Single_Particle_State, state2::Single_Particle_State) = begin
    if state1.statistics > state2.statistics
        return true
    else # state1.statistics <= state2.statistics
        if state1.statistics == state2.statistics
            return state1.dof_indices > state2.dof_indices
        else # state1.statistics < state2.statistics
            return false
        end
    end
end
Base.:<(state1::Single_Particle_State, state2::Single_Particle_State) = begin
    if state1.statistics < state2.statistics
        return true
    else # state1.statistics >= state2.statistics
        if state1.statistics == state2.statistics
            return state1.dof_indices < state2.dof_indices
        else # state1.statistics > state2.statistics
            return false
        end
    end
end


# TODO: add arithmic operations for `Single_Particle_State` to construct (symbolic) polynomials




#=====================================================================================================================================#
#========================================================= Many-Particle State =======================================================#
#=====================================================================================================================================#
"""
General Many Particle State `Many_Particle_State{N<:Integer}`
---
with `N` the number of single-particle states. Note: this representation even support *infinitely* (or even *uncountably infinitely*) occupied many-particle states!
"""
struct Many_Particle_State{N,S} <: Abstract_State where {N<:Integer,S<:Quantum_Statistics}
    single_particle_states::NTuple{N,<:Single_Particle_State{S}}
end
"add a show method for `Many_Particle_State{N}`"
Base.show(io::IO, state::Many_Particle_State{N}) where {N} = begin
    for (i, ψ_single) in enumerate(state.single_particle_states)
        print(io, "$ψ_single")
        if i != N
            print(io, " ⊗ ")
        end
    end
end
"add a method to conveniently construct `Many_Particle_State{S,N}` with a tuple of single-particle states"
Many_Particle_State(single_particle_states::NTuple{N,<:Single_Particle_State{S}}) where {N<:Integer,S<:Quantum_Statistics} = Many_Particle_State{N,S}(single_particle_states)
Many_Particle_State(single_particle_states::Vector{<:Single_Particle_State{S}}) where {S<:Quantum_Statistics} = Many_Particle_State{length(single_particle_states),S}(Tuple(single_particle_states))