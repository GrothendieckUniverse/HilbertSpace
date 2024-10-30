
abstract type Abstract_Hilbert_Space end
abstract type Finite_Hilbert_Space <: Abstract_Hilbert_Space end


#=====================================================================================================================================#
#=================================================== Single-Particle Hilbert Space ===================================================#
#=====================================================================================================================================#
"""
Finite-dimensional Single-particle Hilbert Space `Finite_Single_Particle_Hilbert_Space{D<:Integer} <: Finite_Hilbert_Space`
---
with `D<:Integer` the number of indices to label the single-particle states.
- Fields:
    - `statistics::Quantum_Statistics`: quantum statistics of the single-particle states
    - `ndof::Integer`: number of degrees of freedom for the single-particle state (please distinguish from `dof_ndof` below!)
    - `dof_name::NTuple{<:Integer,String}`: tuple of name for each d.o.f of the single-particle state
    - `dof_ndof::NTuple{<:Integer,<:Integer}`: tuple of the number of degrees of freedom for each d.o.f of the single-particle state. For example, `dof_ndof = (2, 3)` means the first d.o.f has 2 states and the second has 3 states
    - `state_generator::Base.Iterators.ProductIterator`: lazy iterator of single-particle states. Note: to save memory, we only store iterator of the single-particle states, rather than the explicit full list of them
    - `state_to_index_map::Dict{<:Single_Particle_State,<:Integer}`: hashmap `single-particle state -> index`
    - `nstate::Integer`: number of single-particle states
"""
struct Finite_Single_Particle_Hilbert_Space{D} <: Finite_Hilbert_Space where {D<:Integer}
    statistics::Quantum_Statistics # quantum statistics of the single-particle states
    ndof::Integer # number of degrees of freedom for the single-particle state: always explicitly stored as the parametric `D<:Integer` (please distinguish from `dof_ndof` below!)
    dof_name::NTuple{D,String} # tuple of name for each d.o.f of the single-particle state
    dof_ndof::NTuple{D,<:Integer} # tuple of the number of degrees of freedom for each d.o.f of the single-particle state. For example, `dof_ndof = (2, 3)` means the first d.o.f has 2 states and the second has 3 states

    state_generator::Base.Generator # lazy generator of single-particle states. Note: to save memory, we only store the generator of the single-particle states, rather than the explicit full list of them
    nstate::Integer # number of single-particle states
    state_to_index_map::Dict{<:Single_Particle_State,<:Integer} # hashmap `single-particle state -> index`
end


function initialize_finite_single_particle_Hilbert_space(;
    statistics::Quantum_Statistics,
    dof_ndof::NTuple{D,<:Integer},
    dof_name::Tuple=()
) where {D}
    if isempty(dof_name)
        dof_name = Tuple("dof_$i" for i in 1:D)
    else
        @assert eltype(dof_name) == String "The type of `dof_name` must be `String`"
        @assert length(dof_name) == D
    end

    single_particle_state_generator = (Single_Particle_State(statistics, finite_dof_indices) for finite_dof_indices in Iterators.product((1:n for n in dof_ndof)...))
    nstate = length(single_particle_state_generator)
    single_particle_state_to_index_map = Dict(
        state => i
        for (i, state) in enumerate(single_particle_state_generator)
    )

    return Finite_Single_Particle_Hilbert_Space{D}(
        statistics,
        D,
        dof_name,
        dof_ndof,
        single_particle_state_generator,
        nstate,
        single_particle_state_to_index_map,
    )
end


# TODO: `Infinite_Single_Particle_Hilbert_Space` for infinite-dimensional single-particle Hilbert space





#=====================================================================================================================================#
#==================================================== Many-Particle Hilbert Space ====================================================#
#=====================================================================================================================================#

"""
Finite-dimensional Many-Particle Hilbert Space `Finite_Many_Particle_Hilbert_Space <: Finite_Hilbert_Space`
---
- Fields:
    - `statistics::Quantum_Statistics`: quantum statistics of the single-particle states
    - `nparticle::Integer`: number of particles
    - `state_generator::Base.Generator`: lazy generator for the list of many-body states. Here each many-body state is a tuple of single-particle states. Note: we should never try to explicitly store the Full list of many-body states, since it can cost formidablely large amount of memory! On the other hand, it is suitable for the on-the-fly algorithm for ED calculation
    - `state_to_index_map::Dict{<:Many_Particle_State,BigInt}`: dict to get the index of each many-body state
    - `nstate::BigInt`: number of many-body states (it can be a large so it's better to store it explicitly)
"""
struct Finite_Many_Particle_Hilbert_Space <: Finite_Hilbert_Space
    statistics::Quantum_Statistics # quantum statistics of the single-particle states
    nparticle::Integer # number of particles
    state_generator::Base.Generator # lazy generator for the list of many-body states. Here each many-body state is a tuple of single-particle states. Note: we should never try to explicitly store the Full list of many-body states, since it can cost formidablely large amount of memory! On the other hand, it is suitable for the on-the-fly algorithm for ED calculation
    nstate::BigInt # number of many-body states (it can be a large so it's better to store it explicitly)
    state_to_index_map::Dict{<:Many_Particle_State,BigInt} # dict to get the index of each many-body state
end

"""
Struct Initialization of `Finite_Many_Particle_Hilbert_Space`
---
- Named_Args:
    - `single_particle_Hilbert_space::Finite_Single_Particle_Hilbert_Space`: single-particle Hilbert space
    - `nfermion::Integer`: number of fermions
"""
function initialize_finite_fermionic_many_particle_Hilbert_space(;
    single_particle_Hilbert_space::Finite_Single_Particle_Hilbert_Space,
    nfermion::Integer
)
    statistics = single_particle_Hilbert_space.statistics
    @assert statistics == Fermionic() "Check Single-particle Statistics: method `initialize_finite_fermionic_many_particle_Hilbert_space()` is defined for Fermionic particles!"
    many_particle_state_generator = Combinatorics.combinations(collect(single_particle_Hilbert_space.state_generator), nfermion)

    nstate::BigInt = binomial(BigInt(length(single_particle_Hilbert_space.nstate)), nfermion)
    many_particle_state_to_index_map = Dict(
        Many_Particle_State(many_particle_state) => BigInt(i)
        for (i, many_particle_state) in enumerate(many_particle_state_generator)
    )
    # many_particle_state_to_index_map = Dict{Many_Particle_State,BigInt}()

    return Finite_Many_Particle_Hilbert_Space(
        statistics,
        nfermion,
        many_particle_state_generator,
        nstate,
        many_particle_state_to_index_map
    )
end


# TODO: `Infinite_Many_Particle_Hilbert_Space` for infinite-dimensional many-particle Hilbert space
