MLStyle.@data Quantum_Statistics begin
    Fermionic()
    Bosonic()
end


# add comparison between `Quantum_Statistics`: we always want to write bosonic components before fermionic components"
Base.:(==)(statistics1::Quantum_Statistics, statistics2::Quantum_Statistics) = begin
    if statistics1 isa Bosonic && statistics2 isa Bosonic ||
       statistics1 isa Fermionic && statistics2 isa Fermionic
        return true
    else
        return false
    end
end
Base.:<(statistics1::Quantum_Statistics, statistics2::Quantum_Statistics) = begin
    if !(statistics1 == statistics2)
        if statistics1 isa Bosonic
            return true
        else
            return false
        end
    else # statistics1 == statistics2
        return false
    end
end
Base.:>(statistics1::Quantum_Statistics, statistics2::Quantum_Statistics) = begin
    if !(statistics1 == statistics2)
        if statistics1 isa Fermionic
            return true
        else
            return false
        end
    else # statistics1 == statistics2
        return false
    end
end