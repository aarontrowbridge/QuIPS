# Quantum Virtual Machine Library

module QVM

export QCircuit, run!, step!, reset!, operate!

using Operators

const C = Complex{Float32}

mutable struct QCircuit
    wfn::Vector{C}
    ops::Vector{Operator}
    out::BitVector
    pos::Int
    N::Int

    function QCircuit(quip::Vector, N::Int)
        wfn = zeros(C, 2^N); wfn[1] = 1
        ops = compile(quip)
        out = BitVector(undef, N)
        new(wfn, ops, out, 0, N)
    end
end

function compile(quip::Vector, verbose=true)
    println("\ncompiling quip to circuit\n")
    Vs = []
    for tag in quip
        if tag[1] == :MEASURE || tag[1] == :M
            push!(Vs, Measurement(tag[2]))
        elseif String(tag[1])[1] == 'C'
            push!(Vs, Control(tag))
        else
            push!(Vs, Gate(tag))
        end
        if verbose
            V = Vs[end]
            if typeof(V) == Control
                println("  compiling ", V.name, " ", V.indx...)
            elseif typeof(V) == Measurement
                println("  compiling M ", V.k)
            else
                println("  compiling ", V.name, " ", V.k)
            end
        end
    end
    println("\ncompiled\n")
    Vs
end

function run!(QC::QCircuit, verbose=true)
    if verbose
        println("starting run now\n")
    end
    for V in QC.ops
        QC.pos += 1
        evolve!(QC, V)
        if verbose
            if typeof(V) == Control
                println(V.name, " ", V.indx...)
            elseif typeof(V) == Measurement
                println("M ", V.k, " -> ", QC.out[V.k])
            else
                println(V.name, " ", V.k)
            end
        end
    end
end

function step!(QC::QCircuit)
    if QC.pos < length(QC.ops)
        QC.pos += 1
        V = QC.ops[QC.pos]
        evolve!(QC, V)
    end
end

function reset!(QC::QCircuit)
    QC.wfn = zero.(QC.wfn); QC.wfn[1] = 1
    QC.out = BitVector(undef, QC.N)
    QC.pos = 0
end

function evolve!(QC::QCircuit, V::Operator)
    if typeof(V) == Measurement
        ψ, outcome = V(QC.wfn, QC.N)
        QC.out[V.k] = outcome
        QC.wfn .= ψ
    else
        ψ = V(QC.wfn, QC.N)
        QC.wfn .= ψ
    end
end

function operate!(QC::QCircuit, Vs::Vector{T} where {T<:Operator})
    for (i, V) in enumerate(Vs)
        insert!(QC.ops, QC.pos + i, V)
    end
    step!(QC)
end

end
