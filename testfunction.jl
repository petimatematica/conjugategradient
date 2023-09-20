

# Extended Powell singular function (the dimension of this function must be a multiple of four)

function singx(x)
    s = 0
    for _ in 1:length(x)
        if i%4 == 1
            s += (x[i]+10x[i+1])^2
        elseif i%4 == 2
            s += 5*(x[i+1]-x[i+2])^2
        elseif i%4 == 3
            s += (x[i-1]-2*x[i])^4
        else
            s += 10*(x[i-3]-x[i])^4
        end
    end
    return s
end

function gradsingx(x)
    v = Float64[]
    for _ in 1:length(x)
        if i%4 == 1
            push!(v, 2*(x[i]+10*x[i+1])+40*(x[i]-x[i+3])^3)
        elseif i%4 == 2
            push!(v, 20*(x[i-1]+10*x[i])+4*(x[i]-2*x[i+1])^3)
        elseif i%4 == 3
            push!(v, 10*(x[i]-x[i+1])-8*(x[i-1]-2*x[i])^3)
        else
            push!(v, -10*(x[i-1]-x[i])-40*(x[i-3]-x[i])^3)
        end
    end
    return v
end
