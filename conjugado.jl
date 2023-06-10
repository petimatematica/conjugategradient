# 
# Strong Wolfe's Linesearch
# α comprimento de passo inicial
#β= c1 e σ =c2 (nocedal)
#y0=f(x)
#g0= gradiente 
# NaN não é um número
function wolfe(f, ∇, x, d; α=1, β=1e-4, σ=0.1)
    y0, g0, y_prev, α_prev = f(x), ∇(x)⋅d, NaN, 0
    αlo, αhi = NaN, NaN

    #bracket phase   #COMPRIMENTO DE PASSOS DESEJAVÉIS
    while true
        y = f(x + α*d)
        if y > y0 + β*α*g0 || (!isnan(y_prev) && y ≥ y_prev)
            αlo, αhi = α_prev, α
            break
        end
        g = ∇(x + α*d)⋅d
        if abs(g) ≤ -σ*g0
            return α
        elseif g ≥ 0
            αlo, αhi = α, α_prev
            break
        end
        y_prev, α_prev, α = y, α, 2α
    end

    # zoom phase
    ylo = f(x + αlo*d)
    while true
        α = (αlo + αhi)/2
        y = f(x + α*d)
        if y > y0 + β*α*g0 || y ≥ ylo
            αhi = α
        else
            g = ∇(x + α*d)⋅d
            if abs(g) ≤ -σ*g0
                return α

            elseif g*(αhi - αlo) ≥ 0
                αhi = αlo
            end
            αlo = α
        end
    end
end      # wolfe atá aqui




#   Armijo's Linesearch (for Gradient Method)
#
function armijo(f, ∇f, x, d; γ=0.9, c1=1e-4, maxiter=1000)
    t = 1;
    iter = 0;
    gradfx = ∇f(x); 
    k = gradfx'*d
    k = k[1]
    while (f(x+t*d) > f(x) + c1*t*k) && (iter < maxiter)
        t = γ*t
        iter = iter + 1
    end

    if iter > maxiter
        return("número máximo de iteradas excedido")
    else
        return t
    end
end


#
#   Goldstein's Linesearch (for Gradient Method)
#
function goldstein(f, ∇, x, d; α=1, β=1e-4, σ=0.1)
    y0, g0, y_prev, α_prev = f(x), ∇(x)⋅d, NaN, 0
    αlo, αhi = NaN, NaN

    #bracket phase   #COMPRIMENTO DE PASSOS DESEJAVÉIS
    while true
        y = f(x + α*d)
        if y0 + (1-β)*α*g0 > y || (!isnan(y_prev) && y ≥ y_prev)
            αlo, αhi = α_prev, α
            break
        end
        g = ∇(x + α*d)⋅d
        if y ≤ y0 + β*α*g0
            return α
        elseif g ≥ 0
            αlo, αhi = α, α_prev
            break
        end
        y_prev, α_prev, α = y, α, 2α
    end

    # zoom phase
    ylo = f(x + αlo*d)
    while true
        α = (αlo + αhi)/2
        y = f(x + α*d)
        if y0 + (1-β)*α*g0 > y || y0 + (1-β)*α*g0 ≥ ylo
            αhi = α
        else
            g = ∇(x + α*d)⋅d
            if y ≤ y0 + β*α*g0
                return α

            elseif g*(αhi - αlo) ≥ 0
                αhi = αlo
            end
            αlo = α
        end
    end
end


function conjugado(x, f, ∇f, n,lineseach,method; maxiter = 1000, ϵ = 1e-6)
    d = -∇f(x)
    gradfx = -d
    iter = 0
    while norm(∇f(x)) > ϵ && iter < maxiter 
        t = lineseach(f, ∇f, x, d)
        x = x + t*d
        newgradfx = ∇f(x)
        if method == 0
            βk=0
        else
            if (iter+1)%n != 0 
                βk = (newgradfx'⋅newgradfx)/(gradfx'⋅gradfx)  
            else
                βk = 0
            end
        end
        d = -newgradfx + βk*d
        iter = iter + 1
        gradfx = newgradfx
    end
    if iter < maxiter
        return x, iter
    else
        return ("não deu certo")
    end
end

