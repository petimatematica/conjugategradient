function cgzhang(x,f,gf,ϵ,ϵ1;maxiter=20)
    ierror = 0
    iter = 0
    gradfxp = 0.0
    d = 0.0
    while true
    
        gradfx = gf(x)
        normgfx = norm(gradfx,2)

        println("$iter   $normgfx   $(f(x))")

        if normgfx < ϵ
            return(x,iter,ierror)
        end

        if iter > 0
            y = gradfx - gradfxp
            dirtest = (d' * y)[1] - ϵ1 * norm(d,2) * normgfx
            if dirtest < 0.0
                d = - gradfx
            else
                β = norm(gradfx)^2 / (d' * y)
                d = -gradfx + β * d
            end 
        else
            d = -gradfx
        end

        (α,ierror) = armijo(f, gradfx, x, d)

        iter += 1

        if iter > maxiter
            ierror = 1
            return(x,iter,ierror)
        end

        x = x + α * d

        gradfxp = gradfx

    end


end



function armijo(f, gradfx, x, d; γ=0.5, c1=0.1, maxk=1000)
    t = 1.0;
    k = 0;
    #gradfx = ∇f(x); 
    gd = gradfx'*d
    gd = c1 * gd[1]

    ierror = 0

    while (f(x+t*d) > f(x) + t*gd) && (k < maxk)
        t = γ*t
        k = k + 1
    end

    if k < maxk
        return t, ierror
    else
       # ierror = 2
        return t, 2
    end
end