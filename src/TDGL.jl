module TDGL
    using Plots
    
    struct parameters
        L::Array{Float64,1}
        N::Array{Int64,1}
        dr::Array{Float64,1}               
        Hz::Float64
        δt :: Float64
        Maxstep ::Int64
        dirname ::String
        κ::Float64
        calc_A::Bool
    end

    function init_system(L,N,Hz;δt = 0.01,Maxstep=1000,name="output",κ=2.0,calc_A=false)
        dr = L ./ (N .-1)
        return parameters(L,N,dr,Hz,δt,Maxstep,name,κ,calc_A)
    end

    function tdGLsimulation(Δ0,param,T,U)
        Nx = param.N[1]
        Ny = param.N[2]
        δt = param.δt
        name = param.dirname
        Δ = copy(Δ0)        
        dΔdt = zeros(ComplexF64,Nx,Ny)
        dUdt = zeros(Float64,2,Nx,Ny)
        count = 0
        if param.calc_A
            U[:,2:Nx+1,2:Ny+1] .= 1
        end
        println("TDGL simulation starts!")
        println("Num | ","\t"," Δ at the center | ","\t","|Δ| at the center | ","\t","dΔdt")
        maps = @animate  for i=1:param.Maxstep            
            calc_dΔdt!(dΔdt,U,Δ,param,T)   
            if param.calc_A                  
                calc_dUdt!(dUdt,U,Δ,param,T)
            end
            Δ += δt*dΔdt   
            if param.calc_A
                U[:,2:Nx+1,2:Ny+1] = U[:,2:Nx+1,2:Ny+1].*exp.(-im*δt.*dUdt[:,:,:]) #-δt*dUdt 
            end
              
            p = heatmap(1:Nx, 1:Ny, abs.(Δ),aspect_ratio=:equal);

            if i % 40 ==1
                count += 1                
                savefig(p,"./"*name*"/amp_"*string(count))     
                println(i,"\t",Δ[div(Nx,2),div(Ny,2)],"\t",abs(Δ[div(Nx,2),div(Ny,2)]),"\t",dΔdt[div(Nx,2),div(Ny,2)])             
            end
                   
        end every 40 
        gif(maps, "./"*name*"/GL.gif", fps = 15)  
        return Δ,U       
    end

    function init_U(param)
        Nx = param.N[1]
        Ny = param.N[2]
        U = zeros(ComplexF64,2,Nx+2,Ny+2) #which includes boundaries       
        dr = param.dr
        for ix=0:Nx+1
            for iy=0:Ny+1
                dx = 1
                jx = ix + dx/2
                jy = iy 
                r = ixiy2r(jx,jy,param)              
                  U[1,ix+1,iy+1] = exp(im*dr[1]*Ax(r,param))
                dy = 1
                jx = ix 
                jy = iy +dy/2
                r = ixiy2r(jx,jy,param)
                U[2,ix+1,iy+1] = exp(im*dr[2]*Ay(r,param))
            end
        end
        return U
    end

    Ax(r,param) = -param.Hz*r[2]/2
    Ay(r,param) = param.Hz*r[1]/2

    function ixiy2r(ix,iy,param)
        rx = (ix-1)*param.dr[1]-param.L[1]/2
        ry = (iy-1)*param.dr[2]-param.L[2]/2
        return rx,ry
    end

    function get_Ux(ix,iy,U)
        return get_U(1,ix,iy,U)
    end

    function get_Uy(ix,iy,U)
        return get_U(2,ix,iy,U)
    end

    function get_U(xy,ix,iy,U)        
        return U[xy,ix+1,iy+1]
    end

    function plaq(ix,iy,U)
        U1 = get_Ux(ix,iy,U)
        U4 = get_Uy(ix,iy,U)

        jx = ix + 1
        jy = iy        
        U2 = get_Uy(jx,jy,U)

        jx = ix
        jy = iy+1
        U3 = get_Ux(jx,jy,U)

        Lzij = U1*U2*conj(U3)*conj(U4)
        return Lzij
    end

    function calc_dΔdt!(dΔdt,U,Δ,param,T)
        Nx = param.N[1]
        Ny = param.N[2]
        dr = param.dr
        for ix=1:Nx
            for iy=1:Ny
                dΔdt[ix,iy] = (1-T)*(1-abs(Δ[ix,iy])^2)*Δ[ix,iy]                
                for dx =-1:1
                    jx = ix + dx
                    jy = iy
                    if jx > Nx || jx < 1
                        t = 0
                    else
                        if dx == 0
                            t = -2*Δ[ix,iy]/dr[1]^2
                        else
                            r = ixiy2r(jx,jy,param)
                            U1 = get_Ux(jx,jy,U)
                            t = Δ[jx,jy]*U1^(dx)/dr[1]^2                            
                        end
                    end
                    dΔdt[ix,iy] += t
                end

                for dy =-1:1
                    jx = ix 
                    jy = iy + dy
                    if jy > Ny || jy < 1
                        t = 0
                    else
                        if dy == 0
                            t = -2*Δ[ix,iy]/dr[2]^2
                        else
                            r = ixiy2r(jx,jy,param)
                            U1 = get_Uy(jx,jy,U)
                            t = Δ[jx,jy]*U1^(dy)/dr[2]^2                            
                        end
                    end
                    dΔdt[ix,iy] += t
                end                
            end
        end        
    end

    function calc_dUdt!(dUdt,U,Δ,param,T)
        Nx = param.N[1]
        Ny = param.N[2]
        dr = param.dr
        for ix=1:Nx
            for iy=1:Ny
                Lzij = plaq(ix,iy,U)
                dx = 1
                jx = ix + 1
                jy = iy 
                if jx <= Nx        
                           
                    U1 = get_Ux(ix,iy,U)                        
                    Jx = imag(conj(Δ[ix,iy])*U1*Δ[jx,jy])
                    #Jx = 0.0 
                else
                    Jx = 0.0
                end
                Lzijm = plaq(ix,iy-1,U)
                
                dUdt[1,ix,iy] = real((1-T)*Jx+((param.κ^2/dr[2]^2)*(Lzij*conj(Lzijm)-1)/im))
                
                dy = 1
                jx = ix 
                jy = iy  + dy
                if jy <= Ny 
                      
                    U1 = get_Uy(ix,iy,U)                                   
                    Jy = imag(conj(Δ[ix,iy])*U1*Δ[jx,jy])
                    #Jy = 0.0   
                else
                    Jy = 0.0
                end                
                Lzimj = plaq(ix-1,iy,U)
                dUdt[2,ix,iy] = real((1-T)*Jy+((param.κ^2/dr[1]^2)*(conj(Lzij)*Lzimj-1)/im))                                
            end
        end        
    end

    function test2()
        L = [25.0,25.0]
        N = [51,51]
        Hz = 0.1
        param = init_system(L,N,Hz,Maxstep = 40000)
        Δ0 = 0.9*ones(ComplexF64,N[1],N[2])
#        Δ0[5,5] += im
        T = 0.2
        U = init_U(param)
        Δ,U = tdGLsimulation(Δ0,param,T,U)
    end
end # module
