
using Plots, Printf, LinearAlgebra, GeoData, NCDatasets, JLD


# enable plotting & saving by default
if !@isdefined do_visu; do_visu = true end
if !@isdefined do_save; do_save = true end

# finite difference stencil operation support functions
@views av(A) = 0.5.*(A[1:end-1].+A[2:end]) # average x-dir
@views inn(A)   = A[2:end-1] # inner points

function mass_balance_constants(xc, yc)
    b_max    = 0.15            # max. Mass balance rate
    lat_min, lat_max = 60, 80
    Xc, Yc   = [Float32(x) for x=xc,y=yc], [Float32(y) for x=xc,y=yc]
    Yc2      = Yc .- minimum(Yc); Yc2 .= Yc2/maximum(Yc2)
    grad_b   = (1.3517 .- 0.014158.*(lat_min.+Yc2*(lat_max-lat_min)))./100.0.*0.91 # Mass Bal. gradient, from doi: 10.1017/jog.2016.75
    z_ELA    = 1300.0 .- Yc2*300.0                                 # Educated guess for ELA altitude
    return grad_b, z_ELA, b_max
end

@views function diffusion_1D(; do_visu=true)

    # Physics
    lx     = 96      # domain size
    #D      = 1.0        # diffusion coefficient
    #ttot   = 0.6        # total simulation time
    #dt     = 0.1        # physical time step
    # physics
    s2y      = 3600*24*365.25  # seconds to years
    rho_i    = 910.0           # ice density
    g        = 9.81            # gravity acceleration
    npow     = 3.0             # Glen's power law exponent
    a0       = 1.5e-24         # Glen's law enhancement term
    dt       = 0.1        # physical time step
    # Numerics
    # Derived numerics
    # Valeurs arbitaires
        #nx     = size(Zbed,1) # numerical grid resolution #A changer avec les donnees extrait de Zbed 
    #ny     = size(Zbed,2)
    #@assert (nx, ny) == size(Zbed) == size(Hice) == size(Mask) "Size doesn't match"
    itMax    = 1e5             # number of iteration (max)
    nout     = 1 #200         # error check frequency
    tolnl    = 1e-6            # nonlinear tolerance
    epsi     = 1e-4            # small number
    damp     = 0.85            # convergence accelerator (this is a tuning parameter, dependent on e.g. grid resolution)
    dtausc   = 1.0/3.0         # iterative dtau scaling
    # derived physics
    a      = 2.0*a0/(npow+2)*(rho_i*g)^npow*s2y  #Viscosité de la glace

    #dtau   = (1.0/(dx^2/D/2.1) + 1.0/dt)^-1 # iterative "timestep" #A changer dans chaque boucle 
    #xc     = x
    nx = 96
    # Array allocation
    qH     = zeros(nx-2)
    dHdtau = zeros(nx-3)
    dtau   = zeros(nx-3) #A initialiser à chaque tour de boucle 
    ResH   = zeros(nx-3)
    Err    = zeros(nx)
    dSdx   = zeros(nx-1)
    M      = ones(nx)
    B      = zeros(nx) #Initialisation d'un tableau de 0 de taille nx 
    H      = zeros(nx)
    D      = zeros(nx-1) #Initialisation d'un tableau de 0 de taile nx pour le coeff de diffusion pour chaque point 
    # Initial condition
    # On charge notre fichier GeoData...
    # Puis on extrait les valeurs dans les variables...
    var = load("BedMachineGreenland_96_176.jld"); # ultra low res data
    Zbed = var["Zbed"]
    xc = var["xc"]
    yc = var["yc"]
    Hice = var["Hice"]
    Mask1 = var["Mask"]
    dx =15600;
    cfl    = dx^2/4.1
    grad_b1, z_ELA1, b_max = mass_balance_constants(xc,yc)
    grad_b=grad_b1[:,80]
    z_ELA = z_ELA1[:,80]
    S = zeros(nx) # Épaisseur de de la glace
    B1 = Zbed.data # Niveau de la mer
    H1 = Hice.data # Hauteur de la glace
    B = B1[:,80]
    H = H1[:,80]
    Mask=Mask1[:,80]
    S .= B .+ H
    #display(plot(xc,S))
    #display(plot!(xc,B))
    #; it = 0; ittot = 0
    # iteration loop
    println(" starting iteration loop:")
    iter = 1; err = 2*tolnl
    # Physical time loop
    #while iter<itMax
        #iter = 0; 
        # Pseudo-transient iteration
        while err>tolnl && iter<itMax #Quelle autre  condition
            Err .= H
            M     .= min.(grad_b.*(S .- z_ELA), b_max) 
            dSdx .= diff(S)/dx
            D     .= a*av(H).^(npow+2) .*dSdx.^(npow-1) #Devient une variable 
            qH         .= .-av(D).*diff(S[2:end])/dx  # flux -> négatif lorsqu'on grimpe la courbe
            ResH  .= .-(diff(qH)/dx) .+ inn(M[1:end-1]) #a quoi sert le m
            dtau  .= dtausc*min.(10.0, cfl./(epsi .+ av(D[1:end-1])))
            dHdtau     .= ResH + damp*dHdtau         # damped rate of change
            H[2:end-2] .= max.(0.0,H[2:end-2] .+ dtau.*dHdtau)   # update rule, sets the BC as H[1]=H[end]=0
            # apply mask (a very poor-man's calving law)
            H[Mask.==0] .= 0.0
            # update surface
            S     .= B .+ H
            # error check
            if mod(iter, nout)==0
                @printf("diff(S) = %f ,dSdx = %f, S = %f, D = %f, qH = %f, ResH = %f, H = %f, dtau = %f", (S[50]-S[51])/dx, dSdx[50], S[50], D[50], qH[50], ResH[50], H[50], dHdtau[50])
                Err .= Err .- H
                err = norm(Err)/length(Err)
                @printf(" Err = %f, iter = %d, error = %1.2e \n", Err[50], iter, err)
                if isnan(err)
                    error("""NaNs encountered.  Try a combination of:
                                decreasing `damp` and/or `dtausc`, more smoothing steps""")
                end
            end
            iter += 1
        end
        
        #ittot += iter; it += 1; t += dt
        #iter += 1
        #t += dt  #Necessaire ? 
   # end
    # compute velocities
    Vx .= -D./(av(H) .+ epsi).*av(dSdx)
    #Vy .= -D./(av(H) .+ epsi).*av_xa(dSdy)
    # return as GeoArrays
    return  as_geoarray(H,  Zbed, name=:thickness),
            as_geoarray(S,  Zbed, name=:surface),
            as_geoarray(M,  Zbed, name=:smb),
            as_geoarray(Vx, Zbed, name=:vel_x, staggerd=true),
            as_geoarray(Vy, Zbed, name=:vel_y, staggerd=true)
end

diffusion_1D(; do_visu=do_visu);

# ------------------------------------------------------------------------------

