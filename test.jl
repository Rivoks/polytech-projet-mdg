
using Plots, Printf, LinearAlgebra, GeoData, NCDatasets, JLD


# enable plotting & saving by default
if !@isdefined do_visu; do_visu = true end
if !@isdefined do_save; do_save = true end

# finite difference stencil operation support functions
@views av(A) = 0.5.*(A[1:end-1].+A[2:end]) # average x-dir
@views inn(A)   = A[2:end-1] # inner points

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
    nout     = 50 #200         # error check frequency
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
    Err    = zeros(nx  )
    dSdx   = zeros(nx-1)
    M      = zeros(nx-1)
    B      = zeros(nx) #Initialisation d'un tableau de 0 de taille nx 
    H      = zeros(nx)
    D      = zeros(nx-1) #Initialisation d'un tableau de 0 de taile nx pour le coeff de diffusion pour chaque point 
    # Initial condition
    # On charge notre fichier GeoData...
    # Puis on extrait les valeurs dans les variables...
    var = load("BedMachineGreenland_96_176.jld"); # ultra low res data
    Zbed = var["Zbed"]
    xc = var["xc"]
    Hice = var["Hice"]
    dx =15600;
    S = zeros(96) # Épaisseur de de la glace
    B1 = Zbed.data # Niveau de la mer
    H1 = Hice.data # Hauteur de la glace
    B = B1[:,80]
    H = H1[:,80]
    S .= B .+ H
    plot(xc,S)
    plot!(xc,B1)
    t = 0.0
    #; it = 0; ittot = 0
    # iteration loop
    println(" starting iteration loop:")
    iter = 1; err = 2*tolnl
    # Physical time loop
    while iter<itMax
        #iter = 0; 
        # Pseudo-transient iteration
        while err>tolnl  #Quelle autre  condition 
            dSdx .= diff(S)/dx
            D     .= a*av(H).^(npow+2) .*dSdx.^(npow-1) #Devient une variable 
            qH         .= .-av(D).*diff(S[2:end])/dx  # flux
            ResH  .= .-(diff(qH)/dx) .+ inn(M) #a quoi sert le m
            #ResH       .= -(H[2:end-1] - Hold[2:end-1])/dt - diff(qH)/dx # residual of the PDE
            dHdtau     .= ResH + damp*dHdtau         # damped rate of change
            H[2:end-2] .= H[2:end-2] + dtau.*dHdtau   # update rule, sets the BC as H[1]=H[end]=0

            # error check
            if mod(iter, nout)==0
                Err .= Err .- H
                err = norm(Err)/length(Err)
                @printf(" iter = %d, error = %1.2e \n", iter, err)
                if isnan(err)
                    error("""NaNs encountered.  Try a combination of:
                                decreasing `damp` and/or `dtausc`, more smoothing steps""")
                end
            end
        end
        
        #ittot += iter; it += 1; t += dt
        iter += 1
        #t += dt  #Necessaire ? 
    end
    # compute velocities
    Vx .= -D./(av(H) .+ epsi).*av_ya(dSdx)
    Vy .= -D./(av(H) .+ epsi).*av_xa(dSdy)
    # return as GeoArrays
    return  as_geoarray(H,  Zbed, name=:thickness),
            as_geoarray(S,  Zbed, name=:surface),
            as_geoarray(M,  Zbed, name=:smb),
            as_geoarray(Vx, Zbed, name=:vel_x, staggerd=true),
            as_geoarray(Vy, Zbed, name=:vel_y, staggerd=true)
end

diffusion_1D(; do_visu=do_visu);

# ------------------------------------------------------------------------------

