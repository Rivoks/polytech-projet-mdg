
using Plots, Printf, LinearAlgebra, GeoData, NCDatasets, JLD


# enable plotting & saving by default
if !@isdefined do_visu; do_visu = true end
if !@isdefined do_save; do_save = true end

# finite difference stencil operation support functions
@views av(A) = 0.5.*(A[1:end-1].+A[2:end]) # moyenne
@views inn(A)   = A[2:end-1] # inner points

function mass_balance_constants(xc, yc) # Permet de récupérer toutes les valeurs nécessaires pour l'ablation
    b_max    = 0.15            # max. Mass balance rate
    lat_min, lat_max = 60, 80
    Xc, Yc   = [Float32(x) for x=xc,y=yc], [Float32(y) for x=xc,y=yc]
    Yc2      = Yc .- minimum(Yc); Yc2 .= Yc2/maximum(Yc2)
    grad_b   = (1.3517 .- 0.014158.*(lat_min.+Yc2*(lat_max-lat_min)))./100.0.*0.91 # Mass Bal. gradient, from doi: 10.1017/jog.2016.75
    z_ELA    = 1300.0 .- Yc2*300.0                                 # Educated guess for ELA altitude
    return grad_b, z_ELA, b_max
end

function as_geoarray(A, template::AbstractGeoArray; name=nothing, staggerd=false)
    @assert length(size(A))==2
    dd = if staggerd
        x, y = dims(template)
        (X( 0.5*(val(x)[2:end] + val(x)[1:end-1]), mode=x.mode, metadata=x.metadata),
         Y( 0.5*(val(y)[2:end] + val(y)[1:end-1]), mode=y.mode, metadata=y.metadata))
    else
        dims(template)
    end
    return GeoArray(A, dims=dd, name=name)
end

@views function diffusion_1D(; do_visu=true)

    # Physics
    # Domaine de taille 96,176
    s2y      = 3600*24*365.25  # seconds to years
    rho_i    = 910.0           # ice density
    g        = 9.81            # gravity acceleration
    npow     = 3.0             # Glen's power law exponent
    a0       = 1.5e-24         # Glen's law enhancement term
    # Numerics
    itMax    = 1e5             # number of iteration (max)
    nout     = 50 #200         # error check frequency
    tolnl    = 1e-6            # nonlinear tolerance
    epsi     = 1e-4            # small number
    damp     = 0.65            # Pour converger plus vite (on l'a arbitrairement mis à 0.65)
    dtausc   = 0.75/3.0        # iterative dtau scaling (on l'a arbitairement mis à 0.75/3, peut changer)
    # derived physics
    cfl    = 0                 # Confition de stzbilite
    a      = 2.0*a0/(npow+2)*(rho_i*g)^npow*s2y  #Viscosité de la glace
    nx = 96 # Nombre de valeurs en x
    # Array allocation (Chaque tableau à une taille dépendante du calcul qu'on y fait)
    Vx     = zeros(nx-1)
    qH     = zeros(nx-2) # Flux
    dHdt   = zeros(nx-3)
    dtau   = zeros(nx-3)
    ResH   = zeros(nx-3)
    Err    = zeros(nx) # Erreur à chaque pas de temps
    dSdx   = zeros(nx-1)
    M      = ones(nx) # Ablation
    B      = zeros(nx) # Niveau du sol 
    H      = zeros(nx) # Epaisseur de la glace
    D      = zeros(nx-1) # Coefficient de diffusion
    S      = zeros(nx) # Hauteur de la glace
   
   
    # Initial condition
    # On charge notre fichier GeoData...
    # Puis on extrait les valeurs dans les variables...
    # On extrait toutes les valeurs nécessaires depuis helpers.jl
    var = load("BedMachineGreenland_96_176.jld"); # ultra low res data
    Zbed = var["Zbed"]
    xc = var["xc"]
    yc = var["yc"]
    Hice = var["Hice"]
    Mask1 = var["Mask"]
    dx =15600;
    grad_b1, z_ELA1, b_max = mass_balance_constants(xc,yc)
    cfl    = dx^2/4.1 
    B1 = Zbed.data # Niveau de la terre
    H1 = Hice.data # Epaisseur de la glace

    # On ne récupère qu'une colonne qu'on étudiera ensuite (Arbitrairement choisi la colonne 80)
    B = B1[:,80]
    H = H1[:,80]
    H0 = copy(H)
    grad_b=grad_b1[:,80] # Help  gradient du bilan massique
    z_ELA = z_ELA1[:,80] # Help latitude de la ligne d'équilibre
    Mask=Mask1[:,80] # Booléen qui nous permettra d'appliquer le masque pour le niveau de la mer
    S .= B .+ H # Evaluation de la hauteur de la surface (Base + Glace)
  
    # iteration loop
    println(" starting iteration loop:")
    iter = 1; err = 2*tolnl
    while err>tolnl && iter<itMax # On sort lorsque l'erreur est en dessous de la tolérance (ou trop d'itérations)
        Err .= H # On remet à jour l'erreur
        M     .= min.(grad_b.*(S .- z_ELA), b_max) # Abblation
        dSdx .= diff(S)/dx 
        D     .= a*av(H).^(npow+2) .*dSdx.^(npow-1) # Coefficient de diffusion NON LINEAIRE
        qH         .= .-av(D).*diff(S[1:end-1])/dx  # flux -> négatif lorsqu'on grimpe la courbe
        ResH  .= .-(diff(qH)/dx .+ inn(M[1:end-1])) 
        dtau  .= dtausc*min.(10.0, cfl./(epsi .+ av(D[2:end])))
        dHdt     .= ResH + damp.*dHdt         # Ajout du damp pour la convergence
        H[2:end-2] .= max.(0.0,H[2:end-2] .+ dtau.*dHdt)  # Permet de mettre le niveau de la mer à 0
        # apply mask (a very poor-man's calving law)
        H[Mask.==0] .= 0.0
        # update surface
        S     .= B .+ H # Mise à jour de la surface
        # error check
        if mod(iter, nout)==0 # Calcul de l'erreur tous les nout temps
            @printf("dSdx = %f, S = %f, D = %f, qH = %f, ResH = %f, H = %f, dtau = %f", dSdx[50], S[50], D[50], qH[50], ResH[50], H[50], dHdt[50])
            Err .= Err .- H
            err = norm(Err)/length(Err) # Length = 96
            @printf(" Err = %f, iter = %d, error = %1.2e \n", Err[50], iter, err)
            if isnan(err)
                error("""NaNs encountered.  Try a combination of:
                            decreasing `damp` and/or `dtausc`, more smoothing steps""")
            end
        end
        iter += 1
    end

    # compute velocities
    Vx .= -D./(av(H) .+ epsi).*dSdx
    # return as GeoArrays
    if do_visu
        display(plot(xc, H0, linewidth=3, label="H initial")); display(plot!(xc, H, legend=true, framestyle=:box, linewidth=3,label="H apres diffusion", xlabel="lx", ylabel="H", title="iceflow 1D (nt=$nout, iters=$iter)"))
        display(plot!(xc, B, linewidth=3, label="lit rocheux"));
    end
end

diffusion_1D(; do_visu = true)




# ------------------------------------------------------------------------------