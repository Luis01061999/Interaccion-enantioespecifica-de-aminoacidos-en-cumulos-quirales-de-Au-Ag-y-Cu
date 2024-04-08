using LinearAlgebra, Statistics, Random#, Plots

## Funciones para abrir, escribir y extrar posiciones de archivos .xyz (y formato FHI-aims)

"""
The function open_xyz_file() recive a path of XYZ file and return an array átoms symbols, array atóm coords and comment 
of the XYZ file
"""
function open_xyz_file(path::AbstractString)
    file = open(path)
    num_atoms = parse(Int, readline(file))
     comment = readline(file)
    symbols = Vector{String}(undef, num_atoms)
    coords = zeros(num_atoms,3)
    for i in 1:num_atoms
        line = readline(file)
        data = split(line)
        symbols[i] = data[1]
        coords[i, :] = [parse(Float64, d) for d in data[2:4]]
    end
    close(file)
    return symbols, coords, comment
end

"""
The function write_xyz recive a name 
Define a function to write an XYZ file with element symbols
"""
function write_xyz(filename, coords, symbols)
    # Get the number of atoms and dimensionality of the coordinates
    natoms, ndim = size(coords)

    # Open the output file for writing
    output = open(filename, "w")

    # Write the header to the output file
    write(output, "$natoms\n")
    write(output, "XYZ file written with Julia\n")

    # Write the coordinates and element symbols to the output file
    for i in 1:natoms
        write(output, "$(symbols[i])")
        for j in 1:ndim
            write(output, " $(coords[i,j])")
        end
        write(output, "\n")
    end

    # Close the output file
    close(output)
end

"""
La función write_multiple_xyz
"""
function write_multiple_xyz(filename, coords, symbols)
    # Open the output file for writing
    output = open(filename, "w")
    write(output, "total estructures $(length(coords))\n")
    # Write each structure to the output file
    for i in 1:length(coords)
        structure = coords[i]
        natoms, ndim = size(structure)

        # Write the header for each structure
        write(output, "estructura $i\n")
        write(output, "$natoms\n")
        write(output, "XYZ file written with Julia\n")

        # Write the coordinates and element symbols for each atom
        for j in 1:natoms
            write(output, "$(symbols[j])")
            for k in 1:ndim
                write(output, " $(structure[j, k])")
            end
            write(output, "\n")
        end

        write(output, "\n")
    end

    # Close the output file
    close(output)
end

"""
La función write_multiple_xyz_FHI recive el nombre del archivo, las coordenadas de las diferentes configuraciones  y un arreglo
de los symbolos de los  de la configuración
"""

function write_multiple_xyz_FHI(filename, coords, symbols)
    # Open the output file for writing
    output = open(filename, "w")
    write(output, "total estructures $(length(coords))\n")
    # Write each structure to the output file
    for i in 1:length(coords)
        structure = coords[i]
        natoms, ndim = size(structure)

        # Write the header for each structure
        write(output, "estructura $i\n")
        write(output, "$natoms\n")
        write(output, "XYZ file written with Julia\n")

        
        # Write the coordinates and element symbols to the output file
        for j in 1:natoms
            write(output, "atom")
        
            for k in 1:ndim
                write(output, " $(structure[j,k])")
            
            end
            write(output, " $(symbols[j])")
            write(output, "\n")
        end


        write(output, "\n")
    end

    # Close the output file
    close(output)
end

"""
La función extraer_posiciones  recibe la ruta de archivos de las funciones write_multiple_xyz_FHI y write_multiple_xyz
y devuelve un arreglo con las poiciones  de  estos archivos
"""
function extraer_posiciones(archivo)
    estructuras = []
    posiciones_totales = []
    posiciones_estructura = []

    for linea in eachline(archivo)
        if occursin("total estructures", linea)
            estructuras = parse(Int, split(linea)[3])
            continue
        end

        if occursin("estructura", linea)
            if !isempty(posiciones_estructura)
                push!(posiciones_totales, posiciones_estructura)
            end
            posiciones_estructura = []
            continue
        end

        if occursin("XYZ file written with Julia", linea)
            continue
        end

        palabras = split(linea)
        if length(palabras) > 1
            coordenadas = [parse(Float64, palabra) for palabra in palabras[2:end]]
            push!(posiciones_estructura, coordenadas)
        end
    end

    if !isempty(posiciones_estructura)
        push!(posiciones_totales, posiciones_estructura)
    end
    return posiciones_totales
end


"""
Extrae posciones  de la evolución de un isomero
"""
function posiciones_iso(archivo)
    posiciones_totales = []
    posiciones_estructura = []


    for linea in eachline("archivo")
        if occursin("Relaxation step number ", linea)
            if !isempty(posiciones_estructura)
                push!(posiciones_totales, posiciones_estructura)
            end
            posiciones_estructura = []
            continue
        end

        palabras = split(linea)
        if length(palabras) > 1
            coordenadas = [parse(Float64, palabra) for palabra in palabras[2:4]]
            push!(posiciones_estructura, coordenadas)
        end
    end

    if !isempty(posiciones_estructura)
        push!(posiciones_totales, posiciones_estructura)
    end
    return posiciones_totales
end

"""
La función convertir_a_adjunto de convierte un arrays de arrays a una
matrix 
"""
function convertir_a_adjunto(vectores)
    matriz = Matrix{Float64}(undef, length(vectores), length(vectores[1]))
    
    for i in 1:length(vectores)
        matriz[i, :] = convert(Vector{Float64}, vectores[i])
    end
    
    return adjoint(matriz)
end


## Funciones para buscar  configuraciones de adsorción  de manera aleatoria y las reflexiones de los enantiomeros

"""
La función Compute_of_CM_R_CCM recive las coordenadas de la estructura
y devulve el centro de masa de la estructura, un radio  de una esfera 
que   cubre nuestra estructura y una matrix con las posiciones 
trasladadas  respecto al centro de masas
"""
function Compute_of_CM_R_CCM(coords)
    # Compute the center of mass of the atoms
    n = size(coords, 1)
    center_of_mass = sum(coords, dims=1) / n
    # Compute the distances between the center of mass and each atom
    distances = zeros(n)
    for i in 1:n
        distances[i] = norm((center_of_mass.-coords)[i,:])
    end
    # Compute the radius of the smallest circle that encloses the atoms
    radius = maximum(distances)
    #Compute the coordinates translated to the center of mass
    coordinates_origin = coords .- center_of_mass
    
    return center_of_mass, radius, coordinates_origin
end

"""
La función rotate_system  recibe  estructuras tipo esfericas  como 
Au34  y una estructura que se desee rotar alrededor de la estructura
esferica.
"""
function rotate_system(MO, MR ; Dmin = 2, Dmax = 3, Δ = 0.1,n = 10)
    #compute systems at the origin
    CM_MO, R_MO, O_MO = Compute_of_CM_R_CCM(MO)
    CM_MR, R_MR, O_MR = Compute_of_CM_R_CCM(MR)
    #compute  rotation 
    Rs = R_MO .+ [i for i in Dmin:Δ:Dmax]
    Rrand = rand(Rs)
    θ = rand(0:π/n:π) 
    ϕ = rand(0:π/n:2*π)
    x = Rrand * sin(θ) * cos(ϕ)
    y = Rrand * sin(θ) * sin(ϕ)
    z = Rrand * cos(θ) 

    T = [x y z]

    traslation = O_MR .+ T
    #generate random rotation matrices
    # The first rotation matrix is for the change of the condition inicial  the cysteine
    rot_matrix1 = qr(randn(3, 3)).Q 
    # The second rotation matrix is for the rotation of the cysteinee
    rot_matrix2 = qr(randn(3, 3)).Q
    #The third rotation matrix joins the previous matrices.
    rot_matrix3 = rot_matrix2*rot_matrix1
    #Apply the rotation matrix to the translated  MR
    RotMR = (rot_matrix3 * traslation')'
    
    return RotMR, O_MO, O_MR
end

"""
La función multiple_rotation recibe los elementos de la función 
rotate_system además de un n número de veces que quiera repetir el 
de la función rotate_system y te devuelve  un arreglo con las
diferentes rotaciones  de  las rotaciones que se desee rotar
"""
function multiple_rotation(MO,MR,n; Dmin = 2, Dmax = 3, Δ = 0.1, r=10)
    # Perform multiple rotations
    MR_otations = []
    for i in 1:n
        Rot = rotate_system(MO, MR ; Dmin =Dmin,Dmax = Dmax, Δ = Δ,n = r)[1]
    push!(MR_otations, Rot)
    end
    return MR_otations
end


"""
La función interaction_an_atom recibe las coordenadas de la estrcutura 
como el Au34  en el centro, un array con los simbolos de las posiciones 
 de los atomos  de la estructura que se rota, las coordenadas de lo que se 
desee rotar, además de 
"""
function interaction_an_atom(coordinate,symbols,rotations,index_atom; err = 0.3,interaction=2.5)
    
    coordinate_origin = Compute_of_CM_R_CCM(coordinate)[3]
    index = findfirst(symbols .== index_atom)
    tolerance = err
    selected_positions = []
    distances =[]
    closet_atoms = []
    for i in 1:length(rotations)
        rotate_selec = rotations[i]
        dis = norm.([vec(row) for row  in eachrow(rotate_selec[index,:]' .- coordinate_origin)])
        closest_atoms=  findall(isapprox.(dis, minimum(dis), atol=tolerance, rtol=tolerance))
        push!(distances, dis)
        push!(closet_atoms,closest_atoms)
        if length(closest_atoms) >= 2
            for j in 1:length(closest_atoms)
                for k in j+1:length(closest_atoms)
                    if isapprox(norm(rotate_selec[index,:].- coordinate_origin[closest_atoms[j], :]),interaction ,atol=tolerance, rtol=tolerance) && 
                        isapprox(norm(rotate_selec[index,:] .- coordinate_origin[closest_atoms[k], :]), interaction,atol=tolerance, rtol=tolerance)
                        push!(selected_positions, i)
                        break
                    end
                end
                if i in selected_positions
                    break
                end
            end
        end
    end
    return [rotations[i] for i in selected_positions]
end
       

function positions_symbols(index_atoms, symbols)
    positions = []
    for i in index_atoms
        p = findall(x -> x == i, symbols)
        push!(positions,p)
    end
    return reduce((x, y) -> union(x, y), positions)
end

function eliminate_interactions_atoms(positions,coordinates,symbols,index_atoms; interaction=2.5, tolerance = 0.3)
    
    coordinates_origin = Compute_of_CM_R_CCM(coordinates)[3]
    positions2 = copy(positions)

    distances = []
    closet_atoms = []
    s_p = []
    for i in 1:length(positions)
        posiciones_selec = positions2[i]
        for p in positions_symbols(index_atoms, symbols)
            dis = [norm(row) for row  in eachrow(posiciones_selec[p,:]' .- coordinates_origin)]
            closest_atoms =  findall(isapprox.(dis, minimum(dis), atol=tolerance, rtol=tolerance))
            push!(distances,dis)
            push!(closet_atoms,closest_atoms)
            for j in 1:length(closest_atoms)
                if  norm(posiciones_selec[p,:].- coordinates_origin[closest_atoms[j],:]) <interaction#isapprox(norm(posiciones_selec[p,:].- coordinates_origin[closest_atoms[j],:]),interaction,atol=tolerance,rtol=tolerance)
                    push!(s_p,i)
                end
            end
       end
    end
    return deleteat!(positions2, unique(s_p))
end
    


function Plano_3p(p1,p2,p3)
    v1 = vec(p2.-p1)
    v2 = vec(p3.-p1)
    A, B , C = cross(v1[1],v2[1])/norm(cross(v1[1],v2[1]))
    D = -A*p1[1][1]-B*p1[1][2]-C*p1[1][3]
    return A, B, C, D
end

function pmax_alejadas_element(p_elementos,simbolos,elementos)
    E = []
    for element in elementos
        el = findall(x -> x == element, simbolos)
        if length(el) == 1
            push!(E,el)
        else
            i = argmax([norm(p_elementos[E[1],:].-p_elementos[i,:]) for i in el])
            push!(E,[el[i]])
        end
    end
    return E
end

# Función para reflejar un punto a través de un plano
function reflect_point(point, plane_parameters)
    A, B, C, D = plane_parameters

    # Calcula la distancia entre el punto y el plano
    distance = (A * point[1] + B * point[2] + C * point[3] + D) / sqrt(A^2 + B^2 + C^2)

    # Calcula la imagen reflejada del punto
    reflected_point = [
        point[1] - 2 * A * distance,
        point[2] - 2 * B * distance,
        point[3] - 2 * C * distance
    ]

    return reflected_point
end


function get_enantiomers(positions,p_elementos,simbolos,elementos)
    r_pos = []
    Pcluster = positions[1][1:34]
    atoms_plane = pmax_alejadas_element(p_elementos,simbolos,elementos)
    for p in positions
        P = p[35:end]
        plane = Plano_3p(P[atoms_plane[1]],P[atoms_plane[2]],P[atoms_plane[3]])
        enantiomers = [reflect_point(point, plane) for point in P]
        push!(r_pos,vcat(Pcluster,enantiomers))
    end
    return r_pos
end

##Funciones para realizar el análisis

function distancias_(coordenadas)
    n = size(coordenadas,1)
    distancias_entre_atómos = []
    for i in 1:n
        for j in i+1:n
            push!(distancias_entre_atómos,norm(coordenadas[i,:]-coordenadas[j,:]))
         end
    end
    
    centro_de_masas = sum(coordenadas[1:34,:], dims=1) / n
    co = coordenadas[1:34,:].-centro_de_masas
    distancias = []
    for i in 1:n
        push!(distancias,norm(co[i,:]))
    end
    return distancias_entre_atómos,distancias,co
end

"""
Devuelve la distancia  entre los átomos de un cumulo
"""
function distancias_entre_atomos(coordenadas)
    n = size(coordenadas,1)
    distancias_entre_atómos = []
    for i in 1:n
        for j in i+1:n
            push!(distancias_entre_atómos,norm(coordenadas[i,:]-coordenadas[j,:]))
        end
    end
    return distancias_entre_atómos
end

"""
Devuelve los átomos  que tienen el npumero de vecinos selecionado
"""
function átomos_con_n_vecinos(coordenadas,r,n_vecinos)
    D = []
    for i in 1:length(coordenadas[:,1]) 
        distan = []
        for j in 1:length(coordenadas[:,1])
            d = norm(coordenadas[i,:]-coordenadas[j,:])
            push!(distan,d)
        end
        push!(D,distan)
    end
    closet_atoms=[]
    for i in 1:length(coordenadas[:,1])
        push!(closet_atoms, findall(x -> r > x, D[i]))
    end
    return findall(x -> n_vecinos == x, length.(closet_atoms))
end

"""
Traslada las coordenadas al origen
"""
function coordendas_origen(coordenadas)
    centro_de_masas = sum(coordenadas, dims=1) / n
    co = coordenadas.-centro_de_masas
    return co
end

"""
Distancias de los atómos al  origen
"""
function distancias_al_origen(coordenadas)
    distancias = []
    for i in 1:length(coordenadas)
        push!(distancias,norm(coordenadas[i,:]))
    end
    return distancias
end







##Forma de calcular la medida de quiralidad de hausdorff con julia

function ρ(Q,Qp) 
    n = size(Q,1)
    sup = Float64[]
    for i in 1:n
        inf = Float64[]
        for j in 1:n
            dis = norm(Q[i,:].-Qp[j,:])
            push!(inf,dis)
        end
        push!(sup,minimum(inf))
    end
    return maximum(sup)    
end

function dQ(coordenadas)
    n = size(coordenadas,1)
    dmax = 0.0
    for i in 1:n 
        for j in i+1:n
             # faster: d = @views norm(coordenadas[:,i]-coordenadas[:,j])
            d = norm(coordenadas[i,:]-coordenadas[j,:]) 
            dmax = max(dmax, d)
        end
    end
    return dmax
end

function HCM(q;Δ =0.025, Δ2 = 0.5,I=0.2)
    A = Compute_of_CM_R_CCM(q)[3]
    δ = (2*π)/360
    B = -1 *A
    R= [RotXYZ(i, j, k) for i in -π:π*Δ:π, j in -π:π*Δ:π, k in -π/2:π*Δ:π/2]
    angulos =[(i,j,k) for i in -π:π*Δ:π, j in -π:π*Δ:π, k in -π/2:π*Δ:π/2]
    HH = Float64[]
    for i in 1:length(R)
        BB = B*R[i]
        HCM = max(ρ(A,BB),ρ(BB,A))/dQ(A)
        #HH = min(HCM,HH)
        push!(HH,HCM)
    end
    
 
        min_b1 = argmin(HH)
        θs = angulos[min_b1]
   
        HH2 = Float64[]

        R2 =[RotXYZ(i, j, k) for i in θs[1]-I*π:Δ2*δ:θs[1]+I*π, 
            j in θs[2]-I*π:Δ2*δ:θs[2]+I*π, k in θs[3]-I/2*π:Δ2*δ:θs[3]+I/2*π]
        angulos2 =[(i,j,k) for i in θs[1]-I*π:Δ2*δ:θs[1]+I*π, 
            j in θs[2]-I*π:Δ2*δ:θs[2]+I*π, k in θs[3]-I/2*π:Δ2*δ:θs[3]+I/2*π]
    
        HH2 = Float64[]
        for i in 1:length(R2)
            BB = B*R2[i]
            HCM = max(ρ(A,BB),ρ(BB,A))/dQ(A)
            push!(HH2,HCM)
        end
        return minimum(HH2), HH2, angulos2, minimum(HH), HH, angulos
 
end




