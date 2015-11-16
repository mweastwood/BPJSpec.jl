function xyz2uvw(x,y,z)
    Nant = length(x)
    Nbase = div(Nant*(Nant-1),2) + Nant
    u = Array{Float64}(Nbase)
    v = Array{Float64}(Nbase)
    w = Array{Float64}(Nbase)
    α = 1
    for i = 1:Nant, j = i:Nant
        u[α] = x[j]-x[i]
        v[α] = y[j]-y[i]
        w[α] = z[j]-z[i]
        α += 1
    end
    u,v,w
end

function ant1ant2(Nant)
    Nbase = div(Nant*(Nant-1),2) + Nant
    ant1 = Array(Int32,Nbase)
    ant2 = Array(Int32,Nbase)
    α = 1
    for i = 1:Nant, j = i:Nant
        ant1[α] = i
        ant2[α] = j
        α += 1
    end
    ant1,ant2
end

function createms(Nant,Nfreq)
    Nbase = div(Nant*(Nant-1),2) + Nant

    # Let's put all of our antennas at the north pole
    # with the phase center at the north pole because
    # this makes the conversion from antenna positions
    # to uvw coordinates very straight forward.

    # Once TTCal has the ability to compute more general
    # baseline coordinates, we can change this to be
    # a little more realistic.

    Rearth = 6.3674447e6
    x = 20*randn(Nant)
    y = 20*randn(Nant)
    z =  5*randn(Nant) + Rearth

    u,v,w = xyz2uvw(x,y,z)
    ν = linspace(40e6,60e6,Nfreq) |> collect
    t = (2015.-1858.)*365.*24.*60.*60. # a rough current Julian date (in seconds)
    ant1,ant2 = ant1ant2(Nant)

    frame = ReferenceFrame()
    pos = Position(pos"ITRF",0.0,0.0,Rearth) # north pole
    set!(frame,Epoch(epoch"UTC",Quantity(t,"s")))
    set!(frame,pos)
    north = Direction(dir"APP",q"0.0deg",q"90.0deg")
    phase_dir = measure(frame,north,dir"J2000")

    name  = tempname()*".ms"
    table = Table(name)

    subtable = Table("$name/SPECTRAL_WINDOW")
    Tables.addRows!(subtable,1)
    subtable["CHAN_FREQ"] = reshape(ν,length(ν),1)
    unlock(subtable)

    subtable = Table("$name/ANTENNA")
    Tables.addRows!(subtable,Nant)
    subtable["POSITION"] = [x y z]'
    unlock(subtable)

    subtable = Table("$name/FIELD")
    Tables.addRows!(subtable,1)
    subtable["PHASE_DIR"] = reshape([longitude(phase_dir);latitude(phase_dir)],2,1)
    unlock(subtable)

    Tables.addRows!(table,Nbase)
    table[kw"SPECTRAL_WINDOW"] = "Table: $name/SPECTRAL_WINDOW"
    table[kw"ANTENNA"] = "Table: $name/ANTENNA"
    table[kw"FIELD"] = "Table: $name/FIELD"
    table["ANTENNA1"] = ant1-1
    table["ANTENNA2"] = ant2-1
    table["UVW"] = [u v w]'
    table["TIME"] = fill(float(t),Nbase)
    table["FLAG_ROW"] = zeros(Bool,Nbase)
    table["FLAG"] = zeros(Bool,4,Nfreq,Nbase)
    table["DATA"] = zeros(Complex64,4,Nfreq,Nbase)
    unlock(table)

    name,MeasurementSet(name)
end

