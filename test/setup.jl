function metadata(Nant, Nfreq)
    frame = ReferenceFrame()
    obs = measure(frame, observatory("OVRO_MMA"), pos"ITRF")
    time = Epoch(epoch"UTC", (2015-1858)*365*24*3600*seconds)
    set!(frame, obs)
    set!(frame, time)

    antennas = TTCal.Antenna[]
    for ant = 1:Nant
        position = Position(pos"ITRF", obs.x+100randn(), obs.y+100randn(), obs.z+100randn())
        push!(antennas, TTCal.Antenna(position))
    end

    baselines = TTCal.Baseline[]
    for ant1 = 1:Nant, ant2 = ant1:Nant
        push!(baselines, TTCal.Baseline(ant1, ant2))
    end

    channels = Nfreq == 1? [45e6] : collect(linspace(45e6, 70e6, Nfreq))
    phase_center = measure(frame, Direction(dir"AZEL", 0degrees, 90degrees), dir"ITRF")
    beam = SineBeam()

    Metadata(antennas, baselines, channels, phase_center, time, beam)
end

