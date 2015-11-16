let
    frame = ReferenceFrame()
    set!(frame,Position(pos"ITRF",q"0.0m",q"0.0deg",q"90.0deg"))
    set!(frame,Epoch(epoch"UTC",Quantity(50237.29,"d")))
    beam = itrf_beam(frame,SineBeam(1.0),45e6)
    for i = 1:length(beam)
        θ,ϕ = LibHealpix.pix2ang_ring(512,i)
        if θ < π/2
            @test abs(beam[i] - cos(θ)) < 1e-5
        else
            @test abs(beam[i]) < 1e-5
        end
    end
end

