pro calculations
; jm01sep7ufoa
; calculate the conversion from Jansky to photons (counts) in each of
; the SIRTF bands (calculate the gain)

    lambda0 = [3.56,4.50,5.69,7.95,23.8,72.5,157.0] ; effective wavelength [micron]
    dlambda = [0.7,1.0,1.4,2.9,4.7,19.0,35.0]       ; bandwidth [micron]

    throughput = [0.8886,0.9021,0.9081,0.9140,1.0,1.0,1.0] ; throughput [percent]

    aperture = 0.85 ; mirror size [m]

    time = [1.0,1.0,1.0,1.0,1.0,1.0,1.0]

    gain = 1.51E7 * (dlambda/lambda0) * throughput * aperture^2.0 * time ; [counts/Jy]
    

return
end
