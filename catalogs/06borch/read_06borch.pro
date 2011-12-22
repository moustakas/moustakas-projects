function read_06borch, zbin, blue=blue, red=red, nolog=nolog, $
  h100=h100, salpeter=salpeter
; jm11dec20ucsd - read the Borch et al. 2006 stellar mass functions;
; h=0.7, Omega0=0.3, OmegaL=0.7, Chabrier IMF

    if n_elements(h100) eq 0 then h100 = 0.7
    
    file = getenv('CATALOGS_DIR')+'/06borch/c17_mfs.txt'

    data = rsex(file)
    ndata = n_elements(data)

; mass limits eyeballed from Fig 9    
    if n_elements(zbin) ne 0 then begin
       case zbin of
          1: begin
             these = where(data.zbin eq 0.3D,ndata)
             minmass = 10.0
          end
          2: begin
             these = where(data.zbin eq 0.5D,ndata)
             minmass = 10.25
          end
          3: begin
             these = where(data.zbin eq 0.7D,ndata)
             minmass = 10.5
          end
          4: begin
             these = where(data.zbin eq 0.9D,ndata)
             minmass = 10.75
          end
          else: message, 'ZBIN must be 1-4!'
       endcase
    endif else these = lindgen(ndata)

    mf = replicate({mass: 0.0, phi: 0.0, phierr: 0.0},ndata)
    mf.mass = data[these].mass

    mf.phi = data[these].phitot
    mf.phierr = (data[these].phitot_up-data[these].phitot_lo)/2.0
    if keyword_set(blue) then begin
       mf.phi = data[these].phiblue
       mf.phierr = (data[these].phiblue_up-data[these].phiblue_lo)/2.0
    endif
    if keyword_set(red) then begin
       mf.phi = data[these].phired
       mf.phierr = (data[these].phired_up-data[these].phired_lo)/2.0
    endif
    good = where(mf.phi ne 0D and mf.mass ge minmass,ngood)
    mf = mf[good]
    
    if keyword_set(salpeter) then mf.mass += 0.26 ; Chabrier-->Salpeter
    
    if (keyword_set(nolog) eq 0) then begin
       mf.phierr = mf.phierr/mf.phi/alog(10)
       mf.phi = alog10(mf.phi)
    endif

return, mf
end
