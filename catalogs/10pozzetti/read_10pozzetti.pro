function read_10pozzetti, zbin=zbin, nolog=nolog, salpeter=salpeter
; jm11aug31ucsd - read the Pozzetti+10 stellar mass functions; h=0.7,
; Omega0=0.25, OmegaL=0.75, Chabrier IMF

    path = getenv('CATALOGS_DIR')+'/10pozzetti/'
    globalfile = path+'global_mass_function.dat'

    data = rsex(globalfile)
    ndata = n_elements(data)
    
    case zbin of
       1: these = where(data.zlo eq 0.10D and data.meanmass gt 8.2  and data.ngal gt 0,ndata)
       2: these = where(data.zlo eq 0.35D and data.meanmass gt 9.4  and data.ngal gt 0,ndata)
       3: these = where(data.zlo eq 0.55D and data.meanmass gt 10.1 and data.ngal gt 0,ndata)
       4: these = where(data.zlo eq 0.75D and data.meanmass gt 10.6 and data.ngal gt 0,ndata)
       else: these = lindgen(ndata)
    endcase

    mf = replicate({ngal: 0, mass: 0.0, phi: 0.0, phierr: 0.0},ndata)
    mf.ngal = data[these].ngal
    mf.mass = data[these].wmeanmass
    mf.phi = data[these].phi
    mf.phierr = data[these].phierr

    if keyword_set(salpeter) then mf.mass += 0.26 ; Chabrier-->Salpeter
    
    if (keyword_set(nolog) eq 0) then begin
       mf.phierr = mf.phierr/mf.phi/alog(10)
       mf.phi = alog10(mf.phi)
    endif

return, mf
end
