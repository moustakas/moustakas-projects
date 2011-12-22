function read_10pozzetti, zbin, type234=type234, type1=type1, nolog=nolog, salpeter=salpeter
; jm11aug31ucsd - read the Pozzetti+10 stellar mass functions; h=0.7,
; Omega0=0.25, OmegaL=0.75, Chabrier IMF

    path = getenv('CATALOGS_DIR')+'/10pozzetti/'
    file = path+'global_mass_function.dat'
    if keyword_set(type1) then file = path+'MF_type1.dat'
    if keyword_set(type234) then file = path+'MF_type234.dat'

    data = rsex(file)
    ndata = n_elements(data)
    
    if n_elements(zbin) ne 0 then begin
       case zbin of
          1: begin
             masslim = 8.3
             if keyword_set(type1) then masslim = 9.0
             if keyword_set(type234) then masslim = 8.3
             these = where(data.zlo eq 0.10D and data.mup ge masslim and data.ngal gt 0,ndata)
          end
          2: begin
             masslim = 9.4
             if keyword_set(type1) then masslim = 9.8
             if keyword_set(type234) then masslim = 9.4
             these = where(data.zlo eq 0.35D and data.mup ge masslim and data.ngal gt 0,ndata)
          end
          3: begin
             masslim = 10.1
             if keyword_set(type1) then masslim = 10.3
             if keyword_set(type234) then masslim = 10.0
             these = where(data.zlo eq 0.55D and data.mup ge masslim and data.ngal gt 0,ndata)
          end
          4: begin
             masslim = 10.6
             if keyword_set(type1) then masslim = 10.6
             if keyword_set(type234) then masslim = 10.4
             these = where(data.zlo eq 0.75D and data.mup ge masslim and data.ngal gt 0,ndata)
          end
          else: message, 'ZBIN must be 1-4!'
       endcase
    endif else these = lindgen(ndata)

;   if n_elements(zbin) ne 0 then begin
;      case zbin of
;         1: these = where(data.zlo eq 0.10D and data.meanmass gt 8.2  and data.ngal gt 0,ndata)
;         2: these = where(data.zlo eq 0.35D and data.meanmass gt 9.4  and data.ngal gt 0,ndata)
;         3: these = where(data.zlo eq 0.55D and data.meanmass gt 10.1 and data.ngal gt 0,ndata)
;         4: these = where(data.zlo eq 0.75D and data.meanmass gt 10.6 and data.ngal gt 0,ndata)
;         else: message, 'ZBIN must be 1-4!'
;      endcase
;   endif else these = lindgen(ndata)

    mf = replicate({ngal: 0, mass: 0.0, phi: 0.0, phierr: 0.0},ndata)
    mf.ngal = data[these].ngal
    mf.mass = (data[these].mup-data[these].mlo)/2.0+data[these].mlo
;   mf.mass = data[these].wmeanmass
    mf.phi = data[these].phi
    mf.phierr = data[these].phierr

    if keyword_set(salpeter) then mf.mass += 0.26 ; Chabrier-->Salpeter
    
    if (keyword_set(nolog) eq 0) then begin
       mf.phierr = mf.phierr/mf.phi/alog(10)
       mf.phi = alog10(mf.phi)
    endif

return, mf
end
