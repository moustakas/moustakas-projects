function read_chaos_sample
; jm13apr01siena - read the sample and pack into a useful format

    file = getenv('CHAOS_DATA')+'/extract/n628_f1_cal.fits'
    data = mrdfits(file,0,hdr,/silent)
    data = data[*,1:n_elements(data[0,*])-1] ; get rid of the response function

    errfile = getenv('CHAOS_DATA')+'/extract/n628_f1_err.fits'
    errdata = mrdfits(errfile,0,/silent)
    errdata = errdata[*,1:n_elements(errdata[0,*])-1]
       
    wave = make_wave(hdr)
;   ploterror, wave, data[*,0], errdata[*,0], /trad
    
    dim = size(data,/dim)
    npix = dim[0]
    ngal = dim[1]
    
    spec1d = replicate({id: 0, z: 0.0, wave: wave, $
      flux: fltarr(npix), ferr: fltarr(npix)},ngal)
    spec1d.id = lindgen(ngal)
    spec1d.z = 0.002192 ; hack!
    spec1d.flux = data
    spec1d.ferr = errdata

; total hack!!!    
    for ii = 0, ngal-1 do begin
       gd = where(errdata[*,ii] gt 0.0,comp=crap,ncomp=ncrap)
       if ncrap ne 0 then spec1d[ii].ferr[crap] = max(errdata[gd,ii])
    endfor
    
return, spec1d
end
