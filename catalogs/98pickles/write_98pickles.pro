pro write_98pickles
;jm08jul25nyu - parse the Pickles stellar library (which I got from
;               C. Tremonti)

    path = filepath('',root_dir=getenv('CATALOGS_DIR'),$
      subdirectory='98pickles')

    star = file_search(path+'*.flam',count=nstar)

    readcol, path+'pickles.wavelength', wave, $
      format='F', /silent
    npix = n_elements(wave)
    
    readcol, path+'pickles_tbl2.dat', id, type, teff, $
      uv, bv, vr, vi, vj, vh, vk, feh, mbol, bcv, $
      format='I,A,F,F,F,F,F,F,F,F,F,F,F', /silent

    p98 = {$
      id:     0, $
      type:  '', $
      teff: 0.0, $ 
      uv:   0.0, $
      bv:   0.0, $
      vr:   0.0, $
      vi:   0.0, $
      vj:   0.0, $
      vh:   0.0, $
      vk:   0.0, $
      feh:  0.0, $
      mbol: 0.0, $
      bcv:  0.0, $
      wave: wave, $
      flux: fltarr(npix)}
    p98 = replicate(p98,nstar)

    p98.id = id
    p98.type = repstr(type,"'",'')
    p98.teff = teff
    p98.uv = uv
    p98.bv = bv
    p98.vr = vr
    p98.vi = vi
    p98.vj = vj
    p98.vh = vh
    p98.vk = vk
    p98.feh = feh
    p98.mbol = mbol
    p98.bcv = bcv
    struct_print, struct_trimtags(p98,except=['wave','flux'])

    for ii = 0L, nstar-1L do begin
       readcol, path+'uk'+strlowcase(p98[ii].type)+'.flam', $
         format='F', flux, /silent
       p98[ii].flux = flux
    endfor
    
    srt = sort(p98.id)
    p98 = p98[srt]

    splog, 'Writing '+path+'98pickles.fits'
    mwrfits, p98, path+'98pickles.fits', /create

return
end
    
