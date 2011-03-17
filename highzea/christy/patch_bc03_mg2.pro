
; Patch Mg II region in BC03 SSP file
outfile = '/home/tremonti/mmt/models/' + $
          'bc03_pad94_chab_z05_ssp_uvpatch_mmtres.fit'

; BC03 Model
bc03file = '/home/tremonti/models/bc03/bc03_padova1994_chab_z05_ssp.fit'

; Read UVBLUE SSPs
mdir = '/home/tremonti/models/UVBLUE/p05/'
ps = mrdfits('/home/tremonti/mmt/models/uvblue_synth.fit', 1)
; Convert to air??
alam = ps.wave
vactoair, alam
alam = alam[sort(alam)]

;------------------------------------------------------------
; Resample BC03 models 

m = mrdfits(bc03file, 1)
nm = n_elements(m.age)

; Make output wavelength grid
npix = (alog10(10000) - alog10(1100))/ 1e-4
log_mwave = findgen(npix)*1e-4 + alog10(1100)
n = {wave: 10.0^log_mwave, flux: fltarr(npix, nm), age: m.age}
;normwl = where(mwave gt 5400 and mwave lt 5600)

bc03_res = 75.0 ; BC03 spectral resolution??
mmt_res = 92.0 ; 
sig_pix = sqrt((mmt_res^2 - bc03_res^2))/ 69.0
uv = where(n.wave lt 3200)
for ii = 0, nm - 1 do begin
  linterp, m.wave, m.flux[*,ii], n.wave, mfluxi
  n.flux[*,ii] = gconv(mfluxi, sig_pix)
  n.flux[uv,ii] = mfluxi[uv] ; don't convolve UV wavelengths!
endfor

;-------------------------------------------------------------
; Patch models more than 2 Gyr old 

for ii = 0, nm - 1 do begin

  if n.age[ii] gt 2e9 then continue
  
  plot, n.wave, n.flux[*,ii], xr=[2500, 3500], $
         title=strtrim(n.age[ii]/1e6, 2) + ' Myr'

  nmod1 = median(n.flux[where(n.wave gt 3000 and n.wave lt 3200), ii])
  scont1 = smooth(medsmooth(n.flux[*,ii], 250), 30)
  
  ; We are missing models less than 10 Myr
  if n.age[ii] lt 10.0e6 then begin
    if n.age[ii] lt 1.0e6 then uv_file = 't35000g45p05k2.flx'
    if n.age[ii] gt 1.0e6 and n.age[ii] lt 2.0e6 then $
       uv_file = 't32000g40p05k2.flx' 
    if n.age[ii] gt 2.0e6 and n.age[ii] lt 3.0e6 then $
       uv_file = 't28000g40p05k2.flx'
    if n.age[ii] gt 3.0e6 and n.age[ii] lt 4.0e6 then $
       uv_file = 't23000g35p05k2.flx'    
    if n.age[ii] gt 4.0e6 and n.age[ii] lt 7.0e6 then $
       uv_file = 't22000g35p05k2.flx'    
    if n.age[ii] gt 7.0e6 and n.age[ii] lt 10.0e6 then $
       uv_file = 't21000g35p05k2.flx'    
    
    s = readUVBLUE(mdir + uv_file)
    bflam1 = broadenUVBLUE(s.lam, s.flam, vel=mmt_res)
    bflam2 = broadenUVBLUE(ps.wave, ps.flux[*,0], vel=mmt_res)

    bflam = bflam1 * 100 + bflam2
    
  endif else begin
    junk = min(abs(m.age[ii] - ps.age), indx) 
    bflam = broadenUVBLUE(ps.wave, ps.flux[*,indx], vel=mmt_res)
  endelse

  linterp, alam, bflam, n.wave, uvfluxi

  ;nmod2 = median(uvfluxi[where(n.wave gt 3000 and n.wave lt 3200)]) 
  ;oplot, n.wave, smooth(uvfluxi/nmod2*nmod1, 10), color=!blue

  scont2 = smooth(medsmooth(uvfluxi, 250), 30)
  new_mflux = uvfluxi / scont2 * scont1
  oplot, n.wave, new_mflux, color=!red
  
  ; Patch in Mg II
  ok = where(n.wave gt 2700 and n.wave lt 3000)  
  n.flux[ok,ii] = new_mflux[ok]
  oplot, n.wave[ok], n.flux[ok,ii], color=!green  

  wait, 2
endfor

mwrfits, n, outfile, /create

end
