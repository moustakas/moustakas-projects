pro build_cosmicimf_mass, mass
; jm10mar22ucsd - build the stellar mass estimates for the
; AGES/COSMICIMF sample

    cosmicimfpath = ages_path(/projects)+'cosmicimf/'
    isedpath = cosmicimfpath+'isedfit/'

    splog, '#########################'
    splog, 'Building the AGES stellar masses'
    sample = read_cosmicimf_sample()
    ngal = n_elements(sample)

; IMF grid: range of high-mass slopes plus Salpeter
    peg = read_cosmicimf_sample(/pegase)
    imf = strtrim(peg.imf,2)
    nimf = n_elements(peg)
    
    mass = {$
      ages_id:            0L,$
      z:                 0.0,$
      sfhgrid02:  fltarr(nimf),$
      sfhgrid03:  fltarr(nimf)}
    mass = replicate(mass,ngal)
    mass.ages_id = sample.ages_id
    mass.z = sample.z

; read and parse all the isedfit tables    
    for ii = 0, nimf-1 do begin
       imfpath = isedpath+imf[ii]+'/'
; from build_pegase_basemodels
       case imf[ii] of 
          'Salpeter': imfstr = 'salp'
          'Kroupa02': imfstr = 'kroupa'
          else: if strmatch(imf[ii],'*cosmicimf*') then $
            imfstr = imf[ii] else message, 'Fix me'
       endcase 
; sfhgrid02
       isedfile = 'BwRIJHKsirac_pegase_'+imfstr+'_calzetti_sfhgrid02.fits.gz'
       ised = mrdfits(imfpath+isedfile,1,rows=mass.ages_id)
       mass.sfhgrid02[ii] = ised.mass
; sfhgrid03
       isedfile = 'BwRIJHKsirac_pegase_'+imfstr+'_sfhgrid03.fits.gz'
       ised = mrdfits(imfpath+isedfile,1,rows=mass.ages_id)
       mass.sfhgrid03[ii] = ised.mass
    endfor

    im_mwrfits, mass, cosmicimfpath+'ages_cosmicimf_mass.fits', /clobber

; for convenience write out the offsets with respect to Salpeter
    off = replicate({imf: '', slope: 0.0, dmsalp_sfhgrid02: 0.0, $
      dmsalp_sfhgrid02_err: 0.0, dmsalp_sfhgrid03: 0.0, $
      dmsalp_sfhgrid03_err: 0.0},nimf)
    off.imf = imf
    off.slope = peg.slope
    for ii = 0, nimf-1 do begin
       off[ii].dmsalp_sfhgrid02 = djs_median(mass.sfhgrid02[0]-mass.sfhgrid02[ii])
       off[ii].dmsalp_sfhgrid02_err = djsig(mass.sfhgrid02[0]-mass.sfhgrid02[ii],sigrej=5.0)
       off[ii].dmsalp_sfhgrid03 = djs_median(mass.sfhgrid03[0]-mass.sfhgrid03[ii])
       off[ii].dmsalp_sfhgrid03_err = djsig(mass.sfhgrid03[0]-mass.sfhgrid03[ii],sigrej=5.0)
    endfor
    im_mwrfits, off, cosmicimfpath+'dmsalp.fits', /clobber

return
end
