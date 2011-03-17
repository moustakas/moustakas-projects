pro synth_tests, result, eigeninfo, snratio=snratio, nmonte=nmonte, $
  postscript=postscript
; jm03sep17uofa
; take the BC best-fitting coefficients for an atlas galaxy

    nsnr = n_elements(snratio)

; call this routine recursively    
    
    if nsnr gt 1L then begin

       for isnr = 0L, nsnr-1L do begin
          synth_tests, result1, snratio=snratio[isnr], nmonte=nmonte, $
            postscript=postscript, write=write
          if isnr eq 0L then result = result1 else result = [ [result], [result1] ]
       endfor

    endif else snratio = 30.0 ; default S/N ratio
    
    if n_elements(nmonte) eq 0L then nmonte = 25L

    stime0 = systime(1)
    splog, 'Signal-to-noise ratio '+string(snratio,format='(G0.0)')+'.'
    
; read the eigen-templates

    splog, 'Reading the eigen-templates.'
    eigendir = filepath('',root_dir=getenv('ISPEC_DIR'),subdirectory='templates')
    eigenfile = 'bc_020_templates.fits' ; BC03 solar metallicity
    
    eigenflux = mrdfits(eigendir+eigenfile,0,eigenhead,/silent)
    eigeninfo = mrdfits(eigendir+eigenfile,1,/silent)
    eigenwave = make_wave(eigenhead)
    ntemplate = eigeninfo.ntemplates

; restore the results of fitting the atlas

    atlas = read_integrated()
    natlas = n_elements(atlas)

; initialize the output data structure

    splog, 'Initializing the output data structure.'
    result = {$
      galaxy:      '', $
      specfile:    '', $
      snratio:    0.0, $
      coeff_true: fltarr(ntemplate), $
      ebv_true:   0.0, $
      chi2_true:  0.0, $
      chi2_fit:   fltarr(nmonte), $
      coeff:      fltarr(ntemplate,nmonte) $
      }
    result = replicate(result,natlas)
    
    result.galaxy = atlas.galaxy
    result.specfile = atlas.driftfile
    result.snratio = snratio
    result.chi2_true = atlas.continuum_chi2
    result.ebv_true = atlas.ebv[0] ; <--- this will need to be changed!!!
    result.coeff_true = atlas.starcoeff

; loop on each object in the atlas

    minwave = 3700.0
    maxwave = 6800.0
    dwave = 2.75
    wave = findgen((maxwave-minwave)/dwave+1)*dwave+minwave
    npix = n_elements(wave)

    kl = k_lambda(wave,/charlot) ; <--- change this later!!!!
;   kl = k_lambda(wave,/odonnell)
       
    for j = 0L, natlas-1L do begin
    
; form the best-fitting model

       splog, 'Preparing spectrum '+strn(j+1)+'/'+strn(natlas)+' ['+strn(result[j].galaxy)+'].'
       
       bigflux = eigenflux # result[j].coeff_true
       bigfluxerr = bigflux/result[j].snratio
       
; degrade the number of Angstroms per pixel to a more manageable size

       combine1fiber, alog10(eigenwave), bigflux, bigfluxerr^(-2.0), $
         newloglam=alog10(wave), newflux=flux, newivar=invvar
       ferr = invvar^(-1./2)
       
; match the templates to the model; redden each template by the "true"
; amount of reddening

       starflux = imatch_templates(eigenflux,eigenwave,wave,$
         starres=3.0,newres=3.0)
       for k = 0L, ntemplate-1L do starflux[*,k] = starflux[*,k]*10^(-0.4*result[j].ebv_true*kl)

; Monte Carlo this bad boy

       for i = 0L, nmonte-1L do begin

          print, format='("Monte Carlo iteration #",I0,"/",I0,".",A1,$)', i+1, nmonte, string(13b)
          
          if i gt 0L then pflux = flux+randomn(seed,npix)*ferr else pflux = flux
          backfit = ibackfit(pflux,wave,starflux,invvar=invvar,nmonte=0,ndust=0,/silent)

          result[j].coeff[*,i] = backfit.starcoeff
          result[j].chi2_fit[i] = backfit.continuum_chi2

       endfor

    endfor 

    splog, format='("Total time = ",G0," minutes.")', (systime(1)-stime0)/60.0

; generate a save set and generate postscript output if requested

    cmsave, result, eigeninfo, filename='synth_tests.idlsave'
    plot_synth_tests, postscript=postscript

stop
    
return
end    
