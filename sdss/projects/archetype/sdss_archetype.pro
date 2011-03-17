pro sdss_archetype
; jm08may15nyu

    platefit_init, 'archetype_parfile.par', /mpa
    
    mpainfo = read_sdss_mpa(/info)
    mpagalindx = read_sdss_mpa(/indx)

    zval = [0.004,0.020,0.050]

    allplates = mpainfo.plateid
    plate = allplates[uniq(allplates,sort(allplates))]
    nplate = n_elements(plate)

    fiberfit1 = struct_addtags({plateid: 0L, fiberid: 0L, mjd: 0L, $
      z: 0.0, v_disp: 0.0, e_bv_sfd: 0.0, $
      spectro_mag: fltarr(3), kcor_mag: fltarr(3)},$
      im_empty_structure(struct_trimtags($
      mpagalindx[0],select=['MODEL_*','TAUV_*','BEST_*'])))
     
    iplate = 0L & ispec = 1L & maxspec = 10L ; n_elements(mpainfo)
    while (iplate lt nplate) and (ispec le maxspec) do begin

       these = where(plate[iplate] eq allplates,nspec)
       fiberfit = replicate(fiberfit1,nspec)
       struct_assign, mpainfo[these], fiberfit
       struct_assign, mpagalindx[these], fiberfit, /nozero
       fiberfit.z = mpainfo[these].z
       fiberfit.v_disp = mpainfo[these].v_disp
       fiberfit.e_bv_sfd = mpainfo[these].e_bv_sfd

       mjd = mpainfo[these].mjd
       fibers = mpainfo[these].fiberid ; fiber numbers
       
; read the 1D spectra; correct for foreground Galactic reddening 

       splog, 'Reading '+string(nspec,format='(I0)')+' fibers '+$
         'from plate '+string(plate[iplate],format='(I4.4)')
       readspec, fiberfit.plateid, fiberfit.fiberid, wave=vacwave, $
         flux=flux, invvar=invvar, /align
;      readspec, plate[iplate], fibers, wave=vacwave, /align, $
;        flux=flux, invvar=invvar

;      vacwave = 10.0^logwave ; vacuum!
;      wave = vacwave & vactoair, wave ; air! (observed frame)

       restflux = flux*0.0
       restinvvar = invvar*0.0
       continuum = flux*0.0
       for ii = 0L, nspec-1L do begin

          fiber_dered, fiberfit[ii], alog10(vacwave), flux[*,ii], $
            invvar[*,ii], restwl=restwave, restflux=restflux1, $
            err=restferr1
          restflux[*,ii] = restflux1
          restinvvar[*,ii] = 1.0/(restferr1+(restferr1 ge 0.0))^2.0 * (restferr1 ne 0.0)

          match_z = where(abs(fiberfit[ii].best_model_z-zval) lt 0.0001)
          resample_model, alog10(vacwave), fiberfit[ii].z, fiberfit[ii].v_disp, match_z
          cont = bc_model_combine(restwave,[fiberfit[ii].tauv_cont,fiberfit[ii].model_coef])
          continuum[*,ii] = cont

          plot, restwave, restflux[*,ii], xsty=3, ysty=3, ps=10
          djs_oplot, restwave, continuum[*,ii], ps=10, color='green'

       endfor

stop       
       
       splog, 'Correcting for foreground Galactic reddening'
       flux = flux_dust*0.0
       invvar = invvar_dust*0.0
       for ii = 0L, nspec-1L do begin
          dustcor = 10.0^(-0.4*k_lambda(wave,/odon)*fiberfit[ii].e_bv_sfd)
          flux[*,ii] = flux_dust[*,ii]*dustcor
          invvar[*,ii] = invvar_dust[*,ii]/dustcor^2.0
       endfor
       
; loop through each fiber and rebuild the continuum fit       
       
       init = 1
       metallicities = fiberfit.best_model_z
       for ii = 0L, nspec-1L do begin

          if fiberfit[ii].best_model_z eq 0.0 then continue
          cont = bc_model_restore(fiberfit[ii],logwave) ; vacuum!

          if (init) then begin
             dims = size(cont,/dimen)
             continuum = fltarr(nspec,dims[0])
             restwave = fltarr(nspec,dims[0])
             init = 0
          endif 
          
          continuum[ii,*] = cont
          restwave[ii,*] = wave/(1.0+fiberfit[ii].z)

          djs_plot, wave, flux[*,ii], ps=10, xsty=3, ysty=3
          djs_oplot, wave, cont, color='green', ps=10

stop
          
       endfor
    
; advance the indices       
       
       iplate = iplate+1L
       ispec = ispec+nspec
       
    endwhile
    
   ;; Use the restore_fiber function to read in again the plate.
   plate = restore_fiber(platenum)
   init = 1
   metallicities = plate.best_model_z
   for i = 0, n_elements(fibers) - 1 do begin
      j = fibers[i]

      if plate[ii].best_model_z eq 0.0 then continue

      cont = bc_model_restore(plate[ii], logwl, restwl)

      if (init) then begin
         dims = size(cont, /dimen)
         continuum = fltarr(n_elements(fibers), dims[0])
         lambda = fltarr(n_elements(fibers), dims[0])
         init = 0
      endif 

      continuum[i,*] = cont
      lambda[i,*] = airwl/(1.0+plate[ii].z)
   endfor
    
    
    restore_continuum_fits, plate[0], fibers=fibers, continuum=cont, $
      lambda=lambda, err=err, flux=flux, logz=logz
    
    
    
    modelspath = getenv('PLATEFIT_DIR')+'/bc03/fits/'
    bnorm = mrdfits(modelspath+'bc_models_subset_v4_0c_bursts.fit',4,/silent)

    z0 = where(mpagalindx.best_model_z eq 0.004)
    z2 = where(mpagalindx.best_model_z eq 0.02)
    z3 = where(mpagalindx.best_model_z eq 0.05)

    fib_mass = fltarr(n_elements(mpainfo))

    fib_mass[z0] = mpagalindx[z0].model_coef ## (1.0/bnorm[0,*])
    fib_mass[z2] = mpagalindx[z2].model_coef ## (1.0/bnorm[2,*])
    fib_mass[z3] = mpagalindx[z3].model_coef ## (1.0/bnorm[3,*])
    fib_mass = 1D-17*fib_mass*(4*!pi*dluminosity(out_postlss.z,/cm))^2/lsun
    
    

return
end
    
