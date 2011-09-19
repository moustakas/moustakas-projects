pro mz_fix_posterior, sdss=sdss
; jm11sep19ucsd - fix a bug in how the posterior quantites are
; calculated (one-time fix!)

    isedpath = mz_path(/isedfit)
    ndraw = isedfit_ndraw() ; number of random draws

    chunkfile = file_search(isedpath+'sfhgrid01/charlot/ages_'+$
      'bc03_chab_chunk_????.fits.gz',count=nchunk)
    for jj = 0, nchunk-1 do begin
       modelgrid1 = mrdfits(chunkfile[jj],1)
       modelgrid1 = struct_trimtags(temporary(modelgrid1),except='modelmaggies')
       if (jj eq 0) then modelgrid = temporary(modelgrid1) else $
         modelgrid = [temporary(modelgrid),temporary(modelgrid1)]
    endfor
    nmodel = n_elements(modelgrid)
    nage = n_elements(modelgrid[0].age)
    nallmodel = nmodel*nage
    
    bigage    = reform(modelgrid.age,nallmodel)
    bigmass   = reform(modelgrid.mstar,nallmodel)
    bigtau    = reform(rebin(reform(modelgrid.tau,1,nmodel),nage,nmodel),nallmodel)
    bigZ      = reform(rebin(reform(modelgrid.Z,1,nmodel),nage,nmodel),nallmodel)
    bigav     = reform(rebin(reform(modelgrid.mu*modelgrid.av,1,nmodel),nage,nmodel),nallmodel)
    bignburst = reform(rebin(reform(modelgrid.nburst,1,nmodel),nage,nmodel),nallmodel)
    
    splog, 'Reconstructing SFHs'
    bigsfr = bigage*0D
    bigsfr100 = bigage*0D       ; average over the previous 100 Myr
    bigb100 = bigage*0D         ; birthrate parameter
    bigmgal = bigage*0D         ; galaxy mass ignoring mass loss 
    for imod = 0L, nmodel-1 do begin
       tindx = lindgen(nage)+imod*nage
       sfr = isedfit_reconstruct_sfh(modelgrid[imod],outage=bigage[tindx],$
         sfr100=sfr100,b100=b100,mgalaxy=mgal)
       bigsfr[tindx] = sfr      ; alog10(sfr)
       bigsfr100[tindx] = sfr100 ; alog10(sfr100) 
       bigb100[tindx] = b100
       bigmgal[tindx] = mgal
    endfor

; read the old results
    isedfile = isedpath+'ages_bc03_chab_charlot_sfhgrid01.fits.gz'
    old = mrdfits(isedfile,1)
    post = mrdfits(repstr(isedfile,'.fits.gz','_post.fits.gz'),1)
    new = old
    good = where(old.modelindx gt -1,ngal)

; see isedfit_compute_posterior       
    for gg = 0L, ngal-1 do begin
       if ((gg mod 1000) eq 0) then splog, gg, ngal
       logscale_err = post[good[gg]].scale_err/post[good[gg]].scale/alog(10)
       logscale = alog10(post[good[gg]].scale) + randomn(seed,ndraw)*logscale_err

       new[good[gg]] = isedfit_packit(new[good[gg]],alog10(bigmass[post[good[gg]].draws])+logscale,type='mass')
       new[good[gg]] = isedfit_packit(new[good[gg]],alog10(bigsfr[post[good[gg]].draws])+logscale,type='sfr')
       new[good[gg]] = isedfit_packit(new[good[gg]],alog10(bigsfr100[post[good[gg]].draws])+logscale,type='sfr100')

       new[good[gg]] = isedfit_packit(new[good[gg]],bigage[post[good[gg]].draws],type='age')
       new[good[gg]] = isedfit_packit(new[good[gg]],bigtau[post[good[gg]].draws],type='tau')
       new[good[gg]] = isedfit_packit(new[good[gg]],bigZ[post[good[gg]].draws],type='Z')
       new[good[gg]] = isedfit_packit(new[good[gg]],bigav[post[good[gg]].draws],type='av')
       new[good[gg]] = isedfit_packit(new[good[gg]],bigb100[post[good[gg]].draws],type='b100')
    endfor

    outfile = repstr(isedfile,'.gz','')
    im_mwrfits, new, outfile, /clobber

stop    
    
return
end
    
