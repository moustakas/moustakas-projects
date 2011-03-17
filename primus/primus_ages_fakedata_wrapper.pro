pro fakedata_filenames, sim, fakefile, photofile, onedfile, $
  psfile=psfile, simtype=simtype, filters=filters, $
  photosuffix=photosuffix, test=test, refband=refband
; sim - simulation number [1-6]    

    outpath = im_primus_path()+'fakedata/'
    fakepath = getenv('PRIMUS_DATA')+'/fakedata/'

    if (n_elements(sim) eq 0L) then sim = 1
    case sim of
       1: begin
          simtype = 1 & photosuffix = 'BRI'
          filters = 'bessell_'+['B','R','I']+'.par'
          refband = 1L
       end
       2: begin
          simtype = 1 & photosuffix = 'UBVRI'
          filters = 'bessell_'+['U','B','V','R','I']+'.par'
          refband = 3L
       end
       3: begin
          simtype = 1 & photosuffix = 'ugriz'
          filters = 'sdss_'+['u','g','r','i','z']+'0.par'
          refband = 2L
       end
       4: begin
          simtype = 2 & photosuffix = 'BRI'
          filters = 'bessell_'+['B','R','I']+'.par'
          refband = 1L
       end
       5: begin
          simtype = 2 & photosuffix = 'UBVRI'
          filters = 'bessell_'+['U','B','V','R','I']+'.par'
          refband = 3L
       end
       6: begin
          simtype = 2 & photosuffix = 'ugriz'
          filters = 'sdss_'+['u','g','r','i','z']+'0.par'
          refband = 2L
       end

    endcase

; build the filenames    

    if keyword_set(test) then prefix = 'test' else prefix = 'fakedata'
    fakefile = fakepath+prefix+'_sim'+$
      string(simtype,format='(I2.2)')+'_'+$
      string(sim,format='(I0)')+'.fits.gz'
    psfile = outpath+repstr(file_basename(fakefile),'.fits.gz','.'+photosuffix+'.ps')
    photofile = repstr(fakefile,'.fits.gz','.'+photosuffix+'.fits.gz')
    onedfile = outpath+prefix+'_sim'+string(simtype,format='(I2.2)')+$
      '_'+string(sim,format='(I0)')+'-zAll.fits.gz'
    
return
end
    
pro primus_ages_fakedata_wrapper, sims, nmock=nmock, paramplots=paramplots, $
  zvzplots=zvzplots, writesim=writesim, fitredshift=fitredshift, test=test
; jm08jul25nyu - generate various PRIMUS/AGES fakedata sets to test
;                the PHOTOZ code
; echo "primus_ages_fakedata_wrapper, /writemock" | idl > & log &
    
    outpath = im_primus_path()+'fakedata/'
    fakepath = getenv('PRIMUS_DATA')+'/fakedata/'

    if (n_elements(nmock) eq 0L) then nmock = 1000L ; number of mock galaxies

; specify which simulations to run

    if (n_elements(sims) eq 0L) then sims = [1,2,3,4,5,6]
    nsim = n_elements(sims)
    
; test the filename routine    
;   for isim = 1, 6 do begin
;      fakedata_filenames, isim, fakefile, photofile, onedfile, $
;        onedfile_photoz, simtype=simtype, filters=filters, $
;        photosuffix=photosuffix
;      print, isim, ' ', file_basename(fakefile)+' '+file_basename(photofile)+' '+$
;        file_basename(onedfile)+' '+file_basename(onedfile_photoz)
;   endfor

; generate the fakedata samples

    if keyword_set(writesim) then begin
       for isim = 0L, nsim-1L do begin ; iterate on each simulation
          fakedata_filenames, sims[isim], fakefile, photofile, onedfile, $
            simtype=simtype, filters=filters, photosuffix=photosuffix, test=test
          make_ages_fakedata, repstr(fakefile,'.gz',''), simtype=simtype, $
            nmock=nmock, /photoz, photosuffix=photosuffix, filterlist=filters
       endfor
    endif

; fit the simulated samples
    
    if keyword_set(fitredshift) then begin
       for isim = 0L, nsim-1L do begin ; iterate on each simulation
          fakedata_filenames, sims[isim], fakefile, photofile, onedfile, $
            simtype=simtype, filters=filters, photosuffix=photosuffix, test=test
          primus_fit_redshift, fakefile, photofile=photofile, outdir=outpath, $
            /exactfile, /nointrachip, /nostar, /noagn, /nopowerlaw, ftweak=0.0
       endfor
    endif

; make the z versus z plots    
    
    if keyword_set(zvzplots) then begin

       template_file = getenv('PRIMUS_DATA')+$
         '/template_basis/ages_basis_set.fits.gz'
       splog, '  Reading '+template_file
       template_data = mrdfits(template_file,1,/silent)
       
       xpage = 8.5 & ypage = 11.0
;      psfile = outpath+'fakedata_test_photoz.ps'

       for isim = 0L, nsim-1L do begin ; iterate on each simulation

          fakedata_filenames, sims[isim], fakefile, photofile, onedfile, $
            psfile=psfile, simtype=simtype, filters=filters, $
            photosuffix=photosuffix, refband=refband, test=test
          splog, 'Generating '+file_basename(psfile)
          dfpsplot, psfile, xsize=xpage, ysize=ypage, /color

          nophotozlabel = 'sim'+string(simtype,format='(I2.2)')
          photozlabel = nophotozlabel+'/'+photosuffix
          
          extract = mrdfits(fakefile,1)
          photoinfo = mrdfits(photofile,1)
          oned = mrdfits(onedfile,1)

          primus_zvz_qaplot, extract, photoinfo, oned, $
            template_data, nophotozlabel, photozlabel, $
            photosuffix, refband

          dfpsclose
          spawn, 'ps2pdf13 '+psfile+' '+repstr(psfile,'.ps','.pdf')
          spawn, 'gzip -f '+psfile
       
       endfor

    endif
    
; make some plot to set the simulation parameters    

    if keyword_set(paramplots) then begin

; magnitude-S/N relation for a representative mask        

       psfile = outpath+'fakedata_params.ps'
       dfpsplot, psfile, /square, /color
       
       thismask = 'd02b0061' & rerun = '0215'
       primus_read_1dinputs, thismask, rerun, $
         extract=ee, slits=ss

       good = where((ee.sn1 gt 0.1) and (ee.sn2 gt 0.1))
       xx = [ss[good].magr,ss[good].magr]
       yy = alog10([ee[good].sn1,ee[good].sn2])

; S/N distribution       
       
       stats = im_stats(yy,sigrej=3.0)
       xstats = [$
         '\Delta = '+strtrim(string(stats.mean_rej,format='(F12.3)'),2),$
         '\sigma = '+strtrim(string(stats.sigma_rej,format='(F12.3)'),2),$
         '\Delta_{med} = '+strtrim(string(stats.median_rej,format='(F12.3)'),2)]
;      xstats = [strtrim(string(stats.minrej,format='(F12.3)'),2)+'-'+$
;        strtrim(string(stats.maxrej,format='(F12.3)'),2),$
;        strtrim(string(stats.median_rej,format='(F12.3)'),2)+' ('+$
;        strtrim(string(stats.mean_rej,format='(F12.3)'),2)+$
;        '+/-'+strtrim(string(stats.sigma_rej,format='(F12.3)'),2)+')']

       im_plothist, yy, bin=0.1, /noplot, xx1, yy1
       im_plothist, yy, bin=0.1, yrange=[0.0,1.05], xrange=alog10([0.5,600]), $
         xsty=1, ysty=1, charsize=1.9, xthick=4.0, ythick=4.0, $
         charthick=3.0, ytitle='Fraction', xtitle='log (S/N)', $
         thick=4.0, norm=max(yy1)
       legend, textoidl(xstats), /left, /top, box=0, $
         charsize=1.6, charthick=3.0, /clear
       legend, thismask+'/'+rerun, /right, /top, box=0, $
         charsize=1.6, charthick=3.0, /clear
       
; S/N vs magnitude relation       
       
       xrange = [17.5,25] & yrange = alog10([0.5,600])
       djs_plot, [0], [0], /nodata, xrange=xrange, $
         yrange=yrange, xsty=1, ysty=1, charsize=1.9, $
         xthick=4.0, ythick=4.0, charthick=3.0, $
         xtitle='R (AB mag)', ytitle='log (S/N)'
       djs_oplot, xx, yy, ps=6, color='grey', sym=0.5, thick=2.0
       legend, thismask+'/'+rerun, /right, /top, box=0, $
         charsize=1.6, charthick=3.0, /clear

       minx = 18.5 & maxx = 23.9
       med = im_medxbin(xx,yy,0.5,/ver,minx=minx,maxx=maxx)
       djs_oplot, med.binctr, med.medy, ps=4, symsize=3.0, $
         color='navy', thick=5.0
       coeff = linfit(med.binctr,med.medy) ; R vs S/N
       coeff2 = [-coeff[0]/coeff[1],1.0/coeff[1]] ; S/N vs R
       raxis = findgen((24.5-18.0)/0.01+1)*0.01+18.0
       djs_oplot, raxis, poly(raxis,coeff), thick=4.0

       thisfit = 'log (S/N) = '+strtrim(string(coeff[0],format='(F12.2)'),2)+$
         strtrim(string(coeff[1],format='(F12.3)'),2)+'*R'
       thisfit2 = 'R = '+strtrim(string(coeff2[0],format='(F12.1)'),2)+$
         strtrim(string(coeff2[1],format='(F12.3)'),2)+'*log (S/N)'
       legend, textoidl([thisfit,thisfit2]), /left, /bottom, box=0, $
         charsize=1.6, charthick=3.0, /clear
       
; S/N (photometry) vs magnitude relation       

       mbase = 24.0 & sigmbase = 0.05 & minerr = 0.03
       mags = findgen((26.0-17.0)/0.01+1)*0.01+17.0
       nmags = n_elements(mags)
       maggies = 10^(-0.4*mags)
       maggiesivar = maggies*0.0

       maggiesivar = 1.0/(minerr^2+(sigmbase*10.0^(-0.4*(mbase-mags)))^2)/$
         (0.4*alog(10.0))^2/maggies^2
       errmaggies = 1.0/sqrt(maggiesivar)
       maggies = maggies + randomn(seed,nmags)*errmaggies
       
       xrange = [16.5,26.5] & yrange = [0.0,45.0]
       djs_plot, [0], [0], /nodata, xrange=xrange, $
         yrange=yrange, xsty=1, ysty=1, charsize=1.9, $
         xthick=4.0, ythick=4.0, charthick=3.0, $
         xtitle='Magnitude (AB)', ytitle='S/N (photometry)'
;      djs_oplot, mags, maggies/errmaggies, line=0, thick=5.0
       djs_oplot, mags, maggies/errmaggies, ps=4, thick=2.0
       legend, thismask+'/'+rerun, /right, /top, box=0, $
         charsize=1.6, charthick=3.0, /clear

       dfpsclose
       spawn, 'ps2pdf13 '+psfile+' '+repstr(psfile,'.ps','.pdf')
       spawn, 'gzip -f '+psfile

    endif
    
return
end
    
