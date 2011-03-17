;; AGES.pro
;
; AIM: remove OH sky residuals from AGES spectra
; 
; INPUT: data - a plate of spectra
;        info - structure containing redshifts of objects
;
; OUTPUT: ps file of PCA process 
;         new data array of skysubtracted spectra
;
;****************************************************************************

FUNCTION vw_RESISTANT_Mean,Y,CUT,Num_Rej
;+
; NAME:
;    Vw_Resistant_Mean
;
; PURPOSE:
;    Outlier-resistant determination of the mean and standard deviation.
;
; EXPLANATION:
;    Vw_Resistant_Mean trims away outliers using the median and the median
;    absolute deviation.    An approximation formula is used to correct for
;    the truncation caused by trimming away outliers
;
; CALLING SEQUENCE:
;    Vw_Resistant_Mean, VECTOR, Sigma_CUT, Mean, Sigma_Mean, Num_RejECTED
;
; INPUT ARGUMENT:
;       VECTOR    = Vector to average
;       Sigma_CUT = Data more than this number of standard deviations from the
;               median is ignored. Suggested values: 2.0 and up.
;
; OUTPUT ARGUMENT:
;       Mean  = the mean of the input vector, numeric scalar
; OPTIONAL OUTPUTS:
;       Sigma_Mean = the approximate standard deviation of the mean, numeric
;            scalar.  This is the Sigma of the distribution divided by sqrt(N-1)
;            where N is the number of unrejected points. The larger
;            SIGMA_CUT, the more accurate. It will tend to underestimate the
;            true uncertainty of the mean, and this may become significant for
;            cuts of 2.0 or less.
;       Num_RejECTED = the number of points trimmed, integer scalar
;
; EXAMPLE:
;       IDL> a = randomn(seed, 10000)    ;Normal distribution with 10000 pts
;       IDL> Vw_Resistant_Mean,a, 3, mean, meansig, num    ;3 Sigma clipping
;       IDL> print, mean, meansig,num
;
;       The mean should be near 0, and meansig should be near 0.01 ( =
;        1/sqrt(10000) ).
; PROCEDURES USED:
;       AVG() - compute simple mean
; REVISION HISTORY:
;       Written, H. Freudenreich, STX, 1989; Second iteration added 5/91.
;       Use MEDIAN(/EVEN)    W. Landsman   April 2002
;       Correct conditional test, higher order truncation correction formula
;                R. Arendt/W. Landsman   June 2002
;       New truncation formula for sigma H. Freudenriech  July 2002
;       Divide Sigma_mean by Num_good rather than Npts W. Landsman/A. Conley
;                          January 2006
;       vw 30th Jan 2006 - change to function and remove division by
;                          No. pixels
;
;-

 On_Error,2
 if N_params() LT 2 then begin
    print,'Syntax - result=Vw_Resistant_Mean(Vector, Sigma_cut, [Num_Rejected]'
    return,-1.
 endif

 Npts    = N_Elements(Y)
 YMed    = MEDIAN(Y,/EVEN,/double)
 AbsDev  = ABS(Y-YMed)
 MedAbsDev = MEDIAN(AbsDev,/EVEN,/double)/0.6745
 IF MedAbsDev LT 1.0E-24 THEN MedAbsDev = AVG(AbsDev)/.8

 Cutoff    = Cut*MedAbsDev

 GoodPts = Y[ WHERE( AbsDev LE Cutoff, Num_Good ) ]
 Mean    = AVG( GoodPts )
 Sigma   = SQRT( TOTAL((GoodPts-Mean)^2)/Num_Good )
 Num_Rej = Npts - Num_Good

; Compenate Sigma for truncation (formula by HF):
 SC = Cut > 1.0
 IF SC LE 4.50 THEN $
   SIGMA=SIGMA/(-0.15405+0.90723*SC-0.23584*SC^2+0.020142*SC^3)

 Cutoff = Cut*Sigma

 GoodPts = Y[ WHERE( AbsDev LE Cutoff, Num_Good ) ]
 mean    = AVG( GoodPts )
 Sigma   = SQRT( TOTAL((GoodPts-mean)^2)/Num_Good )
 Num_Rej = Npts - Num_Good

; Fixed bug (should check for SC not Sigma) & add higher order correction
 SC = Cut > 1.0
 IF SC LE 4.50 THEN $
   SIGMA=SIGMA/(-0.15405+0.90723*SC-0.23584*SC^2+0.020142*SC^3)

; Now the standard deviation of the mean:
; Sigma = Sigma/SQRT(Num_Good-1.)



 RETURN, [mean,sigma]

 END

FUNCTION rms5,x
;calculate 67th percentile
    y = MEAN(x,/nan)
    index = SORT(ABS(x-y))
    a = N_ELEMENTS(x)*.67
    percentile = ABS((x-y)[index[a]])

    RETURN,percentile
 END

;------------------------------------------------------------------

;for plotting bad sky pixels - I'm sure there's much more efficient
;                              ways of writing these bits of code!

PRO plotbits, lambda, spec, filter, _extra=extra

    fil2 = filter

    nbin = n_elements(fil2)
;print, nbin

    if (fil2[0]) then fil2[0] = 2
    n=2

    for i=1, nbin -1 do begin

       if NOT(fil2[i]) then continue ;=1 leave alone

       if (fil2[i-1]) then fil2[i] = fil2[i-1] ;previous=0, this is in same chunk
       if NOT(fil2[i-1])  then fil2[i] = n ;previous=1, this is start of next chunk
       
       if (i eq nbin-1) then break

       if NOT(fil2[i+1]) then n=n+1 ;next=0, this is end of chunk
    endfor


    nchunk = max(fil2)-1
;print, nchunk
;help, lambda, spec

;plot each chunk seperately to avoid joining lines
    for i=0,nchunk-1 do begin
       index = where(fil2 eq i+2)
       djs_oplot, lambda[index],spec[index], _extra=extra
    endfor


 end
;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;for plotting none sky pixels
PRO plotbits2,lambda,spec,filter

    fil2 = filter
    nbin = n_elements(fil2)

    if NOT(fil2[0]) then fil2[0] = 2
    n=2

    for i=1, nbin -1 do begin

       if (fil2[i]) then continue ;=1 leave alone

       if NOT(fil2[i-1]) then fil2[i] = fil2[i-1] ;previous=0, this is in same chunk
       if (fil2[i-1])  then fil2[i] = n ;previous=1, this is start of next chunk
       
       if (i eq nbin-1) then break

       if (fil2[i+1]) then n=n+1 ;next=0, this is end of chunk
    endfor

    nchunk = max(fil2)-1
;plot each chunk seperately to avoid joining lines
    for i=0,nchunk-1 do begin
       index = where(fil2 eq i+2)
       oplot, lambda[index],spec[index],linestyle = 1
    endfor

 END

;------------------------------------------------------------------
;------------------------------------------------------------------

function ages_skysubpca, data, info, pcainfo, qaplotname=qaplotname, scale=scale

    if (n_elements(scale) eq 0L) then scale = 1.0D ; jm07feb28nyu
    
;20 appears reasonable. Increase if there are some plates which need
;it. 
    maxnrecon=20                    

;;*** if ngalaxy is too low PCA may not work. Just keep an eye on
;;*** it. Could add the sky fibers in to help?
    ngalaxy = n_elements(info)
    splog, 'No of galaxies on plate ',ngalaxy

    if (n_elements(data) eq 0L) or (ngalaxy eq 0L) then begin
       print, 'Syntax - newdata = ages_skysubpca(data,info,qaplotname=)'
       return, -1L
    endif

    if (n_elements(qaplotname) eq 0L) then begin
       splog, 'QAPLOTNAME input required.'
       return, -1L
    endif

    emgal_lines_file = getenv('SUBTRACTOH_DIR')+'/emgal_lines.dat'
;   emgal_lines_file = '~/idl/pro/SDSS/SKY/PUBLIC/emgal_lines.dat'
    if (file_test(emgal_lines_file,/regular) eq 0L) then begin
       splog, 'Unable to find line-list file '+emgal_lines_file+'.'
       return, -1L
    endif

    splog, 'Opening '+qaplotname+'.'
    postthick = 2.1
    !p.thick = postthick & !x.thick = postthick & !y.thick = postthick
    !p.charthick = postthick & !p.charsize = 2.0

    dfpsplot, qaplotname, /color

    newdata = data              ; initialize the output data structure 
    pcainfo = replicate({nrecon:0L, RMS_nosky:0.0, RMS_b4sky:0.0, $
      RMS_afsky:0.0},ngalaxy)   ;useful info that should be kept

; --------------------------------------------------
; only consider wavelengths redward of WAVE_CUT
; --------------------------------------------------

    wave_cut = 6700.0
    wave_indx = where(data.wave gt wave_cut,npix)

;; vw31jan06 - multiplied by 1e17, no reason other than got fed up
;;             with dealing with very small numbers
    skyflux = data.skyflux[wave_indx,*]*scale
    skyerr = data.skyferr[wave_indx,*]*scale
    wave = data.wave[wave_indx]
    nsky = (size(data.skyflux,/dimensions))[1]

    flux = data.flux[wave_indx,*]*scale
    ferr = data.ferr[wave_indx,*]*scale

; --------------------------------------------------
; distinguish sky/non-sky pixels
; --------------------------------------------------

; each sky fiber is effected differently by exponential slope in far
; red... median filter them first - additive effect.
    skyflux2 = skyflux
    for i=0L, nsky-1L do skyflux2[*,i] = $
      skyflux[*,i]-medfilt((skyflux[*,i]),findgen(npix),75) ;vw06jan30

    rms = fltarr(npix)
    
;for i=0L, npix-1L do rms[i] = djsig(skyflux[i,*],sigrej=4.0) ; jm06jan22uofa
    for i=0L, npix-1L do rms[i] =  (vw_resistant_mean(skyflux2[i,*],4.0))[1] ;vw06jan30

    !p.multi=[0,1,3]
    plot, wave, rms, title='identifying bad sky pixels: sky-fibre rms/pixel'

    med2 = medfilt(rms,findgen(npix),75)
    plot, wave, rms/med2, title='same as above, with continuum taken out'
    djs_oplot, !x.crange, [1,1], line=0, color='green'


; these levels can be tweaked to get more/fewer pixels identified, but
; pushing the limit too low won't in any way improve things (PCA can't
; do anything with Poisson noise) and might hinder.

    nonskypix = where(rms/med2 lt 1.2,compl=skypix)

    splog,'skypix,noskypix',n_elements(skypix),n_elements(nonskypix)

; plot it out to check it's sensible
    djs_plot, wave[nonskypix], rms[nonskypix]/med2[nonskypix], yrange=[0.2,6], $
      title='sky pixels (red), nonsky pixels (black)'
    djs_oplot, wave[skypix], rms[skypix]/med2[skypix], color='red', ps=4, syms=0.5

    array = fltarr(npix,ngalaxy) ; array for PCA
    med_gal = fltarr(npix,ngalaxy) ; continuum of each galaxy
    mask = lonarr(npix,ngalaxy) ; line and error mask 

; --------------------------------------------------
; now prepare each galaxy spectrum
; --------------------------------------------------
    !p.multi=[0,1,5]
    for i=0L, ngalaxy-1L do begin

; create mask array containing pixels likely to have emission lines
;;*** For Broad Emission line objects you need to change to the
;;*** correct mask file
       tmpflux = fltarr(npix)
       mask[*,i] = masklines(flux[*,i],wave/(1+info[i].z),/silent,$
         emgal_lines_file,tmpflux)

; estimate galaxy continuum with a median filter (emission lines masked)
       med_gal[*,i] = medfilt(tmpflux,findgen(npix),75)
       
; subtract galaxy continuum
       array[*,i] = tmpflux - med_gal[*,i]

; plot a few - is the redshift right?
;    if i le 9 then plot, wave/(1+info[i].z),flux[*,i],/xstyle

; normalise by errors
       ind = where(ferr[*,i] gt 0.0 and flux[*,i] gt -10.0,compl=compl)
       array[ind,i] = array[ind,i]/ferr[ind,i]

; add bad pixels into emission line mask array
       if compl[0] ne -1L then mask[compl,i]=2
       
    endfor
    
; replace the bad pixels and pixels likely to have emission lines with
; the mean/median/weighted mean (for that same pixel) of all the other
; spectra with good values

    for i=0L,npix-1L do begin
       ind = where(mask[i,*] ne 0,compl=compl)
       if ind[0] ne -1L and compl[0] ne -1L then $
         array[i,ind] = median(array[i,compl])
       
;;*** need to sort this out if it becomes a problem
;;*** would be worth removing this pixel from the PCA i.e. if it is a sky
;;*** pixel then demote it to a non-sky pixel
;    if n_elements(compl) lt 10 then $
;      splog,'This plate has a pixel with <10 good values...'
       
    endfor

;------------------------------------------------------------------
; create eigenspectra
;------------------------------------------------------------------
    
; PCA; PCOMP is a standard IDL routine, with a few adaptations I made

    goodgal = indgen(ngalaxy)
    array2 = array[skypix,*]

;; potential here to add a loop to throw out outliers and recreate
;; eigenspectra... actually it didn't seem necessary.
    for j=0,0 do begin
       
       array3 = array2[*,goodgal]
       pcs = PCOMP3(array3,coeff=espec,/covariance,/nodiv,variances=var,/double)

; look for outliers in PC amplitudes
       pcsd = fltarr(n_elements(skypix),2) ;mean and standard deviation
       outlier = fltarr(n_elements(goodgal))
       for i=0L,n_elements(skypix)-1L do pcsd[i,*] = median(pcs[i,*])
;      for i=0L,n_elements(skypix)-1L do pcsd[i,*] = vw_resistant_mean(pcs[i,*],4.)
       
; how many pcs to use to throw out outliers? certainly no more than
; maxnrecon and probably less. Only becomes an issue if need to use
; it. 
       for k=0L,n_elements(goodgal)-1L do $
         outlier[k] = max(abs(pcs[0:10,k]-pcsd[0:10,0])/pcsd[0:10,1])
       ind = where(outlier lt 5,count)
       splog, 'no. good spectra ',count,' out of ',ngalaxy

       goodgal = ind

       !p.multi = [0,2,2]

; PC distributions: if the distribution is a very narrow spike, with a
; single 'outlying' object with a very large value (negative or
; positive), that one object is dominating the eigenspectrum. This
; isn't any good for reconstructing general properties of the
; dataset - the loop is there to remove the worst outliers.

; why are only the first 8 sky pixels checked?  jm06jan22uofa
; the pcs are "principal component amplitudes" = the amount of each
; eigenvector present in each spectral reconstruction. These plots check the
; first 8 eigenspectra (not pixels) are not dominated by outlieing spectra.

       for i=0L, 7L do plothist, pcs[i,*], $
         title='distribution of PC amplitudes, PC '+$
         string(i,form='(I2)'), charsize=1, /halfbin, xbin, ybin
       
    endfor

    
; these variances are a check to see how much of the variance of the
; data you are picking up in the top eigenspectra.
; look at the plot to see how it falls of. In general you don't want
; to use eigenspectra too far beyond the knee as they contain only noise
; or elements from single unusual spectra.

    splog, 'variances: ', var[0:4]

    !p.multi = [0,1,5]
    plot, var,/xstyle,title='variances of eigenspectra'
    oplot, [maxnrecon,maxnrecon],minmax(var),linestyle=1


; plot the eigenspectra: hard work to spot what they are actually
; telling you. Make sure there's nothing obviously wacky.
    
    filter=uintarr(npix)
    filter[skypix]=1
    
    for i=0L, 8L do begin

       plot, wave[skypix],espec[*,i],/nodata,$
         title='eigenspectrum '+string(i,form='(I2)')
       tmp = fltarr(npix)
       tmp[skypix] = espec[*,i]
       plotbits, wave, tmp, filter
       plotbits2, wave, fltarr(npix), filter
       
    endfor

    
; --------------------------------------------------
; reconstruct sky: this would be a useful sanity check, but I didn't
; get around to it.
; --------------------------------------------------
    
    


; --------------------------------------------------
; reconstruct objects
; --------------------------------------------------

    flux2 = flux

    !p.multi=[0,1,4]
    nn_noskysub = 0L
    nrecon = lonarr(ngalaxy)+1
    for i=0L, ngalaxy-1L do begin
       
       newflux = flux[*,i]-med_gal[*,i]

; divide by errors
       ind = where(mask[*,i] ne 2.)
       if ind[0] ne -1 then newflux[ind] = newflux[ind]/ferr[ind,i]

; non emission line, good pixels
       ind_good_sky = where(mask[skypix,i] eq 0) 
       ind_good_nonsky = where(mask[nonskypix,i] eq 0) 

; check whether PCA sky subtraction useful
       RMS_nosky = (vw_resistant_mean(newflux[nonskypix[ind_good_nonsky]],4))[1]
       RMS_b4sky = (vw_resistant_mean(newflux[skypix[ind_good_sky]],4))[1]

       pcainfo.RMS_nosky = RMS_nosky
       pcainfo.RMS_b4sky = RMS_b4sky

       if RMS_b4sky le RMS_nosky then begin
          newdata.flux[wave_indx,i] = data.flux[wave_indx,i]
          nn_noskysub = nn_noskysub+1
          recon =fltarr(n_elements(skypix))
          nrecon[i]=0

;yes - I know. Told you the code needed tidying up!
          goto,plots
       endif

; --------------------------------------------------
; reconstruct the sky noise
; --------------------------------------------------

; There will be a range of numbers of eigenspectra used, depending on
; the quality of each individual spectrum. If something goes wrong,
; the number of eigenspectra used will shoot up, for example if
; there's a misplaced strong galaxy emission line its trying really
; hard to reconstruct. PCA is simply a "least squares" procedure.
       
       pcs2 = newflux[skypix[ind_good_sky]] ## (TRANSPOSE(espec[ind_good_sky,0:maxnrecon-1]))

       recon =fltarr(n_elements(skypix))

;; here we incrementally reconstruct the sky signal, each time
;; checking the new rms against that of the non-sky pixels. Of course
;; this is only an approximation of what is wanted and feel free to
;; change it. But its the best thing we came up with for the SDSS
;; spectra. 

       while (1) do begin

          recon = recon+pcs2[nrecon[i]-1]*espec[*,nrecon[i]-1]
          
;; below is the optimal solution for gappy data, but it's giving a few
;; spurious results and I don't have time to look into why right
;; now. It may improve things a little though in spectra with many bad
;; pixels.  Gappy-PCA: Connolly+Szalay, 1999, AJ, 117. ask for the
;; code if you want to try it.  
;        weights = fltarr(n_elements(skypix))+1.
;        weights[ind_good_sky]=1.
;        pcs2 = gappy_ages(weights,newflux[skypix],espec[*,0:nrecon[i]])
;        recon = pca_reconstruct(pcs2,espec[*,0:nrecon[i]])


          RMS_afsky = (vw_resistant_mean(newflux[skypix[ind_good_sky]]-$
            recon[ind_good_sky],4))[1]
          
          if (RMS_afsky lt RMS_nosky) or (nrecon[i] eq maxnrecon-1) then begin

             flux2[skypix,i] = flux[skypix,i]-recon*ferr[skypix,i] ; this is it!
             
             newdata.flux[wave_indx,i] = flux2[*,i]/scale
             
             pcainfo.RMS_afsky = RMS_afsky
             
             break
          endif
          nrecon[i]=nrecon[i]+1
       endwhile


; plot some of the results
       plots:
       if i/10 eq i/10. then begin
          
          plot, wave,flux[*,i],title='original spectrum'+string(i,format='(I4)')+' + masked pixels (red)'
          ind = where(mask[*,i] eq 0,compl=compl)
          ind1 = where(mask[*,i] eq 1) ;emission lines
          ind2 = where(mask[*,i] eq 2) ;bad pixels
          
          if ind[0] ne -1 then mask[ind,i]=-10
          if compl[0] ne -1 then  mask[compl,i]=0
          
          djs_oplot, wave,mask[*,i],color='green' ;bad pixels+emission lines
          if ind2[0] ne -1 then mask[ind2,i] = -10
          djs_oplot, wave,mask[*,i],color='red' ;emission lines only
          
          plot, wave,array[*,i],title='input spec to PCA: med sub, err norm, em lines/bad pixels replaced',ytick_get = yrange
          djs_oplot, wave[skypix],array[skypix,i],psym=1,color='red'
          
          plot, wave[skypix],recon,title='PCA reconstruction '+$
            string(nrecon[i],form='(I0)'),yrange=minmax(yrange),$
            /ystyle,psym=1
          recon2 = fltarr(npix)
          recon2[skypix]=recon
          
          plot, wave,flux[*,i],title='sky residual subtracted spectrum (red)'
          djs_oplot, wave,flux2[*,i],color='red'
          
       endif
       
    endfor

    pcainfo.nrecon=nrecon

;; this gives you a quick idea of how bad your plate is for sky
;; residual noise. 
    splog, 'number not in need of sky subtraction',nn_noskysub

; analysis plots: make sure you don't get "upside down hedgehogs"
; anywhere, it means you've created a high pass filter and are
; removing Poisson noise.

    rms1 = fltarr(npix)
    rms2 = rms1 & rms3 = rms1 & rms4 = rms1

    for i=0L, npix-1L do rms1[i] = (vw_resistant_mean((flux[i,*]-med_gal[i,*]),4.0))[1]
    for i=0L, npix-1L do rms2[i] = (vw_resistant_mean((flux2[i,*]-med_gal[i,*]),4.0))[1] 

    for i=0L, npix-1L do begin
       ind = where(ferr[i,*] gt 0.0)
       rms3[i] = (vw_resistant_mean((flux[i,ind]-med_gal[i,ind])/ferr[i,ind],4.0))[1]
    endfor

    for i=0L, npix-1L do begin
       ind = where(ferr[i,*] gt 0.0)
       rms4[i] = (vw_resistant_mean((flux2[i,ind]-med_gal[i,ind])/ferr[i,ind],4.0))[1]
    endfor

    !p.multi=[0,1,2]
    plot, wave, rms1,title='rms/pixel,67th percentile, before residual subtraction',$
      yrange=minmax(rms1), charsize=1.0,ytick_get=yrange
    plot, wave, rms2,/ystyle,title='rms/pixel,67th percentile, after residual subtraction',$
      yrange=minmax(yrange), charsize=1.0

;; these are the main ones - rms4 should be flat.
;; I don't quite understand why it's tailing away in the red.
;; if rms4 isn't yet flat, suggests you may need to increase the
;; maximum number of eigenspectra used (maxnrecon), but be a bit
;; careful of this.

    plot, wave, rms3,title='rms/pixel,error normalised,67th percentile, before',$
      yrange=minmax(rms3), charsize=1.0,ytick_get=yrange
    plot, wave, rms4,/ystyle,title='rms/pixel,error normalised,67th percentile, after',$
      yrange=minmax(yrange), charsize=1.0

    skysubflux = flux2 

    cleanplot, /silent
    dfpsclose
    spawn, ['gzip -f '+qaplotname], /sh

return, newdata
end

