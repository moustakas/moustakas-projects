pro get_bcgs_photometry
; jm10jun05ucsd - retrieve the multiband photometry for
;   Anthony's BCG sample

    bcgspath = ages_path(/projects)+'bcgs/'
    datapath = getenv('DATA_DIR')+'/bootes/'
    ss = rsex(bcgspath+'bcgs_sample_v3.sex')
    iband = mrdfits(datapath+'bootes_I.fits.gz',1)

    spherematch, iband.alpha_j2000, iband.delta_j2000, $
      ss.ra, ss.dec, 3.0/3600.0, m1, m2
;   m1 = im_spherematch(irac,ss,ratagname1='alpha_j2000',dectagname1='delta_j2000',$
;     match2=m2,radius=10.0,raoff=raoff,decoff=decoff)
    srt = sort(m2) & m1 = m1[srt] & m2 = m2[srt]
    help, ss, m1
    
    bands = ['Bw','R','I','J','H','Ks','ch1','ch2','ch3','ch4']
    for jj = 0, n_elements(bands)-1 do begin
       phot1 = mrdfits(datapath+'bootes_'+bands[jj]+$
         '.fits.gz',rows=m1,1)
       tags = tag_names(phot1)
       phot1 = im_struct_trimtags(phot1,select=tags,$
         newtags=bands[jj]+'_'+tags)
       if (jj eq 0) then phot = phot1 else phot = $
         struct_addtags(temporary(phot),phot1)
    endfor

    im_mwrfits, phot, bcgspath+'bcgs_photometry_v3.fits', /clobber
    
return
end
