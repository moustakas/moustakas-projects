; jm00june5ucb

pro trgb_manualpsf, image

        masterdir = rdtxt('master.dir')
        masterdir = masterdir[0]

; parameter file

        spawn, ['find '+masterdir+' -name "'+image+'.pars" -print'], parfile
        nold = where(strpos(parfile,'OLD') eq -1,count)
        if count ne 0 then parfile = parfile[nold]
        parfile = parfile[0]

; psf header        
        
        spawn, ['find '+masterdir+' -name "'+image+'.psf.fits" -print'], irafpsf
        nold = where(strpos(irafpsf,'OLD') eq -1,count)
        if count ne 0 then irafpsf = irafpsf[nold]
        irafpsf = irafpsf[0]

        psfhead = headfits(irafpsf)

        trgb_makeopt, parfile, psfhead

;       setenv, 'imdir='+masterdir
;       spawn, ['alias imdir cd $imdir']
        
return
end
