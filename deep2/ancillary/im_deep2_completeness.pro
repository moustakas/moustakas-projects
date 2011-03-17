;+
; NAME:
;   deep2_completeness
; REVISION HISTORY:
;   2005-04-07 MRB, NYU
;-
;------------------------------------------------------------------------------
pro im_deep2_completeness

pcatm= mrdfits(deep2_path(/analysis)+'deep2_pcat.fits', 1)
zcat=mrdfits(deep2_path(/analysis)+'zcat.dr3.uniq.good.fits.gz',1)

; apply color/magnitude cuts and star-galaxy cut of 0.2 

ikeep=where(deep2_color_cuts(pcatm.magb, pcatm.magr, pcatm.magi) gt 0 AND $
            pcatm.pgal gt 0.2, nkeep)
pcatm=pcatm[ikeep]

spherematch, zcat.ra, zcat.dec, pcatm.ra, pcatm.dec, 1./3600., m1, m2

; NUM is the number of grid points for completeness as a function of R
; and R-I; initialize the output arrays

num=15L 
rmilim=[0., 2.]
rlim=[19., 24.1]
npcat=fltarr(num,num)
nattempt=fltarr(num,num)
ngot=fltarr(num,num)

; photometric catalog

rmi=pcatm.magr-pcatm.magi
r=pcatm.magr
ii=where(rmi gt rmilim[0] and $
         rmi lt rmilim[1] and $
         r gt rlim[0] and $
         r lt rlim[1])
irmi=long((rmi[ii]-rmilim[0])/(rmilim[1]-rmilim[0])*float(num)) ; R-I index number ([0-14])
ir=long((r[ii]-rlim[0])/(rlim[1]-rlim[0])*float(num))           ; R-band index number ([0-14])
igrid=ir+irmi*num
isort=sort(igrid)
iuniq=uniq(igrid[isort])
indx=igrid[isort[iuniq]]
nindx=iuniq-[-1L, iuniq]
npcat[indx]=float(nindx)

; full spectroscopic catalog (no redshift quality cut)

rmi=zcat.magr-zcat.magi
r=zcat.magr
ii=where(rmi gt rmilim[0] and $
         rmi lt rmilim[1] and $
         r gt rlim[0] and $
         r lt rlim[1])
irmi=long((rmi[ii]-rmilim[0])/(rmilim[1]-rmilim[0])*float(num))
ir=long((r[ii]-rlim[0])/(rlim[1]-rlim[0])*float(num))
igrid=ir+irmi*num
isort=sort(igrid)
iuniq=uniq(igrid[isort])
indx=igrid[isort[iuniq]]
nindx=iuniq-[-1L, iuniq]
nattempt[indx]=float(nindx)

; "successful" spectroscopic catalog (redshift Q>3)

rmi=zcat.magr-zcat.magi
r=zcat.magr
ii=where(rmi gt rmilim[0] and $
         rmi lt rmilim[1] and $
         r gt rlim[0] and $
         r lt rlim[1] and $
         zcat.zquality ge 3) ; note!
irmi=long((rmi[ii]-rmilim[0])/(rmilim[1]-rmilim[0])*float(num))
ir=long((r[ii]-rlim[0])/(rlim[1]-rlim[0])*float(num))
igrid=ir+irmi*num
isort=sort(igrid)
iuniq=uniq(igrid[isort])
indx=igrid[isort[iuniq]]
nindx=iuniq-[-1L, iuniq]
ngot[indx]=float(nindx)

; FATTEMPT is given by the fraction of photometric targets that
; were placed on a slit

fattempt=((nattempt+float(npcat eq 0.))/(npcat>1.) ) < 1.
fgot=((ngot+float(nattempt eq 0.))/(nattempt>1.)) < 1.
outfile=deep2_path(/analysis)+'deep2_completeness_full.fits'
mwrfits, fattempt, outfile, /create
mwrfits, fgot, outfile

; FGOT is given by the fraction of photometric targets for which a
; successful redshift (Q>3) was obtained; following Willmer et
; al. (2006), blue galaxies without successful redshifts are at z>1.4,
; so set FGOT=1 for R-I<0.9

rmidiv=0.9
irmidiv=long((rmidiv-rmilim[0])/(rmilim[1]-rmilim[0])*float(num))
fgot[*,0:irmidiv]=1.
outfile = deep2_path(/analysis)+'deep2_completeness.fits'
mwrfits, fattempt, outfile, /create
mwrfits, fgot, outfile

stop

return
end

