; zptfile = '/Users/dey/DECaLS/ZeroPoints/decals-zpt-all-2015july13.fits'
zptfile = 'decals-zpt-all-2015oct19.fits'
;zptfile = '/Users/dey/DECaLS/ZeroPoints/decals-zpt-all-2015oct19.fits'
;savefile='decals-obs-all-2015oct19.sav'
savefitsfile='decals-obs-2015oct19.fits'

z = mrdfits(zptfile,1,/silent)
nz = n_elements(z.ra)

; Select g/r/z photometric data with good photometry
ipho = where(abs(z.zpt-z.ccdzpt) lt 0.05 $
                and z.ccdnmatch ge 20 and z.ccdphrms le 0.2 $
		and z.ccdnum ne 31 and z.dec ge -20 $
		and ((z.filter eq 'g' and abs(z.ccdzpt-25.08) le 0.20) $
		  or (z.filter eq 'r' and abs(z.ccdzpt-25.29) le 0.15) $
		  or (z.filter eq 'z' and abs(z.ccdzpt-24.92) le 0.15)) $
                ,np)
photflag = intarr(nz)
photflag[ipho] = 1

ig = where(z.filter eq 'g' and z.avsky gt 0 and z.ccdzpt lt 30 and z.ccdzpt gt 0 and photflag eq 1)
ir = where(z.filter eq 'r' and z.avsky gt 0 and z.ccdzpt lt 30 and z.ccdzpt gt 0 and photflag eq 1)
iz = where(z.filter eq 'z' and z.avsky gt 0 and z.ccdzpt lt 30 and z.ccdzpt gt 0 and photflag eq 1)

; strip the filenames of the path to the data
fz = z.filename
n = n_elements(fz)
fznew = strarr(n)
for j=0,n-1 do begin
        fj = fz[j]
        strs = strsplit(fj,'/',/extract)
        fj2 = strs[-1]
        fznew[j] = fj2
endfor
z.filename = fznew
files = z.filename

ff = files[sort(files)]
ff = ff[uniq(ff)]
nf = n_elements(ff)
print,'Number of files=',nf
next = 62

ct4mlat = -(30.d0+10./60.+10.78/3600.d0)
ct4mlon = -(70.d0+48./60.+23.49/3600.d0)
ct4malt = 2206.8

nf = n_elements(z.ccdra)
filter = z.filter
ra = z.ccdra
dec = z.ccddec
mjd = z.mjd_obs
jd = mjd + 2400000.5d0
airmass = z.airmass

; Calculate zenith distance
eq2hor,ra,dec,jd,alt,az,lat=ct4mlat,lon=ct4mlon,altitude=ct4malt
zd=90.-alt

; Calculate the sun position
sunpos,jd,sunra,sundec
eq2hor,sunra,sundec,jd,altsun,azsun,lat=ct4mlat,lon=ct4mlon,altitude=ct4malt
sunsep = sphdist(ra,dec,sunra,sundec,/degrees)

; Calculate the moon position, phase and separation
moonpos,jd,moonra,moondec
mphase,jd,moonphase
eq2hor,moonra,moondec,jd,altmoon,azmoon,lat=ct4mlat,lon=ct4mlon,altitude=ct4malt
moonsep = sphdist(ra,dec,moonra,moondec,/degrees)
; ... and whether it is up
moonup = intarr(nf)
ii = where(altmoon gt 0)
moonup[ii] = 1

exptime = z.exptime
expnum = z.expnum
zpt = z.ccdzpt
seeing = z.seeing
ccdnum = z.ccdnum
ccdname = z.ccdname

; Calculate sky brightness and 5-sigma point source depth
skybr=-2.5*alog10(z.avsky/0.262^2/exptime)+zpt
skysig = z.ccdskyrms
npix = !pi*(1.35*seeing/2./0.262)^2
depth= -2.5*alog10(5.*skysig*sqrt(npix)/exptime) + zpt

; create the structure into which the header stuff will go
ccddat = create_struct( $
        'FILENAME',files, $
        'EXPNUM',expnum, $
        'EXPTIME',exptime, $
	'ZD',zd, $
        'FILTER',filter, $
        'SEEING',seeing, $
        'RA',ra, $
        'DEC',dec, $
	'DATE',jd, $
        'CCDNUM',ccdnum, $
        'CCDNAME',ccdname, $
        'CCDZPT',zpt, $
	'SKYBR',skybr, $
	'SKYSIG',skysig, $
	'DEPTH',depth, $
	'SUNSEP',sunsep, $
	'MOONSEP',moonsep, $
	'MOONPHASE',moonphase, $
	'MOONUP',moonup, $
	'PHOTFLAG',photflag)

; write it out
;save,ccddat,filename=savefile
mwrfits,ccddat,savefitsfile,/create

; make plots

glim=24.77-2.5*alog10(sqrt(3.))
rlim=24.07-2.5*alog10(sqrt(3.))
zlim=23.10-2.5*alog10(sqrt(3.))

gmed =median(depth[ig])
rmed =median(depth[ir])
zmed =median(depth[iz])

yg=histogram(depth[ig],min=18,max=30,binsiz=0.1,loca=xx)
yr=histogram(depth[ir],min=18,max=30,binsiz=0.1,loca=xx)
yz=histogram(depth[iz],min=18,max=30,binsiz=0.1,loca=xx)

ps_open,'Depths_DR2',/ps,/encap
device,/schoolbook,xsiz=10,ysiz=5,/inch

!p.multi=[0,3,1,0]

plot,xx,yg,psym=10,xr=[21.9,26.1],/xsty,/ysty,title='g-band' $
	,charsize=2,thick=2,charthick=2
oplot,[glim,glim],[0,100000],lines=2
oplot,[gmed,gmed],[0,100000]

plot,xx,yr,psym=10,xr=[20.9,25.1],/xsty,/ysty,title='r-band' $
	,xtitle='5sigma Point Source Depth (apdia = 1.35xFWHM)' $
	,charsize=2,thick=2,charthick=2
oplot,[rlim,rlim],[0,100000],lines=2
oplot,[rmed,rmed],[0,100000]

plot,xx,yz,psym=10,xr=[19.9,24.1],/xsty,/ysty,title='z-band' $
	,charsize=2,thick=2,charthick=2
oplot,[zlim,zlim],[0,100000],lines=2
oplot,[zmed,zmed],[0,100000]

!p.multi=0
ps_close

stop
end

stop
end
