;+
; NAME:   
;            edge_detect
;
; PURPOSE: 
;    To detect via maximum likelihood analysis the edge of the 
;    TRGB, based on 
;
; CATEGORY:
;    Analysis routine
;
;
; CALLING SEQUENCE:
;    edge_detect, maglist, lnlikelihood
;
; INPUTS: 
;    maglist -- list of magnitudes for input stars 
;
; OUTPUT:
;    lnlike -- ln of likelihood curve as a function of m_c, the cutoff
;                  magnitude
;
; KEYWORDS:
;    maxmag -- maximum magnitude for search
;    minmag -- minimum magnitude for search
;
; MODIFICATION HISTORY:
; md- 21Jun00
;-



pro edge_detect, maglist, lnlike, minmag=minmag, maxmag=maxmag

; first input data

if not Keyword_set(minmag) then minmag=20.
if not Keyword_set(maxmag) then maxmag=25.
 

bmag=fix((maglist-minmag)*100.) ;binned magnitudes in range minmag,maxmag
numberbins=100*(maxmag-minmag)
ii=where(bmag ge 0 and bmag lt  numberbins)
bmag=bmag(ii)  ; subset of stars in selected magnitude range
Ntotal=n_elements(bmag)
lnlike =fltarr(50,41)
print, Ntotal,' stars within range', minmag, ' to ',maxmag

for k= 0, 40 do begin
cutfactor=k/20.
for i= 0, 49 do begin           ;loop over cutoff magnitude
    mc=i/10.  ; magnitude cutoff in range minmag, minmag+5
    g=form_gg(mc,cutfactor)     ;note that this produces 500 points
    lnlike[i,k] = -Ntotal*alog(total(g[0:numberbins-1])) + total(alog(g[bmag]))
end     
end

lnlike=lnlike-max(lnlike)
;print, lnlike
return
end











