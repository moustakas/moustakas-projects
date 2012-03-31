pro envsfrm_zgroup, mlimit=mlimit, mproj=mproj, mz=mz

;; lf_distmod(0.05)=35.96
;; 17.77 - 35.96 = -18.19 +5.*alog10(h10) = ~-19.0
;;  ~meanseparation (V/N)^0.3333 or (V/N/4)^0.3333
;;  V ~ 40Mpc^3 * 1/6, N with Mr < -18 = 74668.
;;  meanseparation ~ 3-4 Mpc
;;  linking length 0.2 gives 0.6-0.8 Mpc, in redshift units: 70/3E+5 * 0.6-0.8 ~ 0.0002
;;  linking length 0.05 gives 0.15-0.2 Mpc, in redshift units: 70/3E+5 * 0.15 -0.2 ~ 0.00005

    if (not keyword_set(mlimit)) then mlimit = -19.
    if (not keyword_set(mproj)) then mproj = 1.0E-4
    if (not keyword_set(mz)) then mz = 3.5E-4

    infile0 = getenv('IM_DATA_DIR')+'/nsa/nsa_v0_1_2.fits.gz'
    outfile = './envsfrm_zgroup_all_group_mlimit'+string(abs(mlimit), format='(f4.1)') + $
      '_mp' + string(mproj, format='(f7.5)') + '_mz'+ string(mz, format='(f7.5)') + '.fits'

    nsa = mrdfits(infile0, 1) ; ,range=[0,1000])
    group0 = {ra:0.0D, dec:0.0D, zdist:0.0, mult:0L, first:0L, $
              x:0.0D, y:0.0D, z:0.0D }
    member0 = {ra:0.0D, dec:0.0D, objno:0L, zdist:0.0, next:0L, ing:0L, $
               mult:0L, first:0L, $
               neargroup2:0L, distgroup2:0.0D, $
               neargroup3:0L, distgroup3:0.0D, $
               neargroup5:0L, distgroup5:0.0D, $
               neargroup8:0L, distgroup8:0.0D, $
               neargroup10:0L, distgroup10:0.0D, $
               ellipneargroup2:0L, ellipdistgroup2:0.0D, ellipneargroupmult2:0.0D, $
               ellipneargroup3:0L, ellipdistgroup3:0.0D, $
               ellipneargroup5:0L, ellipdistgroup5:0.0D, $
               ellipneargroup8:0L, ellipdistgroup8:0.0D, $
               ellipneargroup10:0L, ellipdistgroup10:0.0D, $
               x:0.0D, y:0.0D, z:0.0D $
               }

; run group finder, FoF
    ii = where(nsa.absmag[4] le mlimit, nn)
    ing = zgroup(nsa[ii].zdist, nsa[ii].ra, nsa[ii].dec, mproj, mz,  $
            mult=mult, first=first, next=next)

    group = replicate(group0, n_elements(mult)) 
    group.mult = mult
    group.first = first
    member = replicate(member0, n_elements(nsa))
    member.ra = nsa.ra
    member.dec = nsa.dec
    member.objno = nsa.objno
    member.zdist = nsa.zdist
    member[ii].next = next
    member[ii].ing = ing
    member[ii].mult = mult[ing]
    member[ii].first = first[ing]

;; convert to redshif unit
    print, "convert to redshif unit"
    for i=0L, n_elements(member)-1L do begin
        tmp_ra = member[i].ra
        tmp_dec = member[i].dec
        tmp_zdist = member[i].zdist
        tmp_x = angles_to_x(tmp_ra, tmp_dec)
        member[i].x = tmp_x[0]*tmp_zdist
        member[i].y = tmp_x[1]*tmp_zdist
        member[i].z = tmp_x[2]*tmp_zdist
    endfor
;; calculate group centers
    print, "calculate group centers."
    for i=0L, n_elements(group)-1 do begin
        tt = where(member.ing eq i, ss)
        if (ss gt 1) then begin
            tmp_ra = mean(member[tt].ra)
            tmp_dec = mean(member[tt].dec)
            tmp_zdist = mean(member[tt].zdist)
        endif else begin
            tmp_ra = member[tt].ra
            tmp_dec = member[tt].dec
            tmp_zdist = member[tt].zdist
        endelse
        group[i].ra = tmp_ra
        group[i].dec = tmp_dec
        group[i].zdist = tmp_zdist
        tmp_x = angles_to_x(tmp_ra, tmp_dec)
        group[i].x = tmp_x[0]*tmp_zdist
        group[i].y = tmp_x[1]*tmp_zdist
        group[i].z = tmp_x[2]*tmp_zdist
    endfor

;; calculate distance from the closest group (2 members, 3, 5, 8, 10)
    print, "calculate distance from the closest group (2 members, 3, 5, 8, 10)"
    for i=0L, 4L do begin
        print, 'nmembers=', i
        case i of
            0: nmember=3
            1: nmember=5
            2: nmember=8
            3: nmember=10
            4: nmember=2
        endcase
        ghigh = where(group.mult ge nmember, nghigh)
        for j=0L, n_elements(member)-1 do begin
;           print, j, n_elements(member)-1
            all_dist = (member[j].x-group[ghigh].x)^2 + $
                       (member[j].y-group[ghigh].y)^2 + $
                       (member[j].z-group[ghigh].z)^2
            dist_min = min(all_dist, imin)
;; find parallel difference, dot product
            midx = (member[j].x+group[ghigh].x)*0.5
            midy = (member[j].y+group[ghigh].y)*0.5
            midz = (member[j].z+group[ghigh].z)*0.5
            midnorm = sqrt(midx^2+midy^2+midz^2)
            middirx = midx/midnorm
            middiry = midy/midnorm
            middirz = midz/midnorm
            diffx = member[j].x - group[ghigh].x
            diffy = member[j].y - group[ghigh].y
            diffz = member[j].z - group[ghigh].z
            diffpar = abs(middirx*diffx+middiry*diffy+middirz*diffz)
;; find projected difference
            diffproj = sqrt(diffx^2+diffy^2+diffz^2 - diffpar^2)

            ellipdist_min = min((diffpar/12.E-4)^2+(diffproj/2.E-4)^2, ellipimin)
            ellipdist_min2 = min((diffpar/8.E-4)^2+(diffproj/2.E-4)^2, ellipimin2)

            case i of 
                0: begin
                   member[j].neargroup3 = ghigh[imin]
                   member[j].distgroup3 = sqrt(dist_min)
                   member[j].ellipneargroup3 = ghigh[ellipimin]
                   member[j].ellipdistgroup3 = sqrt(ellipdist_min)
                   end
                1: begin
                   member[j].neargroup5 = ghigh[imin]
                   member[j].distgroup5 = sqrt(dist_min)
                   member[j].ellipneargroup5 = ghigh[ellipimin]
                   member[j].ellipdistgroup5 = sqrt(ellipdist_min)
                   end
                2: begin
                   member[j].neargroup8 = ghigh[imin]
                   member[j].distgroup8 = sqrt(dist_min)
                   member[j].ellipneargroup8 = ghigh[ellipimin]
                   member[j].ellipdistgroup8 = sqrt(ellipdist_min)
                   end
                3: begin
                   member[j].neargroup10 = ghigh[imin]
                   member[j].distgroup10 = sqrt(dist_min)
                   member[j].ellipneargroup10 = ghigh[ellipimin]
                   member[j].ellipdistgroup10 = sqrt(ellipdist_min)
                   end
                4: begin
                   member[j].neargroup2 = ghigh[imin]
                   member[j].distgroup2 = sqrt(dist_min)
                   member[j].ellipneargroup2 = ghigh[ellipimin2]
                   member[j].ellipdistgroup2 = sqrt(ellipdist_min2)
                   member[j].ellipneargroupmult2 = group[ghigh[ellipimin2]].mult
                   end
           endcase
       endfor
    endfor

    im_mwrfits, member, outfile, /nogzip, /clobber
    im_mwrfits, group, outfile, /append, /gzip

    stop
    ii = where(member.mult gt 50L, nn)
    jj = where(member.mult eq 1L, mm)
    print, nn, mm
    djs_plot, member.ra, member.dec, yra=[-2,65], xra=[100,300], psym=3
    djs_oplot, member[ii].ra, member[ii].dec, color='red', psym=3
;   djs_oplot, member[jj].ra, member[jj].dec, color='green', psym=3
    djs_plot, member.ra, member.dec, yra=[25,35], xra=[190,205], psym=3
    djs_oplot, member[ii].ra, member[ii].dec, color='red', psym=3

return
end
