function pmfspecid, s, plate=plate, mjd=mjd, fiberid = fiberid

    if (n_elements(plate) eq 0L) then plate = s.plateid
    if (n_elements(mjd) eq 0L) then mjd = s.mjd
    if (n_elements(fiberid) eq 0L) then fiberid = s.fiberid

    sid = long64(string(plate, format='(I4.4)') + $
      string(mjd, format='(I5)') + $
      string(fiberid, format='(I3.3)'))

return, sid
end
