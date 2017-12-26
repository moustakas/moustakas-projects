pro phys310_f17_grades, alldata, test=test, sendit=sendit, final=final, drop=drop
; jm15sep20 - parse the grades for this class 

    path = getenv('TEACHING_DIR')+'/310-F17/grades/'
    
    date = '17dec17' ; update this
    semester = 'Fall 2017'
    class = 'PHYS310 - Mechanics I'

; read the grade spreadsheet downloaded from GoogleDocs
    gradefile = path+'phys310_grades_'+date+'.csv'
    data = read_gradefile(gradefile,unique_assignments=assign)

; specify the complete list of *possible* assignments and their
; relative weights
    allassign = ['Homework','Lab','Midterm','Final']
    weight = [0.20,0.20,0.45,0.15]
    if keyword_set(drop) then begin
       droplowest = [1,0,1,0]
    endif else begin
       droplowest = [0,0,0,0]
    endelse

    data.final *= 1.1

;   keep = where(data.final_exam gt 0.0)
;   keep = where(strmatch(data.first_name,'**'))
;   data = data[keep]

    process_grades, data, assign=assign, allassign=allassign, $
      weight=weight, class=class, semester=semester, test=test, $
      sendit=sendit, alldata=alldata, final=final, droplowest=droplowest
    srt = sort(alldata.current_grade)
    struct_print, alldata[srt]

return
end
