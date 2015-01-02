pro phys130_f14_grades, alldata, test=test, sendit=sendit, final=final
; jm14feb02siena - parse the grades for this class 

    path = getenv('TEACHING_DIR')+'/130-F14/grades/'
;   path = getenv('TEACHING_DIR')+'/shared/Phys130/130-F14/moustakas/grades/'
    
    date = '14dec13' ; update this
    semester = 'Fall 2014'
    class = 'Physics 130 - General Physics I'

; read the grade spreadsheet downloaded from GoogleDocs
    gradefile = path+'phys130_grades_'+date+'.csv'
    data = read_gradefile(gradefile,unique_assignments=assign)

; specify the complete list of *possible* assignments and their
; relative weights
    allassign = ['Homework','Lab','In-Class Reading Quizzes','Midterm 1','Midterm 2','Final Exam']
    weight = [0.25,0.15,0.10,0.15,0.15,0.20]
    droplowest = [0,0,1,0,0,0]

;   keep = where(strmatch(data.last_name,'*Nap*'))
;   keep = where(data.final_exam gt 0.0)
;   data = data[keep]
    
    process_grades, data, assign=assign, allassign=allassign, $
      weight=weight, class=class, semester=semester, test=test, $
      sendit=sendit, alldata=alldata, final=final, droplowest=droplowest
    srt = sort(alldata.current_grade)
    struct_print, alldata[srt]

return
end
