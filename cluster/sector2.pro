pro sector2,flnm,angmin,angmax,r_min,r_width,nsectors, $
              speeds,t0s,nomspeeds,allgood=allgood, $
              prob=prob

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;                                                            ;;
;; Purpose: Divides the area enclosed by angmin,angmax,       ;;
;;          r_width into sectors and runs cluster_analyze on  ;;
;;          each sector.                                      ;;
;;                                                            ;;
;; Inputs: flnm - name of diff file to be sectored            ;;
;;         angmin - min angle (degrees) of region of interest ;;
;;         angmax - max angle (degrees) of region of interest ;;
;;         r_min - minimum alt (pixels) of height-time ims    ;;
;;         r_width - range of alts (pixels) of height-time ims;;
;;         nsectors - number of angular bins to divide region ;;
;;           into                                             ;;
;;                                                            ;;
;; Outputs: speeds, t0s - speed, onset time array for all     ;;
;;           detected events, concatenated at each angular bin;; 
;;                                                            ;;
;; Dependencies: cluster_analyze                              ;;
;;                                                            ;;
;; Created: 11/02/07                                          ;;
;;                                                            ;;
;; Created by: Shaela Jones                                   ;;
;;                                                            ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


; parameters
nmean=14.


; input assurance
angmin=FLOAT(angmin)
angmax=FLOAT(angmax)
nsectors=FLOOR(nsectors)


;initializations
angrange=angmax-angmin
angwidth=angrange/FLOAT(nsectors)
maxlen=0



; put problysis to work
for i=0,nsectors-1 do begin
  print,'Processing sector ',i+1,' of ',nsectors,'.'
  angle=angmin+(i+0.5)*angwidth
  if KEYWORD_SET(prob) then begin
    if KEYWORD_SET(allgood) then begin
      problysis,flnm,angle,angwidth,r_min,r_width, $
              speeds,t0s,nomspeeds,/allgood
    endif else problysis,flnm,angle,angwidth,r_min,r_width, $
              speeds,t0s,nomspeeds
  endif else begin
    cluster_analyze,flnm,angle,angwidth,nmean,r_min,r_width, $
         speeds,t0s,/srad,/savefile
  endelse
  if t0s[0] ne -1 then begin
    nt0s=N_ELEMENTS(allt0s)
    if nt0s ne 0 then begin
      tmp=allt0s
      length=N_ELEMENTS(t0s)
      if length le maxlen then begin
        allt0s=STRARR(i+1,maxlen)
        allt0s[0:i-1,*]=tmp
        allt0s[i,0:length-1]=t0s
        tmp=allspeeds
        allspeeds=FLTARR(i+1,maxlen)
        allspeeds[0:i-1,*]=tmp
        allspeeds[i,0:length-1]=speeds
      endif else begin
        allt0s=STRARR(i+1,length)
        allt0s[0:i-1,0:maxlen-1]=tmp
        allt0s[i,*]=t0s
        tmp=allspeeds
        allspeeds=FLTARR(i+1,length)
        allspeeds[0:i-1,0:maxlen-1]=tmp
        allspeeds[i,*]=speeds
        maxlen=length
      endelse
;        allt0s=STRARR(nt0s+length)
;        allt0s[0:nt0s-1]=tmp
;        allt0s[nt0s:*]=t0s
;        tmp=allspeeds
;        allspeeds=FLTARR(nt0s+length)
;        allspeeds[0:nt0s-1]=tmp
;        allspeeds[nt0s:*]=speeds
    endif else begin
      allspeeds=speeds
      allt0s=t0s
      maxlen=N_ELEMENTS(t0s)
    endelse
  endif
endfor
print,'maxlen= ',maxlen
if N_ELEMENTS(allt0s) ne 0 then begin
  speeds=TEMPORARY(allspeeds)
  t0s=TEMPORARY(allt0s)
endif else begin
  speeds=-1
  t0s=-1
endelse



save,filename='/home/shaela/Desktop/jic0504.sav',speeds, $
        t0s,nomspeeds,flnm,angmin,angmax,r_min,r_width,nsectors
end




; obsolete

;  if N_ELEMENTS(t0s) eq 0 OR N_ELEMENTS(speeds) eq 0 then begin
;    print,'Trouble with problysis result.'
;    stop
;  endif
