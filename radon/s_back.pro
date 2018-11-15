function s_back,trans,thetas,rs,xcoords,ycoords,accels

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;                                                            ;;
;; Purpose: Backprojects the radon transform into an image    ;;
;;          space whose spacing is given by the x,y coords.   ;;
;;                                                            ;;
;; Inputs: trans - the Radon transform to be backprojected    ;;
;;         thetas,rs - the theta,rho values corresponding to  ;;
;;           each column,row in the transform space           ;;
;;         xcoords,ycoords - x,y values corresponding to each ;;
;;           column,row in the image space                    ;;
;;                                                            ;;
;; Optional Inputs: accels - acceleration values correspond-  ;;
;;                    ing to the optional third dimension in  ;;
;;                    the transform space                     ;;
;;                                                            ;;
;; Outputs: back - the backprojected image                    ;;
;;                                                            ;;
;; Note: This function uses a very simple backprojection meth-;;
;;         od that does not try to solve all of the projec-   ;;
;;         tion equations simultaneously.                     ;;
;;                                                            ;;
;; Created: 01/24/08                                          ;;
;;                                                            ;;
;; Created by: Shaela Jones                                   ;;
;;                                                            ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; parameters

; input checks
sztrans=SIZE(trans)
if sztrans[0] ne 2 and sztrans[0] ne 3 then begin
  print,'Radon transform must be a 2D or 3D array.'
  return,-1
endif
nthet=N_ELEMENTS(thetas)
nr=N_ELEMENTS(rs)
if nthet ne sztrans[1] or nr ne sztrans[2] then begin
  print,'Variables thetas and rs must have the same length ', $
           'as transform dimensions 1 and 2.'
  return,-1
endif
accelson=0
if N_PARAMS() eq 6 then begin
  if sztrans[0] ne 3 then begin
    print,'Warning: accelerations provided with 2D transform!  '
    print,'Ignoring accelerations.'
  endif else begin
    accelson=1
    nacc=N_ELEMENTS(accels)
    if nacc ne sztrans[3] then begin
      print,'Variable accels must have length equal to third ', $
              'transform dimension length.'
      return,-1
    endif
  endelse
endif




; initializations
nx=N_ELEMENTS(xcoords)
ny=N_ELEMENTS(ycoords)
back=FLTARR(nx,ny)
maxx=MAX(xcoords,min=minx)
maxy=MAX(ycoords,min=miny)
rexs=(xcoords-minx)/(maxx-minx)-0.5
reys=(ycoords-miny)/(maxy-miny)-0.5
thetasp=thetas
rsp=rs
if accelson eq 1 then accelsp=accels
if nx gt ny then deltr=MEAN(rexs[1:*]-rexs) else $
        deltr=MEAN(reys[1:*]-reys)
nts=SQRT(nx^2+ny^2)



; main
for i1=0,nthet-1 do begin
  thet=thetasp[i1]
  if accelson eq 0 then begin
    t=rsp
    if MAX(trans[i1,*]) gt 0 then begin
      for i2=0,nr-1 do begin
        temp=FLTARR(nx,ny)
        r0=rsp[i2]
        xy=LINEAR(thet,r0,t)
        list=WHERE(xy[0,*] ge -0.5 and xy[0,*] lt 0.5 and $
                xy[1,*] ge -0.5 and xy[1,*] lt 0.5)
        if list[0] ne -1 then begin
          num=N_ELEMENTS(list)
          for point=0,num-1 do begin
            xp=xy[0,list[point]]
            yp=xy[1,list[point]]
            garbage=MIN(ABS(xp-rexs),nnxp)
            garbage=MIN(ABS(yp-reys),nnyp)
            temp[nnxp,nnyp]++
          endfor
          temp=temp*trans[i1,i2]/num
          back=back+temp
        endif
      endfor
    endif
  endif else begin
    t=(INDGEN(nts)+1.)/nts
    sthet=SIN(thet)
    cthet=COS(thet)
    for i2=0,nr-1 do begin
      dont=0
      acc=WHERE(trans[i1,i2,*] gt 0)
      if acc[0] ne -1 then begin       
        temp=FLTARR(nx,ny)
        rho=rsp[i2]
        x1=0.5*sthet/cthet-rho/cthet              ; x at y=-0.5
        y1=-0.5
        if x1 lt -0.5 then begin
          x1=-0.5
          y1=0.5*cthet/sthet-rho/sthet            ; y at x=-0.5
          if y1 lt -0.5 or y1 gt 0.5 then dont=1
        endif
        x2=-0.5*sthet/cthet-rho/cthet             ; x at y=0.5
        y2=0.5
        if x2 gt 0.5 then begin
          x2=0.5
          y2=-0.5*cthet/sthet-rho/sthet           ; y at x=0.5
          if y2 lt -0.5 or y2 gt 0.5 then dont=1
        endif
        maxacc=ABS(y2-y1)/(x1-x2)^2       ; max accel possible for given r,theta
        inds=WHERE(accels[acc] gt maxacc or accels[acc] lt -maxacc)
        if inds[0] ne -1 then begin
          print,'WTF?'
          stop
        endif
        if dont ne 1 then begin
          for i3=0,N_ELEMENTS(acc)-1 do begin
            xy=PARAB(accels[acc[i3]],[x1,x2],[y1,y2],t)
            list=WHERE(xy[0,*] gt -0.5 and xy[0,*] le 0.5 and $
                         xy[1,*] gt -0.5 and xy[1,*] le 0.5)
            if list[0] ne -1 then begin
              num=N_ELEMENTS(list)
              for point=0,num-1 do begin
                xp=xy[0,list[point]]
                yp=xy[1,list[point]]
                garbage=MIN(ABS(xp-rexs),nnxp)
                garbage=MIN(ABS(yp-reys),nnyp)
                temp[nnxp,nnyp]++
              endfor
              temp=temp*trans[i1,i2,acc[i3]]/num
              back=back+temp 
            endif else begin
              print,'WTF?'
              stop
            endelse
          endfor
        endif
      endif
    endfor
  endelse
endfor




;stop
return,back
end


; obsolete
