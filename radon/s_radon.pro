pro s_radon,image,xcoords,ycoords,nparams,trans,L,param1= $
        param1,param2=param2,param3=param3


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;                                                            ;;
;; Purpose: Calculates the Radon transform of the given image.;;
;;                                                            ;;
;; Inputs: image - original image to be transformed           ;;
;;         xcoords,ycoords - coordinate values associated     ;;
;;           with each column,row of pixels                   ;;
;;         nparams - number of dimensions in transform space  ;; 
;;                                                            ;;
;; Outputs: trans - the Radon transform of image              ;;
;;          L - the number of parametrized points in image    ;;
;;            along each line in the transform space          ;;
;;                                                            ;;
;; Keywords: param1,param2,param3 - allowed values of the     ;;
;;             parameters in transform space; lines or proj-  ;;
;;             ection are either straight for nparams=2, or   ;;
;;             parabolic:                                     ;;
;;                  f(x)=param3*x^2+param2*x+param3           ;;
;;             for nparams=3.  If nparams=2 then param1 is    ;;
;;             thetas, param2 is rs                           ;;
;;                                                            ;;
;; Created: 01/04/08                                          ;;
;;                                                            ;;
;; Created by: Shaela Jones                                   ;;
;;                                                            ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; parameters
minpts=100



; input checks
if nparams ne 2 and nparams ne 3 then begin
  print,'Variable nparams must be either 2 or 3.'
  trans=-1
  return
endif
szxcoords=SIZE(xcoords)
szycoords=SIZE(ycoords)
if szxcoords[0] ne 1 then begin
  print,'Variable xcoords must be a 1D array'
  trans=-1
  return
endif
if szycoords[0] ne 1 then begin
  print,'Variable ycoords must be a 1D array'
  trans=-1
  return
endif
n1=N_ELEMENTS(param1)
n2=N_ELEMENTS(param2)
n3=N_ELEMENTS(param3)
if n1 eq 0 or n2 eq 0 then begin
  print,'Keywords param1 and param2 must be set.'
  trans=-1
  return
endif
if nparams eq 3 and n3 eq 0 then begin
  print,'Keyword param3 must be set.'
  trans=-1
  return
endif
locs=WHERE(xcoords[1:*]-xcoords le 0)  ; want strictly increasing x,y
if locs[0] ne -1 then begin
  print,'Problem with x coordinate vector'
  trans=-1
  return
endif
locs=WHERE(ycoords[1:*]-ycoords le 0)
if locs[0] ne -1 then begin
  print,'Problem with y coordinate vector'
  trans=-1
  return
endif
mn=MEAN(image)
if mn gt .01 then begin
  print,'Variable image should have zero mean,' $
            ,' running adj_cntrst on image.'
  nwimg=adj_cntrst(image,/zeromean)
  if MEAN(nwimg) gt .01 then stop
endif else nwimg=image



; initializations
nx=szxcoords[1]
ny=szycoords[1]
maxx=MAX(xcoords,min=minx)
maxy=MAX(ycoords,min=miny)
rexs=(xcoords-minx)/(maxx-minx)-0.5
reys=(ycoords-miny)/(maxy-miny)-0.5
xrange=[minx,maxx]
yrange=[miny,maxy]
if nx gt ny then deltr=MEAN(rexs[1:*]-rexs) else $
        deltr=MEAN(reys[1:*]-reys)
nts=SQRT(FLOAT(nx)^2+ny^2)



; main
if nparams eq 2 then begin
  L=INTARR(n1,n2)
  rad=FLTARR(n1,n2)
  t=param2
  for i1=0,n1-1 do begin
    thet=param1[i1]
    for i2=0,n2-1 do begin
      xy=LINEAR(thet,param2[i2],t)
      points=BILINEAR_UNEVEN(nwimg,rexs,reys,xy)
      list=WHERE(points ne 0.)
      if N_ELEMENTS(list) gt minpts then begin
        rad[i1,i2]=TOTAL(points)/N_ELEMENTS(list)
        L[i1,i2]=N_ELEMENTS(list)
      endif     ;else print,'nothing in range'
    endfor
  endfor
endif else begin   
  L=INTARR(n1,n2,n3+1)
  rad=FLTARR(n1,n2,n3+1)
;  t=(INDGEN(nts)+1.)/nts
  for i1=0,n1-1 do begin
    print,'i1= ',i1
    thet=param1[i1]
    cthet=COS(thet)
    sthet=SIN(thet)
    for i2=0,n2-1 do begin
      t=param2
      dont=0
      rho=param2[i2]
;      tempts=(INDGEN(100)+1.)/100.-0.5
;      tryxs=-cthet*rho+sthet*tempts
;      tryys=-sthet*rho-cthet*tempts
;      if i1 eq 0 then plot,tryxs,tryys,xrange=[-0.5,0], $
;          yrange=[0,0.5],ytitle='Altitude',xtitle='Time', $
;          title='Image Plane coverage' else oplot,tryxs,tryys
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
      if dont ne 1 then begin
;        oplot,[x1,x2],[y1,y2],psym=2
        maxacc=ABS(y2-y1)/(x1-x2)^2       ; max accel possible for given r,theta
        inds=WHERE(param3 le maxacc and param3 ge -maxacc)
        if inds[0] ne -1 then begin
          xy=LINEAR(thet,param2[i2],t)
          locs=WHERE(xy[0,*] gt -0.5 and xy[0,*] lt 0.5 $
                 and xy[1,*] gt -0.5 and xy[1,*] lt 0.5)
          num=N_ELEMENTS(locs)       
          t=(INDGEN(num)+1.)/num
          accel=param3[inds]
          for i3=0,N_ELEMENTS(inds)-1 do begin
            xy=PARAB_PTS(accel[i3],[x1,x2],[y1,y2],t)
;            oplot,xy[0,*],xy[1,*]
            points=BILINEAR_UNEVEN(nwimg,rexs,reys,xy)
            list=WHERE(points ne 0)
;stop
            if N_ELEMENTS(list) gt minpts then begin
              rad[i1,i2,inds[i3]]=TOTAL(points)/N_ELEMENTS(list)
              L[i1,i2,inds[i3]]=N_ELEMENTS(list)
            endif
          endfor
        endif else begin
          tlin=param2
          xy=LINEAR(thet,param2[i2],tlin)
          points=BILINEAR_UNEVEN(nwimg,rexs,reys,xy)
          list=WHERE(points ne 0.)
          if N_ELEMENTS(list) gt minpts then begin
            rad[i1,i2,n3]=TOTAL(points)/N_ELEMENTS(list)
            L[i1,i2,n3]=N_ELEMENTS(list)
          endif
        endelse
      endif
    endfor
  endfor
endelse



trans=TEMPORARY(rad)
save,filename='/home/shaela/Desktop/trans3.sav',trans,L,image, $
          xcoords,ycoords,param1,param2,param3,nwimg
;stop
return
end


; obsolete
