pro trans,file,rmin,rwidth,ang,angwidth,nmean,image,thetas,rhos, $
         fun2=fun2,fun3=fun3,fun4=fun4,fun5=fun5,fun6=fun6, $
         fun7=fun7,srad=srad

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;                                                            ;;
;; Purpose: Makes a height-time image from the images in file ;;
;;          and uses the Radon Transform to enhance linear    ;;
;;          features in one of several ways.                  ;;
;;                                                            ;;
;; Inputs: file - name of diff file to create height-time     ;;
;;           image for                                        ;;
;;         rmin,rwidth - min alt, range of alts (pixels) of   ;;
;;           height-time image                                ;;
;;         ang,angwidth - min angle, size of angular bin (in  ;;
;;           degrees) of height-time images                   ;;
;;         nmean - multiple of mean value at which to thresh- ;;
;;           old transform                                    ;;
;;                                                            ;;
;; Optional Inputs: image - original height-time image        ;;
;;                  thetas,rhos - radon transform parameters  ;;
;;                                                            ;;
;; Keywords: fun2 - radon transform thresholded at median val ;;
;;           fun3 - fun2 thresholded at nmean*mean(fun2)      ;;
;;           fun4 - fun3 closed using morph_close idl function;;
;;           fun5 - fun4 smoothed                             ;;
;;           fun6 - fun3 closed using morph_close idl function;;
;;             with gray keyword set                          ;;
;;           fun7 - fun6 smoothed                             ;;
;;           thetas,rhos - radon transform                    ;;
;;                                                            ;;
;; Created: 10/12/07                                          ;;
;;                                                            ;;
;; Created by: Shaela Jones                                   ;;
;;                                                            ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


; parameters
resize=4
sizestruc=3
outsz=1024



; get original image
image=' '
r_strp,resize,rmin,rwidth,ang,angwidth,diffs_file=file,radpic=image
szimage=SIZE(image)
if szimage[0] le 1 then begin
  print,'Trouble with r_strp result.'
  return
endif



;get info for transform
ave=MEAN(image,/nan)
image=image-ave




; get transform difference image
if keyword_set(srad) then begin
  restore,file
  diffs=0
  axes,dttms,rmin,rwidth,filelist[0],ctr,ang,outsz,tms,alts
  filelist=0
  dttms=0
  find_params,image,tms,alts,thetas,rhos,as
  s_radon,image,tms,alts,2,try,L,param1=thetas,param2=rhos
;  s_radon,flat,tms,alts,thetas,rhos,tryflat
endif else begin 
  nthets=CEIL(!dpi*SQRT((szimage[1]^2+szimage[2]^2)/2))/2
  thetas=(FINDGEN(nthets)+1)*!dpi/2/nthets+!dpi/2
  thetas=find_thetas(image)
  nthets=N_ELEMENTS(thetas)
  try=radon(image,/double,theta=thetas,rho=rhos,nrho=nthets, $
         ntheta=nthets)
;  tryflat=radon(flat,/double,theta=thetas,rho=rhos,nrho= $
;         nthets,ntheta=nthets)
endelse




; threshold and close image
threshold,try,nmean,sizestruc,fun2=fun2,fun3=fun3,fun4=fun4, $
     fun5=fun5,fun6=fun6,fun7=fun7




end


; obsolete

;nthets=CEIL(!dpi*SQRT((szimage[1]^2+szimage[2]^2)/2))/2
;thetas=(FINDGEN(nthets)+1)*!dpi/2/nthets+!dpi/2
