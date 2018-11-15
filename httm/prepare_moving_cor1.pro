pro prepare_moving_cor1,date,spcrft,windw,file, $
         quarter=quarter,files=files,fivetwelve=fivetwelve, $
         entire=entire,rotate=rotate,sigfilt=sigfilt,nodiff= $
         nodiff,calibrate=calibrate,no_bkg=no_bkg, $
         smoother=smoother,mass=mass,rescale=rescale

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;                                                            ;;
;; Purpose: Makes difference images who base image is a       ;;
;;           moving minimum of the local (in time) images.    ;;
;;                                                            ;;
;; Inputs: date - date of images to be processed              ;;
;;         spcrft - spacecraft (a or b) of interest           ;;
;;         windw - number of images to combine to make base   ;;
;;                 image (must be even)                       ;;
;;                                                            ;;
;;; Optional Output: file - name of file to which difference   ;;
;;         movie was saved                                    ;;
;;                                                            ;;
;; Keywords: quarter - only processes the specified (1-4)     ;;
;;                     quarter of the available images for    ;;
;;                     date.                                  ;;
;;           files - list of files to be processed by name    ;;
;;                   instead of by date or quarter            ;;
;;           fivetwelve - bins images to 512x512 resolution   ;;
;;           entire - uses entire filelist to make base image ;;
;;                    (ignores windw variable)                ;;
;;           rotate - sets rotate keyword when calling secchi ;;
;;                    prep; images are prepped with solar     ;;
;;                    north up                                ;;
;;           sigfilt - uses sigma_filter to remove stars,     ;;
;;                     cosmic rays                            ;;
;;           nodiff - make a regular image sequence; no diff- ;;
;;                    erencing                                ;;
;;           calibrate - if set calfac_off and calimg_off     ;;
;;                    secchi_prep keywords are not used       ;;
;;           no_bkg - if set the monthly median background is ;;
;;                    not subtracted by secchi_prep           ;;
;;           smoother - smooths each pixel's value as a funct-;;
;;                    ion of time                             ;;
;;           mass - converts images to mass in units of g     ;;
;;           rescale - scale image intensity between param-   ;;
;;                     eters cor1max and cor1min              ;;
;;                                                            ;;
;; Created: 08/09/04                                          ;;
;;                                                            ;;
;; Modified: 11/20/07 - modified to work with get_flnms2,     ;;
;;                      removed obsolete keyword bckgrnd, set ;;
;;                      rotate_on keyword to secchi_prep      ;;
;;           04/13/08 - modified to include nodiff keyword    ;;
;;                                                            ;;
;; Created by: Shaela Jones                                   ;;
;;                                                            ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


; parameters
cor1max=60.0
cor1min=-10.0



; input checking
spcrft=STRLOWCASE(spcrft)
if (keyword_set(quarter)) and (keyword_set(files)) then begin
  print,'files and quarter keywords both set; setting ',$
         'quarter to zero.'
  quarter=0
endif
if not(keyword_set(entire)) then begin
  if (windw mod 2) ne 0 then begin
    windw++
    print,'changing windw to even number: ',windw
  endif
endif




; initializations
outsz= (keyword_set(fivetwelve)) ? 512:1024
rot_on= (keyword_set(rotate)) ? 1:0
calib= KEYWORD_SET(calibrate) ? 0:1
back_off= KEYWORD_SET(no_bkg) ? 1:0
instrument='COR1'



; get filenames
if N_ELEMENTS(files) eq 0 then begin
  filesset=0
  prob=0
  get_flnms2,date,spcrft,instrument,files
  if N_ELEMENTS(files) eq 1 then begin
    print,'Problem with get_flnms2 result.'
    return
  endif
endif else begin
  filesset=1
  print,'Warning: Ends of window not trustworthy.'
endelse
szfiles=SIZE(files)
nframes=szfiles[2]
;stop



; find files during (quarter of) day, previous and subsequent
if keyword_set(quarter) then begin
  ends=FLTARR(5)
  ends[0]=-1
  for i=1,3 do begin
    ends[i]=i*n/4-1
    span=ends[i]-ends[i-1]
    if span mod 3 ne 0 then ends[i]=ends[i]-(span mod 3)
  endfor
  ends[4]=n-1
  beg=ends[quarter-1]+1
  en=ends[quarter]
  if quarter eq 1 then alpha='A' $
  else if quarter eq 2 then alpha='B' $
  else if quarter eq 3 then alpha='C' $
  else if quarter eq 4 then alpha='D'
  if NOT(keyword_set(entire)) then begin
    if alpha ne 'A' then previous=files[*,beg-windw/2-1:beg-1] $
    else begin
      get_flnms2,YESTERDAY(date),spcrft,instrument,oldfiles
      if oldfiles[0] eq -1 then begin
        print,"Problem with yesterday's get_flnms2 result."
        return
      endif
      szold=SIZE(oldfiles)
      previous=oldfiles[*,szold[2]-windw/2-1:szold[2]-1]
    endelse
    if alpha ne 'D' then subseq=files[*,en+1:en+windw/2] $
    else begin
      get_flnms2,TOMORROW(date),spcrft,instrument,nextfiles
      if nextfiles[0] eq -1 then begin
        print,"Problem with tomorrow's get_flnms2 result."
        stop
      endif
      subseq=nextfiles[*,0:windw/2-1]
    endelse
  endif
  files=files[*,beg:en]
endif else if NOT(KEYWORD_SET(nodiff)) and $
   filesset ne 1 then begin    ; if processing entire day
  get_flnms2,YESTERDAY(date),spcrft,instrument,oldfiles
  szold=SIZE(oldfiles)
  if szold[0] lt 2 then begin
    print,"Problem with yesterday's get_flnms2 result."
    return
  endif
  previous=oldfiles[*,szold[2]-windw/2:szold[2]-1]
  oldfiles=0
  get_flnms2,TOMORROW(date),spcrft,instrument,nextfiles
  sznext=SIZE(nextfiles)
  if sznext[0] lt 2 then begin
    print,"Problem with tomorrow's get_flnms2 result."
    return
  endif 
  subseq=nextfiles[*,0:windw/2-1]
  nextfiles=0
endif
szfiles=SIZE(files)
nframes=szfiles[2]




; more error checking
szfiles=SIZE(files)
nframes=szfiles[2]
if nframes ge 1000 then begin
  print,'Unable to process.  Too many images.'
  return
endif
if not(keyword_set(entire)) then begin
  if (nframes lt windw+1) then begin
    print,'Variable windw is too large for number of available frames.'
    return
  endif
endif



; more initializations
diffs=FLTARR(outsz,outsz,nframes)
dttms=STRARR(nframes)
ctr=FLTARR(2)


;stop
; read in original images
for i=0,nframes-1 do begin
  print,'processing image number: ',i
  secchi_prep,files[*,i],hdrs,temp,calfac_off=calib,calimg_off= $
            calib,/silent,/smask_on,outsize=outsz,/polariz_on, $
            rotate_on=rot_on,rotinterp_on=rot_on,bkgimg_off= $
            back_off,/calroll,/interp  
  diffs[*,*,i]=temp
  if i eq 0 then begin
    wcshead=FITSHEAD2WCS(hdrs)
    ctr=WCS_GET_PIXEL(wcshead,[0,0])
  endif
  dttms[i]=hdrs[0].date_obs
endfor


;stop
; smooth pixels in time if desired
if KEYWORD_SET(smoother) then begin
  for i=0,outsz-1 do begin
    for j=0,outsz-1 do begin
;      diffs[i,j,*]=SMOOTH(diffs[i,j,*],smoother)
      diffs[i,j,*]=MEDIAN(REFORM(diffs[i,j,*]),smoother)
    endfor
  endfor
endif




if NOT(KEYWORD_SET(nodiff)) then begin
;  read images to finish off ends
  if NOT(keyword_set(entire)) and filesset ne 1then begin
    previms=FLTARR(outsz,outsz,windw/2)
    nxtimgs=FLTARR(outsz,outsz,windw/2)
    for i=0,windw/2-1 do begin
      secchi_prep,previous[*,i],hdrs,temp,calfac_off=calib,calimg_off=calib, $
          /silent,outsize=outsz,/polariz_on,/smask_on,rotate_on= $
          rot_on,rotinterp_on=rot_on,bkgimg_off=back_off, $
          /calroll,/interp
      previms[*,*,i]=temp
      secchi_prep,subseq[*,i],hdrs,temp,calfac_off=calib,calimg_off=calib, $
          /silent,outsize=outsz,/polariz_on,/smask_on,rotate_on= $
          rot_on,rotinterp_on=rot_on,bkgimg_off=back_off, $
          /calroll,/interp
      nxtimgs[*,*,i]=temp
    endfor
  endif




;  make background images and subtract from diffs
  if KEYWORD_SET(entire) then begin   ; makes strange base difference
    moving=MIN(diffs,dimension=3)
    for i=0,nframes-1 do diffs[*,*,i]=diffs[*,*,i]-moving
  endif else if filesset ne 1 then begin    ; ideal case
    moving=FLTARR(outsz,outsz)
    temp=FLTARR(outsz,outsz,windw+1)
    for i=0,nframes-1 do begin
      if i eq 0 then begin    ; must use previous images in window
        temp[*,*,0:windw/2-1]=previms
        temp[*,*,windw/2:*]=diffs[*,*,0:windw/2]
      endif else if i ge nframes-windw/2 then begin    ; must use subsequent images in window
        temp[*,*,0:windw-1]=temp[*,*,1:*]
        temp[*,*,windw]=nxtimgs[*,*,i-(nframes-windw/2)]
      endif else begin
        temp[*,*,0:windw-1]=temp[*,*,1:*]
        temp[*,*,windw]=diffs[*,*,i+windw/2]
      endelse 
      moving=MIN(temp,dimension=3)
      diffs[*,*,i]=diffs[*,*,i]-moving
    endfor
  endif else begin                               ; ends untrustworthy
    moving=FLTARR(outsz,outsz)
    for i=0,nframes-1 do begin
      if i lt windw/2 then begin
        temp=diffs[*,*,0:windw]
      endif else if i lt nframes-windw/2 then begin
        temp[*,*,0:windw-1]=temp[*,*,1:*]
        temp[*,*,windw]=diffs[*,*,i+windw/2]
      endif
      moving=MIN(temp,dimension=3)
      diffs[*,*,i]=diffs[*,*,i]-moving
    endfor
  endelse

endif          ; if making a difference movie



; sigma filter
if keyword_set(sigfilt) then begin
  szdiffs=SIZE(diffs)
  for i=0,szdiffs[3]-1 do diffs[*,*,i]= $
         sigma_filter(diffs[*,*,i],box, n_sigma=nsig,/iterate)
endif



;stop
; convert images to mass
if KEYWORD_SET(mass) then begin
  secchi_prep,files[*,0],hdr,img,calfac_off=calib,calimg_off= $
            calib,/silent,/smask_on,outsize=outsz,/polariz_on, $
            rotate_on=rot_on,rotinterp_on=rot_on,bkgimg_off= $
            back_off,/calroll,/interp  
  for i=0,szdiffs[3]-1 do diffs[*,*,i]= $
          SCC_CALC_CME_MASS(diffs[*,*,i],hdr,/all)
endif



if KEYWORD_SET(rescale) then begin
  ; normalize diffs and scale for screen
  locs=where(diffs ge cor1max)
  if locs[0] ne -1 then diffs[locs]=cor1max
  locs=where(diffs le cor1min)
  if locs[0] ne -1 then diffs[locs]=cor1min
  diffs=255*(TEMPORARY(diffs)-cor1min)/(cor1max-cor1min)
endif



;stop
; save difference movie images
if NOT(keyword_set(rotate)) then begin
  rot=' '
  rotate=0 
endif else rot='rot'
if NOT(keyword_set(sigfilt)) then begin
  sig=' '
  sigfilt=0
endif else sig='sig'
if KEYWORD_SET(calibrate) then cal='cal' else cal=' '
if KEYWORD_SET(smoother) then sm='sm'+STRING(smoother) $
           else sm=' '
if KEYWORD_SET(mass) then mass='mass' else mass=' '
filelist=files
if keyword_set(entire) then windw=nframes
instrument=STRLOWCASE(instrument)
if KEYWORD_SET(nodiff) then dir='/home/shaela/stereo/nodiffs/'+ $
   instrument+'/' else dir='/home/shaela/stereo/diffs/moving/'+ $
   instrument+'/'
if N_ELEMENTS(alpha) eq 0 then alpha=' '
if KEYWORD_SET(nodiff) then phrase='nodiff' else phrase='moving'
if KEYWORD_SET(no_bkg) then bkg='nobkg' else bkg=' '
file=dir+date+spcrft+STRING(outsz)+alpha+phrase+ $
    STRING(windw)+instrument+rot+sig+cal+bkg+sm+mass+'.sav'
file=STRCOMPRESS(file,/remove_all)
save,filename=file,diffs,spcrft,dttms,ctr,filelist, $
        rotate,instrument,calib,wcshead,smoother


;stop
end






; obsolete


;;           bckgrnd - subrtracts monthly background image    ;;
;;                     before finding window minimums         ;;

;; optionally subtract monthly background image
;if keyword_set(bckgrnd) then begin
;  ;fils=scc_getbkgimg(h,/total)
;  str='/service/stereo2/cor1/background/b/monthly_min/'+ $
;        '200703/mc1B*070321*.fts'
;  fils=file_search(str)
;  backfiles=sccreadfits(fils,hds)
;  backfiles=(2./3.)*TOTAL(backfiles,3)
;  sz=SIZE(diffs)
;  for i=0,sz[3]-1 do diffs[*,*,i]=diffs[*,*,i]-backfiles
;  if NOT(keyword_set(entire)) then begin
;    for i=0,windw/2-1 do previms[*,*,i]=previms[*,*,i]-backfiles
;    for i=0,windw/2-1 do nxtimgs[*,*,i]=nxtimgs[*,*,i]-backfiles
;  endif
;endif

;if keyword_set(bckgrnd) then begin
;  len=STRLEN(file)
;  file=STRMID(file,0,len-4)+'bkg.sav'
;endif



;  temp=(2.0/3.0)*(TOTAL(imgs,3))
;  cor1_quickpol,imgs,temp
;  if keyword_set(fivetwelve) then begin
;    diffs[*,*,i]=REBIN(temp,512,512)
;  endif else begin
;    sz=SIZE(temp)
;    if sz[1] gt 1024 then begin
;      print,'Image resolution is too high.  Rebinning to 1024x1024.'
;      diffs[*,*,i]=REBIN(temp,1024,1024)
;    endif else if sz[1] lt 1024 then begin
;      print,'Original images have a resolution less than 1024.'
;      stop
;    endif else diffs[*,*,i]=temp
;  endelse



;  ctr[0]=ctr[0]+TOTAL(hdrs.crpix1)
;  ctr[1]=ctr[1]+TOTAL(hdrs.crpix2)
;ctr=ctr/(3.*nframes)
;if keyword_set(fivetwelve) then ctr=ctr/2
;    cor1_quickpol,imgs,temp
;    if keyword_set(fivetwelve) then nxtimgs[*,*,i]=REBIN(temp,512,512) $
;          else nxtimgs[*,*,i]=temp



;    cor1_quickpol,imgs,temp
;    if keyword_set(fivetwelve) then previms[*,*,i]=REBIN(temp,512,512) $
;          else previms[*,*,i]=temp



;sz=SIZE(diffs)
;for i=0,sz[3]-1 do begin
;  for j=0,sz[1]-1 do begin
;    for k=0,sz[2]-1 do begin
;      temp2=diffs[j,k,i]
;      if temp2 ge 60.0 then temp2=60.0
;      if temp2 le -10.0 then temp2=-10.0
;      diffs[j,k,i]=temp2
;    endfor
;  endfor
;endfor
