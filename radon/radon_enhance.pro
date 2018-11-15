function radon_enhance,httm,bottom_val=bottom_val,mints=mints


szhttm=SIZE(httm)
try=httm-mean(httm)
result=radon(try,theta=ts,rho=rs,/gray)
result2=MAKE_ARRAY(size=SIZE(result))
if N_ELEMENTS(bottom_val) eq 0 then begin
  histup,result
  stop
  read,'bottom_val= ',bottom_val
endif
list=WHERE(result gt bottom_val)
result2[list]=result[list]
stop
if N_ELEMENTS(mints) eq 0 then begin
  print,SIZE(ts)
  read,'minimum ts= ',mints
endif
ts2=ts[mints:*]
result3=radon(result2[mints:*,*],/backproject,nx= $
        szhttm[1],ny=szhttm[2],theta=ts2,rho=rs)
result4=MAKE_ARRAY(size=SIZE(result3))
list2=WHERE(result3 ne 0)
result4[list2]=MEAN(result3[list2])


return,result4
end
