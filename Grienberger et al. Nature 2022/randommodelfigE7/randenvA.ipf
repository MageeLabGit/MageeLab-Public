#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#include <Waves Average>



//make chains


Macro makemchain()
	variable j,c, jj, m, n, o, p,l,a, state, nextpeekloc, cell, first,last,thresh,starts,index,all,fix,peek,offset,check, rise, tau_gauss,V_FitError,V_FitOptions,g=0,f=0,b=0,d=0, taua,spontsub,chi1,chi2,counter,refnumber,skipper, past, interval, 
	string	cha, cst, jst, cstavg, smthwave2, smthwave3, fase, rate2, Vm2, startname, start, wavname, Twave, Gwave, Rwave, rate, Vm, diffwave, kwave, Vwave, Iwave, Pwave, Vname, Iname, Pname, curname, poslist, Ext2, ExtR2, Ext, ExtR, current, current2
	Silent 1;
	
c=0
do
m=0
j=0
cha="mch"+num2istr(c)

make/o /n=5000 $cha
•SetScale/P x 0,0.1,"", $cha

a=abs(round(enoise (999)))+1
if (a>=996)
	state=1
	print c
else
	state=0
endif


do
	
		if (state==0 )
			m=abs(round(enoise (1000)))
		
			jj=(1000-996)
		
		if(m>=1000-jj)
			state=1
			$cha[j]=1
		else
			state=0
			$cha[j]=0
		endif
	else 
		m=abs(round(enoise (1000)))
			jj=(1000-50)

		if(m>=1000-jj)
			state=1
			$cha[j]=1
		else
			state=0
			$cha[j]=0
		endif
	endif
	
	j+=1
	
	while (j<5000)
	
	c+=1
	
	while (c<2000)

end

//calculate dwell times
Macro onoffchains()
	variable j,c, k, m, n, total, lstpnt, nextpeekloc, cell, first,last,thresh,starts,index,all,fix,peek,offset,check, rise, tau_gauss,V_FitError,V_FitOptions,g=0,f=0,b=0,l=0,d=0, taua,spontsub,chi1,chi2,counter,refnumber,skipper, past, interval, 
	string	cha, fase2, smthwave, smthwave4, smthwave2, smthwave3, fase, rate2, Vm2, startname, start, wavname, Twave, Gwave, Rwave, rate, Vm, diffwave, kwave, Vwave, Iwave, Pwave, Vname, Iname, Pname, curname, poslist, Ext2, ExtR2, Ext, ExtR, current, current2
	Silent 1;

	make /o/n=0 opentime
	make /o/n=0 closetime

//$wavname+=2


c=0
m=0
n=0
do
	cha="mch"+num2istr(c)
	if ($cha[0]==0)
	findlevels /q/dest=trans $cha 0.5
	j=0
	
	do	
//		print "j", m,j, (trans[j+1]-trans[j])
		insertpoints m,1, opentime
		
		opentime[m]=(trans[j+1]-trans[j])
		m+=1
		j+=2
		
	while(j<=V_LevelsFound-2)
	
	k=1
	
	do
		insertpoints n,1, closetime
//		print "k", n,k, (trans[k+1]-trans[k])
		closetime[n]=(trans[k+1]-trans[k])
		
		n+=1
		k+=2
	while(k<=V_LevelsFound-2)
	else
	print c
	endif
	
	c+=1
	
	while (c<2000)
	
	
		Make/N=50/O opentime_Hist
		Histogram/B={0,1,50} opentime,opentime_Hist;
		Display opentime_Hist
		
		Make/N=100/O closetime_Hist
		Histogram/B={0,1,200} closetime,closetime_Hist;
		Display closetime_Hist


	
end



//cut chains into laps 

Macro chains2laps()
	variable j, jj, m, n, o, p,l,jjj, lstpnt, endy, cell, first,last,thresh,starts,index,all,fix,peek,offset,check, rise, tau_gauss,V_FitError,V_FitOptions,g=0,f=0,b=0,d=0, taua,spontsub,chi1,chi2,counter,refnumber,skipper, past, interval, 
	string	cst, jst, ast, smthwave4, smthwave2, smthwave3, fase, rate2, Vm2, startname,  wavname, Twave, Gwave, Rwave, rate, Vm, diffwave, kwave, Vwave, Iwave, Pwave, Vname, Iname, Pname, curname, poslist, Ext2, ExtR2, Ext, ExtR, current, current2
	Silent 1;
	
make /o/n=0 activelist	

cell=00
do

ast="lapsa"+num2istr(cell)

jj=0	
j=0
jst="mch"+num2istr(cell)
SetScale/P x 0,0.1,"", $jst
	do
		
		if (j==00)
			duplicate /o /r=[0+j,99+j] $jst $ast
			SetScale/P x 0,0.1,"", $ast
			
			findlevel /q $ast, 0.99
			make /o/n=0 active
			if (v_flag==0)
				insertpoints jj,1,active
				active[jj]=1
			else
				insertpoints jj,1,active
				active[jj]=0
			endif

			
		else
			duplicate /o /r=[0+j,99+j] $jst temp1
			SetScale/P x 0,0.1,"", temp1
			$ast+=temp1
			smooth 3, temp1
			findlevel /q temp1, 0.99
			
			if (v_flag==0)
				insertpoints jj,1,active
				active[jj]=1
			else
				insertpoints jj,1,active
				active[jj]=0
			endif

		endif
		
		j+=100
		
		while (j<5000)
	insertpoints cell, 1, activelist	
	wavestats /q active
	activelist[cell]=V_avg
		
	smooth 3, $ast	
	cell+=1
while (cell<2000)
	
end


//calculate properties of chains
	
	Macro chaincorr()
	variable j, jj, m, n, o, p,l,a, state, nextpeekloc, cells, first,last,thresh,starts,index,all,fix,peek,offset,check, rise, tau_gauss,V_FitError,V_FitOptions,g=0,f=0,b=0,d=0, taua,spontsub,chi1,chi2,counter,refnumber,skipper, past, interval, endy
	string	cst, jst, cstavg, smthwave2, smthwave3, fase, rate2, Vm2, startname, start, wavname, Twave, Gwave, Rwave, rate, Vm, diffwave, kwave, Vwave, Iwave, Pwave, Vname, Iname, Pname, curname, poslist, Ext2, ExtR2, Ext, ExtR, current, current2
	Silent 1;
	
	
	make /o /n=0 runcor
	make /o /n=0 runprob
	make /o /n=0 runsig
	make /o /n=0 runsel
	make /o /n=0 cellid
	make /o /n=0 cellpl
	make /o/n=0 celltyp

	j=0
	m=0
	jj=0
	do
			insertpoints m, 1, celltyp
			celltyp[m]=endy

			jst="laps"+num2istr(j)
			
//			•Smooth 3,$jst
			wavestats /q $jst
			insertpoints m, 1, runsel
			runsel[m]=V_max/v_avg
			insertpoints m, 1, cellpl
			cellpl[m]=V_maxloc

			insertpoints m, 1, cellid
			cellid[m]=m
	
			statslinearcorrelationTest /q $jst, run3_l
			insertpoints m, 1, runcor
			runcor[m]=W_StatsLinearCorrelationTest[1]
			insertpoints m, 1, runprob
			runprob[m]=W_StatsLinearCorrelationTest[7]
			insertpoints m, 1, runsig
			if (W_StatsLinearCorrelationTest[7]<=0.049)
				runsig[m]=1
			else
				runsig[m]=0
			endif
			
			
	if (j==0)
		duplicate /o $jst lapavg
	else
		lapavg+=$jst
	endif


			
	m+=1
	j+=1
	while (m<2000)


// sort cellpl,	runcor, runprob,runsig,runsel,cellid, cellpl, celltyp

display lapavg
end






//heatmap chains

Macro chainspace(last)
	variable cell, first, last, total, r, c, shift, maxa, f, j
	string	rate, allrate, rateb
	silent 1
	
	make /o/n=1 cshw
	make /o/n=100 temp
	make /o/n=(100,last) rateplotfa2
	make /o/n=100 ratepavg
	make /o/n=1000 maxloc

	c=0
	r=0

	do
		j=cellid2[c]

		rate = "lapsa" + num2istr(j)
		
		r=0
		temp=0
		do
			wavestats /q /r=(ttime1[r],ttime2[r]) $rate
			rateplotfa2[r][c]=V_avg
			temp[r]+=V_avg
			
			r+=1
		while (r<100)
		
		
		wavestats /q temp 
		maxloc[c]=V_maxloc
		maxa=V_max
		
		r=0
		do
		
			rateplotfa2[r][c]/=maxa
			
			r+=1
		while (r<100)


			
		c+=1
	while(c<last)
	

	display; appendimage rateplotfa2
	ModifyImage rateplotfa2 ctab= {*,*,Rainbow,1}
 	SetAxis left 999.5,-0.5

	end
	
	
	
//heatmap single chains per lap

Macro chains2indilaps(cell)
	variable r, j, jj, m, n, o, p,l,jjj, lstpnt, endy, cell, first,last,thresh,starts,index,all,fix,peek,offset,check, rise, tau_gauss,V_FitError,V_FitOptions,g=0,f=0,b=0,d=0, taua,spontsub,chi1,chi2,counter,refnumber,skipper, past, interval, 
	string	rp, cst, jst, ast, smthwave4, smthwave2, smthwave3, fase, rate2, Vm2, startname,  wavname, Twave, Gwave, Rwave, rate, Vm, diffwave, kwave, Vwave, Iwave, Pwave, Vname, Iname, Pname, curname, poslist, Ext2, ExtR2, Ext, ExtR, current, current2
	Silent 1;
	
j=0
jj=0
jst="mch"+num2istr(cell)
SetScale/P x 0,0.1,"", $jst

	do
			ast="vm"+num2istr(cell)+"_"+num2istr(jj)
			duplicate /o /r=[0+j,99+j] $jst $ast
			SetScale/P x 0,0.1,"", $ast
		
		j+=100
		jj+=1
	while (j<5000)
		
//	smooth 3, $ast	
  rp="rateplot"+num2istr(cell)

	make /o/n=(100,50) $rp

	jj=0
	do
		ast="vm"+num2istr(cell)+"_"+num2istr(jj)		
		r=0

		do
			wavestats /q /r=(ttime1[r],ttime2[r]) $ast
			$rp[r][jj]=V_avg

			r+=1
		while (r<100)
			
		jj+=1
	while(jj<50)
	

	display; appendimage $rp
	ModifyImage $rp ctab= {*,*,Rainbow,1}
 	SetAxis left 49.5,-0.5

	
end	


///---------------------------------------
//postsynaptic summation of chains


//make posysynaptic summations

Macro combochains(psc,cell)
	variable j, jj, m, n, o, p,l, jjj, psc, laps, endy, cell, first,last,thresh,starts,index,all,fix,peek,offset,check, rise, tau_gauss,V_FitError,V_FitOptions,g=0,f=0,b=0,d=0, taua,spontsub,chi1,chi2,counter,refnumber,skipper, past, interval, 
	string	cha, cst, jst, smthwave4, smthwave2, smthwave3, fase, rate2, Vm2, startname,  wavname, Twave, Gwave, Rwave, rate, Vm, diffwave, kwave, Vwave, Iwave, Pwave, Vname, Iname, Pname, curname, poslist, Ext2, ExtR2, Ext, ExtR, current, current2
	Silent 1;


n=1
m=0
j=0
p=0	

	do
	
		make /o/n=(cell) rancell
		l=0
		do
			rancell[l]=abs(round(enoise (1999)))
		
			l+=1
		while (l<cell)
		
			jj=0
			do
			
				jst="cavg"+num2istr(cell)+"_"+num2istr(p)
				
				j=(rancell[jj])
				
				cha="mch"+num2istr(j)
			
				if (jj==0)				

					duplicate /o $cha $jst
					SetScale/P x 0,0.1,"", $jst
				
				else
					$jst+=$cha
					SetScale/P x 0,0.1,"", $jst
				endif
			
			jj+=1
			m+=1
		while (jj<cell)
		
		Smooth 10,$jst
		
		
		wavestats /q $jst
		n+=1
		p+=1
	while (n<=psc)


	
end

//make inhibtion trace 

Macro makeinhib()
	variable j, jj, m, n, o, p,l,jjj, lstpnt, endy, cell, first,last,thresh,starts,index,all,fix,peek,offset,check, rise, tau_gauss,V_FitError,V_FitOptions,g=0,f=0,b=0,d=0, taua,spontsub,chi1,chi2,counter,refnumber,skipper, past, interval, 
	string	ast, jst, smthwave4, smthwave2, smthwave3, fase, rate2, Vm2, startname,  wavname, Twave, Gwave, Rwave, rate, Vm, diffwave, kwave, Vwave, Iwave, Pwave, Vname, Iname, Pname, curname, poslist, Ext2, ExtR2, Ext, ExtR, current, current2
	Silent 1;
	
	


	
j=0
n=1
jj=0

	do
		ast="mch"+num2istr(j)
		if (j==0)
			
			duplicate /o $ast, inhib2
		else
			inhib2+=$ast
		endif


	j+=1
	while (j<2000)
	
	Smooth 10,inhib2
	
	inhib2/=(2000/100)  // normalize to summations
	display inhib2
end



//Find threshold crossings

Macro find2(pcs,cell)
	variable j, jj, m, n, o, p,l, jjj, pcs, endy, cell, first,last,thresh,starts,index,all,fix,peek,offset,check, rise, tau_gauss,V_FitError,V_FitOptions,g=0,f=0,b=0,d=0, taua,spontsub,chi1,chi2,counter,refnumber,skipper, past, laps 
	string	cst, jst, smthwave4, smthwave2, smthwave3, fase, rate2, Vm2, startname,  wavname, Twave, Gwave, Rwave, rate, Vm, diffwave, kwave, Vwave, Iwave, Pwave, Vname, Iname, Pname, curname, poslist, Ext2, ExtR2, Ext, ExtR, current, current2
	Silent 1;


make /o /n=0 ps1
make /o /n=0 ts1
make /o /n=0 ps
make /o /n=0 ts
make /o /n=0 pks
make /o /n=0 dt
make /o /n=0 pkt

n=1
m=0
j=0
p=0	
	do
	
		jst="cavg"+num2istr(cell)+"_"+num2istr(p)
		$jst-=inhib2
	
		findlevels /q /edge=1 /d=desty $jst, 10.25
		
		
	
			wavestats /q $jst
			
			
			if (V_max>=13)
				insertpoints 0,1,pks
				pks[0]= runs(V_maxloc)
				insertpoints 0,1,pkt
				pkt[0]= V_maxloc
				
				print p,pks[0], v_maxloc
			endif

	
	
	
		jjj=0	
		if (V_flag==1)
			insertpoints jjj,1,ts1
			ts1[jjj]= timess(desty[0])
			insertpoints jjj,1,ps1
			ps1[jjj]= runs(desty[0])


		do
			insertpoints jjj,1,ts
			ts[jjj]= timess(desty[jjj])
			insertpoints jjj,1,ps
			ps[jjj]= runs(desty[jjj])
			insertpoints jjj,1,dt
			dt[jjj]= desty[jjj]
		
			jjj+=1
		
		while (jjj<V_LevelsFound)
		endif
		
		$jst+=inhib2
		
		p+=1
	while (p<pcs)

make /o /n=30 ps_Hist
make /o /n=30 ts_Hist
make /o /n=30 pks_Hist
•Histogram/B={0,6,30} ps,ps_Hist;DelayUpdate
•Histogram/B={0,6,30} ts,ts_Hist;DelayUpdate
Histogram/B={0,6,30} pks,pks_Hist;
make /o /n=30 ps1_Hist
make /o /n=30 ts1_Hist
•Histogram/B={0,6,30} ps1,ps1_Hist;DelayUpdate
•Histogram/B={0,6,30} ts1,ts1_Hist;DelayUpdate
wavestats /q ps
print V_npnts

end






///---------------------------------------
//cell-cell correaltions


//make laps for 1-25 laps (chains500-1500)


Macro chains2laps2()
	variable j, jj, m, n, o, p,l,jjj, lstpnt, endy, cell, first,last,thresh,starts,index,all,fix,peek,offset,check, rise, tau_gauss,V_FitError,V_FitOptions,g=0,f=0,b=0,d=0, taua,spontsub,chi1,chi2,counter,refnumber,skipper, past, interval, 
	string	cst, jst, ast, smthwave4, smthwave2, smthwave3, fase, rate2, Vm2, startname,  wavname, Twave, Gwave, Rwave, rate, Vm, diffwave, kwave, Vwave, Iwave, Pwave, Vname, Iname, Pname, curname, poslist, Ext2, ExtR2, Ext, ExtR, current, current2
	Silent 1;
	
make /o/n=0 activelist	

cell=00
do

ast="lapsa"+num2istr(cell)

jj=0	
j=00
jst="mch"+num2istr(cell+500)
SetScale/P x 0,0.1,"", $jst
	do
		
		if (j==00)
			duplicate /o /r=[0+j,99+j] $jst $ast
			SetScale/P x 0,0.1,"", $ast
			

			
		else
			duplicate /o /r=[0+j,99+j] $jst temp1
			SetScale/P x 0,0.1,"", temp1
			$ast+=temp1
			

		endif
		
		j+=100
		
		while (j<2500)
	insertpoints cell, 1, activelist	
	wavestats /q active
	activelist[cell]=V_avg
		
	smooth 3, $ast	
	cell+=1
while (cell<1000)
	
end

//make laps for 26-50 laps (chains500-1500)

Macro chains2laps3()
	variable j, jj, m, n, o, p,l,jjj, lstpnt, endy, cell, first,last,thresh,starts,index,all,fix,peek,offset,check, rise, tau_gauss,V_FitError,V_FitOptions,g=0,f=0,b=0,d=0, taua,spontsub,chi1,chi2,counter,refnumber,skipper, past, interval, 
	string	cst, jst, ast, smthwave4, smthwave2, smthwave3, fase, rate2, Vm2, startname,  wavname, Twave, Gwave, Rwave, rate, Vm, diffwave, kwave, Vwave, Iwave, Pwave, Vname, Iname, Pname, curname, poslist, Ext2, ExtR2, Ext, ExtR, current, current2
	Silent 1;
	
make /o/n=0 activelist	

cell=00
do

ast="lapsb"+num2istr(cell)

jj=2500	
j=00
jst="mch"+num2istr(cell+500)
SetScale/P x 0,0.1,"", $jst
	do
		
		if (j==2500)
			duplicate /o /r=[0+j,99+j] $jst $ast
			SetScale/P x 0,0.1,"", $ast
			

			
		else
			duplicate /o /r=[0+j,99+j] $jst temp1
			SetScale/P x 0,0.1,"", temp1
			$ast+=temp1
			

		endif
		
		j+=100
		
		while (j<5000)
	insertpoints cell, 1, activelist	
	wavestats /q active
	activelist[cell]=V_avg
		
	smooth 3, $ast	
	cell+=1
while (cell<1000)
	
end

//calculate properties of chains for lap1-25

Macro chaincorr2()
	variable j, jj, m, n, o, p,l,a, state, nextpeekloc, cells, first,last,thresh,starts,index,all,fix,peek,offset,check, rise, tau_gauss,V_FitError,V_FitOptions,g=0,f=0,b=0,d=0, taua,spontsub,chi1,chi2,counter,refnumber,skipper, past, interval, endy
	string	cst, jst, cstavg, smthwave2, smthwave3, fase, rate2, Vm2, startname, start, wavname, Twave, Gwave, Rwave, rate, Vm, diffwave, kwave, Vwave, Iwave, Pwave, Vname, Iname, Pname, curname, poslist, Ext2, ExtR2, Ext, ExtR, current, current2
	Silent 1;
	
	
	make /o /n=0 runcor2
	make /o /n=0 runprob2
	make /o /n=0 runsig2
	make /o /n=0 runsel2
	make /o /n=0 cellid2
	make /o /n=0 cellpl2
	make /o/n=0 celltyp2

	j=00
	m=0
	jj=0
	do
			insertpoints m, 1, celltyp2
			celltyp2[m]=endy

			jst="lapsa"+num2istr(j)
			
//			•Smooth 3,$jst
			wavestats /q $jst
			insertpoints m, 1, runsel2
			runsel2[m]=V_max/v_avg
			insertpoints m, 1, cellpl2
			cellpl2[m]=V_maxloc

			insertpoints m, 1, cellid2
			cellid2[m]=j
	
			statslinearcorrelationTest /q $jst, run3_l
			insertpoints m, 1, runcor2
			runcor2[m]=W_StatsLinearCorrelationTest[1]
			insertpoints m, 1, runprob2
			runprob2[m]=W_StatsLinearCorrelationTest[7]
			insertpoints m, 1, runsig2
			if (W_StatsLinearCorrelationTest[7]<=0.049)
				runsig2[m]=1
			else
				runsig2[m]=0
			endif
			
			
	if (j==00)
		duplicate /o $jst lapavg2
	else
		lapavg2+=$jst
	endif


			
	m+=1
	j+=1
	while (m<1000)


// sort cellpl,	runcor, runprob,runsig,runsel,cellid, cellpl, celltyp

display lapavg3
end


//calculate properties of chains for lap26-50

	Macro chaincorr3()
	variable j, jj, m, n, o, p,l,a, state, nextpeekloc, cells, first,last,thresh,starts,index,all,fix,peek,offset,check, rise, tau_gauss,V_FitError,V_FitOptions,g=0,f=0,b=0,d=0, taua,spontsub,chi1,chi2,counter,refnumber,skipper, past, interval, endy
	string	cst, jst, cstavg, smthwave2, smthwave3, fase, rate2, Vm2, startname, start, wavname, Twave, Gwave, Rwave, rate, Vm, diffwave, kwave, Vwave, Iwave, Pwave, Vname, Iname, Pname, curname, poslist, Ext2, ExtR2, Ext, ExtR, current, current2
	Silent 1;
	
	
	make /o /n=0 runcor3
	make /o /n=0 runprob3
	make /o /n=0 runsig3
	make /o /n=0 runsel3
	make /o /n=0 cellid3
	make /o /n=0 cellpl3
	make /o/n=0 celltyp3

	j=00
	m=0
	jj=0
	do
			insertpoints m, 1, celltyp3
			celltyp3[m]=endy

			jst="lapsb"+num2istr(j)
			
//			•Smooth 3,$jst
			wavestats /q $jst
			insertpoints m, 1, runsel3
			runsel3[m]=V_max/v_avg
			insertpoints m, 1, cellpl3
			cellpl3[m]=V_maxloc

			insertpoints m, 1, cellid3
			cellid3[m]=j
	
			statslinearcorrelationTest /q $jst, run3_l
			insertpoints m, 1, runcor3
			runcor3[m]=W_StatsLinearCorrelationTest[1]
			insertpoints m, 1, runprob3
			runprob3[m]=W_StatsLinearCorrelationTest[7]
			insertpoints m, 1, runsig3
			if (W_StatsLinearCorrelationTest[7]<=0.049)
				runsig3[m]=1
			else
				runsig3[m]=0
			endif
			
			
	if (j==00)
		duplicate /o $jst lapavg3
	else
		lapavg3+=$jst
	endif


			
	m+=1
	j+=1
	while (m<1000)


// sort cellpl,	runcor, runprob,runsig,runsel,cellid, cellpl, celltyp

display lapavg3
end


//calculate cell-cell correlations
	
Macro meancor()
	variable j, jj, m, n, o, p,l,a, state, nextpeekloc, cells, first,last,thresh,starts,index,all,fix,peek,offset,check, rise, tau_gauss,V_FitError,V_FitOptions,g=0,f=0,b=0,d=0, taua,spontsub,chi1,chi2,counter,refnumber,skipper, past, interval, endy
	string	ratea, rateb, cst, jst, cstavg, smthwave2, smthwave3, fase, rate2, Vm2, startname, start, wavname, Twave, Gwave, Rwave, rate, Vm, diffwave, kwave, Vwave, Iwave, Pwave, Vname, Iname, Pname, curname, poslist, Ext2, ExtR2, Ext, ExtR, current, current2
	Silent 1;
	
	
	make /o /n=0 cellcor1

	j=0
	m=0
	do
			m=cellid2[j]

			ratea = "lapsa" + num2istr(m)
			rateb = "lapsb" + num2istr(m)
			
			statslinearcorrelationTest /q $ratea, $rateb
			
			insertpoints j, 1, cellcor1
			cellcor1[j]=W_StatsLinearCorrelationTest[1]

		
			j+=1
		while(j<1000)
		
		print median(cellcor1)
	end

//find high correlation chains

Macro unity()
	variable j, jj, m, n, o, p,l,a, s, state, nextpeekloc, cells, first,last,thresh,starts,index,all,fix,peek,offset,check, rise, tau_gauss,V_FitError,V_FitOptions,g=0,f=0,b=0,d=0, taua,spontsub,chi1,chi2,counter,refnumber,skipper, past, interval, endy
	string	ratea, rateb, cst, jst, cstavg, smthwave2, smthwave3, fase, rate2, Vm2, startname, start, wavname, Twave, Gwave, Rwave, rate, Vm, diffwave, kwave, Vwave, Iwave, Pwave, Vname, Iname, Pname, curname, poslist, Ext2, ExtR2, Ext, ExtR, current, current2
	Silent 1;
	
	
	make /o /n=0 union2
	make /o /n=0 union3
	
	j=0
	m=0
	n=0
	do
			m=cellpl2[j]
			n=cellpl3[j]	
			
			if ((n<=(m+0.3)) && (n>=(m-0.3)))
			
				insertpoints j,1, union2
				union2[j]=cellpl2[j]
				insertpoints j,1, union3
				union3[j]=cellpl3[j]
				
			endif
			j+=1

		while	(j<1000)

display union2 vs union3
wavestats /q union2 
print V_npnts/1000
		
end
