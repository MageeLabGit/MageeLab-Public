#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#include <Waves Average>



// make chains w/ elevated density but uniform p01

Macro makemchain3()
	variable ss, s,j,c, jj, m, n, o, p,l,a,ll, state, nextpeekloc, cell, first,last,thresh,starts,index,all,fix,peek,offset,check, rise, tau_gauss,V_FitError,V_FitOptions,g=0,f=0,b=0,d=0, taua,spontsub,chi1,chi2,counter,refnumber,skipper, past, interval, 
	string	cha, cst, jst, cstavg, smthwave2, smthwave3, fase, rate2, Vm2, startname, start, wavname, Twave, Gwave, Rwave, rate, Vm, diffwave, kwave, Vwave, Iwave, Pwave, Vname, Iname, Pname, curname, poslist, Ext2, ExtR2, Ext, ExtR, current, current2
	Silent 1;
	
c=00
ss=20
trun0=0
trun0[14]=1
s=16
duplicate /o trun0 trunat
do
	m=0
	j=0
	cha="mch"+num2istr(c)

	make/o /n=5100 $cha
	•SetScale/P x 0,0.1,"", $cha


	a=abs(round(enoise (999)))+1
	if (a>=996)
		state=1
		print c
	else
		state=0
	endif

		
		if(c<=1599)
			ss=16
		endif
		if(c>=1600 && c<=1619)
			ss=20
		endif
		if(c>=1620 && c<=1724)
			ss=21
		endif
		if(c>=1725 && c<=1850)
			ss=8
		endif


		if (c==1600)
				ss=20
		endif
		if (c==1620)
				ss=21
			
		endif
		if (c==1725)
				ss=8
		endif
//	print c, ss
		if(c==00+s)
			rotate 1, run0
			SetScale/P x 0,0.1,"", run0
			rotate 1, trun0
			SetScale/P x 0,0.1,"", trun0

			s+=ss
		print c,"jjjj"
			if(c==1620 || c==1725)
				run0=1
				run0[14,23]=7
				print "tttttt"
				trun0=0
				trun0[14]=1

			endif


			duplicate /o run0, runmod0
			
			n=0	
			do
				concatenate /np{run0}, runmod0
				n+=1
			while (n<50)
	
		endif
		
		
		trunat+=trun0
		


do				
		
		if(c<=1850)
//				if (C<=1325 && C>=1275) //use if nonuniform 
//				ll=(runmod0[j]*4)-3
			ll=runmod0[j]
			l=1
		else
			ll=1
			l=1
		endif

	
	if (state==0 )
		m=abs(round(enoise (1000)))
		
		jj=(1000-996)*ll
		
		if(m>=1000-jj)
			state=1
			$cha[j]=1
		else
			state=0
			$cha[j]=0
		endif
	else 
		m=abs(round(enoise (1000)))
		

		jj=(50*l)+900
		
		if(m>=1000-jj)
		
			state=1
			$cha[j]=1
		else
			state=0
			$cha[j]=0
		endif
	endif
	
	j+=1
	
	while (j<5100)
	
	truna+=trun0	
	
	c+=1
	
	while (c<2000)

end


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
p=00	
	do
	
		jst="cavg"+num2istr(cell)+"_"+num2istr(p)
		$jst-=inhib2
	
		findlevels /r=[100,5100]/q /edge=1 /d=desty $jst, 11.5
		wavestats /q $jst
			
			if (V_max>=11)
//				insertpoints 0,1,pks
//	   		pks[0]= runs1(V_maxloc)
//				insertpoints p,1,pkt
//				pkt[0]= V_maxloc
			endif

	
	
	
		jjj=0	
		if (V_flag==1)
			insertpoints jjj,1,ts1
			ts1[jjj]= timess1(desty[0])
			insertpoints jjj,1,ps1
			ps1[jjj]= runs2(desty[0])
			
			insertpoints jjj,1,pks
			pks[jjj]= V_max
			insertpoints jjj,1,pkt
			pkt[jjj]= p
	


		do
			insertpoints jjj,1,ts
			ts[jjj]= timess1(desty[jjj])
			insertpoints jjj,1,ps
			ps[jjj]= runs2(desty[jjj])
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
wavestats /q ps1
print "1", V_npnts

end