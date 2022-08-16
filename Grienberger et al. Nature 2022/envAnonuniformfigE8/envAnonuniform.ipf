#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#include <Waves Average>

/// all other analysis code is the same as found in envAuniform
/// make chains 000-999

Macro makemchain1()
	variable ss, s,j,c, jj, m, n, o, p,l,a,ll, state, nextpeekloc, cell, first,last,thresh,starts,index,all,fix,peek,offset,check, rise, tau_gauss,V_FitError,V_FitOptions,g=0,f=0,b=0,d=0, taua,spontsub,chi1,chi2,counter,refnumber,skipper, past, interval, 
	string	cha, cst, jst, cstavg, smthwave2, smthwave3, fase, rate2, Vm2, startname, start, wavname, Twave, Gwave, Rwave, rate, Vm, diffwave, kwave, Vwave, Iwave, Pwave, Vname, Iname, Pname, curname, poslist, Ext2, ExtR2, Ext, ExtR, current, current2
	Silent 1;
	
c=000
s=7
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


do
		
		
		if(c==000+s)
			rotate 1, run0
			SetScale/P x 0,0.1,"", run0

			s+=7

			duplicate /o run0, runmod0
			
			n=0	
			do
				concatenate /np{run0}, runmod0
				n+=1
			while (n<50)
	
		endif

			
		
		if(c<=699)
			if (C<=325 && C>=275)
				ll=(runmod0[j]*5)-4

				l=1
			else
				ll=runmod0[j]

				l=1
			endif
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
	
	c+=1
	
	while (c<1000)

end




/// make chains 1000-1999

Macro makemchain2()
	variable ss, s,j,c, jj, m, n, o, p,l,a,ll, state, nextpeekloc, cell, first,last,thresh,starts,index,all,fix,peek,offset,check, rise, tau_gauss,V_FitError,V_FitOptions,g=0,f=0,b=0,d=0, taua,spontsub,chi1,chi2,counter,refnumber,skipper, past, interval, 
	string	cha, cst, jst, cstavg, smthwave2, smthwave3, fase, rate2, Vm2, startname, start, wavname, Twave, Gwave, Rwave, rate, Vm, diffwave, kwave, Vwave, Iwave, Pwave, Vname, Iname, Pname, curname, poslist, Ext2, ExtR2, Ext, ExtR, current, current2
	Silent 1;
	
c=1000
s=7
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


do
		
		
		if(c==1000+s)
			rotate 1, run0
			SetScale/P x 0,0.1,"", run0

			s+=7

			duplicate /o run0, runmod0
			
			n=0	
			do
				concatenate /np{run0}, runmod0
				n+=1
			while (n<50)
	
		endif

			
		
		if(c<=1699)
			if (C<=1325 && C>=1275)
				ll=(runmod0[j]*5)-4

				l=1
			else
				ll=runmod0[j]

				l=1
			endif
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
	
	c+=1
	
	while (c<2000)

end

