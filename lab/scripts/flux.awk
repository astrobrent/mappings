BEGIN {
    Energy[1-1802]=0.0;Flux[1-1802]=0;iFlux[1-1802]=0;i=0;files=0;
    cls=2.99792458E10;
    plk = 6.6260755E-27;
    eV = 1.60217733E-12;
    pi=3.1415926536;
    evplk=eV/plk;
    Ryd=13.59844

}
	   	   

{	if ($2 == "PHOTON") {files+=1};
        if ($2 == "RUN:") {Runtitle=$3};
        if ($1 != "%" && NF == 2  ) {i+=1;
	   Energy[i]=$1*evplk;Flux[i]=$2;
	   if ($1 >=Ryd) {iFlux[i]=$2}
                       }
}

END{Total=0.0;iTotal=0.0;
    for (j=1; j<i; j+=1){
        Total +=Flux[j]*(Energy[j+1]-Energy[j]);
	iTotal+=iFlux[j]*(Energy[j+1]-Energy[j])
	    };
    printf ("Total Flux\t%12.6e erg s-1 cm-2\n",Total*4.0*pi);
    printf ("Total Ionizing Flux\t%12.6e erg s-1 cm-2\n",iTotal*4.0*pi);
    exit
}
