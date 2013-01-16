BEGIN {
  total[0-1000]=0.0;
  T_e[0-1000]=0.0;
  dr[0-1000]=0.0;
  i=0
} 
	   	   

{    if ($1 == "Tphot") {temp=$6};
 if ($1 != "Tphot" && $1 != "Dphot" && $1 != "#" && $1 != "Step") {temp2=$2;temp3=$5};
     if ($1 == "Step") {total[i]=temp;
     T_e[i]=temp2;
     dr[i]=temp3
			     i+=1
                       }
}

END{    for (j=1; j<i; j+=1){
  printf ("Temp: %12.6e, dr: %12.6e, Flux: %12.6e erg s-1 cm-2\n",T_e[j],dr[j],total[j])}
    exit
}
