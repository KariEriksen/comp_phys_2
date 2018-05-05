N_mc="17";
for N_d in "3"
	do
		for N_p in "10"
			do
				for dt in "0.0001" "0.001" "0.01" "0.1" "1.0"
					do
						echo "simulating with Nd: $N_d and N_p: $N_p and Nmc: $N_mc and dt: $dt" 
						./app_c.x $N_p $N_d $N_mc $dt 
		done
	done
done
