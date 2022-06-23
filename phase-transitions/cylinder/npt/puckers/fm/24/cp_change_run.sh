file0=in.cyl_npt.py	#file name 
kappa=2.0		#bending 		
kstretch=100  		#harmonic spring coefficient  
strain=0	     	#pre strain before relaxation  
wrap=1			#planar=0, cylinder=1
noise=1			#add noise to position
#kin=0.100		#kin energy/temperature
spin=0   		#puckers=0, stitches=1, initialize spins height in the right locations
order=1 		#AFM=-1, FM=1
height=-0.45		#initialize "spin" height
frac=0.1		#dilated bonds 
L=24
for kin in 0.100 0.250 
do 
mkdir $kin
cd $kin 
for run in 1 2  
do 
mkdir r$run
cd r$run
cp ../../$file0 .
sed -i 's/kappa=2.0/kappa='$kappa'/g' $file0
sed -i 's/kstretch=100/kstretch='$kstretch'/g' $file0 
sed -i 's/nx=48/nx='$L'/g' $file0
sed -i 's/ny=48/ny='$L'/g' $file0
sed -i 's/nz=48/nz='$L'/g' $file0
sed -i 's/strain=0/strain='$strain'/g' $file0  
sed -i 's/wrap=0/wrap='$wrap'/g' $file0
sed -i 's/noise=0/noise='$noise'/g' $file0
sed -i 's/kin=0.100/kin='$kin'/g' $file0
sed -i 's/spin=0/spin='$spin'/g' $file0  
sed -i 's/order=-1/order='$order'/g' $file0  
sed -i 's/height=0.45/height='$height'/g' $file0
sed -i 's/frac=0.1/frac='$frac'/g' $file0
sed -i 's/Run=1/Run='$run'/g' $file0
#use your script to submit the job to your cluster
#sbatch ../../singularity_cpu.sh
cd ..
done
cd ..
done

#NOTE 
#the DEFAULT input uses parameter for AFM puckered sheets of size 48x48 
