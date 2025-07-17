#!/bin/bash
set -e -o pipefail

version () {
	echo "version: 1.0.0"
	exit
}

usage (){
	cat << EOF
Usage: `basename $0` -b stat_k-mer_size -e end_k-mer_size -i step_size -r reads
-b <INT>                                      K-mer length for beginning of loop forecasts (Bp)                                      [required]
-e <INT>                                      K-mer length for ending of loop forecasts (Bp)                                         [required]
-s <INT>                                      K-mer length iteration step size for loop forecasts (Bp)                               [required]
-r  <.cram|.[bs]am|.db|.dam|.f[ast][aq][.gz]> datasets for generating K-mers                                                         [required] 
-p  <INT>                                     ploidy (1, 2, 3, 4, 5, or 6) for GenomeScope2 to use                                   [default: 2]
-f                                            full mode to keep intermediate files                                                   [default: disabled]
-c                                            re-correct the HiFi consensus reads using hifiasm                                      [default: disabled]
-o  <output_prefix>                           prefix of the output files                                                             [default: out]
-t  <INT>                                     the number of threads to use                                                           [default: 1]
-M  <INT>                                     Use -M GB of memory in K-mer counting steps of FastK                                   [default: 12]
-W  <workpdir>                                directory to execute loop forecasts                                                    [default: none]
-k <-v|-t|-p|-bc|-c>                          other options for KMC, details see https://github.com/refresh-bio/KMC             [default: none]
-g <-l|-m......>                              other options for GenomeScope2, details see https://github.com/tbenavi1/genomescope2.0 [default: none]
-V   								          show version number
-h                                            display this help and exit
EOF
	exit
}

export bg_len=""
export ed_len=""
export step=""
export reads=""
export pl="2"
export fullmode=""
export correct=""
export outputpref="out"
export threads="1"
export mem="12"
export workdir=""
export kmcoption=""
export genomescopeoption=""

if [ $# -lt 1 ]
then
    usage
    exit 1
fi

while getopts ":b:e:s:r:p:fco:t:M:W:k:g:hV" opt; do
	case $opt in
		b)
			export bg_len="$OPTARG"
			;;
		e)
			export ed_len="$OPTARG"
			;;
		s)
			export step="$OPTARG"
			;;
		r)
			export reads="$OPTARG"
			;;
		p)
			export pl="$OPTARG"
			;;
		f)
			export fullmode="1"
			;;
		c)
			export correct="1"
			;;
		o)
			export outputpref="$OPTARG"
			;;
		t)
			export threads="$OPTARG"
			;;
		M)
			export mem="$OPTARG"
			;;
		W)
			export workdir="$OPTARG"
			;;
		k)
			export kmcoption="$OPTARG"
			;;
		g)
			export genomescopeoption="$OPTARG"
			;;
		h)
			usage
			;;
		V)
			version
			;;
		\?) 
			echo "Invalid option -$OPTARG" >&2
			usage
			exit 1
			;;
		:)
			echo "Option -$OPTARG requires arguments." >&2
			usage
			exit 1
			;;
	esac
done

if [ "$bg_len" == "" ]; then
    echo "Option -b needed."
    exit 1
else
    echo "Loop forecasts beginning at K: $bg_len"
fi

if [ "$ed_len" == "" ]; then
    echo "Option -e needed."
    exit 1
else
    echo "Loop forecasts ending at K: $ed_len"
fi

if [ "$step" == "" ]; then
    echo "Option -s needed."
    exit 1
else
    echo "Iteration step K: $step"
fi

if [ "$reads" == "" ]; then
    echo "Option -r needed."
    exit 1
else
    echo "reads file: $reads"
fi

echo "Using $threads threads"

if [ "$workdir" == "" ]; then
	export workdir=$PWD
else
	mkdir -p $workdir
fi

echo "Using \"$workdir\" as work space."

if ! [ -x "$(command -v FastK)" ]; then
    echo 'Error: FastK is not available. details see https://github.com/thegenemyers/FASTK'
    exit 1
fi

if ! [ -x "$(command -v Histex)" ]; then
    echo 'Error: Histex is not available. details see https://github.com/thegenemyers/FASTK'
    exit 1
fi

if ! [ -x "$(command -v genomescope.R)" ]; then
    echo 'Error: genomescope.R is not available. details see https://github.com/tbenavi1/genomescope2.0'
    exit 1
fi



cd $workdir
if [ "$correct" == "1" ]; then
	if ! [ -x "$(command -v hifiasm)" ]; then
		echo 'Error: hifiasm is not available. details see https://github.com/chhylp123/hifiasm'
		exit 1
	else
		echo "Re-correcting HiFi reads"
		hifiasm -o cor --write-ec -e -t $threads ${reads}
		rm cor.*.bin
		reads="$PWD/cor.ec.fa"
	fi
fi

loop_forecast_kmc.sh
report.R
export DFtest=`grep "Dickey-Fuller test:" ${outputpref}.report | awk '{split($0,A,";");print A[2]}'`
mkdir -p ${outputpref}_figure
fig.R


limiK=`grep "Limiting K:" ${outputpref}.report | cut -d " " -f 3 `
cp -r intermediate/k${limiK} ./
mv k${limiK} ${outputpref}_limiV

if [ "$fullmode" != "1" ]; then
	rm -r intermediate
fi






