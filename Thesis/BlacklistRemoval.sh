

mkdir Blacklists

wget http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/mm10-mouse/mm10.blacklist.bed.gz

# https://www.encodeproject.org/annotations/ENCSR636HFF/

MouseBlacklist=/rds/general/user/ah3918/ephemeral/MARREFRESH/chipseq/Blacklists/mm10.blacklist.bed

#Adds "chr" prefix, bed files currently missing

for folders in /rds/general/user/ah3918/ephemeral/MARREFRESH/chipseq/GSE*; do
  for files in $folders/MACS2FINALOUTPUT; do
    for beds in $files/*.bed; do
      bname=`basename $beds`
      awk '{print "chr" $0}' $beds > $files/${bname}.prefixed.bed
    done
  done
done

#Filters out blacklisted regions based on MouseBlacklist (no human datasets left)

for folders in /rds/general/user/ah3918/ephemeral/MARREFRESH/chipseq/GSE*; do
  for files in $folders/MACS2FINALOUTPUT; do
    for beds in $files/*.prefixed.bed; do
      bname=`basename $beds`
      bedtools intersect -v -a $beds -b $MouseBlacklist > $folders/BLACKLISTFILTERED/${bname}.filtered.bed
    done
  done
done


for folders in /rds/general/user/ah3918/ephemeral/MARREFRESH/chipseq/GSE*; do
  for files in $folders/BLACKLISTFILTERED; do
    rm $files/*.filtered
  done
done

for i in
do Pipeline.sh
done

Peakfiles=/rdsgpfs/general/ephemeral/user/ah3918/ephemeral/MARREFRESH/GRANGES/PeakFilesOnly

for folders in /rds/general/user/ah3918/ephemeral/MARREFRESH/chipseq/GSE*; do
  for files in $folders/MACS2FINALOUTPUT; do
    cp $files/*_peaks.xls $Peakfiles
  done
done

ls > SRRList.txt
