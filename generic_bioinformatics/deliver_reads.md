Transfer data to EPA
Dry run!
aws s3 sync --profile epa ./ s3://edap-ncct-external/HTTr/<CallNumber>/<PlateNumber>/ --exclude "*bspc*" --exclude "*bsnc*" --dryrun | wc -l
- Call number is the number of the bucket (68HER...) - Plate number are the 2 folders
- Plate number i


s3://dmap-stage-ccte-httr/BioSpyder/


Dry run!
aws s3 cp --profile biospyder ./ s3://dmap-stage-ccte-httr/BioSpyder/68HERC22F0496/HiSeq/ --exclude "*bspc*" --exclude "*bsnc*" --dryrun --recursive | wc -l


aws s3 cp --profile biospyder ./ s3://dmap-stage-ccte-httr/BioSpyder/68HERC22F0496/NovaSeq/ --exclude "*bspc*" --exclude "*bsnc*" --dryrun --recursive | wc -l



BIOS2557 - HiSeq
BIOS2837 - NovaSeq











