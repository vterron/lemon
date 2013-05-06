#!/bin/bash

INITIAL_NUMBER=$1;	# will be zero if no argument is provided
PREFIX="ferM_";
PREFIX_LENGTH=${#PREFIX};
EXTENSION=".fits";
NUMBER_OF_DIGITS=4;	


echo "Initial counting number: $INITIAL_NUMBER";
echo "Filename prefix: \"$PREFIX\""
echo "Filename prefix length: $PREFIX_LENGTH"
echo "File extension: \"$EXTENSION\""
echo "Number of digits in the filename:  $NUMBER_OF_DIGITS"

LAST_NUMBER=`ls *$EXTENSION | sort| tail -n 1`;
echo "Last file in destination directory is: $LAST_NUMBER";

for filename in *.fits;
do
	originalName=$filename;
	echo "Filename: $filename";
	filename=`basename $filename $EXTENSION`;
	#echo "Stripping file extension: $filename"; 
	filename=${filename:$PREFIX_LENGTH};
	#echo "Stripping filename prefix: $filename";

	fileNumber=$INITIAL_NUMBER;
	let INITIAL_NUMBER++;
	echo "Assigned file number: $fileNumber"
	
	until [ ${#fileNumber} -ge $NUMBER_OF_DIGITS ];
	do
		fileNumber="0$fileNumber";
	done;

	echo "Adding the necessary number of zeros: $fileNumber"	
	DEST_NAME="$PREFIX$fileNumber$EXTENSION"
	echo "New filename: $DEST_NAME"
	echo "Renaming $originalName to $DEST_NAME"	
	mv $originalName $DEST_NAME

done

