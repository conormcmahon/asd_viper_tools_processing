# Iterate through all EMIT .nc files and reprocess to orthorectified .ENVI format

input_dir="/mnt/g/EMIT_MESMA/nc/"
output_dir="/mnt/g/EMIT_MESMA/envi/"
echo "Processing all files in directory ${input_dir} and saving outputs to the directory ${output_dir}"
echo ""
for d in ${input_dir}/* ; do
	echo ""
	echo "Processing EMIT file ${d}"
	python ~/src/emit_utils/emit_utils/reformat.py "${d}" "${output_dir}" --orthorectify
done
