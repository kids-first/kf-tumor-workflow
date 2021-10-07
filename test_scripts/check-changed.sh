for cwl_file in $(git diff --name-only | grep '.*\.cwl$')
do
  cwltool --validate "$cwl_file"
done
